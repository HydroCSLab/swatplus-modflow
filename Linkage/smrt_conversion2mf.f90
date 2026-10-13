      subroutine smrt_conversion2mf

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the necessary SWAT+ variables (per hru) into SMRT variables
!!    then converts these into the proper units for MODFLOW
      
!     Import variables
      use hru_module
      use hydrograph_module, only: res,sp_ob,sp_ob1,ob
			use water_body_module, only: res_wat_d
      use basin_module
			use soil_module, only:soil
			use sd_channel_module, only : sd_ch
			use output_landscape_module, only: hwb_d
			use time_module
      use GLOBAL, only:LENUNI,ITMUNI,ncol,nrow,nlay,botm,ibound,DELC,DELR,IUNIT,HNEW !MODFLOW
	    use GWFRCHMODULE, only:RECH
	    use GWFEVTMODULE, only:EVTR,EXDP,SURF
			use GWFUPWMODULE, only:SC2UPW  !MODFLOW
      use GWFLPFMODULE, only:SC2  !MODFLOW
			use GWFRIVMODULE, only:RIVR,NRIVER !MODFLOW
			use GWFWELMODULE, only:NWELLS !modflow well package
			use GWFDRNMODULE, only:NDRAIN !modflow drain package
			use GWFRESMODULE, only:BRES,HRES !modflow reservoir package
      use smrt_parm !SWAT-MODFLOW
			use mf_rt_link, only:rt_active
	    use rt_global, only:cnew
			
      implicit none

	    integer i,j,jj,k,kk,n
      integer cell_row,cell_col,ob_num,ob_cha_num,cha_num,chan_id,res_id
	    integer num_channels         !          |number of channels connected to the river cell
			integer riv_count            !          |number of river cells connected to channels
			integer dum
	    real :: recharge             !mm        |previous day's recharge 
	    real :: hru_rech_volume      !m3        |recharge volume from the hru
	    real :: hru_gwet_volume      !m3        |remaining ET volume from the hru
	    real :: cell_rech_volume     !m3        |recharge given to the cell
	    real :: cell_gwet_volume     !m3        |remaining ET given to the cell
	    real :: chan_depth           !m         |channel depth (weighted)
			real :: cha_depth            !m         |channel depth (single channel)
	    real :: fraction
	    real :: hru_Q                      !m3     |volume transferred from cell to the soil profile
      real :: hru_soilz                  !m      |thickness of HRU soil profile
      real :: vadose_z    							 !m			 |thickness of cell vadose zone
      real :: poly_area                  !m2     |area of cell within the HRU
      real :: tran_mass_no3              !g      |nitrate mass transferred from cell to hru soil
	    real :: tran_mass_p                !g      |phosphorus mass transferred from cell to hru soil
			real :: gw_mass                    !g      |current solute mass in cell
      real :: water_depth(100)           !m      |depth of groundwater in each soil layer
      real :: water_depth_tot            !m      |total depth of groundwater in the soil profile
      real :: sol_thick                  !m      |thickness of soil layer
      real :: layer_fraction             !       |fraction of saturated soil profile within the layer
      real :: layer_transfer             !       |amount of water and solute transferred to the soil layer
      real :: hru_area_m2                !m2     |surface area of the hru
			real :: cell_area                  !m2     |surface area of grid cell
			real :: thickness									 !m      |aquifer saturated thickness
			real :: sy                         !       |aquifer specific yield
			real :: gwvol_avail                !m3     |current groundwater storage in cell
			real :: flow_area                  !m2     |area of water flow between aquifer and channel
	    real :: res_depth									 !m      |depth of water in reservoir
			
	
	
			write(smrt_log,*) '     converting variables to modflow'
			
      !calculate recharge to the water table, based on the percolation from the soil profile --------------------------
      do i=1,sp_ob%hru
        recharge = gw_rech(i)
        gw_rech(i) = 0.
        gw_rech(i) = ((1.-gw_delaye(i))*sepbtm(i)) + (gw_delaye(i)*recharge)
        if (gw_rech(i) < 1.e-6) gw_rech(i) = 0.
      enddo

	    
	    
	    !map recharge to MODFLOW grid cell arrays ------------------------------------------------------------------------
	    RECH = 0.
	    ob_num = sp_ob1%hru  !object number of first HRU
      do k=1,sp_ob%hru
        !recharge volume from hru
	      hru_rech_volume = (gw_rech(k)/1000.) * (ob(ob_num)%area_ha * 10000.) !m * m2 = m3
        do i=1,mf_hru_num_cells(k)
		      cell_row = mf_hru_cell_row(k,i)
			    cell_col = mf_hru_cell_col(k,i)
		      !recharge volume to add to cell
          cell_rech_volume = hru_rech_volume * mf_hru_cell_fract(k,i)
			    RECH(cell_col,cell_row) = RECH(cell_col,cell_row) + cell_rech_volume !m3
        enddo
        ob_num = ob_num + 1
      enddo
	    
			
			
			!map groundwater ET to MODFLOW grid cell arrays --------------------------------------------------------------------
	    if(mf_evt_flag == 1) then
				EVTR = 0.
				ob_num = sp_ob1%hru  !object number of first HRU
				do k=1,sp_ob%hru
					!ET volume from hru
					hru_gwet_volume = (etremain(k)/1000.) * (ob(ob_num)%area_ha * 10000.) !m * m2 = m3
					do i=1,mf_hru_num_cells(k)
						cell_row = mf_hru_cell_row(k,i)
						cell_col = mf_hru_cell_col(k,i)
						!groundwater ET to remove from the cell
						cell_gwet_volume = hru_gwet_volume * mf_hru_cell_fract(k,i)
						EVTR(cell_col,cell_row) = EVTR(cell_col,cell_row) + cell_gwet_volume !m3
					enddo
					ob_num = ob_num + 1
				enddo
			endif
	    
			
	     
	    !calculate groundwater to be transferred to the HRU soil profile ------------------------------------------------
	    RECH_down = RECH !store original recharge values; for solute leaching
			gw_transfer = 0.
			mf_hru_soil = 0.
			if(rt_active == 1) rt_soil_load = 0.
			if(mf_soil_flag == 1) then !flag is turned on; proceed
	    do k=1,sp_ob%hru
			  
	 	    hru_area_m2 = ob(k)%area_ha * 10000. !area (m2)    
        hru_soilz = soil(k)%phys(soil(k)%nly)%d / 1000. !soil thickness (m)
        !loop through the modflow cells connected to the HRU
		    hru_Q = 0.
        do n=1,mf_hru_num_cells(k)
          
          !row and column of the current cell
			    cell_row = mf_hru_cell_row(k,n)
			    cell_col = mf_hru_cell_col(k,n)
          
          !only proceed if the cell is active
          if(ibound(cell_col,cell_row,1) > 0) then
          
					!cell properties
			    cell_area = DELR(cell_col) * DELC(cell_row) !m2
			    if(IUNIT(23).gt.0 .or. IUNIT(1).gt.0) then
            sy = SC2(cell_col,cell_row,1) / cell_area !specific yield from cell storage capacity
          else
            sy = SC2UPW(cell_col,cell_row,1) / cell_area !specific yield from cell storage capacity
          endif
								
          !calculate the current thickness of the vadose zone (m)
          vadose_z = 5000.
					thickness = 0.
					do kk=1,nlay
				    !vadose zone (m)
			      if(HNEW(cell_col,cell_row,kk) < BOTM(cell_col,cell_row,kk-1) & !BOTM(J,cell_row,0) contains the top of layer1 
      					    .and. &
               HNEW(cell_col,cell_row,kk) >= BOTM(cell_col,cell_row,kk)) then 
              vadose_z = BOTM(cell_col,cell_row,0) - HNEW(cell_col,cell_row,kk)
							!for the cell, calculate the saturated thickness (m)
							thickness = HNEW(cell_col,cell_row,kk) - BOTM(cell_col,cell_row,kk)
            endif			
				  enddo !go to next layer
				  gwvol_avail = thickness * cell_area * sy !m3
            
          !if water table is within the soil profile --> calculate groundwater volume (Q) to transfer; then transfer to soil layer
          hru_Q = 0.
          if(vadose_z < hru_soilz) then !water table is within the soil profile
				    !volume of groundwater to transfer to the soil profile
				    poly_area = cell_area * mf_cell_fract(k,n) !area (m2) of cell within HRU    
            hru_Q = (hru_soilz - vadose_z) * poly_area * sy !m3 of groundwater to transfer
					  gw_transfer(cell_col,cell_row) = gw_transfer(cell_col,cell_row) + hru_Q
					
						!solutes
            if (rt_active == 1) then
						  !no3
					    gw_mass = gwvol_avail * cnew(cell_col,cell_row,1,1) !g
					    tran_mass_no3 = hru_Q * cnew(cell_col,cell_row,1,1) !g
					    if(tran_mass_no3 > gw_mass) then
						    tran_mass_no3 = gw_mass
						  endif
              !p
						  gw_mass = gwvol_avail * cnew(cell_col,cell_row,1,2) !g
					    tran_mass_p = hru_Q * cnew(cell_col,cell_row,1,2) !g
					    if(tran_mass_p > gw_mass) then
						    tran_mass_p = gw_mass
						  endif
            endif !end solutes
							
            !determine which HRU soil layers are to receive groundwater
            water_depth = 0.
            water_depth_tot = 0.
            do jj=1,soil(k)%nly
              sol_thick = soil(k)%phys(jj)%thick / 1000. !mm --> m
              if((soil(k)%phys(jj)%d/1000.) > vadose_z) then
                if(jj == 1) then !top layer
                  water_depth(jj) = (soil(k)%phys(jj)%d/1000.) - vadose_z
                else
                  if(vadose_z > (soil(k)%phys(jj-1)%d/1000.)) then !water table is within this layer
                    water_depth(jj) = (soil(k)%phys(jj)%d/1000.) - vadose_z
                  else !water table is above this layer
                    water_depth(jj) = sol_thick
                  endif
                endif
                water_depth_tot = water_depth_tot + water_depth(jj)
              endif  
            enddo !next soil layer
							
					  !transfer the groundwater and the solute mass to the HRU soil layers
            do jj=1,soil(k)%nly
              if(water_depth_tot > 0) then
                layer_fraction = water_depth(jj) / water_depth_tot
              else
                layer_fraction = 0.
              endif
              layer_transfer = (hru_Q*layer_fraction) / hru_area_m2 * 1000. !m3 --> mm
              soil(k)%phys(jj)%st = soil(k)%phys(jj)%st + layer_transfer !mm
              gwsoilq(k) = gwsoilq(k) + layer_transfer !store for hru output
						  if (rt_active == 1) then
						    !no3
						    layer_transfer = (tran_mass_no3*layer_fraction) / 1000. / ob(k)%area_ha !g --> kg/ha  
							  mf_hru_soil(k,jj,1) = mf_hru_soil(k,jj,1) + layer_transfer !kg/ha (mass added to soil profile in nut_nlch)
								rt_soil_load(cell_col,cell_row,1) = (tran_mass_no3*layer_fraction) / 1000.
						    !p
						    layer_transfer = (tran_mass_p*layer_fraction) / 1000. / ob(k)%area_ha !g --> kg/ha  
							  mf_hru_soil(k,jj,2) = mf_hru_soil(k,jj,2) + layer_transfer !kg/ha (mass added to soil profile in nut_solp)
								rt_soil_load(cell_col,cell_row,2) = (tran_mass_p*layer_fraction) / 1000.
              endif
            enddo  

          endif !check for water table in soil profile
          endif !check if cell is active
        enddo !go to next connected cell
				hwb_d(k)%gwsoil = gwsoilq(k)
	    enddo !go to next hru
			RECH = RECH - gw_transfer !remove groundwater (add to RECH array as a negative)
	    endif
	    
			
	    
	    !prepare data for MODFLOW River package (channels) --------------------------------------------------------------
	    ob_num = sp_ob1%modflow !object number of first River Cell
	    do i=1,sp_ob%modflow
        !object and channel number  
        ob_cha_num = ob(ob_num)%obj_out(1)
        chan_id = ob(ob_cha_num)%num
				!store channel stage in RIVR array
				RIVR(4,i) = RIVR(6,i) + sd_ch(chan_id)%chd !channel bed elev + channel depth
				!calculate conductance using current channel width		
				flow_area = sd_ch(chan_id)%chw * mf_rivr_lngt(i) !width * length
				RIVR(5,i) = (flow_area * mf_rivr_hydc(i)) / mf_rivr_thck(i)
        ob_num = ob_num + 1 !increment the object number
	    enddo !go to next river cell
	    NRIVER = sp_ob%modflow
	     
			
			
			!prepare data for MODFLOW River package (canals) ----------------------------------------------------------------
			!turn on cells if the canals have water
			if(mf_cnl_flag == 1) then
				riv_count = sp_ob%modflow
				do i=1,mf_num_canl
					riv_count = riv_count + 1
					if(time%day.ge.mf_canl_beg(i) .and. time%day.le.mf_canl_end(i)) then
						RIVR(5,riv_count) = mf_canl_cond(i)
					else
						RIVR(5,riv_count) = 0.
					endif
				enddo
			endif
			NRIVER = NRIVER + mf_num_canl
			
			
			
			!prepare data for MODFLOW Well package --------------------------------------------------------------------------
			!well data are processed in smrt_well
			if(iunit(2) > 0) then
			  NWELLS = mf_num_well  
			endif
			
			
			
			!prepare data for MODFLOW Drain package -------------------------------------------------------------------------
			!everything is already stored in DRAI - nothing changes during the simulation
			if(mf_drn_flag == 1) then
			  NDRAIN = mf_num_drng
			endif
			
			

			!prepare data for MODFLOW Reservoir package ---------------------------------------------------------------------
			if(IUNIT(17) > 0 .and. mf_res_flag == 1) then
			  do i=1,mf_num_resv
			    !re-compute reservoir stage
				  res_id = mf_resv_reservoir(i)
				  res_depth = res(res_id)%flo / res_wat_d(res_id)%area_ha / 10000. !m
					cell_row = mf_resv_row(i)
					cell_col = mf_resv_col(i)
			    HRES(i) = BRES(cell_col,cell_row) + res_depth
			  enddo
			endif
			
			
			
      return
      end subroutine smrt_conversion2mf
	