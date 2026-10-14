      subroutine smrt_conversion2swat

	
!!    This subroutine stores discharge rates and nutrient loads into MODFLOW River Cell objects
!!    (these rates/loads are then passed to the SWAT+ Channel objects in the "command" subroutine (sum all receiving hydrographs)

      
      !import variables
      use hru_module
      use hydrograph_module, only:res,sp_ob,sp_ob1,ob,ch_stor,ch_gwex_d
      use basin_module
      use time_module,  only:time
	use water_allocation_module, only:div_volume_total
      use GLOBAL,       only:LENUNI,ITMUNI,NROW,NCOL,NLAY,HNEW,BOTM, !modflow
     &                       DELR,DELC,IBOUND,IUNIT
      use GWFRIVMODULE, only:RIVR,MXRIVR !modflow river package
	use GWFDRNMODULE !modflow drain package
      use GWFUPWMODULE, only:SC2UPW  !modflow upw package
      use GWFLPFMODULE, only:SC2  !modflow lpw package
      USE GWFEVTMODULE, only:EVTvol !modflow evt package
	use GWFRESMODULE, only:BRES,HRES,CRES,BBRES !modflow reservoir package
      
      use smrt_parm !swat-modflow linkage
      use mf_rt_link,   only:rt_active,rt_rivmass,rt_drnmass !MODFLOW-RT3D Linkage
      use rt_global,    only:CNEW !RT3D
      implicit none
      
      !local variables
	integer :: i,j,k,m,n
	integer :: mhru
	integer :: ob_num
	integer :: nriver_cells
	integer :: cell_lay,cell_col,cell_row
	integer :: mf_row,mf_col
	integer :: mf_cellID                 !              |ID of the MODFLOW river cell, in river package
	integer :: riv_cellID                !              |ID of the MODFLOW river cell, in modflow.con
	integer :: riv_count                 !              |number of river cells that are in connection with channels
	integer :: irec											 !              |recall object ID
	real :: basin_area                   !m2						|watershed area
	real :: rate,hriv,criv,rbot,gw_head
	real :: no3_gw_loading,p_gw_loading
	real :: gw_volume,gw_change,gw_conhead
	real :: thickness,cell_area,sy,cell_gw_volume,gw_et
	real :: hru_area,gw_sw_discharge,sw_gw_seepage,gw_recharge,soil_perc
      real :: const
	!canal exchange
	integer :: canl_flow
	!drainage
	integer :: drn_lay,drn_row,drn_col,chan_id
	real :: drn_elev,cell_head,drn_cond,drn_flow
	!reservoir exchange
	integer :: res_id
	real :: HR,CR,HHNEW	
	
	
	write(smrt_log,*) '     converting to swat+ variables'
	
      !get the total area of the basin
      basin_area = bsn%area_tot_ha * 10000. !m2
      
	
	!calculate the water exchange between the aquifer and the channels ----------------------------------------------
      !(only proceed if River package is active)
	if(rt_active == 1) rt_gwsw_load = 0.
      if(IUNIT(4) > 0) then      
      nriver_cells = sp_ob%modflow !number of MODFLOW River Cell objects 
      !loop through all MODFLOW River Cells. If there is a match with a River Cell object, then calculate and
      !store the discharge rate.
      do n=1,nriver_cells
        !determine if the River cell is linked with a SWAT+ sub-basin
        mf_row = RIVR(2,n) !row of MODFLOW river cell
        mf_col = RIVR(3,n) !column of MODFLOW river cell
        mf_cellID = ((mf_row-1)*NCOL) + mf_col !grid ID of MODFLOW river cell
        !loop through all of the River Cell objects
        ob_num = sp_ob1%modflow !object number of first River Cell
        do m=1,nriver_cells
          riv_cellID = ob(ob_num)%props !MODFLOW cell ID of the current river cell object
          !only proceed if there is a match
          if(mf_cellID == riv_cellID) then
            !Use original MODFLOW code to calculate the rate of loss/gain from/to rivers (LENUNI^3/ITMUNI)
            !the below code is borrowed from the MODFLOW subroutine LMT7RIV7 to calculate RATE
            cell_lay = RIVR(1,n) !layer
            cell_row = RIVR(2,n) !row
            cell_col = RIVR(3,n) !column
            hriv = RIVR(4,n) !head in the river
            criv = RIVR(5,n) !river bed conductance
            rbot = RIVR(6,n) !river bottom elevation
            gw_head = HNEW(cell_col,cell_row,cell_lay) !groundwater head in the cell
            !calculate discharge by comparing head in aquifer to head in river
				rate = 0.
            if(gw_head > rbot) then
              rate = (criv*hriv) - (criv*gw_head)
            elseif(gw_head <= rbot) then
              rate = criv * (hriv-rbot)
            endif
				!store for smrt output
				mf_chan_rates(n) = rate
				!update the water volume in the channel
				chan_id = mf_chan_channel(n)
				if(rate > 0) then !channel seepage
				  if(rate > ch_stor(chan_id)%flo) then !can only remove what is there
                rate = ch_stor(chan_id)%flo
              endif
				endif
				ch_stor(chan_id)%flo = ch_stor(chan_id)%flo + (rate*(-1))	  
				!store for channel output
				ch_gwex_d(chan_id)%gwsw = ch_gwex_d(chan_id)%gwsw + (rate*(-1))
            !if rt3d is active, add groundwater nutrient mass loading to the River Cell object
            !(this will be added to the channel objects in the "command" subroutine)
            if(rt_active == 1) then
				  !no3
              no3_gw_loading = rt_rivmass(n,1) / 1000. !kg
					if(rate > 0) then !channel seepage
					  if(no3_gw_loading > ch_stor(chan_id)%no3) then
						  no3_gw_loading = ch_stor(chan_id)%no3
					  endif
					endif
					ch_stor(chan_id)%no3 = ch_stor(chan_id)%no3 + (no3_gw_loading*(-1)) !kg
					rt_gwsw_load(cell_col,cell_row,1) = no3_gw_loading
					!p
              p_gw_loading = rt_rivmass(n,2) / 1000. !kg
					if(rate > 0) then !channel seepage
					  if(p_gw_loading > ch_stor(chan_id)%solp) then
						  p_gw_loading = ch_stor(chan_id)%solp
					  endif
					endif
					ch_stor(chan_id)%solp = ch_stor(chan_id)%solp + (p_gw_loading*(-1)) !kg
					rt_gwsw_load(cell_col,cell_row,2) = p_gw_loading
					!store for channel output
					ch_gwex_d(chan_id)%gwsw_n = ch_gwex_d(chan_id)%gwsw_n + 
     &																	 (no3_gw_loading*(-1))
					ch_gwex_d(chan_id)%gwsw_p = ch_gwex_d(chan_id)%gwsw_p + 
     &                                   (p_gw_loading *(-1))
            endif
            goto 10 !step out of loop - go to the next MODFLOW River Cell
          endif
          ob_num = ob_num + 1 !increment the object number
        enddo !go to the next River Cell object 
 10   enddo  !go to the next MODFLOW River Cell
      endif !check if River Package is active
      
	
	!calculate the water exchange between the canals and the channels -----------------------------------------------
      !(only proceed if River package is active)
	if(rt_active == 1) rt_canl_load = 0.
	if(mf_cnl_flag == 1) then
	  !loop through the canal cells
	  if(IUNIT(4) > 0) then
		  riv_count = sp_ob%modflow !number of MODFLOW River Cell objects 
			do i=1,mf_num_canl
			  riv_count = riv_count + 1
			  cell_lay = RIVR(1,riv_count) !layer
            cell_row = RIVR(2,riv_count) !row
            cell_col = RIVR(3,riv_count) !column
            hriv = RIVR(4,riv_count) !head in the canal
            criv = RIVR(5,riv_count) !canal bed conductance
            rbot = RIVR(6,riv_count) !canal bottom elevation
				gw_head = HNEW(cell_col,cell_row,cell_lay) !groundwater head in the cell
            !calculate discharge by comparing head in aquifer to head in canal
				canl_flow = 0.
            if(gw_head > rbot) then
              canl_flow = (criv*hriv) - (criv*gw_head)
            elseif(gw_head <= rbot) then
              canl_flow = criv * (hriv-rbot)
            endif
				!store for groundwater budget
				mf_canl_rates(i) = canl_flow
				!remove the canal seepage water from the canal water source; if the canal is attached 
				!to a canal diversion, then the seepage water is removed from the volume of diverted water.
				!if not, then the seepage water is removed from the source channel.
				irec = mf_canl_div(i)
				if(irec > 0) then
				  if(div_volume_total(irec) > 0) then
				    div_volume_total(irec) = div_volume_total(irec) - canl_flow
					endif
				else
				  chan_id = mf_canl_channel(i)
					if(chan_id > 0) then
						if(canl_flow > ch_stor(chan_id)%flo) then
							canl_flow = ch_stor(chan_id)%flo
						endif
						ch_stor(chan_id)%flo = ch_stor(chan_id)%flo - canl_flow
					endif
					!store for channel output
				  ch_gwex_d(chan_id)%canl = ch_gwex_d(chan_id)%canl - canl_flow
				endif
            !if rt3d is active, add to nutrient mass in the original channel
            if(rt_active == 1) then
					no3_gw_loading = rt_rivmass(riv_count,1) / 1000. !kg
					p_gw_loading = rt_rivmass(riv_count,2) / 1000. !kg
					if(chan_id > 0) then
                ch_stor(chan_id)%no3 = ch_stor(chan_id)%no3 + 
     &                                 no3_gw_loading
                ch_stor(chan_id)%solp = ch_stor(chan_id)%solp +
     &                                  p_gw_loading
              endif
					rt_canl_load(cell_col,cell_row,1) = no3_gw_loading
					rt_canl_load(cell_col,cell_row,2) = p_gw_loading
            endif
			enddo !go to next canal cell
		endif
	endif !check if gw-canal exchange is active

	
      !calculate the water exchange between the drains and the channels -----------------------------------------------
      !(only proceed if Drain package is active)
	if(mf_drn_flag == 1) then
	  mf_drng_rates = 0.
        if(rt_active == 1) rt_drng_load = 0.
		if(IUNIT(3) > 0) then   
          do i=1,mf_num_drng
	      !cell indices
	      drn_lay = DRAI(1,i)
            drn_row = DRAI(2,i)
            drn_col = DRAI(3,i)
	      !only proceed if modflow cell is active
		    if(IBOUND(drn_col,drn_row,drn_lay) > 0) then
		      drn_elev = DRAI(4,i) !elevation of drain
				  !only proceed if groundwater head is above the drain
				  cell_head = HNEW(drn_col,drn_row,drn_lay)
		      if(cell_head > drn_elev) then
		        drn_cond = DRAI(5,i) !conductance between drain and aquifer
					  drn_flow = drn_cond * (cell_head - drn_elev) !L3/day
					  !update the channel flow volume
					  chan_id = mf_drng_channel(i)
					  ch_stor(chan_id)%flo = ch_stor(chan_id)%flo + drn_flow
						!store for channel output
				    ch_gwex_d(chan_id)%drng = ch_gwex_d(chan_id)%drng + drn_flow
				    !store for groundwater budget
					  mf_drng_rates(i) = drn_flow * (-1) !change to negative (removed from aquifer)
					  !nutrient loading from the drain; add to channel
					  if(rt_active == 1) then
                  no3_gw_loading = rt_drnmass(i,1) / 1000. !kg
						  p_gw_loading = rt_drnmass(i,2) / 1000. !kg
                  ch_stor(chan_id)%no3 = ch_stor(chan_id)%no3 + 
     &                                   no3_gw_loading
                  ch_stor(chan_id)%solp = ch_stor(chan_id)%solp +
     &                                    p_gw_loading
							!store for channel output
                  ch_gwex_d(chan_id)%drng_n = ch_gwex_d(chan_id)%drng_n
     &																			+ no3_gw_loading 
                  ch_gwex_d(chan_id)%drng_p = ch_gwex_d(chan_id)%drng_p
     &																			+ p_gw_loading 
                  !store for groundwater budget
							rt_drng_load(cell_col,cell_row,1) = no3_gw_loading
							rt_drng_load(cell_col,cell_row,2) = p_gw_loading
					  endif
				  endif
		    endif
	    enddo
	  endif
	endif
	
	
	!calculate the water exchange between the cells and the lakes and the channels ----------------------------------
      !(only proceed if Reservoir package is active)
	if(mf_res_flag == 1) then
	  mf_resv_rates = 0.
		if(rt_active == 1) rt_resv_load = 0.
	  if(IUNIT(17) > 0) then   
		  do i=1,mf_num_resv
			  cell_row = mf_resv_row(i)
				cell_col = mf_resv_col(i)
				res_id = mf_resv_reservoir(i)
				rate = 0.
				HR = HRES(i) !reservoir stage
		    if(HR > BRES(cell_col,cell_row)) then !only proceed if there is water in the reservoir
				  CR = CRES(cell_col,cell_row) !conductance (m2/day)
              RBOT = BBRES(cell_col,cell_row) !bottom elevation of bed sediments (m)
              HHNEW = HNEW(cell_col,cell_row,1) !groundwater head in cell (m)
				  if(HHNEW <= RBOT) then 
					  rate = CR * (HR-RBOT) !reservoir --> aquifer
						if(res(res_id)%flo > rate) then
						  res(res_id)%flo = res(res_id)%flo - rate
						endif
						mf_resv_rates(i) = rate 
					else
					  rate = CR * (HHNEW-RBOT) !aquifer --> reservoir
						res(res_id)%flo = res(res_id)%flo + rate
						mf_resv_rates(i) = rate * (-1)
					endif
				endif
			enddo
		endif
	endif
	
	
 100  format(4i6,20f12.4)

      return
      end