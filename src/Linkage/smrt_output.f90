      subroutine smrt_output

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calculates and writes out groundwater budget terms
	
			use hydrograph_module, only:sp_ob1,sp_ob,ob
			use basin_module
			use time_module
			use smrt_parm
			use GLOBAL,       only:NROW,NCOL,NLAY,HNEW,BOTM,DELR,DELC,IBOUND,IUNIT
			use GWFBASMODULE, only:VBVL,HNOFLO
			use GWFUPWMODULE, only:SC2UPW  !MODFLOW
      use GWFLPFMODULE, only:SC2  !MODFLOW
      use GWFEVTMODULE, only:EVTvol
			use GWFWELMODULE, only:NWELLS,WELL !modflow well package
			USE GWFMNW2MODULE,only:MNWMAX,MNW2,MNWINT !modflow mnw2 package
			use GWFRIVMODULE, only:RIVR !modflow
			use GWFDRNMODULE, only:DRAI !modflow
			use mf_rt_link, only:rt_active
			use rt_global,    only:CNEW !RT3D
			
			
			implicit none
			
			integer :: i,j,k,n,dum  	!none					|counters
			integer :: ob_num         !none         |object number
			integer :: cell_row
			integer :: cell_col
			integer :: riv_count      !             |number of river cells in connection with channels
			real :: basin_area        !m2           |watershed area
			!groundwater volume
			real :: gw_volume         !m3           |total groundwater volume in model domain
			real :: thickness         !m            |saturated thickness of MODFLOW cell
			real :: cell_area         !m2           |surface area of MODFLOW cell
			real :: cell_gw_volume    !m3           |groundwater volume in MODFLOW cell
			real :: gw_change         !m3           |total groundwater change in model domain
			real :: sy                !-            !aquifer specific yield
			!groundwater fluxes (domain-wide)
			real :: rech              !m3						|recharge to water table
			real :: soil              !m3           |groundwater transfer to soil profile
			real :: gwet              !m3						|groundwater ET
			real :: gwsw              !m3						|groundwater discharge to channels
			real :: swgw              !m3           |channel seepage to groundwater
			real :: bndr              !m3           |constant head boundary flow
			real :: welq              !m3           |groundwater pumping - Well package
			real :: mnwq							!m3           |groundwater pumping - MNW2 package
			real :: canl							!m3           |groundwater exchange with canals
			real :: drng              !m3           |groundwater discharge to drains
			real :: resv              !m3           |groundwater exchange with reservoirs
			!groundwater nutrient loads (domain-wide)
			real :: rech_no3,rech_p   !kg						|recharge to water table
			real :: soil_no3,soil_p   !kg           |groundwater transfer to soil profile
			real :: gwsw_no3,gwsw_p   !kg						|groundwater discharge to channels
			real :: swgw_no3,swgw_p   !kg           |channel seepage to groundwater
			real :: well_no3,well_p   !kg           |groundwater pumping - Well package
			real :: canl_no3,canl_p   !kg           |groundwater exchange with canals
			real :: drng_no3,drng_p   !kg           |groundwater discharge to drains
			real :: resv_no3,resv_p   !kg           |groundwater exchange with reservoirs
			real :: rctn_no3					!kg           |groundwater denitrification
			!other
			real :: rate              !m3           |groundwater-channel flow rate
			real :: const             !             |constant used to compute temporal averages
			real :: mnw2_rate         !m3           |pumping rate for mnw2 wells
			
			write(smrt_log,*) '     writing output to smrt files'
			
			!get the total area of the basin
      basin_area = bsn%area_tot_ha * 10000. !m2
			
			
			
			!domain-wide values ---------------------------------------------------------------------------------------------
      
			!groundwater volume -------------------------------------------------------------------------
      gw_volume = 0.
      do i=1,NROW
        do j=1,NCOL
          do k=1,NLAY
            if(IBOUND(j,i,k).ne.0) then
            !Determine saturated thickness of the current layer
            if(HNEW(j,i,k) < BOTM(j,i,k-1) .and. HNEW(j,i,k) > BOTM(j,i,k)) then 
              thickness = HNEW(j,i,k) - BOTM(j,i,k)
            else
              thickness = BOTM(j,i,k-1) - BOTM(j,i,k)
            endif
            cell_area = DELR(j) * DELC(i)
            if(IUNIT(23).gt.0) then
              sy = SC2(j,i,k) / cell_area !specific yield from cell storage capacity
            else
              sy = SC2UPW(j,i,k) / cell_area !specific yield from cell storage capacity
            endif
            cell_gw_volume = thickness * cell_area * sy
            gw_volume = gw_volume + cell_gw_volume
            endif
          enddo
        enddo
      enddo
      gw_volume = gw_volume / basin_area !convert to m
      gw_volume = gw_volume * 1000. !convert m to mm of water
			
			!change in groundwater volume (m3) ----------------------------------------------------------
      gw_change = (VBVL(3,1) - VBVL(4,1)) * (-1)
      gw_change = gw_change / basin_area !convert to m
      gw_change = gw_change * 1000. !convert m to mm of water
			
      !recharge to water table --------------------------------------------------------------------
			rech = 0.
			rech_no3 = 0.
			rech_p = 0.
			if(IUNIT(8) > 0) then
			  !flow
				do i=1,NROW
					do j=1,NCOL
						rech = rech + RECH_down(j,i)
					enddo
				enddo
				rech = rech / basin_area !convert to m
				rech = rech * 1000. !convert m to mm of water
				!totals
				smrt_rech_tot_mo = smrt_rech_tot_mo + RECH_down
				smrt_rech_tot_yr = smrt_rech_tot_yr + RECH_down
				!loads
				if(rt_active == 1) then
				  do i=1,NROW
						do j=1,NCOL
							rech_no3 = rech_no3 + rt_rech_load(j,i,1) !kg
							rech_p = rech_p + rt_rech_load(j,i,2) !kg
							!write(28944,*) rt_rech_load(j,i,1),rt_rech_load(j,i,2)
						enddo
					enddo
					!totals
				  smrt_rt_rech_tot_mo = smrt_rt_rech_tot_mo + rt_rech_load
					smrt_rt_rech_tot_yr = smrt_rt_rech_tot_yr + rt_rech_load
				endif
			endif
			
      !groundwater transfer to soil profile -------------------------------------------------------
      soil = 0.
			soil_no3 = 0.
			soil_p = 0.
      if(IUNIT(8) > 0) then
				do i=1,NROW
					do j=1,NCOL
						soil = soil + (gw_transfer(j,i)*(-1))
					enddo
				enddo
				soil = soil / basin_area !convert to m
				soil = soil * 1000. !convert m to mm of water
				!totals
				smrt_soil_tot_mo = smrt_soil_tot_mo + (gw_transfer*(-1))
				smrt_soil_tot_yr = smrt_soil_tot_yr + (gw_transfer*(-1))
				if(rt_active == 1) then
				  do i=1,NROW
						do j=1,NCOL
							soil_no3 = soil_no3 + (rt_soil_load(j,i,1)*(-1)) !kg
							soil_p = soil_p + (rt_soil_load(j,i,2)*(-1)) !kg
						enddo
					enddo
					!totals
				  smrt_rt_soil_tot_mo = smrt_rt_soil_tot_mo + (rt_soil_load*(-1))
					smrt_rt_soil_tot_yr = smrt_rt_soil_tot_yr + (rt_soil_load*(-1))
				endif
      endif
			
      !total groundwater removed via ET (EVT package) ---------------------------------------------
      gwet = 0.
			if(IUNIT(5) > 0 .and. mf_evt_flag == 1) then
				do i=1,NROW
					do j=1,NCOL
						do k=1,NLAY
							gwet = gwet + EVTvol(j,i,k)
							smrt_gwet_tot_mo(j,i) = smrt_gwet_tot_mo(j,i) + EVTvol(j,i,k)
							smrt_gwet_tot_yr(j,i) = smrt_gwet_tot_yr(j,i) + EVTvol(j,i,k)
						enddo
					enddo
				enddo
				gwet = gwet / basin_area !convert to m
				gwet = gwet * 1000. !convert m to mm of water
			endif
			
			!groundwater to channels; channels to groundwater -------------------------------------------
      ob_num = sp_ob1%modflow
      gwsw = 0.
      swgw = 0.
			gwsw_no3 = 0.
			gwsw_p = 0.
			swgw_no3 = 0.
			swgw_p = 0.
			if(IUNIT(4) > 0) then
				do i=1,sp_ob%modflow
					rate = mf_chan_rates(i) !flow rate
					if(rate < 0) then !discharge
						gwsw = gwsw + rate
					else !seepage
						swgw = swgw + rate
					endif
					cell_row = RIVR(2,i) !row
          cell_col = RIVR(3,i) !column
					smrt_gwsw_tot_mo(cell_col,cell_row) = smrt_gwsw_tot_mo(cell_col,cell_row) + rate
					smrt_gwsw_tot_yr(cell_col,cell_row) = smrt_gwsw_tot_yr(cell_col,cell_row) + rate
					ob_num = ob_num + 1
				enddo
				gwsw = gwsw / basin_area !convert m3 to m
				gwsw = gwsw * 1000. !m to mm of water
				swgw = swgw/ basin_area
				swgw = swgw * 1000.
				if(rt_active == 1) then
				  ob_num = sp_ob1%modflow
				  do i=1,sp_ob%modflow
					  cell_row = RIVR(2,i) !row
						cell_col = RIVR(3,i) !column
						rate = ob(ob_num)%hd(1)%flo !flow rate
						if(rate > 0) then !discharge
							gwsw_no3 = gwsw_no3 + rt_gwsw_load(cell_col,cell_row,1)
							gwsw_p = gwsw_p + rt_gwsw_load(cell_col,cell_row,2)
						else !seepage
						  swgw_no3 = swgw_no3 + rt_gwsw_load(cell_col,cell_row,1)
							swgw_p = swgw_p + rt_gwsw_load(cell_col,cell_row,2)
						endif
						ob_num = ob_num + 1
					enddo
				  !totals
				  smrt_rt_gwsw_tot_mo = smrt_rt_gwsw_tot_mo + rt_gwsw_load
					smrt_rt_gwsw_tot_yr = smrt_rt_gwsw_tot_yr + rt_gwsw_load
				endif				
			endif
			
      !total groundwater added/removed via constant head boundaries -------------------------------
      bndr = (VBVL(3,2) - VBVL(4,2))
      bndr = bndr / basin_area !convert to m
      bndr = bndr * 1000. !convert m to mm of water

      !total groundwater removed via pumping (WELL package) ---------------------------------------
			welq = 0.
			well_no3 = 0.
			well_p = 0.
			if(IUNIT(2) > 0) then
			  do n=1,mf_num_well
			    welq = welq + WELL(4,n)
					cell_row = WELL(2,n)
          cell_col = WELL(3,n)
					smrt_well_tot_mo(cell_col,cell_row) = smrt_well_tot_mo(cell_col,cell_row) + WELL(4,n)
					smrt_well_tot_yr(cell_col,cell_row) = smrt_well_tot_yr(cell_col,cell_row) + WELL(4,n)
			  enddo
				!totals
        welq = welq / basin_area !convert to m
        welq = welq * 1000. !convert m to mm of water
				if(rt_active == 1) then
				  do i=1,NROW
						do j=1,NCOL
							well_no3 = well_no3 + (rt_well_load(j,i,1)*-1) !kg
							well_p = well_p + (rt_well_load(j,i,2)*-1) !kg
						enddo
					enddo  
					!totals
				  smrt_rt_well_tot_mo = smrt_rt_well_tot_mo + rt_well_load
					smrt_rt_well_tot_yr = smrt_rt_well_tot_yr + rt_well_load
				endif
			endif
	
			!total groundwater removed via pumping (MNW2 package) ---------------------------------------
	    mnwq = 0.
			if(IUNIT(50) > 0) then
			  do n=1,MNWMAX
			    mnw2_rate = MNW2(18,n) !actual pumping rate
				  mnwq = mnwq + mnw2_rate 
					cell_row = MNWINT(3,n)
          cell_col = MNWINT(4,n)
					smrt_mnw2_tot_mo(cell_col,cell_row) = smrt_mnw2_tot_mo(cell_col,cell_row) + mnw2_rate
					smrt_mnw2_tot_yr(cell_col,cell_row) = smrt_mnw2_tot_yr(cell_col,cell_row) + mnw2_rate
			  enddo
				!totals
        mnwq = mnwq / basin_area !convert to m
        mnwq = mnwq * 1000. !convert m to mm of water
			endif
			
			!total groundwater exchanged with canals (RIV package) --------------------------------------
			canl = 0.
			canl_no3 = 0.
			canl_p = 0.
			if(mf_cnl_flag == 1) then
				if(IUNIT(4) > 0) then
					riv_count = sp_ob%modflow
				  do n=1,mf_num_canl
						riv_count = riv_count + 1
					  canl = canl + mf_canl_rates(n)
						cell_row = RIVR(2,riv_count) !row
						cell_col = RIVR(3,riv_count) !column
						smrt_canl_tot_mo(cell_col,cell_row) = smrt_canl_tot_mo(cell_col,cell_row) + mf_canl_rates(n)
						smrt_canl_tot_yr(cell_col,cell_row) = smrt_canl_tot_yr(cell_col,cell_row) + mf_canl_rates(n)
					enddo
					canl = canl / basin_area !convert to m
					canl = canl * 1000. !convert m to mm of water
				endif
				if(rt_active == 1) then
				  do i=1,NROW
						do j=1,NCOL
							canl_no3 = canl_no3 + rt_canl_load(j,i,1) !kg
							canl_p = canl_p + rt_canl_load(j,i,2) !kg
						enddo
					enddo 
				  !totals
				  smrt_rt_canl_tot_mo = smrt_rt_canl_tot_mo + rt_canl_load
					smrt_rt_canl_tot_yr = smrt_rt_canl_tot_yr + rt_canl_load
				endif
			endif
			
			!total groundwater removed via drains (DRN package) -----------------------------------------
			drng = 0.
			drng_no3 = 0.
			drng_p = 0.
			if(mf_drn_flag == 1) then
				if(IUNIT(3) > 0) then
					do n=1,mf_num_drng
						drng = drng + mf_drng_rates(n)
						cell_row = DRAI(2,n)
						cell_col = DRAI(3,n)
						smrt_drng_tot_mo(cell_col,cell_row) = smrt_drng_tot_mo(cell_col,cell_row) + mf_drng_rates(n)
						smrt_drng_tot_yr(cell_col,cell_row) = smrt_drng_tot_yr(cell_col,cell_row) + mf_drng_rates(n)
					enddo
					drng = drng / basin_area !convert to m
					drng = drng * 1000. !convert m to mm of water
				endif
				if(rt_active == 1) then
				  !totals
				  do i=1,NROW
						do j=1,NCOL
							drng_no3 = drng_no3 + rt_drng_load(j,i,1) !kg
							drng_p = drng_p + rt_drng_load(j,i,2) !kg
						enddo
					enddo 
				  smrt_rt_drng_tot_mo = smrt_rt_drng_tot_mo + rt_drng_load
					smrt_rt_drng_tot_yr = smrt_rt_drng_tot_yr + rt_drng_load
				endif
			endif
			
			!total groundwater exchange with reservoirs (RES package) -----------------------------------
			resv = 0.
			resv_no3 = 0.
			resv_p = 0.
			if(mf_res_flag == 1) then
				if(IUNIT(17) > 0) then
					do n=1,mf_num_resv
						resv = resv + mf_resv_rates(n)
						cell_row = mf_resv_row(n)
				    cell_col = mf_resv_col(n)
						smrt_resv_tot_mo(cell_col,cell_row) = smrt_resv_tot_mo(cell_col,cell_row) + mf_resv_rates(n)
						smrt_resv_tot_yr(cell_col,cell_row) = smrt_resv_tot_yr(cell_col,cell_row) + mf_resv_rates(n)
					enddo
					resv = resv / basin_area !convert to m
					resv = resv * 1000. !convert m to mm of water
					if(rt_active == 1) then
					  do i=1,NROW
							do j=1,NCOL
								resv_no3 = resv_no3 + rt_resv_load(j,i,1) !kg
								resv_p = resv_p + rt_resv_load(j,i,2) !kg
							enddo
						enddo 
					  !totals
						smrt_rt_resv_tot_mo = smrt_rt_resv_tot_mo + rt_resv_load
						smrt_rt_resv_tot_yr = smrt_rt_resv_tot_yr + rt_resv_load
					endif
				endif
			endif
			
			!total groundwater denitrification (RT3D) ---------------------------------------------------
			rctn_no3 = 0.
			if(rt_active == 1) then
			  do i=1,NROW
					do j=1,NCOL
					  do k=1,NLAY
						  rctn_no3 = rctn_no3 + (rt_rctn_load(j,i,k)*-1) !kg
						enddo
					enddo
				enddo  
			endif
			
			
      !write to water balance files -----------------------------------------------------------------------------------
      
			!daily
			write(out_gwbal_dy,120) time%day,time%mo,time%day_mo,time%yrc, &
			                        gw_volume,gw_change, &
															rech,soil,gwet,gwsw,swgw,bndr,welq,mnwq,canl,drng,resv
											
      !monthly
      gw_volume_mo = gw_volume_mo + gw_volume
      gw_change_mo = gw_change_mo + gw_change
			rech_mo = rech_mo + rech
			soil_mo = soil_mo + soil
			gwet_mo = gwet_mo + gwet
			gwsw_mo = gwsw_mo + gwsw
			swgw_mo = swgw_mo + swgw
			bndr_mo = bndr_mo + bndr
			well_mo = well_mo + welq
			mnw2_mo = mnw2_mo + mnwq
			canl_mo = canl_mo + canl
			drng_mo = drng_mo + drng
			resv_mo = resv_mo + resv
      if(time%end_mo == 1) then
        const = float(ndays(time%mo + 1) - ndays(time%mo))
        gw_volume_mo = gw_volume_mo / const
        write(out_gwbal_mo,120) time%day,time%mo,time%day_mo,time%yrc, &
			                          gw_volume_mo,gw_change_mo, &
															  rech_mo,soil_mo,gwet_mo,gwsw_mo,swgw_mo,bndr_mo,well_mo,mnw2_mo,canl_mo,drng_mo,resv_mo
        !zero out for next month
        gw_volume_mo = 0.
        gw_change_mo = 0.
			  rech_mo = 0.
			  soil_mo = 0.
			  gwet_mo = 0.
			  gwsw_mo = 0.
			  swgw_mo = 0.
			  bndr_mo = 0.
			  well_mo = 0.
				mnw2_mo = 0.
				canl_mo = 0.
			  drng_mo = 0.
			  resv_mo = 0.
      endif

      !yearly
      gw_volume_yr = gw_volume_yr + gw_volume
      gw_change_yr = gw_change_yr + gw_change
			rech_yr = rech_yr + rech
			soil_yr = soil_yr + soil
			gwet_yr = gwet_yr + gwet
			gwsw_yr = gwsw_yr + gwsw
			swgw_yr = swgw_yr + swgw
			bndr_yr = bndr_yr + bndr
			well_yr = well_yr + welq
			mnw2_yr = mnw2_yr + mnwq
			canl_yr = canl_yr + canl
			drng_yr = drng_yr + drng
			resv_yr = resv_yr + resv
      if (time%end_yr == 1) then
        const = time%day_end_yr
        gw_volume_yr = gw_volume_yr / const
        write(out_gwbal_yr,120) time%day,time%mo,time%day_mo,time%yrc, &
			                          gw_volume_yr,gw_change_yr, &
															  rech_yr,soil_yr,gwet_yr,gwsw_yr,swgw_yr,bndr_yr,well_yr,mnw2_yr,canl_yr,drng_yr,resv_yr										
        !zero out for next year
        gw_volume_yr = 0.
        gw_change_yr = 0.
			  rech_yr = 0.
			  soil_yr = 0.
			  gwet_yr = 0.
			  gwsw_yr = 0.
			  swgw_yr = 0.
			  bndr_yr = 0.
			  well_yr = 0.
				mnw2_yr = 0.
				canl_yr = 0.
			  drng_yr = 0.
			  resv_yr = 0.
      endif

      !average annual
      gw_volume_aa = gw_volume_aa + gw_volume
      gw_change_aa = gw_change_aa + gw_change
			rech_aa = rech_aa + rech
			soil_aa = soil_aa + soil
			gwet_aa = gwet_aa + gwet
			gwsw_aa = gwsw_aa + gwsw
			swgw_aa = swgw_aa + swgw
			bndr_aa = bndr_aa + bndr
			well_aa = well_aa + welq
			mnw2_aa = mnw2_aa + mnwq
			canl_aa = canl_aa + canl
			drng_aa = drng_aa + drng
			resv_aa = resv_aa + resv
			if (time%end_sim == 1) then
        gw_volume_aa = gw_volume_aa / time%nbyr
        write(out_gwbal_aa,120) time%day,time%mo,time%day_mo,time%yrc, &
			                          gw_volume_aa,gw_change_aa, &
															  rech_aa,soil_aa,gwet_aa,gwsw_aa,swgw_aa,bndr_aa,well_aa,mnw2_aa,canl_aa,drng_aa,resv_aa	
      endif

			
			
			
			!groundwater head -----------------------------------------------------------------------------------------------
			
			!totals
			do k=1,nlay
			  do i=1,nrow
				  do j=1,ncol
					  if(ibound(j,i,k).ne.0) then
						  smrt_head_tot_mo(j,i,k) = smrt_head_tot_mo(j,i,k) + HNEW(j,i,k)
							smrt_head_tot_yr(j,i,k) = smrt_head_tot_yr(j,i,k) + HNEW(j,i,k)
						else
						  smrt_head_tot_mo(j,i,k) = -999
							smrt_head_tot_yr(j,i,k) = -999
						endif
					enddo
				enddo
			enddo

			!monthly output
			if (time%end_mo == 1) then
			  do k=1,nlay
					do i=1,nrow
						do j=1,ncol
							if(ibound(j,i,k).ne.0) then
								smrt_head_tot_mo(j,i,k) = smrt_head_tot_mo(j,i,k) / time%day_mo 
							else
								smrt_head_tot_mo(j,i,k) = HNOFLO
							endif
						enddo
					enddo
				enddo
			  write(mf_head_mo,*), time%yrc,time%mo
			  do k=1,NLAY
          write(mf_head_mo,*),'layer:',k
          do i=1,NROW
            write(mf_head_mo,101) (smrt_head_tot_mo(j,i,k),j=1,NCOL)
          enddo
        enddo
				write(mf_head_mo,*)
				smrt_head_tot_mo = 0.
			endif
			
			!yearly output
			if(time%end_yr == 1) then
			  do k=1,nlay
					do i=1,nrow
						do j=1,ncol
							if(ibound(j,i,k).ne.0) then
								smrt_head_tot_yr(j,i,k) = smrt_head_tot_yr(j,i,k) / time%day 
							else
								smrt_head_tot_yr(j,i,k) = HNOFLO
							endif
						enddo
					enddo
				enddo
			  write(mf_head_yr,*), time%yrc
			  do k=1,NLAY
          write(mf_head_yr,*),'layer:',k
          do i=1,NROW
            write(mf_head_yr,101) (smrt_head_tot_yr(j,i,k),j=1,NCOL)
          enddo
        enddo
				write(mf_head_yr,*)
				smrt_head_tot_yr = 0.
			endif
			
			
			!groundwater fluxes ---------------------------------------------------------------------------------------------
	    !monthly output
			if (time%end_mo == 1) then
			  !recharge
			  if(IUNIT(8) > 0) then
				  smrt_rech_tot_mo = smrt_rech_tot_mo / time%day_mo
					write(mf_rech_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_rech_mo,101) (smrt_rech_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_rech_mo,*)
					smrt_rech_tot_mo = 0.
				endif
				!gw-->soil
				if(IUNIT(8) > 0) then
				  smrt_soil_tot_mo = smrt_soil_tot_mo / time%day_mo
					write(mf_soil_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_soil_mo,101) (smrt_soil_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_soil_mo,*)
					smrt_soil_tot_mo = 0.
				endif
				!gwet
				if(IUNIT(5) > 0) then
				  smrt_gwet_tot_mo = smrt_gwet_tot_mo / time%day_mo
					write(mf_gwet_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_gwet_mo,101) (smrt_gwet_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_gwet_mo,*)
					smrt_gwet_tot_mo = 0.
				endif
				!gwsw
				if(IUNIT(4) > 0) then
				  smrt_gwsw_tot_mo = smrt_gwsw_tot_mo / time%day_mo
					write(mf_gwsw_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_gwsw_mo,101) (smrt_gwsw_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_gwsw_mo,*)
					smrt_gwsw_tot_mo = 0.
				endif
				!pumping - well package
				if(IUNIT(2) > 0) then
				  smrt_well_tot_mo = smrt_well_tot_mo / time%day_mo
					write(mf_well_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_well_mo,101) (smrt_well_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_well_mo,*)
					smrt_well_tot_mo = 0.
				endif
				!pumping - mnw2 package
				if(IUNIT(50) > 0) then
				  smrt_mnw2_tot_mo = smrt_mnw2_tot_mo / time%day_mo
					write(mf_mnw2_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_mnw2_mo,101) (smrt_mnw2_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_mnw2_mo,*)
					smrt_mnw2_tot_mo = 0.
				endif
				!canal
				if(IUNIT(4) > 0) then
				  smrt_canl_tot_mo = smrt_canl_tot_mo / time%day_mo
					write(mf_canl_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_canl_mo,101) (smrt_canl_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_canl_mo,*)
					smrt_canl_tot_mo = 0.
				endif
				!drainage
				if(IUNIT(3) > 0) then
				  smrt_drng_tot_mo = smrt_drng_tot_mo / time%day_mo
					write(mf_drng_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_drng_mo,101) (smrt_drng_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_drng_mo,*)
					smrt_drng_tot_mo = 0.
				endif
				!reservoir
			  if(IUNIT(17) > 0) then
				  smrt_resv_tot_mo = smrt_resv_tot_mo / time%day_mo
					write(mf_resv_mo,*), time%yrc,time%mo
          do i=1,NROW
            write(mf_resv_mo,101) (smrt_resv_tot_mo(j,i),j=1,NCOL)
          enddo
					write(mf_resv_mo,*)
					smrt_resv_tot_mo = 0.
				endif
			endif
			
			!yearly output
			if(time%end_yr == 1) then
			!recharge
			  if(IUNIT(8) > 0) then
				  smrt_rech_tot_yr = smrt_rech_tot_yr / time%day
					write(mf_rech_yr,*), time%yrc
          do i=1,NROW
            write(mf_rech_yr,101) (smrt_rech_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_rech_yr,*)
					smrt_rech_tot_yr = 0.
				endif
				!gw-->soil
				if(IUNIT(8) > 0) then
				  smrt_soil_tot_yr = smrt_soil_tot_yr / time%day
					write(mf_soil_yr,*), time%yrc
          do i=1,NROW
            write(mf_soil_yr,101) (smrt_soil_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_soil_yr,*)
					smrt_soil_tot_yr = 0.
				endif
				!gwet
				if(IUNIT(5) > 0) then
				  smrt_gwet_tot_yr = smrt_gwet_tot_yr / time%day
					write(mf_gwet_yr,*), time%yrc
          do i=1,NROW
            write(mf_gwet_yr,101) (smrt_gwet_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_gwet_yr,*)
					smrt_gwet_tot_yr = 0.
				endif
				!gwsw
				if(IUNIT(4) > 0) then
				  smrt_gwsw_tot_yr = smrt_gwsw_tot_yr / time%day
					write(mf_gwsw_yr,*), time%yrc
          do i=1,NROW
            write(mf_gwsw_yr,101) (smrt_gwsw_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_gwsw_yr,*)
					smrt_gwsw_tot_yr = 0.
				endif
				!pumping - well package
				if(IUNIT(2) > 0) then
				  smrt_well_tot_yr = smrt_well_tot_yr / time%day
					write(mf_well_yr,*), time%yrc
          do i=1,NROW
            write(mf_well_yr,101) (smrt_well_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_well_yr,*)
					smrt_well_tot_yr = 0.
				endif
				!pumping - mnw2 package
				if(IUNIT(50) > 0) then
				  smrt_mnw2_tot_yr = smrt_mnw2_tot_yr / time%day
					write(mf_mnw2_yr,*), time%yrc
          do i=1,NROW
            write(mf_mnw2_yr,101) (smrt_mnw2_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_mnw2_yr,*)
					smrt_mnw2_tot_yr = 0.
				endif
				!canal
				if(IUNIT(4) > 0) then
				  smrt_canl_tot_yr = smrt_canl_tot_yr / time%day
					write(mf_canl_yr,*), time%yrc
          do i=1,NROW
            write(mf_canl_yr,101) (smrt_canl_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_canl_yr,*)
					smrt_canl_tot_yr = 0.
				endif
				!drainage
				if(IUNIT(3) > 0) then
				  smrt_drng_tot_yr = smrt_drng_tot_yr / time%day
					write(mf_drng_yr,*), time%yrc
          do i=1,NROW
            write(mf_drng_yr,101) (smrt_drng_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_drng_yr,*)
					smrt_drng_tot_yr = 0.
				endif
				!reservoir
			  if(IUNIT(17) > 0) then
				  smrt_resv_tot_yr = smrt_resv_tot_yr / time%day
					write(mf_resv_yr,*), time%yrc
          do i=1,NROW
            write(mf_resv_yr,101) (smrt_resv_tot_yr(j,i),j=1,NCOL)
          enddo
					write(mf_resv_yr,*)
					smrt_resv_tot_yr = 0.
				endif  
			endif
			
			
			
			!groundwater nutrient output-------------------------------------------------------------------------------------
			if(rt_active == 1) then
			
			!mass balance files
			!no3-n
			!daily
			write(out_rtbal_dy_n,121) time%day,time%mo,time%day_mo,time%yrc, &
														  	rech_no3,soil_no3,gwsw_no3,swgw_no3,well_no3,canl_no3,drng_no3,resv_no3,rctn_no3		
      !monthly
			rech_mo_no3 = rech_mo_no3 + rech_no3
			soil_mo_no3 = soil_mo_no3 + soil_no3
			gwsw_mo_no3 = gwsw_mo_no3 + gwsw_no3
			swgw_mo_no3 = swgw_mo_no3 + swgw_no3
			well_mo_no3 = well_mo_no3 + well_no3
			canl_mo_no3 = canl_mo_no3 + canl_no3
			drng_mo_no3 = drng_mo_no3 + drng_no3
			resv_mo_no3 = resv_mo_no3 + resv_no3
			rctn_mo_no3 = rctn_mo_no3 + rctn_no3
      if(time%end_mo == 1) then
        write(out_rtbal_mo_n,121) time%day,time%mo,time%day_mo,time%yrc, &
															    rech_mo_no3,soil_mo_no3,gwsw_mo_no3,swgw_mo_no3,well_mo_no3,canl_mo_no3,drng_mo_no3,resv_mo_no3,rctn_mo_no3
        !zero out for next month
			  rech_mo_no3 = 0.
			  soil_mo_no3 = 0.
			  gwsw_mo_no3 = 0.
			  swgw_mo_no3 = 0.
			  well_mo_no3 = 0.
				canl_mo_no3 = 0.
			  drng_mo_no3 = 0.
			  resv_mo_no3 = 0.
				rctn_mo_no3 = 0.
      endif
      !yearly
			rech_yr_no3 = rech_yr_no3 + rech_no3
			soil_yr_no3 = soil_yr_no3 + soil_no3
			gwsw_yr_no3 = gwsw_yr_no3 + gwsw_no3
			swgw_yr_no3 = swgw_yr_no3 + swgw_no3
			well_yr_no3 = well_yr_no3 + well_no3
			canl_yr_no3 = canl_yr_no3 + canl_no3
			drng_yr_no3 = drng_yr_no3 + drng_no3
			resv_yr_no3 = resv_yr_no3 + resv_no3
			rctn_yr_no3 = rctn_yr_no3 + rctn_no3
      if (time%end_yr == 1) then
        const = time%day_end_yr
        write(out_rtbal_yr_n,121) time%day,time%mo,time%day_mo,time%yrc, &
															    rech_yr_no3,soil_yr_no3,gwsw_yr_no3,swgw_yr_no3,well_yr_no3,canl_yr_no3,drng_yr_no3,resv_yr_no3,rctn_yr_no3										
        !zero out for next year
			  rech_yr_no3 = 0.
			  soil_yr_no3 = 0.
			  gwsw_yr_no3 = 0.
			  swgw_yr_no3 = 0.
			  well_yr_no3 = 0.
				canl_yr_no3 = 0.
			  drng_yr_no3 = 0.
			  resv_yr_no3 = 0.
				rctn_yr_no3 = 0.
      endif
      !average annual
			rech_aa_no3 = rech_aa_no3 + rech_no3
			soil_aa_no3 = soil_aa_no3 + soil_no3
			gwsw_aa_no3 = gwsw_aa_no3 + gwsw_no3
			swgw_aa_no3 = swgw_aa_no3 + swgw_no3
			well_aa_no3 = well_aa_no3 + well_no3
			canl_aa_no3 = canl_aa_no3 + canl_no3
			drng_aa_no3 = drng_aa_no3 + drng_no3
			resv_aa_no3 = resv_aa_no3 + resv_no3
			rctn_aa_no3 = rctn_aa_no3 + rctn_no3
			if (time%end_sim == 1) then
        write(out_rtbal_aa_n,121) time%day,time%mo,time%day_mo,time%yrc, &
															    rech_aa_no3,soil_aa_no3,gwsw_aa_no3,swgw_aa_no3,well_aa_no3,canl_aa_no3,drng_aa_no3,resv_aa_no3,rctn_aa_no3
      endif
			!p
			!daily
			write(out_rtbal_dy_p,121) time%day,time%mo,time%day_mo,time%yrc, &
														  	rech_p,soil_p,gwsw_p,swgw_p,well_p,canl_p,drng_p,resv_p						
      !monthly
			rech_mo_p = rech_mo_p + rech_p
			soil_mo_p = soil_mo_p + soil_p
			gwsw_mo_p = gwsw_mo_p + gwsw_p
			swgw_mo_p = swgw_mo_p + swgw_p
			well_mo_p = well_mo_p + well_p
			canl_mo_p = canl_mo_p + canl_p
			drng_mo_p = drng_mo_p + drng_p
			resv_mo_p = resv_mo_p + resv_p
      if(time%end_mo == 1) then
        write(out_rtbal_mo_p,121) time%day,time%mo,time%day_mo,time%yrc, &
															    rech_mo_p,soil_mo_p,gwsw_mo_p,swgw_mo_p,well_mo_p,canl_mo_p,drng_mo_p,resv_mo_p
        !zero out for next month
			  rech_mo_p = 0.
			  soil_mo_p = 0.
			  gwsw_mo_p = 0.
			  swgw_mo_p = 0.
			  well_mo_p = 0.
				canl_mo_p = 0.
			  drng_mo_p = 0.
			  resv_mo_p = 0.
      endif
      !yearly
			rech_yr_p = rech_yr_p + rech_p
			soil_yr_p = soil_yr_p + soil_p
			gwsw_yr_p = gwsw_yr_p + gwsw_p
			swgw_yr_p = swgw_yr_p + swgw_p
			well_yr_p = well_yr_p + well_p
			canl_yr_p = canl_yr_p + canl_p
			drng_yr_p = drng_yr_p + drng_p
			resv_yr_p = resv_yr_p + resv_p
      if (time%end_yr == 1) then
        const = time%day_end_yr
        write(out_rtbal_yr_p,121) time%day,time%mo,time%day_mo,time%yrc, &
															    rech_yr_p,soil_yr_p,gwsw_yr_p,swgw_yr_p,well_yr_p,canl_yr_p,drng_yr_p,resv_yr_p										
        !zero out for next year
			  rech_yr_p = 0.
			  soil_yr_p = 0.
			  gwsw_yr_p = 0.
			  swgw_yr_p = 0.
			  well_yr_p = 0.
				canl_yr_p = 0.
			  drng_yr_p = 0.
			  resv_yr_p = 0.
      endif
      !average annual
			rech_aa_p = rech_aa_p + rech_p
			soil_aa_p = soil_aa_p + soil_p
			gwsw_aa_p = gwsw_aa_p + gwsw_p
			swgw_aa_p = swgw_aa_p + swgw_p
			well_aa_p = well_aa_p + well_p
			canl_aa_p = canl_aa_p + canl_p
			drng_aa_p = drng_aa_p + drng_p
			resv_aa_p = resv_aa_p + resv_p
			if (time%end_sim == 1) then
        write(out_rtbal_aa_p,121) time%day,time%mo,time%day_mo,time%yrc, &
															    rech_aa_p,soil_aa_p,gwsw_aa_p,swgw_aa_p,well_aa_p,canl_aa_p,drng_aa_p,resv_aa_p	
      endif
			
			!groundwater concentration
			!totals
			do k=1,nlay
			  do i=1,nrow
				  do j=1,ncol
					  if(ibound(j,i,k).ne.0) then
						  smrt_conc_tot_mo(j,i,k,1) = smrt_conc_tot_mo(j,i,k,1) + cnew(j,i,k,1)
							smrt_conc_tot_yr(j,i,k,1) = smrt_conc_tot_yr(j,i,k,1) + cnew(j,i,k,1)
							smrt_conc_tot_mo(j,i,k,2) = smrt_conc_tot_mo(j,i,k,1) + cnew(j,i,k,2)
							smrt_conc_tot_yr(j,i,k,2) = smrt_conc_tot_yr(j,i,k,1) + cnew(j,i,k,2)
						else
						  smrt_conc_tot_mo(j,i,k,1) = -5
							smrt_conc_tot_yr(j,i,k,1) = -5
							smrt_conc_tot_mo(j,i,k,2) = -5
							smrt_conc_tot_yr(j,i,k,2) = -5
						endif
					enddo
				enddo
			enddo

			!monthly output
			if (time%end_mo == 1) then
			  do k=1,nlay
					do i=1,nrow
						do j=1,ncol
							if(ibound(j,i,k).ne.0) then
								smrt_conc_tot_mo(j,i,k,1) = smrt_conc_tot_mo(j,i,k,1) / time%day_mo 
								smrt_conc_tot_mo(j,i,k,2) = smrt_conc_tot_mo(j,i,k,2) / time%day_mo 
							else
								smrt_conc_tot_mo(j,i,k,1) = -5
								smrt_conc_tot_mo(j,i,k,2) = -5
							endif
						enddo
					enddo
				enddo
			  write(rt_conc_mo,*), time%yrc,time%mo
			  do k=1,NLAY
          write(rt_conc_mo,*),'layer:',k
					!no3
					write(rt_conc_mo,*),'no3-n'
          do i=1,NROW
            write(rt_conc_mo,101) (smrt_conc_tot_mo(j,i,k,1),j=1,NCOL)
          enddo
					!p
					write(rt_conc_mo,*),'p'
          do i=1,NROW
            write(rt_conc_mo,101) (smrt_conc_tot_mo(j,i,k,2),j=1,NCOL)
          enddo
        enddo
				write(rt_conc_mo,*)
				smrt_conc_tot_mo = 0.
			endif
			
			!yearly output
			if (time%end_yr == 1) then
			  do k=1,nlay
					do i=1,nrow
						do j=1,ncol
							if(ibound(j,i,k).ne.0) then
								smrt_conc_tot_yr(j,i,k,1) = smrt_conc_tot_yr(j,i,k,1) / time%day
								smrt_conc_tot_yr(j,i,k,2) = smrt_conc_tot_yr(j,i,k,2) / time%day
							else
								smrt_conc_tot_yr(j,i,k,1) = -5
								smrt_conc_tot_yr(j,i,k,2) = -5
							endif
						enddo
					enddo
				enddo
			  write(rt_conc_yr,*), time%yrc
			  do k=1,NLAY
          write(rt_conc_yr,*),'layer:',k
					!no3
					write(rt_conc_yr,*),'no3-n'
          do i=1,NROW
            write(rt_conc_yr,101) (smrt_conc_tot_yr(j,i,k,1),j=1,NCOL)
          enddo
					!p
					write(rt_conc_yr,*),'p'
          do i=1,NROW
            write(rt_conc_yr,101) (smrt_conc_tot_yr(j,i,k,2),j=1,NCOL)
          enddo
        enddo
				write(rt_conc_yr,*)
				smrt_conc_tot_yr = 0.
			endif
			endif
			
			!groundwater nutrient loads ----------------------------------------------------------------------------------------
	    if(rt_active == 1) then
			!monthly output
			if (time%end_mo == 1) then
			  !recharge
			  if(IUNIT(8) > 0) then
				  smrt_rt_rech_tot_mo = smrt_rt_rech_tot_mo / time%day_mo
					write(rt_rech_mo,*), time%yrc,time%mo
					!no3
					write(rt_rech_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_rech_mo,101) (smrt_rt_rech_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_rech_mo,*), 'p'
          do i=1,NROW
            write(rt_rech_mo,101) (smrt_rt_rech_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_rech_mo,*)
					smrt_rt_rech_tot_mo = 0.
				endif
				!gw-->soil
				if(IUNIT(8) > 0) then
				  smrt_rt_soil_tot_mo = smrt_rt_soil_tot_mo / time%day_mo
					write(rt_soil_mo,*), time%yrc,time%mo
					!no3
					write(rt_soil_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_soil_mo,101) (smrt_rt_soil_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_soil_mo,*), 'p'
          do i=1,NROW
            write(rt_soil_mo,101) (smrt_rt_soil_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_soil_mo,*)
					smrt_rt_soil_tot_mo = 0.
				endif
				!gwsw
				if(IUNIT(4) > 0) then
				  smrt_rt_gwsw_tot_mo = smrt_rt_gwsw_tot_mo / time%day_mo
					write(rt_gwsw_mo,*), time%yrc,time%mo
					!no3
					write(rt_gwsw_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_gwsw_mo,101) (smrt_rt_gwsw_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_gwsw_mo,*), 'p'
          do i=1,NROW
            write(rt_gwsw_mo,101) (smrt_rt_gwsw_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_gwsw_mo,*)
					smrt_rt_gwsw_tot_mo = 0.   
				endif
				!pumping - well package
				if(IUNIT(2) > 0) then
				  smrt_rt_well_tot_mo = smrt_rt_well_tot_mo / time%day_mo
					write(rt_well_mo,*), time%yrc,time%mo
					!no3
					write(rt_well_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_well_mo,101) (smrt_rt_well_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_well_mo,*), 'p'
          do i=1,NROW
            write(rt_well_mo,101) (smrt_rt_well_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_well_mo,*)
					smrt_rt_well_tot_mo = 0.  
				endif
				!canal
				if(IUNIT(4) > 0) then
				  smrt_rt_canl_tot_mo = smrt_rt_canl_tot_mo / time%day_mo
					write(rt_canl_mo,*), time%yrc,time%mo
					!no3
					write(rt_canl_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_canl_mo,101) (smrt_rt_canl_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_canl_mo,*), 'p'
          do i=1,NROW
            write(rt_canl_mo,101) (smrt_rt_canl_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_canl_mo,*)
					smrt_rt_canl_tot_mo = 0.  
				endif
				!drainage
				if(IUNIT(3) > 0) then
				  smrt_rt_drng_tot_mo = smrt_rt_drng_tot_mo / time%day_mo
					write(rt_drng_mo,*), time%yrc,time%mo
					!no3
					write(rt_drng_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_drng_mo,101) (smrt_rt_drng_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_drng_mo,*), 'p'
          do i=1,NROW
            write(rt_drng_mo,101) (smrt_rt_drng_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_drng_mo,*)
					smrt_rt_drng_tot_mo = 0.  
				endif
				!reservoir
			  if(IUNIT(17) > 0) then
				  smrt_rt_resv_tot_mo = smrt_rt_resv_tot_mo / time%day_mo
					write(rt_resv_mo,*), time%yrc,time%mo
					!no3
					write(rt_resv_mo,*) 'no3-n'
          do i=1,NROW
            write(rt_resv_mo,101) (smrt_rt_resv_tot_mo(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_resv_mo,*), 'p'
          do i=1,NROW
            write(rt_resv_mo,101) (smrt_rt_resv_tot_mo(j,i,2),j=1,NCOL)
          enddo
					write(rt_resv_mo,*)
					smrt_rt_resv_tot_mo = 0.  
				endif 
			endif
			
			!yearly output
			if(time%end_yr == 1) then
			  !recharge
			  if(IUNIT(8) > 0) then
				  smrt_rt_rech_tot_yr = smrt_rt_rech_tot_yr / time%day
					write(rt_rech_yr,*), time%yrc
					!no3
					write(rt_rech_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_rech_yr,101) (smrt_rt_rech_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_rech_yr,*), 'p'
          do i=1,NROW
            write(rt_rech_yr,101) (smrt_rt_rech_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_rech_yr,*)
					smrt_rt_rech_tot_yr = 0.
				endif
				!gw-->soil
				if(IUNIT(8) > 0) then
				  smrt_rt_soil_tot_yr = smrt_rt_soil_tot_yr / time%day
					write(rt_soil_yr,*), time%yrc
					!no3
					write(rt_soil_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_soil_yr,101) (smrt_rt_soil_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_soil_yr,*), 'p'
          do i=1,NROW
            write(rt_soil_yr,101) (smrt_rt_soil_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_soil_yr,*)
					smrt_rt_soil_tot_yr = 0.
				endif
				!gwsw
				if(IUNIT(4) > 0) then
				  smrt_rt_gwsw_tot_yr = smrt_rt_gwsw_tot_yr / time%day
					write(rt_gwsw_yr,*), time%yrc
					!no3
					write(rt_gwsw_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_gwsw_yr,101) (smrt_rt_gwsw_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_gwsw_yr,*), 'p'
          do i=1,NROW
            write(rt_gwsw_yr,101) (smrt_rt_gwsw_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_gwsw_yr,*)
					smrt_rt_gwsw_tot_yr = 0.   
				endif
				!pumping - well package
				if(IUNIT(2) > 0) then
				  smrt_rt_well_tot_yr = smrt_rt_well_tot_yr / time%day
					write(rt_well_yr,*), time%yrc
					!no3
					write(rt_well_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_well_yr,101) (smrt_rt_well_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_well_yr,*), 'p'
          do i=1,NROW
            write(rt_well_yr,101) (smrt_rt_well_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_well_yr,*)
					smrt_rt_well_tot_yr = 0.  
				endif
				!canal
				if(IUNIT(4) > 0) then
				  smrt_rt_canl_tot_yr = smrt_rt_canl_tot_yr / time%day
					write(rt_canl_yr,*), time%yrc
					!no3
					write(rt_canl_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_canl_yr,101) (smrt_rt_canl_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_canl_yr,*), 'p'
          do i=1,NROW
            write(rt_canl_yr,101) (smrt_rt_canl_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_canl_yr,*)
					smrt_rt_canl_tot_yr = 0.  
				endif
				!drainage
				if(IUNIT(3) > 0) then
				  smrt_rt_drng_tot_yr = smrt_rt_drng_tot_yr / time%day
					write(rt_drng_yr,*), time%yrc
					!no3
					write(rt_drng_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_drng_yr,101) (smrt_rt_drng_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_drng_yr,*), 'p'
          do i=1,NROW
            write(rt_drng_yr,101) (smrt_rt_drng_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_drng_yr,*)
					smrt_rt_drng_tot_yr = 0.  
				endif
				!reservoir
			  if(IUNIT(17) > 0) then
				  smrt_rt_resv_tot_yr = smrt_rt_resv_tot_yr / time%day
					write(rt_resv_yr,*), time%yrc
					!no3
					write(rt_resv_yr,*) 'no3-n'
          do i=1,NROW
            write(rt_resv_yr,101) (smrt_rt_resv_tot_yr(j,i,1),j=1,NCOL)
          enddo
					!p
					write(rt_resv_yr,*), 'p'
          do i=1,NROW
            write(rt_resv_yr,101) (smrt_rt_resv_tot_yr(j,i,2),j=1,NCOL)
          enddo
					write(rt_resv_yr,*)
					smrt_rt_resv_tot_yr = 0.  
				endif  
			endif
			endif
			
			
101   format(10000(e13.5))			  
120   format(i8,i8,i8,i8,e16.7,e16.7,1000(e13.4))
121   format(i8,i8,i8,i8,1000(e13.4))
			
      return
			end
