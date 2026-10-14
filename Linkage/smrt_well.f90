      subroutine smrt_well(ob_id_src,ob_id_dmd,demand_vol,extracted,dmd_unmet,demand_type) !rtb gwflow

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine determines the volume of groundwater that is extracted
!!    from modflow grid cells
!!    (pumping volume are used in mf_run, in modflow storage equations)
      
      use gwflow_module
      use organic_mineral_mass_module, only : soil1
      use hru_module, only : hru,irrn,irrp
      use hydrograph_module, only : wet,ob
      use soil_module, only : soil
      use GWFWELMODULE, only:WELL !modflow
			use GLOBAL, only:HNEW,BOTM,ncol,nrow,DELC,DELR,IUNIT !modflow
			use GWFUPWMODULE, only:SC2UPW,LAYTYPUPW,SC1 !modflow
      use GWFLPFMODULE, only:SC2  !modflow
			use smrt_parm
			use mf_rt_link, only:rt_active
			use rt_global, only:cnew
			
      
      
      implicit none

			integer, intent (in) :: ob_id_src	    !       |object designated as the water source
			integer, intent (in) :: ob_id_dmd	    !       |object that has a water demand
      real, intent (in) :: demand_vol       !m3     |volume of water demand
      real, intent (inout) :: extracted     !m3     |volume of groundwater extracted from aquifer for irrigation
      real, intent (inout) :: dmd_unmet     !m3     |volume of irrigation demand not met by aquifer
			character*4, intent (in) :: demand_type
      integer :: i                          !       |counter
      integer :: s                          !       |solute counter
			integer :: cell_row,cell_col          !       |row and column of modflow cell
      integer :: cell_id                    !       |modflow cell id
      integer :: wetland                    !       |wetland flag
			integer :: hru_id
      real :: thickness                     !m      |groundwater saturated thickness
			real :: sy                            !       |aquifer specific yield
			real :: cell_area											!m2     |surface area of modflow cell
			real :: cell_demand                   !m3     |volume of irrigation water demand per gwflow cell
      real :: gwvol_avail                   !m3     |groundwater storage in the gwflow cell
      real :: gwvol_removed									!m3     |groundwater removed from the gwflow cell for irrigation demand
      real :: gwvol_unmet										!m3     |groundwater not available to meet irrigation demand
      real :: gw_mass                       !kg     |mass of solute in groundwater of the gwflow cell
      real :: irr_mass_no3                  !kg     |mass of nitrate removed from aquifer for irrigation
			real :: irr_mass_p                    !kg     |mass of phosphorus removed from aquifer for irrigation
      real :: mass_diff                     !kg     |difference between irrigation mass and actual groundwater mass
      real :: hru_area_m2                   !m2     |HRU area
			real :: stor_coeff                    !-      |groundwater storage coefficient (either specific yield or storativity)
			real :: remainder
      
      
      !zero out the total met and unmet demand
      extracted = 0.
      dmd_unmet = 0.

			!first option: hru irrigation from all cells geographically connected to the hru
			!hru demand --> irrigation ------------------------------------------------------------------------------------------------
			if(demand_type == "hru" .and. ob_id_src == 0) then
			
      !only proceed if the HRU is connected (intersected with) gwflow grid cells
      hru_id = ob_id_dmd
			if(mf_hru_num_cells(hru_id) > 0) then
          
        !groundwater volume to remove from each cell connected to the HRU
        cell_demand = demand_vol / mf_hru_num_cells(hru_id) !groundwater to remove from each cell connected to the HRU
        
        !loop through the cells that are connected to the HRU
        do i=1,mf_hru_num_cells(hru_id)
          
          !row and column of the current cell
				  cell_row = mf_hru_cell_row(hru_id,i)
			    cell_col = mf_hru_cell_col(hru_id,i)
				  
					!check available groundwater in the cell
					!determine saturated thickness (m) of the current layer
          thickness = 0.
					if(HNEW(cell_col,cell_row,mf_well_lay) < BOTM(cell_col,cell_row,0) .and. HNEW(cell_col,cell_row,mf_well_lay) > BOTM(cell_col,cell_row,mf_well_lay)) then 
            thickness = HNEW(cell_col,cell_row,mf_well_lay) - BOTM(cell_col,cell_row,mf_well_lay)
          else
            thickness = BOTM(cell_col,cell_row,0) - BOTM(cell_col,cell_row,1)
          endif
					cell_area = DELR(cell_col) * DELC(cell_row)
					stor_coeff = 0.
					if(LAYTYPUPW(mf_well_lay) > 0) then !not confined - use specific yield
					  if(IUNIT(23).gt.0) then
              stor_coeff = SC2(cell_col,cell_row,mf_well_lay) / cell_area
            else
              stor_coeff = SC2UPW(cell_col,cell_row,mf_well_lay) / cell_area
            endif
					elseif(LAYTYPUPW(mf_well_lay) == 0) then !confined - use specific storage
					  stor_coeff = (SC1(cell_col,cell_row,mf_well_lay) / cell_area) * thickness
					endif
					gwvol_avail = thickness * cell_area * stor_coeff !m3
					
					!remove irrigation demand from storage; if storage is less than demand, remove all groundwater
          gwvol_removed = 0.
          gwvol_unmet = 0.
          if(gwvol_avail < cell_demand) then
            gwvol_removed = gwvol_avail
            gwvol_unmet = cell_demand - gwvol_avail !amount that is not available for irrigation
          else
            gwvol_removed = cell_demand
          endif
          extracted = extracted + gwvol_removed
          dmd_unmet = dmd_unmet + gwvol_unmet
					
					!add pumping rate to well package array
					mf_num_well = mf_num_well + 1
					mf_well_id = mf_well_id + 1
					WELL(1,mf_well_id) = mf_well_lay
          WELL(2,mf_well_id) = cell_row
          WELL(3,mf_well_id) = cell_col
          WELL(4,mf_well_id) = gwvol_removed * (-1) !make it extraction
					
          !add solute mass in irrigation water to HRU soil profile
          if (rt_active == 1) then
					  !first: determine the nutrient mass in irrigation water
					  !no3
					  gw_mass = cnew(cell_col,cell_row,mf_well_lay,1) * gwvol_avail !g
					  irr_mass_no3 = cnew(cell_col,cell_row,mf_well_lay,1) * gwvol_removed !g
						mass_diff = 0.
            if(irr_mass_no3 > gw_mass) then
              mass_diff = irr_mass_no3 - gw_mass
            endif 
						irr_mass_no3 = irr_mass_no3 - mass_diff
						if(irr_mass_no3 < 0) irr_mass_no3 = 0.
						irr_mass_no3 = irr_mass_no3 / 1000. !kg
						rt_well_load(cell_col,cell_row,mf_well_lay) = irr_mass_no3
					  !p
						gw_mass = cnew(cell_col,cell_row,mf_well_lay,2) * gwvol_avail !g
					  irr_mass_p = cnew(cell_col,cell_row,mf_well_lay,2) * gwvol_removed !g
						mass_diff = 0.
            if(irr_mass_p > gw_mass) then
              mass_diff = irr_mass_p - gw_mass
            endif 
						irr_mass_p = irr_mass_p - mass_diff
						if(irr_mass_p < 0) irr_mass_p = 0. 
						irr_mass_p = irr_mass_p / 1000. !kg
						rt_well_load(cell_col,cell_row,2) = irr_mass_p
						!second: add the nutrient mass to the soil profile (hru)
            wetland = hru(hru_id)%dbs%surf_stor !check if HRU is a wetland
            if(wetland > 0) then !add to wetland
              !nutrients
              wet(hru_id)%no3 = wet(hru_id)%no3 + irr_mass_no3 !kg
              wet(hru_id)%solp = wet(hru_id)%solp + irr_mass_p !kg
            else !add to soil profile
              !nutrients
              soil1(hru_id)%mn(1)%no3 = soil1(hru_id)%mn(1)%no3 + (irr_mass_no3/hru(hru_id)%area_ha) !kg/ha - add to soil layer
              soil1(hru_id)%mp(1)%lab = soil1(hru_id)%mp(1)%lab + (irr_mass_p/hru(hru_id)%area_ha) !kg/ha - add to soil layer
            endif
          endif
          
        enddo !go to next modflow cell
        
			endif
      			
			
			
			!irrigation pumping from a single cell ------------------------------------------------------------------------------------
			elseif(demand_type == "hru" .and. ob_id_src > 0) then
        
			  hru_id = ob_id_dmd
        cell_demand = demand_vol !m3 groundwater volume to remove from the cell  
				
				!determine the row and column of the cell
			  remainder = mod(ob_id_src,ncol)
        if(remainder == 0) then
          cell_row = cell_id / ncol !row
        else
				  cell_row = int(cell_id/ncol) + 1 !row
        endif
			  cell_col = cell_id - ((cell_row-1)*ncol) !column
				
				!check available groundwater in the cell
				!determine saturated thickness (m) of the current layer
        thickness = 0.
				if(HNEW(cell_col,cell_row,mf_well_lay) < BOTM(cell_col,cell_row,0) .and. HNEW(cell_col,cell_row,mf_well_lay) > BOTM(cell_col,cell_row,mf_well_lay)) then 
          thickness = HNEW(cell_col,cell_row,mf_well_lay) - BOTM(cell_col,cell_row,mf_well_lay)
        else
          thickness = BOTM(cell_col,cell_row,0) - BOTM(cell_col,cell_row,1)
        endif
				cell_area = DELR(cell_col) * DELC(cell_row)
				stor_coeff = 0.
				if(LAYTYPUPW(mf_well_lay) > 0) then !not confined - use specific yield
					if(IUNIT(23).gt.0) then
            stor_coeff = SC2(cell_col,cell_row,mf_well_lay) / cell_area
          else
            stor_coeff = SC2UPW(cell_col,cell_row,mf_well_lay) / cell_area
          endif
				elseif(LAYTYPUPW(mf_well_lay) == 0) then !confined - use specific storage
					stor_coeff = (SC1(cell_col,cell_row,mf_well_lay) / cell_area) * thickness
				endif
				gwvol_avail = thickness * cell_area * stor_coeff !m3
				
        !remove irrigation demand from storage; if storage is less than demand, remove all groundwater
        gwvol_removed = 0.
        gwvol_unmet = 0.
        if(gwvol_avail < cell_demand) then
          gwvol_removed = gwvol_avail
          gwvol_unmet = cell_demand - gwvol_avail !amount that is not available for irrigation
        else
          gwvol_removed = cell_demand
        endif
        extracted = extracted + gwvol_removed
        dmd_unmet = dmd_unmet + gwvol_unmet

				!add pumping rate to well package array
				mf_num_well = mf_num_well + 1
				mf_well_id = mf_well_id + 1
				WELL(1,mf_well_id) = mf_well_lay
        WELL(2,mf_well_id) = cell_row
        WELL(3,mf_well_id) = cell_col
        WELL(4,mf_well_id) = gwvol_removed * (-1) !make it extraction
				  
        !add solute mass in irrigation water to HRU soil profile
        if (rt_active == 1) then
					!first: determine the nutrient mass in irrigation water
					!no3
					gw_mass = cnew(cell_col,cell_row,mf_well_lay,1) * gwvol_avail !g
					irr_mass_no3 = cnew(cell_col,cell_row,mf_well_lay,1) * gwvol_removed !g
					mass_diff = 0.
          if(irr_mass_no3 > gw_mass) then
            mass_diff = irr_mass_no3 - gw_mass
          endif 
					irr_mass_no3 = irr_mass_no3 - mass_diff
					if(irr_mass_no3 < 0) irr_mass_no3 = 0.
					irr_mass_no3 = irr_mass_no3 / 1000. !kg
					rt_well_load(cell_col,cell_row,mf_well_lay) = irr_mass_no3
					!p
					gw_mass = cnew(cell_col,cell_row,mf_well_lay,2) * gwvol_avail !g
					irr_mass_p = cnew(cell_col,cell_row,mf_well_lay,2) * gwvol_removed !g
					mass_diff = 0.
          if(irr_mass_p > gw_mass) then
            mass_diff = irr_mass_p - gw_mass
          endif 
					irr_mass_p = irr_mass_p - mass_diff
					if(irr_mass_p < 0) irr_mass_p = 0. 
					irr_mass_p = irr_mass_p / 1000. !kg
					rt_well_load(cell_col,cell_row,2) = irr_mass_p
					!second: add the nutrient mass to the soil profile (hru)
          wetland = hru(hru_id)%dbs%surf_stor !check if HRU is a wetland
          if(wetland > 0) then !add to wetland
            !nutrients
            wet(hru_id)%no3 = wet(hru_id)%no3 + irr_mass_no3 !kg
            wet(hru_id)%solp = wet(hru_id)%solp + irr_mass_p !kg
          else !add to soil profile
            !nutrients
            soil1(hru_id)%mn(1)%no3 = soil1(hru_id)%mn(1)%no3 + (irr_mass_no3/hru(hru_id)%area_ha) !kg/ha - add to soil layer
            soil1(hru_id)%mp(1)%lab = soil1(hru_id)%mp(1)%lab + (irr_mass_p/hru(hru_id)%area_ha) !kg/ha - add to soil layer
          endif
        endif
				
			  
				
			!municipal demand from all cells connected geographically to the hru ------------------------------------------------------
			elseif(demand_type == "muni" .and. ob_id_src == 0) then

			hru_id = ob_id_dmd
			if(mf_hru_num_cells(hru_id) > 0) then
          
        !groundwater volume to remove from each cell connected to the HRU
        cell_demand = demand_vol / mf_hru_num_cells(hru_id) !groundwater to remove from each cell connected to the HRU
        
        !loop through the cells that are connected to the HRU
        do i=1,mf_hru_num_cells(hru_id)
          
          !row and column of the current cell
				  cell_row = mf_hru_cell_row(hru_id,i)
			    cell_col = mf_hru_cell_col(hru_id,i)
				  
					!check available groundwater in the cell
					!determine saturated thickness (m) of the current layer
          thickness = 0.
					if(HNEW(cell_col,cell_row,mf_well_lay) < BOTM(cell_col,cell_row,0) .and. HNEW(cell_col,cell_row,mf_well_lay) > BOTM(cell_col,cell_row,mf_well_lay)) then 
            thickness = HNEW(cell_col,cell_row,mf_well_lay) - BOTM(cell_col,cell_row,mf_well_lay)
          else
            thickness = BOTM(cell_col,cell_row,0) - BOTM(cell_col,cell_row,1)
          endif
					cell_area = DELR(cell_col) * DELC(cell_row)
					stor_coeff = 0.
					if(LAYTYPUPW(mf_well_lay) > 0) then !not confined - use specific yield
					  if(IUNIT(23).gt.0) then
              stor_coeff = SC2(cell_col,cell_row,mf_well_lay) / cell_area
            else
              stor_coeff = SC2UPW(cell_col,cell_row,mf_well_lay) / cell_area
            endif
					elseif(LAYTYPUPW(mf_well_lay) == 0) then !confined - use specific storage
					  stor_coeff = (SC1(cell_col,cell_row,mf_well_lay) / cell_area) * thickness
					endif
					gwvol_avail = thickness * cell_area * stor_coeff !m3
					
					!remove municipal demand from storage; if storage is less than demand, remove all groundwater
          gwvol_removed = 0.
          gwvol_unmet = 0.
          if(gwvol_avail < cell_demand) then
            gwvol_removed = gwvol_avail
            gwvol_unmet = cell_demand - gwvol_avail !amount that is not available for irrigation
          else
            gwvol_removed = cell_demand
          endif
          extracted = extracted + gwvol_removed
          dmd_unmet = dmd_unmet + gwvol_unmet
					
					!add pumping rate to well package array
					mf_num_well = mf_num_well + 1
					mf_well_id = mf_well_id + 1
					WELL(1,mf_well_id) = mf_well_lay
          WELL(2,mf_well_id) = cell_row
          WELL(3,mf_well_id) = cell_col
          WELL(4,mf_well_id) = gwvol_removed * (-1) !make it extraction
					
        enddo !go to next modflow cell
        
			endif
				
				
				
			!municipal demand from a single cell --------------------------------------------------------------------------------------
			elseif(demand_type == "muni" .and. ob_id_src > 0) then

			  !m3 groundwater volume to remove from the cell  
        cell_demand = demand_vol 
        
				!determine the row and column of the cell
			  remainder = mod(ob_id_src,ncol)
        if(remainder == 0) then
          cell_row = cell_id / ncol !row
        else
				  cell_row = int(cell_id/ncol) + 1 !row
        endif
			  cell_col = cell_id - ((cell_row-1)*ncol) !column
				
        !check available groundwater in the cell
				!determine saturated thickness (m) of the current layer
        thickness = 0.
				if(HNEW(cell_col,cell_row,mf_well_lay) < BOTM(cell_col,cell_row,0) .and. HNEW(cell_col,cell_row,mf_well_lay) > BOTM(cell_col,cell_row,mf_well_lay)) then 
          thickness = HNEW(cell_col,cell_row,mf_well_lay) - BOTM(cell_col,cell_row,mf_well_lay)
        else
          thickness = BOTM(cell_col,cell_row,0) - BOTM(cell_col,cell_row,1)
        endif
				cell_area = DELR(cell_col) * DELC(cell_row)
				stor_coeff = 0.
				if(LAYTYPUPW(mf_well_lay) > 0) then !not confined - use specific yield
					if(IUNIT(23).gt.0) then
            stor_coeff = SC2(cell_col,cell_row,mf_well_lay) / cell_area
          else
            stor_coeff = SC2UPW(cell_col,cell_row,mf_well_lay) / cell_area
          endif
				elseif(LAYTYPUPW(mf_well_lay) == 0) then !confined - use specific storage
					stor_coeff = (SC1(cell_col,cell_row,mf_well_lay) / cell_area) * thickness
				endif
				gwvol_avail = thickness * cell_area * stor_coeff !m3
					
				!remove municipal demand from storage; if storage is less than demand, remove all groundwater
        gwvol_removed = 0.
        gwvol_unmet = 0.
        if(gwvol_avail < cell_demand) then
          gwvol_removed = gwvol_avail
          gwvol_unmet = cell_demand - gwvol_avail !amount that is not available for irrigation
        else
          gwvol_removed = cell_demand
        endif
        extracted = extracted + gwvol_removed
        dmd_unmet = dmd_unmet + gwvol_unmet
					
				!add pumping rate to well package array
				mf_num_well = mf_num_well + 1
				mf_well_id = mf_well_id + 1
				WELL(1,mf_well_id) = mf_well_lay
        WELL(2,mf_well_id) = cell_row
        WELL(3,mf_well_id) = cell_col
        WELL(4,mf_well_id) = gwvol_removed * (-1) !make it extraction
				
			endif !check: hru or muni
			

			return
      end subroutine smrt_well
			