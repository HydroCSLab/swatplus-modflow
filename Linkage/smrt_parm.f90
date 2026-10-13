      module smrt_parm

	    implicit none
	
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This module contains variables used for linking SWAT variables to MODFLOW and RT3D

      !Flag for whether modflow is active -----------------------------------------------------------------------------
      integer :: mf_active
      integer :: leapyr
      
	    !package flags
	    integer :: mf_obs_flag = 0
	    integer :: mf_evt_flag = 0
			integer :: mf_cnl_flag = 0
	    integer :: mf_drn_flag = 0
			integer :: mf_res_flag = 0
	
			!groundwater balance files
			integer :: out_gwbal_dy = 29996
			integer :: out_gwbal_mo = 29997
			integer :: out_gwbal_yr = 29998
			integer :: out_gwbal_aa = 29999
			
			!groundwater nutrient balance files
			integer :: out_rtbal_dy_n = 29980
			integer :: out_rtbal_mo_n = 29981
			integer :: out_rtbal_yr_n = 29982
			integer :: out_rtbal_aa_n = 29983
			integer :: out_rtbal_dy_p = 29984
			integer :: out_rtbal_mo_p = 29985
			integer :: out_rtbal_yr_p = 29986
			integer :: out_rtbal_aa_p = 29987
			
	    !HRU-cell linkage
	    integer, dimension (:), allocatable :: mf_hrus_connected			 !        |flag for each HRU: connected (1) or not connected (0) to cells
	    integer, dimension (:), allocatable :: mf_hru_num_cells				 !        |number of cells connected to each HRU
	    integer, dimension (:,:), allocatable :: mf_hru_cell_row		   !        |row of each cell connected to the HRU
	    integer, dimension (:,:), allocatable :: mf_hru_cell_col		   !        |column of each cell connected to the HRU
	    real, dimension (:,:), allocatable :: mf_hru_cell_fract   		 !        |fraction of HRU intersected with the cell
	    real, dimension (:,:), allocatable :: mf_cell_fract				     !        |fraction of cell intersected with the HRU
	     
			!channel-cell linkage
			integer, dimension (:), allocatable :: mf_chan_channel         !				|ID of channel connected to modflow cell
			real, dimension (:), allocatable :: mf_rivr_lngt               !m				|length of channel in the modflow cell
			real, dimension (:), allocatable :: mf_rivr_hydc               !m/day   |hydraulic conductivity of channel bed
			real, dimension (:), allocatable :: mf_rivr_thck               !m       |thickness of channel bed
			real, dimension (:), allocatable :: mf_chan_rates              !m3      |daily exchange rate between modflow cell and channel
			
	    !arrays for transferring water and solute mass from MODFLOW cells to SWAT HRUs (soil profile)
      integer :: mf_soil_flag = 0
			real, dimension (:,:), allocatable :: gw_transfer
      real, dimension (:,:), allocatable :: RECH_down
      real, dimension (:,:,:), allocatable :: mass_cell
	    real, dimension (:,:,:), allocatable :: mf_hru_soil
			
      !groundwater ET (EVT package)
      real, dimension (:), allocatable :: etremain
      
			!groundwater pumping (WELL package)
			integer :: mf_num_well = 0
			integer :: mf_well_id = 0
			integer :: NWELLS_og																					 !        |number of well cells besides smrt linkage wells
			integer :: mf_well_lay = 0																		 !        |layer from which agricultural pumping occurs
			
			!groundwater-canal exchange (RIV package)
			integer :: mf_num_canl = 0
			integer, dimension (:), allocatable :: mf_canl_channel
			integer, dimension (:), allocatable :: mf_canl_div
			integer, dimension (:), allocatable :: mf_canl_beg
			integer, dimension (:), allocatable :: mf_canl_end
			real, dimension (:), allocatable :: mf_canl_cond
			real, dimension (:), allocatable :: mf_canl_rates              !m3      |daily exchange rate between modflow cell and canal
			
			!groundwater drainage (DRN package)
			integer :: mf_num_drng
			integer, dimension (:), allocatable :: mf_drng_channel
			real, dimension (:), allocatable :: mf_drng_rates
			
			!groundwater-reservoir exchange (RES package)
			integer :: mf_num_resv
			integer, dimension (:), allocatable :: mf_resv_reservoir
			integer, dimension (:), allocatable :: mf_resv_row,mf_resv_col
			real, dimension (:), allocatable :: mf_resv_rates
			
      !flags for writing output data (SWAT, MODFLOW)
      integer :: out_SWAT_recharge,out_MF_recharge,out_SWAT_channel,out_MF_channel,out_MODFLOW_gwsw,out_SWAT_gwsw

      !array for storing MODFLOW observation cells
      integer :: num_mf_obs
      integer, dimension (:,:), allocatable :: mf_obs

      !variables for calculating groundwater balance
			!monthly
      real :: gw_volume_mo = 0.
			real :: gw_change_mo = 0.
			real :: rech_mo = 0.
			real :: soil_mo = 0.
			real :: gwet_mo = 0.
			real :: gwsw_mo = 0.
			real :: swgw_mo = 0.
			real :: bndr_mo = 0.
			real :: well_mo = 0.
			real :: mnw2_mo = 0.
			real :: canl_mo = 0.
			real :: drng_mo = 0.
			real :: resv_mo = 0.
			!yearly
      real :: gw_volume_yr = 0.
			real :: gw_change_yr = 0.
			real :: rech_yr = 0.
			real :: soil_yr = 0.
			real :: gwet_yr = 0.
			real :: gwsw_yr = 0.
			real :: swgw_yr = 0.
			real :: bndr_yr = 0.
			real :: well_yr = 0.
			real :: mnw2_yr = 0.
			real :: canl_yr = 0.
			real :: drng_yr = 0.
			real :: resv_yr = 0.
			!average annual
      real :: gw_volume_aa = 0.
			real :: gw_change_aa = 0.
			real :: rech_aa = 0.
			real :: soil_aa = 0.
			real :: gwet_aa = 0.
			real :: gwsw_aa = 0.
			real :: swgw_aa = 0.
			real :: bndr_aa = 0.
			real :: well_aa = 0.
			real :: mnw2_aa = 0.
			real :: canl_aa = 0.
			real :: drng_aa = 0.
			real :: resv_aa = 0.
			
			!variables for calculating groundwater nutrient balance
			!no3-n
			!monthly
			real :: rech_mo_no3 = 0.
			real :: soil_mo_no3 = 0.
			real :: gwsw_mo_no3 = 0.
			real :: swgw_mo_no3 = 0.
			real :: well_mo_no3 = 0.
			real :: canl_mo_no3 = 0.
			real :: drng_mo_no3 = 0.
			real :: resv_mo_no3 = 0.
			real :: rctn_mo_no3 = 0.
			!yearly
			real :: rech_yr_no3 = 0.
			real :: soil_yr_no3 = 0.
			real :: gwsw_yr_no3 = 0.
			real :: swgw_yr_no3 = 0.
			real :: well_yr_no3 = 0.
			real :: canl_yr_no3 = 0.
			real :: drng_yr_no3 = 0.
			real :: resv_yr_no3 = 0.
			real :: rctn_yr_no3 = 0.
			!average annual
			real :: rech_aa_no3 = 0.
			real :: soil_aa_no3 = 0.
			real :: gwsw_aa_no3 = 0.
			real :: swgw_aa_no3 = 0.
			real :: well_aa_no3 = 0.
			real :: canl_aa_no3 = 0.
			real :: drng_aa_no3 = 0.
			real :: resv_aa_no3 = 0.
			real :: rctn_aa_no3 = 0.
			!p
			!monthly
			real :: rech_mo_p = 0.
			real :: soil_mo_p = 0.
			real :: gwsw_mo_p = 0.
			real :: swgw_mo_p = 0.
			real :: well_mo_p = 0.
			real :: canl_mo_p = 0.
			real :: drng_mo_p = 0.
			real :: resv_mo_p = 0.
			!yearly
			real :: rech_yr_p = 0.
			real :: soil_yr_p = 0.
			real :: gwsw_yr_p = 0.
			real :: swgw_yr_p = 0.
			real :: well_yr_p = 0.
			real :: canl_yr_p = 0.
			real :: drng_yr_p = 0.
			real :: resv_yr_p = 0.
			!average annual
			real :: rech_aa_p = 0.
			real :: soil_aa_p = 0.
			real :: gwsw_aa_p = 0.
			real :: swgw_aa_p = 0.
			real :: well_aa_p = 0.
			real :: canl_aa_p = 0.
			real :: drng_aa_p = 0.
			real :: resv_aa_p = 0.
			
      !routing soil percolation water to the water table
      real, dimension (:), allocatable :: gw_delaye,gw_rech
      real, dimension (:,:), allocatable :: mf_percsol               !kg/ha    |solute mass leaving the soil profile
      real, dimension (:,:), allocatable :: mf_rechsol               !kg/ha    |solute mass in daily recharge (reaching water table)
			
			!other flags
			logical :: mf_nat_model

			!arrays for storing monthly and annual head and flux values
			real, dimension (:,:,:), allocatable :: smrt_head_tot_mo, &
			                                        smrt_head_tot_yr
			real, dimension (:,:), allocatable :: smrt_rech_tot_mo, &
																						smrt_rech_tot_yr
			real, dimension (:,:), allocatable :: smrt_soil_tot_mo, &
																						smrt_soil_tot_yr
			real, dimension (:,:), allocatable :: smrt_gwet_tot_mo, &
																						smrt_gwet_tot_yr
			real, dimension (:,:), allocatable :: smrt_gwsw_tot_mo, &
																						smrt_gwsw_tot_yr
			real, dimension (:,:), allocatable :: smrt_well_tot_mo, &
																						smrt_well_tot_yr
			real, dimension (:,:), allocatable :: smrt_mnw2_tot_mo, &
																						smrt_mnw2_tot_yr																																											
      real, dimension (:,:), allocatable :: smrt_canl_tot_mo, &
																						smrt_canl_tot_yr
      real, dimension (:,:), allocatable :: smrt_drng_tot_mo, &
																						smrt_drng_tot_yr
			real, dimension (:,:), allocatable :: smrt_resv_tot_mo, &
																						smrt_resv_tot_yr																	
			
			!arrays for storing monthly and annual conc and nutrient flux values
			real, dimension (:,:,:,:), allocatable :: smrt_conc_tot_mo, &
			                                          smrt_conc_tot_yr											
			real, dimension (:,:,:), allocatable :: smrt_rt_rech_tot_mo, &
															  							smrt_rt_rech_tot_yr
			real, dimension (:,:,:), allocatable :: smrt_rt_soil_tot_mo, &
																		  				smrt_rt_soil_tot_yr
			real, dimension (:,:,:), allocatable :: smrt_rt_gwsw_tot_mo, &
																	  					smrt_rt_gwsw_tot_yr
			real, dimension (:,:,:), allocatable :: smrt_rt_well_tot_mo, &
																  						smrt_rt_well_tot_yr
		  real, dimension (:,:,:), allocatable :: smrt_rt_canl_tot_mo, &
															  							smrt_rt_canl_tot_yr
			real, dimension (:,:,:), allocatable :: smrt_rt_drng_tot_mo, &
															  							smrt_rt_drng_tot_yr
			real, dimension (:,:,:), allocatable :: smrt_rt_resv_tot_mo, &
														  								smrt_rt_resv_tot_yr	
						
			!nutrient load arrays
			real, dimension (:,:,:), allocatable :: rt_rech_load, &
			                                        rt_soil_load, &
															  							rt_gwsw_load, &
																  						rt_well_load, &
																							rt_canl_load, &
																	  					rt_drng_load, &
																						  rt_resv_load, &
			                                        rt_rctn_load
      																						
			!output files
			!smrt log file
			integer :: smrt_log = 30400		
			integer :: mf_head_ob_in = 30401
			integer :: mf_head_ob_out = 30402
			!groundwater head and solute concentration
			integer :: mf_head_mo = 30500
			integer :: mf_head_yr = 30501
			integer :: rt_conc_mo = 30502
			integer :: rt_conc_yr = 30503
			!groundwater flux	- monthly
			integer :: mf_rech_mo = 30600
			integer :: mf_soil_mo = 30601
			integer :: mf_gwet_mo = 30602
			integer :: mf_gwsw_mo = 30603
			integer :: mf_well_mo = 30604
			integer :: mf_mnw2_mo = 30605
			integer :: mf_canl_mo = 30606
			integer :: mf_drng_mo = 30607
			integer :: mf_resv_mo = 30608
			integer :: rt_rech_mo = 30650
			integer :: rt_soil_mo = 30651
			integer :: rt_gwsw_mo = 30652
			integer :: rt_well_mo = 30653
			integer :: rt_canl_mo = 30654
			integer :: rt_drng_mo = 30655
			integer :: rt_resv_mo = 30656
			!groundwater flux	- yearly
			integer :: mf_rech_yr = 30700
			integer :: mf_soil_yr = 30701
			integer :: mf_gwet_yr = 30702
			integer :: mf_gwsw_yr = 30703
			integer :: mf_well_yr = 30704
			integer :: mf_mnw2_yr = 30705
			integer :: mf_canl_yr = 30706
			integer :: mf_drng_yr = 30707
			integer :: mf_resv_yr = 30708
			integer :: rt_rech_yr = 30750
			integer :: rt_soil_yr = 30751
			integer :: rt_gwsw_yr = 30752
			integer :: rt_well_yr = 30753
			integer :: rt_canl_yr = 30754
			integer :: rt_drng_yr = 30755
			integer :: rt_resv_yr = 30756
			
			
      end module smrt_parm
