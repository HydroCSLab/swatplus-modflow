      subroutine smrt_initialize

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine prepares information for a linked SWAT, MODFLOW, RT3D
!!    simulation.
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mf_active    |none          |index whether or not to use MODFLOW
!!    rt_active    |none          |index whether or not to use RT3D
!!    mf_interval  |days          |number of days between MODFLOW calls
!!    mf_obs_flag  |none          |index whether or not MODFLOW observation cells are
!!                                 read in (for detailed output)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use hru_module
      use hydrograph_module
      use basin_module
      use GLOBAL, only: mf_interval,mf_ran,NROW,NCOL,NLAY,IBOUND !MODFLOW
	use GWFBASMODULE, only:HNOFLO
      use mf_rt_link, only: rt_active !MODFLOW-(UZF-RT3D) Linkage
      use smrt_parm !SWAT-MODFLOW

      implicit none

	character*30 :: header
	character(len=13) :: hdr_bal(17),hdr_bal_rt(13),hdr_bal_rt_p(12)
	logical :: i_exist
      integer :: i,j,k,n
      
	
	!write message to log file
	write(smrt_log,*)
	write(smrt_log,*) 'hyd_connect --> smrt_initialize'
	write(smrt_log,*) 'reading all smrt connections...'

	!set default values -----------------------------------------------------------------------------------
	mf_interval = 1.0
	mf_ran = .false.
	print *, 'SWAT+/modflow: modflow is being used' !rtb
	
	
	!set up MODFLOW data and allocate arrays --------------------------------------------------------------
      write(smrt_log,*) '     reading modflow input files'
	call mf_read
      write(smrt_log,*) '     completed reading modflow input files'
	
	!read SWAT-MODFLOW connection information -------------------------------------------------------------
      !note: connections between MODFLOW river cells and SWAT+ channels occur smrt_chan_read
      call smrt_hrucell   !hru-cell connections
	call smrt_canl_read !canal-cell connections
	call smrt_drng_read !drain-cell connections
	call smrt_resv_read !reservoir-cell connections
	
      !allocate other necessary arrays ----------------------------------------------------------------------
      allocate(etremain(sp_ob%hru)) !groundwater ET
	allocate(gw_transfer(NCOL,NROW)) !groundwater transfer to soil profile
      etremain = 0.
	gw_transfer = 0.
	
      !open up files that will contain water balance output -------------------------------------------------
      write(smrt_log,*) '     opening groundwater balance output files'
	
	!daily
	open(out_gwbal_dy,file='smrt_out_mf_balance_dy')
      write(out_gwbal_dy,*) 'Daily Groundwater Balance (mm)'
      write(out_gwbal_dy,*) 'Note: gwvolume values are approximate'
      write(out_gwbal_dy,*)
	write(out_gwbal_dy,*) 'watershed area (m2):',(bsn%area_tot_ha*10000.)
	write(out_gwbal_dy,*)
	write(out_gwbal_dy,*) 'gwvolume:      mm  total groundwater volume'
	write(out_gwbal_dy,*) 'gwchange:      mm  groundwater volume change'
	write(out_gwbal_dy,*) 'rech:          mm  recharge to water table'
	write(out_gwbal_dy,*) 'soil:          mm  transfer to soil profile'
	write(out_gwbal_dy,*) 'gwet:          mm  groundwater ET'
	write(out_gwbal_dy,*) 'gwsw:          mm  groundwater --> streams'
	write(out_gwbal_dy,*) 'swgw:          mm  streams --> groundwater'
	write(out_gwbal_dy,*) 'bndr:          mm  constant head boundary'
	write(out_gwbal_dy,*) 'well:          mm  groundwater pumping'
	write(out_gwbal_dy,*) 'mnw2:          mm  groundwater pumping'
	write(out_gwbal_dy,*) 'canl:          mm  canal exhange'
	write(out_gwbal_dy,*) 'drng:          mm  groundwater --> drains'
	write(out_gwbal_dy,*) 'resv:          mm  reservoir exchange'
	write(out_gwbal_dy,*)
	
	hdr_bal = (/character(len=13) :: "jday","mon","day","yr","gwvolume","gwchange","rech",
     &						                                              "soil",
     &                                                          "gwet",
     &                                                          "gwsw",
     &                                                          "swgw",
     &                                                          "bndr",
     &                                                          "well",
     &																													"mnw2",
     &																													"canl",
     &                                                          "drng",
     &                                                          "resv"/)
      write(out_gwbal_dy,120) (hdr_bal(j),j=1,17)
	
	!monthly
	open(out_gwbal_mo,file='smrt_out_mf_balance_mo')
      write(out_gwbal_mo,*) 'Daily Groundwater Balance (mm)'
      write(out_gwbal_mo,*) 'Note: gwvolume values are approximate'
      write(out_gwbal_mo,*)
	write(out_gwbal_mo,*) 'watershed area (m2):',(bsn%area_tot_ha*10000.)
	write(out_gwbal_mo,*)
	write(out_gwbal_mo,*) 'gwvolume:      mm  total groundwater volume'
	write(out_gwbal_mo,*) 'gwchange:      mm  groundwater volume change'
	write(out_gwbal_mo,*) 'rech:          mm  recharge to water table'
	write(out_gwbal_mo,*) 'soil:          mm  transfer to soil profile'
	write(out_gwbal_mo,*) 'gwet:          mm  groundwater ET'
	write(out_gwbal_mo,*) 'gwsw:          mm  groundwater --> streams'
	write(out_gwbal_mo,*) 'swgw:          mm  streams --> groundwater'
	write(out_gwbal_mo,*) 'bndr:          mm  constant head boundary'
	write(out_gwbal_mo,*) 'well:          mm  groundwater pumping'
	write(out_gwbal_mo,*) 'mnw2:          mm  groundwater pumping'
	write(out_gwbal_mo,*) 'canl:          mm  canal exchange'
	write(out_gwbal_mo,*) 'drng:          mm  groundwater --> drains'
	write(out_gwbal_mo,*) 'resv:          mm  reservoir exchange'
	write(out_gwbal_mo,*)
	hdr_bal = (/character(len=13) :: "jday","mon","day","yr","gwvolume","gwchange","rech",
     &						                                              "soil",
     &                                                          "gwet",
     &                                                          "gwsw",
     &                                                          "swgw",
     &                                                          "bndr",
     &                                                          "well",
     &																													"mnw2",
     &                                                          "canl",
     &                                                          "drng",
     &                                                          "resv"/)
      write(out_gwbal_mo,120) (hdr_bal(j),j=1,17)

	!yearly
	open(out_gwbal_yr,file='smrt_out_mf_balance_yr')
      write(out_gwbal_yr,*) 'Daily Groundwater Balance (mm)'
      write(out_gwbal_yr,*) 'Note: gwvolume values are approximate'
      write(out_gwbal_yr,*)
	write(out_gwbal_yr,*) 'watershed area (m2):',(bsn%area_tot_ha*10000.)
	write(out_gwbal_yr,*)
	write(out_gwbal_yr,*) 'gwvolume:      mm  total groundwater volume'
	write(out_gwbal_yr,*) 'gwchange:      mm  groundwater volume change'
	write(out_gwbal_yr,*) 'rech:          mm  recharge to water table'
	write(out_gwbal_yr,*) 'soil:          mm  transfer to soil profile'
	write(out_gwbal_yr,*) 'gwet:          mm  groundwater ET'
	write(out_gwbal_yr,*) 'gwsw:          mm  groundwater --> streams'
	write(out_gwbal_yr,*) 'swgw:          mm  streams --> groundwater'
	write(out_gwbal_yr,*) 'bndr:          mm  constant head boundary'
	write(out_gwbal_yr,*) 'well:          mm  groundwater pumping'
	write(out_gwbal_yr,*) 'mnw2:          mm  groundwater pumping'
	write(out_gwbal_yr,*) 'canl:          mm  canal exchange'
	write(out_gwbal_yr,*) 'drng:          mm  groundwater --> drains'
	write(out_gwbal_yr,*) 'resv:          mm  reservoir exchange'
	write(out_gwbal_yr,*)
	hdr_bal = (/character(len=13) :: "jday","mon","day","yr","gwvolume","gwchange","rech",
     &						                                              "soil",
     &                                                          "gwet",
     &                                                          "gwsw",
     &                                                          "swgw",
     &                                                          "bndr",
     &                                                          "well",
     &																													"mnw2",
     &																													"canl",
     &                                                          "drng",
     &                                                          "resv"/)
      write(out_gwbal_yr,120) (hdr_bal(j),j=1,16)

	!average annual
	open(out_gwbal_aa,file='smrt_out_mf_balance_aa')
      write(out_gwbal_aa,*) 'Daily Groundwater Balance (mm)'
      write(out_gwbal_aa,*) 'Note: gwvolume values are approximate'
      write(out_gwbal_aa,*)
	write(out_gwbal_aa,*) 'watershed area (m2):',(bsn%area_tot_ha*10000.)
	write(out_gwbal_aa,*)
	write(out_gwbal_aa,*) 'gwvolume:      mm  total groundwater volume'
	write(out_gwbal_aa,*) 'gwchange:      mm  groundwater volume change'
	write(out_gwbal_aa,*) 'rech:          mm  recharge to water table'
	write(out_gwbal_aa,*) 'soil:          mm  transfer to soil profile'
	write(out_gwbal_aa,*) 'gwet:          mm  groundwater ET'
	write(out_gwbal_aa,*) 'gwsw:          mm  groundwater --> streams'
	write(out_gwbal_aa,*) 'swgw:          mm  streams --> groundwater'
	write(out_gwbal_aa,*) 'bndr:          mm  constant head boundary'
	write(out_gwbal_aa,*) 'well:          mm  groundwater pumping'
	write(out_gwbal_aa,*) 'mnw2:          mm  groundwater pumping'
	write(out_gwbal_aa,*) 'canl:          mm  canal exchange'
	write(out_gwbal_aa,*) 'drng:          mm  groundwater --> drains'
	write(out_gwbal_aa,*) 'resv:          mm  reservoir exchange'
	write(out_gwbal_aa,*)
	hdr_bal = (/character(len=13) :: "jday","mon","day","yr","gwvolume","gwchange","rech",
     &						                                              "soil",
     &                                                          "gwet",
     &                                                          "gwsw",
     &                                                          "swgw",
     &                                                          "bndr",
     &                                                          "well",
     &																													"mnw2",
     &                                                          "canl",
     &                                                          "drng",
     &                                                          "resv"/)
      write(out_gwbal_aa,120) (hdr_bal(j),j=1,16)
	
	
      !open up files that will contain head and flux output -------------------------------------------------
      write(smrt_log,*) '     opening groundwater flux output files' 
	
	!groundwater head
      open(mf_head_mo,file='smrt_out_mf_gwhead_mo')
      write(mf_head_mo,*) 'Monthly Averaged Head Values for MODFLOW'
	write(mf_head_mo,*) 'Units: meter'
      write(mf_head_mo,*)
      open(mf_head_yr,file='smrt_out_mf_gwhead_yr')
      write(mf_head_yr,*) 'Yearly Averaged Head Values (m) - MODFLOW'
	write(mf_head_yr,*) 'Units: meter'
      write(mf_head_yr,*)
	
      !gw recharge
      open(mf_rech_mo,file='smrt_out_mf_rech_mo')
      write(mf_rech_mo,*) 'Monthly Averaged Recharge'
	write(mf_rech_mo,*) 'Units: m3/day'
      write(mf_rech_mo,*)
      open(mf_rech_yr,file='smrt_out_mf_rech_yr')
      write(mf_rech_yr,*) 'Yearly Averaged Recharge'
	write(mf_rech_yr,*) 'Units: m3/day'
      write(mf_rech_yr,*)

      !gw-->soil transfer
      open(mf_soil_mo,file='smrt_out_mf_soil_mo')
      write(mf_soil_mo,*) 'Monthly Averaged GW-->Soil'
	write(mf_soil_mo,*) 'Units: m3/day'
      write(mf_soil_mo,*)
      open(mf_soil_yr,file='smrt_out_mf_soil_yr')
      write(mf_soil_yr,*) 'Yearly Averaged GW-->Soil'
	write(mf_soil_yr,*) 'Units: m3/day'
      write(mf_soil_yr,*)
	
	!gw ET
      open(mf_gwet_mo,file='smrt_out_mf_gwet_mo')
      write(mf_gwet_mo,*) 'Monthly Averaged GW ET'
	write(mf_gwet_mo,*) 'Units: m3/day'
      write(mf_gwet_mo,*)
      open(mf_gwet_yr,file='smrt_out_mf_gwet_yr')
      write(mf_gwet_yr,*) 'Yearly Averaged GW ET'
	write(mf_gwet_yr,*) 'Units: m3/day'
      write(mf_gwet_yr,*)
	
	!gw-channel exchange
      open(mf_gwsw_mo,file='smrt_out_mf_gwsw_mo')
      write(mf_gwsw_mo,*) 'Monthly Averaged GW-Channel Exchange'
	write(mf_gwsw_mo,*) 'Units: m3/day'
	write(mf_gwsw_mo,*) 'Positive: channel seepage to aquifer'
	write(mf_gwsw_mo,*) 'Negative: groundwater discharge to channel'
      write(mf_gwsw_mo,*)
      open(mf_gwsw_yr,file='smrt_out_mf_gwsw_yr')
      write(mf_gwsw_yr,*) 'Yearly Averaged GW-Channel Exchange'
	write(mf_gwsw_yr,*) 'Units: m3/day'
	write(mf_gwsw_yr,*) 'Positive: channel seepage to aquifer'
	write(mf_gwsw_yr,*) 'Negative: groundwater discharge to channel'
      write(mf_gwsw_yr,*)
	
	!gw pumping - well package
      open(mf_well_mo,file='smrt_out_mf_well_mo')
      write(mf_well_mo,*) 'Monthly Averaged Pumping - Well Package'
	write(mf_well_mo,*) 'Units: m3/day'
      write(mf_well_mo,*)
      open(mf_well_yr,file='smrt_out_mf_well_yr')
      write(mf_well_yr,*) 'Yearly Averaged Pumping - Well Package'
	write(mf_well_yr,*) 'Units: m3/day'
      write(mf_well_yr,*)
	
	!gw pumping - mnw2 package
      open(mf_mnw2_mo,file='smrt_out_mf_mnw2_mo')
      write(mf_mnw2_mo,*) 'Monthly Averaged Pumping - MNW2 Package'
	write(mf_mnw2_mo,*) 'Units: m3/day'
      write(mf_mnw2_mo,*)
      open(mf_mnw2_yr,file='smrt_out_mf_mnw2_yr')
      write(mf_mnw2_yr,*) 'Yearly Averaged Pumping - MNW2 Package'
	write(mf_mnw2_yr,*) 'Units: m3/day'
      write(mf_mnw2_yr,*)
	
	!gw-canal exchange
      open(mf_canl_mo,file='smrt_out_mf_canl_mo')
      write(mf_canl_mo,*) 'Monthly Averaged GW-Canal Exchange'
	write(mf_canl_mo,*) 'Units: m3/day'
      write(mf_canl_mo,*)
      open(mf_canl_yr,file='smrt_out_mf_canl_yr')
      write(mf_canl_yr,*) 'Yearly Averaged GW-Canal Exchange'
	write(mf_canl_yr,*) 'Units: m3/day'
      write(mf_canl_yr,*)
	
	!gw drainage
      open(mf_drng_mo,file='smrt_out_mf_drng_mo')
      write(mf_drng_mo,*) 'Monthly Averaged Drainage'
	write(mf_drng_mo,*) 'Units: m3/day'
      write(mf_drng_mo,*)
      open(mf_drng_yr,file='smrt_out_mf_drng_yr')
      write(mf_drng_yr,*) 'Yearly Averaged Drainage'
	write(mf_drng_yr,*) 'Units: m3/day'
      write(mf_drng_yr,*)
	
	!gw-reservoir exchange
      open(mf_resv_mo,file='smrt_out_mf_resv_mo')
      write(mf_resv_mo,*) 'Monthly Averaged GW-Reservoir Exchange'
	write(mf_resv_mo,*) 'Units: m3/day'
	write(mf_resv_mo,*) 'Positive: reservoir seepage to aquifer'
	write(mf_resv_mo,*) 'Negative: groundwater discharge to reservoir'
      write(mf_resv_mo,*)
      open(mf_resv_yr,file='smrt_out_mf_resv_yr')
      write(mf_resv_yr,*) 'Yearly Averaged GW-Reservoir Exchange'
	write(mf_resv_yr,*) 'Units: m3/day'
	write(mf_resv_yr,*) 'Positive: reservoir seepage to aquifer'
	write(mf_resv_yr,*) 'Negative: groundwater discharge to reservoir'
      write(mf_resv_yr,*)
      

	!allocate arrays for writing out monthly and annual values --------------------------------------------
	!head
	allocate(smrt_head_tot_mo(ncol,nrow,nlay))
      allocate(smrt_head_tot_yr(ncol,nrow,nlay))
	smrt_head_tot_mo = 0.
	smrt_head_tot_yr = 0.
	!fluxes
	allocate(smrt_rech_tot_mo(ncol,nrow))
      allocate(smrt_rech_tot_yr(ncol,nrow))
	smrt_rech_tot_mo = 0.
	smrt_rech_tot_yr = 0.
	allocate(smrt_soil_tot_mo(ncol,nrow))
      allocate(smrt_soil_tot_yr(ncol,nrow))
	smrt_soil_tot_mo = 0.
	smrt_soil_tot_yr = 0.
	allocate(smrt_gwet_tot_mo(ncol,nrow))
      allocate(smrt_gwet_tot_yr(ncol,nrow))
	smrt_gwet_tot_mo = 0.
	smrt_gwet_tot_yr = 0.
	allocate(smrt_gwsw_tot_mo(ncol,nrow))
      allocate(smrt_gwsw_tot_yr(ncol,nrow))
	smrt_gwsw_tot_mo = 0.
	smrt_gwsw_tot_yr = 0.
	allocate(smrt_well_tot_mo(ncol,nrow))
      allocate(smrt_well_tot_yr(ncol,nrow))
	smrt_well_tot_mo = 0.
	smrt_well_tot_yr = 0.
	allocate(smrt_mnw2_tot_mo(ncol,nrow))
      allocate(smrt_mnw2_tot_yr(ncol,nrow))
	smrt_mnw2_tot_mo = 0.
	smrt_mnw2_tot_yr = 0.
	allocate(smrt_canl_tot_mo(ncol,nrow))
      allocate(smrt_canl_tot_yr(ncol,nrow))
	smrt_canl_tot_mo = 0.
	smrt_canl_tot_yr = 0.
	allocate(smrt_drng_tot_mo(ncol,nrow))
      allocate(smrt_drng_tot_yr(ncol,nrow))
	smrt_drng_tot_mo = 0.
	smrt_drng_tot_yr = 0.
      allocate(smrt_resv_tot_mo(ncol,nrow))
      allocate(smrt_resv_tot_yr(ncol,nrow))
	smrt_resv_tot_mo = 0.
	smrt_resv_tot_yr = 0.
	
	
	
	
      !set up RT3D data and allocate arrays -----------------------------------------------------------------
      !determine if RT3D is active
	write(smrt_log,*)
	write(smrt_log,*) 'check if rt3d is active for nutrient transport'
	inquire (file='modflow_rt.nam',exist=i_exist)
      if(i_exist) then
	  rt_active = 1
		print *, 'SWAT+/modflow: rt3d is being used for nutrients' !rtb
		print *
		write(smrt_log,*) '     modflow_rt.nam found: rt3d is assumed active'
	else
	  rt_active = 0
	endif
	if(rt_active.eq.1) then
        	   
	  !read data and allocation arrays ----------------------------------------------------------
	  write(smrt_log,*) '     read rt3d input files'
	  call rt_read
		write(smrt_log,*) '     finished reading rt3d input files'
		write(smrt_log,*) '     open concentration and loading output files'
		
		!NO3-N basin wide fluxes ------------------------------------------------------------------
		!daily
		open(out_rtbal_dy_n,file='smrt_out_rt_balance_dy_no3')
		write(out_rtbal_dy_n,*) 'Daily NO3-N Loads (kg)'
		write(out_rtbal_dy_n,*)
		write(out_rtbal_dy_n,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_dy_n,*)
		write(out_rtbal_dy_n,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_dy_n,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_dy_n,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_dy_n,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_dy_n,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_dy_n,*) 'canl:          kg  canal exhange'
		write(out_rtbal_dy_n,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_dy_n,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_dy_n,*) 'rctn:          kg  denitrification'
		write(out_rtbal_dy_n,*)
		hdr_bal_rt = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                             "soil",
     &                                         "gwsw",
     &                                         "swgw",
     &                                         "well",
     &																         "canl",
     &                                         "drng",
     &                                         "resv",
     &                                         "rctn"/)
		write(out_rtbal_dy_n,120) (hdr_bal_rt(j),j=1,13)
		!monthly
		open(out_rtbal_mo_n,file='smrt_out_rt_balance_mo_no3')
		write(out_rtbal_mo_n,*) 'Monthly NO3-N Loads (kg)'
		write(out_rtbal_mo_n,*)
		write(out_rtbal_mo_n,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_mo_n,*)
		write(out_rtbal_mo_n,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_mo_n,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_mo_n,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_mo_n,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_mo_n,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_mo_n,*) 'canl:          kg  canal exhange'
		write(out_rtbal_mo_n,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_mo_n,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_mo_n,*) 'rctn:          kg  denitrification'
		write(out_rtbal_mo_n,*)
		hdr_bal_rt = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                             "soil",
     &                                         "gwsw",
     &                                         "swgw",
     &                                         "well",
     &																         "canl",
     &                                         "drng",
     &                                         "resv",
     &                                         "rctn"/)
		write(out_rtbal_mo_n,120) (hdr_bal_rt(j),j=1,13)
		!yearly
		open(out_rtbal_yr_n,file='smrt_out_rt_balance_yr_no3')
		write(out_rtbal_yr_n,*) 'Yearly NO3-N Loads (kg)'
		write(out_rtbal_yr_n,*)
		write(out_rtbal_yr_n,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_yr_n,*)
		write(out_rtbal_yr_n,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_yr_n,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_yr_n,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_yr_n,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_yr_n,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_yr_n,*) 'canl:          kg  canal exhange'
		write(out_rtbal_yr_n,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_yr_n,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_yr_n,*) 'rctn:          kg  denitrification'
		write(out_rtbal_yr_n,*)
		hdr_bal_rt = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                             "soil",
     &                                         "gwsw",
     &                                         "swgw",
     &                                         "well",
     &																         "canl",
     &                                         "drng",
     &                                         "resv",
     &                                         "rctn"/)
		write(out_rtbal_yr_n,120) (hdr_bal_rt(j),j=1,13)
		!average annual
		open(out_rtbal_aa_n,file='smrt_out_rt_balance_aa_no3')
		write(out_rtbal_aa_n,*) 'Average Annual NO3-N Loads (kg)'
		write(out_rtbal_aa_n,*)
		write(out_rtbal_aa_n,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_aa_n,*)
		write(out_rtbal_aa_n,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_aa_n,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_aa_n,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_aa_n,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_aa_n,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_aa_n,*) 'canl:          kg  canal exhange'
		write(out_rtbal_aa_n,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_aa_n,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_aa_n,*) 'rctn:          kg  denitrification'
		write(out_rtbal_aa_n,*)
		hdr_bal_rt = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                             "soil",
     &                                         "gwsw",
     &                                         "swgw",
     &                                         "well",
     &																         "canl",
     &                                         "drng",
     &                                         "resv",
     &                                         "rctn"/)
		write(out_rtbal_aa_n,120) (hdr_bal_rt(j),j=1,13)
		
		!P basin wide fluxes ----------------------------------------------------------------------
		!daily
		open(out_rtbal_dy_p,file='smrt_out_rt_balance_dy_p')
		write(out_rtbal_dy_p,*) 'Daily P Loads (kg)'
		write(out_rtbal_dy_p,*)
		write(out_rtbal_dy_p,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_dy_p,*)
		write(out_rtbal_dy_p,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_dy_p,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_dy_p,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_dy_p,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_dy_p,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_dy_p,*) 'canl:          kg  canal exhange'
		write(out_rtbal_dy_p,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_dy_p,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_dy_p,*)
		hdr_bal_rt_p = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                               "soil",
     &                                           "gwsw",
     &                                           "swgw",
     &                                           "well",
     &																           "canl",
     &                                           "drng",
     &                                           "resv"/)
		write(out_rtbal_dy_p,120) (hdr_bal_rt_p(j),j=1,12)
		!monthly
		open(out_rtbal_mo_p,file='smrt_out_rt_balance_mo_p')
		write(out_rtbal_mo_p,*) 'Monthly P Loads (kg)'
		write(out_rtbal_mo_p,*)
		write(out_rtbal_mo_p,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_mo_p,*)
		write(out_rtbal_mo_p,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_mo_p,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_mo_p,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_mo_p,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_mo_p,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_mo_p,*) 'canl:          kg  canal exhange'
		write(out_rtbal_mo_p,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_mo_p,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_mo_p,*)
		hdr_bal_rt_p = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                               "soil",
     &                                           "gwsw",
     &                                           "swgw",
     &                                           "well",
     &																           "canl",
     &                                           "drng",
     &                                           "resv"/)
		write(out_rtbal_mo_p,120) (hdr_bal_rt_p(j),j=1,12)
		!yearly
		open(out_rtbal_yr_p,file='smrt_out_rt_balance_yr_p')
		write(out_rtbal_yr_p,*) 'Yearly P Loads (kg)'
		write(out_rtbal_yr_p,*)
		write(out_rtbal_yr_p,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_yr_p,*)
		write(out_rtbal_yr_p,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_yr_p,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_yr_p,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_yr_p,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_yr_p,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_yr_p,*) 'canl:          kg  canal exhange'
		write(out_rtbal_yr_p,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_yr_p,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_yr_p,*)
		hdr_bal_rt_p = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                               "soil",
     &                                           "gwsw",
     &                                           "swgw",
     &                                           "well",
     &																           "canl",
     &                                           "drng",
     &                                           "resv"/)
		write(out_rtbal_yr_p,120) (hdr_bal_rt_p(j),j=1,12)
		!average annual
		open(out_rtbal_aa_p,file='smrt_out_rt_balance_aa_p')
		write(out_rtbal_aa_p,*) 'Average Annual P Loads (kg)'
		write(out_rtbal_aa_p,*)
		write(out_rtbal_aa_p,*) 'watershed area (m2):',(bsn%area_tot_ha*10000)
		write(out_rtbal_aa_p,*)
		write(out_rtbal_aa_p,*) 'rech:          kg  recharge to water table'
		write(out_rtbal_aa_p,*) 'soil:          kg  transfer to soil profile'
		write(out_rtbal_aa_p,*) 'gwsw:          kg  groundwater --> streams'
		write(out_rtbal_aa_p,*) 'swgw:          kg  streams --> groundwater'
		write(out_rtbal_aa_p,*) 'well:          kg  groundwater pumping'
		write(out_rtbal_aa_p,*) 'canl:          kg  canal exhange'
		write(out_rtbal_aa_p,*) 'drng:          kg  groundwater --> drains'
		write(out_rtbal_aa_p,*) 'resv:          kg  reservoir exchange'
		write(out_rtbal_aa_p,*)
		hdr_bal_rt_p = (/character(len=13) :: "jday","mon","day","yr","rech",
     &						                               "soil",
     &                                           "gwsw",
     &                                           "swgw",
     &                                           "well",
     &																           "canl",
     &                                           "drng",
     &                                           "resv"/)
		write(out_rtbal_aa_p,120) (hdr_bal_rt_p(j),j=1,12)
		
		!groundwater concentration files ----------------------------------------------------------
		open(rt_conc_mo,file='smrt_out_rt_conc_mo')
		write(rt_conc_mo,*) 'Monthly Averaged Concentration'
		write(rt_conc_mo,*) 'Units: g/m3 = mg/L'
		write(rt_conc_mo,*)
		open(rt_conc_yr,file='smrt_out_rt_conc_yr')
		write(rt_conc_yr,*) 'Yearly Averaged Concentration'
		write(rt_conc_yr,*) 'Units: g/m3 = mg/L'
		write(rt_conc_yr,*)
	  
		!nutrient flux files ----------------------------------------------------------------------
		!gw recharge
		open(rt_rech_mo,file='smrt_out_rt_rech_mo')
		write(rt_rech_mo,*) 'Monthly Averaged Recharge Load'
		write(rt_rech_mo,*) 'Units: kg/day'
		write(rt_rech_mo,*)
		open(rt_rech_yr,file='smrt_out_rt_rech_yr')
		write(rt_rech_yr,*) 'Yearly Averaged Recharge Load'
		write(rt_rech_yr,*) 'Units: kg/day'
		write(rt_rech_yr,*)

		!gw-->soil transfer
		open(rt_soil_mo,file='smrt_out_rt_soil_mo')
		write(rt_soil_mo,*) 'Monthly Averaged GW-->Soil Load'
		write(rt_soil_mo,*) 'Units: kg/day'
		write(rt_soil_mo,*)
		open(rt_soil_yr,file='smrt_out_rt_soil_yr')
		write(rt_soil_yr,*) 'Yearly Averaged GW-->Soil Load'
		write(rt_soil_yr,*) 'Units: kg/day'
		write(rt_soil_yr,*)
	
		!gw-channel exchange
		open(rt_gwsw_mo,file='smrt_out_rt_gwsw_mo')
		write(rt_gwsw_mo,*) 'Monthly Averaged GW-Channel Exchange Load'
		write(rt_gwsw_mo,*) 'Units: kg/day'
		write(rt_gwsw_mo,*) 'Positive: channel seepage to aquifer'
		write(rt_gwsw_mo,*) 'Negative: groundwater discharge to channel'
		write(rt_gwsw_mo,*)
		open(rt_gwsw_yr,file='smrt_out_rt_gwsw_yr')
		write(rt_gwsw_yr,*) 'Yearly Averaged GW-Channel Exchange Load'
		write(rt_gwsw_yr,*) 'Units: kg/day'
		write(rt_gwsw_yr,*) 'Positive: channel seepage to aquifer'
		write(rt_gwsw_yr,*) 'Negative: groundwater discharge to channel'
		write(rt_gwsw_yr,*)
	
		!gw pumping
		open(rt_well_mo,file='smrt_out_rt_well_mo')
		write(rt_well_mo,*) 'Monthly Averaged Pumping Load - Well Package'
		write(rt_well_mo,*) 'Units: kg/day'
		write(rt_well_mo,*)
		open(rt_well_yr,file='smrt_out_rt_well_yr')
		write(rt_well_yr,*) 'Yearly Averaged Pumping Load - Well Package'
		write(rt_well_yr,*) 'Units: kg/day'
		write(rt_well_yr,*)
	
		!gw-canal exchange
		open(rt_canl_mo,file='smrt_out_rt_canl_mo')
		write(rt_canl_mo,*) 'Monthly Averaged GW-Canal Load'
		write(rt_canl_mo,*) 'Units: kg/day'
		write(rt_canl_mo,*)
		open(rt_canl_yr,file='smrt_out_rt_canl_yr')
		write(rt_canl_yr,*) 'Yearly Averaged GW-Canal Load'
		write(rt_canl_yr,*) 'Units: kg/day'
		write(rt_canl_yr,*)
		
		!gw drainage
		open(rt_drng_mo,file='smrt_out_rt_drng_mo')
		write(rt_drng_mo,*) 'Monthly Averaged Drainage Load'
		write(rt_drng_mo,*) 'Units: kg/day'
		write(rt_drng_mo,*)
		open(rt_drng_yr,file='smrt_out_rt_drng_yr')
		write(rt_drng_yr,*) 'Yearly Averaged Drainage Load'
		write(rt_drng_yr,*) 'Units: kg/day'
		write(rt_drng_yr,*)
	
		!gw-reservoir exchange
		open(rt_resv_mo,file='smrt_out_rt_resv_mo')
		write(rt_resv_mo,*) 'Monthly Averaged GW-Reservoir Exchange Load'
		write(rt_resv_mo,*) 'Units: kg/day'
		write(rt_resv_mo,*) 'Positive: reservoir seepage to aquifer'
		write(rt_resv_mo,*) 'Negative: groundwater discharge to reservoir'
		write(rt_resv_mo,*)
		open(rt_resv_yr,file='smrt_out_rt_resv_yr')
		write(rt_resv_yr,*) 'Yearly Averaged GW-Reservoir Exchange Load'
		write(rt_resv_yr,*) 'Units: kg/day'
		write(rt_resv_yr,*) 'Positive: reservoir seepage to aquifer'
		write(rt_resv_yr,*) 'Negative: groundwater discharge to reservoir'
		write(rt_resv_yr,*)
		
		!concentration
		allocate(smrt_conc_tot_mo(ncol,nrow,nlay,2))
		allocate(smrt_conc_tot_yr(ncol,nrow,nlay,2))
		smrt_conc_tot_mo = 0.
		smrt_conc_tot_yr = 0.
		!loads
		allocate(smrt_rt_rech_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_rech_tot_yr(ncol,nrow,2))
		smrt_rt_rech_tot_mo = 0.
		smrt_rt_rech_tot_yr = 0.
		allocate(smrt_rt_soil_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_soil_tot_yr(ncol,nrow,2))
		smrt_rt_soil_tot_mo = 0.
		smrt_rt_soil_tot_yr = 0.
		allocate(smrt_rt_gwsw_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_gwsw_tot_yr(ncol,nrow,2))
		smrt_rt_gwsw_tot_mo = 0.
		smrt_rt_gwsw_tot_yr = 0.
		allocate(smrt_rt_well_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_well_tot_yr(ncol,nrow,2))
		smrt_rt_well_tot_mo = 0.
		smrt_rt_well_tot_yr = 0.
		allocate(smrt_rt_canl_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_canl_tot_yr(ncol,nrow,2))
		smrt_rt_canl_tot_mo = 0.
		smrt_rt_canl_tot_yr = 0.
		allocate(smrt_rt_drng_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_drng_tot_yr(ncol,nrow,2))
		smrt_rt_drng_tot_mo = 0.
		smrt_rt_drng_tot_yr = 0.
		allocate(smrt_rt_resv_tot_mo(ncol,nrow,2))
		allocate(smrt_rt_resv_tot_yr(ncol,nrow,2))
		smrt_rt_resv_tot_mo = 0.
		smrt_rt_resv_tot_yr = 0.
		
		!allocate other necessary arrays
		allocate(rt_rech_load(ncol,nrow,2))
		allocate(rt_soil_load(ncol,nrow,2))
		allocate(rt_gwsw_load(ncol,nrow,2))
		allocate(rt_well_load(ncol,nrow,2))
		allocate(rt_canl_load(ncol,nrow,2))
		allocate(rt_drng_load(ncol,nrow,2))
		allocate(rt_resv_load(ncol,nrow,2))
		allocate(rt_rctn_load(ncol,nrow,nlay))
		rt_rech_load = 0.
		rt_soil_load = 0.
		rt_gwsw_load = 0.
		rt_well_load = 0.
		rt_canl_load = 0.
		rt_drng_load = 0.
		rt_resv_load = 0.
		rt_rctn_load = 0.
		allocate(mf_hru_soil(sp_ob%hru,20,2))
		allocate(mf_percsol(sp_ob%hru,2))
		allocate(mf_rechsol(sp_ob%hru,2))
		mf_hru_soil = 0.
		mf_percsol = 0.
		mf_rechsol = 0.
		
      endif



120   format(4x,a8,a8,a8,a4,a16,a19,50(a13))	

      return
      end
