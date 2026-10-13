      subroutine smrt_run

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls conversions from SWAT+ to MODFLOW, then runs MODFLOW
!!    for one day, then converts back the groundwater results from MODFLOW to SWAT+

        use hru_module
        use hydrograph_module
        use basin_module
        use GLOBAL, only: mf_interval,mf_ran !MODFLOW
        use mf_rt_link, only: rt_active !MODFLOW-(UZF-RT3D) Linkage
        use smrt_parm !SWAT-MODFLOW
		use time_module
        
        implicit none
        
        integer m
        
		
		!message to log file
        write(smrt_log,*)
		write(smrt_log,*) 'running smrt for time:',time%yrc,time%day
		
        !convert SWAT units into MODFLOW units
        call smrt_conversion2mf
          
        !run MODFLOW (water table, river discharge/seepage, pumping); run RT3D if selected
        call mf_run
        
        !convert MODFLOW units back to SWAT units
        call smrt_conversion2swat


      end subroutine smrt_run
