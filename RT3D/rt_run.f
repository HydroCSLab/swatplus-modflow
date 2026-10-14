*     This subroutine runs the RT3D simulator for the time period indicated. Typically,
*     the stress period and flow time step are the same as for SWAT (1 day).

*     This subroutine is called for each flow time step of MODFLOW

      subroutine rt_run

      use rt_global
      use mf_rt_link
      use GLOBAL,only:ITMUNI,NROW,NCOL,NLAY,MF_leapyr !MODFLOW

      implicit none
      integer   i,j,k,N,L,dum
      common   /GCGIDX/L(19) !for implicit shceme
           


*     This subroutine is called within the MODFLOW flow time step loop. Thus, it
*     calculates the reactive transport of chemical species for 1 day (to match
*     the SWAT step size)
      
      !Calculate the dispersion coefficients for the current flow time step
      if(TRNOP(2)) call dsp_coeff

      dt = rt_delt !(calculated in mf_run)
      t1 = t2
      t2 = t2 + dt
      
	!zero out mass transfer arrays
      rt_rivmass = 0.
      rt_drnmass = 0.
	rt_resmass = 0.
	
      !loop through the transport time steps for the flow time step
      do ntrans=1,5000

        !advance once transport step
        call btn_ts

        !for each mobile species, compute change in concentration due to advection, 
        !dispersion, and source-sink mixing. this is done implicitly
        do icomp=1,ncomp
          if(phase(icomp).eq.1) then !species is mobile
            call rt_implicit
          endif
        enddo
          
        !compute change in concentration due to chemical reactions
        if(trnop(4)) call rct_solve

        !compute mass budgets and output results for the transport step
        do icomp=1,ncomp
          call btn_budget !mass budgets  
        enddo 

        !if end of flow time step, exit
        if(time2.ge.t2) goto 100

      enddo !go to next transport time step


100   return
	end
