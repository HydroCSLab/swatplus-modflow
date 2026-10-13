      subroutine sd_channel_output (ichan)
    
      use sd_channel_module
      use basin_module
      use time_module
      use hydrograph_module
      use water_body_module
      
      implicit none
      
      integer, intent (in) :: ichan         !             |
      integer :: iob                        !             |
      integer :: ii                         !             
      real :: const                         !             |

      iob = sp_ob1%chandeg + ichan - 1

      !ch_stor_m(ichan) = ch_stor_m(ichan) + ch_stor(ichan)
      ch_in_m(ichan) = ch_in_m(ichan) + ch_in_d(ichan)
      ch_out_m(ichan) = ch_out_m(ichan) + ch_out_d(ichan)
      ch_wat_m(ichan) = ch_wat_m(ichan) + ch_wat_d(ichan)
			
			!rtb modflow
			!monthy flow rate
			ch_gwex_m(ichan)%gwsw = ch_gwex_m(ichan)%gwsw + ch_gwex_d(ichan)%gwsw
			ch_gwex_m(ichan)%gwsw_n = ch_gwex_m(ichan)%gwsw_n + ch_gwex_d(ichan)%gwsw_n
			ch_gwex_m(ichan)%gwsw_p = ch_gwex_m(ichan)%gwsw_p + ch_gwex_d(ichan)%gwsw_p
			ch_gwex_m(ichan)%canl = ch_gwex_m(ichan)%canl + ch_gwex_d(ichan)%canl
			ch_gwex_m(ichan)%canl_n = ch_gwex_m(ichan)%canl_n + ch_gwex_d(ichan)%canl_n
			ch_gwex_m(ichan)%canl_p = ch_gwex_m(ichan)%canl_p + ch_gwex_d(ichan)%canl_p
			ch_gwex_m(ichan)%drng = ch_gwex_m(ichan)%drng + ch_gwex_d(ichan)%drng
			ch_gwex_m(ichan)%drng_n = ch_gwex_m(ichan)%drng_n + ch_gwex_d(ichan)%drng_n
      ch_gwex_m(ichan)%drng_p = ch_gwex_m(ichan)%drng_p + ch_gwex_d(ichan)%drng_p
			
         if (pco%day_print == "y" .and. time%step > 1 .and. pco%int_day_cur == pco%int_day) then
          do ii = 1, time%step 
            write (2508,101)  time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ii, ob(iob)%name,  &
                     ob(iob)%hyd_flo(1,ii)
          end do
         end if
           
!!!!! daily print
       if (pco%day_print == "y" .and. pco%int_day_cur == pco%int_day) then
        if (pco%sd_chan%d == "y") then
          write (2500,100) time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,       &
            ch_wat_d(ichan)%area_ha, ch_wat_d(ichan)%precip, ch_wat_d(ichan)%evap, ch_wat_d(ichan)%seep,        &                
            ch_stor(ichan), ch_in_d(ichan), ch_out_d(ichan), wtemp, ch_gwex_d(ichan)
           if (pco%csvout == "y") then
             write (2504,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,  &
               ch_wat_d(ichan)%area_ha, ch_wat_d(ichan)%precip, ch_wat_d(ichan)%evap, ch_wat_d(ichan)%seep,     &
               ch_stor(ichan), ch_in_d(ichan), ch_out_d(ichan), wtemp, ch_gwex_d(ichan)
           end if
        end if
      end if

			!rtb modflow
			ch_gwex_d(ichan)%gwsw = 0.
			ch_gwex_d(ichan)%gwsw_n = 0.
			ch_gwex_d(ichan)%gwsw_p = 0.
			ch_gwex_d(ichan)%canl = 0.
			ch_gwex_d(ichan)%canl_n = 0.
			ch_gwex_d(ichan)%canl_p = 0.
			ch_gwex_d(ichan)%drng = 0.
			ch_gwex_d(ichan)%drng_n = 0.
      ch_gwex_d(ichan)%drng_p = 0.
			
!!!!! monthly print
        if (time%end_mo == 1) then
          !ch_stor_y(ichan) = ch_stor_y(ichan) + ch_stor_m(ichan)
          const = float (ndays(time%mo + 1) - ndays(time%mo))
          ch_in_y(ichan) = ch_in_y(ichan) + ch_in_m(ichan)
          ch_out_y(ichan) = ch_out_y(ichan) + ch_out_m(ichan)
          ch_wat_y(ichan) = ch_wat_y(ichan) + ch_wat_m(ichan)
					
					!rtb modflow
					ch_gwex_y(ichan)%gwsw = ch_gwex_y(ichan)%gwsw + ch_gwex_m(ichan)%gwsw
			    ch_gwex_y(ichan)%gwsw_n = ch_gwex_y(ichan)%gwsw_n + ch_gwex_m(ichan)%gwsw_n
			    ch_gwex_y(ichan)%gwsw_p = ch_gwex_y(ichan)%gwsw_p + ch_gwex_m(ichan)%gwsw_p
			    ch_gwex_y(ichan)%canl = ch_gwex_y(ichan)%canl + ch_gwex_m(ichan)%canl
			    ch_gwex_y(ichan)%canl_n = ch_gwex_y(ichan)%canl_n + ch_gwex_m(ichan)%canl_n
			    ch_gwex_y(ichan)%canl_p = ch_gwex_y(ichan)%canl_p + ch_gwex_m(ichan)%canl_p
			    ch_gwex_y(ichan)%drng = ch_gwex_y(ichan)%drng + ch_gwex_m(ichan)%drng
			    ch_gwex_y(ichan)%drng_n = ch_gwex_y(ichan)%drng_n + ch_gwex_m(ichan)%drng_n
          ch_gwex_y(ichan)%drng_p = ch_gwex_y(ichan)%drng_p + ch_gwex_m(ichan)%drng_p
					
          ch_in_m(ichan)%flo = ch_in_m(ichan)%flo / const
          ch_out_m(ichan)%flo = ch_out_m(ichan)%flo / const
          ch_wat_m(ichan) = ch_wat_m(ichan) // const            !! // only divides area (daily average values)
					
          if (pco%sd_chan%m == "y") then
          write (2501,100) time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,  &
            ch_wat_m(ichan)%area_ha, ch_wat_m(ichan)%precip, ch_wat_m(ichan)%evap, ch_wat_m(ichan)%seep,   &
            ch_stor(ichan), ch_in_m(ichan), ch_out_m(ichan), wtemp, ch_gwex_m(ichan)

          if (pco%csvout == "y") then
            write (2505,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,   &
            ch_wat_m(ichan)%area_ha, ch_wat_m(ichan)%precip, ch_wat_m(ichan)%evap, ch_wat_m(ichan)%seep,   &
            ch_stor(ichan), ch_in_m(ichan), ch_out_m(ichan), wtemp, ch_gwex_m(ichan)
          end if
        end if
        ch_in_m(ichan) = chaz
        ch_out_m(ichan) = chaz
        ch_wat_m(ichan) = wbodz
				!rtb modflow
				ch_gwex_m(ichan)%gwsw = 0.
			  ch_gwex_m(ichan)%gwsw_n = 0.
			  ch_gwex_m(ichan)%gwsw_p = 0.
			  ch_gwex_m(ichan)%canl = 0.
			  ch_gwex_m(ichan)%canl_n = 0.
			  ch_gwex_m(ichan)%canl_p = 0.
			  ch_gwex_m(ichan)%drng = 0.
			  ch_gwex_m(ichan)%drng_n = 0.
        ch_gwex_m(ichan)%drng_p = 0.
        end if

!!!!! yearly print
      if (time%end_yr == 1) then
        !ch_stor_a(ichan) = ch_stor_a(ichan) + ch_stor_y(ichan)
        const = time%day_end_yr
        ch_in_y(ichan)%flo = ch_in_y(ichan)%flo / const
        ch_out_y(ichan)%flo = ch_out_y(ichan)%flo / const
        ch_in_a(ichan) = ch_in_a(ichan) + ch_in_y(ichan)
        ch_out_a(ichan) = ch_out_a(ichan) + ch_out_y(ichan)
        ch_wat_a(ichan) = ch_wat_a(ichan) + ch_wat_y(ichan)
        
				!rtb modflow
				ch_gwex_a(ichan)%gwsw = ch_gwex_a(ichan)%gwsw + ch_gwex_y(ichan)%gwsw
			  ch_gwex_a(ichan)%gwsw_n = ch_gwex_a(ichan)%gwsw_n + ch_gwex_y(ichan)%gwsw_n
			  ch_gwex_a(ichan)%gwsw_p = ch_gwex_a(ichan)%gwsw_p + ch_gwex_y(ichan)%gwsw_p
			  ch_gwex_a(ichan)%canl = ch_gwex_a(ichan)%canl + ch_gwex_y(ichan)%canl
			  ch_gwex_a(ichan)%canl_n = ch_gwex_a(ichan)%canl_n + ch_gwex_y(ichan)%canl_n
			  ch_gwex_a(ichan)%canl_p = ch_gwex_a(ichan)%canl_p + ch_gwex_y(ichan)%canl_p
			  ch_gwex_a(ichan)%drng = ch_gwex_a(ichan)%drng + ch_gwex_y(ichan)%drng
			  ch_gwex_a(ichan)%drng_n = ch_gwex_a(ichan)%drng_n + ch_gwex_y(ichan)%drng_n
        ch_gwex_a(ichan)%drng_p = ch_gwex_a(ichan)%drng_p + ch_gwex_y(ichan)%drng_p
				
        !ch_stor_y(ichan) = ch_stor_y(ichan) / const     !! all storage variables are averages
        ch_wat_y(ichan) = ch_wat_y(ichan) // const      !! // only divides area (daily average values)
          
        if (pco%sd_chan%y == "y") then 
          write (2502,100) time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,   &
            ch_wat_y(ichan)%area_ha, ch_wat_y(ichan)%precip, ch_wat_y(ichan)%evap, ch_wat_y(ichan)%seep,    &
            ch_stor(ichan), ch_in_y(ichan), ch_out_y(ichan), wtemp, ch_gwex_y(ichan)
          if (pco%csvout == "y") then
           write (2506,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,    &
            ch_wat_y(ichan)%area_ha, ch_wat_y(ichan)%precip, ch_wat_y(ichan)%evap, ch_wat_y(ichan)%seep,    &
            ch_stor(ichan), ch_in_y(ichan), ch_out_y(ichan), wtemp, ch_gwex_y(ichan)
          end if
        end if
      end if

!!!!! average annual print
      if (time%end_sim == 1) then
        !ch_stor_a(ichan) = ch_stor_a(ichan) / time%yrs_prt      !! all storage variables (averaged) must be divided by years
        ch_in_a(ichan) = ch_in_a(ichan) / time%yrs_prt          !! all inflow and outflow varaibles (summed) are divided by years
        ch_out_a(ichan) = ch_out_a(ichan) / time%yrs_prt
        ch_wat_a(ichan) = ch_wat_a(ichan) / time%yrs_prt        !! all summed variables divided by years
        ch_wat_a(ichan) = ch_wat_a(ichan) // time%yrs_prt       !! all averaged variables divided by years
        
				!rtb modflow
				ch_gwex_a(ichan)%gwsw = ch_gwex_a(ichan)%gwsw / time%yrs_prt  
			  ch_gwex_a(ichan)%gwsw_n = ch_gwex_a(ichan)%gwsw_n / time%yrs_prt  
			  ch_gwex_a(ichan)%gwsw_p = ch_gwex_a(ichan)%gwsw_p / time%yrs_prt  
			  ch_gwex_a(ichan)%canl = ch_gwex_a(ichan)%canl / time%yrs_prt  
			  ch_gwex_a(ichan)%canl_n = ch_gwex_a(ichan)%canl_n / time%yrs_prt  
			  ch_gwex_a(ichan)%canl_p = ch_gwex_a(ichan)%canl_p / time%yrs_prt  
			  ch_gwex_a(ichan)%drng = ch_gwex_a(ichan)%drng / time%yrs_prt  
			  ch_gwex_a(ichan)%drng_n = ch_gwex_a(ichan)%drng_n / time%yrs_prt  
        ch_gwex_a(ichan)%drng_p = ch_gwex_a(ichan)%drng_p / time%yrs_prt  
				
        if (pco%sd_chan%a == "y") then
        write (2503,100) time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,  &
          ch_wat_a(ichan)%area_ha, ch_wat_a(ichan)%precip, ch_wat_a(ichan)%evap, ch_wat_a(ichan)%seep,   &
          ch_stor(ichan), ch_in_a(ichan), ch_out_a(ichan), wtemp, ch_gwex_a(ichan)
        if (pco%csvout == "y") then
          write (2507,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ichan, ob(iob)%gis_id, ob(iob)%name,   &
          ch_wat_a(ichan)%area_ha, ch_wat_a(ichan)%precip, ch_wat_a(ichan)%evap, ch_wat_a(ichan)%seep,   &
          ch_stor(ichan), ch_in_a(ichan), ch_out_a(ichan), wtemp, ch_gwex_a(ichan)
        end if
       end if
     end if 
      
      return

100   format (4i6,2i8,2x,a,73e15.4)
101   format (4i6,3i8,2x,a,e15.4)  
       
      end subroutine sd_channel_output