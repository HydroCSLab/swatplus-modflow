      subroutine smrt_conversion2rt3d

!!    This subroutine converts the necessary SWAT variables (per hru) into RT3D variables
      
      !Import variables
      use GLOBAL, only:LENUNI,ITMUNI,NCOL,NROW !MODFLOW
      use GWFRIVMODULE, only:RIVR,MXRIVR !MODFLOW
	use GWFRCHMODULE, only: RECH
      use mf_rt_link, only: rt_criv,crch !MODFLOW-RT3D linkage
      use hydrograph_module
      use channel_module
      use smrt_parm !SWAT-MODFLOW
      
      implicit none
      
      !Define local variables
      integer mf_lengthUnit,mf_timeUnit,nriver_cells
      integer i,j,k
	integer	cell_col,cell_row,cha_num,ob_num,ob_cha_num,num_channels
      integer mf_gridID,riv_cellID,chan_id,dum
      real river_no3,river_p,cha_no3,cha_p,fraction,chan_volume
      real :: rech_conc             !kg/ha       |mass of solute in recharge water
      real :: hru_mass_no3					!g           |nitrate mass in leaching water (hru total)
	real :: hru_mass_p            !g           |p mass in leaching water (hru total)
	real :: cell_mass_no3         !g           |nitrate mass in leaching water to cell
	real :: cell_mass_p           !g					 |p mass in leaching water to cell
	
	
	
	   
      !leaching mass (kg/ha) in recharge water --------------------------------------------------------------
	do k=1,sp_ob%hru
	  !no3-n
	  rech_conc = mf_rechsol(k,1)
	  mf_rechsol(k,1) = ((1.-gw_delaye(k))*mf_percsol(k,1)) + 
     &                    (gw_delaye(k)*rech_conc)
	  !p
        rech_conc = mf_rechsol(k,2)
	  mf_rechsol(k,2) = ((1.-gw_delaye(k))*mf_percsol(k,2)) + 
     &                    (gw_delaye(k)*rech_conc)
	enddo
	
	
	!convert leaching mass (kg/ha) to cell concentration for RT3D arrays ----------------------------------
	open(28944,file='smrt_out_cell_leach')
	crch = 0.
	rt_rech_load = 0.
	ob_num = sp_ob1%hru  !object number of first HRU
      do k=1,sp_ob%hru
        !leaching mass (g)
	  hru_mass_no3 = mf_rechsol(k,1) * ob(ob_num)%area_ha * 1000. !g 
	  hru_mass_p = mf_rechsol(k,2) * ob(ob_num)%area_ha * 1000. !g 
        do i=1,mf_hru_num_cells(k)
		  cell_row = mf_hru_cell_row(k,i)
			cell_col = mf_hru_cell_col(k,i)
			!no3-n
			cell_mass_no3 = hru_mass_no3 * mf_hru_cell_fract(k,i) !g
			if(RECH_down(cell_col,cell_row) > 0) then !only divided if there is recharge
			  crch(cell_col,cell_row,1) = cell_mass_no3 / 
     &                                  RECH_down(cell_col,cell_row)
          else
			  crch(cell_col,cell_row,1) = 0.
			endif
			rt_rech_load(cell_col,cell_row,1) = cell_mass_no3 / 1000. !kg
			!p
			cell_mass_p = hru_mass_p * mf_hru_cell_fract(k,i) !g
			if(RECH_down(cell_col,cell_row) > 0) then !only divided if there is recharge
			  crch(cell_col,cell_row,2) = cell_mass_p / 
     &                                  RECH_down(cell_col,cell_row)
          else
			  crch(cell_col,cell_row,2) = 0.
			endif
			rt_rech_load(cell_col,cell_row,2) = cell_mass_p / 1000. !kg
			write(28944,*) cell_row,cell_col,
     &                   rt_rech_load(cell_col,cell_row,1),
     &                   rt_rech_load(cell_col,cell_row,2)
		enddo
        ob_num = ob_num + 1
	enddo
	  
	
	!calculate concentration (g/m3) in RT3D river cells ---------------------------------------------------
	ob_num = sp_ob1%modflow !object number of first River Cell
	do i=1,sp_ob%modflow
        !object and channel number  
        ob_cha_num = ob(ob_num)%obj_out(1)   
	  chan_id = ob(ob_cha_num)%num
		chan_volume = ch_stor(chan_id)%flo !m3
		if(chan_volume > 0) then
		  rt_criv(i,1) = (ch_stor(chan_id)%no3 * 1000.) / chan_volume !g/m3 no3-n
		  rt_criv(i,2) = (ch_stor(chan_id)%solp * 1000.) / chan_volume !g/m3 p
		else
		  rt_criv(i,1) = 0.
			rt_criv(i,2) = 0.
		endif
	  ob_num = ob_num + 1 !increment the object number
	enddo
	
	!zero our denitrification mass array ------------------------------------------------------------------
	rt_rctn_load = 0.
	
		
      return
      end