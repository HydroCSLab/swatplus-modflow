
      !This subroutine reads in canal connection information
  
      !  Prepared by: Ryan Bailey, Colorado State University
      !  August 2024
	
      subroutine smrt_canl_read
      
			use smrt_parm
      use hydrograph_module
			use GLOBAL, only:NROW,NCOL !modflow
			use GWFRIVMODULE, only:RIVR !modflow river package
      
      implicit none

			character*30 header
			logical :: i_exist
      integer  i,j,k,cell_id,canal_id
			integer  dumid,num_canals,riv_count
			integer	 cell_row,cell_col,cell_lay
			real     remainder
			real     cond,canl_length,canl_stage
			integer, dimension (:), allocatable :: canl_cha,canl_div,canl_beg,canl_end
			real, dimension (:), allocatable :: canl_wid,canl_dep,canl_thk,canl_hdc
			
			
			
			!determine if feature is active
			inquire (file='smrt.canalcells',exist=i_exist)
      if(i_exist) then
			
			write(smrt_log,*) '     smrt.canalcells: begin reading'
			mf_cnl_flag = 1 !set flag
			
      !open file
      open(1280,file='smrt.canalcells')
      
			!read in canal information ---------------------------------------------------------------------
			do i=1,11
			  read(1280,*) header
			enddo
			read(1280,*) num_canals
			allocate(canl_cha(num_canals))
			allocate(canl_div(num_canals))
			allocate(canl_wid(num_canals))
			allocate(canl_dep(num_canals))
			allocate(canl_thk(num_canals))
			allocate(canl_hdc(num_canals))
			allocate(canl_beg(num_canals))
			allocate(canl_end(num_canals))
			read(1280,*)
			do i=1,num_canals
			  read(1280,*) dumid,canl_cha(i),canl_div(i),canl_wid(i),canl_dep(i), &
				                   canl_thk(i),canl_hdc(i),canl_beg(i),canl_end(i)
			enddo
			
			!read in cell-canal connection information --------------------------------------------------
			read(1280,*) header
			read(1280,*) mf_num_canl
			do i=1,6
			  read(1280,*) header
			enddo
			allocate(mf_canl_div(mf_num_canl))
			allocate(mf_canl_beg(mf_num_canl))
			allocate(mf_canl_end(mf_num_canl))
			allocate(mf_canl_cond(mf_num_canl))
			allocate(mf_canl_channel(mf_num_canl))
			mf_canl_div = 0
			mf_canl_beg = 0
			mf_canl_end = 0
			mf_canl_cond = 0.
			mf_canl_channel = 0
			riv_count = sp_ob%modflow
			read(1280,*) header
			do i=1,mf_num_canl
			  read(1280,*) cell_id,cell_lay,canal_id,canl_length,canl_stage
			  !determine the row and column of cell
				remainder = mod(cell_id,NCOL)
        if(remainder == 0) then
          cell_row = cell_id / NCOL !row
        else
				  cell_row = int(cell_id/NCOL) + 1 !row
        endif
			  cell_col = cell_id - ((cell_row-1)*NCOL) !column
			  !calculate conductance
				cond = (canl_wid(canal_id)*canl_length) * canl_hdc(canal_id) / canl_thk(canal_id) !m2/day
				!populate RIVR array
				riv_count = riv_count + 1
				RIVR(1,riv_count) = cell_lay
				RIVR(2,riv_count) = cell_row
				RIVR(3,riv_count) = cell_col
				RIVR(4,riv_count) = canl_stage
				!RIVR(5,riv_count) conductance is turned on/off in smrt_conversion2mf
				RIVR(6,riv_count) = canl_stage - canl_dep(canal_id)
				mf_canl_div(i) = canl_div(canal_id)
				mf_canl_beg(i) = canl_beg(canal_id)
				mf_canl_end(i) = canl_end(canal_id)
				mf_canl_cond(i) = cond
				mf_canl_channel(i) = canl_cha(canal_id)
			enddo !go to next canal cell
      close(1280)
			
			!allocate for later use
			allocate(mf_canl_rates(mf_num_canl))
			mf_canl_rates = 0.
			
			write(smrt_log,*) '     smrt.canalcells: finished reading'
			endif !check if canal connection is active
			
			
			
      end subroutine smrt_canl_read      
			