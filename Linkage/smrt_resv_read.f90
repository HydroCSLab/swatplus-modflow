
      !This subroutine reads in cell-lake connection information
  
      !  Prepared by: Ryan Bailey, Colorado State University
      !  August 2024
	
      subroutine smrt_resv_read
      
			use smrt_parm
			use GLOBAL, only:NROW,NCOL,BOTM,IUNIT,DELR,DELC !modflow
			use GWFRESMODULE, only:NRES,NRESOP,IRES,BRES,CRES,BBRES,HRES,IRESCB !modflow reservoir package
      
      implicit none

			character*30 header
			logical :: i_exist
      integer  i,j,k,cell_id,res_id
			integer	 cell_row,cell_col,cell_lay
			real     remainder,res_stage,res_hydc,res_thick,res_depth
			real     cell_area
			
			
      !determine if feature is active
			inquire (file='smrt.resvcells',exist=i_exist)
      if(i_exist) then
			
			write(smrt_log,*) '     smrt.resvcells: being reading'
			mf_res_flag = 1 !set flag
			
			!open file
      open(1280,file='smrt.resvcells')
      
			!read number of drain cells -----------------------------------------------------------------
			read(1280,*) header
			read(1280,*) header
			read(1280,*) mf_num_resv
			
			!turn on the reservoir package, set defaults, and allocate arrays
			IUNIT(17) = 54321
			NRES = mf_num_resv
			NRESOP = 1 !layer
			IRESCB = 0
			allocate(IRES(NCOL,NROW))		!reservoir cell id (from list)
			allocate(HRES(mf_num_resv)) !reservoir stage (m)
			allocate(BRES(NCOL,NROW))   !reservoir bottom elevation (m)
			allocate(BBRES(NCOL,NROW))  !reservoir bottom elevation, including sediments (m)
			allocate(CRES(NCOL,NROW))   !reservoir bottom conductance (m2/day)
			IRES = 0
			HRES = 0.
			BRES = 0.
			BBRES = 0.
			CRES = 0.
			
			!read information for each reservoir cell ---------------------------------------------------
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			allocate(mf_resv_reservoir(mf_num_resv))
			allocate(mf_resv_row(mf_num_resv))
			allocate(mf_resv_col(mf_num_resv))
			mf_resv_reservoir = 0
			!loop through the drain cells
			do i=1,mf_num_resv
			  read(1280,*) cell_id,res_id,res_stage,res_depth,res_hydc,res_thick
			  !determine the row and column of cell
				remainder = mod(cell_id,NCOL)
        if(remainder == 0) then
          cell_row = cell_id / NCOL !row
        else
				  cell_row = int(cell_id/NCOL) + 1 !row
        endif
			  cell_col = cell_id - ((cell_row-1)*NCOL) !column
				cell_area = DELR(cell_col) * DELC(cell_row) !m2
				mf_resv_row(i) = cell_row
				mf_resv_col(i) = cell_col
				!populate RES package arrays
				IRES(cell_col,cell_row) = i
				HRES(i) = res_stage 
				BRES(cell_col,cell_row) = res_stage - res_depth
				BBRES(cell_col,cell_row) = BRES(cell_col,cell_row) - res_thick 
				CRES(cell_col,cell_row) = (res_hydc*DELC(cell_row)*DELR(cell_col)) / res_thick
				!store reservoir to which cell is connected (i.e., gives outflow)
				mf_resv_reservoir(i) = res_id
			enddo !go to next cell
      close(1280)
			
			!allocate for later use
			allocate(mf_resv_rates(mf_num_resv))
			mf_resv_rates = 0.
			
			write(smrt_log,*) '     smrt.resvcells: finished reading'
			endif !check if reservoir connection is active
			
			
      end subroutine smrt_resv_read      
			