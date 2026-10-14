
      !This subroutine reads in drainage connection information
  
      !  Prepared by: Ryan Bailey, Colorado State University
      !  August 2024
	
      subroutine smrt_drng_read
      
			use smrt_parm
      use hydrograph_module
			use GLOBAL, only:NROW,NCOL,BOTM !modflow
			use GWFDRNMODULE !modflow drain package
      
      implicit none

			character*30 header
			logical :: i_exist
      integer  i,j,k,cell_id,cond_zone
			integer	 drain_chan
			integer  nzone,dumid
			integer	 cell_row,cell_col,cell_lay
      real     drain_depth
			real     remainder
			real, dimension (:), allocatable :: zone_cond
			
	
			
			!determine if feature is active
			inquire (file='smrt.drngcells',exist=i_exist)
      if(i_exist) then
			
			write(smrt_log,*) '     smrt.drngcells: begin reading'
			mf_drn_flag = 1 !set flag
			
			
      !open file
      open(1280,file='smrt.drngcells')
      
			!read in drain conductance zones ---------------------------------------------------------------
			read(1280,*) header
			read(1280,*) header
			read(1280,*) nzone
			allocate(zone_cond(nzone))
			read(1280,*) header
			do i=1,nzone
			  read(1280,*) dumid,zone_cond(i)
			enddo
			
			!read number of drain cells --------------------------------------------------------------------
			read(1280,*) header
			read(1280,*) mf_num_drng
			
			!read cell-channel connection information ------------------------------------------------------
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			allocate(mf_drng_channel(mf_num_drng))
			mf_drng_channel = 0
			!loop through the drain cells
			do i=1,mf_num_drng
			  read(1280,*) cell_id,cell_lay,drain_depth,drain_chan,cond_zone
			  !determine the row and column of cell
				remainder = mod(cell_id,NCOL)
        if(remainder == 0) then
          cell_row = cell_id / NCOL !row
        else
				  cell_row = int(cell_id/NCOL) + 1 !row
        endif
			  cell_col = cell_id - ((cell_row-1)*NCOL) !column
				!populate DRAI array
				DRAI(1,i) = cell_lay
				DRAI(2,i) = cell_row
				DRAI(3,i) = cell_col
				DRAI(4,i) = BOTM(cell_col,cell_row,0) - drain_depth !drain elevation (m)
				DRAI(5,i) = zone_cond(cond_zone) !conductance (m2/day)
				!store channel to which drainage cell is connected (i.e., gives outflow)
				mf_drng_channel(i) = drain_chan
			enddo !go to next cell
      close(1280)
			
			!set the modflow number (used in modflow subroutines)
			NDRAIN = mf_num_drng
			
			!allocate for later use
			allocate(mf_drng_rates(mf_num_drng))
			mf_drng_rates = 0.
			
			write(smrt_log,*) '     smrt.drngcells: finished reading'
			endif !check if drainage connection is active
			
			
			
      end subroutine smrt_drng_read      
			