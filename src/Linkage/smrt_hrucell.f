      subroutine smrt_hrucell

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine reads in the geographic connection information between HRUs and MODFLOW cells

!     Import variables
	use hydrograph_module, only:sp_ob
      use GLOBAL, only: ncol,nrow,nlay,delr,delc !modflow
      use smrt_parm !modflow linkage
	
      implicit none
      
      !Initialize local variables
	character*30 :: header
	integer :: i,k,n                !          |
      integer :: hru_id               !          |
      integer :: nhru_connected       !          |number of HRUs connected to grid cells
	integer :: num_hru							!          |
	integer :: cell_count           !          |
	integer :: cell_row
	integer :: cell_col
	integer :: cell_id    
	integer :: hru_read
      real ::    poly_area            !m2        |
	real ::    hru_area             !m2        |surface area of the hru
	real ::    cell_area            !m2        |surface area of the cell
	real ::    remainder
      real ::    delay
	
	
	write(smrt_log,*) '     smrt.hrucells: begin reading'
	
	!retrieve number of HRUs in the model
	num_hru = sp_ob%hru
	
	
	!open file
      open(6004,file='smrt.hrucells')
      read(6004,*) header
	
	
	!read in groundwater delay term -----------------------------------------------------------------------
      allocate(gw_delaye(sp_ob%hru))
      allocate(gw_rech(sp_ob%hru))
      read(6004,*) header
      read(6004,*) hru_read
      if(hru_read.eq.0) then
        read(6004,*) delay
        gw_delaye = Exp(-1./(delay + 1.e-6))
      else
        do n=1,sp_ob%hru
          read(6004,*) delay
          gw_delaye(n) = Exp(-1./(delay + 1.e-6))
        enddo
      endif
      gw_rech = 0.
	
	
	!read in flag for simulating groundater-->soil transfer -----------------------------------------------
	read(6004,*) header
	read(6004,*) mf_soil_flag
	
	
      !read in list of HRUs that are spatially connected to grid cells --------------------------------------
      read(6004,*) header
	read(6004,*) nhru_connected !number of HRUs spatially connected to grid cells
      allocate(mf_hrus_connected(num_hru))
      mf_hrus_connected = 0
      do i=1,nhru_connected
        read(6004,*) hru_id
        mf_hrus_connected(hru_id) = 1
      enddo
      
	!read in the HRU-cell connection information ----------------------------------------------------------
      read(6004,*) header
      read(6004,*) header
      allocate(mf_hru_num_cells(num_hru))
      allocate(mf_hru_cell_row(num_hru,5000))
	allocate(mf_hru_cell_col(num_hru,5000))
      allocate(mf_hru_cell_fract(num_hru,5000))
      allocate(mf_cell_fract(num_hru,5000))
      mf_hru_num_cells = 0
      mf_hru_cell_fract = 0.
      do k=1,num_hru !loop through the HRUs: read in the connected cells for each HRU
        if(mf_hrus_connected(k) == 1) then
        hru_id = k
        cell_count = 0
        do while (hru_id == k)
          !read in the information for the cell connected to the HRU
          read(6004,*) hru_id,hru_area,cell_id,poly_area
			cell_count = cell_count + 1
			!determine the row and column of the cell
			remainder = mod(cell_id,ncol)
          if(remainder == 0) then
            cell_row = cell_id / ncol !row
          else
				cell_row = int(cell_id/ncol) + 1 !row
          endif
			cell_col = cell_id - ((cell_row-1)*ncol) !column
			mf_hru_cell_row(k,cell_count) = cell_row
			mf_hru_cell_col(k,cell_count) = cell_col
			!calculate the fraction of HRU area that is occupied by the cell
			mf_hru_cell_fract(k,cell_count) = poly_area / hru_area
			!calculate the fraction of the cell area that is occupied by the HRU
			cell_area = delr(cell_col) * delc(cell_row) 
			mf_cell_fract(k,cell_count) = poly_area / cell_area !fraction of cell area
          read(6004,*,end=10) hru_id
          backspace(6004)
        enddo
10      mf_hru_num_cells(k) = cell_count
        endif
      enddo !go to the next HRU
	close(6004)
      
	
	write(smrt_log,*) '     smrt.hrucells: finished reading'
	
      return
	
      end subroutine smrt_hrucell
	
	
	