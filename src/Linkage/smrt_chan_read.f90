
      !This subroutine reads in channel-cell connection information and writes the modflow.con input file
  
      !  Prepared by: Ryan Bailey, Colorado State University
      !  August 2024
	
      subroutine smrt_chan_read
      
			use smrt_parm
      use hydrograph_module
			use GLOBAL, only:NROW,NCOL !modflow
			use GWFRIVMODULE, only:RIVR !modflow
      
      implicit none

			character*30 header
      character(len=8) :: col_head_con(17)
      integer  i,j,k,in_chan,in_con,cell_id,chan_id,chan_zone,dum1,dum2,chan_cell
			integer	 num_chancells
			integer  nzone,dumid
			integer	 cell_row,cell_col,cell_lay
      real     bed_elev,chan_length,dum3
			real     remainder
			real, dimension (:), allocatable :: zone_hydc,zone_thck
			
			
			write(smrt_log,*)
			write(smrt_log,*) 'smrt.chancells: begin reading'
			write(smrt_log,*) 'smrt.chancells: create modflow.con'
			
      !integers for input and output files
      open(1280,file='smrt.chancells')
      open(1281,file='modflow.con')
      
			!prepare headers for modflow.con file ----------------------------------------------------------
      write(1281,*) 'modflow.con: channel-cell spatial connections'
      col_head_con = [character(len=17) :: "NUMB","NAME","GISID","AREA","LAT","LONG","ELEV","CELL","WST","CONST","OVER",  &
            "RULE","SRC_TOT ","OBTYPE_OUT1 ","OBTYPNO_OUT1 ","HTYPE_OUT1 ","FRAC_OUT1 "]
      write(1281,103) (col_head_con(j),j=1,17)
      dum1 = 1
      dum2 = 0
      dum3 = 1.00
			
			!read in the channel bed properties, by zone ---------------------------------------------------
			read(1280,*) header
			read(1280,*) header
			read(1280,*) nzone
			allocate(zone_hydc(nzone))
			allocate(zone_thck(nzone))
			read(1280,*) header
			do i=1,nzone
			  read(1280,*) dumid,zone_hydc(i),zone_thck(i)
			enddo
			
			!read cell-channel connection information ------------------------------------------------------
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			read(1280,*) header
			allocate(mf_chan_channel(sp_ob%modflow))
			allocate(mf_rivr_lngt(sp_ob%modflow))
			allocate(mf_rivr_hydc(sp_ob%modflow))
			allocate(mf_rivr_thck(sp_ob%modflow))
			mf_chan_channel = 0
			mf_rivr_lngt = 0.
			mf_rivr_hydc = 0.
			mf_rivr_thck = 0.
			!loop through the cells connected to channels
			do i=1,sp_ob%modflow
			  read(1280,*) cell_id,cell_lay,bed_elev,chan_id,chan_length,chan_zone
				mf_chan_channel(i) = chan_id
			  !determine the row and column of cell
				remainder = mod(cell_id,NCOL)
        if(remainder == 0) then
          cell_row = cell_id / NCOL !row
        else
				  cell_row = int(cell_id/NCOL) + 1 !row
        endif
			  cell_col = cell_id - ((cell_row-1)*NCOL) !column
				!populate RIVR array
				RIVR(1,i) = cell_lay
				RIVR(2,i) = cell_row
				RIVR(3,i) = cell_col
				!RIVR(4,i): channel stage is populated in smrt_conversion2mf
				!RIVR(5,i): channel bed conductance is populated in smrt_conversion2mf		
				RIVR(6,i) = bed_elev
				!store length, K, and thickness for conductance calculations in smrt_conversion2mf
				mf_rivr_lngt(i) = chan_length
				mf_rivr_hydc(i) = zone_hydc(chan_zone)
				mf_rivr_thck(i) = zone_thck(chan_zone)
				!write channel connection to modflow.con
				write(1281,107) dum1,dum1,dum1,dum2,dum2,dum2,dum2,cell_id,dum1,dum2,dum2,dum2,dum1,'sdc',chan_id,'tot',dum3
			enddo !go to next cell
      close(1280)
			close(1281)
			
			!allocate for later use
			allocate(mf_chan_rates(sp_ob%modflow))
			mf_chan_rates = 0.
			
			write(smrt_log,*) 'smrt.chancells: finished reading'
			write(smrt_log,*)
			write(smrt_log,*)
			
      
 103  format(30(a10))   
 107  format(i5,i5,i6,i5,i4,i5,i5,i10,i4,i6,i5,i5,i8,a12,i13,a11,f6.2)
 
      end subroutine smrt_chan_read      
			