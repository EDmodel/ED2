!==========================================================================================!
!==========================================================================================!
subroutine ed_masterput_processid(nproc,headnode_num,masterworks,par_run)
!------------------------------------------------------------------------------------------!
!    This routine gives basic processor ID info to the nodes.                              !
!------------------------------------------------------------------------------------------!

  use ed_para_coms, only: mainnum,nmachs,machsize,machnum
  use ed_node_coms, only: mynum,nnodetot,sendnum,recvnum,master_num,machs

  implicit none
  integer :: headnode_num,nproc
  include 'mpif.h'
  integer :: nm
  integer :: ierr
  integer :: par_run
  logical :: masterworks

  mainnum=headnode_num
  master_num=headnode_num
  nmachs=nproc

  if (masterworks) then
    mynum=nmachs+1
    nnodetot=machsize
    sendnum=1
    recvnum=nmachs
  else
    mynum=-1
    nnodetot=nmachs
    sendnum=-1
    recvnum=-1
  end if


  if (par_run == 0) return

  do nm=1,nmachs
     machnum(nm)=nm
     machs(nm)=nm
  enddo

  machs(machsize)=0  !Thats me!!

  do nm=1,nmachs
    call MPI_Send(mainnum,1,MPI_INTEGER,machnum(nm),11,MPI_COMM_WORLD,ierr)
    call MPI_Send(machnum(nm),1,MPI_INTEGER,machnum(nm),12,MPI_COMM_WORLD,ierr)
    call MPI_Send(nm,1,MPI_INTEGER,machnum(nm),13,MPI_COMM_WORLD,ierr)
    call MPI_Send(nmachs,1,MPI_INTEGER,machnum(nm),14,MPI_COMM_WORLD,ierr)
    call MPI_Send(machnum,nmachs,MPI_INTEGER,machnum(nm),15,MPI_COMM_WORLD,ierr)
    call MPI_Send(machsize,1,MPI_INTEGER,machnum(nm),16,MPI_COMM_WORLD,ierr)
  enddo



  return
end subroutine ed_masterput_processid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_masterput_nl(par_run)
!------------------------------------------------------------------------------------------!
!   This subroutine is responsible for sending all the namelist-related information to all !
! the nodes. This is done by the head node, because the head node has read this            !
! information to define the polygon distribution in regional runs (and in the future, the  !
! patches for the SOI runs).                                                               !
!------------------------------------------------------------------------------------------!
   use ed_para_coms,    only: mainnum
   use max_dims,        only: str_len,max_soi,max_ed_regions,nzgmax,n_pft,maxgrds,maxpvars
   use misc_coms,       only: expnme, runtype,itimea,iyeara,imontha,idatea ,itimez,iyearz  &
                             ,imonthz,idatez,dtlsm,radfrq,ifoutput,idoutput,imoutput,iyoutput &
                             ,iclobber,frqfast,sfilin,ffilout,ied_init_mode,ed_inputs_dir   &
                             ,integration_scheme,end_time,current_time,sfilout,frqstate     &
                             ,isoutput,iprintpolys,printvars,pfmtstr,ipmin,ipmax,iedcnfgf   &
                             ,outfast,outstate,out_time_fast,out_time_state,nrec_fast       &
                             ,nrec_state,irec_fast,irec_state,unitfast,unitstate,event_file

   use ed_misc_coms,only: attach_metadata
   use grid_coms,       only: nzg,nzs,ngrids,nnxp,nnyp,deltax,deltay,polelat,polelon       &
                             ,centlat,centlon,time,timmax,nstratx,nstraty
   use soil_coms,       only: isoilflg,nslcon,slz,slmstr,stgoff,veg_database,soil_database &
                             ,soilstate_db,soildepth_db,isoilstateinit,isoildepthflg       &
                             ,isoilbc,runoff_time,zrough,layer_index,nlon_lyr,nlat_lyr
   use met_driver_coms, only: ed_met_driver_db,imettype,metcyc1,metcycf,initial_co2, lapse_scheme
   use mem_sites,       only: n_soi,n_ed_region,grid_type,grid_res,soi_lat,soi_lon         &
                             ,ed_reg_latmin,ed_reg_latmax,ed_reg_lonmin,ed_reg_lonmax      &
                             ,edres,maxpatch,maxcohort
   use physiology_coms, only: istoma_scheme,n_plant_lim
   use phenology_coms , only: iphen_scheme,iphenys1,iphenysf,iphenyf1,iphenyff,phenpath,repro_scheme
   use decomp_coms,     only: n_decomp_lim
   use pft_coms,        only: include_these_pft,pft_1st_check
   use disturb_coms,    only: include_fire,ianth_disturb, treefall_disturbance_rate
   use optimiz_coms,    only: ioptinpt
  
   implicit none
   include 'mpif.h'
   integer :: ierr
   integer :: par_run
   integer :: n
   if (par_run == 0 ) return

!----- First, the namelist-derived type, before I forget... -------------------------------!
   call MPI_Bcast(ngrids,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(timmax,1,MPI_DOUBLE_PRECISION,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(current_time%year,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%month,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%date,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%ifirst,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(end_time%year,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%month,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%date,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%ifirst,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

!----- Now the namelist -------------------------------------------------------------------!
   call MPI_Bcast(expnme,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runtype,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imontha,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatea,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyeara,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimea,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imonthz,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatez,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyearz,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimez,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(dtlsm,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radfrq,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ifoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(imoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoutput,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iclobber,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqfast,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outfast,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outstate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitfast,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitstate,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
                      
   call MPI_Bcast(sfilin,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ffilout,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ied_init_mode,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(sfilout,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqstate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nzg ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nzs ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilflg,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nslcon,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(slz ,nzgmax,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(stgoff,nzgmax,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slmstr,nzgmax,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
 
   call MPI_Bcast(veg_database,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soil_database,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_inputs_dir,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_met_driver_db,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(isoilstateinit,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoildepthflg,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilbc,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(n_soi,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_ed_region,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(grid_type,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(grid_res,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soi_lat,max_soi,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soi_lon,max_soi,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmin,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmax,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmin,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmax,max_ed_regions,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nnxp,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nnyp,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstratx,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstraty,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltax,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltay,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelat,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelon,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlat,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlon,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(integration_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(istoma_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphen_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(repro_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(lapse_scheme,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_plant_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_decomp_lim,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_fire,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ianth_disturb,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_these_pft,n_pft,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pft_1st_check,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(treefall_disturbance_rate,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runoff_time,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iprintpolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   do n=1,maxpvars
      call MPI_Bcast(printvars(n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(pfmtstr(n),str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   end do
   call MPI_Bcast(ipmin,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipmax,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(imettype,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcyc1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcycf,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(initial_co2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iphenys1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenysf,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyf1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyff,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iedcnfgf,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(phenpath,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(event_file,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(maxpatch,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(maxcohort,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

  
   call MPI_Bcast(ioptinpt,str_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(zrough,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(edres,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(attach_metadata,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
!------------------------------------------------------------------------------------------!
!   One last thing to send is the layer index based on the soil_depth. It is not really a  !
! namelist thing, but it is still a setup variable.                                        !
!------------------------------------------------------------------------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Just to wait until the matrix is allocated
   call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_masterput_nl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_masterput_met_header(par_run)
!------------------------------------------------------------------------------------------!
!    This subroutine sends the met driver information to the nodes, which can then read    !
! the hdf5 in parallel                                                                     !
!------------------------------------------------------------------------------------------!
   use ed_para_coms, only: mainnum
   use max_dims, only: max_met_vars,str_len
   use met_driver_coms, only: nformats, met_names, met_nlon,   &
        met_nlat, met_dx, met_dy, met_xmin, met_ymin, met_nv,   &
        met_vars, met_frq, met_interp, ed_met_driver_db, no_ll,  &
        metname_len,metvars_len
   
   implicit none
   include 'mpif.h'
   integer, intent(in) :: par_run
   integer             :: ierr, nsize,f,v

   if (par_run == 0) return
   
   nsize=nformats*max_met_vars

!----- First I send the scalars -----------------------------------------------------------!
   call MPI_Bcast (nformats,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast (no_ll,1,MPI_LOGICAL,mainnum,MPI_COMM_WORLD,ierr)
   
!------------------------------------------------------------------------------------------!
!   Here I need a MPI Barrier. The master has the variables already allocated, but I need  !
! the nodes with their structures already allocated before I proceed sending the info      !
!------------------------------------------------------------------------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   do f=1,nformats
     call MPI_Bcast(met_names(f),metname_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   end do
   
   call MPI_Bcast(met_nlon,nformats,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nlat,nformats,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dx,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dy,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_xmin,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_ymin,nformats,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nv,nformats,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
  
   do f=1,nformats
      do v=1,max_met_vars
         call MPI_Bcast(met_vars(f,v),metvars_len,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
      end do
   end do
   
   call MPI_Bcast(met_frq,nsize,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_interp,nsize,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

!! Just double checking:
!   write(unit=60,fmt='(a,1x,i5)') 'nformats=',nformats
!   write(unit=60,fmt='(a,1x,l1)') 'no_ll=',no_ll
!   do f=1,nformats
!     write (unit=60,fmt='(a,i5,a)')     ' + Format: ',f,':'
!     write (unit=60,fmt='(2a)')         '   - Name    : ',trim(met_names(f))
!     write (unit=60,fmt='(a,i5)')       '   - Nlon    : ',met_nlon(f) 
!     write (unit=60,fmt='(a,i5)')       '   - Nlat    : ',met_nlat(f) 
!     write (unit=60,fmt='(a,es12.5)')   '   - dx      : ',met_dx(f) 
!     write (unit=60,fmt='(a,es12.5)')   '   - dy      : ',met_dy(f) 
!     write (unit=60,fmt='(a,es12.5)')   '   - xmin    : ',met_xmin(f) 
!     write (unit=60,fmt='(a,es12.5)')   '   - ymin    : ',met_ymin(f) 
!     write (unit=60,fmt='(a,i5)')       '   - nv      : ',met_nv(f) 
!     do v=1,met_nv(f)
!        write (unit=60,fmt='(a,i5)')       '   - Variable: ',v
!        write (unit=60,fmt='(2a)')         '     ~ vars    : ',trim(met_vars(f,v))
!        write (unit=60,fmt='(a,es12.5)')   '     ~ frq     : ',met_frq(f,v)
!        write (unit=60,fmt='(a,i5)')       '     ~ interp  : ',met_interp(f,v)
!     end do
!   end do
!!
   return
end subroutine ed_masterput_met_header
!==========================================================================================!
!==========================================================================================!

subroutine ed_masterput_poly_dims(par_run)
   use ed_state_vars , only : gdpy,py_off
   use grid_coms     , only : nnxp,nnyp,ngrids
   use ed_work_vars  , only : work_e
   use ed_para_coms  , only : mainnum,nmachs
   use mem_sites     , only : n_ed_region, n_soi
   implicit none

   include 'mpif.h'
   integer             :: ierr
   integer, intent(in) :: par_run
   integer             :: ifm,nm,npolys,offset
   integer             :: polys_per_node,extras,extra1

   if (par_run == 1 .and. n_soi > 0 ) then
      write (unit=*,fmt='(a)') '--------------------------------------------------------------------'
      write (unit=*,fmt='(a)') '                             FATAL ERROR!                           ' 
      write (unit=*,fmt='(a)') '--------------------------------------------------------------------'
      write (unit=*,fmt='(a)')
      write (unit=*,fmt='(a)') '--------------------------------------------------------------------'
      write (unit=*,fmt='(a)') '                                                                    '
      write (unit=*,fmt='(a)') '  Dear ED user,                                                     '
      write (unit=*,fmt='(a)') '                                                                    '
      write (unit=*,fmt='(a)') '  Thank you for choosing ED, the Ecosystem Demography Model!        '
      write (unit=*,fmt='(a)') 'We know that choosing a model is more than looking for a            '
      write (unit=*,fmt='(a)') 'numeric solver, it is an investment in your academic career.        '
      write (unit=*,fmt='(a)') 'Therefore, it is our commitment to provide our valuable researchers ' 
      write (unit=*,fmt='(a)') 'with the very best in ecosystem dynamics modelling and we want to   '
      write (unit=*,fmt='(a)') 'reaffirm that your satisfaction is our number-one priority.         '
      write (unit=*,fmt='(a)') '  Unfortunately, the option of using parallelism under SOI is not   '
      write (unit=*,fmt='(a)') 'available yet, and we would like to offer our most sincere apologies'
      write (unit=*,fmt='(a)') 'for any inconvenience that this may have caused. Our team is        '
      write (unit=*,fmt='(a)') 'currently working on implementing this capability, and we hope to   '
      write (unit=*,fmt='(a)') 'achieve this goal in a very near future, thus meeting your needs of '
      write (unit=*,fmt='(a)') 'a fast and reliable model. Meanwhile you may find useful to run the '
      write (unit=*,fmt='(a)') 'SOI runs in serial mode instead.                                    '
      write (unit=*,fmt='(a)') '                                                                    '
      write (unit=*,fmt='(a)') '--------------------------------------------------------------------' 
      stop '====> STOP! Create_polygon_map (ed_para_init.f90)'
   end if
  
   ! Now I am just setting the number of polygons that each node will deal with. The polygon ID 
   if (par_run == 0) then
     do ifm=1,n_ed_region
        gdpy(1,ifm)=count(work_e(ifm)%land)
        py_off(1,ifm)=0
     end do
     do ifm=n_ed_region+1,ngrids
        gdpy(1,ifm)=1
        py_off(1,ifm)=0
     end do
   else
      do ifm=1,ngrids

         ! Essentially, the domain decomposition is 
         ! determined in the next sequence of lines
         ! ----------------------------------------

         polys_per_node = floor(real(count(work_e(ifm)%land))/real(nmachs+1))
         extras = count(work_e(ifm)%land)-polys_per_node*real(nmachs+1)

         if (extras>0) then
            extra1=1
         else
            extra1=0
         endif
         ! ----------------------------------------
         
         npolys=0
         offset=0
         do nm=1,nmachs 
            offset=offset+npolys
            npolys = polys_per_node + extra1

            gdpy(nm,ifm)   = npolys
            py_off(nm,ifm) = offset
            call MPI_Bcast(npolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(offset,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

            extras = extras - 1
            if (extras>0) then
               extra1=1
            else
               extra1=0
            endif
            
           
         end do
         nm=nmachs+1
         offset=offset+npolys
         npolys = polys_per_node + extra1

         gdpy(nm,ifm)   = npolys
         py_off(nm,ifm) = offset

      end do
   end if

   return
end subroutine ed_masterput_poly_dims
!==========================================================================================!
!==========================================================================================!

subroutine ed_masterput_worklist_info(par_run)

  use max_dims, only: maxmach
  use grid_coms, only: ngrids
  use ed_work_vars, only: work_e,work_vars,ed_alloc_work,ed_nullify_work,ed_dealloc_work
  use ed_para_coms, only: nmachs,mainnum,machnum
  use soil_coms, only: isoilflg
  use ed_node_coms, only: mxp,myp,i0,j0
  use ed_state_vars, only: gdpy,py_off

  implicit none
  include 'mpif.h'
  integer, intent(in) :: par_run
  integer :: npoly
  integer :: offset
  integer :: nm
  integer :: ifm
  integer :: mpiid
  integer :: ierr

  type(work_vars), allocatable, dimension(:)   :: sc_work

  real,            allocatable, dimension(:) :: rscratch
  integer,         allocatable, dimension(:) :: iscratch

  

  if (par_run == 1) then
     
     do nm=1,nmachs
        do ifm=1,ngrids
           
           npoly  = gdpy(nm,ifm)
           offset = py_off(nm,ifm)
           
           allocate(rscratch(npoly),iscratch(npoly))

           rscratch(1:npoly) = work_e(ifm)%vec_glon(offset+1:offset+npoly)
           mpiid=0000+maxmach*(ifm-1)+nm
           call MPI_Send(rscratch,npoly,MPI_REAL,machnum(nm),mpiid,MPI_COMM_WORLD,ierr)
  
           rscratch(1:npoly) = work_e(ifm)%vec_glat(offset+1:offset+npoly)
           mpiid=10000+maxmach*(ifm-1)+nm
           call MPI_Send(rscratch,npoly,MPI_REAL,machnum(nm),mpiid,MPI_COMM_WORLD,ierr)
           
           rscratch(1:npoly) = work_e(ifm)%vec_landfrac(offset+1:offset+npoly)
           mpiid=20000+maxmach*(ifm-1)+nm
           call MPI_Send(rscratch,npoly,MPI_REAL,machnum(nm),mpiid,MPI_COMM_WORLD,ierr)
           
           iscratch(1:npoly) = work_e(ifm)%vec_ntext(offset+1:offset+npoly)
           mpiid=30000+maxmach*(ifm-1)+nm
           call MPI_Send(iscratch,npoly,MPI_INTEGER,machnum(nm),mpiid,MPI_COMM_WORLD,ierr)

           deallocate(rscratch,iscratch)

        enddo
     enddo
     
  endif
  
  ! Deallocate each of the global vectors, and reallocate (and fill) as local vectors

  nm=nmachs+1
  allocate(sc_work(ngrids))

  do ifm=1,ngrids
     
     npoly  = gdpy(nm,ifm)
     offset = py_off(nm,ifm)

     allocate(sc_work(ifm)%vec_glon(npoly))
     allocate(sc_work(ifm)%vec_glat(npoly))
     allocate(sc_work(ifm)%vec_landfrac(npoly))
     allocate(sc_work(ifm)%vec_ntext(npoly))

     sc_work(ifm)%vec_glon(1:npoly) = work_e(ifm)%vec_glon(offset+1:offset+npoly)
     sc_work(ifm)%vec_glat(1:npoly) = work_e(ifm)%vec_glat(offset+1:offset+npoly)
     sc_work(ifm)%vec_landfrac(1:npoly) = work_e(ifm)%vec_landfrac(offset+1:offset+npoly)
     sc_work(ifm)%vec_ntext(1:npoly) = work_e(ifm)%vec_ntext(offset+1:offset+npoly)

     call ed_dealloc_work(work_e(ifm))
     

  enddo

  deallocate(work_e)

  ! So here it will be 2 (node-style) if it is a parallel run, and 0 if it is a serial run
  call ed_mem_alloc(2*par_run) 

  do ifm=1,ngrids

     npoly  = gdpy(nm,ifm)
     
     ! Longitudes
     allocate(work_e(ifm)%vec_glon(npoly))
     work_e(ifm)%vec_glon(1:npoly) = sc_work(ifm)%vec_glon(1:npoly)
     
     ! Latitudes
     allocate(work_e(ifm)%vec_glat(npoly))
     work_e(ifm)%vec_glat(1:npoly) = sc_work(ifm)%vec_glat(1:npoly)
     
     ! Land fractions
     allocate(work_e(ifm)%vec_landfrac(npoly))
     work_e(ifm)%vec_landfrac(1:npoly) = sc_work(ifm)%vec_landfrac(1:npoly)
     
     ! Soil Textures
     allocate(work_e(ifm)%vec_ntext(npoly))
     work_e(ifm)%vec_ntext(1:npoly) = sc_work(ifm)%vec_ntext(1:npoly)

     ! Need we say more
     deallocate(sc_work(ifm)%vec_glon)
     deallocate(sc_work(ifm)%vec_glat)
     deallocate(sc_work(ifm)%vec_landfrac)
     deallocate(sc_work(ifm)%vec_ntext)

   end do
   
   deallocate(sc_work)


  
  return
end subroutine ed_masterput_worklist_info


!==========================================================================================!
!==========================================================================================!


subroutine ed_nodeget_processid(init)

  use max_dims
  use ed_node_coms

  implicit none
  integer :: init

  include 'mpif.h'
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr

  if(init == 1) then
     
     call MPI_Recv(master_num,1,MPI_INTEGER,0,11,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(mchnum,1,MPI_INTEGER,0,12,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(mynum,1,MPI_INTEGER,0,13,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(nmachs,1,MPI_INTEGER,0,14,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(machs,nmachs,MPI_INTEGER,0,15,MPI_COMM_WORLD,status,ierr)
     call MPI_Recv(nnodetot,1,MPI_INTEGER,0,16,MPI_COMM_WORLD,status,ierr)
     
     recvnum = mynum-1
     sendnum = mynum+1
     if (mynum == nmachs) sendnum=0
  endif
  write(unit=*,fmt='(a,1x,i5,1x,a)') '---> Node',mynum,'got first message!'

  return
end subroutine ed_nodeget_processid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_nl

  use ed_node_coms, only: master_num,mynum
!------------------------------------------------------------------------------------------!
!   This subroutine is responsible for getting all the namelist-related information in     !
! every node.                                                                              !
!------------------------------------------------------------------------------------------!
   use max_dims,        only: str_len,max_soi,max_ed_regions,nzgmax,n_pft,maxgrds,maxpvars
   use misc_coms,       only: expnme, runtype,itimea,iyeara,imontha,idatea ,itimez,iyearz  &
                             ,imonthz,idatez,dtlsm,radfrq,ifoutput,idoutput,imoutput, iyoutput &
                             ,iclobber,frqfast,sfilin,ffilout,ied_init_mode,ed_inputs_dir   &
                             ,integration_scheme,end_time,current_time,isoutput,sfilout    &
                             ,frqstate,iprintpolys,printvars,pfmtstr,ipmin,ipmax,iedcnfgf  &
                             ,outfast,outstate,out_time_fast,out_time_state,nrec_fast      &
                             ,nrec_state,irec_fast,irec_state,unitfast,unitstate,event_file

   use grid_coms,       only: nzg,nzs,ngrids,nnxp,nnyp,deltax,deltay,polelat,polelon       &
                             ,centlat,centlon,time,timmax,nstratx,nstraty
   use soil_coms,       only: isoilflg,nslcon,slz,slmstr,stgoff,veg_database,soil_database &
                             ,soilstate_db,soildepth_db,isoilstateinit,isoildepthflg       &
                             ,isoilbc,runoff_time,zrough,layer_index,nlon_lyr,nlat_lyr
   use met_driver_coms, only: ed_met_driver_db,imettype,metcyc1,metcycf,initial_co2,lapse_scheme
   use mem_sites,       only: n_soi,n_ed_region,grid_type,grid_res,soi_lat,soi_lon         &
                             ,ed_reg_latmin,ed_reg_latmax,ed_reg_lonmin,ed_reg_lonmax      &
                             ,edres,maxpatch,maxcohort
   use physiology_coms, only: istoma_scheme,n_plant_lim
   use phenology_coms , only: iphen_scheme,iphenys1,iphenysf,iphenyf1,iphenyff,phenpath,repro_scheme
   use decomp_coms,     only: n_decomp_lim
   use disturb_coms,    only: include_fire,ianth_disturb, treefall_disturbance_rate
   use optimiz_coms,    only: ioptinpt
   use ed_misc_coms,only: attach_metadata
   use pft_coms,        only: include_these_pft,pft_1st_check

   implicit none
   include 'mpif.h'
   integer :: ierr
   integer :: n

!----- First, the namelist-derived type, before I forget... -------------------------------!
   call MPI_Bcast(ngrids,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(timmax,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(current_time%year,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%month,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%date,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(current_time%ifirst,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(end_time%year,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%month,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%date,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(end_time%ifirst,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

!----- Now the namelist -------------------------------------------------------------------!
   call MPI_Bcast(expnme,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runtype,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imontha,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatea,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyeara,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimea,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(imonthz,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idatez,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyearz,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itimez,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(dtlsm,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(radfrq,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ifoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(idoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(imoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iyoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoutput,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iclobber,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqfast,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outfast,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(outstate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitfast,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(unitstate,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
                      
   call MPI_Bcast(sfilin,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ffilout,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ied_init_mode,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(sfilout,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(frqstate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nzg ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nzs ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilflg,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nslcon,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(slz ,nzgmax,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(stgoff,nzgmax,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(slmstr,nzgmax,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
 
   call MPI_Bcast(veg_database,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soil_database,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_inputs_dir,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_met_driver_db,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(isoilstateinit,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoildepthflg,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isoilbc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(n_soi,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_ed_region,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(grid_type,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(grid_res,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soi_lat,max_soi,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(soi_lon,max_soi,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmin,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_latmax,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmin,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ed_reg_lonmax,max_ed_regions,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(nnxp,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nnyp,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstratx,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nstraty,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltax,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltay,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelat,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(polelon,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlat,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(centlon,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(integration_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(istoma_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphen_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(repro_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(lapse_scheme,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_plant_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(n_decomp_lim,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_fire,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ianth_disturb,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(include_these_pft,n_pft,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pft_1st_check,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(treefall_disturbance_rate,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(runoff_time,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iprintpolys,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   do n=1,maxpvars
      call MPI_Bcast(printvars(n),str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(pfmtstr(n),str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   end do
   call MPI_Bcast(ipmin,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipmax,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(imettype,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcyc1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(metcycf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(initial_co2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iphenys1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenysf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyf1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iphenyff,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(iedcnfgf,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(phenpath,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(event_file,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(maxpatch,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(maxcohort,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

  
   call MPI_Bcast(ioptinpt,str_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(zrough,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(edres,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(attach_metadata,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

!------------------------------------------------------------------------------------------!
!     Receiving the layer index based on soil_depth. This is allocatable, so I first       !
! allocate, then let the master know that it is safe to send to me and I reveive the data. !
!------------------------------------------------------------------------------------------!
   if (allocated(layer_index)) deallocate(layer_index)
   allocate(layer_index(nlat_lyr,nlon_lyr))
   call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Safe to receive the data.
   call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_nodeget_nl
!=========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_met_header()
!------------------------------------------------------------------------------------------!
!    This subroutine sends the met driver information to the nodes, which can then read    !
! the hdf5 in parallel                                                                     !
!------------------------------------------------------------------------------------------!
   use ed_node_coms, only: master_num,mynum
   use max_dims, only: max_met_vars,str_len
   use met_driver_coms, only: nformats, met_names, met_nlon,   &
        met_nlat, met_dx, met_dy, met_xmin, met_ymin, met_nv,   &
        met_vars, met_frq, met_interp, ed_met_driver_db, no_ll,  &
        metname_len,metvars_len

   implicit none
   include 'mpif.h'
   integer             :: ierr, nsize,f,v

   

!----- First I get the scalars ------------------------------------------------------------!
   call MPI_Bcast (nformats,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast (no_ll,1,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)

   nsize=nformats*max_met_vars

!----- Allocate the vectors and matrices --------------------------------------------------!
   allocate(met_names(nformats))
   allocate(met_nlon(nformats))
   allocate(met_nlat(nformats))
   allocate(met_dx(nformats))
   allocate(met_dy(nformats))
   allocate(met_xmin(nformats))
   allocate(met_ymin(nformats))
   allocate(met_nv(nformats))
   allocate(met_vars(nformats, max_met_vars))
   allocate(met_frq(nformats, max_met_vars))
   allocate(met_interp(nformats, max_met_vars))
       
!------------------------------------------------------------------------------------------!
!   Here I need a MPI Barrier. I don't want the master sending information before the      !
! variables are allocated in this node.                                                    !
!------------------------------------------------------------------------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   do f=1,nformats
     call MPI_Bcast(met_names(f),metname_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   end do
   
   call MPI_Bcast(met_nlon,nformats,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nlat,nformats,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dx,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_dy,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_xmin,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_ymin,nformats,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_nv,nformats,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
   do f=1,nformats
      do v=1,max_met_vars
         call MPI_Bcast(met_vars(f,v),metvars_len,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
      end do
   end do
   
   call MPI_Bcast(met_frq,nsize,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(met_interp,nsize,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine ed_nodeget_met_header
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_poly_dims
   use ed_state_vars, only: gdpy,py_off
   use ed_node_coms, only: master_num,nmachs,mynum
   use grid_coms, only: ngrids
   implicit none
   include 'mpif.h'
   integer :: ierr
   integer :: ifm,nm
  
   do ifm=1,ngrids
      do nm=1,nmachs
         call MPI_Bcast(gdpy(nm,ifm),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(py_off(nm,ifm),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
      end do
   end do
   return
end subroutine ed_nodeget_poly_dims
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
subroutine ed_nodeget_worklist_info

  use max_dims, only: maxmach
  use grid_coms,  only: ngrids
  use ed_work_vars,  only: work_e
  use ed_node_coms,  only: mynum,nmachs,master_num
  use ed_state_vars, only: gdpy

  implicit none
  include 'mpif.h'
  
  integer,dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr
  integer :: npolygons
  integer :: mpiid
  integer :: ifm
  
  do ifm=1,ngrids
     
     npolygons = gdpy(mynum,ifm)
     
     ! Allocate the work vectors
     
     allocate(work_e(ifm)%vec_glon(npolygons))
     allocate(work_e(ifm)%vec_glat(npolygons))
     allocate(work_e(ifm)%vec_landfrac(npolygons))
     allocate(work_e(ifm)%vec_ntext(npolygons))

     mpiid=0000+maxmach*(ifm-1)+mynum
     call MPI_Recv(work_e(ifm)%vec_glon(1:npolygons),npolygons, &
          MPI_REAL,0,mpiid,MPI_COMM_WORLD,status,ierr)
     
     mpiid=10000+maxmach*(ifm-1)+mynum
     call MPI_Recv(work_e(ifm)%vec_glat(1:npolygons),npolygons, &
          MPI_REAL,0,mpiid,MPI_COMM_WORLD,status,ierr)
     
     mpiid=20000+maxmach*(ifm-1)+mynum
     call MPI_Recv(work_e(ifm)%vec_landfrac(1:npolygons),npolygons, &
          MPI_REAL,0,mpiid,MPI_COMM_WORLD,status,ierr)
     
     mpiid=30000+maxmach*(ifm-1)+mynum
     call MPI_Recv(work_e(ifm)%vec_ntext(1:npolygons),npolygons, &
          MPI_INTEGER,0,mpiid,MPI_COMM_WORLD,status,ierr)

     
  end do

  return
end subroutine ed_nodeget_worklist_info
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
subroutine mk_2_buff(a,b,n1,n2,m1,m2,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,m1,m2,i1,i2,j1,j2
  real :: a(n1,n2),b(m1,m2)
     
     b(1:m1,1:m2)=a(i1:i2,j1:j2)

  return
end subroutine mk_2_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mk_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)
     
     b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

  return
end subroutine mk_2p_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mk_3_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)
     
     b(1:m1,1:m2,1:m3)=a(1:n1,i1:i2,j1:j2)

  return
end subroutine mk_3_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mk_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2
  real :: a(n1,n2,n3,m4),b(m1,m2,m3,m4)
     
     b(1:m1,1:m2,1:m3,1:m4)=a(1:n1,i1:i2,j1:j2,1:n4)

  return
end subroutine mk_4_buff
!==========================================================================================!
!==========================================================================================!
