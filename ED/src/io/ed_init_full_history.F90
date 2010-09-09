!==========================================================================================!
!==========================================================================================!
subroutine init_full_history_restart()


  use ed_max_dims, only: n_pft
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, phenology
  use ed_misc_coms, only: sfilin, ied_init_mode,current_time
  use mem_polygons, only: grid_res,edres
  use consts_coms, only: pio180
  use ed_misc_coms, only: use_target_year, restart_target_year,ied_init_mode,runtype
  use ed_state_vars,only: polygontype,sitetype,patchtype,edtype, &
       edgrid_g,allocate_sitetype,allocate_patchtype,allocate_polygontype
  use soil_coms, only: alloc_soilgrid
  use grid_coms,only:ngrids
  use ed_node_coms, only: mynum,nmachs,nnodetot,mchnum,machs
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs
  use allometry, only: dbh2h

  implicit none
  
  character(len=1)  :: vnam
  character(len=3)  :: cgr
  character(len=128) :: hnamel
  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  logical :: exists ! File existence
  real,allocatable :: file_lats(:),file_lons(:)
  integer,allocatable :: pysi_n(:),pysi_id(:)
  integer,allocatable :: sipa_n(:),sipa_id(:)
  integer,allocatable :: paco_n(:),paco_id(:)
  
  real :: ll_tolerance
  real :: minrad, currad

  integer :: ngr,ifpy,ipft
  integer :: ipy,isi,ipa,ico
  integer :: py_index,si_index,pa_index

  ! HDF5 types are defined here
  integer :: hdferr
  include 'mpif.h'
  real(kind=8) :: dbletime

 

  ! ------------------------------------------------------------------------------
  ! There are two types of history restarts that can be done.  An exact restart
  ! Assumes that you are continuing with the exact same configuration as a given
  ! model simulation that wrote the file in which you are using.  In this
  ! case, each node starts with a list of polygons, and searches the HDF history
  ! file for these polygons.  It expects to find at least one polygon in the file 
  ! within 250 meters of proximity.  If it does not find this it stops. It then fills
  ! each of these polygons with data, and traverses the data hierarchical tree
  ! that roots from each polygon, initializing the model to the exact same state
  ! as the end of the previous run.  A search based restart does not expect to 
  ! find exact matches, and may not use all of the polygons in the file. But none
  ! the less, it will traverse the tree from these polygons and populate the model 
  ! states with what is found in the files tree.
  ! -------------------------------------------------------------------------------
  ! Currently, we are not doing collective reads from the dataset, this may not
  ! be a feasible option, because at least during the polygon read-in, we are not
  ! sure which chunks to read in, so we take the whole vector of polygons for each
  ! node. Note that all nodes read the history restart separately.
  ! -------------------------------------------------------------------------------

  
  
  ! Set the tolerance on a matched latitude or longitude (100 meters)
  ! If this is a true history restart, the tolerance should be very small
  ! If this is an initialization, perhaps a new grid, using the closest,
  ! then there is no tolerance, ie 90.0 degrees, or some value which indicates
  ! the destination location is no-where close to the donor location

!!!  if (trim(runtype) == 'HISTORY' ) then

     ll_tolerance = (1.0/115.0)*(1.0/10.0)     ! 1/10th km
!!     ll_tolerance = (1.0/115.0)*(2.0/1.0)      ! 2km

     print*,"====================================================="
     print*,"         Entering Full State Initialization      "

!!!  else if (trim(runtype)=='INITIAL' .and. ied_init_mode==4 ) then

!!!     ll_tolerance = 20.0                       ! 20 km

!!     print*,"====================================================="
!!     print*," Entering Nearest Neighbor State File Initialization "


!!!  else
!!     call fatal_error ('Innapropriate run type encountered here'         &
 !!!         ,'init_full_history_restart','ed_history_io.f90')
!!!  end if


   
  ! at equator: (1 degree / 115 kilometers)  (1 km / 10 100-meter invervals)

  ! Open the HDF environment

  call h5open_f(hdferr)


  ! Turn off automatic error printing. This is done because there
  ! may be datasets that are not in the file, but it is OK. If
  ! data is missing that should be there, ED2 error reporting
  ! will detect it. If something is truly missing, the following
  ! call can be bypassed. Note, that automatic error reporting
  ! is turned back on at the end.
  
  call h5eset_auto_f(0,hdferr)


  ! Construct the file name for reinitiatlizing from
  ! The history file
  
  vnam = 'S'
  
  do ngr=1,ngrids
     
     cgrid => edgrid_g(ngr)
     

     !=======================================
     ! 1) Open the HDF5 HISTORY FILE
     !=======================================

     write(cgr,'(a1,i2.2)') 'g',ngr

     dbletime=dble(current_time%time)
     
     call makefnam(hnamel,sfilin(1),dbletime,current_time%year, &
          current_time%month,current_time%date,0,vnam,cgr,'h5 ')

     inquire(file=trim(hnamel),exist=exists)

     if (.not.exists) then
        call fatal_error ('File '//trim(hnamel)//' not found.'         &
                         ,'init_full_history_restart','ed_history_io.f90')
     else
        call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr < 0) then
           print *, 'Error opening HDF5 file - error - ',hdferr
           print *, '   Filename: ',trim(hnamel)
           call fatal_error('Error opening HDF5 file - error - '//trim(hnamel) &
                           ,'init_full_history_restart','ed_history_io.f90')
        end if
     end if


     !=======================================
     ! 2) Retrieve global vector sizes
     !=======================================

     !
     ! TO DO!!!!  INCLUDE NZG,NPFT,NBBH AS PART OF
     !            THE HDF5 HEADER, AND COMPARE
     !            AS SANITY CHECK TO THOSE VALUES
     !            READ IN THE NAMELIST. IF THEY
     !            DONT MATCH THEN THE NAMELIST
     !            ENTRY SHOULD BE CHANGED....
     !


     ! Only read in global grid data if this is a history restart
     ! otherwise, this data is not correct
     ! ==========================================================
     
     globdims = 0_8
     chnkdims = 0_8
     chnkoffs = 0_8
     
     globdims(1) = 1_8
     
     call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%npolygons_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%nsites_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%npatches_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%ncohorts_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     
     !=======================================
     ! 3) Retrieve the mapping of the data tree
     !=======================================
     
     globdims = 0_8
     globdims(1) = int(cgrid%npolygons_global,8)
     
     allocate(pysi_n(cgrid%npolygons_global))
     allocate(pysi_id(cgrid%npolygons_global))
     
     call h5dopen_f(file_id,'PYSI_N', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,pysi_n,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'PYSI_ID', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,pysi_id,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     globdims(1) = int(cgrid%nsites_global,8)
     
     allocate(sipa_n(cgrid%nsites_global))
     allocate(sipa_id(cgrid%nsites_global))
     
     call h5dopen_f(file_id,'SIPA_N', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,sipa_n,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'SIPA_ID', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,sipa_id,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     globdims(1) = int(cgrid%npatches_global,8)
     
     allocate(paco_n(cgrid%npatches_global))
     allocate(paco_id(cgrid%npatches_global))
     
     call h5dopen_f(file_id,'PACO_N', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,paco_n,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'PACO_ID', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,paco_id,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     
     ! ======================================
     ! 4) Retrieve the polygon coordinates data
     
     globdims(1) = int(cgrid%npolygons_global,8)
     allocate(file_lats(cgrid%npolygons_global))
     allocate(file_lons(cgrid%npolygons_global))
     
     call h5dopen_f(file_id,'LATITUDE', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_REAL,file_lats,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     call h5dopen_f(file_id,'LONGITUDE', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_REAL,file_lons,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     ! ======================================
     ! 5) Loop the polygons in the model state
     !    and match them with those int he file
     !    After the match. Walk through the 
     !    data from that polygon and initialize.
     !    A polygon match must have both latitudes
     !    and longitudes within 100 meters

     do ipy = 1,cgrid%npolygons
        
        py_index = 0

        cpoly => cgrid%polygon(ipy)

!!!        minrad = sqrt(2*(ll_tolerance**2))

        do ifpy = 1,cgrid%npolygons_global
           
!!!           currad = sqrt( (file_lats(ifpy)-cgrid%lat(ipy))**2 + (file_lons(ifpy)-cgrid%lon(ipy))**2 )

           if ( abs(file_lats(ifpy)-cgrid%lat(ipy)) < ll_tolerance .and. &
                abs(file_lons(ifpy)-cgrid%lon(ipy)) < ll_tolerance ) then !!! .and. &
!!!                (currad <  minrad) ) then
              py_index = ifpy
!!!              minrad   = currad
           end if
           
        enddo

        if (py_index==0) then
           print*,"COULD NOT MATCH A POLYGON WITH THE DATASET"
           print*,"STOPPING"
           print*,"THIS IS THE ",ipy,"th POLYGON"
           print*,"GRID LATS: ",cgrid%lat(ipy)
           print*,"GRID LONS: ",cgrid%lon(ipy)
!           print*,"FILE LATS: ",file_lats
!           print*,"FILE LONS: ",file_lons
           call fatal_error('Mismatch between polygon and dataset'         &
                           ,'init_full_history_restart','ed_history_io.f90')
        endif

        
        ! ========================================
        ! Get all necessary polygon variables
        ! associated with this index for the
        ! current polygon, scalar reads

        call fill_history_grid(cgrid,ipy,py_index)

        if (pysi_n(py_index) > 0) then
           
  !         print*,"Allocating: ",pysi_n(py_index)," sites from polygon",py_index
           
           call allocate_polygontype(cpoly,pysi_n(py_index))
           
           ! ========================================
           ! Get all necessary site variables
           ! associated with this index for the
           ! current polygon, vector reads
           
           call fill_history_polygon(cpoly,pysi_id(py_index),cgrid%nsites_global)
           
           do isi = 1,cpoly%nsites
              csite => cpoly%site(isi)
              
              ! Calculate the index of this site's data in the HDF
              si_index = pysi_id(py_index) + isi - 1
              
              if (sipa_n(si_index) > 0) then
                 
 !                print*,"Allocating: ",sipa_n(si_index)," patches from site ",si_index
                 
                 call allocate_sitetype(csite,sipa_n(si_index))

                 ! ========================================
                 ! Get all necessary patch variables
                 ! associated with this index for the
                 ! current site
                 
                 call fill_history_site(csite,sipa_id(si_index),cgrid%npatches_global)

                 csite%hcapveg = 0.

                 do ipa = 1,csite%npatches
                    cpatch => csite%patch(ipa)
                    
                    pa_index = sipa_id(si_index) + ipa - 1

                    if (paco_n(pa_index) > 0) then

!                       print*,"Allocating: ",paco_n(pa_index)," cohorts from patch ",pa_index
                       
                       call allocate_patchtype(cpatch,paco_n(pa_index))
                       
                       ! ========================================
                       ! Get all necessary patch variables
                       ! associated with this index for the
                       ! current site
                       
                       call fill_history_patch(cpatch,paco_id(pa_index),cgrid%ncohorts_global &
                                              ,cpoly%green_leaf_factor(:,isi))
                       
!                       do ipft = 1,n_pft
!                          csite%old_stoma_data_max(ipft,ipa)%recalc = 1
!                       enddo
                       do ico = 1,cpatch%ncohorts
                          csite%hcapveg(ipa) = csite%hcapveg(ipa) + cpatch%hcapveg(ico)
                       end do
                       

                    else

                       cpatch%ncohorts = 0
                       
                    endif
                 
                 enddo

              else

                 print*,"ATTEMPTING TO FILL SITE WITH PATCH VECTOR DATA"
                 print*,"NO PATCHES WERE FOUND in SIPA_N(SI_INDEX)"
                 print*,"THIS IS UNLIKELY AND MORALLY QUESTIONABLE."
                 stop
                 
              endif

           enddo
           
        else

           print*,"ATTEMPTING TO FILL A POLYGON WITH SITE VECTOR DATA"
           print*,"NO SITES WERE FOUND AT PYSI_N(PY_INDEX)"
           print*,"THIS IS EVEN MORE WRONG THAN HAVING NO PATCHES"
           print*,"AND THAT WAS REALLY REALLY WRONG"
           print*,"THIS IS WORSE THAN LETTING GIL BUY LUNCH FOR YOU"
           stop

        endif

        
     enddo


     call h5fclose_f(file_id, hdferr)
     if (hdferr.ne.0) then
         print*,"COULD NOT CLOSE THE HDF FILE"
         print*,hdferr
         stop	
     endif

     deallocate(file_lats,file_lons)
     deallocate(paco_n,paco_id)
     deallocate(sipa_n,sipa_id)
     deallocate(pysi_n,pysi_id )

  enddo

  ! Turn automatic error reporting back on.
  ! This is probably unecessary, because the environment
  ! is about to be flushed.

  call h5eset_auto_f(1,hdferr)


  ! Close the HDF environment
  
  call h5close_f(hdferr)

  ! Initialize the disturbance transition rates
  
  write(*,'(a,i2.2)')'    Initializing anthropogenic disturbance forcing. Node: ',mynum
  call landuse_init


  return
end subroutine init_full_history_restart
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid(cgrid,ipy,py_index)

  use ed_state_vars,only: edtype,polygontype
  use grid_coms,only : nzg
  use ed_max_dims,only : n_pft,n_dbh,n_age,n_dist_types
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize,datatype_id

  implicit none

  
#if USE_INTERF
  interface
     subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_r
     subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_d
     subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_i
  end interface
#endif

  type(edtype),target ::       cgrid

  integer,intent(in) :: ipy,py_index
  integer :: iparallel,dsetrank
  integer(SIZE_T) :: sz
  integer :: hdferr

  iparallel = 0
 

  globdims = 0_8
  chnkdims = 0_8
  chnkoffs = 0_8
  memoffs  = 0_8
  memdims  = 0_8
  memsize  = 1_8

  dsetrank = 1

  ! These are the dimensions in the filespace
  ! itself. Global is the size of the dataset,
  ! chnkoffs is the offset of the chunk we
  ! are going to read.  Chnkdims is the size
  ! of the slab that is to be read.
  
  globdims(1) = int(cgrid%npolygons_global,8)
  chnkdims(1) = 1_8
  chnkoffs(1) = int(py_index - 1,8)

  ! These are the dimensions for the memory space
  ! this should essentially be the same dimensioning
  ! as the buffer that we are filling. This routine
  ! is just filling a scalar point in a vector
  ! of polygons.

  memdims(1)  = 1_8
  memoffs(1)  = 0_8
  memsize(1)  = 1_8


  call hdf_getslab_d(cgrid%walltime_py(ipy:ipy),'WALLTIME_PY ',dsetrank,iparallel,.false.)

  call hdf_getslab_i(cgrid%lsl(ipy:ipy),'LSL ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.)
  
  call hdf_getslab_r(cgrid%Te(ipy:ipy),'TE ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%zbar(ipy:ipy),'ZBAR ',dsetrank,iparallel,.true.)

!!  call hdf_getslab_r(cgrid%tau(ipy:ipy),'TAU ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%sheat(ipy:ipy),'SHEAT ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%baseflow(ipy:ipy),'BASEFLOW ',dsetrank,iparallel,.true.)

  call hdf_getslab_i(cgrid%load_adjacency(ipy:ipy),'LOAD_ADJACENCY ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%swliq(ipy:ipy),'SWLIQ ',dsetrank,iparallel,.true.)
  
  ! All daily and monthly variables need to be retrieved if you are loading there...
  
  if (associated(cgrid%dmean_pcpg           ))                                             &
     call hdf_getslab_r(cgrid%dmean_pcpg           (ipy:ipy) ,'DMEAN_PCPG            '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_runoff         ))                                             &
     call hdf_getslab_r(cgrid%dmean_runoff         (ipy:ipy) ,'DMEAN_RUNOFF          '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_drainage       ))                                             &
     call hdf_getslab_r(cgrid%dmean_drainage       (ipy:ipy) ,'DMEAN_DRAINAGE        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_gpp            ))                                             &
     call hdf_getslab_r(cgrid%dmean_gpp            (ipy:ipy) ,'DMEAN_GPP             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_evap           ))                                             &
     call hdf_getslab_r(cgrid%dmean_evap           (ipy:ipy) ,'DMEAN_EVAP            '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_transp         ))                                             &
     call hdf_getslab_r(cgrid%dmean_transp         (ipy:ipy) ,'DMEAN_TRANSP          '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_sensible_vc    ))                                             &
     call hdf_getslab_r(cgrid%dmean_sensible_vc    (ipy:ipy) ,'DMEAN_SENSIBLE_VC     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_sensible_gc    ))                                             &
     call hdf_getslab_r(cgrid%dmean_sensible_gc    (ipy:ipy) ,'DMEAN_SENSIBLE_GC     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_sensible_ac    ))                                             &
     call hdf_getslab_r(cgrid%dmean_sensible_ac    (ipy:ipy) ,'DMEAN_SENSIBLE_AC     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vapor_vc       ))                                             &
     call hdf_getslab_r(cgrid%dmean_vapor_vc       (ipy:ipy) ,'DMEAN_VAPOR_VC        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vapor_gc       ))                                             &
     call hdf_getslab_r(cgrid%dmean_vapor_gc       (ipy:ipy) ,'DMEAN_VAPOR_GC        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vapor_ac       ))                                             &
     call hdf_getslab_r(cgrid%dmean_vapor_ac       (ipy:ipy) ,'DMEAN_VAPOR_AC        '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_nep            ))                                             &
     call hdf_getslab_r(cgrid%dmean_nep            (ipy:ipy) ,'DMEAN_NEP             '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_plresp         ))                                             &
     call hdf_getslab_r(cgrid%dmean_plresp         (ipy:ipy) ,'DMEAN_PLRESP          '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_rh             ))                                             &
     call hdf_getslab_r(cgrid%dmean_rh             (ipy:ipy) ,'DMEAN_RH              '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_leaf_resp      ))                                             &
     call hdf_getslab_r(cgrid%dmean_leaf_resp      (ipy:ipy) ,'DMEAN_LEAF_RESP       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_root_resp      ))                                             &
     call hdf_getslab_r(cgrid%dmean_root_resp      (ipy:ipy) ,'DMEAN_ROOT_RESP       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_growth_resp    ))                                             &
     call hdf_getslab_r(cgrid%dmean_growth_resp    (ipy:ipy) ,'DMEAN_GROWTH_RESP     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_storage_resp   ))                                             &
     call hdf_getslab_r(cgrid%dmean_storage_resp   (ipy:ipy) ,'DMEAN_STORAGE_RESP    '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vleaf_resp     ))                                             &
     call hdf_getslab_r(cgrid%dmean_vleaf_resp     (ipy:ipy) ,'DMEAN_VLEAF_RESP      '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_fs_open        ))                                             &
     call hdf_getslab_r(cgrid%dmean_fs_open        (ipy:ipy) ,'DMEAN_FS_OPEN         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_fsw            ))                                             &
     call hdf_getslab_r(cgrid%dmean_fsw            (ipy:ipy) ,'DMEAN_FSW             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_fsn            ))                                             &
     call hdf_getslab_r(cgrid%dmean_fsn            (ipy:ipy) ,'DMEAN_FSN             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_temp       ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_temp       (ipy:ipy) ,'DMEAN_CAN_TEMP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_shv        ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_shv        (ipy:ipy) ,'DMEAN_CAN_SHV         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_prss       ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_prss       (ipy:ipy) ,'DMEAN_CAN_PRSS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_theta      ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_theta      (ipy:ipy) ,'DMEAN_CAN_THETA       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_theiv      ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_theiv      (ipy:ipy) ,'DMEAN_CAN_THEIV       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_co2        ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_co2        (ipy:ipy) ,'DMEAN_CAN_CO2         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_rhos       ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_rhos       (ipy:ipy) ,'DMEAN_CAN_RHOS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_energy     ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_energy     (ipy:ipy) ,'DMEAN_VEG_ENERGY      '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_water      ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_water      (ipy:ipy) ,'DMEAN_VEG_WATER       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_hcap       ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_hcap       (ipy:ipy) ,'DMEAN_VEG_HCAP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_temp       ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_temp       (ipy:ipy) ,'DMEAN_VEG_TEMP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_temp       ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_temp       (ipy:ipy) ,'DMEAN_ATM_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
  
  if (associated(cgrid%dmean_rshort       ))                                               &
       call hdf_getslab_r(cgrid%dmean_rshort       (ipy:ipy) ,'DMEAN_RSHORT        '       &
       ,dsetrank,iparallel,.false.)
  
  if (associated(cgrid%dmean_rlong       ))                                                &
       call hdf_getslab_r(cgrid%dmean_rlong       (ipy:ipy) ,'DMEAN_RLONG        '         &
       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_shv        ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_shv        (ipy:ipy) ,'DMEAN_ATM_SHV         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_prss       ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_prss       (ipy:ipy) ,'DMEAN_ATM_PRSS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_vels       ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_vels       (ipy:ipy) ,'DMEAN_ATM_VELS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_co2_residual   ))                                             &
     call hdf_getslab_r(cgrid%dmean_co2_residual   (ipy:ipy) ,'DMEAN_CO2_RESIDUAL    '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_energy_residual))                                             &
     call hdf_getslab_r(cgrid%dmean_energy_residual(ipy:ipy) ,'DMEAN_ENERGY_RESIDUAL '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_water_residual ))                                             &
     call hdf_getslab_r(cgrid%dmean_water_residual (ipy:ipy) ,'DMEAN_WATER_RESIDUAL  '     &
                       ,dsetrank,iparallel,.false.)


  if (associated(cgrid%mmean_co2_residual   ))                                             &
     call hdf_getslab_r(cgrid%mmean_co2_residual   (ipy:ipy) ,'MMEAN_CO2_RESIDUAL    '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_energy_residual))                                             &
     call hdf_getslab_r(cgrid%mmean_energy_residual(ipy:ipy) ,'MMEAN_ENERGY_RESIDUAL '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_water_residual ))                                             &
     call hdf_getslab_r(cgrid%mmean_water_residual (ipy:ipy) ,'MMEAN_WATER_RESIDUAL  '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_gpp            ))                                             &
     call hdf_getslab_r(cgrid%mmean_gpp            (ipy:ipy) ,'MMEAN_GPP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_evap           ))                                             &
     call hdf_getslab_r(cgrid%mmean_evap           (ipy:ipy) ,'MMEAN_EVAP            '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_transp         ))                                             &
     call hdf_getslab_r(cgrid%mmean_transp         (ipy:ipy) ,'MMEAN_TRANSP          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_sensible_vc    ))                                             &
     call hdf_getslab_r(cgrid%mmean_sensible_vc    (ipy:ipy) ,'MMEAN_SENSIBLE_VC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_sensible_gc    ))                                             &
     call hdf_getslab_r(cgrid%mmean_sensible_gc    (ipy:ipy) ,'MMEAN_SENSIBLE_GC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_sensible_ac    ))                                             &
     call hdf_getslab_r(cgrid%mmean_sensible_ac    (ipy:ipy) ,'MMEAN_SENSIBLE_AC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vapor_vc       ))                                             &
     call hdf_getslab_r(cgrid%mmean_vapor_vc       (ipy:ipy) ,'MMEAN_VAPOR_VC        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vapor_gc       ))                                             &
     call hdf_getslab_r(cgrid%mmean_vapor_gc       (ipy:ipy) ,'MMEAN_VAPOR_GC        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vapor_ac       ))                                             &
     call hdf_getslab_r(cgrid%mmean_vapor_ac       (ipy:ipy) ,'MMEAN_VAPOR_AC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_nep            ))                                             &
     call hdf_getslab_r(cgrid%mmean_nep            (ipy:ipy) ,'MMEAN_NEP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_plresp         ))                                             &
     call hdf_getslab_r(cgrid%mmean_plresp         (ipy:ipy) ,'MMEAN_PLRESP          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_rh             ))                                             &
     call hdf_getslab_r(cgrid%mmean_rh             (ipy:ipy) ,'MMEAN_RH              '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_leaf_resp      ))                                             &
     call hdf_getslab_r(cgrid%mmean_leaf_resp      (ipy:ipy) ,'MMEAN_LEAF_RESP       '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_root_resp      ))                                             &
     call hdf_getslab_r(cgrid%mmean_root_resp      (ipy:ipy) ,'MMEAN_ROOT_RESP       '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_growth_resp    ))                                             &
     call hdf_getslab_r(cgrid%mmean_growth_resp    (ipy:ipy) ,'MMEAN_GROWTH_RESP     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_storage_resp   ))                                             &
     call hdf_getslab_r(cgrid%mmean_storage_resp   (ipy:ipy) ,'MMEAN_STORAGE_RESP    '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vleaf_resp     ))                                             &
     call hdf_getslab_r(cgrid%mmean_vleaf_resp     (ipy:ipy) ,'MMEAN_VLEAF_RESP      '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_temp       ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_temp       (ipy:ipy) ,'MMEAN_CAN_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_shv        ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_shv        (ipy:ipy) ,'MMEAN_CAN_SHV         '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_co2        ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_co2        (ipy:ipy) ,'MMEAN_CAN_CO2         '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_rhos       ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_rhos       (ipy:ipy) ,'MMEAN_CAN_RHOS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_can_prss       ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_prss       (ipy:ipy) ,'MMEAN_CAN_PRSS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_can_theta      ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_theta      (ipy:ipy) ,'MMEAN_CAN_THETA       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_can_theiv      ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_theiv      (ipy:ipy) ,'MMEAN_CAN_THEIV       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_veg_energy     ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_energy     (ipy:ipy) ,'MMEAN_VEG_ENERGY      '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_veg_water      ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_water      (ipy:ipy) ,'MMEAN_VEG_WATER       '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_veg_temp       ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_temp       (ipy:ipy) ,'MMEAN_VEG_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_veg_hcap       ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_hcap       (ipy:ipy) ,'MMEAN_VEG_HCAP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_rshort       ))                                               &
       call hdf_getslab_r(cgrid%mmean_rshort       (ipy:ipy) ,'MMEAN_RSHORT        '       &
       ,dsetrank,iparallel,.false.)
  
  if (associated(cgrid%mmean_rlong       ))                                                &
       call hdf_getslab_r(cgrid%mmean_rlong       (ipy:ipy) ,'MMEAN_RLONG        '         &
       ,dsetrank,iparallel,.false.)
  
  if (associated(cgrid%mmean_atm_temp       ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_temp       (ipy:ipy) ,'MMEAN_ATM_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_shv        ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_shv        (ipy:ipy) ,'MMEAN_ATM_SHV         '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_prss       ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_prss       (ipy:ipy) ,'MMEAN_ATM_PRSS        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_vels       ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_vels       (ipy:ipy) ,'MMEAN_ATM_VELS        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_pcpg           ))                                             &
     call hdf_getslab_r(cgrid%mmean_pcpg           (ipy:ipy) ,'MMEAN_PCPG            '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_runoff         ))                                             &
     call hdf_getslab_r(cgrid%mmean_runoff         (ipy:ipy) ,'MMEAN_RUNOFF          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_drainage       ))                                             &
     call hdf_getslab_r(cgrid%mmean_drainage       (ipy:ipy) ,'MMEAN_DRAINAGE        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_fs_open        ))                                             &
     call hdf_getslab_r(cgrid%mmean_fs_open        (ipy:ipy) ,'MMEAN_FS_OPEN         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_fsw            ))                                             &
     call hdf_getslab_r(cgrid%mmean_fsw            (ipy:ipy) ,'MMEAN_FSW             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_fsn            ))                                             &
     call hdf_getslab_r(cgrid%mmean_fsn            (ipy:ipy) ,'MMEAN_FSN             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%stdev_gpp            ))                                             &
     call hdf_getslab_r(cgrid%stdev_gpp            (ipy:ipy) ,'STDEV_GPP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_evap           ))                                             &
     call hdf_getslab_r(cgrid%stdev_evap           (ipy:ipy) ,'STDEV_EVAP            '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_transp         ))                                             &
     call hdf_getslab_r(cgrid%stdev_transp         (ipy:ipy) ,'STDEV_TRANSP          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_sensible       ))                                             &
     call hdf_getslab_r(cgrid%stdev_sensible       (ipy:ipy) ,'STDEV_SENSIBLE        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_nep            ))                                             &
     call hdf_getslab_r(cgrid%stdev_nep            (ipy:ipy) ,'STDEV_NEP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_rh             ))                                             &
     call hdf_getslab_r(cgrid%stdev_rh             (ipy:ipy) ,'STDEV_RH              '     &
                       ,dsetrank,iparallel,.false.)

   ! Variables with 2 dimensions (nzg,npolygons)
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims(1)  = int(nzg,8)
   memsize(1)  = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   call hdf_getslab_i(cgrid%ntext_soil(:,ipy)          ,'NTEXT_SOIL '       ,&
        dsetrank,iparallel,.false.)

   if(associated(cgrid%dmean_soil_temp)) &
      call hdf_getslab_r(cgrid%dmean_soil_temp(:,ipy)  ,'DMEAN_SOIL_TEMP '  ,&
      dsetrank,iparallel,.false.)

   if(associated(cgrid%dmean_soil_water)) &
      call hdf_getslab_r(cgrid%dmean_soil_water(:,ipy) ,'DMEAN_SOIL_WATER ' ,&
      dsetrank,iparallel,.false.)

   if(associated(cgrid%mmean_soil_temp)) &
      call hdf_getslab_r(cgrid%mmean_soil_temp(:,ipy)  ,'MMEAN_SOIL_TEMP '  ,&
      dsetrank,iparallel,.false.)

   if(associated(cgrid%mmean_soil_water)) &
      call hdf_getslab_r(cgrid%mmean_soil_water(:,ipy) ,'MMEAN_SOIL_WATER ' ,&
      dsetrank,iparallel,.false.)


   ! Variables with 2 dimensions (n_pft,npolygons)
   dsetrank    = 2
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims(1)  = int(n_pft,8)
   memsize(1)  = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   if(associated(cgrid%bseeds_pft)) call hdf_getslab_r(cgrid%bseeds_pft(:,ipy) ,'BSEEDS_PFT '    , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%lai_pft)) call hdf_getslab_r(cgrid%lai_pft(:,ipy) ,'LAI_PFT '       , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%wpa_pft)) call hdf_getslab_r(cgrid%wpa_pft(:,ipy) ,'WPA_PFT '       , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%wai_pft)) call hdf_getslab_r(cgrid%wai_pft(:,ipy) ,'WAI_PFT '       , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_lai_pft)) call hdf_getslab_r(cgrid%mmean_lai_pft(:,ipy) ,'MMEAN_LAI_PFT ' , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_wpa_pft)) call hdf_getslab_r(cgrid%mmean_wpa_pft(:,ipy) ,'MMEAN_WPA_PFT ' , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_wai_pft)) call hdf_getslab_r(cgrid%mmean_wai_pft(:,ipy) ,'MMEAN_WAI_PFT ' , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%agb_pft)) call hdf_getslab_r(cgrid%agb_pft(:,ipy) ,'AGB_PFT '       , &
        dsetrank,iparallel,.true.)
   if(associated(cgrid%ba_pft)) call hdf_getslab_r(cgrid%ba_pft(:,ipy) ,'BA_PFT '        ,   &
        dsetrank,iparallel,.true.)


   ! Variables with 2 dimensions (n_dbh,npolygons)
   dsetrank    = 2
   globdims(1) = int(n_dbh,8)
   chnkdims(1) = int(n_dbh,8)
   memdims(1)  = int(n_dbh,8)
   memsize(1)  = int(n_dbh,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   if(associated(cgrid%dmean_gpp_dbh)) call hdf_getslab_r(cgrid%dmean_gpp_dbh(:,ipy) , &
        'DMEAN_GPP_DBH ' ,dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_gpp_dbh)) call hdf_getslab_r(cgrid%mmean_gpp_dbh(:,ipy) , &
        'MMEAN_GPP_DBH ' ,dsetrank,iparallel,.false.)

   ! Variables with 2 dimensions (13,npolygons)
   dsetrank    = 2
   globdims(1) = int(13,8)
   chnkdims(1) = int(13,8)
   memdims(1)  = int(13,8)
   memsize(1)  = int(13,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   if(associated(cgrid%workload)) call hdf_getslab_r(cgrid%workload(:,ipy) , &
        'WORKLOAD ' ,dsetrank,iparallel,.true.)


   ! Variables with three dimensions(n_dist_types,n_dist_types,npolygons)
   dsetrank    = 3
   globdims(1) = int(n_dist_types,8)
   chnkdims(1) = int(n_dist_types,8)
   memdims(1)  = int(n_dist_types,8)
   memsize(1)  = int(n_dist_types,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2) = int(n_dist_types,8)
   chnkdims(2) = int(n_dist_types,8)
   memdims(2)  = int(n_dist_types,8)
   memsize(2)  = int(n_dist_types,8)
   chnkoffs(2) = 0_8
   memoffs(2)  = 0_8

   globdims(3)  = int(cgrid%npolygons_global,8)
   chnkdims(3)  = 1_8
   chnkoffs(3)  = int(py_index - 1,8)
   memdims(3)   = 1_8
   memsize(3)   = 1_8
   memoffs(3)   = 0_8
   if (associated(cgrid%disturbance_rates))                                                &
       call hdf_getslab_r(cgrid%disturbance_rates(:,:,ipy),'DISTURBANCE_RATES '            &
                         ,dsetrank,iparallel,.false.)

   return
 end subroutine fill_history_grid
 !=========================================================================================!
 !=========================================================================================!






 !==========================================================================================!
 !==========================================================================================!
 subroutine fill_history_polygon(cpoly,pysi_index,nsites_global)

   use ed_state_vars,only: polygontype
   use hdf5
   use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
        globdims,chnkdims,chnkoffs,cnt,stride, &
        memdims,memoffs,memsize

   use grid_coms,only : nzg
   use ed_max_dims,only : n_pft,n_dbh,n_dist_types

   implicit none

#if USE_INTERF
     interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in):: required
      end subroutine hdf_getslab_r
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_d
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_i
   end interface
#endif

   type(polygontype),target :: cpoly
   integer,intent(in) :: pysi_index
   integer,intent(in) :: nsites_global
   integer :: iparallel
   integer :: dsetrank

   iparallel = 0

   dsetrank = 1_8
   globdims = 0_8
   chnkdims = 0_8
   chnkoffs = 0_8
   memoffs  = 0_8
   memdims  = 0_8
   memsize  = 1_8

   globdims(1) = int(nsites_global,8)
   chnkdims(1) = int(cpoly%nsites,8)
   chnkoffs(1) = int(pysi_index - 1,8)
   memdims(1)  = int(cpoly%nsites,8)
   memsize(1)  = int(cpoly%nsites,8)
   memoffs(1)  = 0_8

   call hdf_getslab_i(cpoly%patch_count,'PATCH_COUNT ',dsetrank,iparallel,.true.)  
   call hdf_getslab_i(cpoly%sitenum,'SITENUM ',dsetrank,iparallel,.true.)

   call hdf_getslab_i(cpoly%lsl,'LSL_SI ',dsetrank,iparallel,.true.)   
   call hdf_getslab_r(cpoly%area,'AREA_SI ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%patch_area,'PATCH_AREA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%elevation,'ELEVATION ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%slope,'SLOPE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%aspect,'ASPECT ',dsetrank,iparallel,.true.)

   call hdf_getslab_i(cpoly%num_landuse_years,'NUM_LANDUSE_YEARS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%TCI,'TCI ',dsetrank,iparallel,.true.)      
   call hdf_getslab_r(cpoly%pptweight,'pptweight ',dsetrank,iparallel,.true.)      
   call hdf_getslab_i(cpoly%hydro_next,'HYDRO_NEXT ',dsetrank,iparallel,.true.)
   call hdf_getslab_i(cpoly%hydro_prev,'HYDRO_PREV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%moist_W,'MOIST_W ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%moist_f,'MOIST_F ',dsetrank,iparallel,.true.)  
   call hdf_getslab_r(cpoly%moist_tau,'MOIST_TAU ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%moist_zi,'MOIST_ZI ',dsetrank,iparallel,.true.) 
   call hdf_getslab_r(cpoly%baseflow,'BASEFLOW_SI ',dsetrank,iparallel,.true.) 
!   call hdf_getslab_i(cpoly%metplex_beg_month,'METPLEX_BEG_MONTH ',dsetrank,iparallel,.true.)
!   call hdf_getslab_i(cpoly%metplex_beg_year,'METPLEX_BEG_YEAR ',dsetrank,iparallel,.true.)
!   call hdf_getslab_i(cpoly%metplex_end_year,'METPLEX_END_YEAR ',dsetrank,iparallel,.true.)

   call hdf_getslab_r(cpoly%min_monthly_temp,'MIN_MONTHLY_TEMP ',dsetrank,iparallel,.true.)
!   call hdf_getslab_r(cpoly%removed_biomass,'REMOVED_BIOMASS ',dsetrank,iparallel,.true.) 
!   call hdf_getslab_r(cpoly%harvested_biomass,'HARVESTED_BIOMASS ', &
!        dsetrank,iparallel,.true.) 
   call hdf_getslab_i(cpoly%plantation,'PLANTATION_SI ',dsetrank,iparallel,.true.) 
   call hdf_getslab_i(cpoly%agri_stocking_pft,'AGRI_STOCKING_PFT ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agri_stocking_density,'AGRI_STOCKING_DENSITY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_i(cpoly%plantation_stocking_pft,'PLANTATION_STOCKING_PFT ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%plantation_stocking_density,'PLANTATION_STOCKING_DENSITY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%primary_harvest_memory,'PRIMARY_HARVEST_MEMORY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%secondary_harvest_memory,'SECONDARY_HARVEST_MEMORY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%fire_disturbance_rate,'FIRE_DISTURBANCE_RATE ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%ignition_rate,'IGNITION_RATE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%treefall_disturbance_rate,'TREEFALL_DISTURBANCE_RATE ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%nat_disturbance_rate,'NAT_DISTURBANCE_RATE ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_i(cpoly%nat_dist_type,'NAT_DIST_TYPE ',dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims(1)  = int(n_pft,8)
   memsize(1)  = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   if (associated(cpoly%lai_pft)) call hdf_getslab_r(cpoly%lai_pft,'LAI_PFT_SI ', &
        dsetrank,iparallel,.false.)
   if (associated(cpoly%wpa_pft)) call hdf_getslab_r(cpoly%wpa_pft,'WPA_PFT_SI ', &
        dsetrank,iparallel,.false.)
   if (associated(cpoly%wai_pft)) call hdf_getslab_r(cpoly%wai_pft,'WAI_PFT_SI ', &
        dsetrank,iparallel,.false.)
   call hdf_getslab_r(cpoly%green_leaf_factor,'GREEN_LEAF_FACTOR ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%leaf_aging_factor,'LEAF_AGING_FACTOR ', &
        dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims(1)  = int(nzg,8)
   memsize(1)  = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   call hdf_getslab_i(cpoly%ntext_soil,'NTEXT_SOIL_SI ',dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = 12_8
   chnkdims(1) = 12_8
   memdims(1)  = 12_8
   memsize(1)  = 12_8
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   call hdf_getslab_r(cpoly%lambda_fire,'LAMBDA_FIRE ',dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = int(n_dist_types,8)
   chnkdims(1) = int(n_dist_types,8)
   memdims(1)  = int(n_dist_types,8)
   memsize(1)  = int(n_dist_types,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   call hdf_getslab_r(cpoly%loss_fraction,'LOSS_FRACTION ',dsetrank,iparallel,.true.)

   dsetrank    = 3_8
   globdims(1:2) = int(n_dist_types,8)
   chnkdims(1:2) = int(n_dist_types,8)
   memdims(1:2)  = int(n_dist_types,8)
   memsize(1:2)  = int(n_dist_types,8)
   chnkoffs(1:2) = 0
   memoffs(1:2)  = 0
   globdims(3)  = int(nsites_global,8)
   chnkdims(3)  = int(cpoly%nsites,8)
   chnkoffs(3)  = int(pysi_index - 1,8)
   memdims(3)   = int(cpoly%nsites,8)
   memsize(3)   = int(cpoly%nsites,8)
   memoffs(3)   = 0

   call hdf_getslab_r(cpoly%disturbance_memory,'DISTURBANCE_MEMORY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%disturbance_rates,'DISTURBANCE_RATES_SI ', &
        dsetrank,iparallel,.true.)

   dsetrank    = 3
   globdims(3) = int(nsites_global,8)
   chnkdims(3) = int(cpoly%nsites,8)
   chnkoffs(3) = int(pysi_index - 1,8)
   memdims(3)  = int(cpoly%nsites,8)
   memsize(3)  = int(cpoly%nsites,8)
   memoffs(3)  = 0
   globdims(2) = int(n_dbh,8)
   chnkdims(2) = int(n_dbh,8)
   memdims(2)  = int(n_dbh,8)
   memsize(2)  = int(n_dbh,8)
   chnkoffs(2) = 0
   memoffs(2)  = 0
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims(1)  = int(n_pft,8)
   memsize(1)  = int(n_pft,8)
   chnkoffs(1) = 0
   memoffs(1)  = 0

   call hdf_getslab_r(cpoly%basal_area,'BASAL_AREA_SI ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb,'AGB_SI ',dsetrank,iparallel,.true.)
!   call hdf_getslab_r(cpoly%pldens,'PLDENS_SI ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(cpoly%basal_area_growth,'BASAL_AREA_GROWTH ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb_growth,'AGB_GROWTH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%basal_area_mort,'BASAL_AREA_MORT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%basal_area_cut,'BASAL_AREA_CUT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb_mort,'AGB_MORT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb_cut,'AGB_CUT ',dsetrank,iparallel,.true.)


   return
 end subroutine fill_history_polygon
 !==========================================================================================!
 !==========================================================================================!






 !==========================================================================================!
 !==========================================================================================!
 subroutine fill_history_site(csite,sipa_index,npatches_global)

   use ed_state_vars,only: sitetype
   use grid_coms,only : nzg,nzs
   use c34constants,only:n_stoma_atts
   use ed_max_dims,only : n_pft,n_dbh
   use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
        globdims,chnkdims,chnkoffs,cnt,stride, &
        memdims,memoffs,memsize,datatype_id,setsize
   use fusion_fission_coms, only: ff_ndbh
   use hdf5

   implicit none

#if USE_INTERF
     interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_r
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_d
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_i
   end interface
#endif

   type(sitetype),target :: csite
   integer,intent(in) :: sipa_index
   integer,intent(in) :: npatches_global
   integer :: iparallel
   integer :: dsetrank
   integer :: hdferr
   integer :: ipa,ipft
   real(kind=8),allocatable, dimension(:,:) ::  buff

   iparallel = 0

   dsetrank = 1
   globdims = 0
   chnkdims = 0
   chnkoffs = 0
   memoffs  = 0
   memdims  = 0
   memsize  = 1

   ! These are the dimensions in the filespace
   ! itself. Global is the size of the dataset,
   ! chnkoffs is the offset of the chunk we
   ! are going to read.  Chnkdims is the size
   ! of the slab that is to be read.

   globdims(1) = int(npatches_global,8)
   chnkdims(1) = int(csite%npatches,8)
   chnkoffs(1) = int(sipa_index - 1,8)

   memdims(1)  = int(csite%npatches,8)
   memsize(1)  = int(csite%npatches,8)
   memoffs(1)  = 0

   call hdf_getslab_i(csite%dist_type,'DIST_TYPE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%age,'AGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%area,'AREA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fast_soil_C,'FAST_SOIL_C ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%slow_soil_C,'SLOW_SOIL_C ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%structural_soil_C,'STRUCTURAL_SOIL_C ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%structural_soil_L,'STRUCTURAL_SOIL_L ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fast_soil_N,'FAST_SOIL_N ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sum_dgd,'SUM_DGD ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sum_chd,'SUM_CHD ',dsetrank,iparallel,.true.)
   call hdf_getslab_i(csite%plantation,'PLANTATION ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_i(csite%cohort_count,'COHORT_COUNT ',dsetrank,iparallel)
   call hdf_getslab_r(csite%can_theiv,'CAN_THEIV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_prss,'CAN_PRSS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_theta,'CAN_THETA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_temp,'CAN_TEMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_shv,'CAN_SHV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_co2,'CAN_CO2 ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_rhos,'CAN_RHOS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_depth,'CAN_DEPTH ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_i(csite%pname,'PNAME ',dsetrank,iparallel)
   call hdf_getslab_r(csite%lai,'LAI_PA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wpa,'WPA_PA ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%wai,'WAI_PA ',dsetrank,iparallel,.false.)
   call hdf_getslab_i(csite%nlev_sfcwater,'NLEV_SFCWATER ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ground_shv,'GROUND_SHV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%surface_ssh,'SURFACE_SSH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rough,'ROUGH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%avg_daily_temp,'AVG_DAILY_TEMP ',dsetrank,iparallel,.true.)  
   call hdf_getslab_r(csite%mean_rh,'MEAN_RH ',dsetrank,iparallel,.true.)

   if (associated(csite%dmean_rh       )) &
        call hdf_getslab_r(csite%dmean_rh,'DMEAN_RH_PA ',dsetrank,iparallel,.false.)
   if (associated(csite%mmean_rh       )) &
        call hdf_getslab_r(csite%mmean_rh,'MMEAN_RH_PA ',dsetrank,iparallel,.false.)

   call hdf_getslab_r(csite%lambda_light,'LAMBDA_LIGHT ',dsetrank,iparallel,.true.)

   if (associated(csite%dmean_lambda_light       )) &
        call hdf_getslab_r(csite%dmean_lambda_light,'DMEAN_LAMBDA_LIGHT ',dsetrank,iparallel,.false.)

   if (associated(csite%mmean_lambda_light       )) &
        call hdf_getslab_r(csite%mmean_lambda_light,'MMEAN_LAMBDA_LIGHT ',dsetrank,iparallel,.false.)
   
   call hdf_getslab_r(csite%mean_nep,'MEAN_NEP ',dsetrank,iparallel,.true.)

   call hdf_getslab_r(csite%wbudget_loss2atm,'WBUDGET_LOSS2ATM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_denseffect,'WBUDGET_DENSEFFECT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_precipgain,'WBUDGET_PRECIPGAIN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_loss2runoff,'WBUDGET_LOSS2RUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_initialstorage,'WBUDGET_INITIALSTORAGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_loss2atm,'EBUDGET_LOSS2ATM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_denseffect,'EBUDGET_DENSEFFECT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_loss2runoff,'EBUDGET_LOSS2RUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_netrad,'EBUDGET_NETRAD ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_precipgain,'EBUDGET_PRECIPGAIN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_initialstorage,'EBUDGET_INITIALSTORAGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_initialstorage,'CO2BUDGET_INITIALSTORAGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_loss2atm,'CO2BUDGET_LOSS2ATM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_denseffect,'CO2BUDGET_DENSEFFECT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_gpp,'CO2BUDGET_GPP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_plresp,'CO2BUDGET_PLRESP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_rh,'CO2BUDGET_RH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%today_A_decomp,'TODAY_A_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%today_Af_decomp,'TODAY_AF_DECOMP ',dsetrank,iparallel,.true.)
   if (associated(csite%dmean_A_decomp       )) &
        call hdf_getslab_r(csite%dmean_A_decomp,'DMEAN_A_DECOMP ',dsetrank,iparallel,.false.)
   if (associated(csite%dmean_Af_decomp       )) &
        call hdf_getslab_r(csite%dmean_Af_decomp,'DMEAN_AF_DECOMP ',dsetrank,iparallel,.false.)
   if (associated(csite%mmean_A_decomp       )) &
        call hdf_getslab_r(csite%mmean_A_decomp,'MMEAN_A_DECOMP ',dsetrank,iparallel,.false.)
   if (associated(csite%mmean_Af_decomp       )) &
   call hdf_getslab_r(csite%mmean_Af_decomp,'MMEAN_AF_DECOMP ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%veg_rough,'VEG_ROUGH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%veg_height ,'VEG_HEIGHT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fsc_in,'FSC_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ssc_in,'SSC_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ssl_in,'SSL_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fsn_in,'FSN_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%total_plant_nitrogen_uptake,'TOTAL_PLANT_NITROGEN_UPTAKE ',dsetrank,iparallel,.true.)
   
   !  call hdf_getslab_r(csite%rshort_g,'RSHORT_G ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rshort_g_beam,'RSHORT_G_BEAM ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rshort_g_diffuse,'RSHORT_G_DIFFUSE ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rlong_g,'RLONG_G ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rlong_g_surf,'RLONG_G_SURF ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_g_incid,'RLONG_G_INCID ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_s,'RLONG_S ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_s_surf,'RLONG_S_SURF ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_s_incid,'RLONG_S_INCID ',dsetrank,iparallel)
   
   !  call hdf_getslab_r(csite%albedt,'ALBEDT ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%albedo_beam,'ALBEDO_BEAM ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%albedo_diffuse,'ALBEDO_DIFFUSE ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlongup,'RLONGUP ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_albedo,'RLONGUP_ALBEDO ',dsetrank,iparallel)
   call hdf_getslab_r(csite%total_snow_depth,'TOTAL_SNOW_DEPTH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%snowfac,'SNOWFAC ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%A_decomp,'A_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%f_decomp,'F_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rh,'RH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%cwd_rh,'CWD_RH ',dsetrank,iparallel,.true.)
   call hdf_getslab_i(csite%fuse_flag,'FUSE_FLAG ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%plant_ag_biomass,'PLANT_AG_BIOMASS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_wflux,'MEAN_WFLUX ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_latflux,'MEAN_LATFLUX ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_hflux,'MEAN_HFLUX ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_runoff,'MEAN_RUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_qrunoff,'MEAN_QRUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%htry,'HTRY ',dsetrank,iparallel,.true.)
   if (associated(csite%dmean_rk4step)) &
        call hdf_getslab_r(csite%dmean_rk4step,'DMEAN_RK4STEP ',dsetrank,iparallel,.false.)
   if (associated(csite%mmean_rk4step)) &
        call hdf_getslab_r(csite%mmean_rk4step,'MMEAN_RK4STEP ',dsetrank,iparallel,.false.)
   
   dsetrank    = 2
   globdims(1) = int(nzs,8)
   chnkdims(1) = int(nzs,8)
   memdims(1)  = int(nzs,8)
   memsize(1)  = int(nzs,8)
   chnkoffs(1) = 0
   memoffs(1)  = 0
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches,8)
   chnkoffs(2) = int(sipa_index - 1,8)
   memdims(2)  = int(csite%npatches,8)
   memsize(2)  = int(csite%npatches,8)
   memoffs(2)  = 0
   
   call hdf_getslab_r(csite%sfcwater_mass,'SFCWATER_MASS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_energy,'SFCWATER_ENERGY ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_depth,'SFCWATER_DEPTH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rshort_s,'RSHORT_S ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rshort_s_beam,'RSHORT_S_BEAM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rshort_s_diffuse,'RSHORT_S_DIFFUSE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_tempk,'SFCWATER_TEMPK ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_fracliq,'SFCWATER_FRACLIQ ',dsetrank,iparallel,.true.)
   
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims(1)  = int(nzg,8)
   memsize(1)  = int(nzg,8)
   chnkoffs(1) = 0
   memoffs(1)  = 0
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches,8)
   chnkoffs(2) = int(sipa_index - 1,8)
   memdims(2)  = int(csite%npatches,8)
   memsize(2)  = int(csite%npatches,8)
   memoffs(2)  = 0
   
   call hdf_getslab_i(csite%ntext_soil,'NTEXT_SOIL_PA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%soil_energy,'SOIL_ENERGY_PA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%soil_water,'SOIL_WATER_PA ',dsetrank,iparallel,.true.)

   !-----------------------------------------------------------------------------------!
   !  Soil water is double precision, although it may not be DP in the dataset
   !  The following lines make provisions for this by testing the dataset.
   
!   call h5dopen_f(file_id,'SOIL_WATER_PA ', dset_id, hdferr)
!   if (hdferr /= 0 ) then
!      call fatal_error('Dataset did not have soil water?' &
!           ,'fill_history_site','ed_history_io.f90')
!   endif
   ! ---------------------------------------------------------------------------------!
   ! THESE LINES ARE USEFULL FOR DETERMINING DATA SIZE OF ANY GIVEN OBJECT IN A SET   !
   ! ---------------------------------------------------------------------------------!
   
   !call h5dget_type_f(dset_id,datatype_id,hdferr)
   !call h5tget_size_f(datatype_id,setsize,hdferr)
   !call h5dclose_f(dset_id  , hdferr)
 
! =============================================================================================
! KEEP THIS CODE AS A TEMPLATE IN CASE WE NEED TO DO SOMETHING LIKE THIS IN THE FUTURE
!  HELPFUL IF WE CHANGE DATA TYPES
! ---------------------------------------------------------------------------------------------
!   if (setsize==4_8) then  !Old precision
!      call hdf_getslab_r(csite%soil_water,'SOIL_WATER_PA ',dsetrank,iparallel,.true.)
!   else if (setsize==8_8) then ! Newer precision
!      allocate(buff(nzg,csite%npatches))
!      write (unit=*,fmt='(a)') '-------------------------------------------------------------------'
!      write (unit=*,fmt='(a)') '  Loading 8-byte precision soil water and converting to 4-byte'
!      write (unit=*,fmt='(a)') '-------------------------------------------------------------------'
!      call hdf_getslab_d(buff,'SOIL_WATER_PA ',dsetrank,iparallel,.true.)
!      csite%soil_water(1:nzg,1:csite%npatches) = sngl(buff(1:nzg,1:csite%npatches))
!      deallocate(buff)
!  else
!     call fatal_error('Soil water dataset is not real nor double?'                         &
!                     ,'fill_history_site','ed_history_io.f90')
!  end if

  !--------------------------------------------------------------------------------------------  
  
  call hdf_getslab_r(csite%soil_tempk,'SOIL_TEMPK_PA ',dsetrank,iparallel,.true.)
  call hdf_getslab_r(csite%soil_fracliq,'SOIL_FRACLIQ_PA ',dsetrank,iparallel,.true.)

  dsetrank    = 2
  globdims(1) = int(n_pft,8)
  chnkdims(1) = int(n_pft,8)
  memdims(1)  = int(n_pft,8)
  memsize(1)  = int(n_pft,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8
  globdims(2) = int(npatches_global,8)
  chnkdims(2) = int(csite%npatches,8)
  chnkoffs(2) = int(sipa_index - 1,8)
  memdims(2)  = int(csite%npatches,8)
  memsize(2)  = int(csite%npatches,8)
  memoffs(2)  = 0_8
  
  call hdf_getslab_r(csite%A_o_max,'A_O_MAX ',dsetrank,iparallel,.true.) 
  call hdf_getslab_r(csite%A_c_max,'A_C_MAX ',dsetrank,iparallel,.true.) 
  call hdf_getslab_r(csite%repro,'REPRO_PA ',dsetrank,iparallel,.true.)


  dsetrank    = 2
  globdims(1) = int(n_dbh,8)
  chnkdims(1) = int(n_dbh,8)
  memdims(1)  = int(n_dbh,8)
  memsize(1)  = int(n_dbh,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8
  globdims(2) = int(npatches_global,8)
  chnkdims(2) = int(csite%npatches,8)
  chnkoffs(2) = int(sipa_index - 1,8)
  memdims(2)  = int(csite%npatches,8)
  memsize(2)  = int(csite%npatches,8)
  memoffs(2)  = 0_8
  call hdf_getslab_r(csite%co2budget_gpp_dbh,'CO2BUDGET_GPP_DBH ',dsetrank,iparallel,.true.)

!!!! MAY NEED TO ADD THIS ONE
!  call hdf_getslab_r(csite%old_stoma_data_max,'OLD ',dsetrank,iparallel,.true.)

  dsetrank    = 3
  globdims(3) = int(npatches_global,8)
  chnkdims(3) = int(csite%npatches,8)
  chnkoffs(3) = int(sipa_index - 1,8)

  memdims(3)  = int(csite%npatches,8)
  memsize(3)  = int(csite%npatches,8)
  memoffs(3)  = 0_8
  
  globdims(2) = int(ff_ndbh,8)
  chnkdims(2) = int(ff_ndbh,8)
  memdims(2)  = int(ff_ndbh,8)
  memsize(2)  = int(ff_ndbh,8)
  chnkoffs(2) = 0_8
  memoffs(2)  = 0_8

  globdims(1) = int(n_pft,8)
  chnkdims(1) = int(n_pft,8)
  memdims(1)  = int(n_pft,8)
  memsize(1)  = int(n_pft,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8

  call hdf_getslab_r(csite%pft_density_profile,'PFT_DENSITY_PROFILE ',dsetrank,iparallel,.true.)


  dsetrank    = 3
  globdims(3) = int(npatches_global,8)
  chnkdims(3) = int(csite%npatches,8)
  chnkoffs(3) = int(sipa_index - 1,8)

  memdims(3)  = int(csite%npatches,8)
  memsize(3)  = int(csite%npatches,8)
  memoffs(3)  = 0_8
  
  globdims(2) = int(n_pft,8)
  chnkdims(2) = int(n_pft,8)
  memdims(2)  = int(n_pft,8)
  memsize(2)  = int(n_pft,8)
  chnkoffs(2) = 0_8
  memoffs(2)  = 0_8

  globdims(1) = int(n_stoma_atts,8)
  chnkdims(1) = int(n_stoma_atts,8)
  memdims(1)  = int(n_stoma_atts,8)
  memsize(1)  = int(n_stoma_atts,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8

  call hdf_getslab_r(csite%old_stoma_vector_max,'OLD_STOMA_VECTOR_MAX ',dsetrank,iparallel,.true.)

  patchloop: do ipa=1,csite%npatches
     pftloop: do ipft = 1,n_pft
        csite%old_stoma_data_max(ipft,ipa)%recalc = int(csite%old_stoma_vector_max(1,ipft,ipa))
        csite%old_stoma_data_max(ipft,ipa)%T_L    = csite%old_stoma_vector_max(2,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%e_A    = csite%old_stoma_vector_max(3,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%PAR    = csite%old_stoma_vector_max(4,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%rb_factor = csite%old_stoma_vector_max(5,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%prss   = csite%old_stoma_vector_max(6,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%phenology_factor = csite%old_stoma_vector_max(7,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%gsw_open = csite%old_stoma_vector_max(8,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%ilimit = int(csite%old_stoma_vector_max(9,ipft,ipa))
        
        csite%old_stoma_data_max(ipft,ipa)%T_L_residual = csite%old_stoma_vector_max(10,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%e_a_residual = csite%old_stoma_vector_max(11,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%par_residual = csite%old_stoma_vector_max(12,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%rb_residual  = csite%old_stoma_vector_max(13,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%prss_residual= csite%old_stoma_vector_max(14,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%leaf_residual= csite%old_stoma_vector_max(15,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%gsw_residual = csite%old_stoma_vector_max(16,ipft,ipa)
     end do pftloop
  end do patchloop

  return
end subroutine fill_history_site

!==========================================================================================!
!==========================================================================================!

subroutine fill_history_patch(cpatch,paco_index,ncohorts_global,green_leaf_factor)
  
  use ed_state_vars,only: patchtype
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  use consts_coms, only: cliq,cice,t3ple,tsupercool
  use c34constants,only: n_stoma_atts
  use ed_max_dims,only: n_pft, n_mort
  use ed_therm_lib, only : calc_hcapveg
  use allometry, only : area_indices
  use therm_lib, only : qwtk
  implicit none

#if USE_INTERF
  interface
     subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_r
     subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_d
     subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_i
  end interface
#endif

  type(patchtype),target :: cpatch
  integer,intent(in) :: paco_index
  integer,intent(in) :: ncohorts_global
  real, dimension(n_pft), intent(in) :: green_leaf_factor
  integer :: iparallel,dsetrank
  
  ! Needed for reconstructing veg_energy if using an old restart
  ! ------------------------------------------------------------
  real :: plai
  integer :: ico
  ! ------------------------------------------------------------

  iparallel = 0
  
  dsetrank = 1
  globdims = 0_8
  chnkdims = 0_8
  chnkoffs = 0_8
  memoffs  = 0_8
  memdims  = 0_8
  memsize  = 1_8
  
  globdims(1) = int(ncohorts_global,8)
  chnkdims(1) = int(cpatch%ncohorts,8)
  chnkoffs(1) = int(paco_index - 1,8)

  memdims(1)  = int(cpatch%ncohorts,8)
  memsize(1)  = int(cpatch%ncohorts,8)
  memoffs(1)  = 0_8

  if(cpatch%ncohorts>0) then

     call hdf_getslab_i(cpatch%pft,'PFT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%nplant,'NPLANT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%hite,'HITE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%dbh,'DBH ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%agb,'AGB_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%basarea,'BA_CO',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%dagb_dt,'DAGB_DT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%dba_dt,'DBA_DT',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%ddbh_dt,'DDBH_DT',dsetrank,iparallel,.true.)

     call hdf_getslab_r(cpatch%bdead,'BDEAD ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%bleaf,'BLEAF ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%phenology_status,'PHENOLOGY_STATUS ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%balive,'BALIVE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%broot,'BROOT  ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%bsapwood,'BSAPWOOD ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%lai,'LAI_CO ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%llspan,'LLSPAN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%turnover_amp,'TURNOVER_AMP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%vm_bar,'VM_BAR ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%sla,'SLA ',dsetrank,iparallel,.true.)

     call hdf_getslab_r(cpatch%bstorage,'BSTORAGE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%cbr_bar,'CBR_BAR ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%veg_temp,'VEG_TEMP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%veg_water,'VEG_WATER ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%wpa,'WPA_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%wai,'WAI_CO ',dsetrank,iparallel,.true.)

     call hdf_getslab_r(cpatch%veg_energy,'VEG_ENERGY ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%hcapveg,'HCAPVEG ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%veg_fliq,'VEG_FLIQ ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%mean_gpp,'MEAN_GPP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%mean_leaf_resp,'MEAN_LEAF_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%mean_root_resp,'MEAN_ROOT_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_leaf_resp,'TODAY_LEAF_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_root_resp,'TODAY_ROOT_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_gpp,'TODAY_GPP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_gpp_pot,'TODAY_GPP_POT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_gpp_max,'TODAY_GPP_MAX ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%growth_respiration,'GROWTH_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%storage_respiration,'STORAGE_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%vleaf_respiration,'VLEAF_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%fsn,'FSN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%monthly_dndt,'MONTHLY_DNDT ',dsetrank,iparallel,.true.)

     if (associated(cpatch%mmean_gpp       )) &
          call hdf_getslab_r(cpatch%mmean_gpp,'MMEAN_GPP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_leaf_resp       )) &
          call hdf_getslab_r(cpatch%mmean_leaf_resp,'MMEAN_LEAF_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_root_resp       )) &
          call hdf_getslab_r(cpatch%mmean_root_resp,'MMEAN_ROOT_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_growth_resp       )) &
          call hdf_getslab_r(cpatch%mmean_growth_resp,'MMEAN_GROWTH_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_storage_resp       )) &
          call hdf_getslab_r(cpatch%mmean_storage_resp,'MMEAN_STORAGE_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_vleaf_resp       )) &
          call hdf_getslab_r(cpatch%mmean_vleaf_resp,'MMEAN_VLEAF_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_leaf_resp       )) &
          call hdf_getslab_r(cpatch%dmean_leaf_resp,'DMEAN_LEAF_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_root_resp       )) &
          call hdf_getslab_r(cpatch%dmean_root_resp,'DMEAN_ROOT_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_gpp       )) &
          call hdf_getslab_r(cpatch%dmean_gpp,'DMEAN_GPP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_fs_open       )) &
     call hdf_getslab_r(cpatch%dmean_fs_open,'DMEAN_FS_OPEN_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_fs_open       )) &
     call hdf_getslab_r(cpatch%mmean_fs_open,'MMEAN_FS_OPEN_CO ',dsetrank,iparallel,.false.) 
     if (associated(cpatch%dmean_fsw       )) &
     call hdf_getslab_r(cpatch%dmean_fsw,'DMEAN_FSW_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_fsw       )) &
     call hdf_getslab_r(cpatch%mmean_fsw,'MMEAN_FSW_CO ',dsetrank,iparallel,.false.) 
     if (associated(cpatch%dmean_fsn       )) &
     call hdf_getslab_r(cpatch%dmean_fsn,'DMEAN_FSN_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_fsn       )) &
     call hdf_getslab_r(cpatch%mmean_fsn,'MMEAN_FSN_CO ',dsetrank,iparallel,.false.) 
     if (associated(cpatch%mmean_leaf_maintenance )) &
     call hdf_getslab_r(cpatch%mmean_leaf_maintenance,'MMEAN_LEAF_MAINTENANCE ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%mmean_root_maintenance )) &
     call hdf_getslab_r(cpatch%mmean_root_maintenance,'MMEAN_ROOT_MAINTENANCE ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%mmean_leaf_drop       )) &
     call hdf_getslab_r(cpatch%mmean_leaf_drop,'MMEAN_LEAF_DROP_CO ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%mmean_cb       )) &
     call hdf_getslab_r(cpatch%mmean_cb,'MMEAN_CB ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%dmean_light_level       )) &
     call hdf_getslab_r(cpatch%dmean_light_level,'DMEAN_LIGHT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level       )) &
     call hdf_getslab_r(cpatch%mmean_light_level,'MMEAN_LIGHT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level       )) &
     call hdf_getslab_r(cpatch%dmean_light_level_beam,'DMEAN_LIGHT_LEVEL_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level_beam       )) &
     call hdf_getslab_r(cpatch%mmean_light_level_beam,'MMEAN_LIGHT_LEVEL_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_light_level_diff       )) &
     call hdf_getslab_r(cpatch%dmean_light_level_diff,'DMEAN_LIGHT_LEVEL_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level_diff       )) &
     call hdf_getslab_r(cpatch%mmean_light_level_diff,'MMEAN_LIGHT_LEVEL_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_par_v       )) &
     call hdf_getslab_r(cpatch%dmean_par_v,'DMEAN_PAR_V ',dsetrank,iparallel,.false.)

     if (associated(cpatch%dmean_par_v_beam       )) &
     call hdf_getslab_r(cpatch%dmean_par_v_beam,'DMEAN_PAR_V_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_par_v_diff       )) &
     call hdf_getslab_r(cpatch%dmean_par_v_diff,'DMEAN_PAR_V_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_par_v       )) &
     call hdf_getslab_r(cpatch%mmean_par_v,'MMEAN_PAR_V ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_par_v_beam       )) &
     call hdf_getslab_r(cpatch%mmean_par_v_beam,'MMEAN_PAR_V_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_par_v_diff       )) &
     call hdf_getslab_r(cpatch%mmean_par_v_diff,'MMEAN_PAR_V_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_beamext_level       )) &
     call hdf_getslab_r(cpatch%dmean_beamext_level,'DMEAN_BEAMEXT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_beamext_level      )) &
     call hdf_getslab_r(cpatch%mmean_beamext_level,'MMEAN_BEAMEXT_LEVEL ',dsetrank,iparallel,.false.)

     if (associated(cpatch%dmean_diffext_level       )) &
     call hdf_getslab_r(cpatch%dmean_diffext_level,'DMEAN_DIFFEXT_LEVEL ',dsetrank,iparallel,.false.)

     if (associated(cpatch%mmean_diffext_level       )) &
     call hdf_getslab_r(cpatch%mmean_diffext_level,'MMEAN_DIFFEXT_LEVEL ',dsetrank,iparallel,.false.)

     if (associated(cpatch%dmean_norm_par_beam      )) &
     call hdf_getslab_r(cpatch%dmean_norm_par_beam,'DMEAN_NORM_PAR_BEAM ',dsetrank,iparallel,.false.)

     if (associated(cpatch%mmean_norm_par_beam       )) &
     call hdf_getslab_r(cpatch%mmean_norm_par_beam,'MMEAN_NORM_PAR_BEAM ',dsetrank,iparallel,.false.)

     if (associated(cpatch%dmean_norm_par_diff       )) &
     call hdf_getslab_r(cpatch%dmean_norm_par_diff,'DMEAN_NORM_PAR_DIFF ',dsetrank,iparallel,.false.)

     if (associated(cpatch%mmean_norm_par_diff       )) &
     call hdf_getslab_r(cpatch%mmean_norm_par_diff,'MMEAN_NORM_PAR_DIFF ',dsetrank,iparallel,.false.)

     if (associated(cpatch%dmean_lambda_light       )) &
          call hdf_getslab_r(cpatch%dmean_lambda_light,'DMEAN_LAMBDA_LIGHT_CO ',dsetrank,iparallel,.false.)

     if (associated(cpatch%mmean_lambda_light       )) &
     call hdf_getslab_r(cpatch%mmean_lambda_light,'MMEAN_LAMBDA_LIGHT_CO ',dsetrank,iparallel,.false.)

     call hdf_getslab_r(cpatch%Psi_open,'PSI_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%krdepth,'KRDEPTH ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%first_census,'FIRST_CENSUS ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%new_recruit_flag,'NEW_RECRUIT_FLAG ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%light_level,'LIGHT_LEVEL ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%light_level_beam,'LIGHT_LEVEL_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%light_level_diff,'LIGHT_LEVEL_DIFF ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%beamext_level,'BEAMEXT_LEVEL ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%diffext_level,'DIFFEXT_LEVEL ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%norm_par_beam,'NORM_PAR_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%norm_par_diff,'NORM_PAR_DIFF ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%lambda_light,'LAMBDA_LIGHT_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%par_v,'PAR_V ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%par_v_beam,'PAR_V_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%par_v_diffuse,'PAR_V_DIFFUSE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rshort_v,'RSHORT_V ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rshort_v_beam,'RSHORT_V_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rshort_v_diffuse,'RSHORT_V_DIFFUSE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rlong_v,'RLONG_V ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rlong_v_surf,'RLONG_V_SURF ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rlong_v_incid,'RLONG_V_INCID ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rb,'RB ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%A_open,'A_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%A_closed,'A_CLOSED ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%Psi_closed,'PSI_CLOSED ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rsw_open,'RSW_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rsw_closed,'RSW_CLOSED ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%fsw,'FSW ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%fs_open,'FS_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%stomatal_resistance,'STOMATAL_RESISTANCE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%leaf_maintenance,'LEAF_MAINTENANCE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%root_maintenance,'ROOT_MAINTENANCE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%leaf_drop,'LEAF_DROP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%bseeds,'BSEEDS_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%leaf_respiration,'LEAF_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%root_respiration,'ROOT_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%gpp,'GPP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%paw_avg,'PAW_AVG ',dsetrank,iparallel,.true.)
     
     !----- 13-month dimension (12 previous months + current month). ----------------------!
     dsetrank    = 2
     globdims(1) = 13_8
     chnkdims(1) = 13_8
     chnkoffs(1) = 0_8
     memdims(1)  = 13_8
     memsize(1)  = 13_8
     memoffs(2)  = 0_8
     
     globdims(2) = int(ncohorts_global,8)
     chnkdims(2) = int(cpatch%ncohorts,8)
     chnkoffs(2) = int(paco_index - 1,8)
     
     memdims(2)  = int(cpatch%ncohorts,8)
     memsize(2)  = int(cpatch%ncohorts,8)
     memoffs(2)  = 0_8
     
     call hdf_getslab_r(cpatch%cb,'CB ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%cb_max,'CB_MAX ',dsetrank,iparallel,.true.)

     !----- 2-D, dimensioned by the number of mortality rates. ----------------------------!
     dsetrank    = 2
     globdims(1) = int(n_mort,8)
     chnkdims(1) = int(n_mort,8)
     chnkoffs(1) = 0_8
     memdims(1)  = int(n_mort,8)
     memsize(1)  = int(n_mort,8)
     memoffs(2)  = 0_8
     
     globdims(2) = int(ncohorts_global,8)
     chnkdims(2) = int(cpatch%ncohorts,8)
     chnkoffs(2) = int(paco_index - 1,8)
     
     memdims(2)  = int(cpatch%ncohorts,8)
     memsize(2)  = int(cpatch%ncohorts,8)
     memoffs(2)  = 0_8

     call hdf_getslab_r(cpatch%mort_rate,'MORT_RATE_CO ',dsetrank,iparallel,.true.)


     dsetrank    = 2
     globdims(1) = int(n_stoma_atts,8)
     chnkdims(1) = int(n_stoma_atts,8)
     chnkoffs(1) = 0_8
     memdims(1)  = int(n_stoma_atts,8)
     memsize(1)  = int(n_stoma_atts,8)
     memoffs(2)  = 0_8
     
     globdims(2) = int(ncohorts_global,8)
     chnkdims(2) = int(cpatch%ncohorts,8)
     chnkoffs(2) = int(paco_index - 1,8)
     
     memdims(2)  = int(cpatch%ncohorts,8)
     memsize(2)  = int(cpatch%ncohorts,8)
     memoffs(2)  = 0_8
     

     call hdf_getslab_r(cpatch%old_stoma_vector,'OLD_STOMA_VECTOR', &
          dsetrank,iparallel,.true.)

  cohortloop: do ico=1,cpatch%ncohorts
     cpatch%old_stoma_data(ico)%recalc = int(cpatch%old_stoma_vector(1,ico))
     cpatch%old_stoma_data(ico)%T_L    = cpatch%old_stoma_vector(2,ico)
     cpatch%old_stoma_data(ico)%e_A    = cpatch%old_stoma_vector(3,ico)
     cpatch%old_stoma_data(ico)%PAR    = cpatch%old_stoma_vector(4,ico)
     cpatch%old_stoma_data(ico)%rb_factor = cpatch%old_stoma_vector(5,ico)
     cpatch%old_stoma_data(ico)%prss = cpatch%old_stoma_vector(6,ico) 
     cpatch%old_stoma_data(ico)%phenology_factor = cpatch%old_stoma_vector(7,ico)
     cpatch%old_stoma_data(ico)%gsw_open = cpatch%old_stoma_vector(8,ico)
     cpatch%old_stoma_data(ico)%ilimit   = int(cpatch%old_stoma_vector(9,ico))
     cpatch%old_stoma_data(ico)%T_L_residual = cpatch%old_stoma_vector(10,ico)
     cpatch%old_stoma_data(ico)%e_a_residual = cpatch%old_stoma_vector(11,ico)
     cpatch%old_stoma_data(ico)%par_residual = cpatch%old_stoma_vector(12,ico)
     cpatch%old_stoma_data(ico)%rb_residual  = cpatch%old_stoma_vector(13,ico)
     cpatch%old_stoma_data(ico)%prss_residual= cpatch%old_stoma_vector(14,ico) 
     cpatch%old_stoma_data(ico)%leaf_residual= cpatch%old_stoma_vector(15,ico)
     cpatch%old_stoma_data(ico)%gsw_residual = cpatch%old_stoma_vector(16,ico)
  enddo cohortloop

endif


  return
end subroutine fill_history_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
  
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  
  
  implicit none
  
  real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff

  integer :: hdferr,dsetrank
  integer :: iparallel
  character(len=*),intent(in) :: varn
  
  ! Some datasets are required during model initialization, such as the state
  ! variables like soil, canopy and vegetation energy and mass, although some variables
  ! are involved only in diagnostics.  If they are missing from the restart dataset,
  ! that may compromise diagnostic variables that are being averaged over the current time
  ! period, but they will not compromise the transition of the prognostic model state from one
  ! simulation to the next.
  ! If the dataset is not required, pass it in as a .true. argument

  logical,intent(in) :: required

  ! If the the optional argument is not present, take the conservative stance
  ! and make sure that it is "not-not-not-not-required", "not-not-required", or "required".

  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0 .and. required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
          ,'hdf_getslab_r','ed_history_io.f90')
     
  else if (hdferr /= 0 .and. .not.required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
     write (unit=*,fmt='(a)') '                                                           '
     write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
     write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
     write (unit=*,fmt='(a)') ' + his may cause some of your diagnostic output related'
     write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
     write (unit=*,fmt='(a)') ''
     write (unit=*,fmt='(a)') '   This variable has been specified as:'
     write (unit=*,fmt='(a)') '   NOT ABSOUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') ''
     
     buff=0.
     return
     
  else
     
     call h5dget_space_f(dset_id,filespace,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
          chnkdims,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
     if (hdferr /= 0) then
        write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
        write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
        call fatal_error('Could not create the hyperslabs memspace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
          memdims,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     if (iparallel == 1) then
        
        call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace, &
             xfer_prp = plist_id)
        
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if

     else

        call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace )

        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if

     endif
     
     !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
     
     call h5sclose_f(filespace, hdferr)
     call h5sclose_f(memspace , hdferr)
     call h5dclose_f(dset_id  , hdferr)
     
  endif
  
  
  return
end subroutine hdf_getslab_r

!==========================================================================================!
!==========================================================================================!

subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
  
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  
  
  implicit none
  
  real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff

  integer :: hdferr,dsetrank
  integer :: iparallel
  character(len=*),intent(in) :: varn

  ! Some datasets are required during model initialization, such as the state
  ! variables like soil, canopy and vegetation energy and mass, although some variables
  ! are involved only in diagnostics.  If they are missing from the restart dataset,
  ! that may compromise diagnostic variables that are being averaged over the current time
  ! period, but they will not compromise the transition of the prognostic model state from one
  ! simulation to the next.
  ! If the dataset is not required, pass it in as a .true. argument

  logical,intent(in)  :: required


  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0 .and. required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
          ,'hdf_getslab_d','ed_history_io.f90')
     
  else if (hdferr /= 0 .and. .not.required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
     write (unit=*,fmt='(a)') '                                                           '
     write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
     write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
     write (unit=*,fmt='(a)') ' + This may cause some of your diagnostic output related'
     write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
     write (unit=*,fmt='(a)') ''
     write (unit=*,fmt='(a)') '   This variable has been specified as:'
     write (unit=*,fmt='(a)') '   NOT ABSOLUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') ''
     
     buff=0.d0
     return
     
  else
  
     call h5dget_space_f(dset_id,filespace,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
          chnkdims,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
     if (hdferr /= 0) then
        write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
        write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
        call fatal_error('Could not create the hyperslabs memspace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
          memdims,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_r','ed_history_io.f90')
     end if
     
     if (iparallel == 1) then
        
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace, &
             xfer_prp = plist_id)
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if
        
     else
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace )
        
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if
     endif
     
     !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
     
     call h5sclose_f(filespace, hdferr)
     call h5sclose_f(memspace , hdferr)
     call h5dclose_f(dset_id  , hdferr)

  endif
  

  return
end subroutine hdf_getslab_d
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
  
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  
  implicit none
  
  integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff

  integer :: hdferr,dsetrank
  integer :: iparallel
  character(len=*),intent(in) :: varn
  ! Some datasets are required during model initialization, such as the state
  ! variables like soil, canopy and vegetation energy and mass, although some variables
  ! are involved only in diagnostics.  If they are missing from the restart dataset,
  ! that may compromise diagnostic variables that are being averaged over the current time
  ! period, but they will not compromise the transition of the prognostic model state from one
  ! simulation to the next.
  ! If the dataset is not required, pass it in as a .true. argument

  logical,intent(in) :: required
  
  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0 .and. required ) then
     
     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
          ,'hdf_getslab_i','ed_history_io.f90')
     
  else if (hdferr /= 0 .and. .not.required ) then
     
     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
     write (unit=*,fmt='(a)') '                                                           '
     write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
     write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
     write (unit=*,fmt='(a)') ' + This may cause some of your diagnostic output related'
     write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
     write (unit=*,fmt='(a)') ''
     write (unit=*,fmt='(a)') '   This variable has been specified as:'
     write (unit=*,fmt='(a)') '   NOT ABSOLUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') ''
     
     buff=0
     return
     
  else
  
     call h5dget_space_f(dset_id,filespace,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_i','ed_history_io.f90')
     endif
     
     call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
          chnkdims,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_i','ed_history_io.f90')
     endif
     
     call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
     if (hdferr /= 0) then
        write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
        write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
        call fatal_error('Could not create the hyperslabs memspace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_i','ed_history_io.f90')
     endif
     
     call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
          memdims,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
             ,'hdf_getslab_i','ed_history_io.f90')
     end if

     if (iparallel == 1) then
        
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace, &
             xfer_prp = plist_id)
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_i','ed_history_io.f90')
        end if
        
     else
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace )
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if
     end if
     
     !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
     
     call h5sclose_f(filespace, hdferr)
     call h5sclose_f(memspace, hdferr)
     call h5dclose_f(dset_id, hdferr)

  endif
     
  return
end subroutine hdf_getslab_i
!==========================================================================================!
!==========================================================================================!
