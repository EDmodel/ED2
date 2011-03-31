!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine rams_node()

  use mem_grid, only : &
       ideltat,        & ! INTENT(IN)
       iflag,          & ! INTENT(IN) - ALF - Modified in DTSET
       maxsched,       & ! INTENT(IN)
       maxschent,      & ! INTENT(IN)
       nndtrat,        & ! INTENT(IN) - Modified in DTSET
       time,           & ! INTENT(INOUT)
       timmax,         & ! INTENT(IN)
       dtlongn,        & ! INTENT(IN) - a ser Modified in DTSET
       ngrids,         & ! INTENT(IN)
       nsubs,          & ! INTENT(IN) - a ser (OUT) na chamada a modsched
       isstp,          & ! INTENT(OUT)
       isched,         & ! INTENT(IN) - a ser (OUT) na chamada a modsched
       ngrid,          & ! INTENT(OUT)
       dtlt,           & ! INTENT(IN)
       ngbegun,        & ! INTENT(OUT)
       nxtnest,        & ! INTENT(IN)
       nzg,            & ! INTENT(IN)
       nxp, nyp,       & ! INTENT(IN)
       f_thermo_e,     & ! INTENT(OUT)
       f_thermo_w,     & ! INTENT(OUT)
       f_thermo_s,     & ! INTENT(OUT)
       f_thermo_n        ! INTENT(OUT)


  use node_mod, only : &
       ipara,          & ! INTENT(OUT)
       master_num,     & ! INTENT(IN)
       load_bal,       & ! INTENT(IN)
       mzp,            & ! INTENT(IN)
       mxp,            & ! INTENT(IN)
       myp,            & ! INTENT(IN)
       mi0,            & ! INTENT(IN)
       mj0,            & ! INTENT(IN)
       mynum             ! INTENT(IN)

  use io_params, only : & ! Include by Alvaro L.Fazenda
       maxgrds,         & ! INTENT(IN)
       avgtim,          & !INTENT(IN)
       frqmean,         & !INTENT(IN)
       frqboth            !INTENT(IN)

  ! Needed for CATT
  use catt_start, only: CATT           ! intent(in)

  ! Only needed if checking (debugging) reprodutibility
  !use node_mod, only : nmachs ! INTENT(IN)

  ! ALF
  ! Necessary in new advection scheme
  use advect_kit, only :   &
       advect_first_alloc, &  ! Subroutine
       prepare_inv            ! Subroutine

  ! ALF
  ! Necessary for allocation
  use node_mod, only : &
       ibcon,          & ! INTENT(IN)
       ia, iz, izu,    & ! INTENT(IN)
       ja, jz, jzv,    & ! INTENT(IN)
       i0, j0,         & ! INTENT(IN)
       mmxp,           & ! INTENT(IN)
       mmyp,           & ! INTENT(IN)
       mmzp,           & ! INTENT(IN)
       nsend_buff,     & ! intent(out)
       nrecv_buff      ! ! intent(out)
  
  use mem_leaf, only: isfcl ! intent(in)

  use dtset, only: dtset_new ! subroutine

  implicit none

  ! Local Variables:
  include 'interface.h'
  include 'mpif.h'
  integer :: isendflg,isendlite,isendmean,isendboth,nt,npass,icm,ifm,nfeed
  real :: wstart,totcpu,t1,w1,t6,w6
  real(kind=8) :: begtime
  real, external :: walltime
  integer :: nndtflg ! ALF - For local processing
  integer :: isendbackflg ! ALF - For local processing
  !ALF
  real :: dxtmax_local(maxgrds)
  integer :: ng
  !MLO
  integer :: ierr ! For new MPI call style

  ipara=1
  
  ! Reset the send and receive sizes
  nsend_buff(:) = 0
  nrecv_buff(:) = 0
  
  !          Call routine to initialize input parameters
  !               and namelist settings
  !          -----------------------------------

  call init_params(1)

  !          Allocate memory needed on node

  call rams_mem_alloc(2)


  !          Routine to get fields from master and finish initialization
  !          -----------------------------------------------------------

  call init_fields(.true.)

  isendflg=0
  isendlite = 0
  isendmean = 0
  isendboth = 0
  isendbackflg = 0 ! ALF

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !------------------------------------------------
  !     Loop through total number of long timesteps
  !------------------------------------------------

  wstart=walltime(0.)
  nt = 0

  ! ALF
  ! Initialize dtlongn, nndtrat, and nnacoust, and compute the timestep
  ! schedule for all grid operations, locally.
  call node_getdtsched(isendflg,isendlite,isendmean,isendboth)

  ! ALF
  ! Receiving DXTMAX in local domain
  call node_getdxt(dxtmax_local)


  ! ALF
  ! Preparing data for Advection Scheme
  ! Memory allocation for new advection scheme
  call advect_first_alloc(ngrids, mmzp(1:ngrids), mmxp(1:ngrids), &
       mmyp(1:ngrids))
  ! Invariable data
  call prepare_inv(ngrids)

  ! ALF
  ! Checking if the actual node have to run THERMO on the boundaries
  ! Loop through grids
  do ng=1,ngrids
     ! Setting grid properties
     call newgrid(ng)
     call node_index()
     ! Checking Eastern boundary
     f_thermo_e(ng) = i0 == 0
     ! Checking Weastern boundary
     f_thermo_w(ng) = (mxp+i0) == nxp
     ! Checking Southern boundary
     f_thermo_s(ng) = j0 == 0
     ! Checking Northern boundary
     f_thermo_n(ng) = (myp+j0) == nyp
  enddo

  !----- Initialise microphysics tables ---------------------------------------------------!
  call micro_1st()
  !----------------------------------------------------------------------------------------!

  do while (time<timmax)

     totcpu=0
     w1=walltime(wstart)
     call timing(1,t1)

     nt = nt + 1
     begtime=time

!!!! This part is for the dynamic balancing and sending new varfile/sst info

     if (isendflg==1) then

        if (load_bal==1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Dynamic balance experimental for now until OS memory
           ! issues can be solved
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call nodeget_grid_dimens()
           !CALL nodeget_grid_dimens(ngrids) ! Modif.by Craig
           call rams_mem_alloc(3)  ! Dealloc, then allocate new vars
        endif

     endif

     if (isendbackflg==1) then

        call init_fields(.false.) ! Get fields from master

     endif


     !    Receive message from master containing ISENDFLG, new dt's, etc.

     ! ALF
     ! Determines, locally, whether send stuff back at the END of the
     ! timestep
     call comm_time(isendflg,isendlite,isendmean,isendboth,isendbackflg)


     ! ALF
     ! Examine Courant numbers in case model needs to be stopped
     ! or (if ideltat < 0), to update dtlongn, nndtrat,
     ! nnacoust, sspct and isched.
     !
     !call dtset_local(nndtflg, dxtmax_local)
     call dtset_new(mynum, nndtflg, dxtmax_local)

     if (iflag>0) then  !Estourou Courant - modelo ira´ parar. Antes salva tudo
        isendflg = 1
        isendlite = 1
        isendmean = 1
        isendboth = 1
     endif

     if (nndtflg>0) then
        call modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
     endif
     !

     !------------------------------------------------------------------------
     !                  Loop through all grids and advance a 'DTLONG' timestep.
     !-------------------------------------------------------------------------

     !                  Start the timestep schedule

     do npass=1,nsubs

        isstp=isched(npass,3)
        ngrid=isched(npass,1)
        call newgrid(ngrid)

        call node_index()

        !---------------------------------------------------------------------
        !         Advance this grid forward by the appropriate timestep.

        time=begtime + (isched(npass,5)-1) * dtlt

        !      Call main timestep driver
        !      ------------------------------
        call timestep()

        ngbegun(ngrid)=1

        !---------------------------------------------------------------------
        !---------------------------------------------------------------------
        !---------------------------------------------------------------------
        !    Is it time to send the coarse domain points to the nested
        !    nodes to interpolate a nested grid's boundaries?

        if(isched(npass,2) /= 0) then
           icm = ngrid
           do ifm = 1, ngrids
              if (nxtnest(ifm) == icm) then
                 ngrid = ifm            ! JP: is this necessary?
                 isstp=isched(npass,3)  ! JP: is this necessary (maybe incorrect)!
                 call newgrid(ifm)
                 call node_sendnbc(ifm, icm)
                 call node_getnbc(ifm, icm)
              end if
           end do
        endif
        !---------------------------------------------------------------------
        !---------------------------------------------------------------------
        !---------------------------------------------------------------------
        !                  Is it time to feedback fields from fine grids to
        !                     coarser grids?

        if (isched(npass,4)/=0) then
           ngrid=isched(npass,1)
           do nfeed=1,isched(npass,4)

              call newgrid(ngrid)

              call node_sendfeed(ngrid)

              call node_getfeed(nxtnest(ngrid),ngrid)
              ngrid=nxtnest(ngrid)

           enddo

        endif
        !---------------------------------------------------------------------

     enddo

     !------------------------------------------------------------------------
     !        Also, average each of the analysis variables over time
     do ngrid=1,ngrids
        call newgrid(ngrid)

        !          THETAM and RVM have not been updated after nesting feedback
        !             This means that these variables are really a timestep
        !             behind the instantaneous variables.

        !          Calculate the means
        ! Mod. by Alvaro L. Fazenda
        if ((avgtim /= 0.).and.(frqmean /= 0. .or. frqboth /= 0.))  &
             call anlavg(mzp,mxp,myp,nzg)
     enddo

     !------------------------------------------------------------------------
     ! Send timing info/CFL numbers back to master.

     do ifm = 1,ngrids
        call newgrid(ifm)
        call cfl(mzp,mxp,myp,mi0(ifm),mj0(ifm),mynum)
     enddo

     call timing(2,t6)
     w6=walltime(wstart)
     totcpu=totcpu+t6-t1

     call node_putcflcpu(t6-t1,w6-w1)

     ! Receiveing CFL to Recalculate DeltaT if necessary
     if (ideltat < 0) then
        call node_getcflmax()
     endif

     !------------------------------------------------------------------------
     ! Send entire subdomains back to master every now and then.

     if (isendflg==1) then
        call node_sendall()
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     if (isendlite==1) then
        call node_sendanl('LITE')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     if (isendmean==1) then
        call node_sendanl('MEAN')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     if (isendboth==1) then
        call node_sendanl('BOTH')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     !------------------------------------------------------------------------
     !                   Update main time variable by a long timestep.

     time=begtime+dtlongn(1)

  enddo

end subroutine rams_node
!
!     ****************************************************************
!
subroutine init_params(init)

  use mem_grid
  use node_mod
  use mem_oda
  use mem_radiate, only: ISWRTYP, ILWRTYP ! Intent(in)
  use mem_leaf   , only: isfcl ! Intent(in)

  implicit none

  integer, intent(in) :: init

  include 'mpif.h'
  include 'interface.h'
  integer :: ierr

  !          get all initialization info from the master process
  !          -------------------------------------------------------

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call nodeget_processid(init)
  call nodeget_nl
  if (isfcl == 5) call nodeget_ednl(master_num)
  call nodeget_gridinit
  if (ibnd .eq. 4 .or. jbnd .eq. 4) then
     call ipaths_cyc_alloc(nnxp(1),nnyp(1),ibnd,jbnd)
  endif
  call nodeget_grid_dimens()
  call nodeget_gridset
  call nodeget_cofnest()
  call nodeget_micphys
  if (if_oda == 1) call nodeget_oda()
  call nodeget_misc
end subroutine init_params
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_fields(init)

   use mem_grid
   use node_mod

   use var_tables

   use mem_leaf  , only : isfcl
   use mem_cuparm, only : nclouds
   use mem_aerad , only : nwave
   use grid_dims , only : ndim_types

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical, intent(in)            :: init
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(ndim_types) :: npvar
   integer                        :: ierr
   integer                        :: hugedim
   integer                        :: ng
   integer                        :: nm
   integer                        :: itype
   integer                        :: i1
   integer                        :: j1
   integer                        :: i2
   integer                        :: j2
   integer                        :: xlbc
   integer                        :: ylbc
   integer                        :: fdzp
   integer                        :: fdep
   integer                        :: idim
   integer                        :: memf
   integer                        :: nv
   !----- Include modules. ----------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !      Initialise surface constants.                                                    !
   !---------------------------------------------------------------------------------------!
   if (init) then
      select case (isfcl)
      case (1,2,5)
         call sfcdata
      case (3)
         !call sfcdata_sib_driver
      end select
   end if
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !          Get all necessary fields from master.                                        !
   !---------------------------------------------------------------------------------------!
   call node_getinit()
   if (isfcl == 5 .and. init) then
      call node_ed_init()
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !       Check feedback buffer.                                                          !
   !---------------------------------------------------------------------------------------!
   ! itype=6 !Changed for reproducibility - Saulo Barros
   itype=7
   !---------------------------------------------------------------------------------------!



   !----- Find number of lbc variables to be communicated. --------------------------------!
   npvar(:) = 0
   do nv = 1,num_var(1)
      if (vtab_r(nv,1)%impt1 == 1 ) then
         idim        = vtab_r(nv,1)%idim_type
         npvar(idim) = npvar(idim) + 1
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find the size of the lateral boundary condition buffers.                         !
   !---------------------------------------------------------------------------------------!
   nbuff_feed=0
   do ng=1,ngrids
      do nm=1,nmachs
         i1         = ipaths(1,itype,ng,nm)
         i2         = ipaths(2,itype,ng,nm)
         j1         = ipaths(3,itype,ng,nm)
         j2         = ipaths(4,itype,ng,nm)
         xlbc       = i2 - i1 + 1
         ylbc       = j2 - j1 + 1 
         
         memf = 0
         do idim = 2, ndim_types
            call ze_dims(ng,idim,.false.,fdzp,fdep)
            memf = memf + fdzp * xlbc * ylbc * fdep * npvar(idim)
         end do 

         nbuff_feed = max(nbuff_feed,memf)
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Allocate long time step send and receive buffers.  Check the size of real and      !
   ! integer packages, and set up the largest as standard, so nothing is lost...           !
   !---------------------------------------------------------------------------------------!
   call MPI_Pack_size(1,MPI_REAL,MPI_COMM_WORLD,mpi_real_size,ierr)
   call MPI_Pack_size(1,MPI_INTEGER,MPI_COMM_WORLD,mpi_int_size,ierr)
   f_ndmd_size = max(mpi_real_size,mpi_int_size)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     De-allocate the buffers if this is not the first time this sub-routine is called. !
   !---------------------------------------------------------------------------------------!
   if (.not. init) then
      do nm=1,nmachs
         call dealloc_node_buff(node_buffs_lbc(nm))
         call dealloc_node_buff(node_buffs_st(nm))
         call dealloc_node_buff(node_buffs_feed(nm))
         call dealloc_node_buff(node_buffs_nest(nm))
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If the buffer is smaller than the actual needs, then de-allocate and re-allocate  !
   ! the buffers.                                                                          !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (nsend_buff(nm) > 0 .or. nrecv_buff(nm) > 0) then
         nbuff_feed = max(nbuff_feed,nsend_buff(nm),nrecv_buff(nm))
         call dealloc_node_buff(node_buffs_lbc(nm))
         call alloc_node_buff(node_buffs_lbc(nm),nbuff_feed,f_ndmd_size)
         call dealloc_node_buff(node_buffs_st(nm))
         call alloc_node_buff(node_buffs_st(nm),nbuff_feed,f_ndmd_size)
         call dealloc_node_buff(node_buffs_feed(nm))
         call alloc_node_buff(node_buffs_feed(nm),nbuff_feed,f_ndmd_size)
         call dealloc_node_buff(node_buffs_nest(nm))
         call alloc_node_buff(node_buffs_nest(nm),nbuff_feed,f_ndmd_size)
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      In case the boundary conditions are cyclic, initialise the cyclic structure.     !
   !---------------------------------------------------------------------------------------!
   if (ibnd == 4 .or. jbnd == 4) then
      call node_cycinit(nnzp(1),nnxp(1),nnyp(1),npvar,nmachs,ibnd,jbnd,mynum)
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_fields
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine node_index()

  use node_mod

  implicit none

  ia_1=max(ia-1,1)
  ia_2=max(ia-2,1)
  ia_3=max(ia-3,1)
  ia1=ia+1
  ia2=ia+2
  ia3=ia+3
  iz_1=iz-1
  iz_2=iz-2
  iz_3=iz-3
  iz1=min(iz+1,mxp)
  iz2=min(iz+2,mxp)
  iz3=min(iz+3,mxp)

  izu=iz
  if(iand(ibcon,2) /= 0) izu=iz-1

  if(myp > 1) then
     ja_1=max(ja-1,1)
     ja_2=max(ja-2,1)
     ja_3=max(ja-3,1)
     ja1=ja+1
     ja2=ja+2
     ja3=ja+3
     jz_1=jz-1
     jz_2=jz-2
     jz_3=jz-3
     jz1=min(jz+1,myp)
     jz2=min(jz+2,myp)
     jz3=min(jz+3,myp)

     jzv=jz
     if(iand(ibcon,8) /= 0) jzv=jz-1

  else
     print*,'Trying to do 2-dimensional run ??????'
     stop 'no parallel 2d'
     ja_1=1
     ja_2=1
     ja_3=1
     ja1=1
     ja2=1
     ja3=1
     jz_1=1
     jz_2=1
     jz_3=1
     jz1=1
     jz2=1
     jz3=1
     jzv=1
  endif

  return
end subroutine node_index

! For reproducibility - Saulo Barros
subroutine print_field_node(mzp,mxp,myp,a,i0,j0,ia,iz,ja,jz,mynum,lev)
  implicit none
  integer :: mzp, mxp, myp, i0, j0, ia, iz, ja, jz, mynum, lev
  real :: a(mzp,mxp,myp)

  integer :: iu, i, j
  !
  iu = 70 + mynum
  write(iu,*) ' new field '
  write(iu,*) 1+j0,myp+j0,1+i0,mxp+i0
  write(iu,*) ja+j0,jz+j0
  do j=ja,jz
     write(iu,*) ia+i0,iz+i0
  enddo
  do j=1,myp
     write(iu,*) (a(lev,i,j),i=1,mxp)
  enddo
  return
end subroutine print_field_node
! -----------------------------------
