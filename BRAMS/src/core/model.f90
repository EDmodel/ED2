!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine model()

  use grid_dims, only: &
       maxgrds

  use mem_grid, only: &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       nxp,           &
       nyp,           &
       nzp,           &
       nzg,           &
       time,          &
       timmax,        &
       istp,          &
       isstp,         &
       ngbegun,       &
       isched,        &
       maxsched,      &
       maxschent,     &
       ngrid,         &
       ngrids,        &
       nxtnest,       &
       nndtrat,       &
       nsubs,         &
       dtlt,          &
       dtlongn,       &
       f_thermo_e,    & ! INTENT(OUT)
       f_thermo_w,    & ! INTENT(OUT)
       f_thermo_s,    & ! INTENT(OUT)
       f_thermo_n       ! INTENT(OUT)

  ! io_params include by Alvaro L.Fazenda
  use io_params, only : &
       avgtim,          & ! INTENT(IN)
       frqmean,         & ! INTENT(IN)
       frqboth            ! INTENT(IN)

  ! Needed for CATT
  use catt_start, only: CATT ! intent(in)

  ! ALF
  ! Necessary in new advection scheme
  use advect_kit, only :   &
       advect_first_alloc, &  ! Subroutine
       prepare_inv            ! Subroutine

  use node_mod, only : ia, iz, izu, ja, jz, jzv ! INTENT(IN)

  use mem_leaf, only : isfcl

  use dtset, only: dtset_new ! subroutine

  implicit none

  !   +------------------------------------------------------------------
  !   ! This routine drives the entire time integration process
  !   !   for a non-parallel run.
  !   +------------------------------------------------------------------

  ! Local Variables:

  integer :: npass,nndtflg,icm,ifm,nfeed,mynum,i
  real :: wtime_start,t1,wtime1,wtime2,t2,wtime_tot
  real, external :: walltime
  character(len=10) :: c0, c1, c2
  character(len=*), parameter :: h="**(model)**"

  real(kind=8) :: begtime

  !ALF
  real :: dxtmax_local(maxgrds)

  wtime_start=walltime(0.)
  istp = 0

  ! ALF
  ! Preparing data for Advection Scheme
  ! Memory allocation for new advection scheme
  call onenode()
  call advect_first_alloc(ngrids, nnzp(1:ngrids), nnxp(1:ngrids), &
       nnyp(1:ngrids))
  ! Invariable data
  call prepare_inv(ngrids)

  ! ALF
  ! Setting the code to call Thermo on 4 boundaries in Serial runs
  f_thermo_e = .true.
  f_thermo_w = .true.
  f_thermo_n = .true.
  f_thermo_s = .true.

  !----- Initialise microphysics tables ---------------------------------------------------!
  call micro_1st()
  !----------------------------------------------------------------------------------------!

  !         Start the timesteps
  write(*,"(/,a,/)") " === Time integration starts (model) ==="
  do while (time .lt. timmax)

     istp = istp + 1
     begtime=time

     !            CPU timing information

     call timing(1,t1)
     wtime1=walltime(wtime_start)

     ! Examine Courant numbers in case model needs to be stopped
     ! or (if ideltat < 0), to update dtlongn, nndtrat,
     ! nnacoust, sspct and isched.

     !call dtset(nndtflg)
     call dtset_new(0, nndtflg, dxtmax_local)


     if (nndtflg .gt. 0) then
        call dump_dtset(nndtflg)
        call modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
        call dump_modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
     endif

     ! Start the timestep schedule to loop through all grids and advance them
     ! in time an increment equal to dtlongn(1).

     do npass=1,nsubs

        isstp=isched(npass,3)
        ngrid=isched(npass,1)
        call newgrid(ngrid)

        time=begtime + (isched(npass,5)-1) * dtlt

        call onenode()

        ! run a timestep at grid ngrid

        call timestep()

        ! if scheduled, send boundary conditions to 
        ! all direct son grids

        ngbegun(ngrid)=1
        if(isched(npass,2) /= 0) then
           icm = ngrid
           do ifm = 1, ngrids
              if (nxtnest(ifm) == icm) then
                 call newgrid(ifm)
                 call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm)  &
                      ,nnzp(icm),nnxp(icm),nnyp(icm),0,0,ifm,0,mynum)
              end if
           end do
        endif

        ! if scheduled, send feedback fields to 
        ! all parent grids

        if(isched(npass,4).ne.0) then
           ifm = isched(npass,1)
           icm = nxtnest(ifm)
           do nfeed = 1,isched(npass,4)
              call newgrid(ifm)
              call nstfeed(ifm,icm)
              ifm = icm
              icm = nxtnest(ifm)
           enddo
        endif

     enddo

     ! At this point, all grids have been advanced forward by time increment DTLONG,
     ! and all nesting operations (interpolation and feedback) have been carried
     ! out.  If there are two hemispheric grids in this simulation, now is the
     ! time to carry out the communication between them.  Subroutine hemintrp
     ! will return immediately if nhemgrd2 is not greater than 1.

     call hemintrp

     ! Compute Courant numbers cflxy and cflz and do averaging.

     do ifm = 1,ngrids
        call newgrid(ifm)
        call cfl(nzp,nxp,nyp,0,0,1)

        !           THETAM and RVM have not been updated after nesting feedback
        !              This means that these variables are really a timestep
        !              behind the instantaneous variables.

        !           Calculate the means

        ! Mod. by Alvaro L. Fazenda
        if((avgtim /= 0.).and.(frqmean /= 0. .or. frqboth /= 0.))  &
             call anlavg(nzp,nxp,nyp,nzg)

     enddo

     ! New position to update 'TIME' - ALF
     time=begtime+dtlongn(1)

     wtime2=walltime(wtime_start)
     call TIMING(2,T2)

     write(*,"(a,i6,a,f9.1,2(a,f6.2),a)") &
          " Timestep ",istp,&
          "; Sim time",time,&
          "s; Wall",wtime2-wtime1,&
          "s; CPU",t2-t1,&
          "s"

     call rams_output

  enddo

  wtime_tot=walltime(wtime_start)
  write(c0,"(f10.1)") wtime_tot
  write(*,"(/,a,/)") " === Time integration ends; Total elapsed time="//&
       &trim(adjustl(c0))//" ===" 
end subroutine MODEL

!     *****************************************************************

subroutine par_model(master_num)

  use grid_dims, only: &
       maxgrds

  use mem_grid, only : &
       ideltat,        & ! INTENT(IN)
       time,           & ! INTENT(INOUT)
       timmax,         & ! INTENT(IN)
       istp,           & ! INTENT(OUT)
       iflag,          & ! INTENT(IN) - Modifyed in DTSET
       isched,         & ! INTENT(INOUT)
       maxsched,       & ! INTENT(IN) - Modifyed in MODSCHED
       maxschent,      & ! INTENT(IN) - Modifyed in MODSCHED
       ngrids,         & ! INTENT(IN)
       nxtnest,        & ! INTENT(IN)
       nndtrat,        & ! INTENT(IN) - Modifyed in DTSET
       nsubs,          & ! INTENT(OUT)
       dtlongn           ! INTENT(IN) - Modifyed in DTSET

  use rpara, only : &
       load_bal,       & ! INTENT(IN)
       nmachs,         & ! INTENT(IN)
       machnum,        & ! INTENT(IN)
       ptimes            ! INTENT(IN)

  use dtset, only: dtset_new ! subroutine

  implicit none
  !   +------------------------------------------------------------------
  !   ! This routine drives the entire time integration process
  !   !   for a parallel run.
  !   +------------------------------------------------------------------

  ! Local Variables:

  integer :: isendflg,isendlite,isendmean,isendboth,nndtflg,ntsend,nmach,n,i
  real :: wtime_start,t1,wtime1,wtime2,t2,wtime_tot,pcpu,pwall
  real, external :: walltime
  integer :: isendbackflg ! ALF - For local processing
  character(len=10) :: c0
  character(len=*), parameter :: h="**(par_model)**"

  !MLO
  integer :: ierr, master_num
  include 'mpif.h'

  !ALF
  real :: dxtmax_local(maxgrds)

  wtime_start=walltime(0.)

  isendflg=0
  isendlite = 0
  isendmean = 0
  isendboth = 0
  isendbackflg = 0 ! ALF

  istp = 0

  ! ALF
  ! Initialize dtlongn, nndtrat, and nnacoust, and compute the timestep
  ! schedule for all grid operations in all nodes.
  call master_putdtsched(isendflg,isendlite,isendmean,isendboth,1)
  !
  ! Sending DXTMAX to local domains
  call master_putdxt(master_num)


  !         Start the main timestep loop

  write(*,'(/,a,/)') ' === Time integration starts ==='
  do while (time<timmax)

     istp = istp + 1

!!!! Reinitialize subdomains - based on Isendflg computed below
!!!! This part is for the dynamic balancing and sending new varfile/sst info

     if (isendflg==1) then

        if (load_bal==1) then
           call node_decomp(0)
           call dump_Domain_Decomposition()
           call masterput_grid_dimens(master_num)
        endif

     endif

     if (isendbackflg==1) then
        call master_sendinit()
     endif

     !            CPU timing information

     call timing(1,t1)
     wtime1=walltime(wtime_start)

     !            ISENDFLG designates whether nodes should send back
     !               stuff things it normally doesn't have to
     !               at the end of timestep for history/analysis write,
     !               load balancing, etc.

     !               Determines whether nodes send stuff back at the END of the
     !               timestep!!!!!

     call comm_time(isendflg,isendlite,isendmean,isendboth,isendbackflg)

     !            Examine Courant numbers in case model needs to be stopped
     !            or (if ideltat < 0), to update dtlongn, nndtrat,
     !            nnacoust, sspct and isched.

     !call dtset(nndtflg)
     call dtset_new(0, nndtflg, dxtmax_local)

     if (iflag > 0) then   !Estourou Courant - modelo ira´ parar. Antes salva tudo
        isendflg = 1
        isendlite = 1
        isendmean = 1
        isendboth = 1
     endif

     if (nndtflg > 0) then
        call dump_dtset(nndtflg)
        call modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
        call dump_modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
     endif

     ! Send timestep schedule and timesteps to nodes only if they have changed.
     !   Need to do it in first timestep regardless...
     ntsend=0
     if(istp == 1 .or. nndtflg > 0) ntsend=1

     !---------------------------------------------------------------------
     !  Bypass the timestep schedule, then update the main time variable.

     !         do npass=1,nsubs
     !         enddo

     time=time+dtlongn(1)

     ! Wait for cpu time, wallclock time, and Courant numbers
     ! cflxy and cflz from nodes.

     call master_getcflcpu()

     ! Send CFL to Recalculate DeltaT if necessary
     if (ideltat < 0) then
        call master_putcflmax(master_num)
     endif

     if(isendflg.eq.1) then
        !                  Wait for whole subdomains from nodes
        call master_getall
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     if(isendlite.eq.1) then
        call master_getanl('LITE')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     if(isendmean.eq.1) then
        call master_getanl('MEAN')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif

     if(isendboth.eq.1) then
        call master_getanl('BOTH')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     endif


     wtime2=walltime(wtime_start)
     call timing(2,t2)

     pcpu=0
     pwall=0
     do n=1,nmachs
        pcpu=pcpu+ptimes(n,1)
        pwall=pwall+ptimes(n,2)
     enddo

     write(*,'(a,i6,a,f12.1,2(a,f12.2),a)') &
          ' Timestep ',istp,&
          '; Sim Time',time,&
          's; Wall',wtime2-wtime1,&
          's; sum CPU slaves',pcpu,'s'

     call rams_output()

  enddo

  wtime_tot=walltime(wtime_start)
  write(c0,'(f10.1)') wtime_tot
  write(*,'(/,a,/)') ' === Time integration ends; Total elapsed time='//&
                     trim(adjustl(c0))//' ===' 
  return
end subroutine par_MODEL
