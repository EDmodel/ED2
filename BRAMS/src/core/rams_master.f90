!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine rams_master(ipara, nslaves, master_num, name_name)

  use grid_dims, only: &
       maxgrds,        &
       BRAMS_version

  use rpara, only: &
       iparallel
         
  use node_mod, only : nmachs,mynum

  use mem_grid, only: &
       runtype,         &
       expnme,          &
       ngrids,          &
       nzp,             &
       nxp,             &
       nyp,             &
       isched,          &
       maxsched,        &
       maxschent,       &
       nxtnest,         &
       nndtrat,         &
       nsubs,           &
       time,            &
       timmax

  use mem_oda, only:if_oda

  use mem_radiate, only: ISWRTYP, ILWRTYP ! Intent(in)
  use mem_leaf, only : isfcl ! Intent(in)
  use dtset, only: dtset_new ! subroutine

  implicit none

  integer, intent(in) :: ipara              ! 0 iff sequential run; 1 iff parallel run
  integer, intent(in) :: nslaves            ! number of slaves on a parallel run
  integer, intent(in) :: master_num         ! this process rank on a parallel run
  character(len=*), intent(in) :: name_name ! namelist file name

  integer, allocatable :: taskids(:)
  character(len=12) :: c0
  character(len=12) :: c1

  integer :: i,ifm,nndtflg

  real :: w1,w2,w3,wtime_start  ! wall time
  real, external :: walltime    ! wall time
  real :: t1,t2                 ! cpu time

  !ALF
  real :: dxtmax_local(maxgrds)

  ![MLO - For MPI
  integer :: ierr
  include 'interface.h'
  include 'mpif.h'
  !MLO]

  ! wall time 

  wtime_start=walltime(0.)
  w1=walltime(wtime_start)

  ! setup iparallel

  iparallel=ipara

  ! start graphics package (just in case we put in graphical stuff)

  call opngks

  ! read namelist file 

  call read_nl(trim(name_name))

  ! print initial banner

  write(*,"(a)") '+-----------------------------------------------------'
  write(*,"(a)") '!             '//BRAMS_version
  write(*,"(a)") '! input namelist filename is '//trim(name_name)
  write(*,"(a)") '! run name is '//trim(expnme)
  if (ipara == 0) then
     write(*,'(a)') '! single process execution on '//&
           trim(adjustl(runtype))//' run'
  else
     write(c0,'(i12)') nslaves
     write(*,'(a)') '! parallel execution with master and '//&
           trim(adjustl(c0))//' slave processes on '//&
           trim(adjustl(runtype))//' run'
     if (iparallel == 0) then
        write(*,'(a)') '! converted into a sequential run for '//&
           trim(adjustl(runtype))//' run'
     end if
  end if
  write(*,"(a)") '+-----------------------------------------------------'

  ! dump namelist

  call nameout()

  ! Various option settings that should normally not be changed
  ! (previously namelist parameters)

  call eng_params  

  ! Reset parallel flag if necessary

  if ( (runtype(1:4) == 'MAKE') .and. iparallel /= 0) then
     iparallel=0
  endif
  
  ! Initiliazing mynum and nmachs at node_mod. It will be used only to abort runs for
  ! now.
  mynum=0
  if (iparallel == 0) then
     nmachs = 0
  else
     nmachs = nslaves
  end if

  ! First check of options, mainly for numbers of grid points

  call opspec1

  ! Basic grid coordinate setup

  call grid_setup(1)


  ! Additional checks, mainly for nesting

  call opspec2()


  ! Check sfc,sst,ndvi files; remake if needed

  call make_sfcfiles()


  ! Behave accordingly to run type

  if (runtype(1:7) == 'MAKESFC') then

     ! done on a MAKESFC run

     write(*,"(a)") ' MAKESFC run complete'


  else if(runtype(1:9) == 'MAKEVFILE') then

     ! on a "MAKEVFILE" run, call ISAN, then exit.

     call isan_driver(name_name)
     write(*,"(a)") ' MAKEVFILE run complete'


  else

     !-----------------------------------------------------------
     ! If we got here, we are doing an actual
     !    simulation (RAMS or HYPACT)
     !-----------------------------------------------------------

     ! Initialize micro arrays. May need to change some settings which affect memory.

     call jnmbinit()


     ! Allocate main memory

     if (iparallel == 0) then
        call rams_mem_alloc(0)
     else
        call rams_mem_alloc(1)
     endif


     ! Call main initialization driver

     call initlz(name_name)


     ! Compute Courant numbers cflxy and cflz.

     do ifm = 1,ngrids
        call newgrid(ifm)
        call cfl(nzp,nxp,nyp,0,0,1)
     enddo


     ! Initialize dtlongn, nndtrat, and nnacoust, and compute the timestep
     ! schedule for all grid operations.

     !call dtset(nndtflg)
     call dtset_new(0, nndtflg, dxtmax_local)

     call dump_dtset(nndtflg)

     call modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
     call dump_modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)


     !  Initialize HYPACT configuration and sources

     if (index(runtype,'HYP') /= 0)  then
        write (*,'(a)') 'ERROR** HYPACT not implemented'
        stop
     endif

     if(iparallel == 1) then

        !   ---------------------------------------------------------------
        !     Initialize parallel processes with all relevant information
        !   --------------------------------------------------------------

        allocate(taskids(nslaves))
        do i=1,nslaves
           taskids(i)=i
        end do
        write (unit=*,fmt='(a)') ' Waiting for all nodes to reach this barrier'
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        write (unit=*,fmt='(a)') ' - Masterput_processid'
        call masterput_processid(nslaves,taskids,master_num)
        write (unit=*,fmt='(a)') ' - Masterput_nl'
        call masterput_nl(master_num)
        if (isfcl == 5) then
           write (unit=*,fmt='(a)') ' - Masterput_ednl'
           call masterput_ednl(master_num)
        end if
        write (unit=*,fmt='(a)') ' - Masterput_gridinit'
        call masterput_gridinit(master_num)
        write (unit=*,fmt='(a)') ' - Node_decomp'
        call node_decomp(.true.)
        write (unit=*,fmt='(a)') ' - dump_Domain_Decomposition'
        call dump_Domain_Decomposition()
        write (unit=*,fmt='(a)') ' - Masterput_grid_dimens'
        call masterput_grid_dimens(master_num)
        write (unit=*,fmt='(a)') ' - Masterput_grid_gridset'
        call masterput_gridset(master_num)
        write (unit=*,fmt='(a)') ' - Masterput_grid_cofnest'
        call masterput_cofnest(master_num)
        write (unit=*,fmt='(a)') ' - Masterput_grid_micphys'
        call masterput_micphys(master_num)
        if (if_oda == 1) then
           write (unit=*,fmt='(a)') ' - Masterput_oda'
           call masterput_oda(master_num)
        end if
        write (unit=*,fmt='(a)') ' - Masterput_misc'
        call masterput_misc(master_num)
        write (unit=*,fmt='(a)') ' - Master_sendinit'
        call master_sendinit()
        write (unit=*,fmt='(a)') ' - Master_ed_init'
        call master_ed_init(iparallel)
        write (unit=*,fmt='(a)') ' - Second Barrier'
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
     else
        call master_ed_init(iparallel)
     end if

     call timing(1,t1)
     w2=walltime(wtime_start)
     write(c0,'(f12.2)') t1
     write(c1,'(f12.2)') w2-w1
     write(*,'(/,a,/)') ' === Finish initialization; CPU(sec)='//&
                        trim(adjustl(c0))//'; Wall(sec)='//trim(adjustl(c1))//&
                        '; Time integration starts (rams_master) ===' 

     ! Exit if doing a zero time run
     if (time < timmax) then
        !  Call the model time integration driver
        if(iparallel == 1) then
           call par_model(master_num)
        else
           call model()
        end if
     end if
  end if

  !  RAMS finished, clean up some last things...

  if (iparallel == 1) then
     deallocate(taskids)
  end if

  call clsgks     ! call this just in case we put in graphical stuff

  return
end subroutine rams_master

!-------------------------------------------------------------------------

subroutine comm_time(isendflg,isendlite,isendmean,isendboth,isendbackflg)

  use mem_varinit
  use mem_cuparm
  use io_params
  use mem_grid

  ! CATT
  use catt_start, only: CATT           ! intent(in)

  implicit none

  integer, intent(out) :: isendflg,isendlite,isendmean,isendboth,isendbackflg
  real(kind=8) :: timemf
  integer :: ifm
  real :: never
  real, parameter  :: frqqueim=86400.  ! CATT
  ! These variables are lowest common denomenators with model time
  ! This simplifies the application of the double into the mod functions
  ! RGK 5-29-07
  real :: time_fein,time_frql,time_frqm,time_frqb
  real :: time_frqh,time_frqa,time_frqp

  !         ISENDFLG designates whether nodes should send back
  !            stuff things it normally doesn't have to
  !            at the end of timestep for history/analysis write,
  !            load balancing, etc.

  !         isendflg  = the usual RAMS stuff
  !         isendlite = the "lite" variables
  !         isendmean = the "mean" variasbles
  !         isendboth = Both the "mean" and "lite" variables

  !            Determines whether nodes send stuff back at the END of the
  !            timestep!!!!!

  timemf = time + dble(dtlongn(1))
  never  = 1.25 * sngl(timmax)
  if (frqqueim > 0) then
    time_fein = real(dmod( timemf, dble(frqqueim)))
  else
    time_fein = never
  end if
  if (frqlite > 0) then
    time_frql = real(dmod( timemf, dble(frqlite )))
  else
    time_frql = never
  end if
  if (frqmean > 0) then
    time_frqm = real(dmod( timemf, dble(frqmean )))
  else
    time_frqm = never
  end if
  if (frqboth > 0) then
    time_frqb = real(dmod( timemf, dble(frqboth )))
  else
    time_frqb = never
  end if
  time_frqh = real(dmod( timemf, dble(frqhis  )))
  time_frqa = real(dmod( timemf, dble(frqanl  )))
  if (frqprt > 0) then
    time_frqp = real(dmod( timemf, dble(frqprt )))
  else
    time_frqp = never
  end if


  isendflg = 0
  isendlite = 0
  isendmean = 0
  isendboth = 0
  isendbackflg = 0 ! ALF

  if (CATT == 1) then
     !----srf-

     !inicializando as fontes as 00UTC
     if ( mod(time_fein + frqqueim - 0.*3600., frqqueim ) .lt. dtlongn(1) .and. &
          timemf.ge.0.*3600.) then
        isendflg = 1
        isendbackflg = 1 ! ALF
        !        print*,'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
        !        print*,'BURN MAP READING: ',isendflg,timemf/3600.
        !        print*,'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
     endif
     !----srf-
  endif

  if ( time_frql  < dtlongn(1)) isendlite=1
  
  if(avgtim > 0.)then
     if(mod(time_frqm + frqmean - avgtim/2.,frqmean) < dtlongn(1) &
     .and. timemf >= avgtim) isendmean=1    !RGK
  elseif(avgtim < 0.)then
     if( time_frqm < dtlongn(1)) isendmean=1
  endif

  if(runtype(1:7) == 'INITIAL'.and.timemf < dtlongn(1)) isendmean=0
  if(runtype(1:7) == 'HISTORY'.and.timemf <= timstr) isendmean=0

  if(avgtim > 0.)then
      if(mod(time_frqb + frqboth - avgtim/2.,frqboth) < dtlongn(1) &
           .and. timemf >= avgtim) isendboth=1
   elseif(avgtim < 0.)then
      if( time_frqb < dtlongn(1)) isendboth=1
  endif
  
  if(runtype(1:7) == 'INITIAL'.and.timemf < dtlongn(1))isendboth=0
  if(runtype(1:7) == 'HISTORY'.and.timemf <= timstr)isendboth=0

  if (ioutput  /=  0) then
     if ( time_frqa  <  dtlongn(1) .or. time_frqh  <  dtlongn(1)) then
        isendflg = 1
        !return
     endif
  endif

  if( timemf  >=  timmax - .01*dtlongn(1) ) then
     isendflg = 1
     !return
  endif

  if  ( nud_type == 2 .and. timemf  >=  vtime2 .and. timemf  <  timmax) then
     isendflg = 1
     isendbackflg = 1 ! ALF
     return
  endif

  if  ( nud_type == 1 .and. timemf  >=  htime2 .and. timemf  <  timmax) then
     isendflg = 1
     isendbackflg = 1 ! ALF
     !return
  endif

  if  ( nud_cond == 1 .and. timemf  >=  condtime2 .and. timemf  <  timmax) then
     isendflg = 1
     isendbackflg = 1 ! ALF
     !return
  endif

  if ( time_frqp  <  dtlongn(1) ) then
     isendflg = 1
     !return
  endif

  if (iupdsst  ==  1 ) then
     do ifm = 1,ngrids
        if (isstflg(ifm)  ==  1) then
           if (timemf  >=  ssttime2(ifm) .and.  &
                timemf  <  timmax) then
              isendflg = 1
              isendbackflg = 1 ! ALF
              return

           endif
        endif
     enddo
  endif

  if (iupdndvi  ==  1 ) then
     do ifm = 1,ngrids
        if (ndviflg(ifm)==1) then
           if (timemf>=ndvitime2(ifm) .and. timemf<timmax) then
              isendflg = 1
              isendbackflg = 1 ! ALF

              !return
           endif
        endif
     enddo
  endif

  if (if_cuinv  ==  1 ) then
     do ifm = 1,ngrids
        if (timemf  >=  cu_times(ncufl+1) .and.  &
             timemf  <  timmax) then
           isendflg = 1
           isendbackflg = 1 ! ALF
           !return
        endif
     enddo
  endif

  return
end subroutine comm_time

!-------------------------------------------------------------------------

subroutine comm_time_new(isendflg, isendlite, isendmean, isendboth, &
     isendbackflg, isendiv, isendsst, isendndvi)

  use mem_varinit
  use mem_cuparm
  use io_params
  use mem_grid
  ! CATT
  use catt_start, only: CATT           ! intent(in)
  ! ALF
  use node_mod, only : mynum           ! intent(in)
  implicit none

  integer, intent(out) :: isendflg,isendlite,isendmean,isendboth,isendbackflg
  integer, intent(out) :: isendiv, isendsst, isendndvi

  real :: time_fein,time_frql,time_frqm,time_frqb
  real :: time_frqh,time_frqa,time_frqp

  real(kind=8) :: timemf
  real, parameter :: frqqueim=86400.  ! CATT
  real            :: never
  integer :: ifm

  ! ALF
  integer :: iyears,imonths,idates,ihours
  character(len=14)  :: itotdate_current
  real(kind=8) :: lastdate_sec, initialdate_sec

  !         ISENDFLG designates whether nodes should send back
  !            stuff things it normally doesn't have to
  !            at the end of timestep for history/analysis write,
  !            load balancing, etc.

  !         isendflg  = the usual RAMS stuff
  !         isendlite = the "lite" variables
  !         isendmean = the "mean" variasbles
  !         isendboth = Both the "mean" and "lite" variables

  !            Determines whether nodes send stuff back at the END of the
  !            timestep!!!!!

  never=1.25*sngl(timmax)
  timemf = time + dble(dtlongn(1))

  if (frqqueim > 0) then
    time_fein = real(dmod( timemf, dble(frqqueim)))
  else
    time_fein = never
  end if
  if (frqlite > 0) then
    time_frql = real(dmod( timemf, dble(frqlite )))
  else
    time_frql = never
  end if
  if (frqmean > 0) then
    time_frqm = real(dmod( timemf, dble(frqmean )))
  else
    time_frqm = never
  end if
  if (frqboth > 0) then
    time_frqb = real(dmod( timemf, dble(frqboth )))
  else
    time_frqb = never
  end if
  time_frqh = real(dmod( timemf, dble(frqhis  )))
  time_frqa = real(dmod( timemf, dble(frqanl  )))
  if (frqprt > 0) then
    time_frqp = real(dmod( timemf, dble(frqprt )))
  else
    time_frqp = never
  end if



  isendflg = 0
  isendlite = 0
  isendmean = 0
  isendboth = 0
  isendbackflg = 0 ! ALF
  isendiv = 0 ! ALF
  isendsst = 0 ! ALF
  isendndvi = 0 ! ALF


  if (CATT == 1) then
     !----srf-
     !time para leitura dos mapas de queimadas
     !inicializando as fontes as 00UTC
     if ( mod(time_fein + frqqueim -0.*3600.,frqqueim) .lt. dtlongn(1) .and. &
          timemf.ge.0.*3600.) then
        isendflg = 1
        isendbackflg = 1 ! ALF
        !        print*,'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
        !        print*,'BURN MAP READING: ',isendflg,timemf/3600.
        !        print*,'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
     endif
     !----srf-
  endif

  if (time_frql < dtlongn(1)) isendlite=1

  if(avgtim > 0.)then
     if(mod( time_frqm + frqmean - avgtim/2.,frqmean) < dtlongn(1) .and.  &
          timemf >= avgtim) isendmean=1
  elseif(avgtim < 0.)then
     if( time_frqm < dtlongn(1)) isendmean=1
  endif
  if(runtype(1:7) == 'INITIAL'.and.timemf < dtlongn(1)) isendmean=0
  if(runtype(1:7) == 'HISTORY'.and.timemf <= timstr) isendmean=0

  if(avgtim > 0.)then
     if(mod(time_frqb + frqboth - avgtim/2.,frqboth) < dtlongn(1) .and.  &
          timemf >= avgtim) isendboth=1
  elseif(avgtim < 0.)then
     if( time_frqb < dtlongn(1)) isendboth=1
  endif
  if(runtype(1:7) == 'INITIAL'.and.timemf < dtlongn(1))isendboth=0
  if(runtype(1:7) == 'HISTORY'.and.timemf <= timstr)isendboth=0

  if (ioutput  /=  0) then
     if ( time_frqa <  dtlongn(1) .or.  &
          time_frqh <  dtlongn(1)) then
        isendflg = 1
        !return
     endif
  endif

  if( timemf  >=  timmax - .01*dtlongn(1) ) then
     isendflg = 1
     !return
  endif

  if (nud_type==2) then

     if (mynum/=0) then
        if  (timemf>=vtime2 .and. timemf<timmax) then
           isendflg = 1
           isendbackflg = 1 ! ALF
           isendiv = 1 ! ALF
        endif
     else ! Master node
        call date_make_big(iyeara,imontha,idatea,itimea,itotdate_current)
        call date_abs_secs(itotdate_current, initialdate_sec)
        call date_abs_secs(lastdate_iv, lastdate_sec)
        if  (timemf>=vtime1 .and. timemf<timmax .and.         &
             ((lastdate_sec-initialdate_sec)+.0001)<timmax) then
           call date_add_to(iyeara,imontha,idatea,itimea*100  &
                ,timmax,'s',iyears,imonths,idates,ihours)
           call date_make_big(iyears,imonths,idates,ihours,itotdate_current)
           if (lastdate_iv<itotdate_current) then
              isendflg = 1
              isendbackflg = 1 ! ALF
              isendiv = 1 ! ALF
           endif
        endif
        !return
     endif

  endif

  if  ( nud_type == 1 .and. timemf  >=  htime2  &
       .and. timemf  <  timmax) then
     isendflg = 1
     isendbackflg = 1 ! ALF
     !return
  endif

  if  ( nud_cond == 1 .and. timemf  >=  condtime2  &
       .and. timemf  <  timmax) then
     isendflg = 1
     isendbackflg = 1 ! ALF
     !return
  endif

  if ( time_frqp  <  dtlongn(1) ) then
     isendflg = 1
     !return
  endif

  if (iupdsst  ==  1 ) then

     if (mynum/=0) then ! Slave Process
        do ifm = 1,ngrids
           if (isstflg(ifm)==1) then
              if (timemf>=ssttime2(ifm) .and. timemf<timmax) then
                 isendflg = 1
                 isendbackflg = 1 ! ALF
                 isendsst = 1 ! ALF
                 !return
              endif
           endif
        enddo
     else
        ! Master Process
        ! First, checking if all values have been read
        call date_make_big(iyeara,imontha,idatea,itimea*100,itotdate_current)
        call date_abs_secs(itotdate_current, initialdate_sec)
        do ifm = 1,ngrids
           call date_abs_secs(lastdate_sst(ifm), lastdate_sec)
           if (isstflg(ifm)==1 .and. timemf>=ssttime1(ifm) .and. &
                timemf<timmax .and. &
                ((lastdate_sec-initialdate_sec)+.0001)<timmax) then
              call date_add_to(iyeara,imontha,idatea,itimea*100  &
                    ,timmax,'s',iyears,imonths,idates,ihours)
              call date_make_big(iyears,imonths,idates,ihours,itotdate_current)
              if (lastdate_sst(ifm)<itotdate_current) then
                 isendflg = 1
                 isendbackflg = 1 ! ALF
                 isendsst = 1 ! ALF
                 !return
              endif
           endif
        enddo
     endif

  endif

  if (iupdndvi  ==  1 ) then
     do ifm = 1,ngrids
        if (ndviflg(ifm)==1) then
           if (timemf>=ndvitime2(ifm) .and. timemf<timmax) then
              isendflg = 1
              isendbackflg = 1 ! ALF
              isendndvi = 1
              !return
           endif
        endif
     enddo
  endif

  if (if_cuinv  ==  1 ) then
     do ifm = 1,ngrids
        if (timemf  >=  cu_times(ncufl+1) .and.  &
             timemf  <  timmax) then
           isendflg = 1
           isendbackflg = 1 ! ALF
           !return
        endif
     enddo
  endif

  return
end subroutine comm_time_new

!-------------------------------------------------------------------------

subroutine rams_output()

  use mem_varinit, only: &
       nud_type,  & ! INTENT(IN)
       htime2,    & ! INTENT(IN) 
       nud_cond,  & ! INTENT(IN)
       condtime2, & ! INTENT(IN)
       vtime2       ! INTENT(IN)

  use mem_cuparm, only: &
       if_cuinv,  & ! INTENT(IN)
       cutime2      ! INTENT(IN)

  use io_params, only: &
       frqprt,  & ! INTENT(IN)
       INITFLD, & ! INTENT(IN)
       frqhis,  & ! INTENT(IN)
       frqanl,  & ! INTENT(IN)
       frqlite, & ! INTENT(IN)
       frqmean, & ! INTENT(IN)
       frqboth, & ! INTENT(IN)
       avgtim,  & ! INTENT(IN)
       iupdsst, & ! INTENT(IN)
       iupdndvi,& ! INTENT(IN)
       ioutput    ! ADDED [RGK]

  use mem_grid, only: &
       ngrids,        & ! INTENT(IN)
       nzp, nxp, nyp, & ! INTENT(IN)
       npatch,        & ! INTENT(IN)
       time,          & ! INTENT(IN)
       dtlongn,       & ! INTENT(IN)
       iflag,         & ! INTENT(IN)
       timmax           ! INTENT(IN)

  use node_mod, only: &
       ia, iz, ja, jz                  ! INTENT(IN)

  ! CATT
  use emission_source_map, only: read_emission_sources_map ! Subroutine

  ! CATT
  use catt_start, only: CATT           ! intent(in)

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM      ! INTENT(IN)
  use mem_emiss, only : isource         ! INTENT(IN)
  
  ![MLO-ED2
  use mem_mass, only:     &
        mass_g,           & ! intent(inout)
        imassflx,         & ! intent(in)
        frqmassave        ! ! intent(in)
  !MLO]

  implicit none

  ! Local Variables
  integer        :: ierr, ifm      !,ifileok
  ! This simplifies the application of the double into the mod functions
  ! RGK 5-29-07
  real :: time_fein,time_frql,time_frqm,time_frqb
  real :: time_frqh,time_frqa,time_frqp,time_frqz
  real,parameter :: frqqueim=86400.  ! CATT
  real           :: never

  never=1.25*sngl(timmax)

  if (frqqueim > 0) then
    time_fein = real(dmod( time, dble(frqqueim)))
  else
    time_fein = never
  end if
  if (frqlite > 0) then
    time_frql = real(dmod( time, dble(frqlite )))
  else
    time_frql = never
  end if
  if (frqmean > 0) then
    time_frqm = real(dmod( time, dble(frqmean )))
  else
    time_frqm = never
  end if
  if (frqboth > 0) then
    time_frqb = real(dmod( time, dble(frqboth )))
  else
    time_frqb = never
  end if
  time_frqh = real(dmod( time, dble(frqhis  )))
  time_frqa = real(dmod( time, dble(frqanl  )))
  if (frqprt > 0) then
    time_frqp = real(dmod( time, dble(frqprt )))
  else
    time_frqp = never
  end if
  if (imassflx == 1) then ! Time to make averages zero
     time_frqz = real (dmod( time, dble(frqmassave)))
  else
     time_frqz = never
  end if
!MLO]

  if (TEB_SPM==1) then
     !EDF  emission module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !----------------------------------------
     if(isource==1)then
        ifm=ngrids
        call newgrid(ifm)
        CALL le_fontes(ifm, nzp, nxp, nyp, npatch, ia, iz, ja, jz) 
     endif
     !EDF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  endif

  do ifm = 1,ngrids
     call newgrid(ifm)

     if(( time_frqp  <dtlongn(1).and.INITFLD==1) .or.iflag == 1)then
        call prtout()
        !c            call uwcomp(a)
     endif

  enddo


  if( time_frqa < dtlongn(1) .or. time  >=  timmax - .01*dtlongn(1) ) then
     !  ADDED HDF - [RGK] 12-20-07
     if (ioutput == 3) then
        call anlhdf('INST')
     elseif (ioutput .gt. 0 ) then
        call anlwrt('no','INST')
     end if
  end if

  !     Call the analysis writing routine again for the other var types
  if( time_frql < dtlongn(1) ) then
     !  ADDED HDF - [RGK] 12-20-07

     if (ioutput == 3) then
        call anlhdf('LITE')
     elseif (ioutput .gt. 0 ) then
        call anlwrt('no','LITE')
     end if
  end if

  if (frqmean > 0.) then 
    if(avgtim > 0.0.and.mod(time_frqm + frqmean - avgtim/2.,frqmean) < dtlongn(1)  &
         .and.time >= avgtim) call anlwrt('no','MEAN')
    if(avgtim < 0.0.and. time_frqm  < dtlongn(1))  &
         call anlwrt('no','MEAN')
  end if

  if (frqboth > 0.) then
    if(avgtim > 0.0.and.mod(time_frqb + frqboth - avgtim/2.,frqboth) < dtlongn(1)  &
         .and.time >= avgtim)call anlwrt('no','BOTH')
    if(avgtim < 0.0.and. time_frqb < dtlongn(1))  &
         call anlwrt('no','BOTH')
  end if
  
  !Moved history to here, so I can use some ED structures that were
  !already updated
  if( time_frqh < dtlongn(1).or. time  >=  timmax - .01*dtlongn(1) .or.  &
       iflag == 1)then
     call hiswrt('no')
  end if

  if (iupdsst  ==  1 .and. time+0.00001 < timmax) then
     do ifm = 1,ngrids
        call sst_read(3,ifm,ierr)
     enddo
  endif

  if (iupdndvi  ==  1 .and. time+0.00001 < timmax) then
     do ifm = 1,ngrids
        call ndvi_read(3,ifm,ierr)
     enddo
  endif

  if( nud_type == 1 .and. time >= htime2 .and. time+.00001 < timmax) then
     call nud_read(2)
  endif

  if( nud_cond == 1 .and. time >= condtime2 .and.time+.00001 < timmax) then
     call cond_read(2)
  endif

  if(if_cuinv == 1 .and. time >= cutime2 .and. time+.00001 < timmax) then
     call cu_read(2)
  endif

  if(nud_type==2) then
     if (time>=vtime2 .and. time+.0001<timmax) then
        call varf_read(2)
     endif
  endif

  ! CATT
  if (CATT == 1) then
     ! Initiating sources at 00UTC
     if ( mod(time_fein + frqqueim - 0.*3600.,frqqueim) .lt. dtlongn(1) .and. &
          time.ge.0.*3600.) then
        ! Reading emission map for all grids
        !        print*,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
        !        print*,'BURN MAP READING FOR ALL GRIDS AT TIME FREQUENCY: ', time/3600.
        call read_emission_sources_map()
        !        print*,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'

     endif
  endif

  if (iflag==1) stop 'IFLAG'

  return
end subroutine rams_output
