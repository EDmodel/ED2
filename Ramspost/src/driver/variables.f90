!==========================================================================================!
!==========================================================================================!
subroutine RAMS_anal_init(nfile,fnames,file_prefix,dep_zlev,iep_nx,iep_ny,iep_nz,iep_ng    &
                         ,iep_np,iep_nc,iep_ngrids)

   use an_header
   use rpost_dims
   use rpost_coms
   use brams_data
   use misc_coms , only : memsize4  ! ! intent(inout)
   use rconstants, only : day_sec   & ! intent(in)
                        , hr_sec    & ! intent(in)
                        , min_sec   ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                          , intent(out)   :: nfile
   character(len=fnm_len), dimension(maxfiles)      , intent(in)    :: fnames
   real                  , dimension(nzpmax,maxgrds), intent(inout) :: dep_zlev
   integer                                          , intent(inout) :: iep_np
   integer                                          , intent(inout) :: iep_nc
   integer                                          , intent(inout) :: iep_ng
   integer                                          , intent(inout) :: iep_ngrids
   integer               , dimension(       maxgrds), intent(inout) :: iep_nx
   integer               , dimension(       maxgrds), intent(inout) :: iep_ny
   integer               , dimension(       maxgrds), intent(inout) :: iep_nz
   !----- Internal variables. -------------------------------------------------------------!
   character(len=fnm_len)                                           :: file_prefix
   character(len=fnm_len)                                           :: fpref
   integer               , dimension(13)                            :: mondays
   integer                                                          :: maxmem
   integer                                                          :: nc
   integer                                                          :: nfn
   integer                                                          :: nv
   integer                                                          :: istrhrs
   integer                                                          :: n
   integer                                                          :: nn
   integer                                                          :: ng
   integer                                                          :: k
   integer                                                          :: itme
   integer                                                          :: iadddays
   integer                                                          :: izhours
   integer                                                          :: izmin
   integer                                                          :: izsec
   integer                                                          :: iiyear
   integer                                                          :: iidate
   integer                                                          :: iimon
   real                                                             :: strtim
   !----- External functions. -------------------------------------------------------------!
   logical                                          , external      :: isleap
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Initialise the number of days in each month.  February will be updated soon if    !
   ! this is a leap year.                                                                  !
   !---------------------------------------------------------------------------------------!
   mondays(:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31 /)

   !----- Find file name size and bind the header suffix. ---------------------------------!
   maxmem     = 0
   fpref      = file_prefix
   nc         = len_trim(fpref)+1
   nfile      = -1
   fpref      = file_prefix
   fpref(nc:) = '*-head.txt'
   write (unit=*,fmt='(a,1x,a)') ' => RAMS_filelist searching for: ',fpref(1:nc+10)

   call RAMS_filelist(fnames,fpref,nfile)

   !---------------------------------------------------------------------------------------!
   !    Read the analysis header and extract information about every variable, such as the !
   ! type of variable, the array size and location in the binary file.                     !
   !---------------------------------------------------------------------------------------!
   do nfn=1,nfile
      open (unit=10,file=fnames(nfn),status='old',action='read',form='formatted')

      read (unit=10,fmt=*) nvbtab

      allocate (anal_table(nvbtab))
      do nv=1,nvbtab
         read(unit=10,fmt=*)  anal_table(nv)%string   , anal_table(nv)%npointer            &
                            , anal_table(nv)%idim_type, anal_table(nv)%ngrid               &
                            , anal_table(nv)%nvalues
      end do
      !----- Read the BRAMS settings and save it. -----------------------------------------!
      call commio('ANAL','READ',10)
      close(unit=10,status='keep')

      ftimes(nfn) = time
      istrhrs     = nint(float(itime1)/100.+0.0001)
      strtim      = istrhrs+float(itime1-istrhrs*100)/60.
      startutc    = strtim
      deallocate (anal_table)
   end do

   call RAMS_fltsort(nfile,ftimes,fnames)

   !----- Loop over all files. ------------------------------------------------------------!
   do nfn=1,nfile
      open(unit=10,file=fnames(nfn),status='old',action='read',form='formatted')
      read(unit=10,fmt=*) nvbtab
      allocate (anal_table(nvbtab))
      do nv=1,nvbtab
         read(unit=10,fmt=*)  anal_table(nv)%string   , anal_table(nv)%npointer            &
                            , anal_table(nv)%idim_type, anal_table(nv)%ngrid               &
                            , anal_table(nv)%nvalues
      end do
      call commio('ANAL','READ',10)
      close(unit=10,status='keep')

      ftimes(nfn)=time

      if (nfn == 1) then
         iep_ngrids=ngrids
         do n=1,ngrids
            maxmem    = max(maxmem,nnxp(n)*nnyp(n)*max(nnzp(n),npatch*nzg,nnzp(n)*nclouds))
            iep_nx(n) = nnxp(n)
            iep_ny(n) = nnyp(n)
            iep_nz(n) = nnzp(n)
            iep_ng    = nzg
            iep_np    = npatch
            iep_nc    = nclouds
            do nn=1,nnzp(n)
               dep_zlev(nn,n)=ztn(nn,n)
            end do
         end do
      end if

      itme     = nint(strtim*hr_sec + time)
      iadddays = itme / day_sec
      izhours  = (itme - iadddays * day_sec) / hr_sec
      izmin    = (itme - iadddays * day_sec - izhours * hr_sec) / min_sec
      izsec    =  itme - iadddays * day_sec - izhours * hr_sec - izmin * min_sec

      iftimes(nfn) = izhours
      iftimes(nfn) = iftimes(nfn) *100 + izmin
      iftimes(nfn) = iftimes(nfn) *100 + izsec
      select case (iyear1)
      case(:50)
         iyear1 = iyear1 + 2000

      case(51:99)
         iyear1 = iyear1 + 1900

      end select

      if (isleap(iyear1)) mondays(2) = 29

      iiyear = iyear1
      iidate = idate1+iadddays
      iimon  = imonth1
      do while (iidate > mondays(iimon))
         iidate = iidate - mondays(iimon)
         iimon  = iimon + 1
      end do
      do while (iimon > 12)
         iiyear = iiyear + 1
         iimon  = iimon  - 12
      end do
   
      ifdates(nfn) = iiyear       * 100  + iimon
      ifdates(nfn) = ifdates(nfn) * 100  + iidate

      nfgrids(nfn) = ngrids

      do ng=1,ngrids
         nfgpnts(1,ng,nfn) = nnxp(ng)
         nfgpnts(2,ng,nfn) = nnyp(ng)
         nfgpnts(3,ng,nfn) = nnzp(ng)
         nfgpnts(4,ng,nfn) = nzg
         fdelx(ng,nfn)     = DELTAXN(NG)
         fdely(ng,nfn)     = DELTAYN(NG)

         do k=1,nnzp(ng)
            flevels(k,ng,nfn) = ztn(k,ng)
         end do
      end do

      httop = zmn(nnzp(1)-1,1)

      close(unit=10,status='keep')

      if(nfn < nfile) deallocate (anal_table)
   end do

   memsize4=maxmem

   return
end subroutine RAMS_anal_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_get_time_init(nfl,iyear,imonth,idate,ihour,imin)
   use brams_data, only : iftimes ! ! intent(in)
   use rpost_coms, only : iyear1  & ! intent(in)
                        , imonth1 & ! intent(in)
                        , idate1  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)  :: nfl
   integer, intent(out) :: iyear
   integer, intent(out) :: imonth
   integer, intent(out) :: idate
   integer, intent(out) :: ihour
   integer, intent(out) :: imin
   !---------------------------------------------------------------------------------------!
   
   iyear =iyear1
   imonth=imonth1
   idate =idate1
   ihour =int(float(iftimes(nfl))/10000.)
   imin  =int(float(iftimes(nfl)-10000*ihour)/100.)
   return
end subroutine RAMS_get_time_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_get_time_step(iistep,hunit,nfiles)
   use brams_data, only : iftimes
   implicit none 
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)  :: nfiles
   integer, intent(out) :: iistep
   integer, intent(out) :: hunit
   !---------------------------------------------------------------------------------------!
   
   if (nfiles == 1) then
      !----- If only one file exists, assume 1-hr interval, it really doesn't matter... ---!
      iistep = 1
      hunit  = 3
   else
      !----- Compute the difference in time between the two variables. --------------------!
      iistep = iftimes(2) - iftimes(1)

      if (iistep >= 6000) then 
         !----- Time difference is in hours. ----------------------------------------------!
       iistep = int(float(iistep) / 10000.)
       hunit  = 3

      elseif (iistep >= 60) then
         !----- Time difference is in minutes. --------------------------------------------!
         iistep=int(float(iistep)/100.)
         hunit=2

      elseif (iistep < 60) then
         !----- Time difference is in seconds. --------------------------------------------!
         hunit=1

      end if

   end if
   return
end subroutine RAMS_get_time_step
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function will read a variable from the binary file and save it into a temporary  !
! array.  The result of this function will an error message.                               !
!------------------------------------------------------------------------------------------!
integer function RAMS_getvar(stringg,itype,ngrd,a,b,flnm)

   use an_header
   use rpost_dims, only : fnm_len     ! ! intent(inout)
   use misc_coms , only : ifound      ! ! intent(inout)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)              , intent(in)    :: flnm
   character(len=*)              , intent(in)    :: stringg
   real            , dimension(*), intent(inout) :: a
   real            , dimension(*), intent(inout) :: b
   integer                       , intent(out)   :: itype
   integer                       , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   character(len=fnm_len)                        :: flng
   character(len=120)                            :: errmsg
   character(len=20)                             :: string
   integer                                       :: lenf
   integer                                       :: ni
   integer                                       :: npts
   integer                                       :: iword
   logical                                       :: success
   !---------------------------------------------------------------------------------------!

   !----- Define the length of the file name and variable. --------------------------------!
   lenf  = len_trim(flnm)

   !---------------------------------------------------------------------------------------!
   !     Here we check whether we are reading averaged analysis.  If so, append an M to    !
   ! the variable name.                                                                    !
   !---------------------------------------------------------------------------------------!
   select case (flnm(lenf-18:lenf-18))
   case ('M','B')
      string = trim(stringg)//'M'
   case default
      string = trim(stringg)
   end select

   !---------------------------------------------------------------------------------------!
   !      Initialise some variable with the default, error values.  They will be replaced  !
   ! if the variable is found.                                                             !
   !---------------------------------------------------------------------------------------!
   success     = .false.
   RAMS_getvar =  1
   itype       = -1
   errmsg      = '       # Variable not available in this run - '//trim(string)

   !---------------------------------------------------------------------------------------!
   !      Now loop over the variable table, looking for the variable we want.              !
   !---------------------------------------------------------------------------------------!
   searchloop: do ni=1,nvbtab
      success = trim(string) == trim(anal_table(ni)%string) .and.                          &
                ngrd == anal_table(ni)%ngrid
      if (success) then
         !---------------------------------------------------------------------------------!
         !      Found variable, now we build the file name and check whether it exists or  !
         ! not.  If it doesn't, return an error message.                                   !
         !---------------------------------------------------------------------------------!
         write (flng,fmt='(2a,i1.1,2a)') trim(flnm),'-g',ngrd,'.vfm',char(0)
         inquire(file=flng,exist=success)
         if (.not. success) then
            errmsg='       # File not found - '//flng
            call error_mess(errmsg)
            exit searchloop
         end if
         !---------------------------------------------------------------------------------!

         npts  = anal_table(ni)%nvalues
         itype = anal_table(ni)%idim_type
         iword = anal_table(ni)%npointer
         
         write (unit=*,fmt='(a,1x,a,3x,a,1x,i3,3x,2(a,1x,i12,3x))')                        &
                                '       # Reading variable:',trim(anal_table(ni)%string)   &
                              , 'itype =',itype,'Dimension = ',npts,'Pointer =',iword 
         call RAMS_c_open(flng,'r'//char(0))
         call vfirecr(10,a,npts,'LIN',b,iword)
         call RAMS_c_close()


         RAMS_getvar=0
         ifound=ifound+1
         exit searchloop
      end if
   end do searchloop

   !----- Give the bad news to the user if it doesn't work... -----------------------------!
   if (.not. success) call error_mess(errmsg)

   return
end function RAMS_getvar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_varlib(cvar,nx,ny,nz,nsl,npat,ncld,ngrd,flnm,cdname,cdunits,ivar_type &
                      ,a,b,a2,a6)
   use rconstants
   use rpost_coms
   use rpost_dims, only : nwave           ! ! intent(in)
   use leaf_coms , only : ustmin          & ! intent(in)
                        , ubmin           ! ! intent(in)
   use misc_coms , only : memsize4        & ! intent(inout)
                        , ierr_getvar     & ! intent(inout)
                        , ifound          ! ! intent(inout)
   use scratch_coms, only : scr           & ! intent(inout)
                          , alloc_scratch ! ! subroutine
   use micro_coms
   use somevars    , only : co2_on
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: nx
   integer                       , intent(in)    :: ny
   integer                       , intent(in)    :: nz
   integer                       , intent(in)    :: ngrd
   integer                       , intent(in)    :: nsl
   integer                       , intent(in)    :: npat
   integer                       , intent(in)    :: ncld
   integer                       , intent(out)   :: ivar_type
   character(len=*)              , intent(in)    :: cvar
   character(len=*)              , intent(in)    :: flnm
   character(len=*)              , intent(out)   :: cdname
   character(len=*)              , intent(out)   :: cdunits
   real            , dimension(*), intent(inout) :: a
   real            , dimension(*), intent(inout) :: b
   real            , dimension(*), intent(inout) :: a2
   real            , dimension(*), intent(inout) :: a6
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: idim_type
   integer                                       :: irecind
   integer                                       :: irecsize
   integer                                       :: ind
   integer                                       :: ispec
   integer                                       :: ierr
   integer                       , save          :: memsave4 = 0
   !----- External functions. -------------------------------------------------------------!
   integer                       , external      :: RAMS_getvar
   !---------------------------------------------------------------------------------------!

   if (memsave4 == 0) then
      write (unit=*,fmt='(a,1x,i12)') '   - Allocating scratch variables; Size=',memsize4
      call alloc_scratch(scr,memsize4)
      memsave4 = memsize4
   elseif (memsize4 > memsave4) then
      write (unit=*,fmt='(a,1x,i12)') '   - Re-allocating scratch variables; New size='    &
                                     ,memsize4
      call alloc_scratch(scr,memsize4)
      memsave4 = memsize4
   end if

   !----- Initialise some variables. ------------------------------------------------------!
   ivar_type=0
   ierr_getvar=0
   ierr=0
   ifound=0

   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !     The huge select case starts here...                                               !
   !---------------------------------------------------------------------------------------!
   select case (trim(cvar))
   case ('u','ue_avg')
      ivar_type   = 3
      ierr        = RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr        = RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_rotate(nx,ny,nz,a,scr%c,ngrd)
      call RAMS_comp_avgu(nx,ny,nz,a)
      cdname      = 'True zonal wind'
      cdunits     = 'm/s'

   case ('v','ve_avg')
      ivar_type   = 3
      ierr        = RAMS_getvar('VP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr        = RAMS_getvar('UP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_rotate(nx,ny,nz,scr%c,a,ngrd)
      call RAMS_comp_avgv(nx,ny,nz,a)
      cdname      = 'True meridional wind'
      cdunits     = 'm/s'

   case ('w','w_avg')
      ivar_type   = 3
      ierr        = RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_avgw(nx,ny,nz,a)
      cdname      = 'True vertical velocity'
      cdunits     = 'm/s'

   case ('speed')
      ivar_type   = 3
      ierr        = RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr        = RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_speed(nx,ny,nz,a,scr%c)
      cdname      = 'total wind speed'
      cdunits     = 'm/s'


   case('direction')
      ivar_type=3
      ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dir(nx,ny,nz,a,scr%c,ngrd)
      cdname='direction'
      cdunits='deg'

   case ('pvap')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)               ! c is pressure [hPa]
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)    ! a is mixing ratio
      call RAMS_comp_noneg(nx,ny,nz,a)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_pvap(nx,ny,nz,scr%c,a)              ! a is vapour pressure [hPa]
      cdname='water vapour pressure'
      cdunits='hPa'

   case ('alpha')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm) ! c is Exner function
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)  ! a is potential temperature [K]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,a,scr%c)              ! a is temperature [K]
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm) ! d is mixing ratio [kg/kg]
      call RAMS_comp_noneg(nx,ny,nz,scr%d)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)                ! c is pressure [hPa]
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)           ! c is pressure [Pa]
      call RAMS_comp_pvap(nx,ny,nz,scr%c,scr%d)           ! d is vapour pressure [hPa]
      call RAMS_comp_spvol(nx,ny,nz,a,scr%d,scr%c)        ! a is specific volume [m3/kg]
      cdname='specific volume'
      cdunits='m3/kg'

   case ('solenoidx')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)    ! c is Exner function
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm) ! e is pot. temperature [K]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%e,scr%c)             ! e is temperature [K]
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm)    ! d is mixing ratio [kg/kg]
      call RAMS_comp_noneg(nx,ny,nz,scr%d)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)                   ! c is pressure [hPa]
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)              ! c is pressure [Pa]
      call RAMS_comp_pvap(nx,ny,nz,scr%c,scr%d)              ! d is vapour pressure [hPa]
      call RAMS_comp_spvol(nx,ny,nz,scr%e,scr%d,scr%c)       ! e is specific volume [m3/kg]
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)  ! h is topography [m]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_solenoidx(nx,ny,nz,scr%e,scr%c,a,scr%h,ngrd)
      cdname='x-solenoid term'
      cdunits='rad/s2'

   case ('solenoidy')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)    ! c is Exner function
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm) ! e is pot. temperature [K]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%e,scr%c)             ! e is temperature [K]
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm)    ! d is mixing ratio [kg/kg]
      call RAMS_comp_noneg(nx,ny,nz,scr%d)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)                   ! c is pressure [hPa]
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)              ! c is pressure [Pa]
      call RAMS_comp_pvap(nx,ny,nz,scr%c,scr%d)              ! d is vapour pressure [hPa]
      call RAMS_comp_spvol(nx,ny,nz,scr%e,scr%d,scr%c)       ! e is specific volume [m3/kg]
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)  ! h is topography [m]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_solenoidy(nx,ny,nz,scr%e,scr%c,a,scr%h,ngrd)
      cdname='y-solenoid term'
      cdunits='rad/s2'

   case ('relvortx')
      ivar_type=3
      ierr= RAMS_getvar('VP',idim_type,ngrd,a,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,a,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('WP',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%c,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_relvortx(nx,ny,nz,a,scr%c,b,scr%d,ngrd)
      cdname='x-vorticity'
      cdunits='rad/s'

   case ('relvorty')
      ivar_type=3
      ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,a,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('WP',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%c,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_relvorty(nx,ny,nz,a,scr%c,b,scr%d,ngrd)
      cdname='y-vorticity'
      cdunits='rad/s'

   case ('relvortz')
      ivar_type=3 
      ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,a,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%c,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_relvortz(nx,ny,nz,a,scr%c,b,scr%d,ngrd)
      cdname='relative z-vorticity'
      cdunits='rad/s'

   case ('absvortz')
      ivar_type=3
      ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,a,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%c,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_totvortz(nx,ny,nz,a,scr%c,b,scr%d,ngrd)
      cdname='absolute z-vorticity'
      cdunits='rad/s'

   case ('potvortz')
      ivar_type=3
      ierr= RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,a,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%c,1.e-6)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_totvortz(nx,ny,nz,a,scr%c,b,scr%d,ngrd)
      call RAMS_comp_dn0(nx,ny,nz,scr%e,b,scr%c,scr%d,ngrd)

      ierr= RAMS_getvar('THETA',idim_type,ngrd,b,scr%e,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_potvortz(nx,ny,nz,a,b,scr%c,scr%e,scr%d,ngrd)
      cdname='potential z-vorticity'
      cdunits='rad/s'

   case ('pi')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Exner function'
      cdunits='J/(kg K)'

   case ('press')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,a)
      cdname='pressure'
      cdunits='mb'

   case ('theta')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential temp'
      cdunits='K'

   case ('thil')
      ivar_type=3
      ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='ice-liquid potential temp'
      cdunits='K'

   case ('co2')
      ivar_type=3
      if (co2_on) then
         ierr= RAMS_getvar('CO2P',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
      else
         write (unit=*,fmt='(a,1x,es12.5)') '       # Assigning constant CO2P =',co2con(1)
         call ae0(nx*ny*nz,a,co2con(1))
      end if
      cdname='carbon dioxide mixing ratio'
      cdunits='umol/mol'

   case ('dn0')
      ivar_type=3
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%c,b,a,scr%e,ngrd)
      cdname='ref density'

      cdunits='kg/m3'

   case ('pi0')
      ivar_type=3
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,a,b,scr%c,scr%e,ngrd)
      cdname='ref Exner func'
      cdunits='J/(kg K)'

   case ('th0')
      ivar_type=3
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,a,scr%c,scr%e,ngrd)
      cdname='reference virtual potential temp'
      cdunits='K'

   case ('tempk')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,a,scr%c)
      cdname='temperature'
      cdunits='K'

   case ('tempc')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,a,scr%c)
      call RAMS_comp_tempC(nx,ny,nz,1,a)
      cdname='temperature'
      cdunits='C'


   case ('thetae_iv','theiv')
      ivar_type=3
      ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%f,scr%c)
      call RAMS_comp_press(nx,ny,nz,scr%c)
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)
      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)
      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%h,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%h)

      call RAMS_comp_thetaeiv(nx,ny,nz,a,scr%f,scr%c,scr%d,scr%e)
      cdname='Ice-vapour equivt pot temp'
      cdunits='K'

   case ('theta_v')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RTP',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_thetv(nx,ny,nz,a,scr%c,scr%d)
      cdname='virtual pot temp'
      cdunits='K'

   case ('rv')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='vapour mixing ratio'
      cdunits='g/kg'

   case ('cloud')
      ivar_type=3
      ierr= RAMS_getvar('RCP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='cloud mixing ratio'
      cdunits='g/kg'

   case ('rain')
      ivar_type=3
      ierr= RAMS_getvar('RRP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='rain mixing ratio'
      cdunits='g/kg'

   case ('pristine')
      ivar_type=3
      ierr= RAMS_getvar('RPP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='pristine ice mixing ratio'
      cdunits='g/kg'

   case ('snow')
      ivar_type=3
      ierr= RAMS_getvar('RSP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='snow mixing ratio'
      cdunits='g/kg'

   case ('aggregates','agg')
      ivar_type=3
      ierr= RAMS_getvar('RAP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='aggregate mixing ratio'
      cdunits='g/kg'

   case ('graupel')
      ivar_type=3
      ierr= RAMS_getvar('RGP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='graupel mixing ratio'
      cdunits='g/kg'

   case ('hail')
      ivar_type=3
      ierr= RAMS_getvar('RHP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='hail mixing ratio'
      cdunits='g/kg'

   case ('liquid')
      ivar_type=3
      call RAMS_comp_zero(nx,ny,nz,a)
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q6',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracliq(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      endif

      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q7',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracliq(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
       endif

      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='liquid mixing ratio'
      cdunits='g/kg'

   case ('ice')
      ivar_type=3
      call RAMS_comp_zero(nx,ny,nz,a)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q6',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracice(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      endif

      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q7',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracice(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      endif

      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='ice mixing ratio'
      cdunits='g/kg'

   case ('total_cond')
      ivar_type=3
      call RAMS_comp_zero(nx,ny,nz,a)
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='cloud mixing ratio'
      cdunits='g/kg'

   case ('rall')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='vapour + condensed mixing ratio'
      cdunits='g/kg'

   case ('rtp')
      ivar_type=3
      ierr= RAMS_getvar('RTP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='total water mixing ratio'
      cdunits='g/kg'

   case ('dewptk')
      ivar_type=3
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_dewK(nx,ny,nz,a,scr%c,scr%d)
      call RAMS_comp_tempK(nx,ny,nz,a,scr%c)
      cdname='dewpoint temp'
      cdunits='K'

   case ('dewptc')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_dewK(nx,ny,nz,a,scr%c,scr%d)
      call RAMS_comp_tempC(nx,ny,nz,1,a)
      cdname='dewpoint temp'
      cdunits='C'

   case ('rh')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_rh(nx,ny,nz,a,scr%c,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='relative humidity'
      cdunits='pct'

   case ('clear_frac')
      ivar_type=2
      ierr= RAMS_getvar('RV',idim_type,ngrd,b,a,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,a,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,a,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_rh(nx,ny,nz,b,scr%c,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,b)
      
      call cldfraction(nx,ny,nz,a,scr%c,b)
      
      cdname='clear sky fraction'
      cdunits='n/d'

   case ('cloud_concen_mg')
      ivar_type=3
   ! variable 18 is iccp
      ierr= RAMS_getvar('CCP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='cloud concen'
      cdunits='#/mg'

   case ('rain_concen_kg')
      ivar_type=3
      ierr= RAMS_getvar('CRP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='rain concen'
      cdunits='#/kg'

   case ('pris_concen_kg')
      ivar_type=3
      ierr= RAMS_getvar('CPP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='pristine concen'
      cdunits='#/kg'

   case ('snow_concen_kg')
      ivar_type=3
      ierr= RAMS_getvar('CSP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='snow concen'
      cdunits='#/kg'

   case ('agg_concen_kg')
      ivar_type=3
      ierr= RAMS_getvar('CAP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='aggregate concen'
      cdunits='#/kg'

   case ('graup_concen_kg')
      ivar_type=3
      ierr= RAMS_getvar('CGP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='graupel concen'
      cdunits='#/kg'

   case ('hail_concen_kg')
      ivar_type=3
      ierr= RAMS_getvar('CHP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='hail concen'
      cdunits='#/kg'

   case ('cloud_concen_cm3')
      ivar_type=3
      ierr= RAMS_getvar('CCP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='cloud concen'
      cdunits='#/cm3'

   case ('rain_concen_m3')
      ivar_type=3
      ierr= RAMS_getvar('CRP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='rain concen'
      cdunits='#/m3'

   case ('pris_concen_m3')
      ivar_type=3
      ierr= RAMS_getvar('CPP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='pristine concen'
      cdunits='#/m3'

   case ('snow_concen_m3')
      ivar_type=3
      ierr= RAMS_getvar('CSP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='snow concen'
      cdunits='#/m3'

   case ('agg_concen_m3')
      ivar_type=3
      ierr= RAMS_getvar('CAP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='aggregates concen'
      cdunits='#/m3'

   case ('graup_concen_m3')
      ivar_type=3
      ierr= RAMS_getvar('CGP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='graupel concen'
      cdunits='#/m3'

   case ('hail_concen_m3')
      ivar_type=3
      ierr= RAMS_getvar('CHP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='hail concen'
      cdunits='#/m3'

   case ('ccn_concen')
      ivar_type=3
      ierr= RAMS_getvar('CCCNP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)
      cdname='ccnx concen'
      cdunits='#/mg'

   case ('ifn_conc')
      ivar_type=3
      ierr= RAMS_getvar('CIFNP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='CN mixing ratio'
      cdunits='#/kg'

   case ('cloud_diam')
      ivar_type=3
      ierr= RAMS_getvar('RCP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CCP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(1),pwmas(1))
      call RAMS_comp_mults(nx,ny,nz,a,1.e6)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='cloud diam'
      cdunits='microns'

   case ('rain_diam')
      ivar_type=3
      ierr= RAMS_getvar('RRP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CRP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(2),pwmas(2))
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='rain diam'
      cdunits='mm'

   case ('pris_diam')
      ivar_type=3
      ierr= RAMS_getvar('RPP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CPP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      !----- More general case: write habit to anal file for cfmas & pwmas index. ---------!
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(3),pwmas(3))
      call RAMS_comp_mults(nx,ny,nz,a,1.e6)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='pristine diam'
      cdunits='microns'

   case ('snow_diam')
      ivar_type=3
      ierr= RAMS_getvar('RSP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CSP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   ! more general case: write habit to anal file for cfmas & pwmas index
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(4),pwmas(4))
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='snow diam'
      cdunits='mm'

   case ('agg_diam')
      ivar_type=3
      ierr= RAMS_getvar('RAP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CAP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(5),pwmas(5))
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='aggregates diam'
      cdunits='mm'

   case ('graup_diam')
      ivar_type=3
      ierr= RAMS_getvar('RGP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CGP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(6),pwmas(6))
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='graupel diam'
      cdunits='mm'

   case ('hail_diam')
      ivar_type=3
      ierr= RAMS_getvar('RHP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CHP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(7),pwmas(7))
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='hail diam'
      cdunits='mm'

   case ('q2','qrain')
      ivar_type=3
      ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Rain internal energy'
      cdunits='J/kg'

   case ('q6','qgraupel')
      ivar_type=3
      ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Graupel internal energy'
      cdunits='J/kg'

   case ('q7','qhail')
      ivar_type=3
      ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Hail internal energy'
      cdunits='J/kg'

   case ('rain_temp')
      ivar_type=3
      ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_raintemp(nx,ny,nz,a)
      cdname='rain temperature'
      cdunits='K'

   case ('graup_temp')
      ivar_type=3
      ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_qtcpcp(nx,ny,nz,a)
      cdname='graupel temperature'
      cdunits='C'

   case ('hail_temp')
      ivar_type=3
      ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_qtcpcp(nx,ny,nz,a)
      cdname='hail temperature'
      cdunits='C'

   case ('rain_air_tempdif')
      ivar_type=3
      ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_raintemp(nx,ny,nz,a)
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%d,scr%c)
      call RAMS_comp_tempC(nx,ny,nz,1,scr%d)
      call RAMS_comp_subt(nx,ny,nz,a,scr%d)
      cdname='rain-air temp'
      cdunits='K'

   case ('graup_air_tempdf')
      ivar_type=3
      ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_qtcpcp(nx,ny,nz,a)
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%d,scr%c)
      call RAMS_comp_tempC(nx,ny,nz,1,scr%d)
      call RAMS_comp_subt(nx,ny,nz,a,scr%d)
      cdname='graupel-air temp'
      cdunits='K'

   case ('hail_air_tempdif')
      ivar_type=3
      ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_qtcpcp(nx,ny,nz,a)
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%d,scr%c)
      call RAMS_comp_tempC(nx,ny,nz,1,scr%d)
      call RAMS_comp_subt(nx,ny,nz,a,scr%d)
      cdname='hail-air temp'
      cdunits='K'

   case ('graup_fracliq')
      ivar_type=3
      ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_fracliq(nx,ny,nz,a)
      cdname='graupel liq frac'
      cdunits=' '

   case ('hail_fracliq')
      ivar_type=3
      ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_fracliq(nx,ny,nz,a)
      cdname='hail liq frac'
      cdunits=' '

   case ('geo')
      ivar_type=3
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_z(nx,ny,nz,a,scr%c,ngrd)
      cdname='geopotential height'
      cdunits='m'

   case ('tke')
      ivar_type=3
      ierr= RAMS_getvar('TKEP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='turb kinetic energy'
      cdunits='m2/s2'

   case ('thsrc') 
      ivar_type=6
      ierr= RAMS_getvar('THSRC',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz*ncld,a,86400.)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='deep conv heat rate'
      cdunits='K/day'

   case ('rtsrc') 
      ivar_type=6
      ierr= RAMS_getvar('RTSRC',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz*ncld,a,86400.)
      call RAMS_comp_mults(nx,ny,nz*ncld,a,1000.)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='deep conv moist rate'
      cdunits='g/kg/day'

   case ('fthrd') 
      ivar_type=3
      ierr= RAMS_getvar('FTHRD',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,86400.)
      cdname='rad heat rate'
      cdunits='K/day'

   case ('fthrd_lw') 
      ivar_type=3
      ierr= RAMS_getvar('FTHRD_LW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,86400.)
      cdname='rad heat rate (longwave)'
      cdunits='K/day'

   case ('khh') 
      ivar_type=3
      ierr= RAMS_getvar('HKH',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='horiz diffusion coeff'
      cdunits='m2/s'

   case ('khv') 
      ivar_type=3
      ierr= RAMS_getvar('VKH',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='vert diffusion coeff'
      cdunits='m2/s'

   case ('accpr','liqpcp') 
      ivar_type=2
      ierr= RAMS_getvar('ACCPR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      select case (trim(cvar))
      case ('accpr')
         cdname='accum rain'
      case('liqpcp')
         cdname='purely liquid precip'
      end select
      cdunits='kg/m2'

   case ('accpp') 
      ivar_type=2
      ierr= RAMS_getvar('ACCPP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='accum pristine'
      cdunits='kg/m2'

   case ('accps') 
      ivar_type=2
      ierr= RAMS_getvar('ACCPS',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='accum snow'
      cdunits='kg/m2'

   case ('accpa') 
      ivar_type=2
      ierr= RAMS_getvar('ACCPA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='accum aggregates'
      cdunits='kg/m2'

   case ('accpg') 
      ivar_type=2
      ierr= RAMS_getvar('ACCPG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='accum graupel'
      cdunits='kg/m2'

   case ('accph') 
      ivar_type=2
      ierr= RAMS_getvar('ACCPH',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='accum hail'
      cdunits='kg/m2'

   case ('totpcp','totpcp_in','precip','precip_in')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('ACCPR',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPS',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPA',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPG',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPH',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)

      select case (trim(cvar))
      case ('precip','precip_in')
         ierr= RAMS_getvar('ACONPR',idim_type,ngrd,scr%c,b,flnm)
         if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
         cdname='total accum precip'
      case ('totpcp','totpcp_in')
         cdname='total resolved precip'
      end select

      select case (trim(cvar))
      case ('totpcp','precip')
         cdunits='kg/m2'
      case ('precip_in','totpcp_in')
         call RAMS_comp_mults(nx,ny,nz,a,.03937)
         cdunits='in liq'
      end select
      call RAMS_comp_noneg(nx,ny,1,a)

   case ('icepcp')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('ACCPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPS',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPA',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)

      cdname='purely ice precip'
      cdunits='kg/m2'
      call RAMS_comp_noneg(nx,ny,1,a)

   case ('mixpcp')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('ACCPG',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPH',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)

      cdname='mixed (ice/liq) precip'
      cdunits='kg/m2'
      call RAMS_comp_noneg(nx,ny,1,a)

   case ('pcprr')
      ivar_type=2
      ierr= RAMS_getvar('PCPRR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,1,a,3600.)
      cdname='rain precip rate'
      cdunits='mm/hr liq equiv'

   case ('pcprp')
      ivar_type=2
      ierr= RAMS_getvar('PCPRP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,1,a,3600.)
      cdname='pristine precip rate'
      cdunits='mm/hr liq equiv'

   case ('psprs')
      ivar_type=2
      ierr= RAMS_getvar('PCPRS',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,1,a,3600.)
      cdname='snow precip rate'
      cdunits='mm/hr liq equiv'

   case ('pcpra')
      ivar_type=2
      ierr= RAMS_getvar('PCPRA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,1,a,3600.)
      cdname='aggregates precip rate'
      cdunits='mm/hr liq equiv'

   case ('pcprg')
      ivar_type=2
      ierr= RAMS_getvar('PCPRG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='graupel precip rate'
      cdunits='mm/hr liq equiv'

   case ('pcprh')
      ivar_type=2
      ierr= RAMS_getvar('PCPRH',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,1,a,3600.)
      cdname='hail precip rate'
      cdunits='mm/hr liq equiv'

   case ('pcpg')
      ivar_type=2
      ierr= RAMS_getvar('PCPG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='pcpg'
      cdunits='kg/m2'

   case ('qpcpg')
      ivar_type=2
      ierr= RAMS_getvar('QPCPG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='qpcpg'
      cdunits='J/m2'

   case ('dpcpg')
      ivar_type=2
      ierr= RAMS_getvar('DPCPG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='dpdpg'
      cdunits='m'

   case ('pcprate','pcprate_in')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('PCPRR',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRS',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRA',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRG',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRH',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='resolved precip rate'

      select case (trim(cvar))
      case ('pcprate','precipr')
         call RAMS_comp_mults(nx,ny,1,a,3600.)
         cdunits='mm/hr'
      case ('pcprate_in','precipr_in')
         call RAMS_comp_mults(nx,ny,1,a,141.732)
         cdunits='in/hr'
      end select

   case ('conpcp','conprr')
      ivar_type=9
      ierr= RAMS_getvar('CONPRR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,1,a,3600.)
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='convective pcp rate'
      cdunits='mm/hr'

   case ('acccon')
      ivar_type=2
      ierr= RAMS_getvar('ACONPR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='accum convective pcp'
      cdunits='mm'




   case ('cape')
      ivar_type=2
   !- rel hum (e)
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_rh(nx,ny,nz,scr%e,scr%c,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,scr%e)
   !- tempk (d)
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%d,scr%c)
   !- press (c)
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)
   !- cape
      call cape_cine(nx,ny,nz,scr%c,scr%d,scr%e,a,'cape',-9.99e33)

      cdname='cape'
      cdunits='J/kg'

   case ('cine')
      ivar_type=2
   !- rel hum (e)
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_rh(nx,ny,nz,scr%e,scr%c,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,scr%e)
   !- tempk (d)
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%d,scr%c)
   !- press (c)
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)
   !- cape
      call cape_cine(nx,ny,nz,scr%c,scr%d,scr%e,a,'cine',-9.99e33)

      cdname='cine'
      cdunits='J/kg'

   ! Tropopause values
   ![ML - Altura da tropopausa
   case ('ztropop')
      ivar_type=2
      ierr=RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm)  ! e= temperatura potencial
      ierr_getvar = ierr_getvar + ierr
      call acha_ztropop(nx,ny,nz,a,scr%e,ngrd)               ! a= altura da tropopausa
      cdname='Tropopause height'
      cdunits='m - sigmaz'
   !ML]   

   ![ML - Temperatura da tropopausa
   case ('ttropop')
      ivar_type=2
      ierr=RAMS_getvar('THETA',idim_type,ngrd,scr%c,b,flnm)  ! c= temperatura potencial
      ierr_getvar = ierr_getvar + ierr
      ierr=RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm)  ! e= temperatura potencial
      ierr_getvar = ierr_getvar + ierr
      ierr=RAMS_getvar('PI',idim_type,ngrd,scr%d,b,flnm)     ! d= funçao de Exner
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%e,scr%d)                 ! e= temperatura em Kelvin
      call acha_ttropop(nx,ny,nz,a,scr%c,scr%e,ngrd)             ! a= temperatura da tropopausa
      cdname='Tropopause temperature'
      cdunits='K'
   !ML]   

   ![ML - Pressao na tropopausa
   case ('ptropop')
      ivar_type=2
      ierr=RAMS_getvar('THETA',idim_type,ngrd,scr%c,b,flnm)  ! c= temperatura potencial
      ierr_getvar = ierr_getvar + ierr
      ierr=RAMS_getvar('PI',idim_type,ngrd,scr%e,b,flnm)     ! e= funçao de Exner
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%e)                   ! e= pressao
      call acha_ptropop(nx,ny,nz,a,scr%c,scr%e,ngrd)             ! a= temperatura da tropopausa
      cdname='Tropopause pressure'
      cdunits='hPa'
   !ML]   


   ![Marcos Parâmetros da convecção para uso em STILT

   case ('cfxup_deep')
      ivar_type=3
      ierr= RAMS_getvar('CFXUP1',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Conv. upward flux - deep'
      cdunits='kg/m2/s'

   case ('cfxdn_deep')
      ivar_type=3
      ierr= RAMS_getvar('CFXDnx',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_nopos(nx,ny,nz,a)
      cdname='Conv. downward flux - deep'
      cdunits='kg/m2/s'

   case ('cfxup_shal')
      ivar_type=3
      ierr= RAMS_getvar('CFXUP2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Conv. upward flux - shallow'
      cdunits='kg/m2/s'

   case ('efxup_deep')
      ivar_type=3
      ierr= RAMS_getvar('EFXUP1',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Updraft entrainment flux - deep'
      cdunits='kg/m2/s'

   case ('efxdn_deep')
      ivar_type=3
      ierr= RAMS_getvar('EFXDnx',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Downdraft entrainment flux - deep'
      cdunits='kg/m2/s'

   case ('efxup_shal')
      ivar_type=3
      ierr= RAMS_getvar('EFXUP2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Updraft entrainment flux - shallow'
      cdunits='kg/m2/s'

   case ('dfxup_deep')
      ivar_type=3
      ierr= RAMS_getvar('DFXUP1',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Updraft detrainment flux - deep'
      cdunits='kg/m2/s'

   case ('dfxdn_deep')
      ivar_type=3
      ierr= RAMS_getvar('DFXDnx',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Downdraft detrainment flux - deep'
      cdunits='kg/m2/s'

   case ('dfxup_shal')
      ivar_type=3
      ierr= RAMS_getvar('DFXUP2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Updraft detrainment flux - shallow'
      cdunits='kg/m2/s'

   case ('cfxup')
      ivar_type=6
      ierr= RAMS_getvar('CFXUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Convective upward flux'
      cdunits='kg/m2/s'

   case ('cfxdn')
      ivar_type=6
      ierr= RAMS_getvar('CFXDN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_nopos(nx,ny,nz*ncld,a)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Convective downward flux'
      cdunits='kg/m2/s'

   case ('dfxup')
      ivar_type=6
      ierr= RAMS_getvar('DFXUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Detrainment upward flux'
      cdunits='kg/m2/s'

   case ('dfxdn')
      ivar_type=6
      ierr= RAMS_getvar('DFXDN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Detrainment upward flux'
      cdunits='kg/m2/s'

   case ('efxup')
      ivar_type=6
      ierr= RAMS_getvar('EFXUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Entrainment upward flux'
      cdunits='kg/m2/s'

   case ('efxdn')
      ivar_type=6
      ierr= RAMS_getvar('EFXDN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Entrainment upward flux'
      cdunits='kg/m2/s'

   ! Extra turbulence parameters

   case ('tkem')
      ivar_type=3
      ierr= RAMS_getvar('TKEPB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Mean turbulent kinetic energy'
      cdunits='m2/s2'

   case ('ltscale')
      ivar_type=3
      ierr= RAMS_getvar('TL',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Mean Lagrangean time scale'
      cdunits='s'

   case ('sigw')
      ivar_type=3
      ierr= RAMS_getvar('SIGW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Mean vertical velocity standard deviation'
      cdunits='m/s'

   case ('snowdepth')
       ivar_type=2
       ierr= RAMS_getvar('SNOW_DEPTH',idim_type,ngrd,a,b,flnm)
       ierr_getvar = ierr_getvar + ierr
       cdname='Depth of the snow layer'
       cdunits='m'

   case ('pblhgt')
       ivar_type=2
       ierr= RAMS_getvar('PBLHGT',idim_type,ngrd,a,b,flnm)
       ierr_getvar = ierr_getvar + ierr
       cdname='PBL Depth'
       cdunits='m'

   case ('lmo')
       ivar_type=2
       ierr= RAMS_getvar('LMO',idim_type,ngrd,a,b,flnm)
       ierr_getvar = ierr_getvar + ierr
       cdname='Obukhov lenght scale'
       cdunits='m'


   ! Vertically-integrated atmospheric moisture

   case ('vertint_rt','vertint_cond')
      ivar_type=2

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%c,b,scr%d,scr%e,ngrd)

      select case(trim(cvar))
      case ('vertint_rt')
         ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         cdname='vertint total water'
      case ('vertint_cond')
         call RAMS_comp_zero(nx,ny,nz,a)
         cdname='vertint condensate'
      end select

      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_vertint(nx,ny,nz,a,scr%e,ngrd)

      cdunits='mm'


   ! 2D SURFACE HEAT, MOISTURE, MOMENTUM AND RADIATIVE FLUXES

   case ('SFLUX_T')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_T',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='SFLUX_T'
      cdunits='m'

   case ('SFLUX_R')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='SFLUX_R'
      cdunits='m'

   case ('uw')
      ivar_type=2
      ierr= RAMS_getvar('UW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='uw'
      cdunits='m'

   case ('vw')
      ivar_type=2
      ierr= RAMS_getvar('VW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='vw'
      cdunits='m'

   case ('SFLUX_W')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_W',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='SFLUX_W'
      cdunits='m'

   case ('SFLUX_C')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_C',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='SFLUX_C'
      cdunits='m'

   case ('hflxca')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_T',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      call RAMS_comp_mults(nx,ny,1,a,1004.)
      cdname='sfc sens heat flx canopy to atmosphere'
      cdunits='W/m2'

   case ('cflxca')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_C',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      call RAMS_comp_mults(nx,ny,1,a,mmdryi)
      cdname='CO2 flux'
      cdunits='umol/m2/s'

   case ('qwflxca')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      call RAMS_comp_mults(nx,ny,1,a,alvl)
      cdname='water flux from canopy to atmosphere'
      cdunits='W/m2'

   case ('etrans')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
   !                 Divide by water density to get depth and 
   !                   convert units from m/s to mm/hour (3600./1000.)
      call RAMS_comp_mults(nx,ny,1,a,3.6)
      cdname='evapo-transpiration'
      cdunits='mm/hour'

   case ('umom_flx')
      ivar_type=2
      ierr= RAMS_getvar('UW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      cdname='sfc u-momentum flx'
      cdunits='Pa'

   case ('vmom_flx')
      ivar_type=2
      ierr= RAMS_getvar('VW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      cdname='sfc v-momentum flx'
      cdunits='Pa'

   case ('wmom_flx')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_W',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      cdname='sfc w-momentum flx'
      cdunits='Pa'

   case ('bowen')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_T',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('SFLUX_R',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_bowen(nx,ny,1,a,scr%c)
      cdname='bowen ratio'
      cdunits=' '

   case ('cosz')
      ivar_type=2
      ierr= RAMS_getvar('COSZ',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='co-sine of zenith angle'
      cdunits=' '

   case ('zen')
      ivar_type=2
      ierr= RAMS_getvar('COSZ',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_comp_zenith(nx,ny,scr%c,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='zenith angle'
      cdunits='deg'

   case ('rshort')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='sfc SW. rad.'
      cdunits='W/m2'

   case ('rshorttoa')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT_TOP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='TOA SW rad.'
      cdunits='W/m2'

   case ('rshortd')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT_DIFFUSE',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='diffuse sfc. SW rad.'
      cdunits='W/m2'

   case ('rlong')
      ivar_type=2
      ierr= RAMS_getvar('RLONG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='rlong'
      cdunits='W/m2'

   case ('rlongup')
      ivar_type=2
      ierr= RAMS_getvar('RLONGUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='rlongup'
      cdunits='W/m2'

   case ('albedt')
      ivar_type=2
      ierr= RAMS_getvar('ALBEDT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='albedt'
      cdunits=' '

   ! 2D TOPOGRAPHY AND GEOGRAPHIC VALUES

   case ('topoa')
      ivar_type=2
      ierr= RAMS_getvar('TOPTA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='topo'
      cdunits='m'

   case ('topo')
      ivar_type=2
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='topo'
      cdunits='m'

   case ('lon','longitude')
      ivar_type=2
      ierr= RAMS_getvar('GLON',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='longitude'
      cdunits='deg'

   case ('lat','latitude')
      ivar_type=2
      ierr= RAMS_getvar('GLAT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='latitude'
      cdunits='deg'

   case ('mynum')
      ivar_type=2
      ierr= RAMS_getvar('MYNUM',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='node ID'
      cdunits=' '

   ! 2D MISCELLANEOUS FIELDS

   case ('slp_rams')
      ivar_type=2
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_z(nx,ny,nz,scr%c,a,ngrd)

      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_slpress(nx,ny,nz,scr%e,scr%d,scr%c,a)
      cdname='sea level pressure'
      cdunits='mb'

   case ('sea_press')
      ivar_type=2
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_slpmm5(nx,ny,nz,scr%e,scr%d,scr%c,a)
      cdname='sea level pressure;'
      cdunits='mb;'



   case ('sfc_div')
      ivar_type=2
      ierr= RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcdiv(nx,ny,nz,a,ngrd)
      cdname='surface divergence'
      cdunits='1/s'

   ! Special use of sst: acquired for patch #1 even where no water exists

   case ('sst')
      ivar_type=2
      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd   &
           ,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call rams_fill_sst(nx,ny,nsl*npat,nsl,a,scr%c)
      cdname='water temperature'
      cdunits='C'


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! LEAF2 variables section

   ! If want a horiz plot, specify a string like 'tgpatch'; it will
   !   return i,j,ip array.
   ! Specify a new ivar_type, not corresponding to anal file var type.  With
   !   horiz plot, get back into iplt.  If have this var type, don't slice.
   ! Need replacement for rams3to2d because windowing is done in there.
   ! Replacement would window but not slice.
   ! Then, if want xz (vert cross section) have name like tgpatch_vert.
   ! This would return entire 4d array from hvlib.f.
   ! Then we have to slice and window with yet another replacement to rams3to2d.

   ! nkk is the record number, where n is the LEAF field number (1, 2, 3, or 4)
   ! and kk is the k level.
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   case ('pfarea')

      ivar_type = 7
      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      irecind = irecind + irecsize
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='patch fractional area'
      cdunits=''

   case ('soil_z0_p','soil_z0_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case(trim(cvar))
      case ('soil_z0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_Z0',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('soil_z0_p')
         ivar_type = 7
      case ('soil_z0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='soil roughness'
      cdunits='m'

   case ('vtype','vtype_bp')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('vtype_bp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
             ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('LEAF_CLASS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_vegclass(irecsize,1,1,a(irecind))
      select case (trim(cvar))
      case ('vtype')
         ivar_type = 7
      case ('vtype_bp')
         ivar_type = 2
         call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),1,npat  &
            ,a(irecind),a(1),b)
      end select
      cdname='vegetation class'
      cdunits='#'

   case ('ndvi','ndvi_bp')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('ndvi_bp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
             ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_NDVIC',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('ndvi')
         ivar_type = 7
      case ('ndvi_bp')
         ivar_type = 2
         call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),1,npat  &
            ,a(irecind),a(1),b)
      end select

      cdname='ndvi'
      cdunits='#'


   case ('qveg_class_p','qveg_class_bp')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('qveg_class_bp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('DATQ_CLASS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_vegclass(irecsize,1,1,a(irecind))

      select case (trim(cvar))
      case ('qveg_class_p')
         ivar_type = 7
      case ('qveg_class_bp')
         ivar_type = 2
         call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),1,npat  &
            ,a(irecind),a(1),b)
      end select

      cdname='q vegetation class'
      cdunits='#'

   case ('vegfrac','veg_fracarea_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('veg_fracarea_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_FRACAREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegfrac')
         ivar_type = 7
      case ('veg_fracarea_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='vegetation frac area'
      cdunits=''

   case ('land')

      ivar_type = 2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_1minus(nnxp(ngrd),nnyp(ngrd),1,a)
      cdname='land frac area'
      cdunits=''

   case ('agb','vegagb','agb_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('agb_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_AGB',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('agb','vegagb')
         ivar_type = 7
      case ('agb_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='above ground biomass'
      cdunits='kgC/m2'

   case ('lai','veglai','lai_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('lai_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_LAI',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('lai','veglai')
         ivar_type = 7
      case ('lai_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='green leaf area index'
      cdunits=''


   case ('tai','tai_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('tai_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_TAI',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('tai')
         ivar_type = 7
      case ('tai_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname=' total leaf area index'
      cdunits=''

   case ('ustar')
      ivar_type = 7
      ierr = RAMS_getvar('USTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='friction velocity'
      cdunits='m/s'

   case ('tstar')
      ivar_type = 7
      ierr = RAMS_getvar('TSTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential temperature scale'
      cdunits='K'

   case ('rstar')
      ivar_type = 7
      ierr = RAMS_getvar('RSTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a,1000.)
      cdname='water vapour mixing ratio scale'
      cdunits='g/kg'

   case ('cstar')
      ivar_type = 7
      ierr = RAMS_getvar('TSTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='CO2 mixing ratio scale'
      cdunits='umol/mol'

   case ('z0','z0_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('z0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('PATCH_ROUGH',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('z0_p')
         ivar_type = 7
      case ('z0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='roughness'
      cdunits='m'


   case ('net_z0_p','net_z0_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('net_z0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('NET_Z0',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('net_z0_p')
         ivar_type = 7
      case ('net_z0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='net roughness'
      cdunits='m'

   case ('vegz0','vegz0_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('vegz0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select


      ierr = RAMS_getvar('VEG_ROUGH',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegz0')
         ivar_type = 7
      case ('vegz0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='vegetation roughness'
      cdunits='m'

   case ('vegdisp','veg_disp_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat

      select case (trim(cvar))
      case ('veg_disp_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_DISP',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegdisp')
         ivar_type = 7
      case ('veg_disp_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='vegetation displacement height'
      cdunits='m'

   case ('patch_wetind')

      ivar_type = 7
      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      irecind = irecind + irecsize
      ierr = RAMS_getvar('WET_INDEX',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='patch wetness index'
      cdunits=''

   case ('snowlevels')

      ivar_type = 7
      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      irecind = irecind + irecsize
      ierr = RAMS_getvar('KSNOW',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='number of snow levels'
      cdunits='#'

   case ('grnd_mixrat_p','grnd_mixrat_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('grnd_mixrat_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('SFC_RS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),1.e3)

      select case (trim(cvar))
      case ('grnd_mixrat_p')
         ivar_type = 7
      case ('grnd_mixrat_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='ground mixing ratio'
      cdunits='g/kg'

   case ('soil_mixrat_p','soil_mixrat_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('soil_mixrat_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('SOIL_RS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),1.e3)

      select case (trim(cvar))
      case ('soil_mixrat_p')
         ivar_type = 7
      case ('soil_mixrat_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='soil mixing ratio'
      cdunits='g/kg'

   case ('lwater_p','lwater_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('lwater_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('lwater_p')
         ivar_type = 7
      case ('lwater_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

       cdname='leaf water'
      cdunits='kg/m2'



   case ('rvcan','rvcan_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('rvcan_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('CAN_RVAP',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),1.e3)

      select case (trim(cvar))
      case ('rvcan')
         ivar_type = 7
      case ('rvcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='canopy mixing ratio'
      cdunits='g/kg'

   case ('co2can','co2can_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('co2can_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('CAN_CO2',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('co2can')
         ivar_type = 7
      case ('co2can_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='CO2 mixing ratio'
      cdunits='umol/mol'

   case ('gpp_p','gpp')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('gpp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('GPP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('gpp_p')
         ivar_type = 7
      case ('gpp')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Gross Primary Production'
      cdunits='umol/m2/s'

   case ('plresp_p','plresp')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('plresp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('PLRESP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('plresp_p')
         ivar_type = 7
      case ('plresp')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Plant respiration'
      cdunits='umol/m2/s'

   case ('resphet_p','resphet')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('resphet')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('RESPHET',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('resphet_p')
         ivar_type = 7
      case ('resphet')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Heterotrophic respiration'
      cdunits='umol/m2/s'

   case ('hflxgc_p','hflxgc')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('hflxgc')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SENSIBLE_GC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('hflxgc_p')
         ivar_type = 7
      case ('hflxgc')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Sensible heat (gnd -> can)'
      cdunits='W/m2'

   case ('hflxvc_p','hflxvc')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('hflxvc')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SENSIBLE_VC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('hflxvc_p')
         ivar_type = 7
      case ('hflxvc')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Sensible heat (veg -> can)'
      cdunits='W/m2'

   case ('h_p','h')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('h')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SENSIBLE_GC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr = RAMS_getvar('SENSIBLE_VC',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call accum(nnxp(ngrd)*nnyp(ngrd)*npat,a(irecind),scr%c)

      select case (trim(cvar))
      case ('h_p')
         ivar_type = 7
      case ('h')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Sensible heat flux'
      cdunits='W/m2'

   case ('qwflxgc_p','qwflxgc')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('qwflxgc')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('EVAP_GC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qwflxgc_p')
         ivar_type = 7
      case ('qwflxgc')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Latent heat (gnd->can)'
      cdunits='W/m2'

   case ('qwflxvc_p','qwflxvc')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('qwflxvc')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('EVAP_VC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qwflxvc_p')
         ivar_type = 7
      case ('qwflxvc')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Latent heat (veg->can)'
      cdunits='W/m2'

   case ('evap_p','evap')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('evap')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('EVAP_GC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr = RAMS_getvar('EVAP_VC',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call accum(nnxp(ngrd)*nnyp(ngrd)*npat,a(irecind),scr%c)

      select case (trim(cvar))
      case ('evap_p')
         ivar_type = 7
      case ('evap')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Latent heat flux due to evaporation'
      cdunits='W/m2'

   case ('transp_p','transp')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('transp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('TRANSP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('transp_p')
         ivar_type = 7
      case ('transp')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Latent heat flux due to transpiration'
      cdunits='W/m2'

   case ('le_p','le')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('le')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('TRANSP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('EVAP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_accum(nx,ny,npat,a(irecind),scr%c)

      select case (trim(cvar))
      case ('le_p')
         ivar_type = 7
      case ('le')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Latent heat flux'
      cdunits='W/m2'

   case ('tveg','tveg_ps')

      irecind = 1
      irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
      select case (trim(cvar))
      case ('tveg_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select
      
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_theta2temp(nx,ny,npat,scr%e,scr%d)
      call RAMS_comp_tempC(nx,ny,1,npat,scr%e)
      
      ierr = RAMS_getvar('VEG_ENERGY',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('VEG_HCAP',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_tvegc(nx,ny,npat,a(irecind),scr%c,scr%d,scr%e)

      !----- Filling first patch with SST. -------------------------------!
      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call rams_fill_sst(nx,ny,nsl*npat,nsl,a(irecind),scr%e)

      select case (trim(cvar))
      case ('tveg')
         ivar_type = 7
      case ('tveg_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='vegetation temperature'
      cdunits='C'

   case ('tcan','tcan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('tcan_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_comp_theta2temp(nx,ny,npat,a(irecind),scr%c)
      call RAMS_comp_tempC(nx,ny,1,npat,a(irecind))

      select case (trim(cvar))
      case ('tcan')
         ivar_type = 7
      case ('tcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='canopy temperature'
      cdunits='C'


   case ('pcan','pcan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('pcan_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),1.e-2)

      select case (trim(cvar))
      case ('pcan')
         ivar_type = 7
      case ('pcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='canopy press'
      cdunits='hPa'

   case ('thcan','thcan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('thcan_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('thcan')
         ivar_type = 7
      case ('thcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='canopy potential temperature'
      cdunits='K'

   case ('thecan','thecan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('thecan_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_THEIV',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('thecan')
         ivar_type = 7
      case ('thecan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='canopy equiv. pot. temperature'
      cdunits='K'

   case ('rshort_gnd','rshort_gnd_ps')

      irecind = 1
      select case (trim(cvar))
      case ('rshort_gnd_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('RSHORT_GND',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('rshort_gnd')
         ivar_type = 7
      case ('rshort_gnd_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Ground shortwave radiation'
      cdunits='W/m2'

   case ('rlong_gnd','rlong_gnd_ps')

      irecind = 1
      select case (trim(cvar))
      case ('rlong_gnd_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('RLONG_GND',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('rlong_gnd')
         ivar_type = 7
      case ('rlong_gnd_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='Ground longwave radiation'
      cdunits='W/m2'

   case ('snow_depth_p','snow_depth_ps')

      irecind = 1

      select case (trim(cvar))
      case ('snow_depth_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SNOW_DEPTH',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_sum_snowlayers(nnxp(ngrd)*nnyp(ngrd),nzs,npat,a(irecind))

      select case (trim(cvar))
      case ('snow_depth_p')
         ivar_type = 7
      case ('snow_depth_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='snow depth'
      cdunits='m'

   case ('snowcover_p','snowcover_ps')

      irecind = 1
      select case (trim(cvar))
      case ('snowcover_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SNOW_MOIST',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_sum_snowlayers(nnxp(ngrd)*nnyp(ngrd),nzs,npat,a(irecind))

      select case (trim(cvar))
      case ('snow_depth_p')
         ivar_type = 7
      case ('snowcover_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),1,npat,a)
      end select

      cdname='snowcover'
      cdunits='kg/m2'

   case ('sltex','sltex_bp')

      irecind = 1
      select case (trim(cvar))
      case ('sltex_bp')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('sltex')
         ivar_type = 8
      case ('sltex_bp')
         ivar_type = 10
         call RAMS_comp_bigpatch(nnxp(ngrd),nnyp(ngrd),nsl,npat  &
            ,a(irecind),a(1),b)
      end select

      cdname='soil textural class'
      cdunits='#'

   case ('soilq','soilq_ps')

      irecind = 1
      
      select case (trim(cvar))
      case ('soilq_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call get_leaf_soil(nx,ny,nsl,npat,a(irecind),a2)

      select case (trim(cvar))
      case ('soilq')
         ivar_type = 8
      case ('soilq_ps')
         ivar_type = 10
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),nsl,npat,a)
      end select

      cdname='soil q'
      cdunits='J/m3'

   case ('smoist','smoist_ps')

      irecind = 1
      select case (trim(cvar))
      case ('smoist_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr


      select case (trim(cvar))
      case ('smoist')
         ivar_type = 8
         call get_leaf_soil(nx,ny,nsl,npat,a,a2)
      case ('smoist_ps')
         ivar_type = 10
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),nsl,npat,a)
      end select

      cdname='soil moisture'
      cdunits='m3/m3'


   case ('tsoil','tsoil_ps')

      irecind   = 1

      select case (trim(cvar))
      case ('tsoil_ps')
         irecsize  = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_copysst(nx,ny,nsl,a(irecind))

      call RAMS_comp_qwtk(nx,ny,nsl,npat,a(irecind),scr%c,scr%d)


      select case (trim(cvar))
      case ('tsoil')
         call RAMS_comp_tempC(nx,ny,nsl,npat,a)
         call get_leaf_soil(nx,ny,nsl,npat,a,a2)
         ivar_type = 8
      case ('tsoil_ps')
         ivar_type = 10
         call RAMS_comp_patchsum(nnxp(ngrd),nnyp(ngrd),nsl,npat,a)
         call RAMS_comp_tempC(nx,ny,nsl,1,a)
      end select

      cdname='soil/sea temp'
      cdunits='C'

   case ('smfrac','smfrac_ps')

      irecind = 1
      select case (trim(cvar))
      case ('smfrac_ps')
         irecsize = nnxp(ngrd) * nnyp(ngrd) * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select
      irecind = irecind + irecsize
      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      
      
      call rams_comp_slmstf(nx,ny,nsl,npat,a(irecind),scr%c)


      select case (trim(cvar))
      case ('smfrac')
         ivar_type = 8
         call get_leaf_soil(nx,ny,nsl,npat,a,a2)
      case ('smfrac_ps')
         ivar_type = 10
         call RAMS_comp_patchsum_l(nnxp(ngrd),nnyp(ngrd),nsl,npat,a)
      end select

      cdname='soil moisture frac'
      cdunits='m3/m3'

   case ('sfcw_temp')
      ivar_type=2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_ENERGY',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_MASS',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcwmeantemp(nx,ny,nzs,npat,scr%c,scr%d,scr%e,scr%f,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='Pond/snow mean temperature'
      cdunits='C'

   case ('sfcw_mass')
      ivar_type=2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_MASS',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcwinteg(nx,ny,nzs,npat,scr%c,scr%e,scr%f,a)
      call RAMS_comp_noneg(nx,ny,1,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='Pond/snow mass'
      cdunits='kg/m2'

   case ('sfcw_depth')
      ivar_type=2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_DEPTH',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcwinteg(nx,ny,nzs,npat,scr%c,scr%d,scr%f,a)
      call RAMS_comp_noneg(nx,ny,1,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='Pond/snow depth'
      cdunits='m'

   ! CATT

   case ('CO')
      ivar_type=3
      ierr= RAMS_getvar('SCLP001',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
     call RAMS_transf_ppb(nx,ny,nz,a)
      cdname='CO Concentration'
      cdunits='ppb'

      case ('src1')
      ivar_type=3
      ierr= RAMS_getvar('scrsc001',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 1'
      cdunits='kg/m2/day'

      case ('src2')
      ivar_type=3
      ierr= RAMS_getvar('scrsc002',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 2'
      cdunits='kg/m2/day'

      case ('src3')
      ivar_type=3
      ierr= RAMS_getvar('scrsc003',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 3'
      cdunits='kg/m2/day'

      case ('src4')
      ivar_type=3
      ierr= RAMS_getvar('scrsc004',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 4'
      cdunits='kg/m2/day'

      case ('src5')
      ivar_type=3
      ierr= RAMS_getvar('scrsc005',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 5'
      cdunits='kg/m2/day'

      case ('src6')
      ivar_type=3
      ierr= RAMS_getvar('scrsc006',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 6'
      cdunits='kg/m2/day'

       case ('src7')
      ivar_type=3
      ierr= RAMS_getvar('scrsc007',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 7'
      cdunits='kg/m2/day'

      case ('src8')
      ivar_type=3
      ierr= RAMS_getvar('scrsc008',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='emission 5'
      cdunits='kg/m2/day'
     
   case ('COstc')
      ivar_type=3
      ierr= RAMS_getvar('SCLP002',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
     call RAMS_transf_ppb(nx,ny,nz,a)
      cdname='CO Conc. without conv. transp'
      cdunits='ppb'

   case ('COANT')
      ivar_type=3
      ierr= RAMS_getvar('SCLP004',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
      call RAMS_transf_ppb(nx,ny,nz,a)
      cdname='CO Concentration ANTRO'
      cdunits='ppb'

   case ('COTOT')
      ivar_type=3
      ierr= RAMS_getvar('SCLP005',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
     call RAMS_transf_ppb(nx,ny,nz,a)
      cdname='CO Conc ANTRO+BB'
      cdunits='ppb'

   case ('PM25')
      ivar_type=3
      ierr= RAMS_getvar('SCLP003',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
   !air density
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)

      call RAMS_transf_ugm3(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='PM25 Concentration'
      cdunits='ug/m3'


   case ('PMINT')
      ivar_type=2
      ierr= RAMS_getvar('SCLP003',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
   !air density
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d) !Unit: kg[pm25]/m3
      call RAMS_comp_vertint(nx,ny,nz,a,scr%e,ngrd) ! Unit: kg[pm25]/m2
      call RAMS_comp_mults(nx,ny,nz,a,1.e+9)  ! converte de kg/m2 para ug/m2

      cdname='PM25 vert int'
      cdunits='ug/m2'

   ! ------------------ AOT ------------------ 
   ! WAVE / 0.256, 0.280, 0.296, 0.319, 0.335, 0.365, 0.420, 0.482,
   !        0.598, 0.690, 0.762, 0.719, 0.813, 0.862, 0.926, 1.005,
   !        1.111, 1.333, 1.562, 1.770, 2.051, 2.210, 2.584, 3.284,
   !        3.809, 4.292,
   !        4.546, 4.878, 5.128, 5.405, 5.714, 6.061, 6.452, 6.897,
   !        7.407, 8.333, 9.009, 10.309,12.500,13.889,16.667,
   !        20.000, 26.316, 35.714, 62.50                         
   case ('aot256')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,1,a,scr%c)
      cdname='AOT 256nm'
      cdunits=' '

   case ('aot296')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,3,a,scr%c)
      cdname='AOT 296nm'
      cdunits=' '

   case ('aot335')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,5,a,scr%c)
      cdname='AOT 335nm'
      cdunits=' '

   case ('aot420')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,7,a,scr%c)
      cdname='AOT 420nm'
      cdunits=' '

   case ('aot482')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,8,a,scr%c)
      cdname='AOT 482nm'
      cdunits=' '


   case ('aot598')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,9,a,scr%c)
      cdname='AOT 598nm'
      cdunits=' '

   case ('aot690')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,10,a,scr%c)
      cdname='AOT 690nm'
      cdunits=' '

   case ('aot500')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,11,a,scr%c)
      cdname='AOT 500nm'
      cdunits=' '

   case ('aot550')
      ivar_type=2
      ierr= RAMS_getvar('AOT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nwave,12,a,scr%c)
      cdname='AOT 550nm'
      cdunits=' '


   case ('secog')
      ivar_type=2
      ierr= RAMS_getvar('DUM1',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nz,2,a,scr%c)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='GOES-8 ABBA CO emission'
      cdunits='kg/m2/day'


   case ('secod')
      ivar_type=2
      ierr= RAMS_getvar('DUM1',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nz,11,a,scr%c)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='Duncan CO emission'
      cdunits='kg/m2/day'

   case ('secoant')
      ivar_type=2
      ierr= RAMS_getvar('DUM1',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nz,11,a,scr%c)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='Antropogenic CO emission'
      cdunits='kg/m2/day'

   case ('secoe')
      ivar_type=2
      ierr= RAMS_getvar('DUM1',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call D3toD2(nx,ny,nz,14,a,scr%c)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/day para kg/day
      cdname='EDGAR CO emission'
      cdunits='kg/m2/day'


   case ('scco')
      ivar_type=2
      ierr= RAMS_getvar('QSC1',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO emitida'
      cdunits='kg/(m2 day)'

   case ('scpm25')
      ivar_type=2
      ierr= RAMS_getvar('QSC2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de PM25 emitida'
      cdunits='kg/(m2 day)'

   case ('sccofe')
      ivar_type=2
      ierr= RAMS_getvar('QSC3',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO FWB - EDGAR emitida'
      cdunits='kg/(m2 day)'

   case ('sccoae')
      ivar_type=2
      ierr= RAMS_getvar('QSC4',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO AWB - EDGAR emitida'
      cdunits='kg/(m2 day)'

   case ('sccobbe')
      ivar_type=2
      ierr= RAMS_getvar('QSC5',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO BB - EDGAR emitida'
      cdunits='kg/(m2 day)'

   case ('sccod')
      ivar_type=2
      ierr= RAMS_getvar('QSC9',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO Duncan emitida'
      cdunits='kg/(m2 day)'

   case ('sccol')
      ivar_type=2
      ierr= RAMS_getvar('QSC3',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO emitida -logan'
      cdunits='kg/(m2 day)'

   case ('sccoant')
      ivar_type=2
      ierr= RAMS_getvar('QSC9',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Massa de CO emitida -ANTROPO'
      cdunits='kg/(m2 day)'

   case ('pw','pwv')
      ivar_type=2
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !air density
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,b,scr%c,scr%d,scr%e,ngrd) ! d=dens_ar
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)         ! aqui a=rv*dens_ar
      call RAMS_comp_vertint(nx,ny,nz,a,scr%e,ngrd) ! agua em kg/m^2
      call RAMS_comp_mults(nx,ny,nz,a,0.1) !converte para cm = 1 kg/m^2 * 100 cm/m / (1000 kg/m^3 dens_agua)
      cdname='precipitable water vapor'
      cdunits='cm'



   ! ------------------------ Stilt-RAMS coupling------------
   case ('afxu')
      ivar_type=3
      ierr= RAMS_getvar('AFXU',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='adv u flux'
      cdunits='kg/m^2s'

   case ('afxub')
      ivar_type=3
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('AFXUB',idim_type,ngrd,a,b,flnm)
      cdname='averaged adv u flux'
      cdunits='kg/m^2s'

   case ('afxv')
      ivar_type=3
      ierr= RAMS_getvar('AFXV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='adv v flux'
      cdunits='kg/m^2s'

   case ('afxvb')
      ivar_type=3
      ierr= RAMS_getvar('AFXVB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged adv v flux'
      cdunits='kg/m^2s'

   case ('afxw')
      ivar_type=3
      ierr= RAMS_getvar('AFXW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='adv w flux'
      cdunits='kg/m^2s'

   case ('afxwb')
      ivar_type=3
      ierr= RAMS_getvar('AFXWB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged adv W flux'
      cdunits='kg/m^2s'

   case ('sigwb')
      ivar_type=3
      ierr= RAMS_getvar('SIGWB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged sigma W'
      cdunits='m/s'

   case ('tlb')
      ivar_type=3
      ierr= RAMS_getvar('TLB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged Lagr timescale'
      cdunits='s'

   case ('tl')
      ivar_type=3
      ierr= RAMS_getvar('TL',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Lagr timescale'
      cdunits='s'

   case ('tkeb')
      ivar_type=3
      ierr= RAMS_getvar('TKEPB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='average turb kinetic energy'
      cdunits='m2/s2'



   !------------Grell cumulus scheme --------------------------

   case ('wdm1')
      ivar_type=2
      ierr= RAMS_getvar('wetdep001',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
      cdname='Wet deposition mass tracer 1'
      cdunits='kg/m2'


   case ('wdm3')
      ivar_type=2
      ierr= RAMS_getvar('wetdep003',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)  ! converte de mg/kg para kg/kg
      cdname='Wet deposition mass tracer 3'
      cdunits='kg/m2'


   case ('cuprliq')
      ivar_type=6
      ierr= RAMS_getvar('CUPRLIQ',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call RAMS_comp_mults(nx,ny,nz*ncld,a,1000.)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Conv. water mixing ratio'
      cdunits='g/kg'

   case ('cuprice')
      ivar_type=6
      ierr= RAMS_getvar('CUPRICE',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call RAMS_comp_mults(nx,ny,nz*ncld,a,1000.)
      call get_cumulus(nx,ny,nz,ncld,a,a6)
      cdname='Conv. water mixing ratio'
      cdunits='g/kg'

   case ('areadn')
      ivar_type=9
      ierr= RAMS_getvar('AREADN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Downdraft relative area'
      cdunits=''

   case ('areaup')
      ivar_type=9
      ierr= RAMS_getvar('AREAUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Updraft relative area'
      cdunits=''


   case ('ierr')
      ivar_type=9
      ierr= RAMS_getvar('XIERR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Grell''s error flag:'
      cdunits=' '

   case ('upmf')
      ivar_type=9
      ierr= RAMS_getvar('UPMF',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='updraft mass flux'
      cdunits='kg/m2/s'

   case ('dnmf')
      ivar_type=9
      ierr= RAMS_getvar('DNMF',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='downdraft mass flux'
      cdunits='kg/m2/s'

   case ('upmx')
      ivar_type=9
      ierr= RAMS_getvar('UPMX',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential updraft mass flux'
      cdunits='kg/m2/s'

   case ('dnmx')
      ivar_type=9
      ierr= RAMS_getvar('DNMX',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential downdraft mass flux'
      cdunits='kg/m2/s'

   case ('wdndraft')
      ivar_type=9
      ierr= RAMS_getvar('WDNDRAFT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='downdraft velocity at origin of downdrafts'
      cdunits='m/s'

   case ('wupdraft')
      ivar_type=9
      ierr= RAMS_getvar('WUPDRAFT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='updraft velocity at origin of updrafts'
      cdunits='m/s'

   case ('wbuoymin')
      ivar_type=9
      ierr= RAMS_getvar('WBUOYMIN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='minimum velocity to reach LCL'
      cdunits='m/s'

   case ('zklnb')
      ivar_type=9
      ierr= RAMS_getvar('ZKLNB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Level of neutral buoyancy'
      cdunits='m'

   case ('zklfc')
      ivar_type=9
      ierr= RAMS_getvar('ZKLFC',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Level of free convection'
      cdunits='m'

   case ('zklcl')
      ivar_type=9
      ierr= RAMS_getvar('ZKLCL',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Lifting condensation level'
      cdunits='m'

   case ('zklod')
      ivar_type=9
      ierr= RAMS_getvar('ZKLOD',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Height of origin of downdraft'
      cdunits='m'

   case ('zklou')
      ivar_type=9
      ierr= RAMS_getvar('ZKLOU',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Height of origin of updraft'
      cdunits='m'

   case ('zkdet')
      ivar_type=9
      ierr= RAMS_getvar('ZKDET',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Top of downdraft detrainment'
      cdunits='m'

   case ('edt')
      ivar_type=9
      ierr= RAMS_getvar('EDT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Downdraft/updraft ratio'
      cdunits=' '

   case ('aadn')
      ivar_type=9
      ierr= RAMS_getvar('AADN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Dndraft cloud work function'
      cdunits='J/kg'

   case ('aaup')
      ivar_type=9
      ierr= RAMS_getvar('AAUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Updraft cloud work function'
      cdunits='J/kg'

   !---------------------------------------------------------------------------------------!
   !     Variables that are reduced to a given height, using the Monin-Obukhov similarity  !
   ! theory.  We also re-compute the "perceived" stars here in case we are running with    !
   ! ED-2.1.                                                                               !
   !---------------------------------------------------------------------------------------!
   case ('tempc2m','theta2m','rv2m','tdewc2m','rhum2m','co22m','u10m','rib_ps','zeta_ps'    &
        ,'ustar_ps','tstar_ps','rstar_ps','cstar_ps')
      ivar_type = 2

      !----- Topography. ------------------------------------------------------------------!
      ierr = RAMS_getvar('TOPT',idim_type,ngrd,scr%c,b,flnm)  ! c = topography

      !----- Winds. -----------------------------------------------------------------------!
      ierr = RAMS_getvar('UP',idim_type,ngrd,scr%d,b,flnm)    ! d = zonal wind
      ierr_getvar = ierr_getvar + ierr
      ierr=RAMS_getvar('VP',idim_type,ngrd,scr%e,b,flnm)      ! e = meridional wind
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_speed(nx,ny,nz,scr%d,scr%e)              ! d = wind magnitude
      call RAMS_flush_to_zero(nx,ny,nz,npat,scr%e,ubmin)      ! d = wind magnitude

      !----- Atmospheric properties. ------------------------------------------------------!
      ierr = RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm) ! e = potential temperature.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('RV',idim_type,ngrd,scr%f,b,flnm)    ! f = H2O mixing ratio
      ierr_getvar = ierr_getvar + ierr
      if (co2_on) then
         ierr= RAMS_getvar('CO2P',idim_type,ngrd,scr%g,b,flnm)
         ierr_getvar = ierr_getvar + ierr
      else
         write (unit=*,fmt='(a,1x,es12.5)') '       # Assigning constant CO2P =',co2con(1)
         call ae0(nx*ny*nz,scr%g,co2con(1))
      end if

      !----- Roughness. -------------------------------------------------------------------!
      ierr = RAMS_getvar('PATCH_ROUGH', idim_type,ngrd,scr%h,b,flnm) ! h = patch roughness.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('RIBULK',idim_type,ngrd,scr%t,b,flnm)       ! t = bulk Ri
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('ZETA',idim_type,ngrd,scr%u,b,flnm)         ! u = z/L
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%m,b,flnm)   ! m = patch area.
      ierr_getvar = ierr_getvar + ierr

      !----- Canopy air space (CAS) properties. ----------------------------------------------!
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,scr%n,b,flnm) ! n = CAS potential temp.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_RVAP',idim_type,ngrd,scr%o,b,flnm)  ! o = CAS mixing ratio.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,scr%s,b,flnm)  ! s = canopy pressure
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_CO2',idim_type,ngrd,scr%w,b,flnm)   ! w = CAS CO2 mixing r.
      ierr_getvar = ierr_getvar + ierr

      !----- Characteristic scales (aka stars). ----------------------------------------------!
      ierr = RAMS_getvar('USTAR',idim_type,ngrd,scr%p,b,flnm)     ! p = ustar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%p,ustmin)          ! p = ustar
      ierr = RAMS_getvar('TSTAR',idim_type,ngrd,scr%q,b,flnm)     ! q = tstar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%q,1.e-6)           ! q = tstar
      ierr = RAMS_getvar('RSTAR',idim_type,ngrd,scr%r,b,flnm)     ! r = rstar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%r,1.e-6)           ! r = rstar
      ierr = RAMS_getvar('CSTAR',idim_type,ngrd,scr%x,b,flnm)     ! x = cstar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%x,1.e-6)           ! x = cstar
      select case (trim(cvar))
      case ('theta2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'THET',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Potential temperature at 2m AGL'
         cdunits = 'K'

      case ('tempc2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TEMP',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_tempC(nx,ny,1,1,a)

         cdname  = 'Temperature at 2m AGL'
         cdunits = 'C'

      case ('tdewc2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TDEW',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_tempC(nx,ny,1,1,a)

         cdname  = 'Dew/frost point at 2m AGL'
         cdunits = 'C'

      case ('rv2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'RVAP',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_mults(nx,ny,1,a,1000.)

         cdname  = 'Vapour mixing ratio at 2m AGL'
         cdunits = 'g/kg'

      case ('co22m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'CO_2',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'CO2 mixing ratio at 2m AGL'
         cdunits = 'umol/mol'

      case ('rhum2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TDEW',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TEMP',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,scr%v)
         call RAMS_comp_relhum(nx,ny,1,a,scr%v)
         call RAMS_comp_mults(nx,ny,1,a,100.)

         cdname  = 'Relative humidity at 2m AGL'
         cdunits = '%'

      case ('zeta_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'ZETA',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Dimensionless height'
         cdunits = '---'

      case ('rib_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'RICH',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Bulk Richardson number'
         cdunits = '---'

      case ('u10m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'WIND',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,10.,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Wind speed at 10m AGL'
         cdunits = 'm/s'

      case ('ustar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'USTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Friction velocity'
         cdunits = 'm/s'

      case ('tstar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TSTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Potential temperature characteristic scale'
         cdunits = 'K'

      case ('rstar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'RSTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Vapour mixing ratio characteristic scale'
         cdunits = 'g/kg'
         call RAMS_comp_mults(nx,ny,1,a,1.e3)

      case ('cstar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'CSTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'CO2 mixing ratio characteristic scale'
         cdunits = 'umol/mol'

      end select


   case default

      write (unit=*,fmt='(2(a,1x))') '       # Variable name not found in RAMS_varlib -'   &
                                    ,trim(cvar)
      ivar_type = 0

   end select

   if (ierr_getvar > 0) then
     write (unit=*,fmt='(3(a,1x))') '       # WARNING! Not all the variables needed for'   &
                                   ,trim(cvar),' are available...' 
     ivar_type=0
   end if

   return
end subroutine RAMS_varlib
!==========================================================================================!
!==========================================================================================!
