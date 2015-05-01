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
   use somevars
   
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                          , intent(out)   :: nfile
   character(len=str_len), dimension(maxfiles)      , intent(inout) :: fnames
   real                  , dimension(nzpmax,maxgrds), intent(inout) :: dep_zlev
   integer                                          , intent(inout) :: iep_np
   integer                                          , intent(inout) :: iep_nc
   integer                                          , intent(inout) :: iep_ng
   integer                                          , intent(inout) :: iep_ngrids
   integer               , dimension(       maxgrds), intent(inout) :: iep_nx
   integer               , dimension(       maxgrds), intent(inout) :: iep_ny
   integer               , dimension(       maxgrds), intent(inout) :: iep_nz
   !----- Internal variables. -------------------------------------------------------------!
   character(len=str_len)                                           :: file_prefix
   character(len=str_len)                                           :: fpref
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
            nnxp(n)   = mynnxp(n)
            nnyp(n)   = mynnyp(n)
            nnzp(n)   = mynnzp(n)
            maxmem    = max(maxmem, mynnxp(n)*mynnyp(n)                                    &
                                  * max(mynnzp(n),npatch*nzg,mynnzp(n)*nclouds))
            iep_nx(n) = mynnxp(n)
            iep_ny(n) = mynnyp(n)
            iep_nz(n) = mynnzp(n)
            iep_ng    = nzg
            iep_np    = npatch
            iep_nc    = nclouds
            do nn=1,mynnzp(n)
               dep_zlev(nn,n)=myztn(nn,n)
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
         nfgpnts(1,ng,nfn) = mynnxp(ng)
         nfgpnts(2,ng,nfn) = mynnyp(ng)
         nfgpnts(3,ng,nfn) = mynnzp(ng)
         nfgpnts(4,ng,nfn) = nzg
         fdelx(ng,nfn)     = DELTAXN(NG)
         fdely(ng,nfn)     = DELTAYN(NG)

         do k=1,mynnzp(ng)
            flevels(k,ng,nfn) = myztn(k,ng)
         end do
      end do

      httop = myzmn(mynnzp(1)-1,1)

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
   use somevars  , only : myiyear1  & ! intent(in)
                        , myimonth1 & ! intent(in)
                        , myidate1  ! ! intent(in)
   use rpost_coms, only : iyear1    & ! intent(in)
                        , imonth1   & ! intent(in)
                        , idate1    ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)  :: nfl
   integer, intent(out) :: iyear
   integer, intent(out) :: imonth
   integer, intent(out) :: idate
   integer, intent(out) :: ihour
   integer, intent(out) :: imin
   !---------------------------------------------------------------------------------------!
   
   iyear1  = myiyear1
   imonth1 = myimonth1
   idate1  = myidate1
   iyear   = myiyear1
   imonth  = myimonth1
   idate   = myidate1
   ihour   = int(float(iftimes(nfl))/10000.)
   imin    = int(float(iftimes(nfl)-10000*ihour)/100.)
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
   use rpost_dims, only : str_len     ! ! intent(inout)
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
   character(len=str_len)                        :: flng
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
   !case ('M','B')
   !   string = trim(stringg)//'M'
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
subroutine ep_getvar(cvar,nx,ny,nz,ng,fn,cdname,cdunits,itype,npatch,nclouds,nzg)
   use rout_coms, only : rout      & ! intent(inout)
                       , rout_vars ! ! variable type
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)            , intent(in)    :: cvar
   character(len=*)            , intent(in)    :: fn
   character(len=*)            , intent(out)   :: cdname
   character(len=*)            , intent(out)   :: cdunits
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: ng
   integer                     , intent(out)   :: itype
   integer                     , intent(in)    :: npatch
   integer                     , intent(in)    :: nzg
   integer                     , intent(in)    :: nclouds
   !---------------------------------------------------------------------------------------!


   !----- Load the variable. --------------------------------------------------------------!
   call RAMS_varlib(cvar,nx,ny,nz,nzg,npatch,nclouds,ng,fn,cdname,cdunits,itype            &
                   ,rout(ng)%abuff,rout(ng)%bbuff)
   !---------------------------------------------------------------------------------------!


   !----- Copy to the appropriate scratch. ------------------------------------------------!
   select case (itype)
   case (2)
      call atob(nx*ny,rout(ng)%abuff,rout(ng)%r2)
      call RAMS_show_range(nx,ny,1,1,rout(ng)%r2,cvar,'EP_getvar')
   case (3)
      call atob(nx*ny*nz,rout(ng)%abuff,rout(ng)%r3)
      call RAMS_show_range(nx,ny,nz,1,rout(ng)%r3,cvar,'EP_getvar')
   case (6)
      call atob(nx*ny*nz*nclouds,rout(ng)%abuff,rout(ng)%r6)
      call RAMS_show_range(nx,ny,nz,nclouds,rout(ng)%r6,cvar,'EP_getvar')
   case (7)
      call atob(nx*ny*npatch,rout(ng)%abuff,rout(ng)%r7)
      call RAMS_show_range(nx,ny,npatch,1,rout(ng)%r7,cvar,'EP_getvar')
   case (8)
      call atob(nx*ny*nzg*npatch,rout(ng)%abuff,rout(ng)%r8)
      call RAMS_show_range(nx,ny,nzg,npatch,rout(ng)%r8,cvar,'EP_getvar')
   case (9)
      call atob(nx*ny*nclouds,rout(ng)%abuff,rout(ng)%r9)
      call RAMS_show_range(nx,ny,nclouds,1,rout(ng)%r9,cvar,'EP_getvar')
   case (10)
      call atob(nx*ny*nzg,rout(ng)%abuff,rout(ng)%r10)
      call RAMS_show_range(nx,ny,nzg,1,rout(ng)%r10,cvar,'EP_getvar')
   end select
   return
end subroutine ep_getvar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ep_setdate(iyear1,imonth1,idate1,strtim,itrec)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer              , intent(in)  :: iyear1
   integer              , intent(in)  :: imonth1
   integer              , intent(in)  :: idate1
   real                 , intent(in)  :: strtim
   integer, dimension(6), intent(out) :: itrec
   !---------------------------------------------------------------------------------------!

   itrec(1) = iyear1
   itrec(2) = imonth1
   itrec(3) = idate1
   itrec(4) = int(mod(strtim,24.))
   itrec(5) = int(mod(strtim,1.)*60)
   itrec(6) = int(mod( (strtim) *3600.,60.))

   return
end subroutine ep_setdate
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine dumps the array to the output binary file (gra file).               !
!------------------------------------------------------------------------------------------!
subroutine ep_putvar(iunit,nxp,nyp,nzp,xa,xz,ya,yz,za,zz,array3d,irec,cvar)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: iunit
   integer                     , intent(in)    :: nxp
   integer                     , intent(in)    :: nyp
   integer                     , intent(in)    :: nzp
   integer                     , intent(in)    :: xa
   integer                     , intent(in)    :: xz
   integer                     , intent(in)    :: ya
   integer                     , intent(in)    :: yz
   integer                     , intent(in)    :: za
   integer                     , intent(in)    :: zz
   real, dimension(nxp,nyp,nzp), intent(in)    :: array3d
   integer                     , intent(inout) :: irec
   character(len=*)            , intent(in)    :: cvar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real, dimension(nxp,nyp)                    :: mat
   character(len=99)                           :: epv
   !---------------------------------------------------------------------------------------!
   do z=za,zz
      do y=1,nyp
         do x=1,nxp
            mat(x,y) = array3d(x,y,z)
         end do
      end do

      irec=irec+1
      write (unit=iunit,rec=irec) ((mat(x,y),x=xa,xz),y=ya,yz)

      write(epv,fmt='(a,i3.3)') 'EP_putvar_LEV=',z
      call RAMS_show_range(nxp,nyp,1,1,mat,cvar,epv)

   end do


   return
end subroutine ep_putvar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine array_interpol(ng,nxg,nyg,nxr,nyr,iinf,jinf,rmi,proj)
   use rout_coms , only : undefflg   ! ! intent(in)
   use misc_coms , only : glong      & ! intent(in)
                        , glatg      ! & ! intent(in)! intent(in)
   use somevars  , only : myxtn      & ! intent(in)
                        , myytn      & ! intent(in)
                        , mydeltaxn  & ! intent(in)
                        , mydeltayn  & ! intent(in)
                        , mypolelat  & ! intent(in)
                        , mypolelon  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)  :: ng
   integer                           , intent(in)  :: nxg
   integer                           , intent(in)  :: nyg
   integer                           , intent(in)  :: nxr
   integer                           , intent(in)  :: nyr
   integer     , dimension(nxg,nyg)  , intent(out) :: iinf
   integer     , dimension(nxg,nyg)  , intent(out) :: jinf
   real(kind=4), dimension(nxg,nyg,4), intent(out) :: rmi
   character(len=*)                  , intent(in)  :: proj
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: i
   integer                                         :: j
   integer                                         :: l
   integer                                         :: i1
   integer                                         :: i2
   integer                                         :: j1
   integer                                         :: j2
   integer                                         :: ix
   integer                                         :: iy
   real(kind=4)                                    :: x
   real(kind=4)                                    :: y
   !---------------------------------------------------------------------------------------!


   select case (trim(proj))
   case ('no')
      !------ Nothing to be done here. ----------------------------------------------------!
      return
      !------------------------------------------------------------------------------------!

   case ('yes')
      !------------------------------------------------------------------------------------!
      !    Build interpolation array.  First, we assign undefined flag everywhere, later   !
      ! we will overwrite this for points inside the domain.                               !
      !------------------------------------------------------------------------------------!
      do i=1,nxg
         do j=1,nyg
            iinf (i,j) = 1
            jinf (i,j) = 1
            do l=1,4
               rmi(i,j,l) = undefflg
            end do
         end do
      end do
      !------------------------------------------------------------------------------------!


      !----- Find the position of the GraDS grid in RAMS. ---------------------------------!
      xloop: do i=1,nxg
         yloop: do j=1,nyg
            call ge_to_xy(mypolelat,mypolelon,glong(i),glatg(j),x,y)

            !----- Skip points outside the domain. ----------------------------------------!
            if (x < myxtn(1,ng) .or. x > myxtn(nxr,ng)) cycle yloop
            if (y < myytn(1,ng) .or. y > myytn(nyr,ng)) cycle yloop
            !------------------------------------------------------------------------------!


            !----- Find x and indices. ----------------------------------------------------!
            xfind: do ix=1,nxr
               if (x <= myxtn(ix,ng)) exit xfind
            end do xfind
            i1        = max(1,ix - 1)
            i2        = ix
            iinf(i,j) = i1
            !------------------------------------------------------------------------------!

            !----- Find y and indices. ----------------------------------------------------!
            yfind: do iy=1,nyr
               if (y <= myytn(iy,ng)) exit yfind
            end do yfind
            j1        = max(1,iy - 1)
            j2        = iy
            jinf(i,j) = j1
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Find matrix components.                                                   !
            !------------------------------------------------------------------------------!
            rmi (i,j,1) = ( x - myxtn(i1,ng) ) / mydeltaxn(ng)
            rmi (i,j,2) = 1.0 - rmi(i,j,1)
            rmi (i,j,3) = ( y - myytn(j1,ng) ) / mydeltayn(ng)
            rmi (i,j,4) = 1.0 - rmi(i,j,3)
            !------------------------------------------------------------------------------!
         end do yloop
         !---------------------------------------------------------------------------------!
      end do xloop
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!
   return
end Subroutine array_interpol
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine proj_rams_to_grads(vp,n,nxr,nyr,nzz,nxg,nyg,rmi,iinf,jinf,this,thisgrads,proj)
   use rout_coms, only : maxnormal  & ! intent(in)
                       , undefflg   ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)                        , intent(in)    :: vp
   integer                                 , intent(in)    :: n
   integer                                 , intent(in)    :: nxr
   integer                                 , intent(in)    :: nyr
   integer                                 , intent(in)    :: nzz
   integer                                 , intent(in)    :: nxg
   integer                                 , intent(in)    :: nyg
   real            , dimension(nxg,nyg,4)  , intent(in)    :: rmi
   integer         , dimension(nxg,nyg)    , intent(in)    :: iinf
   integer         , dimension(nxg,nyg)    , intent(in)    :: jinf
   real            , dimension(nxr,nyr,nzz), intent(in)    :: this
   real            , dimension(nxg,nyg,nzz), intent(inout) :: thisgrads
   character(len=*)                        , intent(in)    :: proj
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                 :: i
   integer                                                 :: j
   integer                                                 :: k
   integer                                                 :: i1
   integer                                                 :: i2
   integer                                                 :: j1
   integer                                                 :: j2
   real                                                    :: r1
   real                                                    :: r2
   real                                                    :: r3
   real                                                    :: r4
   real                                                    :: vy3
   real                                                    :: vy4
   logical         , dimension(nxr,nyr,nzz)                :: weird
   logical                                                 :: weird3
   logical                                                 :: weird4
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Check whether to use lon-lat projection or not.                                   !
   !---------------------------------------------------------------------------------------!
   select case(trim(proj))
   case ('no')
      !----- Ensure that the grads domain is exactly the same as the RAMS one. ------------!
      if (nxg /= nxr .or. nyg /= nyr) then
         call abort_run  ('Projection with problems...','proj_rams_to_grads'               &
                         ,'rpost_main.f90')
      end if
      call atob(nxr*nyr*nzz,this,thisgrads)
      !------------------------------------------------------------------------------------!

   case default

      weird(:,:,:) = abs(this(:,:,:)) > maxnormal

      !------------------------------------------------------------------------------------!
      !     Decide whether to interpolate or use nearest neighbour.   Numerical fields     !
      ! should be interpolated, whereas categorical variables should use nearest           !
      ! neighbour.                                                                         !
      !------------------------------------------------------------------------------------!
      select case (trim(vp))
      case ('mynum','vtype','vtype_bp','scolour','scolour_bp','snowlevels','sltex'         &
           ,'sltex_bp','ierr')
         !---------------------------------------------------------------------------------!
         !    Categorical variable: use nearest neighbour.                                 !
         !---------------------------------------------------------------------------------!
         do i=1,nxg
            do j=1,nyg
               r1 = rmi(i,j,1)
               r2 = rmi(i,j,2)
               r3 = rmi(i,j,3)
               r4 = rmi(i,j,4)
               i1 = iinf(i,j)
               i2 = min(i1+1,nxr)
               j1 = jinf(i,j)
               j2 = min(j1+1,nyr)


               do k=1,nzz
                  !------------------------------------------------------------------------!
                  !     Look for "row 3" values and weights, and skip undefined values.    !
                  !------------------------------------------------------------------------!
                  if (weird(i1,j1,k) .and. weird(i2,j1,k)) then
                     vy3 = undefflg
                  elseif (weird(i1,j1,k)) then
                     vy3 = this(i2,j1,k)
                  elseif (weird(i2,j1,k)) then
                     vy3 = this(i1,j1,k)
                  elseif (r1 <= r2) then
                     vy3 = this(i1,j1,k)
                  else
                     vy3 = this(i2,j1,k)
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Look for "row 4" values and weights, and skip undefined values.    !
                  !------------------------------------------------------------------------!
                  if (weird(i1,j2,k) .and. weird(i2,j2,k)) then
                     vy4 = undefflg
                  elseif (weird(i1,j2,k)) then
                     vy4 = this(i2,j2,k)
                  elseif (weird(i2,j2,k)) then
                     vy4 = this(i1,j2,k)
                  elseif (r1 <= r2) then
                     vy4 = this(i1,j2,k)
                  else
                     vy4 = this(i2,j2,k)
                  end if
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Scale the values using only the good values.                       !
                  !------------------------------------------------------------------------!
                  weird3 = abs(vy3) > maxnormal
                  weird4 = abs(vy4) > maxnormal
                  if (weird3 .and. weird4) then
                     thisgrads(i,j,k) = undefflg
                  elseif (weird3) then
                     thisgrads(i,j,k) = vy4
                  elseif (weird4) then
                     thisgrads(i,j,k) = vy3
                  elseif (r3 <= r4) then
                     thisgrads(i,j,k) = vy3
                  else
                     thisgrads(i,j,k) = vy4
                  end if
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!


      case default
         !---------------------------------------------------------------------------------!
         !    Numerical field: interpolate.                                                !
         !---------------------------------------------------------------------------------!
         do i=1,nxg
            do j=1,nyg
               r1 = rmi(i,j,1)
               r2 = rmi(i,j,2)
               r3 = rmi(i,j,3)
               r4 = rmi(i,j,4)
               i1 = iinf(i,j)
               i2 = min(i1+1,nxr)
               j1 = jinf(i,j)
               j2 = min(j1+1,nyr)


               do k=1,nzz
                  !------------------------------------------------------------------------!
                  !     Look for "row 3" values and weights, and skip undefined values.    !
                  !------------------------------------------------------------------------!
                  if (weird(i1,j1,k) .and. weird(i2,j1,k)) then
                     vy3 = undefflg
                  elseif (weird(i1,j1,k)) then
                     vy3 = this(i2,j1,k)
                  elseif (weird(i2,j1,k)) then
                     vy3 = this(i1,j1,k)
                  else
                     vy3 = this(i1,j1,k) * (1. - r1) + this(i2,j1,k) * (1. - r2)
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Look for "row 4" values and weights, and skip undefined values.    !
                  !------------------------------------------------------------------------!
                  if (weird(i1,j2,k) .and. weird(i2,j2,k)) then
                     vy4 = undefflg
                  elseif (weird(i1,j2,k)) then
                     vy4 = this(i2,j2,k)
                  elseif (weird(i2,j2,k)) then
                     vy4 = this(i1,j2,k)
                  else
                     vy4 = this(i1,j2,k) * (1. - r1) + this(i2,j2,k) * (1. - r2)
                  end if
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Scale the values using only the good values.                       !
                  !------------------------------------------------------------------------!
                  weird3 = abs(vy3) > maxnormal
                  weird4 = abs(vy4) > maxnormal
                  if (weird3 .and. weird4) then
                     thisgrads(i,j,k) = undefflg
                  elseif (weird3) then
                     thisgrads(i,j,k) = vy4
                  elseif (weird4) then
                     thisgrads(i,j,k) = vy3
                  else
                     thisgrads(i,j,k) = vy3 * (1. - r3) + vy4 * (1. - r4)
                  end if
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!

      !------ Show range after interpolation or regridding. -------------------------------!
      call RAMS_show_range(nxg,nyg,nzz,1,thisgrads,vp,'RAMS_2_GRADS')
      !------------------------------------------------------------------------------------!
   end select
   return
end subroutine proj_rams_to_grads
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ge_to_xy(polelat,polelon,xlon,xlat,x,y)
   use rconstants, only : erad   & ! intent(in)
                        , pio180 ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in)  :: polelat
   real, intent(in)  :: polelon
   real, intent(in)  :: xlon
   real, intent(in)  :: xlat
   real, intent(out) :: x
   real, intent(out) :: y
   !----- Local variables. ----------------------------------------------------------------!
   real              :: b
   real              :: f
   real              :: xlonrad
   real              :: xlatrad
   real              :: plonrad
   real              :: platrad
   !---------------------------------------------------------------------------------------!


   !----- Convert coordinates to radians. -------------------------------------------------!
   xlonrad = pio180 * xlon
   xlatrad = pio180 * xlat
   plonrad = pio180 * polelon
   platrad = pio180 * polelat
   !---------------------------------------------------------------------------------------!


   !----- Horizontal transform. -----------------------------------------------------------!
   b = 1.0 + sin(xlatrad)*sin(platrad) + cos(platrad)*cos(platrad)*cos(xlonrad- plonrad)
   !---------------------------------------------------------------------------------------!

   f = 2.00 * erad /b


   y = f * (cos(platrad)*sin(xlatrad) - sin(platrad)*cos(xlatrad)*cos(xlonrad-plonrad))

   x = f * (cos(xlatrad)*sin(xlonrad - plonrad))

   return
end subroutine ge_to_xy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine geo_grid(nx,ny,rlat,rlon,dep_glon1,dep_glon2,dep_glat1,dep_glat2,rlatmin        &
                   ,rlatmax,rlonmin,rlonmax,nxg,nyg,proj)
   use misc_coms, only : glong & ! intent(inout)
                       , glatg ! ! intent(inout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)  :: nx
   integer                           , intent(in)  :: ny
   real(kind=4)    , dimension(nx,ny), intent(in)  :: rlat
   real(kind=4)    , dimension(nx,ny), intent(in)  :: rlon
   real(kind=4)                      , intent(out) :: dep_glon1
   real(kind=4)                      , intent(out) :: dep_glon2
   real(kind=4)                      , intent(out) :: dep_glat1
   real(kind=4)                      , intent(out) :: dep_glat2
   real(kind=4)                      , intent(out) :: rlatmin
   real(kind=4)                      , intent(out) :: rlatmax
   real(kind=4)                      , intent(out) :: rlonmin
   real(kind=4)                      , intent(out) :: rlonmax
   integer                           , intent(out) :: nxg
   integer                           , intent(out) :: nyg
   character(len=*)                  , intent(in)  :: proj
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: i
   integer                                         :: j
   integer                                         :: n
   real(kind=4)                                    :: x
   real(kind=4)                                    :: xx
   real(kind=4)                                    :: rlon1
   real(kind=4)                                    :: rlat1
   !---------------------------------------------------------------------------------------!


   dep_glon1 = rlon( 1,1)
   dep_glon2 = rlon(nx,1)
   do n=1,ny
      if(rlon(1,n)  > dep_glon1) dep_glon1 = rlon( 1,n)
      if(rlon(nx,n) < dep_glon2) dep_glon2 = rlon(nx,n)
   end do
   dep_glon2 = ( dep_glon2 - dep_glon1 ) / ( nx - 1 )


   dep_glat1 = rlat( 1,1)
   dep_glat2 = rlat(1,ny)
   do n=1,nx
      if(rlat(n, 1) > dep_glat1) dep_glat1 = rlat(n, 1)
      if(rlat(n,ny) < dep_glat2) dep_glat2 = rlat(n,ny)
   end do
   dep_glat2 = ( dep_glat2 - dep_glat1 ) / ( ny - 1 )

   !10/08/98
   x  = 0
   xx = 0
   do n=1,ny
      x  = x+rlon(1,n)
      xx = xx+ (rlon(nx,n)-rlon(1,n)) / (nx-1)
   end do
   dep_glon1 = x  / ny
   dep_glon2 = xx / ny


   x  = 0
   xx = 0
   do n=1,nx
      x  = x  + rlat(n,1)
      xx = xx + ( rlat(n,ny) - rlat(n,1) ) / ( ny - 1 )
   end do
   dep_glat1 = x  / nx
   dep_glat2 = xx / nx



   !------ Find domain range. -------------------------------------------------------------!
   rlatmin = minval(rlat)
   rlatmax = maxval(rlat)
   rlonmin = minval(rlon)
   rlonmax = maxval(rlon)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Select whether to re-project the grid or not.                                     !
   !---------------------------------------------------------------------------------------!
   select case (trim(proj))
   case ('no')
      !----- Use grid as is, no re-projection. --------------------------------------------!
      nxg = nx
      nyg = ny
      !------------------------------------------------------------------------------------!

   case ('yes')

      !----- Re-project the grid to a regular lon-lat grid. -------------------------------!
      nxg       = int( ( rlonmax - rlonmin ) / dep_glon2 + 0.5 ) + 4
      nyg       = int( ( rlatmax - rlatmin ) / dep_glat2 + 0.5 ) + 4
      rlon1     = rlonmin - ( nxg - nx - 1 ) * dep_glon2
      rlat1     = rlatmin - ( nyg - ny - 1 ) * dep_glat2
      rlon1     = rlonmin - dep_glon2
      rlat1     = rlatmin - dep_glat2
      dep_glat1 = rlat1
      dep_glon1 = rlon1
      !------------------------------------------------------------------------------------!

   case default
      !----- "Maybe" is not an option... --------------------------------------------------!
      call abort_run   ('Invalid value for iproj: '//trim(proj)//'...'                     &
                       ,'geo_grid','rpost_utils.f90')
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!




   !----- Save gridded domain. ------------------------------------------------------------!
   do i=1,nxg
      glong(i) = dep_glon1 + float( i - 1 ) * dep_glon2
   end do

   do j=1,nyg
      glatg(j) = dep_glat1 + float( j - 1 ) * dep_glat2
   end do
   !---------------------------------------------------------------------------------------!

   return
end Subroutine geo_grid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Show the variable range.                                                            !
!------------------------------------------------------------------------------------------!
subroutine RAMS_show_range(nx,ny,nz,ne,vnow,cvar,loc)
   use rout_coms, only : undefflg
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                 , intent(in) :: nx
   integer                                 , intent(in) :: ny
   integer                                 , intent(in) :: nz
   integer                                 , intent(in) :: ne
   real(kind=4)    , dimension(nx,ny,nz,ne), intent(in) :: vnow
   character(len=*)                        , intent(in) :: cvar
   character(len=*)                        , intent(in) :: loc
   !----- Local variables. ----------------------------------------------------------------!
   integer                                              :: x
   integer                                              :: y
   integer                                              :: z
   integer                                              :: e
   integer                                              :: cnt
   real(kind=4)                                         :: vlow
   real(kind=4)                                         :: vhigh
   real(kind=4)                                         :: vmean
   integer                                              :: nundef
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the lowest, middle and highest values.                                 !
   !---------------------------------------------------------------------------------------!
   vlow   =  huge(1.)
   vhigh  = -huge(1.)
   vmean  = 0.
   nundef = 0
   cnt    = 0
   do e=1,ne
      do z=1,nz
         do y=1,ny
            do x=1,nx
               if (vnow(x,y,z,e) == undefflg) then
                  nundef = nundef + 1
               else
                  cnt   = cnt + 1
                  if (vnow(x,y,z,e) < vlow ) vlow  = vnow(x,y,z,e)
                  if (vnow(x,y,z,e) > vhigh) vhigh = vnow(x,y,z,e)
                  vmean = vmean + (vnow(x,y,z,e)-vmean) / cnt
               end if
            end do
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Show range in the standard output.                                                !
   !---------------------------------------------------------------------------------------!
   write(unit=*,fmt='(4(a,1x),3(a,1x,es12.5,1x),2(a,1x,i6,1x))')                           &
      ' Location: ',trim(loc),':   Variable ',trim(cvar),'.   Minimum =',vlow              &
                             ,';   Mean =',vmean,';   Maximum = ',vhigh                    &
                             ,';   # Fine = ',cnt,';   # Undef =',nundef
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_show_range
!==========================================================================================!
!==========================================================================================!
