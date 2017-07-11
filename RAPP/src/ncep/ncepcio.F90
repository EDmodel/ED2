!==========================================================================================!
!==========================================================================================!
! MEVI. ncepcio.F90  This is the equivalent of rcio for NCEP, so it reads the global       !
!                    attributes and the dimensions. It will read and store only the        !
!                    variables that are needed for conversion into the output formats.     !
!------------------------------------------------------------------------------------------!
subroutine commio_ncep(ngrid)
   use mod_model      , only : polelon         & ! intent(out)
                             , polelat         & ! intent(out)
                             , l2ndlat         & ! intent(out)
                             , centlon         & ! intent(inout)
                             , centlat         & ! intent(inout)
                             , nnxp            & ! intent(inout)
                             , nnyp            & ! intent(inout)
                             , nnzp            & ! intent(inout)
                             , nntp            & ! intent(inout)
                             , nclouds         & ! intent(inout)
                             , npatch          & ! intent(inout)
                             , nzg             & ! intent(inout)
                             , nzs             & ! intent(inout)
                             , nwave           & ! intent(inout)
                             , ihtran          & ! intent(inout)
                             , deltaxn         & ! intent(inout)
                             , deltayn         & ! intent(inout)
                             , zero_time       & ! intent(inout)
                             , this_time       & ! intent(inout)
                             , xtn             & ! intent(inout)
                             , ytn             ! ! intent(inout)
   use mod_maxdims    , only : maxtimes        & ! intent(in)
                             , maxgrds         ! ! intent(in)
   use mod_ioopts     , only : missflg_int     ! ! intent(in)
#if USE_NCDF
   use netcdf
   use mod_netcdf     , only : idnntp          & ! intent(inout)
                             , idnnxp          & ! intent(inout)
                             , idnnyp          & ! intent(inout)
                             , idnnzp          & ! intent(inout)
                             , idnzg           & ! intent(inout)
                             , idnzs           & ! intent(inout)
                             , idnclouds       & ! intent(inout)
                             , idnpatch        & ! intent(inout)
                             , idnwave         ! ! intent(inout)
   use mod_ncdf_globio, only : ncio_glo        & ! function
                             , ncio_glo_sca    & ! function
                             , ncio_dim        & ! function
                             , ncio_1dnotime   & ! function
                             , ncio_ncep_time  ! ! function
   use rconstants     , only : spcon           ! ! intent(in)
#endif


   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer          , intent(in)       :: ngrid   ! Grid ID regarding this variable.
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: curryear,currdoy,it,s,t,z
   real   , dimension(maxtimes)        :: gmttime
   integer, dimension(maxtimes)        :: years,doys
   integer                             :: timelength
!------ Giving the preprocessor the option to entirely skip this routine ------------------!
#if USE_NCDF
   integer                           :: ierr
   character(len=NF90_MAX_NAME)      :: dumchar
   integer                           :: mapproj
   integer                           :: dummyint
   logical                           :: gotz
   integer, dimension(maxgrds), save :: nnzpmax = 0
   !---------------------------------------------------------------------------------------!


   !----- Getting the dimensions ----------------------------------------------------------!
   ierr   = ncio_dim   ('time'                  , .true.  ,nntp(ngrid)    ,idnntp    )
   ierr   = ncio_dim   ('lon'                   , .true.  ,nnxp(ngrid)    ,idnnxp    )
   ierr   = ncio_dim   ('lat'                   , .true.  ,nnyp(ngrid)    ,idnnyp    )

   !---------------------------------------------------------------------------------------!
   !    Not always will the NCEP file have level. If the file has only sfc variables, they !
   ! won't be present, in which case we will assign 1                                      !
   !---------------------------------------------------------------------------------------!
   ierr   = ncio_dim   ('level'                 , .false. ,nnzp(ngrid)    ,idnnzp    )
   if (ierr == NF90_NOERR) then
      if (nnzpmax(ngrid) < dummyint) then
         gotz           = .true.
         nnzp(ngrid)    = dummyint
         nnzpmax(ngrid) = dummyint 
      else
         gotz           = .false.
      end if
   else 
      gotz = .false.
   end if

   !----- Other dimensions will be set to 1 and their IDs to a non-sense number -----------!
   nzg       = 1
   nzs       = 1
   npatch    = 1
   nclouds   = 1
   nwave     = 1 
   idnzg     = missflg_int
   idnzs     = missflg_int
   idnpatch  = missflg_int
   idnclouds = missflg_int
   idnwave   = missflg_int

   !----- NCEP will always have a lon/lat list (either grid or Gaussian). -----------------!
   ihtran = 0

   !----- Getting longitude, latitude and levels ------------------------------------------!
   ierr   = ncio_1dnotime('lon'           ,.true.,nnxp(ngrid), xtn(1:nnxp(ngrid),ngrid))
   ierr   = ncio_1dnotime('lat'           ,.true.,nnyp(ngrid), ytn(1:nnyp(ngrid),ngrid))
   
   polelon = minval(xtn(1:nnxp(ngrid),ngrid),dim=1)
   polelat = maxval(ytn(1:nnyp(ngrid),ngrid),dim=1)
   l2ndlat = polelat
   
   !----- Finding delta-x and delta-y in metres -------------------------------------------! 
   deltaxn(ngrid) = (xtn(2,ngrid)-xtn(1,ngrid))*spcon
   deltayn(ngrid) = (ytn(2,ngrid)-ytn(1,ngrid))*spcon
   !----- Dumping anything on centlon and centlat -----------------------------------------!
   centlon(ngrid) = (xtn(nnxp(ngrid)/2,ngrid))
   centlat(ngrid) = (ytn(nnyp(ngrid)/2,ngrid))
   !----- Retrieving all times ------------------------------------------------------------!
   ierr   = ncio_ncep_time('time', .true., nntp(ngrid),this_time(1:nntp(ngrid),ngrid))

#else
   call fatal_error ('You can''t run ncep without compiling with netcdf!'                  &
                    ,'commio_ncep','ncepcio.F90')
#endif

   return
end subroutine commio_ncep
!==========================================================================================!
!==========================================================================================!
