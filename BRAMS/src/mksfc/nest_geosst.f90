!==========================================================================================!
!======================================= Change Log =======================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine initialises the topography map.                                      !
!------------------------------------------------------------------------------------------!
subroutine toptnest(ngra,ngrb)

   use mem_mksfc
   use mem_grid
   use io_params

   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer, intent(in) :: ngra
   integer, intent(in) :: ngrb
   !------ Local variables. ---------------------------------------------------------------!
   integer             :: ifm
   integer             :: icm
   integer             :: ipat
   integer             :: i
   integer             :: j
   integer             :: k
   integer             :: indfm
   integer             :: ivtime
   integer             :: nc1
   integer             :: mynum
   !---------------------------------------------------------------------------------------!

   do ifm = ngra,ngrb
      icm = nxtnest(ifm)

      !----- Initialize TOPOGRAPHY in toptinit. -------------------------------------------!
      call toptinit(nnxp(ifm),nnyp(ifm),ifm,sfcfile_p(ifm)%topt,sfcfile_p(ifm)%topzo)
      !------------------------------------------------------------------------------------!

      if (icm >= 1 .and. itoptflg(ifm) == 0) then

         !----- Interpolate TOPO from coarser grid. ---------------------------------------!
         call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1,scr1,sfcfile_p(icm)%topt)
         call eintp(scr1,scr2,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
         call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1,scr2,sfcfile_p(ifm)%topt)
         !---------------------------------------------------------------------------------!


         !----- Interpolate TOPO ZO from coarser grid. ------------------------------------!
         call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1,scr1,sfcfile_p(icm)%topzo)
         call eintp(scr1,scr2,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
         call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1,scr2,sfcfile_p(ifm)%topzo)
         !---------------------------------------------------------------------------------!

      elseif (itoptflg(ifm) == 1) then

         !----- Interpolate TOPO from standard dataset. -----------------------------------!
         call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topt,itoptfn(ifm),itoptfn(ifm)     &
                    ,vt2da,vt2db,ifm,'TOP')
         !---------------------------------------------------------------------------------!


         !----- Interpolate TOPO ZO from standard dataset. --------------------------------!
         call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topzo,itoptfn(ifm),itoptfn(ifm)    &
                    ,vt2da,vt2db,ifm,'ZOT')
         !---------------------------------------------------------------------------------!

      elseif (itoptflg(ifm) == 3) then

         !----- Interpolate TOPO from dted dataset. ---------------------------------------!
         call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topt,itoptfn(ifm),itoptfn(ifm)     &
                    ,vt2da,vt2db,ifm,'TOD')
         !---------------------------------------------------------------------------------!

         !----- Interpolate TOPO ZO from dted dataset. ------------------------------------!
         call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topzo,itoptfn(ifm),itoptfn(ifm)    &
                    ,vt2da,vt2db,ifm,'ZOD')
         !---------------------------------------------------------------------------------!

      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    If desired, replace current values of TOPOGRAPHY by some made up by the user/   !
      ! developer in ruser.f90 files.  The default is to not change anything.              !
      !------------------------------------------------------------------------------------!
      call toptinit_user(nnxp(ifm),nnyp(ifm),ifm,sfcfile_p(ifm)%topt ,sfcfile_p(ifm)%topzo)
      !------------------------------------------------------------------------------------!

   end do
   !---------------------------------------------------------------------------------------!




   if (ngra /= ngrb) then

      !------------------------------------------------------------------------------------!
      !      In case topography data have been independently reassigned on any grid,       !
      ! average fine mesh topography sequentially to the coarser grids.                    !
      !------------------------------------------------------------------------------------!
      do ifm = ngrb,ngra,-1
         if (nxtnest(ifm) > ngridsh .and. ifm >= 2) then
            icm = nxtnest(ifm)
            call fdback(sfcfile_p(icm)%topt,sfcfile_p(ifm)%topt,vt2da,scr2,1,nnxp(icm)     &
                       ,nnyp(icm),1,nnxp(ifm),nnyp(ifm),ifm,'terr',vt2db)
         end if
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      In case terrain heights have been independently reassigned on any grid,       !
      ! interpolate coarse grid terrain heights to a temporary fine mesh array.  Fill the  !
      ! fine mesh boundary terrain heights from the temporary array.                       !
      !------------------------------------------------------------------------------------!
      do ifm = ngra,ngrb
         icm = nxtnest(ifm)
         if (icm >= 1) then
            call fillscr(1,nxpmax,nypmax,1,nnxp(icm),nnyp(icm),1,1,scr1,sfcfile_p(icm)%topt)
            call eintp(scr1,scr2,1,nxpmax,nypmax,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
            call fillvar(1,nxpmax,nypmax,1,nnxp(ifm),nnyp(ifm),1,1,scr2,scr1)

            nc1 = jdim * (nstraty(ifm) + 1)
            call ae2(nnxp(ifm),nnyp(ifm),2+nstratx(ifm),nnxp(ifm)-1-nstratx(ifm),1+nc1     &
                    ,nnyp(ifm)-nc1,scr1,sfcfile_p(ifm)%topt)
            call ae1(nnxp(ifm)*nnyp(ifm),sfcfile_p(ifm)%topt,scr1)
         end if
      end do
   end if

   return
end subroutine toptnest
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine will set up the patch areas, land use classes and soil textural    !
! classes for this run.                                                                    !
!------------------------------------------------------------------------------------------!
subroutine geonest_file(ifm)
   use mem_mksfc
   use mem_grid
   use io_params
   use mem_leaf

   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer, intent(in) :: ifm 
   !------ Local variables. ---------------------------------------------------------------!
   integer             :: icm
   integer             :: ipat
   integer             :: i
   integer             :: j
   integer             :: k
   integer             :: ic
   integer             :: jc
   integer             :: indfm
   integer             :: ivtime
   integer             :: nc1
   integer             :: mynum
   !---------------------------------------------------------------------------------------!


   !----- Alias for coarser mesh. ---------------------------------------------------------!
   icm = nxtnest(ifm)



   !---------------------------------------------------------------------------------------!
   !      Assign initial, horizontally homogeneous values for patch area, land use class,  !
   ! and soil textural class in subroutine sfcinit.                                        !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)')         '-----------------------------------------------------'
   write (unit=*,fmt='(a,1x,i5,a)') ' - Initialising patch properties on grid',ifm,'...'
   write (unit=*,fmt='(a)')         '-----------------------------------------------------'
   write (unit=*,fmt='(a)')         ' '
   call sfcinit_file(nnxp(ifm),nnyp(ifm),nzg,npatch,ifm,sfcfile_p(ifm)%patch_area          &
                    ,sfcfile_p(ifm)%leaf_class,sfcfile_p(ifm)%soil_text)
   !---------------------------------------------------------------------------------------!


   write (unit=*,fmt='(a)')         '-----------------------------------------------------'
   write (unit=*,fmt='(a,1x,i5,a)') ' - Starting land use properties on grid',ifm,'...'
   write (unit=*,fmt='(a)')         '-----------------------------------------------------'
   write (unit=*,fmt='(a)')         ' '

   if (icm >= 1 .and. ivegtflg(ifm) == 0) then

      !----- Assign patch areas and patch classes from coarser grid. ----------------------!
      do ipat = 1,npatch
         do j = 1,nnyp(ifm)
            do i = 1,nnxp(ifm)
               ic = ipm(i,ifm)
               jc = ipm(j,ifm)

               sfcfile_p(ifm)%patch_area(i,j,ipat) = sfcfile_p(icm)%patch_area(ic,jc,ipat)
               sfcfile_p(ifm)%leaf_class(i,j,ipat) = sfcfile_p(icm)%leaf_class(ic,jc,ipat)
            end do
         end do
      end do
      !------------------------------------------------------------------------------------!


   elseif (ivegtflg(ifm) == 1) then

      !----- Assign PATCH AREAS and PATCH CLASSES from standard dataset. ------------------!
      call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat,ivegtflg(ifm),ivegtfn(ifm)  &
                       ,isoilflg(ifm),isoilfn(ifm),ndviflg(ifm),ndvifn(ifm)                &
                       ,vndvifil(1,ifm),'veg',platn(ifm),plonn(ifm)                        &
                       ,sfcfile_p(ifm)%soil_text,sfcfile_p(ifm)%patch_area                 &
                       ,sfcfile_p(ifm)%leaf_class,sfcfile_p(ifm)%veg_ndvif)
      !------------------------------------------------------------------------------------!
   end if


   write (unit=*,fmt='(a)')         '-----------------------------------------------------'
   write (unit=*,fmt='(a,1x,i5,a)') ' - Starting soil texture properties on grid',ifm,'...'
   write (unit=*,fmt='(a)')         '-----------------------------------------------------'
   write (unit=*,fmt='(a)')         ' '

   if (icm >= 1 .and. isoilflg(ifm) == 0) then

      !----- Assign soil texture class from coarser grid. ---------------------------------!
      do ipat = 2,npatch
         do k = 1,nzg
            do j = 1,nnyp(ifm)
               do i = 1,nnxp(ifm)
                  ic = ipm(i,ifm)
                  jc = ipm(j,ifm)
                  sfcfile_p(ifm)%soil_text(k,i,j,ipat) =                                   &
                                                     sfcfile_p(icm)%soil_text(k,ic,ic,ipat)
               end do
            end do
         end do
      end do
      !------------------------------------------------------------------------------------!

   elseif (isoilflg(ifm) == 1) then

      !----- Assign soil texture class from standard dataset. -----------------------------!
      call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat,ivegtflg(ifm),ivegtfn(ifm)  &
                       ,isoilflg(ifm),isoilfn(ifm),ndviflg(ifm),ndvifn(ifm)                &
                       ,vndvifil(1,ifm),'soil',platn(ifm),plonn(ifm)                       &
                       ,sfcfile_p(ifm)%soil_text,sfcfile_p(ifm)%patch_area                 &
                       ,sfcfile_p(ifm)%leaf_class,sfcfile_p(ifm)%veg_ndvif)
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      If desired, override current values of patch area, land use class, soil texture  !
   !  class, and NDVI in ruser.f90 subroutines.                                            !
   !---------------------------------------------------------------------------------------!
   call sfcinit_file_user(nnxp(ifm),nnyp(ifm),nzg,npatch,ifm,sfcfile_p(ifm)%patch_area     &
                         ,sfcfile_p(ifm)%leaf_class,sfcfile_p(ifm)%soil_text)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      As a final initialization step, eliminate any land patch area that is less than  !
   ! 1% of the total grid cell area.  Set its area to zero, and compensate by enlarging    !
   ! areas of remaining patches.                                                           !
   !---------------------------------------------------------------------------------------!
   call patch_minsize(nnxp(ifm),nnyp(ifm),npatch,sfcfile_p(ifm)%patch_area)
   !---------------------------------------------------------------------------------------!

   return
end subroutine geonest_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine either initialises or interpolates various LEAF variables (plus a   !
! couple of radiation variables) for which standard RAMS datasets never exist.             !
!------------------------------------------------------------------------------------------!
subroutine geonest_nofile(ngra,ngrb)

   use mem_leaf
   use mem_basic
   use mem_scratch
   use mem_grid
   use io_params
   use mem_radiate      , only : radiate_g  & ! intent(inout)
                               , iswrtyp    & ! intent(in)
                               , ilwrtyp    ! ! intent(in)
   use mem_soil_moisture, only : soil_moist !  ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ngra
   integer, intent(in) :: ngrb
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: isiz
   integer             :: ifm
   integer             :: icm
   integer             :: ipat
   integer             :: i
   integer             :: j
   integer             :: k
   integer             :: indfm
   integer             :: ivtime
   integer             :: nc1
   integer             :: mynum
   integer             :: ic
   integer             :: jc
   integer             :: can_shv
   integer             :: veg_fliq
   !---------------------------------------------------------------------------------------!


   isiz = maxnxp * maxnyp

   do ifm = ngra,ngrb
      icm = nxtnest(ifm)

      !------------------------------------------------------------------------------------!
      !     Copy CO2 to a scratch array.  This way we can use the same sub-routine even    !
      ! when CO2 is not allocated.                                                         !
      !------------------------------------------------------------------------------------!
      if (co2_on) then
         call atob(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),basic_g(ifm)%co2p,scratch%vt3do)
      else
         call ae0(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),scratch%vt3do,co2con(1))
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Find the cosine of the zenith angle, only if radiation will be called.  Then   !
      ! copy the cosine of the zenith angle to a scratch array; otherwise, fill it with a  !
      ! dummy value of one in case radiation is off.  This way we can initialise LEAF      !
      ! using the same sub-routine in both cases.                                          !
      !------------------------------------------------------------------------------------!
      if (iswrtyp + ilwrtyp > 0) then
         call zen(nnxp(ifm),nnyp(ifm),1,nnxp(ifm),1,nnyp(ifm),iswrtyp,ilwrtyp               &
                 ,grid_g(ifm)%glon,grid_g(ifm)%glat,radiate_g(ifm)%cosz)
         call atob(nnxp(ifm)*nnyp(ifm),radiate_g(ifm)%cosz,scratch%vt2dq)
      else
         call azero(nnxp(ifm)*nnyp(ifm),scratch%vt2dq)
      end if
      !------------------------------------------------------------------------------------!


      !----- Then, fill NOFILE LEAF variables with default values in sfcinit_nofile. ------!
      call sfcinit_nofile( nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,nzs,npatch,ifm                &
                         , basic_g(ifm)%theta             , basic_g(ifm)%pi0               &
                         , basic_g(ifm)%pp                , basic_g(ifm)%rv                &
                         , scratch%vt3do                  , leaf_g(ifm)%seatp              &
                         , leaf_g(ifm)%seatf              , leaf_g(ifm)%soil_water         &
                         , leaf_g(ifm)%soil_energy        , leaf_g(ifm)%soil_text          &
                         , leaf_g(ifm)%sfcwater_mass      , leaf_g(ifm)%sfcwater_energy    &
                         , leaf_g(ifm)%sfcwater_depth     , leaf_g(ifm)%ustar              &
                         , leaf_g(ifm)%tstar              , leaf_g(ifm)%rstar              &
                         , leaf_g(ifm)%cstar              , leaf_g(ifm)%zeta               &
                         , leaf_g(ifm)%ribulk             , leaf_g(ifm)%veg_fracarea       &
                         , leaf_g(ifm)%veg_agb            , leaf_g(ifm)%veg_lai            &
                         , leaf_g(ifm)%veg_tai            , leaf_g(ifm)%veg_rough          &
                         , leaf_g(ifm)%veg_height         , leaf_g(ifm)%veg_displace       &
                         , leaf_g(ifm)%veg_albedo         , leaf_g(ifm)%patch_area         &
                         , leaf_g(ifm)%patch_rough        , leaf_g(ifm)%patch_wetind       &
                         , leaf_g(ifm)%leaf_class         , leaf_g(ifm)%soil_rough         &
                         , leaf_g(ifm)%sfcwater_nlev      , leaf_g(ifm)%stom_condct        &
                         , leaf_g(ifm)%ground_rsat        , leaf_g(ifm)%ground_rvap        &
                         , leaf_g(ifm)%ground_temp        , leaf_g(ifm)%ground_fliq        &
                         , leaf_g(ifm)%veg_water          , leaf_g(ifm)%veg_hcap           &
                         , leaf_g(ifm)%veg_energy         , leaf_g(ifm)%can_prss           &
                         , leaf_g(ifm)%can_theiv          , leaf_g(ifm)%can_theta          &
                         , leaf_g(ifm)%can_rvap           , leaf_g(ifm)%can_co2            &
                         , leaf_g(ifm)%sensible_gc        , leaf_g(ifm)%sensible_vc        &
                         , leaf_g(ifm)%evap_gc            , leaf_g(ifm)%evap_vc            &
                         , leaf_g(ifm)%transp             , leaf_g(ifm)%gpp                &
                         , leaf_g(ifm)%plresp             , leaf_g(ifm)%resphet            &
                         , leaf_g(ifm)%veg_ndvip          , leaf_g(ifm)%veg_ndvic          &
                         , leaf_g(ifm)%veg_ndvif          , leaf_g(ifm)%snow_mass          &
                         , leaf_g(ifm)%snow_depth         , scratch%vt2dq                  &
                         , scratch%vt2dr                  , scratch%vt2ds                  &
                         , scratch%vt2da                  , scratch%vt2db                  &
                         , scratch%vt2dc                  , scratch%vt2dd                  &
                         , scratch%vt2de                  , grid_g(ifm)%glat               &
                         , grid_g(ifm)%glon               , grid_g(ifm)%topzo              &
                         , grid_g(ifm)%flpw               , grid_g(ifm)%rtgt               )
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     In case radiation is on, copy back the albedo and upwelling longwave           !
      ! radiation.                                                                         !
      !------------------------------------------------------------------------------------!
      if (iswrtyp + ilwrtyp > 0) then
         call atob(nnxp(ifm)*nnyp(ifm),scratch%vt2dr,radiate_g(ifm)%rlongup)
         call atob(nnxp(ifm)*nnyp(ifm),scratch%vt2ds,radiate_g(ifm)%albedt )
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find out whether we should interpolate from a coarser grid, or simply copy.     !
      !------------------------------------------------------------------------------------!
      if (icm > 0) then
         select case (nofilflg(ifm))
         case (0)
            !----- Assign values from coarse grid cells and patches. ----------------------!
            call coarse2fine_driver(icm,ifm)

         case (1)

            !------------------------------------------------------------------------------!
            !      Interpolate from coarse grid.  We can interpolate water patch directly. !
            ! For land patches, do this by first averaging all coarse grid land patches,   !
            ! interpolate, then assign back to all fine grid land patches.                 !
            !------------------------------------------------------------------------------!
            call patch_interp_driver(icm,ifm)
         end select
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide whether to run an heterogeneous soil moisture initialisation.           !
      !------------------------------------------------------------------------------------!
      select case(trim(soil_moist))
      case ('i','I','a','A')
         call soil_moisture_init(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,npatch,ifm              &
                                ,leaf_g(ifm)%can_theta       ,leaf_g(ifm)%can_prss         &
                                ,grid_g(ifm)%glat            ,grid_g(ifm)%glon             &
                                ,leaf_g(ifm)%soil_water      ,leaf_g(ifm)%soil_energy      &
                                ,leaf_g(ifm)%soil_text       )
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     In case radiation is on, copy back the albedo and upwelling longwave           !
      ! radiation.                                                                         !
      !------------------------------------------------------------------------------------!
      ! if (iswrtyp + ilwrtyp > 0) then
      !    call atob(nnxp(ifm)*nnyp(ifm),radiate_g(ifm)%rlongup,scratch%vt2dr)
      !    call atob(nnxp(ifm)*nnyp(ifm),radiate_g(ifm)%albedt ,scratch%vt2ds)
      ! else
      !    call azero(nnxp(ifm)*nnyp(ifm),scratch%vt2dr)
      !    call azero(nnxp(ifm)*nnyp(ifm),scratch%vt2ds)
      ! end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Override any of the above variable assignments by user-specified changes to   !
      ! subroutine sfcinit_nofile_user.                                                    !
      !------------------------------------------------------------------------------------!
      ! call sfcinit_nofile_user(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,nzs,npatch,ifm          &
      !            , basic_g(ifm)%theta                , basic_g(ifm)%pi0                  &
      !            , basic_g(ifm)%pp                   , basic_g(ifm)%rv                   &
      !            , scratch%vt3do                     , leaf_g(ifm)%soil_water            &
      !            , leaf_g(ifm)%soil_energy           , leaf_g(ifm)%soil_text             &
      !            , leaf_g(ifm)%sfcwater_mass         , leaf_g(ifm)%sfcwater_energy       &
      !            , leaf_g(ifm)%sfcwater_depth        , leaf_g(ifm)%ustar                 &
      !            , leaf_g(ifm)%tstar                 , leaf_g(ifm)%rstar                 &
      !            , leaf_g(ifm)%cstar                 , leaf_g(ifm)%zeta                  &
      !            , leaf_g(ifm)%ribulk                , leaf_g(ifm)%veg_fracarea          &
      !            , leaf_g(ifm)%veg_agb               , leaf_g(ifm)%veg_lai               &
      !            , leaf_g(ifm)%veg_tai               , leaf_g(ifm)%veg_rough             &
      !            , leaf_g(ifm)%veg_height            , leaf_g(ifm)%veg_displace          &
      !            , leaf_g(ifm)%veg_albedo            , leaf_g(ifm)%patch_area            &
      !            , leaf_g(ifm)%patch_rough           , leaf_g(ifm)%patch_wetind          &
      !            , leaf_g(ifm)%leaf_class            , leaf_g(ifm)%soil_rough            &
      !            , leaf_g(ifm)%sfcwater_nlev         , leaf_g(ifm)%stom_condct           &
      !            , leaf_g(ifm)%ground_rsat           , leaf_g(ifm)%ground_rvap           &
      !            , leaf_g(ifm)%ground_temp           , leaf_g(ifm)%ground_fliq           &
      !            , leaf_g(ifm)%veg_water             , leaf_g(ifm)%veg_hcap              &
      !            , leaf_g(ifm)%veg_energy            , leaf_g(ifm)%can_prss              &
      !            , leaf_g(ifm)%can_theiv             , leaf_g(ifm)%can_theta             &
      !            , leaf_g(ifm)%can_rvap              , leaf_g(ifm)%can_co2               &
      !            , leaf_g(ifm)%sensible_gc           , leaf_g(ifm)%sensible_vc           &
      !            , leaf_g(ifm)%evap_gc               , leaf_g(ifm)%evap_vc               &
      !            , leaf_g(ifm)%transp                , leaf_g(ifm)%gpp                   &
      !            , leaf_g(ifm)%plresp                , leaf_g(ifm)%resphet               &
      !            , leaf_g(ifm)%veg_ndvip             , leaf_g(ifm)%veg_ndvic             &
      !            , leaf_g(ifm)%veg_ndvif             , leaf_g(ifm)%snow_mass             &
      !            , leaf_g(ifm)%snow_depth            , scratch%vt2dq                     &
      !            , scratch%vt2dr                     , scratch%vt2ds                     &
      !            , scratch%vt2da                     , scratch%vt2db                     &
      !            , scratch%vt2dc                     , scratch%vt2dd                     &
      !            , scratch%vt2de                     , grid_g(ifm)%glat                  &
      !            , grid_g(ifm)%glon                  , grid_g(ifm)%topzo                 &
      !            , grid_g(ifm)%flpw                  , grid_g(ifm)%rtgt                  )
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     In case radiation is on, copy back the albedo and upwelling longwave           !
      ! radiation.                                                                         !
      !------------------------------------------------------------------------------------!
      ! if (iswrtyp + ilwrtyp > 0) then
      !    call atob(nnxp(ifm)*nnyp(ifm),scratch%vt2dr,radiate_g(ifm)%rlongup)
      !    call atob(nnxp(ifm)*nnyp(ifm),scratch%vt2ds,radiate_g(ifm)%albedt )
      ! end if
      !------------------------------------------------------------------------------------!
   end do

   return
end subroutine geonest_nofile
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine patch_interp(icm,ifm,nc1,nc2,nc3,nc4,nf1,nf2,nf3,nf4,ac,af,pareac,pareaf,avgc   &
                       ,avgf,slabc,slabf)
   use mem_scratch

   implicit none

   integer :: icm,ifm,nc1,nc2,nc3,nc4,nf1,nf2,nf3,nf4
   real :: ac(nc1,nc2,nc3,nc4), af(nf1,nf2,nf3,nf4)
   real :: pareac(nc2,nc3,nc4), pareaf(nf2,nf3,nf4)
   real :: avgc(nc1,nc2,nc3),   avgf(nf1,nf2,nf3)
   real :: slabc(nc2,nc3), slabf(nf2,nf3)

   integer :: k,i,j

   ! Average coarse grid field over all land patches

   call patch_land_average(nc1,nc2,nc3,nc4,pareac,ac,avgc)

   ! Interpolate patch-averaged to fine grid

   do k=1,nc1    ! nc1 and nf1 are the same

      do j=1,nc3
         do i=1,nc2
            slabc(i,j)=avgc(k,i,j)
         enddo
      enddo

      call fmint2d(icm,ifm,'t',slabc,slabf)

      do j=1,nf3
         do i=1,nf2
            avgf(k,i,j)=slabf(i,j)
         enddo
      enddo

   enddo

   ! Fill fine grid field back into all land patches

   call patch_land_unaverage(nf1,nf2,nf3,nf4,avgf,af)

   return
end subroutine patch_interp


!     *****************************************************************

subroutine fmint5(var1,var2,dn0xc,dn0xf,vt2da,ifm,icm,vpnt,idwt)

  use mem_scratch
  use mem_grid

  implicit none
  integer :: ifm,icm,idwt

  real, dimension(*) :: var1,var2,vt2da,dn0xc,dn0xf
  character(len=*) :: vpnt

  if (icm .eq. 0) return

  call fillscr(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
       ,1,nnzp(icm),scratch%scr1,var1)


  if (idwt .eq. 1) then
     call dnswt2(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,scratch%scr1,dn0xc,vpnt,1)
  endif

  call eintp(scratch%scr1,scratch%scr2,maxnzp,maxnxp,maxnyp  &
       ,nnzp(ifm),nnxp(ifm),nnyp(ifm),ifm,3,vpnt,0,0)

  call fillvar(maxnzp,maxnxp,maxnyp,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
       ,1,nnzp(ifm),scratch%scr2,var2)

  if (idwt .eq. 1) call dnswt2(nnzp(ifm),nnxp(ifm),nnyp(ifm),nnzp(ifm)  &
       ,nnxp(ifm),nnyp(ifm),var2,dn0xf,vpnt,2)

  call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),var2,vt2da  &
       ,grid_g(ifm)%topt,ifm,vpnt)

  return
end subroutine fmint5
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine patch_minsize(nx,ny,np,patch_area)
   use leaf_coms, only : min_patch_area ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nx
   integer                  , intent(in)    :: ny
   integer                  , intent(in)    :: np
   real, dimension(nx,ny,np), intent(inout) :: patch_area
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: x,y,p
   real                                     :: parea_tot
   !---------------------------------------------------------------------------------------!

   yloop: do y=1,ny
      xloop: do x=1,nx
         !----- First ensure that the water patch has at least 0.1% of the area. ----------!
         if (patch_area(x,y,1) < min_patch_area) patch_area(x,y,1) = min_patch_area

         !---------------------------------------------------------------------------------!
         !     Then we loop over the other patches, eliminating those that are tiny or     !
         ! make no sense.                                                                  !
         !---------------------------------------------------------------------------------!
         parea_tot = patch_area(x,y,1)
         elimploop: do p=2,np
            if (patch_area(x,y,p) > 0. .and. patch_area(x,y,p) < min_patch_area) then
               patch_area(x,y,p) = 0.
            elseif (patch_area(x,y,p) < 0. .or. patch_area(x,y,p) > 1.001) then
               patch_area(x,y,p) = 0.
            end if
            parea_tot = parea_tot + patch_area(x,y,p)
         end do elimploop

         !----- Loop over the patches and re-scale them. ----------------------------------!
         rsclploop: do p=1,np
            patch_area(x,y,p) = patch_area(x,y,p) / parea_tot
         end do rsclploop
         !---------------------------------------------------------------------------------!

      end do xloop
   end do yloop

   return
end subroutine patch_minsize
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! TEB
!#############################################################################
subroutine fusonest(ngra,ngrb)
!#############################################################################

  use mem_mksfc
  use mem_grid
  use io_params

  implicit none

  integer :: ngra,ngrb 

  integer :: ifm,icm,ipat,i,j,k,indfm,ivtime,nc1,mynum



  do ifm = ngra,ngrb
     icm = nxtnest(ifm)
     ! Initialize FUSO in fusoinit.
     
!     write(*,*)'glon =',grid_g(ifm)%glon (1,1)
!     call fusoinit(nnxp(ifm),nnyp(ifm),ifm  &
!          ,sfcfile_p(ifm)%fuso(1,1),grid_g(ifm)%glon (1,1))

     if (icm .ge. 1 .and. ifusflg(ifm) .eq. 0) then

        ! Interpolate FUSO from coarser grid:
        call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
             ,scr1,sfcfile_p(icm)%fuso)
        call eintp(scr1,scr2,1,maxnxp,maxnyp  &
             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
        call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
             ,scr2,sfcfile_p(ifm)%fuso)

     elseif (ifusflg(ifm) .eq. 1) then

        ! Interpolate FUSO from standard dataset:
        call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%fuso  &
             ,ifusfn(ifm),ifusfn(ifm),vt2da,vt2db,ifm,'FUS')

     endif


  enddo

  if (ngra .eq. ngrb) return

  return
end subroutine fusonest
!==========================================================================================!
!==========================================================================================!
