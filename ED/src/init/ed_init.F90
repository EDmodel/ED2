module ed_init
  contains

!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign the longitude, latitude, and soil class for all non-empty !
! polygons.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine set_polygon_coordinates()
   use grid_coms     , only : ngrids    ! ! intent(in)
   use ed_work_vars  , only : work_v    ! ! structure
   use ed_node_coms  , only : mynum     ! ! intent(in)
   use ed_state_vars , only : edgrid_g  & ! structure
                            , gdpy      & ! intent(in)
                            , edtype    ! ! intent(in)
   implicit none 
   !----- Local variables -----------------------------------------------------------------!
   integer               :: ifm
   integer               :: ipy
   integer               :: npoly
   type(edtype), pointer :: cgrid
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !  Ifort-11.0.083 complained about using edgrid_g(ifm) directly, replaced by cgrid.     !
   !---------------------------------------------------------------------------------------!
   gridloop: do ifm=1,ngrids
      cgrid => edgrid_g(ifm)

      npoly=gdpy(mynum,ifm)
      !----- Go through every polygon. ----------------------------------------------------!
      polyloop: do ipy=1,npoly
         cgrid%lon (ipy) = work_v(ifm)%glon   (ipy)
         cgrid%lat (ipy) = work_v(ifm)%glat   (ipy)
         cgrid%xatm(ipy) = work_v(ifm)%xid    (ipy)
         cgrid%yatm(ipy) = work_v(ifm)%yid    (ipy)
      end do polyloop
      !------------------------------------------------------------------------------------!
   end do gridloop
   !---------------------------------------------------------------------------------------!


   return
end subroutine set_polygon_coordinates
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine assigns the site areas and the soil types when not running with       !
! ied_init_node 3 or 4.                                                                    !
!------------------------------------------------------------------------------------------!
subroutine set_site_defprops()
   use grid_coms     , only : ngrids               & ! intent(in)
                            , nzg                  ! ! intent(in)
   use ed_work_vars  , only : work_v               ! ! structure
   use ed_state_vars , only : edgrid_g             & ! structure
                            , edtype               & ! structure
                            , polygontype          & ! structure
                            , allocate_polygontype ! ! subroutine
   use soil_coms     , only : soil                 & ! intent(in)
                            , slz                  & ! intent(in)
                            , nslcon
   use mem_polygons  , only : maxsite              ! ! intent(in)
   use ed_misc_coms  , only : min_site_area        ! ! intent(in)
   implicit none 
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)     , pointer :: cgrid
   type(polygontype), pointer :: cpoly
   integer                    :: ifm
   integer                    :: ipy
   integer                    :: isi
   integer                    :: itext
   integer                    :: nsite
   integer                    :: k
   integer                    :: sc
   real                       :: zmin
   real                       :: fa
   real                       :: fb
   real                       :: te
   real                       :: t0
   real                       :: k0
   real                       :: text_area
   real                       :: text_area_i
   !---------------------------------------------------------------------------------------!



   gridloop: do ifm=1,ngrids
      cgrid => edgrid_g(ifm)

      polyloop: do ipy=1,cgrid%npolygons
      
         !----- Initialise load adjacency with dummy value. -------------------------------!
         cgrid%load_adjacency(ipy) = 0

         !----- Alias to current polygon. -------------------------------------------------!
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !    Find the total usable area and the number of sites that will be allocated.   !
         !---------------------------------------------------------------------------------!
         nsite     = min(maxsite, count(work_v(ifm)%soilfrac(:,ipy) > min_site_area))
         text_area = sum(work_v(ifm)%soilfrac(1:nsite,ipy))
         !----- Sanity check. -------------------------------------------------------------!
         if (nsite == 0 .or. text_area <= 0.) then
            write (unit=*,fmt='(a)'          ) '------------------------------------------'
            write (unit=*,fmt='(a,1x,i12)'   ) ' + NSITE     = ',nsite
            write (unit=*,fmt='(a,1x,es12.5)') ' + TEXT_AREA = ',text_area
            write (unit=*,fmt='(a,1x,es12.5)') ' + MIN_VALID = ',min_site_area
            write (unit=*,fmt='(a,1x,es12.5)') ' + MIN_AREA  = '                           &
                                              ,minval(work_v(ifm)%soilfrac(:,ipy))
            write (unit=*,fmt='(a,1x,es12.5)') ' + MAX_AREA  = '                           &
                                              ,maxval(work_v(ifm)%soilfrac(:,ipy))
            write (unit=*,fmt='(a)'          ) '------------------------------------------'
            call fatal_error('Invalid number of sites!','set_site_defprops'                &
                            ,'ed_init.f90')
         else
            text_area_i = 1. / text_area
         end if


         !------ Allocate the number of sites that have enough area. ----------------------!
         call allocate_polygontype(cpoly,nsite)
         call soil_default_fill(cgrid,ifm,ipy)

         !---------------------------------------------------------------------------------!
         !     Populate the sites.                                                         !
         !---------------------------------------------------------------------------------!
         isi = 0
         siteloop: do itext=1,maxsite
            if (work_v(ifm)%soilfrac(itext,ipy) > min_site_area) then
               !----- Update counter. -----------------------------------------------------!
               isi = isi +1

               !---------------------------------------------------------------------------!
               !     Initialise area, using the soil texture area scaled by the total be-  !
               ! ing used.                                                                 !
               !---------------------------------------------------------------------------!
               cpoly%area(isi) = work_v(ifm)%soilfrac(itext,ipy) * text_area_i
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Use the soil type and populate the site-level soil texture.           !
               !---------------------------------------------------------------------------!
               do k=1,nzg
                  cpoly%ntext_soil(k,isi) = nslcon(k) !work_v(ifm)%ntext(itext,ipy)
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !       Set soil moisture decay function, based on second layer's K value.  !
               ! We use the second layer instead of the top in case top is organic/peat.   !
               !---------------------------------------------------------------------------!
               sc = cpoly%ntext_soil(nzg-1,isi)
               cpoly%moist_f(isi) = -log(soil(sc)%slcons / soil(sc)%slcons0) / 2.0
               !---------------------------------------------------------------------------!




               !----- Derive adjustments to f. --------------------------------------------!
               zmin = slz(cpoly%lsl(isi))
               fa   = -1.0/zmin !! should be 1/(depth to bedrock)
               if(cpoly%moist_f(isi)*zmin < 0.0) then
                  fb = cpoly%moist_f(isi)/(1.0-exp(cpoly%moist_f(isi)*zmin))
                  cpoly%moist_f(isi) = max(fa,fb)
               else
                  cpoly%moist_f(isi) = fa
               endif

               cpoly%sitenum(isi)   = 1
               cpoly%elevation(isi) = 0.0
               cpoly%slope(isi)     = 0.0
               cpoly%aspect(isi)    = 0.0
               cpoly%TCI(isi)       = 0.0
            end if
         end do siteloop

         !----- This should not happen in this sub-routine, but in case things change... --!
         if (cgrid%load_adjacency(ipy) /= 0) then
            call calc_flow_routing(cgrid,ipy)
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Not sure what these things do, just copying from hydrology...               !
         !---------------------------------------------------------------------------------!
         !----- Part 1. -------------------------------------------------------------------!
         Te = 0.0 
         do isi = 1,cpoly%nsites
            sc = cpoly%ntext_soil(nzg-1,isi)
            K0 = soil(sc)%slcons0
            T0 = K0 / cpoly%moist_f(isi)
            Te = Te + T0*cpoly%area(isi)
         end do
         cgrid%Te(ipy) = Te
         !----- Part 2. -------------------------------------------------------------------!
         cgrid%wbar(ipy) = 0.0
         do isi = 1,cpoly%nsites
            sc = cpoly%ntext_soil(nzg-1,isi)
            K0 = soil(sc)%slcons0
            T0 = K0 / cpoly%moist_f(isi)
            cpoly%moist_W(isi) = cpoly%TCI(isi) + log(Te) - log(T0)
            cgrid%wbar(ipy)    = cgrid%wbar(ipy) + cpoly%moist_W(isi) * cpoly%area(isi)
         end do
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

   end do gridloop
   !---------------------------------------------------------------------------------------!


   return
end subroutine set_site_defprops
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine fills the lsl, soil colour, and soil texture based on the defaults.  !
! In case isoildepthflg was zero, then the layer_index matrix was filled with ones, so we  !
! do not need to worry about this here.                                                    !
!------------------------------------------------------------------------------------------!
subroutine soil_default_fill(cgrid,ifm,ipy)
   
   use soil_coms     , only : layer_index & ! intent(in)
                            , nslcon 
   use ed_state_vars , only : edtype      & ! structure
                            , polygontype ! ! structure
   use ed_work_vars  , only : work_v      ! ! structure
   use grid_coms     , only : nzg         ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype) , target       :: cgrid
   integer      , intent(in)   :: ifm
   integer      , intent(in)   :: ipy
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   integer                     :: ilat_bin
   integer                     :: ilon_bin
   integer                     :: isi
   integer                     :: k
   !---------------------------------------------------------------------------------------!


   ilat_bin = min(180,int(90.0 - cgrid%lat(ipy)) + 1)
   ilon_bin = int(180.0 + cgrid%lon(ipy)) + 1

   cpoly => cgrid%polygon(ipy)
   do isi=1,cpoly%nsites
      !------------------------------------------------------------------------------------!
      !    Require at least 2 layers.  This requirement was taken in consideration when    !
      ! layer_index was filled at the first initialization, so it is safe to just copy.    !
      !------------------------------------------------------------------------------------!
      cpoly%lsl(isi) =layer_index(ilat_bin,ilon_bin) 
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the commonest soil type and populate the site-level soil texture.          !
      !------------------------------------------------------------------------------------!
      do k=1,nzg
         cpoly%ntext_soil(k,isi) = nslcon(k) !work_v(ifm)%ntext(1,ipy) EJL
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the polygon-level soil colour to populate the site-level.                  !
      !------------------------------------------------------------------------------------!
      cpoly%ncol_soil(isi) = work_v(ifm)%nscol(ipy)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine soil_default_fill
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Several Procedures requiring ASCII reads follow.  If this is a parallel run, then the !
! nodes must queue.  Since we will access sequential format files, each node needs to wait !
! its turn to access it... MPI_File commands won't work with ASCII files, so that's a      !
! bottleneck here. If the run is serial mynum=nnodetot, so I don't need to wait.           !
!------------------------------------------------------------------------------------------!
subroutine load_ecosystem_state()
   use landuse_init_module
   use ed_nbg_init
   use ed_misc_coms      , only : ied_init_mode   & ! intent(in)
                                , ibigleaf        ! ! intent(in)
   use phenology_startup , only : phenology_init  ! ! intent(in)
#if defined(RAMS_MPI)
   use ed_node_coms      , only : mynum           & ! intent(in)
                                , nmachs          & ! intent(in)
                                , nnodetot        & ! intent(in)
                                , mchnum          & ! intent(in)
                                , machs           & ! intent(in)
                                , master_num      & ! intent(in)
                                , sendnum         & ! intent(in)
                                , recvnum         ! ! intent(in)
#else
   use ed_node_coms      , only : mynum           ! ! intent(in)
#endif
   use grid_coms         , only : ngrids,nzl          ! ! intent(in)
   use ed_state_vars     , only : edgrid_g        ! ! structure

   implicit none
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Local variables -----------------------------------------------------------------!
   integer                :: igr
   integer                :: ping 
   !----- Local variables (MPI only). -----------------------------------------------------!
#if defined(RAMS_MPI)
   integer                :: ierr
#endif
   !---------------------------------------------------------------------------------------!

   ping = 741776


   if (mynum == 1) write(unit=*,fmt='(a)') ' + Doing sequential initialization over nodes.'


   !---------------------------------------------------------------------------------------!
   ! STEP 1: Read in Site files and initialize hydrologic adjacencies.                     !
   !---------------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (mynum /= 1) &
      call MPI_Recv(ping,1,MPI_INTEGER,recvnum,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#endif

   select case (ied_init_mode)
   case (3)
      !----- Hydrology run.  Use specific scheme. -----------------------------------------!
      do igr = 1,ngrids
         call read_site_file(edgrid_g(igr),igr)
      end do
   case (4,7)
      continue
   case default
      call set_site_defprops()
   end select

#if defined(RAMS_MPI)
   if (mynum < nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,100,MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   ! STEP 3: Do ASCII type restart initialization of site patch and cohort biophysical     !
   !         states.                                                                       !
   !---------------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (mynum /= 1) &
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#endif

   select case (ied_init_mode)
   case (-8,-1,0)
   
      select case (ibigleaf)
      case (0)
         !----- Initialize everything with near-bare ground -------------------------------!
         if (mynum /= 1) then
            write (unit=*,fmt='(a)') ' + Doing near bare ground initialization...'
         end if
         do igr=1,ngrids
              call near_bare_ground_init(edgrid_g(igr))
         end do

      case (1)
         !----- Initialize everything with near-bare ground -------------------------------!
         if (mynum /= 1) then
            write(unit=*,fmt='(a)') ' + Doing near-bare-ground big-leaf initialization...'
         end if
         do igr=1,ngrids
            call near_bare_ground_big_leaf_init(edgrid_g(igr))
         end do
      end select

   case (1,2,3,6)
      !----- Initialize with ED1-type restart information. --------------------------------!
      write(unit=*,fmt='(a,i3.3)') ' + Initializing from ED restart file. Node: ',mynum
      call read_ed10_ed20_history_file

      select case (ibigleaf)
      case (1)
         do igr=1,ngrids
            call ed_bigleaf_init(edgrid_g(igr))
         end do
      end select

   case (4)   
      write(unit=*,fmt='(a,i3.3)') ' + Initializing from ED2.1 state file. Node: ',mynum
      call read_ed21_history_file
      select case (ibigleaf)
      case (1)
         do igr=1,ngrids
            call ed_bigleaf_init(edgrid_g(igr))
         end do
      end select

   case (5,99)
      write(unit=*,fmt='(a,i3.3)')                                                         &
          ' + Initializing from a collection of ED2.1 state files. Node: ',mynum
      call read_ed21_history_unstruct
      select case (ibigleaf)
      case (1)
         do igr=1,ngrids
            call ed_bigleaf_init(edgrid_g(igr))
         end do
      end select

   case (7)
      write(unit=*,fmt='(a,i3.3)')                                                         &
           ' + Initializing Soils+Veg from nearest neighbor ED2.1 files. Node: ',mynum
      call read_ed21_polyclone
      select case (ibigleaf)
      case (1)
         do igr=1,ngrids
            call ed_bigleaf_init(edgrid_g(igr))
         end do
      end select
      
   end select



#if defined(RAMS_MPI)
   if (mynum < nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,101,MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! STEP 4: Initialize phenology parameters and thermal sums.                             !
   !---------------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (mynum /= 1) &
      call MPI_Recv(ping,1,MPI_INTEGER,recvnum,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#endif

   write(unit=*,fmt='(a,i3.3)') ' + Initializing phenology. Node: ',mynum
   call phenology_init()

#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,102,MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! STEP 5: Initialize anthropogenic disturbance.                                         !
   !---------------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (mynum /= 1) &
     call MPI_Recv(ping,1,MPI_INTEGER,recvnum,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#endif

   write(unit=*,fmt='(a,i3.3)')                                                            &
      ' + Initializing anthropogenic disturbance forcing. Node: ',mynum

   call landuse_init()

#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,103,MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!

   if (mynum == 1) then
      do igr=1,ngrids
         call ed_newgrid(igr)
         call print_soil_info(edgrid_g(igr),igr)
      end do
   end if


   return
end subroutine load_ecosystem_state
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine defines the variables related to the soil layers, and also initial-  !
! ises some RK4 variables that depend on the soil grid.                                    !
!------------------------------------------------------------------------------------------!
subroutine sfcdata_ed()
   use rk4_coms    , only : ipercol           ! ! intent(in)
   use grid_coms   , only : nzg               & ! intent(in)
                          , nzs               & ! intent(in)
                          , nzl               ! ! intent(in)
   use soil_coms   , only : ed_nstyp          & ! intent(in)
                          , slz               & ! intent(in)
                          , dslz              & ! intent(out)
                          , dslzo2            & ! intent(out)
                          , dslzi             & ! intent(out)
                          , dslzidt           & ! intent(out)
                          , slzt              & ! intent(out)
                          , dslzt             & ! intent(out)
                          , dslzti            & ! intent(out)
                          , dslztidt          & ! intent(out)
                          , slz8              & ! intent(in)
                          , dslz8             & ! intent(out)
                          , dslzo28           & ! intent(out)
                          , dslzi8            & ! intent(out)
                          , dslzidt8          & ! intent(out)
                          , slzt8             & ! intent(out)
                          , dslzt8            & ! intent(out)
                          , dslzti8           & ! intent(out)
                          , dslztidt8         & ! intent(out)
                          , fhydraul          & ! intent(out)
                          , slcons1           & ! intent(out)
                          , slcons18          & ! intent(out)
                          , slden             & ! intent(out)
                          , soil              & ! intent(in)
                          , thicknet          & ! intent(out)
                          , thick             & ! intent(out)
                          , olz               & ! intent(in)
                          , dolz              & ! intent(out)
                          , dolzo2            & ! intent(out)
                          , dolzi             & ! intent(out)
                          , dolzidt           & ! intent(out)
                          , olzt              & ! intent(out)
                          , dolzt             & ! intent(out)
                          , dolzti            & ! intent(out)
                          , dolztidt          & ! intent(out)
                          , olz8              & ! intent(in)
                          , dolz8             & ! intent(out)
                          , dolzo28           & ! intent(out)
                          , dolzi8            & ! intent(out)
                          , dolzidt8          & ! intent(out)
                          , olzt8             & ! intent(out)
                          , dolzt8            & ! intent(out)
                          , dolzti8           & ! intent(out)
                          , dolztidt8         & ! intent(out)
                          , ohydraul          & ! intent(out)
                          , olcons1           & ! intent(out)
                          , olcons18          & ! intent(out)
                          , olden             ! ! intent(out)
   use consts_coms , only : wdns              & ! intent(in)
                          , wdnsi8            ! ! intent(in)
   use rk4_coms    , only : rk4min_sfcw_moist & ! intent(in)
                          , rk4min_virt_moist & ! intent(in)
                          , rk4min_sfcw_mass  & ! intent(out)
                          , rk4min_virt_water ! ! intent(out)
   use ed_misc_coms, only : dtlsm             ! ! intent(in)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer                :: k
   integer                :: nnn
   integer                :: kzs
   real                   :: refdepth
   real                   :: thik
   real                   :: stretch
   real                   :: slz0
   real                   :: ezg
   real                   :: olz0
   real                   :: oezg
   !---------------------------------------------------------------------------------------!


   !----- Soil vertical grid spacing arrays (some with timestep info). --------------------!
   slz (nzg+1) = 0.
   slz8(nzg+1) = 0.d0

   do k = 1,nzg
      dslz    (k) = slz(k+1) - slz(k)
      dslzo2  (k) = .5 * dslz(k)
      dslzi   (k) = 1. / dslz(k)
      dslzidt (k) = dslzi(k) * dtlsm
      slzt    (k) = .5 * (slz(k) + slz(k+1))
      slz8    (k) = dble(slz    (k))
      dslz8   (k) = dble(dslz   (k))
      dslzo28 (k) = dble(dslzo2 (k))
      dslzi8  (k) = dble(dslzi  (k))
      dslzidt8(k) = dble(dslzidt(k))
      slzt8   (k) = dble(slzt   (k))
   end do

   !----- Find the exponential increase factor to estimate the bottom boundary condition. -!
   ezg  = log(slz(1)/slz(nzg)) / log(real(nzg))
   slz0 = slz(1) * (real(nzg+1)/real(nzg))**ezg
   !---------------------------------------------------------------------------------------!


   !----- Find the thickness of the bottom boundary condition layer. ----------------------!
   dslz    (0) = slz(1) - slz0
   dslzo2  (0) = .5 * dslz(0)
   dslzi   (0) = 1. / dslz(0)
   dslzidt (0) = dslzi(0) * dtlsm
   dslz8   (0) = dble(dslz   (0))
   dslzo28 (0) = dble(dslzo2 (0))
   dslzi8  (0) = dble(dslzi  (0))
   dslzidt8(0) = dble(dslzidt(0))
   !---------------------------------------------------------------------------------------!


   !----- Find the height at the middle of the bottom boundary condition. -----------------!
   slzt   (0) = .5 * (slz0 + slz(1))
   slzt8  (0) = dble(slzt(0))
   !---------------------------------------------------------------------------------------!

   do k = 1,nzg
      dslzt    (k) = slzt(k) - slzt(k-1)
      dslzti   (k) = 1. / dslzt(k)
      dslztidt (k) = dslzti(k) * dtlsm
      dslzt8   (k) = dble(dslzt   (k))
      dslzti8  (k) = dble(dslzti  (k))
      dslztidt8(k) = dble(dslztidt(k))
   end do


   !----- Soil constants. -----------------------------------------------------------------!
   refdepth = -0.5

   do nnn = 1,ed_nstyp
      if (nnn /= 13) then
         fhydraul(nnn) = log (soil(nnn)%slcons / soil(nnn)%slcons0) / refdepth
      else
         fhydraul(nnn) = 0.
      end if
      do k = 0,nzg
      
         select case (ipercol)
         case (0,1)
            !----- Original form, constant with depth.  -----------------------------------!
            slcons1(k,nnn) = soil(nnn)%slcons
         case (2)
            !------------------------------------------------------------------------------!
            !    TOPMODEL form, similar to CLM.  Here we use the same definition of slcons !
            ! from Cosby et al. (1984) because it has a stronger spread and it accounts    !
            ! for sand and clay contents.                                                  !
            !------------------------------------------------------------------------------!
            slcons1(k,nnn) = soil(nnn)%slcons * exp ( - slzt(k) / refdepth)
         end select

         !------ Find the double precision. -----------------------------------------------!
         slcons18(k,nnn) = dble(slcons1(k,nnn))
      end do

      slden    (nnn) =  soil(nnn)%slden    
   end do

   !----- Defining some snow thickness variables ------------------------------------------!
   stretch = 2.0
   do kzs = 1,nzs
      thik          = 1.0
      thicknet(kzs) = 0.0
      do k = 1,(kzs+1)/2
         thick(k,kzs)       = thik
         thick(kzs+1-k,kzs) = thik
         thicknet(kzs)      = thicknet(kzs) + 2. * thik
         thik               = thik * stretch
      end do
      if (mod(kzs, 2) .ne. 0) thicknet(kzs) = thicknet(kzs) - thik/stretch
      do k = 1,kzs
         thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
      end do
   end do

   !----- Assigning some soil grid-dependent RK4 variables --------------------------------!
   rk4min_sfcw_mass  = rk4min_sfcw_moist * wdns   * dslz(nzg)
   rk4min_virt_water = rk4min_virt_moist * wdns   * dslz(nzg)


   !----- Organic Soil vertical grid spacing arrays (some with timestep info). ------------!
   olz (nzl+1) = 0.
   olz8(nzl+1) = 0.d0

   do k = 1,nzl
      dolz    (k) = olz(k+1) - olz(k)
      dolzo2  (k) = .5 * dolz(k)
      dolzi   (k) = 1. / dolz(k)
      dolzidt (k) = dolzi(k) * dtlsm
      olzt    (k) = .5 * (olz(k) + olz(k+1))
      olz8    (k) = dble(olz    (k))
      dolz8   (k) = dble(dolz   (k))
      dolzo28 (k) = dble(dolzo2 (k))
      dolzi8  (k) = dble(dolzi  (k))
      dolzidt8(k) = dble(dolzidt(k))
      olzt8   (k) = dble(olzt   (k))
   end do

   !----- Find the exponential increase factor to estimate the bottom boundary condition. -!
   oezg  = log(olz(1)/olz(nzl)) / log(real(nzl))
   olz0 = olz(1) * (real(nzl+1)/real(nzl))**oezg
   !---------------------------------------------------------------------------------------!


   !----- Find the thickness of the bottom boundary condition layer. ----------------------!
   dolz    (0) = olz(1) - olz0
   dolzo2  (0) = .5 * dolz(0)
   dolzi   (0) = 1. / dolz(0)
   dolzidt (0) = dolzi(0) * dtlsm
   dolz8   (0) = dble(dolz   (0))
   dolzo28 (0) = dble(dolzo2 (0))
   dolzi8  (0) = dble(dolzi  (0))
   dolzidt8(0) = dble(dolzidt(0))
   !---------------------------------------------------------------------------------------!


   !----- Find the height at the middle of the bottom boundary condition. -----------------!
   olzt   (0) = .5 * (olz0 + olz(1))
   olzt8  (0) = dble(olzt(0))
   !---------------------------------------------------------------------------------------!

   do k = 1,nzl
      dolzt    (k) = olzt(k) - olzt(k-1)
      dolzti   (k) = 1. / dolzt(k)
      dolztidt (k) = dolzti(k) * dtlsm
      dolzt8   (k) = dble(dolzt   (k))
      dolzti8  (k) = dble(dolzti  (k))
      dolztidt8(k) = dble(dolztidt(k))
   end do


   !----- Soil constants. -----------------------------------------------------------------!
   refdepth = -0.5

         ohydraul(:) = log (soil(12)%slcons / soil(12)%slcons0) / refdepth
      do k = 0,nzg
      
         select case (ipercol)
         case (0,1)
            !----- Original form, constant with depth.  -----------------------------------!
            olcons1(k,:) = soil(12)%slcons
         case (2)
            !------------------------------------------------------------------------------!
            !    TOPMODEL form, similar to CLM.  Here we use the same definition of slcons !
            ! from Cosby et al. (1984) because it has a stronger spread and it accounts    !
            ! for sand and clay contents.                                                  !
            !------------------------------------------------------------------------------!
            olcons1(k,:) = soil(12)%slcons * exp ( - olzt(k) / refdepth)
         end select

         !------ Find the double precision. -----------------------------------------------!
         olcons18(k,:) = dble(olcons1(k,:))
      end do

      olden    (:) =  soil(12)%slden    


   return
end subroutine sfcdata_ed
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine prints the soil properties for a single-polygon run.                 !
!------------------------------------------------------------------------------------------!
subroutine print_soil_info(cgrid,ifm)
   use ed_state_vars  , only : edtype            & ! structure
                             , polygontype       ! ! structure
   use mem_polygons   , only : n_poi             ! ! intent(in)
   use soil_coms      , only : isoilflg          & ! intent(in)
                             , soil              & ! intent(in)
                             , slxclay           & ! intent(in)
                             , slxsand           ! ! intent(in)
   use grid_coms      , only : nzg               ! ! intent(in)
   use ed_misc_coms   , only : sfilout           ! ! intent(in)
   use ed_max_dims    , only : str_len
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)          , target     :: cgrid
   integer               , intent(in) :: ifm
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)     , pointer    :: cpoly
   character(len=str_len)             :: polyname
   logical                            :: prescribed
   integer                            :: nsoil
   integer                            :: ipy
   integer                            :: isi
   integer                            :: slash
   integer                            :: endstr
   integer                            :: k
   !---------------------------------------------------------------------------------------!
   
   if (ifm /=1 .or. n_poi /= 1 .or. cgrid%npolygons /= 1) return
   ipy = 1
   cpoly => cgrid%polygon(ipy)

   !----- Find the polygon name. ----------------------------------------------------------!
   slash    = index(sfilout,'/',back=.true.) + 1
   endstr   = len_trim(sfilout)
   polyname = sfilout(slash:endstr)

   !----- Find whether the soil type characteristics were re-defined. ---------------------!
   prescribed = isoilflg(ifm)==2 .and. slxclay > 0. .and. slxsand > 0. .and.               &
                (slxclay + slxsand) <= 1.

   write (unit=*,fmt='(a)')           ' '
   write (unit=*,fmt='(a)') '   --------------------------------------------------------'
   write (unit=*,fmt='(a)')           '    Soil information:'
   write (unit=*,fmt='(a)')           ' '
   write (unit=*,fmt='(a,1x,a)')      '    Polygon name               :',trim(polyname)
   write (unit=*,fmt='(a,1x,f11.3)')  '    Longitude                  :',cgrid%lon(ipy)
   write (unit=*,fmt='(a,1x,f11.3)')  '    Latitude                   :',cgrid%lat(ipy)
   write (unit=*,fmt='(a,11x,l1)')    '    Prescribed sand and clay   :',prescribed
   write (unit=*,fmt='(a,1x,i11)')    '    # of sites                 :',cpoly%nsites

   do isi = 1,cpoly%nsites
      write (unit=*,fmt='(a,1x,i11)')   '    Site :',isi
    do k = 1,nzg
      nsoil = cpoly%ntext_soil(k,isi)
      write (unit=*,fmt='(a,1x,i5)')          '    Level: ',k
      write (unit=*,fmt='(a,1x,i10)')   '      - Type :',nsoil
      write(unit=*,fmt='(a,1x,es12.5)') '      - Clay fraction  =', soil(nsoil)%xclay
      write(unit=*,fmt='(a,1x,es12.5)') '      - Sand fraction  =', soil(nsoil)%xsand
      write(unit=*,fmt='(a,1x,es12.5)') '      - Silt fraction  =', soil(nsoil)%xsilt
      write(unit=*,fmt='(a,1x,es12.5)') '      - SLBS           =', soil(nsoil)%slbs
      write(unit=*,fmt='(a,1x,es12.5)') '      - SLPOTS         =', soil(nsoil)%slpots
      write(unit=*,fmt='(a,1x,es12.5)') '      - SLCONS         =', soil(nsoil)%slcons
      write(unit=*,fmt='(a,1x,es12.5)') '      - Dry air soil   =', soil(nsoil)%soilcp
      write(unit=*,fmt='(a,1x,es12.5)') '      - Wilting point  =', soil(nsoil)%soilwp
      write(unit=*,fmt='(a,1x,es12.5)') '      - Field capacity =', soil(nsoil)%sfldcap
      write(unit=*,fmt='(a,1x,es12.5)') '      - Saturation     =', soil(nsoil)%slmsts
      write(unit=*,fmt='(a,1x,es12.5)') '      - Heat capacity  =', soil(nsoil)%slcpd
    end do
   end do
   write (unit=*,fmt='(a)') '   --------------------------------------------------------'
   write (unit=*,fmt='(a)') ' '


   return
end subroutine print_soil_info
!==========================================================================================!
!==========================================================================================!

subroutine read_site_file(cgrid,igr)
   ! function for loading sites within polygons and initializing polygon parms
   ! call prior to loading pss/css files but after basic polygon established
   use soil_coms, only: soil,slz
   use grid_coms, only: nzg
   use ed_misc_coms, only: sfilin, vary_elev, vary_rad, vary_hyd
   use mem_polygons, only: edres
   use ed_state_vars, only: edtype,polygontype,sitetype,allocate_polygontype
   use ed_max_dims, only: max_site,n_pft,str_len

   implicit none
   integer :: igr
   logical :: no_rad   = .false.  !! true turns effect OFF
   logical :: no_lapse = .false.
   logical :: no_hyd   = .false.

   character(len=str_len) :: site_name,pss_name,css_name
  
   type(edtype) :: cgrid
   type(polygontype),pointer :: cpoly

   character(len=200) :: cdummy,cdummy2
   integer :: i,nsc,nsites=1
   real(kind=8) :: area_sum = 0.0d+0
   integer :: sc             ! soil classa
   integer :: fformat = 0    ! file format
   logical :: fexist         ! file exists
   real :: Te,T0,K0
   real :: fa,fb,zmin=-2.0
   integer :: sitenum,get_site_line,get_mat_line,get_header,found_mat_header=0,lcount=0,mcount=0
   real :: area,TCI,elevation,slope,aspect
   integer,allocatable :: soilclass(:)
   integer :: ipy,isi
   integer :: ierr

   if(vary_elev == 0)  no_lapse = .true.
   if(vary_rad == 0)   no_rad = .true.
   if(vary_hyd == 0) no_hyd = .true.




   ! ASSUMING FOR NOW THAT THERE IS NO WATER SITE

   ! init values
   if (associated(cgrid%load_adjacency)) cgrid%load_adjacency = 0

  
   do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)

      call create_ed10_ed20_fname(cgrid%lat(ipy), edres, cgrid%lon(ipy) &
                                 ,trim(sfilin(igr)),pss_name,css_name,site_name)
      !! check if site file exists
      inquire(file=trim(site_name),exist=fexist)
      if(.not.fexist) then
         print*,"error opening site file ",site_name
         print*,"setting ied_init_mode to 2 and loading as single site"
      end if

      
      ! If there is no terrestrial site
      if(.not. fexist) then  !! have no site file or ied_init_mode is not 3

         ! Allocate single site vector information to
         ! the swap polygon
         
         call allocate_polygontype(cpoly,1)
         call soil_default_fill(cgrid,igr,ipy)

         cpoly%area(1) = 1.0             ! Initialize the area to all

         ! Set soil moisture decay function, based on second layer's K value
         ! use the second layer instead of the top in case top is organic/peat
         sc = cpoly%ntext_soil(nzg-1,1)
         cpoly%moist_f(1) = -log(soil(sc)%slcons / soil(sc)%slcons0) / 2.0

         !! derive adjustments to f
         zmin = slz(cpoly%lsl(1))
         fa = -1.0/zmin !! should be 1/(depth to bedrock)
         if(cpoly%moist_f(1)*zmin < 0.0) then
            fb = cpoly%moist_f(1)/(1.0-exp(cpoly%moist_f(1)*zmin))
            cpoly%moist_f(1) = max(fa,fb)
         else
            cpoly%moist_f(1) = fa
         endif

         cpoly%sitenum(1)   = 1
         cpoly%elevation(1) = 0.0
         cpoly%slope(1)     = 0.0
         cpoly%aspect(1)    = 0.0
         cpoly%TCI(1)       = 0.0

      else

         !! Read data from site file
         open(unit=12,file=trim(site_name),form='formatted',status='old')

         !/* read file format line */
         read(unit=12,fmt=*)cdummy,nsites,cdummy2,fformat
         !/*line format: "nsites <nsites> format <fformat>" */   
         
!         print*,"reading",nsites,"sites using file format",fformat

         call allocate_polygontype(cpoly,nsites)
         call soil_default_fill(cgrid,igr,ipy)
     
         if(fformat <=0 .or. fformat > 3) then
            print*,""
            print*,"ERROR :: unrecognized file format specifier in",site_name
            stop
         endif
         if(fformat == 2) cgrid%load_adjacency(ipy) = 1
         nsc = 1
         if(fformat == 3) nsc = nzg

         !/* discard file header line */
         read(unit=12,fmt=*)
         
         !/* read data*/
         allocate(soilclass(nzg))
         lcount = 0 
         mcount = 0
         found_mat_header = 0
         isi = 0
         area_sum = 0.0d+0
         count_sites: do

            isi = isi + 1  ! The site counter
            
            !           read(12,*,iostat=ierr)time,pname,trk,age,area,fsc,stsc,  &
            !                stsl,ssc,psc,msn,fsn,water(1:nwater)
            !           if(ierr /= 0)exit count_patches
            
            lcount = lcount + 1 !! line counter
            get_site_line=0
            get_mat_line=0
            get_header=0 !! line flags
            
            !/********* decide what type of line to read ***************/
            select case (fformat)
            case (1,3)
               if(lcount <= nsites) then !we're reading data
                  get_site_line=1
                  get_mat_line=0
                  get_header=0
               else                      !we're past data, discard
                  get_site_line=0
                  get_mat_line=0
                  get_header=1             
               endif
            case (2)
               if(lcount <= nsites) then  !know we're in the site section
                  get_site_line=1
                  get_mat_line=0

               elseif(found_mat_header.eq.1)then !know we're in the matrix section
                  if(mcount < nsites)then
                     get_site_line = 0
                     get_mat_line=1
                     get_header = 0
                  else  !//already got all the matrix lines, discard what's left
                     get_site_line = 0
                     get_mat_line = 0
                     get_header = 1
                  endif
               else 
                  
                  !!assume line is header, remove and then assume header found
                  
                  get_site_line = 0
                  get_mat_line = 0
                  get_header = 1
                  found_mat_header=1
                  
               endif
               
            end select

!            print*,"line indicators",get_site_line,get_mat_line,get_header
    
            if(get_site_line == 1)then   !/********** READ SITE LINE ***************/
               
               !/* line format: sitenum, area, TCI, elevation, slope, aspect,ntext_soil(s) */              
               read(unit=12,fmt=*,iostat=ierr)sitenum,area,TCI,elevation,slope,aspect,soilclass(1:nsc)
               if(ierr == 0) then
                  !/*create data object for each new site */
!                  print*,sitenum, area, TCI, elevation,slope,aspect,soilclass(1:nsc)
                  cpoly%area(isi) = area            ! Initialize the area to all

                  area_sum = area_sum + dble(area)
                  cpoly%sitenum(isi)      = sitenum
                  cpoly%elevation(isi)    = elevation
                  cpoly%slope(isi)        = slope
                  cpoly%aspect(isi)       = aspect
                  cpoly%TCI(isi)          = TCI+13.96962

                  !! flags to turn effects off
                  if(no_lapse) cpoly%elevation(isi) = 0.0
                  if(no_hyd)   cpoly%TCI(isi)       = 8.0
                  if(no_rad) then
                               cpoly%slope(isi)     = 0.0
                               cpoly%aspect(isi)    = 0.0
                  end if

!print*,"SITE",cpoly%elevation(isi),cpoly%TCI(isi),cpoly%slope(isi),cpoly%aspect(isi)

                  if(fformat == 3) then
                     do i=1,nzg
                        cpoly%ntext_soil(i,isi) = soilclass(i)
                     end do
                  else
                     do i=1,nzg
                        cpoly%ntext_soil(i,isi) = soilclass(i)
                     end do
                  end if
                  
                  !//Currently do nothing with setting site-level soils

                  sc = cpoly%ntext_soil(nzg-1,1)
                  cpoly%moist_f(isi) = -log(soil(sc)%slcons / soil(sc)%slcons0) / 2.0
                  !! derive adjustments to f
                  zmin = slz(cpoly%lsl(isi))
                  fa = -1.0/zmin !! should be 1/(depth to bedrock)
                  if(cpoly%moist_f(isi)*zmin < 0.0) then
                     fb = cpoly%moist_f(isi)/(1.0-exp(cpoly%moist_f(isi)*zmin))
                     cpoly%moist_f(isi) = max(fa,fb)
                  else
                     cpoly%moist_f(isi) = fa
                  endif
               end if  ! end valid line
            end if        !//********** END READ SITE LINE**************

            ! RGK
            if(get_header == 1) read(unit=12,fmt=*,iostat=ierr) !/* discard file line */
            ! RGK
            if(get_mat_line == 1) then

               !/*read line into site adjacency matrix*/
               
               read(unit=12,fmt=*,iostat=ierr) cgrid%site_adjacency(mcount,1:(nsites+1),ipy)

               mcount = mcount+1
               
            endif
            if(ierr /= 0) exit count_sites
        
         end do count_sites
         deallocate(soilclass)
      end if

      !adjust areas
      !assume that if area_sum ~ 1 need to renormalize terrestrial
      if( area_sum > 0.995d+0) then
         cpoly%area(:) = real(dble(cpoly%area(:))/area_sum)
      end if


      if(cgrid%load_adjacency(ipy) /= 0) then  
         call calc_flow_routing(cgrid,ipy)
      endif

      ! calculate summary stats - pass 1: Te
      ! On terrestrial site still
      Te = 0.0 
      area_sum = 0.0d+0
      do isi = 1,cpoly%nsites
         sc = cpoly%ntext_soil(nzg-1,isi)
         K0 = soil(sc)%slcons0
         T0 = K0/cpoly%moist_f(isi)
         Te = Te + T0*cpoly%area(isi)
         area_sum = area_sum + dble(cpoly%area(isi))
      end do

      Te = Te/real(area_sum)
      cgrid%Te(ipy) = Te

      !pass 2: W, Wbar

      cgrid%wbar(ipy) = 0.0

      do isi = 1,cpoly%nsites
         sc = cpoly%ntext_soil(nzg-1,isi)
         K0 = soil(sc)%slcons0
         T0 = K0/cpoly%moist_f(isi)
         cpoly%moist_W(isi) = cpoly%TCI(isi) + log(Te) - log(T0)
         cgrid%wbar(ipy) = cgrid%wbar(ipy) + real(dble(cpoly%moist_W(isi))*dble(cpoly%area(isi))/area_sum)
      end do

      !call dump_ed(myPolygon)

   end do


end subroutine read_site_file

!==========================================================================================!
!==========================================================================================!
subroutine calc_flow_routing(cgrid,ipy)
   use ed_state_vars, only: polygontype, edtype
   type(edtype), target :: cgrid ! Alias for current grid
   integer, intent(in) :: ipy    ! Current polygon Polygon ID
   integer :: i,ihys,ines ! ihys -> current hydro site ID. ines-> next hydro site ID
   type(polygontype),pointer :: cpoly ! Alias for current polygon
   real :: side = 10.0  !cell side length (m) used for calculating adjacency
   real(kind=8) :: row_sum
   integer, external :: find_rank
   integer, allocatable, dimension(:) :: hyrank
   integer :: nsites

   !if we loaded an adjacency matrix from file, we don't need to calculate flow routing
   if(cgrid%load_adjacency(ipy) == 1) return   

   cpoly => cgrid%polygon(ipy)
   nsites=cpoly%nsites
   
   allocate (hyrank(nsites))

   ! set flow routing -> sort by TCI, it will give the sites the hydro order
   call rank_up(nsites,cpoly%TCI,hyrank)


   !init structures and recalculate sitenums list based on hydro order
   do ihys=1,nsites
      do ines=1,nsites
         cgrid%site_adjacency(ihys,ines,ipy) = 0.0
      end do
   end do
   !routing established, now calc approximate adjacency matrix
   do i=1,nsites-1
      ihys=find_rank(i,nsites,hyrank)
      ines=find_rank(i+1,nsites,hyrank)
      cgrid%site_adjacency(ihys,ihys,ipy) = 2.0*side*(side-1.0)*cpoly%area(ihys) &
               + side/2.0*cpoly%area(ihys)/(cpoly%area(ihys)+cpoly%area(ines))
      cgrid%site_adjacency(ihys,ines,ipy) = side*cpoly%area(ines)/(cpoly%area(ihys)+cpoly%area(ihys))
   end do
   !Now find the adjacency for the last rank site:
   ihys=find_rank(nsites,nsites,hyrank)
   ines=find_rank(1,nsites,hyrank)
   cgrid%site_adjacency(ihys,ihys,ipy) = 2.0*side*(side-1.0)*cpoly%area(ihys) &
            + side/2.0*cpoly%area(ihys)/(cpoly%area(ihys)+cpoly%area(ines))
   cgrid%site_adjacency(ihys,ines,ipy) = side*cpoly%area(ines)/(cpoly%area(ihys)+cpoly%area(ihys))

   !normalize routing (rows sum to 1)  
   do ihys=1,nsites
      row_sum=sum(dble(cgrid%site_adjacency(ihys,:,ipy)))
      cgrid%site_adjacency(ihys,:,ipy)=real(dble(cgrid%site_adjacency(ihys,:,ipy))/row_sum)
   end do

   !check routing
!   print*,"CHECK ROUTING" 
!   do i=1,myPolygon%nsites
!      print*,cgrid%site_adjacency(i,:,ipy)
!   end do

   deallocate(hyrank)
   return
end subroutine calc_flow_routing

end module ed_init
