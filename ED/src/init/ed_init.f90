!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign the longitude, latitude, and soil class for all non-empty !
! polygons.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine set_polygon_coordinates()
   use grid_coms     , only : ngrids    & ! intent(in)
                            , nzg       ! ! intent(in)
   use ed_work_vars  , only : work_v    ! ! structure
   use ed_node_coms  , only : mynum     ! ! intent(in)
   use ed_state_vars , only : edgrid_g  & ! structure
                            , gdpy      ! ! intent(in)
   implicit none 
   !----- Local variables -----------------------------------------------------------------!
   integer                :: ifm
   integer                :: ipy
   integer                :: npoly
   !---------------------------------------------------------------------------------------!

   gloop: do ifm=1,ngrids

      npoly=gdpy(mynum,ifm)

      ploop: do ipy=1,npoly
         edgrid_g(ifm)%lon(ipy)              = work_v(ifm)%glon   (ipy)
         edgrid_g(ifm)%lat(ipy)              = work_v(ifm)%glat   (ipy)
         edgrid_g(ifm)%xatm(ipy)             = work_v(ifm)%xid    (ipy)
         edgrid_g(ifm)%yatm(ipy)             = work_v(ifm)%yid    (ipy)

         !----- Assign the commonest soil type to the polygon. ----------------------------!
         edgrid_g(ifm)%ntext_soil(1:nzg,ipy) = work_v(ifm)%ntext(1,ipy)
         !---------------------------------------------------------------------------------!
      end do ploop

      

   end do gloop

    
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
   use ed_node_coms  , only : mynum                ! ! intent(in)
   use ed_state_vars , only : edgrid_g             & ! structure
                            , edtype               & ! structure
                            , polygontype          & ! structure
                            , allocate_polygontype ! ! subroutine
   use soil_coms     , only : soil                 & ! intent(in)
                            , slz                  ! ! intent(in)
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
   integer                    :: npoly
   integer                    :: nsite
   integer                    :: k
   real                       :: sc
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


               !------ Initialise the lowest soil layer. ----------------------------------!
               cpoly%lsl(isi) = cgrid%lsl(ipy)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Use the soil type and populate the site-level soil texture.           !
               !---------------------------------------------------------------------------!
               do k=1,nzg
                  cpoly%ntext_soil(k,isi) = work_v(ifm)%ntext(itext,ipy)
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
!     This subroutine fills the lsl variables based on the soil_depth file.  In case       !
! isoildepthflg was zero, then the layer_index matrix was filled with zeroes, so we do not !
! need to worry about this here.                                                           !
!------------------------------------------------------------------------------------------!
subroutine soil_depth_fill(cgrid,igr)
   
   use soil_coms     , only : layer_index ! ! intent(in)
   use ed_state_vars , only : edtype      ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype) , target     :: cgrid
   integer      , intent(in) :: igr
   !----- Local variables -----------------------------------------------------------------!
   integer                   :: ilat_bin
   integer                   :: ilon_bin
   integer                   :: ipy
   !---------------------------------------------------------------------------------------!

   do ipy = 1,cgrid%npolygons
      ilat_bin = min(180,int(90.0 - cgrid%lat(ipy)) + 1)
      ilon_bin = int(180.0 + cgrid%lon(ipy)) + 1

      !------------------------------------------------------------------------------------!
      !    Require at least 2 layers.  This requirement was taken in consideration when    !
      ! layer_index was filled at the first initialization, so it is safe to just copy.    !
      !------------------------------------------------------------------------------------!
      cgrid%lsl(ipy) =layer_index(ilat_bin,ilon_bin) 
   end do

   return
end subroutine soil_depth_fill
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
   use phenology_coms    , only : iphen_scheme    ! ! intent(in)
   use ed_misc_coms      , only : ied_init_mode   ! ! intent(in)
   use phenology_startup , only : phenology_init  ! ! intent(in)
   use ed_node_coms      , only : mynum           & ! intent(in)
                                , nmachs          & ! intent(in)
                                , nnodetot        & ! intent(in)
                                , mchnum          & ! intent(in)
                                , machs           & ! intent(in)
                                , master_num      & ! intent(in)
                                , sendnum         & ! intent(in)
                                , recvnum         ! ! intent(in)
   use grid_coms         , only : ngrids          ! ! intent(in)
   use ed_state_vars     , only : edgrid_g        ! ! structure

   implicit none
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer                :: ierr
   integer                :: igr
   integer                :: ping 
   !---------------------------------------------------------------------------------------!

   ping = 741776


   if (mynum == 1) write(unit=*,fmt='(a)') ' + Doing sequential initialization over nodes.'

   !---------------------------------------------------------------------------------------!
   ! STEP 1: Find lowest soil layer for each site (derived from soil depth).               !
   !---------------------------------------------------------------------------------------!
   do igr=1,ngrids
      call ed_newgrid(igr)
      call soil_depth_fill(edgrid_g(igr),igr)
   end do


   !---------------------------------------------------------------------------------------!
   ! STEP 2: Read in Site files and initialize hydrologic adjacencies.                     !
   !---------------------------------------------------------------------------------------!
   if (mynum /= 1) &
      call MPI_Recv(ping,1,MPI_INTEGER,recvnum,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  
   select case (ied_init_mode)
   case (3)
      !----- Hydrology run.  Use specific scheme. -----------------------------------------!
      do igr = 1,ngrids
         call read_site_file(edgrid_g(igr),igr)
      end do
   case (4)
      continue
   case default
         call set_site_defprops()
   end select
  
   if (mynum < nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,100,MPI_COMM_WORLD,ierr)
  
   !---------------------------------------------------------------------------------------!
   ! STEP 3: Do ASCII type restart initialization of site patch and cohort biophysical     !
   !         states.                                                                       !
   !---------------------------------------------------------------------------------------!
   if (mynum /= 1) &
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  

   select case (ied_init_mode)
   case(-8,-1,0)
      !----- Initialize everything with near-bare ground ----------------------------------!
      if (mynum /= 1) write(unit=*,fmt='(a)') ' + Doing near bare ground initialization...'
      do igr=1,ngrids
           call near_bare_ground_init(edgrid_g(igr))
      end do

   case(1,2,3,6)
      !----- Initialize with ED1-type restart information. --------------------------------!
      write(unit=*,fmt='(a,i3.3)') ' + Initializing from ED restart file. Node: ',mynum
      call read_ed10_ed20_history_file

   case(4)   
      write(unit=*,fmt='(a,i3.3)') ' + Initializing from ED2.1 state file. Node: ',mynum
      call read_ed21_history_file

   case(5,99)
      write(unit=*,fmt='(a,i3.3)')                                                         &
          ' + Initializing from a collection of ED2.1 state files. Node: ',mynum
      call read_ed21_history_unstruct

   end select

   if (mynum < nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,101,MPI_COMM_WORLD,ierr)

   !---------------------------------------------------------------------------------------!
   ! STEP 4: Initialize phenology parameters and thermal sums.                             !
   !---------------------------------------------------------------------------------------!
   if (mynum /= 1) &
      call MPI_Recv(ping,1,MPI_INTEGER,recvnum,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

   write(unit=*,fmt='(a,i3.3)') ' + Initializing phenology. Node: ',mynum
   call phenology_init()

   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,102,MPI_COMM_WORLD,ierr)
  
  
   !---------------------------------------------------------------------------------------!
   ! STEP 5: Initialize anthropogenic disturbance.                                         !
   !---------------------------------------------------------------------------------------!
   if (mynum /= 1) &
     call MPI_Recv(ping,1,MPI_INTEGER,recvnum,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

   write(unit=*,fmt='(a,i3.3)')                                                            &
      ' + Initializing anthropogenic disturbance forcing. Node: ',mynum

   call landuse_init()

   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,103,MPI_COMM_WORLD,ierr)

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
                          , nzs               ! ! intent(in)
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
                          , emisg             & ! intent(out)
                          , soil              & ! intent(in)
                          , thicknet          & ! intent(out)
                          , thick             ! ! intent(out)
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
   !---------------------------------------------------------------------------------------!


   !----- Soil vertical grid spacing arrays (some with timestep info). --------------------!
   slz(nzg+1) = 0.

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

   !----- Find the exponential increase factor. -------------------------------------------!
   ezg  = log(slz(1)/slz(nzg)) / log(real(nzg))
   slz0 = slz(1) * (real(nzg+1)/real(nzg))**ezg


   slzt (0) = .5 * (slz0 + slz(1))
   slzt8(0) = dble(slzt(0))

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
      emisg (nnn) = .98
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
      if ((kzs+1)/2 .ne. kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
      do k = 1,kzs
         thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
      end do
   end do

   !----- Assigning some soil grid-dependent RK4 variables --------------------------------!
   rk4min_sfcw_mass  = rk4min_sfcw_moist * wdns   * dslz(nzg)
   rk4min_virt_water = rk4min_virt_moist * wdns   * dslz(nzg)

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
      nsoil = cpoly%ntext_soil(nzg,isi)
      write (unit=*,fmt='(a)')          ' '
      write (unit=*,fmt='(a,1x,i11)')   '    Site :',isi
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
   write (unit=*,fmt='(a)') '   --------------------------------------------------------'
   write (unit=*,fmt='(a)') ' '


   return
end subroutine print_soil_info
!==========================================================================================!
!==========================================================================================!
