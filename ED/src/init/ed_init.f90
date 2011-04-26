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
         edgrid_g(ifm)%lon(ipy)              = work_v(ifm)%glon(ipy)
         edgrid_g(ifm)%lat(ipy)              = work_v(ifm)%glat(ipy)
         edgrid_g(ifm)%ntext_soil(1:nzg,ipy) = work_v(ifm)%ntext(ipy)
         edgrid_g(ifm)%xatm(ipy)             = work_v(ifm)%xid(ipy)
         edgrid_g(ifm)%yatm(ipy)             = work_v(ifm)%yid(ipy)

      end do ploop

      

   end do gloop

    
   return
end subroutine set_polygon_coordinates
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
   case (4)
      continue
   case default
      do igr = 1,ngrids
         call read_site_file(edgrid_g(igr),igr)
      end do
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

   do k = 2,nzg
      dslzt    (k) = slzt(k) - slzt(k-1)
      dslzti   (k) = 1. / dslzt(k)
      dslztidt (k) = dslzti(k) * dtlsm
      dslzt8   (k) = dble(dslzt   (k))
      dslzti8  (k) = dble(dslzti  (k))
      dslztidt8(k) = dble(dslztidt(k))
   end do

   !----- These must be defined for free drainage bc (RGK) --------------------------------!
   dslzt    (1) = 2.0*slz(1) - slzt(1)
   dslzti   (1) = 1./dslzt(1)
   dslztidt (1) = dslzti(1) * dtlsm
   dslzt8   (1) = dble(dslzt   (1))
   dslzti8  (1) = dble(dslzti  (1))
   dslztidt8(1) = dble(dslztidt(1))

   !----- Soil constants. -----------------------------------------------------------------!
   refdepth = -2.0

   do nnn = 1,ed_nstyp
      if (nnn /= 13) then
         fhydraul(nnn) = log (soil(nnn)%slcons / soil(nnn)%slcons0) / refdepth
      else
         fhydraul(nnn) = 0.
      end if
      do k = 1,nzg
         slcons1(k,nnn) = soil(nnn)%slcons     ! ORIGINAL form - const with depth
   !     slcons1(k,nnn) = soilparms(5,nnn)  &  ! TOPMODEL form - large at surface
   !        * exp(slz(k) * fhydraul(nnn))      !    and exp decrease with depth
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
