!==========================================================================================!
!==========================================================================================!
!    Subroutines based on the RAMS node decomposition. The main difference between the     !
! original code and this one is that when we split the domain we need to consider whether  !
! the polygon will fall on land or water. The water ones will be removed, so this should   !
! be taken into account for the standalone version.                                        !
!------------------------------------------------------------------------------------------!
subroutine ed_node_decomp(init,standalone,masterworks)

   use grid_coms   , only : ngrids            & ! intent(in)
                          , nnxp              & ! intent(in)
                          , nnyp              ! ! intent(in)
   use ed_node_coms, only : mmxp              & ! intent(out)
                          , mmyp              & ! intent(out)
                          , mia               & ! intent(in)
                          , miz               & ! intent(in)
                          , mja               & ! intent(in)
                          , mjz               & ! intent(in)
                          , mi0               & ! intent(in)
                          , mj0               & ! intent(in)
                          , mibcon            ! ! intent(in)
   use ed_para_coms, only : nmachs            ! ! intent(in)
   use mem_polygons, only : n_ed_region       & ! intent(in)
                          , maxsite           ! ! intent(in)
   use ed_work_vars, only : work_e            & ! intent(in)
                          , work_v            & ! intent(in)
                          , ed_alloc_work     & ! subroutine
                          , ed_nullify_work   ! ! subroutine
   use soil_coms   , only : isoilflg          ! ! subroutine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: init
   logical, intent(in) :: standalone
   logical, intent(in) :: masterworks
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ngr
   integer             :: nsiz
   integer             :: ntotmachs
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    This is a logical flag to test wheter the master should also simulate polygons.    !
   ! This is true for offline runs, but false for coupled runs.                            !
   !---------------------------------------------------------------------------------------!
   if (masterworks) then
      ntotmachs=nmachs+1
   else
      ntotmachs=nmachs
   end if
     
   allocate(work_e(ngrids),work_v(ngrids))

   !---------------------------------------------------------------------------------------!
   !      Decompose all grids into subdomains.                                             !
   !---------------------------------------------------------------------------------------!
   do ngr = 1,ngrids
      !------------------------------------------------------------------------------------!
      !     POI grids always have one point only. Since the structure will be sent.  I am  !
      ! filling the structures.  Only 4 extra numbers will be sent through MPI if we do    !
      ! this.                                                                              !
      !------------------------------------------------------------------------------------!
      mmxp(ngr) = nnxp(ngr)
      mmyp(ngr) = nnyp(ngr)

      call ed_nullify_work(work_e(ngr))
      call ed_alloc_work(work_e(ngr),nnxp(ngr),nnyp(ngr),maxsite)
   end do

   call get_grid()

   do ngr = 1,ngrids!n_ed_region

      !------------------------------------------------------------------------------------!
      !      Obtain estimates of the fraction of computational time (work) required for    !
      ! each column in the region of the domain.                                           !
      !------------------------------------------------------------------------------------!
      call get_work(ngr,mmxp(ngr),mmyp(ngr))
      call ed_parvec_work(ngr,mmxp(ngr),mmyp(ngr))
   end do

   !----- Check whether we have a good first guess of how different polygons work. --------!
   call ed_load_work_from_history()

   !---------------------------------------------------------------------------------------!
   !     Polygons of interest (POI) are single polygons, so they get the simplest          !
   ! configuration: full work, over land.                                                  !
   !---------------------------------------------------------------------------------------!
   do ngr=n_ed_region+1,ngrids
      call ed_newgrid(ngr)
      work_e(ngr)%work(1,1)=1.
      work_e(ngr)%land(1,1)=.true.
      call ed_parvec_work(ngr,mmxp(ngr),mmyp(ngr))
   end do

   return
end subroutine ed_node_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will assign the longitude and latitude of all polygon candidates in   !
! this run.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine get_grid
   use mem_polygons, only : grid_type      & ! intent(in)
                          , grid_res       & ! intent(in)
                          , n_ed_region    & ! intent(in)
                          , poi_lat        & ! intent(in)
                          , poi_lon        & ! intent(in)
                          , ed_reg_lonmin  & ! intent(in)
                          , ed_reg_latmin  ! ! intent(in)
   use grid_coms   , only : ngrids         & ! intent(in)
                          , nnxp           & ! intent(in)
                          , nnyp           & ! intent(in)
                          , nstratx        & ! intent(in)
                          , nstraty        ! ! intent(in)
   use ed_work_vars, only : work_e         ! ! intent(inout)
   implicit none

   !----- Local variables. ----------------------------------------------------------------!
   integer :: ifm
   integer :: i
   integer :: j
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Here we assign the longitude and latitude depending on how the user wants the     !
   ! grid to be defined (regular longitude/latitude grid, or polar-stereographic).         !
   !---------------------------------------------------------------------------------------!
   select case (grid_type)
   case (0) !----- Regular longitude/latitude grid. ---------------------------------------!
      do ifm=1,n_ed_region
         do i=1,nnxp(ifm)
            do j=1,nnyp(ifm)
               work_e(ifm)%glon(i,j) = ed_reg_lonmin(ifm)                                  &
                                     + (float(i) - 0.5) * grid_res / real(nstratx(ifm))
               work_e(ifm)%glat(i,j) = ed_reg_latmin(ifm)                                  &
                                     + (float(j) - 0.5) * grid_res / real(nstraty(ifm))
            end do
         end do
      end do

   case (1) !----- Polar-stereographic grid. ----------------------------------------------!
      if (n_ed_region > 0) call ed_gridset(1)
      do ifm=1,n_ed_region
         call ed_newgrid(ifm)
         call ed_polarst(nnxp(ifm),nnyp(ifm),work_e(ifm)%glat,work_e(ifm)%glon)
      end do
      
   case default !----- Eventually we'll have polygons, but not yet. -----------------------!
      call fatal_error('Invalid grid_type in ED_grid_setup.','get_grid','ed_para_init.f90')
   end select
   do ifm=n_ed_region+1,ngrids
      do i=1,nnxp(ifm)
         do j=1,nnyp(ifm)
            work_e(ifm)%glon(i,j)=poi_lon(ifm-n_ed_region)
            work_e(ifm)%glat(i,j)=poi_lat(ifm-n_ed_region)
         end do
      end do
   end do

   return
end subroutine get_grid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine get_work(ifm,nxp,nyp)

   use ed_work_vars, only : work_e         ! ! structure
   use soil_coms   , only : veg_database   & ! intent(in)
                          , soil_database  & ! intent(in)
                          , isoilflg       & ! intent(in)
                          , nslcon         ! ! intent(in)
   use mem_polygons, only : n_poi          & ! intent(in)
                          , poi_res        & ! intent(in)
                          , grid_res       & ! intent(in)
                          , grid_type      & ! intent(in)
                          , maxsite        ! ! intent(in)
   use ed_misc_coms, only : min_site_area  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ifm
   integer, intent(in) :: nxp
   integer, intent(in) :: nyp
   !----- Local variables. ----------------------------------------------------------------!
   integer :: npoly
   real   , dimension(:,:), allocatable :: lat_list
   real   , dimension(:,:), allocatable :: lon_list
   integer, dimension(:,:), allocatable :: leaf_class_list
   integer, dimension(:,:), allocatable :: ntext_soil_list
   real   , dimension(:,:), allocatable :: ipcent_land
   real   , dimension(:,:), allocatable :: ipcent_soil
   integer                              :: datsoil
   integer                              :: ipy
   integer                              :: i
   integer                              :: j
   integer                              :: jboff
   integer                              :: jtoff
   integer                              :: iloff
   integer                              :: iroff
   integer                              :: itext
   real                                 :: maxwork
   !---------------------------------------------------------------------------------------!


   !----- Assume that every grid cell will become a polygon. ------------------------------!
   npoly = nxp*nyp
   allocate(lat_list       (      3,npoly))
   allocate(lon_list       (      3,npoly))
   allocate(leaf_class_list(maxsite,npoly))
   allocate(ntext_soil_list(maxsite,npoly))
   allocate(ipcent_land    (maxsite,npoly))
   allocate(ipcent_soil    (maxsite,npoly))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Fill lat/lon lists.  The longitude and latitude are initially assigned as arrays. !
   ! The second index (the npoly one) is used to define each polygon.  The first index is  !
   ! used to define the relative coordinate of important points within that polygon:       !
   ! 1 means grid centre, 2 means Northwestern corner, and 3 means Southeastern corner.    !
   !---------------------------------------------------------------------------------------!
   if (n_poi > 0 .and. ifm <= n_poi) then

      ipy = 0
      do j = 1,nyp
         do i=1,nxp
            ipy = ipy + 1

            !----- Grid mid-point. --------------------------------------------------------!
            lon_list(1,ipy) = work_e(ifm)%glon(i,j)
            lat_list(1,ipy) = work_e(ifm)%glat(i,j)
            !----- Northwestern corner. ---------------------------------------------------!
            lon_list(2,ipy) = work_e(ifm)%glon(i,j) - 0.5 * poi_res(ifm)
            lat_list(2,ipy) = work_e(ifm)%glat(i,j) + 0.5 * poi_res(ifm)
            !----- Southeastern corner. ---------------------------------------------------!
            lon_list(3,ipy) = work_e(ifm)%glon(i,j) + 0.5 * poi_res(ifm)
            lat_list(3,ipy) = work_e(ifm)%glat(i,j) - 0.5 * poi_res(ifm)

            !----- Adjusting the longitudes to be between -180 and 180. -------------------!
            if (lon_list(1,ipy) >=  180.) lon_list(1,ipy) = lon_list(1,ipy) - 360.
            if (lon_list(1,ipy) <= -180.) lon_list(1,ipy) = lon_list(1,ipy) + 360.
            if (lon_list(2,ipy) >=  180.) lon_list(2,ipy) = lon_list(2,ipy) - 360.
            if (lon_list(2,ipy) <= -180.) lon_list(2,ipy) = lon_list(2,ipy) + 360.
            if (lon_list(3,ipy) >=  180.) lon_list(3,ipy) = lon_list(3,ipy) - 360.
            if (lon_list(3,ipy) <= -180.) lon_list(3,ipy) = lon_list(3,ipy) + 360.
            !------------------------------------------------------------------------------!

         end do
      end do
   else
      ipy = 0
      do j = 1,nyp
         do i=1,nxp
            ipy = ipy + 1

            !----- Grid mid-point. --------------------------------------------------------!
            lon_list(1,ipy) = work_e(ifm)%glon(i,j)
            lat_list(1,ipy) = work_e(ifm)%glat(i,j)

            !------------------------------------------------------------------------------!
            !     Here we find the corners differently depending on whether the grid is a  !
            ! regular longitude-latitude or a polar-stereographic...                       !
            !------------------------------------------------------------------------------!
            select case (grid_type)
            case (0) !----- Regular longitude/latitude grid. ------------------------------!
               !----- Northwestern corner. ------------------------------------------------!
               lon_list(2,ipy) = work_e(ifm)%glon(i,j) - 0.5 * grid_res
               lat_list(2,ipy) = work_e(ifm)%glat(i,j) + 0.5 * grid_res
               !----- Southeastern corner. ------------------------------------------------!
               lon_list(3,ipy) = work_e(ifm)%glon(i,j) + 0.5 * grid_res
               lat_list(3,ipy) = work_e(ifm)%glat(i,j) - 0.5 * grid_res
            case (1) !----- Polar-stereographic. ------------------------------------------!
               !---------------------------------------------------------------------------!
               !   Setting the offsets, accounting for the corners.  In case we are in one !
               ! of the edges, we must adjust the offset.                                  !
               !---------------------------------------------------------------------------!
               if (i == 1) then
                  iloff = -1
                  iroff =  1
               elseif (i == nxp) then
                  iloff =  1
                  iroff = -1
               else
                  iloff =  1
                  iroff =  1
               end if
               if (j == 1) then
                  jtoff =  1
                  jboff = -1
               elseif (j == nyp) then
                  jtoff = -1
                  jboff =  1
               else
                  jtoff =  1
                  jboff =  1
               end if
               !---------------------------------------------------------------------------!



               !----- Northwestern corner. ------------------------------------------------!
               lon_list(2,ipy) = work_e(ifm)%glon(i,j) + real(iloff) * 0.5                 &
                               * (work_e(ifm)%glon(i-iloff,j)-work_e(ifm)%glon(i,j))
               lat_list(2,ipy) = work_e(ifm)%glat(i,j) + real(jtoff) * 0.5                 &
                               * (work_e(ifm)%glat(i,j+jtoff) - work_e(ifm)%glat(i,j))
               !----- Southeastern corner. ------------------------------------------------!
               lon_list(3,ipy) = work_e(ifm)%glon(i,j) + real(iroff) * 0.5                 &
                               * (work_e(ifm)%glon(i+iroff,j) - work_e(ifm)%glon(i,j))
               lat_list(3,ipy) = work_e(ifm)%glat(i,j) + real(jboff) * 0.5                 &
                               * (work_e(ifm)%glat(i,j-jboff) - work_e(ifm)%glat(i,j))
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !----- Adjusting the longitudes to be between -180 and 180. -------------------!
            if (lon_list(1,ipy) >=  180.) lon_list(1,ipy) = lon_list(1,ipy) - 360.
            if (lon_list(1,ipy) <= -180.) lon_list(1,ipy) = lon_list(1,ipy) + 360.
            if (lon_list(2,ipy) >=  180.) lon_list(2,ipy) = lon_list(2,ipy) - 360.
            if (lon_list(2,ipy) <= -180.) lon_list(2,ipy) = lon_list(2,ipy) + 360.
            if (lon_list(3,ipy) >=  180.) lon_list(3,ipy) = lon_list(3,ipy) - 360.
            if (lon_list(3,ipy) <= -180.) lon_list(3,ipy) = lon_list(3,ipy) + 360.
            !------------------------------------------------------------------------------!

         end do
      end do
   end if

   !----- Generate the land/sea mask. -----------------------------------------------------!
   write(unit=*,fmt=*) ' => Generating the land/sea mask.'

   call leaf_database(trim(veg_database(ifm)),maxsite,npoly,'leaf_class'                   &
                     ,lat_list,lon_list,leaf_class_list,ipcent_land)

   if (isoilflg(ifm) == 1) then
      call leaf_database(trim(soil_database(ifm)),maxsite,npoly,'soil_text'                &
                        ,lat_list,lon_list,ntext_soil_list,ipcent_soil)
   else
      !------------------------------------------------------------------------------------!
      !   Allow for only one site by making the first site with the default soil type and  !
      ! area 1., and the others with area 0.                                               !
      !------------------------------------------------------------------------------------!
      ntext_soil_list        (:,:) = nslcon
      ipcent_soil            (:,:) = 0.
      ipcent_soil            (1,:) = 1.
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!



   !----- Re-map the land cover classes. --------------------------------------------------!
   ipy     = 0
   maxwork = epsilon(0.0)
   do j = 1,nyp
      do i=1,nxp
         ipy = ipy + 1
         work_e(ifm)%land(i,j) = ipcent_land(1,ipy) > min_site_area

         if (work_e(ifm)%land(i,j)) then
            work_e(ifm)%landfrac(i,j) = ipcent_land(1,ipy)

            work_e(ifm)%work(i,j) = 0.0
            do itext = 1,maxsite
               if (ipcent_soil(itext,ipy) > min_site_area) then
                  work_e(ifm)%work(i,j) = work_e(ifm)%work(i,j) + 1.
               end if
               work_e(ifm)%soilfrac(itext,i,j) = ipcent_soil(itext,ipy)
               work_e(ifm)%ntext   (itext,i,j) = ntext_soil_list (itext,ipy)
            end do
            maxwork = max(maxwork,work_e(ifm)%work(i,j))

         else
            !----- Making this grid point 100% water --------------------------------------!
            work_e(ifm)%landfrac  (i,j) = 0.
            work_e(ifm)%work      (i,j) = epsilon(0.0)
            work_e(ifm)%ntext   (:,i,j) = 0
            work_e(ifm)%soilfrac(:,i,j) = 0.
         end if
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !----- Normalise the workload. ---------------------------------------------------------!
   do j=1,nyp
      do i=1,nxp
         work_e(ifm)%work(i,j) = work_e(ifm)%work(i,j) / maxwork
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !------ De-allocate the temporary arrays. ----------------------------------------------!
   deallocate(lat_list       )
   deallocate(lon_list       )
   deallocate(leaf_class_list)
   deallocate(ntext_soil_list)
   deallocate(ipcent_land    )
   deallocate(ipcent_soil    )
   !---------------------------------------------------------------------------------------!

   return
end subroutine get_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_parvec_work(ifm,nxp,nyp)

   use ed_work_vars,  only : work_e              & ! intent(in)
                           , work_v              & ! intent(out)
                           , npolys_run          & ! intent(out)
                           , ed_alloc_work_vec   & ! subroutine
                           , ed_nullify_work_vec ! ! subroutine
   use soil_coms    , only : nslcon              & ! intent(in)
                           , ed_nstyp            ! ! intent(in)
   use mem_polygons , only : maxsite             ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ifm
   integer, intent(in) :: nxp
   integer, intent(in) :: nyp
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: poly
   integer             :: i
   integer             :: j
   integer             :: itext
   !---------------------------------------------------------------------------------------!



   !----- Compute total work load over each row and over entire domain. -------------------!
   npolys_run(ifm) = 0
   do j = 1,nyp
      do i = 1,nxp
         if(work_e(ifm)%land(i,j)) then
            npolys_run(ifm) = npolys_run(ifm) + 1
         end if
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !----- Allocate the polygon vectors. ---------------------------------------------------!
   call ed_nullify_work_vec(work_v(ifm))
   call ed_alloc_work_vec(work_v(ifm),npolys_run(ifm),maxsite)
   !---------------------------------------------------------------------------------------!




   !----- Copy variables to the vector version of the work structure. ---------------------!
   poly = 0
   do j = 1,nyp
      do i = 1,nxp
         
         if(work_e(ifm)%land(i,j)) then

            poly = poly + 1

            work_v(ifm)%glon    (poly) = work_e(ifm)%glon(i,j)
            work_v(ifm)%glat    (poly) = work_e(ifm)%glat(i,j)
            work_v(ifm)%landfrac(poly) = work_e(ifm)%landfrac(i,j)
            work_v(ifm)%work    (poly) = work_e(ifm)%work(i,j)
            work_v(ifm)%xid     (poly) = i
            work_v(ifm)%yid     (poly) = j

            do itext=1,maxsite
               work_v(ifm)%ntext   (itext,poly) = work_e(ifm)%ntext   (itext,i,j)
               work_v(ifm)%soilfrac(itext,poly) = work_e(ifm)%soilfrac(itext,i,j)
            end do
         end if
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_parvec_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will retrieve the workload information stored in the history file,   !
! and use this information to fill the local work array.                                   !
!------------------------------------------------------------------------------------------!
subroutine ed_load_work_from_history()
   use ed_max_dims , only : str_len             & ! intent(in)
                          , maxfiles            & ! intent(in)
                          , maxlist             & ! intent(in)
                          , huge_polygon        ! ! intent(in)
   use ed_misc_coms, only : runtype             & ! intent(in)
                          , ied_init_mode       & ! intent(in)
                          , sfilin              & ! intent(in)
                          , current_time        ! ! intent(in)
   use grid_coms   , only : ngrids              ! ! intent(in)
   use ed_work_vars, only : work_v              & ! intent(inout)
                          , ed_alloc_work_vec   & ! subroutine
                          , ed_nullify_work_vec & ! subroutine
                          , npolys_run          ! ! intent(out)
   use mem_polygons, only : n_ed_region         ! ! intent(in)
   use hdf5_coms   , only : file_id             & ! intent(inout)
                          , dset_id             & ! intent(inout)
                          , dspace_id           & ! intent(inout)
                          , plist_id            & ! intent(inout)
                          , globdims            & ! intent(inout)
                          , chnkdims            & ! intent(inout)
                          , chnkoffs            & ! intent(inout)
                          , cnt                 & ! intent(inout)
                          , stride              & ! intent(inout)
                          , memdims             & ! intent(inout)
                          , memoffs             & ! intent(inout)
                          , memsize             & ! intent(inout)
                          , datatype_id         ! ! intent(inout)
#if USE_HDF5
   use hdf5
#endif
   implicit none
#if USE_HDF5
   !----- Local variables. ----------------------------------------------------------------!
   character(len=str_len), dimension(maxlist)               :: full_list
   character(len=str_len), dimension(maxfiles)              :: histo_list
   character(len=str_len)                                   :: hnamel
   character(len=3)                                         :: cgr
   integer, dimension(:)                      , allocatable :: pclosest
   integer, dimension(:)                      , allocatable :: psrcfile
   integer                                                  :: hdferr
   integer                                                  :: npolys_histo
   integer                                                  :: nflist
   integer                                                  :: nhisto
   integer                                                  :: nf
   integer                                                  :: ifm
   integer                                                  :: ipr
   integer                                                  :: iph
   integer                                                  :: inn
   integer                                                  :: dsetrank
   integer                                                  :: ipya
   integer                                                  :: ipyz
   logical                                                  :: exists
   logical                                                  :: success
   real(kind=8)                                             :: dbletime
   real                                                     :: maxwork
   real, dimension(13,huge_polygon)                         :: worklastyear
   real, dimension(huge_polygon)                            :: histowork
   real, dimension(huge_polygon)                            :: histolon
   real, dimension(huge_polygon)                            :: histolat
   real, dimension(huge_polygon)                            :: polydist
   !----- Local constants. ----------------------------------------------------------------!
   character(len=1), parameter   :: vnam        = 'S'
   integer         , parameter   :: iparallel   = 0
   real            , parameter   :: one_twelfth = 1./12.
   !----- External functions. -------------------------------------------------------------!
   real            , external    :: dist_gc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Here we decide whether this is a history or an ED-2.1 restart run.  In case none  !
   ! of them are true, we simply do nothing, and assume homogeneous workload.              !
   !---------------------------------------------------------------------------------------!
   if (ied_init_mode /= 4 .and. ied_init_mode /= 5 .and. trim(runtype) /= 'HISTORY') then
      return
   end if

   write (unit=*,fmt='(a)') '-------------------------------------------------------------'
   write (unit=*,fmt='(a)') '    === Reading past work from history (restart) file ===    '
   write (unit=*,fmt='(a)') '-------------------------------------------------------------'


   !---------------------------------------------------------------------------------------!
   !     Open the HDF5 environment and                                                     !
   !---------------------------------------------------------------------------------------!
   call h5open_f(hdferr)

   !---------------------------------------------------------------------------------------!
   !    Loop over all regions.  This is useful only for regional runs, so we skip the      !
   ! POI grids.                                                                            !
   !---------------------------------------------------------------------------------------!
   gridloop: do ifm = 1, n_ed_region

      !------------------------------------------------------------------------------------!
      !     Initialise the file list.                                                      !
      !------------------------------------------------------------------------------------!
      full_list(:)  = ''
      histo_list(:) = '' 
      write(cgr,fmt='(a1,i2.2)') 'g',ifm
      
      if (trim(runtype) == 'HISTORY') then
         !----- Full history restart. -----------------------------------------------------!
         dbletime=dble(current_time%time)     
         call makefnam(hnamel,sfilin(1),dbletime,current_time%year,current_time%month      &
                      ,current_time%date,0,vnam,cgr,'h5 ')
         full_list (1) = hnamel
         histo_list(1) = hnamel
         nflist        = 1
         nhisto        = 1

         !----- Checking whether the only history (restart) file exists. ------------------!
         inquire(file=trim(hnamel),exist=exists)
         if (.not.exists) then
            call fatal_error ('File '//trim(hnamel)//' not found.'                         &
                             ,'ed_load_work_from_history','ed_para_init.f90')
         end if

      elseif (ied_init_mode == 4) then
         !----- Standard restart. ---------------------------------------------------------!
         hnamel = trim(sfilin(ifm))//"-"//trim(cgr)//".h5"
         full_list (1) = hnamel
         histo_list(1) = hnamel
         nflist        = 1
         nhisto        = 1

         !----- Checking whether the only history (restart) file exists. ------------------!
         inquire(file=trim(hnamel),exist=exists)
         if (.not.exists) then
            call fatal_error ('File '//trim(hnamel)//' not found.'                         &
                             ,'ed_load_work_from_history','ed_para_init.f90')
         end if

      elseif (ied_init_mode == 5) then
         !----- Unstructured restart. -----------------------------------------------------!

         !----- Retrieve all files with the specified prefix. -----------------------------!
         call ed_filelist(full_list,sfilin(ifm),nflist)
         !----- Check every file and save only those that are actually history files. -----!
         call ed21_fileinfo(nflist,full_list,nhisto,histo_list)
      end if



      !------------------------------------------------------------------------------------!
      !     Loop over all files.  First we initialise some variables.                      !
      !------------------------------------------------------------------------------------!
      success = .true.
      ipyz    = 0
      h5loop: do nf=1,nhisto
         hnamel = histo_list(nf)

         write (unit=*,fmt='(a,1x,2a)') ' - Reading ',trim(hnamel),'...'
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr /= 0) then
            write(unit=*,fmt='(a,1x,i8)') 'Error opening HDF5 file - error - ',hdferr
            call fatal_error('Error opening HDF5 file '//trim(hnamel)//'...'               &
                            ,'read_ed21_history_fill','ed_history_io.f90')
         end if


         !----- Initialise the dimensional control variables for the H5 slabs. ------------!
         globdims = 0_8
         chnkdims = 0_8
         chnkoffs = 0_8
         memoffs  = 0_8
         memdims  = 0_8
         memsize  = 1_8  

         !---------------------------------------------------------------------------------!
         !     Retrieve the global number of polygons.  Inside the file loop, npolys_histo !
         ! is the number of polygons in this file.  Its definition will change outside the !
         ! loop...                                                                         !
         !---------------------------------------------------------------------------------!
         globdims(1) = 1_8

         call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,npolys_histo,globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
 
         !----- Determine the first and last polygon "global" position. -------------------!
         ipya = ipyz + 1
         ipyz = ipyz + npolys_histo


         !----- Load 1D dataset: longitude and latitude. ----------------------------------!
         dsetrank     = 1
         globdims (1) = int(npolys_histo,8)
         chnkdims (1) = int(npolys_histo,8)
         memdims  (1) = int(npolys_histo,8)
         memsize  (1) = int(npolys_histo,8)
         chnkoffs (1) = 0_8
         memoffs  (1) = 0_8

         call hdf_getslab_r(histolon(ipya:ipyz),'LONGITUDE ',dsetrank,iparallel,.true.)
         call hdf_getslab_r(histolat(ipya:ipyz),'LATITUDE ' ,dsetrank,iparallel,.true.)


         !----- Load the workload into a temporary array. ---------------------------------!
         dsetrank     = 2
         globdims (1) = int(13,8)
         chnkdims (1) = int(13,8)
         memdims  (1) = int(13,8)
         memsize  (1) = int(13,8)
         chnkoffs (1) = 0_8
         memoffs  (1) = 0_8
         globdims (2) = int(npolys_histo,8)
         chnkdims (2) = int(npolys_histo,8)
         memdims  (2) = int(npolys_histo,8)
         memsize  (2) = int(npolys_histo,8)
         chnkoffs (2) = 0_8
         memoffs  (2) = 0_8

         call hdf_getslab_r(worklastyear(:,ipya:ipyz),'WORKLOAD '                          &
                           ,dsetrank,iparallel,.false.)

         !----- Here we close the HDF5 file. ----------------------------------------------!
         call h5fclose_f(file_id, hdferr)
         if (hdferr /= 0) then
            write (unit=*,fmt='(a)'      ) 'Could not close the HDF file...'
            write (unit=*,fmt='(a,1x,i6)') 'HDF error=',hdferr
            call fatal_error('Could not close the HDF file','ed_load_work_from_history'       &
                            ,'ed_para_init.f90')
         end if

         !---------------------------------------------------------------------------------!
         !    Workload is a recent variable and may not be available in a history file     !
         ! generated by an earlier version.  Here we check whether anything was read.  In  !
         ! case not, we assume that all polygons have the same workload.                   !
         !---------------------------------------------------------------------------------!
         if (all(worklastyear(:,ipya:ipyz) == 0.)) then
            write (unit=*,fmt='(a)') ' - Failed reading the workload :( ...'
            !------------------------------------------------------------------------------! 
            !     No workload was found, assign ones to all nodes, even if some files      !
            ! had the information.  Also, we leave the loop and use homogeneous workload   !
            !------------------------------------------------------------------------------!
            histowork(:) = 1.
            success      = .false.
            exit h5loop
         else
            write (unit=*,fmt='(a)') ' - Workload was successfully read :) ...'
            !----- And here we compute the average workload during the last year. ---------!
            foundloop: do iph=ipya,ipyz
               histowork(iph) = one_twelfth * sum(worklastyear(1:12,iph))
            end do foundloop
         end if
      end do h5loop
      !------------------------------------------------------------------------------------!
      !     Here we check whether the workload reading was successful or not.  Our action  !
      ! is going to be different depending on what we provide.                             !
      !------------------------------------------------------------------------------------!
      if (success) then
         !---------------------------------------------------------------------------------!
         !     History work was successfully loaded, now we match the workload with the    !
         ! actual polygons by using the nearest neighbour.  We then divide by the maximum  !
         ! workload that is actually used in the simulation to find the relative workload. !
         !---------------------------------------------------------------------------------!
         !----- Change the definition of npolys_histo to the global number of polygons. ---!
         npolys_histo = ipyz
         maxwork = 0.
         polydist = 1.e20
         runpolyloopyes: do ipr= 1,npolys_run(ifm)
            histopolyloop: do iph = 1, npolys_histo
               polydist(iph) = dist_gc(work_v(ifm)%glon(ipr),histolon(iph)                 &
                                      ,work_v(ifm)%glat(ipr),histolat(iph))
            end do histopolyloop

            !----- The closest polygon is the one with minimum distance... ----------------!
            iph                   = minloc(polydist,dim=1)
            work_v(ifm)%work(ipr) = histowork(iph)
            maxwork               = max(maxwork,work_v(ifm)%work(ipr))

         end do runpolyloopyes
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Fiasco.  Since we didn't properly read the workload, use homogeneous        !
         ! instead.                                                                        !
         !---------------------------------------------------------------------------------!
         runpolyloopno: do ipr= 1,npolys_run(ifm)
            work_v(ifm)%work(ipr) = 1.
         end do runpolyloopno
         maxwork = 1.
      end if



      !------------------------------------------------------------------------------------!
      !    Find the maximum workload, and scale all polygons by this number.   Then we     !
      ! switch the relative workload by the cumulative sum.                                !
      !------------------------------------------------------------------------------------!
      if (maxwork <= 0) then
         write (unit=*,fmt='(a,1x,i6)'    ) ' + Grid:             ',ifm
         write (unit=*,fmt='(a,1x,es12.5)') ' + Maximum workload: ',maxwork
         call fatal_error('The maximum workload can''t be right!!!'                        &
                         ,'ed_load_work_from_history','ed_para_init.f90')
      else
         normpolyloop: do ipr= 1, npolys_run(ifm)
            work_v(ifm)%work(ipr) = work_v(ifm)%work(ipr) / maxwork
         end do normpolyloop
      end if
      !------------------------------------------------------------------------------------!
   end do gridloop

   !----- Close the HDF environment. ------------------------------------------------------!
   call h5close_f(hdferr)
#else
   call fatal_error('ED2 now requires HDF5...','ed_load_work_from_history'                 &
                   ,'ed_para_init.F90')
#endif

   write (unit=*,fmt='(a)') '-------------------------------------------------------------'
   write (unit=*,fmt='(a)') ' '

   return
end subroutine ed_load_work_from_history
!==========================================================================================!
!==========================================================================================!
