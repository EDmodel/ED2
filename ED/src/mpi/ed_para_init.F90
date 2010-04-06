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
   use mem_sites   , only : n_ed_region       ! ! intent(in)
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
      !     SOI grids always have one point only. Since the structure will be sent.  I am  !
      ! filling the structures.  Only 4 extra numbers will be sent through MPI if we do    !
      ! this.                                                                              !
      !------------------------------------------------------------------------------------!
      mmxp(ngr) = nnxp(ngr)
      mmyp(ngr) = nnyp(ngr)

      call ed_nullify_work(work_e(ngr))
      call ed_alloc_work(work_e(ngr),nnxp(ngr),nnyp(ngr))
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
   !     Polygons of interest (SOI/POI) are single polygons, so they get the simplest      !
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
   use mem_sites   , only : grid_type      & ! intent(in)
                          , grid_res       & ! intent(in)
                          , n_ed_region    & ! intent(in)
                          , soi_lat        & ! intent(in)
                          , soi_lon        & ! intent(in)
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
            work_e(ifm)%glon(i,j)=soi_lon(ifm-n_ed_region)
            work_e(ifm)%glat(i,j)=soi_lat(ifm-n_ed_region)
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
   use mem_sites   , only : n_soi          & ! intent(in)
                          , grid_res       & ! intent(in)
                          , grid_type      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ifm
   integer, intent(in) :: nxp
   integer, intent(in) :: nyp
   !----- Local variables. ----------------------------------------------------------------!
   integer :: npoly
   real   , dimension(:,:), allocatable :: lat_list
   real   , dimension(:,:), allocatable :: lon_list
   integer, dimension(:)  , allocatable :: leaf_class_list
   integer, dimension(:)  , allocatable :: ntext_soil_list
   integer, dimension(:)  , allocatable :: ipcent_land
   integer                              :: datsoil
   integer                              :: ipy
   integer                              :: i
   integer                              :: j
   integer                              :: jboff
   integer                              :: jtoff
   integer                              :: iloff
   integer                              :: iroff
   !----- Local constants. ----------------------------------------------------------------!
   integer                , parameter   :: min_land_pcent = 25 ! Mininum percentage of land
                                                               !   that a polygon must have
                                                               !   to be considered in the
                                                               !   regional run.
   real                   , parameter   :: soi_edge_deg = 0.05 ! 100th of a degree, about 
                                                               !   5.5 km at the Equator.
   !---------------------------------------------------------------------------------------!

   npoly = nxp*nyp
   allocate(lat_list(3,npoly))
   allocate(lon_list(3,npoly))
   allocate(leaf_class_list(npoly))
   allocate(ipcent_land(npoly))

   !---------------------------------------------------------------------------------------!
   !     Fill lat/lon lists.  The longitude and latitude are initially assigned as arrays. !
   ! The second index (the npoly one) is used to define each polygon.  The first index is  !
   ! used to define the relative coordinate of important points within that polygon:       !
   ! 1 means grid centre, 2 means Northwestern corner, and 3 means Southeastern corner.    !
   !---------------------------------------------------------------------------------------!
   if (n_soi > 0 .and. ifm <= n_soi) then
      ipy = 0
      do i=1,nxp
         do j = 1,nyp
            ipy = ipy + 1

            !----- Grid mid-point. --------------------------------------------------------!
            lon_list(1,ipy) = work_e(ifm)%glon(i,j)
            lat_list(1,ipy) = work_e(ifm)%glat(i,j)
            !----- Northwestern corner. ---------------------------------------------------!
            lon_list(2,ipy) = work_e(ifm)%glon(i,j) - soi_edge_deg
            lat_list(2,ipy) = work_e(ifm)%glat(i,j) + soi_edge_deg
            !----- Southeastern corner. ---------------------------------------------------!
            lon_list(3,ipy) = work_e(ifm)%glon(i,j) + soi_edge_deg
            lat_list(3,ipy) = work_e(ifm)%glat(i,j) - soi_edge_deg

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
      do i=1,nxp
         do j = 1,nyp
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

   call leaf_database(trim(veg_database),npoly,'leaf_class',lat_list,lon_list,ipcent_land)

   if (isoilflg(ifm) == 1) then
      allocate(ntext_soil_list(npoly))
      call leaf_database(trim(soil_database),npoly,'soil_text',lat_list,lon_list           &
                        ,ntext_soil_list)
   end if
 
   !----- Re-map the land cover classes. --------------------------------------------------!
   ipy = 0
   do i=1,nxp
      do j = 1,nyp
         ipy = ipy + 1
         work_e(ifm)%land(i,j) = ipcent_land(ipy) > min_land_pcent

         if (work_e(ifm)%land(i,j)) then
            work_e(ifm)%work(i,j)     = 1.0
            work_e(ifm)%landfrac(i,j) = 0.01 * real(ipcent_land(ipy))

            select case (isoilflg(ifm))
            case (1)  !----- Set from data base or LEAF-3. --------------------------------!
               datsoil = ntext_soil_list(ipy)

               !---------------------------------------------------------------------------!
               !     This is to prevent datsoil to be zero when the polygon was assumed    !
               ! land.                                                                     !
               !---------------------------------------------------------------------------!
               if (datsoil == 0) datsoil=nslcon
               work_e(ifm)%ntext(i,j) = datsoil
            case (2) !! set from ED2IN/RAMSIN
               work_e(ifm)%ntext(i,j) = nslcon
            end select
         else
            !----- Making this grid point 100% water ---------------------------------------!
            work_e(ifm)%landfrac(i,j)  = 0.
            work_e(ifm)%work(i,j)      = epsilon(0.0)
            work_e(ifm)%ntext(i,j)     = 0
         end if
      end do
   end do

  
  
   deallocate(lat_list)
   deallocate(lon_list)
   deallocate(leaf_class_list)
   deallocate(ipcent_land)
   if (allocated(ntext_soil_list)) deallocate (ntext_soil_list)

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
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ifm
   integer, intent(in) :: nxp
   integer, intent(in) :: nyp
   !----- Local variables. ----------------------------------------------------------------!
   integer :: poly
   integer :: i
   integer :: j
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
  
   !----- Allocate the polygon vectors. ---------------------------------------------------!
   call ed_nullify_work_vec(work_v(ifm))
   call ed_alloc_work_vec(work_v(ifm),npolys_run(ifm))

   !----- Copy variables to the vector version of the work structure. ---------------------!
   poly = 0
   do j = 1,nyp
      do i = 1,nxp
         
         if(work_e(ifm)%land(i,j)) then
            poly = poly + 1
            
            work_v(ifm)%glon(poly)     = work_e(ifm)%glon(i,j)
            work_v(ifm)%glat(poly)     = work_e(ifm)%glat(i,j)
            work_v(ifm)%landfrac(poly) = work_e(ifm)%landfrac(i,j)
            work_v(ifm)%work(poly)     = work_e(ifm)%work(i,j)
            if (work_e(ifm)%ntext(i,j) >= 1 .and. work_e(ifm)%ntext(i,j) <= ed_nstyp) then
               work_v(ifm)%ntext(poly)    = work_e(ifm)%ntext(i,j)
            else
               work_v(ifm)%ntext(poly)    = nslcon
            end if
            work_v(ifm)%xid(poly)      = i
            work_v(ifm)%yid(poly)      = j
         end if
      end do
   end do
  
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
   use ed_max_dims , only : str_len             ! ! intent(in)
   use ed_misc_coms, only : runtype             & ! intent(in)
                          , ied_init_mode       & ! intent(in)
                          , sfilin              & ! intent(in)
                          , current_time        ! ! intent(in)
   use grid_coms   , only : ngrids              ! ! intent(in)
   use ed_work_vars, only : work_v              & ! intent(inout)
                          , ed_alloc_work_vec   & ! subroutine
                          , ed_nullify_work_vec & ! subroutine
                          , npolys_run          ! ! intent(out)
   use mem_sites   , only : n_ed_region         ! ! intent(in)
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
   character(len=str_len)              :: hnamel
   character(len=3)                    :: cgr
   integer                             :: hdferr
   integer                             :: npolys_histo
   integer                             :: ifm
   integer                             :: ipr
   integer                             :: iph
   integer                             :: inn
   integer                             :: dsetrank
   logical                             :: exists
   real(kind=8)                        :: dbletime
   real                                :: maxwork
   real, dimension(:,:)  , allocatable :: worklastyear
   real, dimension(:)    , allocatable :: histowork
   real, dimension(:)    , allocatable :: histolon
   real, dimension(:)    , allocatable :: histolat
   real, dimension(:)    , allocatable :: polydist
   !----- Local constants. ----------------------------------------------------------------!
   character(len=1)      , parameter   :: vnam        = 'S'
   real                  , parameter   :: one_twelfth = 1./12.
   integer               , parameter   :: iparallel   = 0
   !----- External functions. -------------------------------------------------------------!
   real                  , external    :: dist_gc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Here we decide whether this is a history or an ED-2.1 restart run.  In case none  !
   ! of them are true, we simply do nothing, and assume homogeneous workload.              !
   !---------------------------------------------------------------------------------------!
   if (ied_init_mode /= 4 .and. trim(runtype) /= 'HISTORY') then
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
   ! SOI/POI grids.                                                                        !
   !---------------------------------------------------------------------------------------!
   gridloop: do ifm = 1, n_ed_region
      write(cgr,fmt='(a1,i2.2)') 'g',ifm
      
      if (trim(runtype) == 'HISTORY') then
         dbletime=dble(current_time%time)     
         call makefnam(hnamel,sfilin,dbletime,current_time%year,current_time%month         &
                      ,current_time%date,0,vnam,cgr,'h5 ')
         
      else !if (ied_init_mode == 4) then
         hnamel = trim(sfilin)//"-"//trim(cgr)//".h5"
      end if

      !----- Opening the history (restart) file. ------------------------------------------!
      inquire(file=trim(hnamel),exist=exists)

      if (.not.exists) then
         call fatal_error ('File '//trim(hnamel)//' not found.'                            &
                          ,'ed_load_work_from_history','ed_para_init.f90')
      else
         write (unit=*,fmt='(a,1x,2a)') ' - Reading ',trim(hnamel),'...'
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr /= 0) then
            write(unit=*,fmt='(a,1x,i8)') 'Error opening HDF5 file - error - ',hdferr
            call fatal_error('Error opening HDF5 file '//trim(hnamel)//'...'               &
                            ,'read_ed21_history_fill','ed_history_io.f90')
         end if
      end if


      !----- Initialise the dimensional control variables for the H5 slabs. ---------------!
      globdims = 0_8
      chnkdims = 0_8
      chnkoffs = 0_8
      memoffs  = 0_8
      memdims  = 0_8
      memsize  = 1_8  

      !------------------------------------------------------------------------------------!
      !     Retrieve the global number of polygons.                                        !
      !------------------------------------------------------------------------------------!
      globdims(1) = 1_8

      call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,npolys_histo,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
     
      !------------------------------------------------------------------------------------!
      !     Now that we know the number of polygons in the history, we allocate the work   !
      ! array, and the scratch that will receive the temporary workload from last year.    !
      !------------------------------------------------------------------------------------!
      allocate(histolon(npolys_histo),histolat(npolys_histo),polydist(npolys_histo))
      allocate(histowork(npolys_histo),worklastyear(13,npolys_histo))

      !----- Load 1D dataset: longitude and latitude. -------------------------------------!
      dsetrank     = 1
      globdims (1) = int(npolys_histo,8)
      chnkdims (1) = int(npolys_histo,8)
      memdims  (1) = int(npolys_histo,8)
      memsize  (1) = int(npolys_histo,8)
      chnkoffs (1) = 0_8
      memoffs  (1) = 0_8

      call hdf_getslab_r(histolon,'LONGITUDE ',dsetrank,iparallel,.true.)
      call hdf_getslab_r(histolat,'LATITUDE ' ,dsetrank,iparallel,.true.)


      !----- Load the workload into a temporary array. ------------------------------------!
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

      call hdf_getslab_r(worklastyear,'WORKLOAD ' ,dsetrank,iparallel,.false.)

      !----- Here we close the HDF5 file. -------------------------------------------------!
      call h5fclose_f(file_id, hdferr)
      if (hdferr /= 0) then
         write (unit=*,fmt='(a)'      ) 'Could not close the HDF file...'
         write (unit=*,fmt='(a,1x,i6)') 'HDF error=',hdferr
         call fatal_error('Could not close the HDF file','ed_load_work_from_history'       &
                         ,'ed_para_init.f90')
      end if

      !------------------------------------------------------------------------------------!
      !    Workload is a recent variable and may not be available in a history file        !
      ! generated by an earlier version.  Here we check whether anything was read.  In     !
      ! case not, we assume that all polygons have the same workload.                      !
      !------------------------------------------------------------------------------------!
      if (all(worklastyear == 0.)) then
         write (unit=*,fmt='(a,1x,2a)') ' - Failed reading the workload :( ...'
         !----- No workload was found, assign ones to all nodes... ------------------------!
         notfoundloop: do iph=1,npolys_histo
            histowork(iph) = 1.
         end do notfoundloop
      else
         write (unit=*,fmt='(a,1x,2a)') ' - Workload was successfully read...'
         !----- And here we compute the average workload during the last year. ------------!
         foundloop: do iph=1,npolys_histo
            histowork(iph) = one_twelfth * sum(worklastyear(1:12,iph))
         end do foundloop
      end if

      !------------------------------------------------------------------------------------!
      !     Now that the history work is loaded, we match the workload with the actual     !
      ! polygons by using the nearest neighbour.  We then divide by the maximum workload   !
      ! that is actually used in the simulation to find the relative workload.             !
      !------------------------------------------------------------------------------------!
      maxwork = 0.
      runpolyloop: do ipr= 1,npolys_run(ifm)
         histopolyloop: do iph = 1, npolys_histo
            polydist(iph) = dist_gc(work_v(ifm)%glon(ipr),histolon(iph)                    &
                                   ,work_v(ifm)%glat(ipr),histolat(iph))
         end do histopolyloop

         !----- The closest polygon is the one with minimum distance... -------------------!
         iph                   = minloc(polydist,dim=1)
         work_v(ifm)%work(ipr) = histowork(iph)
         maxwork               = max(maxwork,work_v(ifm)%work(ipr))
      end do runpolyloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Find the maximum workload, and scale all polygons by this number.   Then we     !
      ! switch the relative workload by the cumulative sum.                                !
      !------------------------------------------------------------------------------------!
      if (maxwork <= 0) then
         write (unit=*,fmt='(a,1x,i6)') ' + Grid:             ',ifm
         write (unit=*,fmt='(a,1x,i6)') ' + Maximum workload: ',maxwork
         call fatal_error('The maximum workload can''t be right!!!'                        &
                         ,'ed_load_work_from_history','ed_para_init.f90')
      else
         normpolyloop: do ipr= 1, npolys_run(ifm)
            work_v(ifm)%work(ipr) = work_v(ifm)%work(ipr) / maxwork
         end do normpolyloop
      end if
      !------------------------------------------------------------------------------------!



      !----- Free the memory associated with the scratch arrays. --------------------------!
      deallocate(histolon,histolat,polydist,histowork,worklastyear)

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
