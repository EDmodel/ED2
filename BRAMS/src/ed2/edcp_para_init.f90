!==========================================================================================!
!==========================================================================================!
!    Subroutines based on the RAMS node decomposition. The main difference between the     !
! original code and this one is that when we split the domain we need to consider whether  !
! the polygon will fall on land or water. The water ones will be removed, so this should   !
! be taken into account for the standalone version.                                        !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will determine the workload related to ED for each of the grid      !
! points.  A polygon will be assigned to the grid cell whenever the cell has at least 1%   !
! of land in the area.  This number could go to the name list at some point...             !
!------------------------------------------------------------------------------------------!
subroutine edcp_get_work(ifm,nxp,nyp,inode,mxp,myp,ia,iz,i0,ja,jz,j0)

   use ed_work_vars, only : work_e                 ! ! intent(inout)
   use soil_coms   , only : veg_database           & ! intent(in)
                          , soil_database          & ! intent(in)
                          , nslcon                 & ! intent(in)
                          , isoilcol               ! ! intent(in)
   use io_params   , only : b_isoilflg => isoilflg & ! intent(in)
                          , b_ivegtflg => ivegtflg ! ! intent(in)
   use mem_polygons, only : n_poi                  & ! intent(in)
                          , poi_res                & ! intent(in)
                          , maxsite                ! ! intent(in)
   use ed_misc_coms, only : min_site_area          ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in)  :: ifm
   integer                , intent(in)  :: nxp
   integer                , intent(in)  :: nyp
   integer                , intent(in)  :: inode
   integer                , intent(in)  :: mxp
   integer                , intent(in)  :: myp
   integer                , intent(in)  :: ia
   integer                , intent(in)  :: iz
   integer                , intent(in)  :: i0
   integer                , intent(in)  :: ja
   integer                , intent(in)  :: jz
   integer                , intent(in)  :: j0
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(:,:), allocatable :: lat_list
   real   , dimension(:,:), allocatable :: lon_list
   integer, dimension(:,:), allocatable :: leaf_class_list
   integer, dimension(:,:), allocatable :: ntext_soil_list
   integer, dimension(:,:), allocatable :: ncol_soil_list
   real   , dimension(:,:), allocatable :: ipcent_land
   real   , dimension(:,:), allocatable :: ipcent_soil
   integer                              :: npoly
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


   !----- This is the total number of potential polygons. ---------------------------------!
   npoly = (iz-ia+1) * (jz-ja+1)


   !----- Allocate some scratch arrays. ---------------------------------------------------!
   allocate(lat_list       (      3,npoly))
   allocate(lon_list       (      3,npoly))
   allocate(leaf_class_list(maxsite,npoly))
   allocate(ntext_soil_list(maxsite,npoly))
   allocate(ncol_soil_list (maxsite,npoly))
   allocate(ipcent_land    (maxsite,npoly))
   allocate(ipcent_soil    (maxsite,npoly))
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !  Fill in the longitude and latitude lists.                                            !
   !                                                                                       !
   !  i index: West-East (longitude).  It increases eastwards.                             !
   !  j index: South-North (latitude).  It increases northwards.                           !
   !---------------------------------------------------------------------------------------!
   if (n_poi > 0 .and. ifm <= n_poi) then
      ipy = 0
      latmap1: do j = ja,jz
         lonmap1: do i =ia,iz

            !----- Update the polygon count. ----------------------------------------------! 
            ipy = ipy + 1

            !----- Co-ordinates at grid centre. -------------------------------------------!
            lon_list(1,ipy) = work_e(ifm)%glon(i,j)
            lat_list(1,ipy) = work_e(ifm)%glat(i,j)

            !----- Co-ordinates at grid North-eastern corner. -----------------------------!
            lon_list(2,ipy) = work_e(ifm)%glon(i,j) - 0.5 * poi_res(ifm)
            lat_list(2,ipy) = work_e(ifm)%glat(i,j) + 0.5 * poi_res(ifm)

            !----- Co-ordinates at grid South-western corner. -----------------------------!
            lon_list(3,ipy) = work_e(ifm)%glon(i,j) + 0.5 * poi_res(ifm)
            lat_list(3,ipy) = work_e(ifm)%glat(i,j) - 0.5 * poi_res(ifm)

            !----- Make sure all co-ordinates are in the expected range. ------------------!
            if (lon_list(1,ipy) >=  180.) lon_list(1,ipy) = lon_list(1,ipy) - 360.
            if (lon_list(1,ipy) <= -180.) lon_list(1,ipy) = lon_list(1,ipy) + 360.
            if (lon_list(2,ipy) >=  180.) lon_list(2,ipy) = lon_list(2,ipy) - 360.
            if (lon_list(2,ipy) <= -180.) lon_list(2,ipy) = lon_list(2,ipy) + 360.
            if (lon_list(3,ipy) >=  180.) lon_list(3,ipy) = lon_list(3,ipy) - 360.
            if (lon_list(3,ipy) <= -180.) lon_list(3,ipy) = lon_list(3,ipy) + 360.

         end do lonmap1
      end do latmap1

   else
      ipy = 0
      latmap2: do j = ja,jz
         lonmap2: do i = ia,iz

            !----- Update the polygon count. ----------------------------------------------!
            ipy = ipy + 1

            !----- Define the standard longitudinal off-sets for neighbour polygons. ------!
            if (i == 1) then
               !---- Western edge. --------------------------------------------------------!
               iloff = -1
               iroff =  1
            elseif (i == mxp) then
               !---- Eastern edge. --------------------------------------------------------!
               iloff =  1
               iroff = -1
            else
               !---- Middle of the domain. ------------------------------------------------!
               iloff =  1
               iroff =  1
            end if
            !------------------------------------------------------------------------------!


            !----- Define the standard latitudinal off-sets for neighbour polygons. -------!
            if (j == 1) then
               !---- Southern edge. -------------------------------------------------------!
               jboff = -1
               jtoff =  1
            elseif (j == myp) then
               !---- Northern edge. -------------------------------------------------------!
               jboff =  1
               jtoff = -1
            else
               !---- Middle of the domain. ------------------------------------------------!
               jboff =  1
               jtoff =  1
            end if
            !------------------------------------------------------------------------------!

            !----- Co-ordinates at grid centre. -------------------------------------------!
            lon_list(1,ipy) = work_e(ifm)%glon(i,j)
            lat_list(1,ipy) = work_e(ifm)%glat(i,j)

            !----- Co-ordinates at grid North-eastern corner. -----------------------------!
            lon_list(2,ipy) = work_e(ifm)%glon(i,j)                                        &
                            + real(iloff) * 0.5 * ( work_e(ifm)%glon(i-iloff,j      )      &
                                                  - work_e(ifm)%glon(i      ,j      ) )
            lat_list(2,ipy) = work_e(ifm)%glat(i,j)                                        &
                            + real(jtoff) * 0.5 * ( work_e(ifm)%glat(i      ,j+jtoff)      &
                                                  - work_e(ifm)%glat(i      ,j      ) )

            !----- Co-ordinates at grid South-western corner. -----------------------------!
            lon_list(3,ipy) = work_e(ifm)%glon(i,j)                                        &
                            + real(iroff) * 0.5 * ( work_e(ifm)%glon(i+iroff,j      )      &
                                                  - work_e(ifm)%glon(i      ,j      ) )
            lat_list(3,ipy) = work_e(ifm)%glat(i,j)                                        &
                            + real(jboff) * 0.5 * ( work_e(ifm)%glat(i      ,j-jboff)      &
                                                  - work_e(ifm)%glat(i      ,j      ) )


            !----- Make sure all co-ordinates are in the expected range. ------------------!
            if (lon_list(1,ipy) >=  180.) lon_list(1,ipy) = lon_list(1,ipy) - 360.
            if (lon_list(1,ipy) <= -180.) lon_list(1,ipy) = lon_list(1,ipy) + 360.
            if (lon_list(2,ipy) >=  180.) lon_list(2,ipy) = lon_list(2,ipy) - 360.
            if (lon_list(2,ipy) <= -180.) lon_list(2,ipy) = lon_list(2,ipy) + 360.
            if (lon_list(3,ipy) >=  180.) lon_list(3,ipy) = lon_list(3,ipy) - 360.
            if (lon_list(3,ipy) <= -180.) lon_list(3,ipy) = lon_list(3,ipy) + 360.
            
         end do lonmap2
      end do latmap2
   end if


   !----- Generate the land/sea mask. -----------------------------------------------------!

   write(unit=*,fmt='(a)')       '   + Obtaining land/sea mask information.'
   select case (b_ivegtflg(ifm))
   case (0,1,2)
      call leaf3init_overwrite(mxp,myp,ia,iz,i0,ja,jz,j0,maxsite,npoly,ifm,'leaf'          &
                              ,leaf_class_list,ipcent_land)
   case (3)
      call leaf_database(trim(veg_database(ifm)),maxsite,npoly,'leaf_class'                &
                        ,lat_list,lon_list,leaf_class_list,ipcent_land)
   end select
   write(unit=*,fmt='(a)')       '   + Obtaining soil texture class information.'

   select case (b_isoilflg(ifm))
   !---------------------------------------------------------------------------------------!
   !    This allows us to use the vegetation and soil type defined for LEAF in coupled     !
   ! runs. Note that to use the ED dataset in coupled runs we need to use isoilflg =3.     !
   !---------------------------------------------------------------------------------------!
   case (0)
      !------------------------------------------------------------------------------------!
      !   Allow for only one site by making the first site with the default soil type and  !
      ! area 1., and the others with area 0.                                               !
      !------------------------------------------------------------------------------------!
      ntext_soil_list        (:,:) = nslcon
      ipcent_soil            (:,:) = 0.
      ipcent_soil            (1,:) = 1.
      !------------------------------------------------------------------------------------!
   case (1,2)
      call leaf3init_overwrite(mxp,myp,ia,iz,i0,ja,jz,j0,maxsite,npoly,ifm,'soil'          &
                              ,ntext_soil_list,ipcent_soil)
   case (3)
      call leaf_database(trim(soil_database(ifm)),maxsite,npoly,'soil_text'                &
                        ,lat_list,lon_list,ntext_soil_list,ipcent_soil)
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      For the time being, soil colour is constant.  Only if results look promising we  !
   ! will attempt to read a map.                                                           !
   !---------------------------------------------------------------------------------------!
   ncol_soil_list (:,:) = isoilcol
   !---------------------------------------------------------------------------------------!

   write(unit=*,fmt='(a)')       '   + Successfully obtain land/sea mask and soil type.'

   !----- Re-map the land cover classes. --------------------------------------------------!
   ipy     = 0
   maxwork = epsilon(0.)
   latwork: do j = ja,jz
      lonwork: do i = ia,iz

         !---- Update the polygon index. --------------------------------------------------!
         ipy = ipy + 1

         work_e(ifm)%land(i,j)      = ipcent_land(1,ipy) > min_site_area

         if (work_e(ifm)%land(i,j)) then

            work_e(ifm)%landfrac(i,j)  = ipcent_land(1,ipy)

            work_e(ifm)%work(i,j) = 0.0
            do itext = 1,maxsite
               if (ipcent_soil(itext,ipy) > min_site_area) then
                  work_e(ifm)%work(i,j) = work_e(ifm)%work(i,j) + 1.
               end if
               work_e(ifm)%soilfrac(itext,i,j) = ipcent_soil     (itext,ipy)
               work_e(ifm)%ntext   (itext,i,j) = ntext_soil_list (itext,ipy)
            end do
            work_e(ifm)%nscol   (i,j) = ncol_soil_list  (1,ipy)
            maxwork = max(maxwork,work_e(ifm)%work(i,j))
         else
            !----- Making this grid point 100% water --------------------------------------!
            work_e(ifm)%landfrac  (i,j)  = 0.
            work_e(ifm)%work      (i,j)  = epsilon(0.0)
            work_e(ifm)%nscol     (i,j)  = 0
            work_e(ifm)%ntext   (:,i,j)  = 0
            work_e(ifm)%soilfrac(:,i,j)  = 0
         end if

      end do lonwork
   end do latwork
   !---------------------------------------------------------------------------------------!




   !----- Normalise the workload. ---------------------------------------------------------!
   do j=ja,jz
      do i=ia,iz
         work_e(ifm)%work(i,j) = work_e(ifm)%work(i,j) / maxwork
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !----- De-allocate all temporary arrays. -----------------------------------------------!
   deallocate(lat_list        )
   deallocate(lon_list        )
   deallocate(leaf_class_list )
   deallocate(ntext_soil_list )
   deallocate(ncol_soil_list  )
   deallocate(ipcent_land     )
   deallocate(ipcent_soil     )
   !---------------------------------------------------------------------------------------!

   return
end subroutine edcp_get_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will copy the information from the (lon,lat) matrix-shaped work     !
! structure to the vector form.  Here we will also update the patch area in the LEAF-3     !
! structure, to match with ED polygons.  Notice that now the work_e structure has the same !
! dimension as leaf, and the local edges will have the values updated by message-passing   !
! routines from BRAMS.                                                                     !
!------------------------------------------------------------------------------------------!
subroutine edcp_parvec_work(ifm,nxp,nyp,ia,iz,i0,ja,jz,j0)
   use ed_work_vars, only : work_e              & ! intent(in)
                          , work_v              & ! intent(inout)
                          , ed_nullify_work_vec & ! subroutine
                          , ed_alloc_work_vec   ! ! subroutine
   use mem_leaf    , only : leaf_g              ! ! intent(inout)
   use ed_node_coms, only : mynum               ! ! intent(in)
   use mem_polygons, only : maxsite             ! ! intent(in)
   use soil_coms   , only : isoilcol            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                   , intent(in) :: ifm
   integer                   , intent(in) :: nxp
   integer                   , intent(in) :: nyp
   integer                   , intent(in) :: ia
   integer                   , intent(in) :: iz
   integer                   , intent(in) :: i0
   integer                   , intent(in) :: ja
   integer                   , intent(in) :: jz
   integer                   , intent(in) :: j0
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: npolygons
   integer                                :: poly
   integer                                :: i
   integer                                :: j
   integer                                :: itext
   !---------------------------------------------------------------------------------------!

   !----- Compute total work load over each row and over entire domain. -------------------!
   npolygons = count(work_e(ifm)%land)

   !----- Allocate the polygon vectors. ---------------------------------------------------!
   call ed_nullify_work_vec(work_v(ifm))
   call ed_alloc_work_vec(work_v(ifm),npolygons,maxsite)

   !----- Loop only over the polygons that will be actually run in this node. -------------!
   poly = 0
   latloop: do j = ja,jz
      lonloop: do i = ia,iz

         if (work_e(ifm)%land(i,j)) then
            poly = poly + 1
            work_v(ifm)%glon       (poly) = work_e(ifm)%glon    (i,j)
            work_v(ifm)%glat       (poly) = work_e(ifm)%glat    (i,j)
            work_v(ifm)%landfrac   (poly) = work_e(ifm)%landfrac(i,j)
            work_v(ifm)%xid        (poly) = work_e(ifm)%xatm    (i,j)
            work_v(ifm)%yid        (poly) = work_e(ifm)%yatm    (i,j)
            work_v(ifm)%nscol      (poly) = work_e(ifm)%nscol   (i,j)

            do itext=1,maxsite
               work_v(ifm)%ntext   (itext,poly) = work_e(ifm)%ntext   (itext,i,j)
               work_v(ifm)%soilfrac(itext,poly) = work_e(ifm)%soilfrac(itext,i,j)
            end do


            !----- Update the patches area so both LEAF-3 and ED-2 structures agree. ------!
            leaf_g(ifm)%patch_area(i,j,1) = 1.0 - work_e(ifm)%landfrac(i,j)
            do itext=1,maxsite
               leaf_g(ifm)%patch_area(i,j,1+itext) = work_e(ifm)%landfrac(i,j)             &
                                                   * work_e(ifm)%soilfrac(itext,i,j)
            end do
         else
            !----- Update the patches area so both LEAF-3 and ED-2 structures agree. ------!
            leaf_g(ifm)%patch_area(i,j,1) = 1.0
            do itext=1,maxsite
               leaf_g(ifm)%patch_area(i,j,1+itext) = 0.0
            end do
         end if
      end do lonloop
   end do latloop

   !----- Just a sanity check. ------------------------------------------------------------!
   if (npolygons /= poly) then
      write (unit=*,fmt='(a)')       '----------------------------------------------------'
      write (unit=*,fmt='(a,1x,i7)') ' - MYNUM     =',mynum
      write (unit=*,fmt='(a,1x,i7)') ' - POLY      =',poly
      write (unit=*,fmt='(a,1x,i7)') ' - NPOLYGONS =',npolygons
      write (unit=*,fmt='(a)')       '----------------------------------------------------'
      call abort_run  ('POLY and NPOLYGONS should be the same, but they are not...'        &
                      ,'edcp_parvec_work','edcp_para_init.f90')
   end if

   return
end subroutine edcp_parvec_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine will copy LEAF-3 initial condition to ED, so we can use other LEAF     !
! databases to decide whether a polygon is inland or offshore, and aslo other soil texture !
! dataset.                                                                                 !
!------------------------------------------------------------------------------------------!
subroutine leaf3init_overwrite(nxp,nyp,ia,iz,i0,ja,jz,j0,nsite,nlandsea,ifm,varname        &
                              ,classout,pctout)
   use mem_leaf    , only : leaf_g ! ! intent(in)
   use mem_grid    , only : nzg    ! ! intent(in)
   use node_mod    , only : mynum  ! ! intent(in)
   use mem_polygons, only : maxsite ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                           , intent(in)  :: nxp
   integer                           , intent(in)  :: nyp
   integer                           , intent(in)  :: ia
   integer                           , intent(in)  :: iz
   integer                           , intent(in)  :: i0
   integer                           , intent(in)  :: ja
   integer                           , intent(in)  :: jz
   integer                           , intent(in)  :: j0
   integer                           , intent(in)  :: nsite
   integer                           , intent(in)  :: nlandsea
   integer                           , intent(in)  :: ifm
   character(len=4)                  , intent(in)  :: varname
   integer, dimension(nsite,nlandsea), intent(out) :: classout
   real   , dimension(nsite,nlandsea), intent(out) :: pctout
   !----- Local variables -----------------------------------------------------------------!
   integer                                         :: i
   integer                                         :: j
   integer                                         :: iglob
   integer                                         :: jglob
   integer                                         :: ipy
   integer                                         :: itext
   integer                                         :: ilpat
   real                                            :: land_frac
   !---------------------------------------------------------------------------------------!

   select case (varname)
   case ('leaf')
      ipy=0

      do j=ja,jz
         jglob = j + j0
         do i = ia,iz
            iglob = i + i0

            ipy         = ipy + 1
            !------------------------------------------------------------------------------!
            !     Classout shouldn't be used, here we make the first the sum of all land   !
            ! patches in LEAF and attribute.  We don't use any other value other than the  !
            ! first.                                                                       !
            !------------------------------------------------------------------------------!
            classout      (:,ipy) = nint(leaf_g(ifm)%leaf_class(iglob,jglob,2))
            pctout        (:,ipy) = 0.
            pctout        (1,ipy) = 1. - leaf_g(ifm)%patch_area(iglob,jglob,1)
            !------------------------------------------------------------------------------!
         end do
      end do
   case ('soil')
      ipy=0
      do j=ja,jz
         jglob = j + j0
         do i=ia,iz
            iglob = i + i0

            ipy         = ipy + 1
            !------------------------------------------------------------------------------!
            !    Find the water fraction.  Get the information about this polygon only if  !
            ! it is not 100% water.                                                        !
            !------------------------------------------------------------------------------!
            land_frac = 1. - leaf_g(ifm)%patch_area(iglob,jglob,1)
            if (land_frac > 0.) then
               do itext = 1,maxsite
                  ilpat = itext + 1
                  classout(itext,ipy) = nint(leaf_g(ifm)%soil_text(nzg,iglob,jglob,ilpat))
                  pctout  (itext,ipy) = leaf_g(ifm)%patch_area(iglob,jglob,ilpat)          &
                                      / land_frac
               end do
            else
               !----- Assign a dummy texture, it shouldn't be used. -----------------------!
               classout(:,ipy) = nint(leaf_g(ifm)%soil_text(nzg,iglob,jglob,2))
               pctout  (:,ipy) = 0.
               pctout  (1,ipy) = 1.
            end if
            !------------------------------------------------------------------------------!
         end do
      end do
   case default
      call abort_run('Invalid key: '//varname,'leaf3init_overwrite.f90','edcp_met_init.f90')
   end select
   return
end subroutine leaf3init_overwrite
!==========================================================================================!
!==========================================================================================!
