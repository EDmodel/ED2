!============================= Change Log =================================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This routine decomposes grid domains of size (nnxp,nnyp) into a number, specified by !
! nodes, of rectangular subdomains.  The convention is followed that any internal          !
! boundaries (between subdomains) that are parallel to the x-axis run continuously across  !
! the full domain, while boundaries parallel to the y-axis may or may not run the full     !
! distance across the domain.  For convenience, regions of the domain bounded by adjacent  !
! east-west internal boundaries are termed "slabs", while smaller divisions within each    !
! slab are termed "blocks".  Each block is required to have a minimum dimension of 6 by 6  !
! grid cells.  If this cannot be satisfied with the given input parameters, the subroutine !
! stops.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine par_decomp(ngr,nxp,nyp,work)
   use grid_dims, only : maxmach ! ! intent(in)
   use rpara    , only : ixb     & ! intent(in)
                       , ixe     & ! intent(in)
                       , iyb     & ! intent(in)
                       , iye     & ! intent(in)
                       , nmachs  ! ! intent(in)
   implicit none
   !----- Arguments. -----------------------------------------------------------------------!
   integer                    , intent(in)    :: ngr
   integer                    , intent(in)    :: nxp
   integer                    , intent(in)    :: nyp
   real   , dimension(nxp,nyp), intent(in)    :: work
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(nxp+nyp)                :: workrow
   real   , dimension(nxp+nyp)                :: workload
   real   , dimension(nxp+nyp)                :: workblock
   real   , dimension(nxp+nyp)                :: workcol
   integer, dimension(nxp+nyp)                :: nblocks
   integer, dimension(nxp+nyp)                :: jrows
   integer, dimension(nxp+nyp)                :: jrow
   real   , dimension(maxmach)                :: relspeed
   integer                                    :: inode
   integer                                    :: i
   integer                                    :: j
   integer                                    :: islab
   integer                                    :: jnodes
   integer                                    :: nslabs
   integer                                    :: min_blocks
   integer                                    :: nbigslabs
   integer                                    :: iblock
   integer                                    :: jnode
   integer                                    :: knode
   real                                       :: aslabs
   real                                       :: totspeed
   real                                       :: workdom
   real                                       :: workaccum
   real                                       :: worksofar
   real                                       :: slabspeed
   real                                       :: workslab
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the relative speed for all nodes with default value of 1.0.            !
   !---------------------------------------------------------------------------------------!
   relspeed(:) = 1.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Estimate the number of slabs to be used (aslabs), and compute a final nearest    !
   ! integer value (nslabs) which is limited to allowable values.  Zero out array for      !
   ! accumulating number of columns for each node.                                         !
   !---------------------------------------------------------------------------------------!
   aslabs = sqrt(real(nmachs) * real(nyp) / real(nxp))
   nslabs = min(nmachs,max(1,nint(aslabs)))
   !---------------------------------------------------------------------------------------!

   totspeed = 0.
   do inode = 1,nmachs
      ixe(inode,ngr) = 0
      totspeed = totspeed + relspeed(inode)
   enddo

   !----- Compute total work load over each row and over entire domain. -------------------!
   workdom = 0.
   do j = 1,nyp
      workrow(j) = 0.
      do i = 1,nxp
         workrow(j) = workrow(j) + work(i,j)
      end do
      workdom = workdom + workrow(j)
   end do
   workrow(2)     = workrow(2)     + workrow(1)
   workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)
   !---------------------------------------------------------------------------------------!



   !----- Determine number of blocks and the average workload for each slab. --------------!
   min_blocks = nmachs / nslabs
   nbigslabs  = nmachs - min_blocks * nslabs
   inode = 0
   do islab = 1,nslabs
      workload(islab) = 0.
      nblocks (islab) = min_blocks

      if (islab <= nbigslabs) nblocks(islab) = min_blocks + 1

      do iblock = 1,nblocks(islab)
         inode = inode + 1
         workload(islab) = workload(islab) + workdom * relspeed(inode) / totspeed
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Assign all j-rows to their respective slabs in a way that balances the work load  !
   ! among slabs according to their respective numbers of nodes (blocks).  The array jrows !
   ! counts the number of rows in each slab, and the array jrow is the index of the        !
   ! southernmost row in each slab.                                                        !
   !---------------------------------------------------------------------------------------!
   do islab = 1,nslabs
      jrows(islab) = 0
   end do

   workaccum = 0.
   worksofar = 0.
   islab     = 0

   do j = 2,nyp-1
      workaccum = workaccum + workrow(j)
      if (workaccum - .5 * workrow(j) > worksofar .and. islab < nslabs) then
         islab       = islab + 1
         jrow(islab) = j
         worksofar   = worksofar + workload(islab)
      end if
      jrows(islab) = jrows(islab) + 1
   end do

   inode = 0
   jnode = 0
   knode = 0
   do islab = 1,nslabs
      !------------------------------------------------------------------------------------!
      !      Compute the total work load for each slab and for each i-column in the slab.  !
      !------------------------------------------------------------------------------------!
      slabspeed = 0.
      workslab  = 0.
      do i = 1,nxp
         workcol(i) = 0.

         do j = jrow(islab),jrow(islab)+jrows(islab)-1
            workcol(i) = workcol(i) + work(i,j)
         end do

         workslab = workslab + workcol(i)
      end do
      workcol(2)     = workcol(2) + workcol(1)
      workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)

      !----- Determine average workload for each block. -----------------------------------!
      do iblock = 1,nblocks(islab)
         jnode = jnode + 1
         slabspeed = slabspeed + relspeed(jnode)
      end do

      do iblock = 1,nblocks(islab)
         knode = knode + 1
         workblock(iblock) = workslab * relspeed(knode) / slabspeed
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Assign the i-columns of each slab to their respective blocks in a way that     !
      ! balances the work load among the blocks.  The array ncols counts the number of     !
      ! i-columns on each node, and the array ncol is the index of the westernmost         !
      ! i-column on each node.                                                             !
      !------------------------------------------------------------------------------------!
      workaccum = 0.
      worksofar = 0.

      iblock = 0
      do i = 2,nxp-1
         workaccum = workaccum + workcol(i)

         if (workaccum - .5 * workcol(i) > worksofar .and. iblock < nblocks(islab)) then
            iblock = iblock + 1

            !----- Define the node bound variables here -----------------------------------!
            inode          = inode + 1
            iyb(inode,ngr) = jrow(islab)
            ixb(inode,ngr) = i
            iye(inode,ngr) = iyb(inode,ngr) + jrows(islab) - 1
            !------------------------------------------------------------------------------!

            worksofar = worksofar + workblock(iblock)

         end if
         ixe(inode,ngr) = ixe(inode,ngr) + 1
      end do
   end do

   !----- Defining node variable here. ----------------------------------------------------!
   do jnode = 1,nmachs
      ixe(jnode,ngr) = ixb(jnode,ngr) + ixe(jnode,ngr) - 1
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   ! Check to make sure that each subdomain has at least 2 interior rows and columns.      !
   !---------------------------------------------------------------------------------------!
   do jnode = 1,nmachs
      if (iye(jnode,ngr)-iyb(jnode,ngr) < 1 .or. ixe(jnode,ngr)-ixb(jnode,ngr) < 1) then
         write (unit=*,fmt='(a)')       '-------------------------------------------------'
         write (unit=*,fmt='(a,1x,i6)') ' NXP  = ',nxp
         write (unit=*,fmt='(a,1x,i6)') ' NYP  = ',nyp
         write (unit=*,fmt='(a,1x,i6)') ' Node = ',jnode
         write (unit=*,fmt='(a,1x,i6)') ' IXB  = ',ixb(jnode,ngr)
         write (unit=*,fmt='(a,1x,i6)') ' IXE  = ',ixe(jnode,ngr)
         write (unit=*,fmt='(a,1x,i6)') ' IYB  = ',iyb(jnode,ngr)
         write (unit=*,fmt='(a,1x,i6)') ' IYE  = ',iye(jnode,ngr)
         call abort_run('Sub-domain would be too small','par_decomp','par_decomp.f90')
      end if
   end do

   return
end subroutine par_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     Sample routine to fill work elements with values proportional to the time required   !
! to perform model operations.                                                             !
!------------------------------------------------------------------------------------------!
subroutine par_est_time(nxp,nyp,work,cput,init)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)  :: nxp
   integer                 , intent(in)  :: nyp
   logical                 , intent(in)  :: init
   real, dimension(nxp,nyp), intent(out) :: work
   real, dimension(nxp,nyp), intent(in)  :: cput
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: i
   integer                               :: j
   !----- Local constants. ----------------------------------------------------------------!
   real                    , parameter   :: bfact = 0.2
   !---------------------------------------------------------------------------------------!


   if (init) then
      do j = 2,nyp-1
         do i = 2,nxp-1
            work(i,j) = 1.
         end do
      end do
   else
      do j = 2,nyp-1
         do i = 2,nxp-1
            work(i,j) = cput(i,j)
         end do
      end do
   end if

   !----- Fill real boundaries with .2 of interior points. --------------------------------!
   do j = 1,nyp
      work(1,j)   = bfact * work(2,j)
      work(nxp,j) = bfact * work(nxp-1,j)
   end do

   do i = 1,nxp
      work(i,1)   = bfact * work(i,2)
      work(i,nyp) = bfact * work(i,nyp-1)
   end do

   call ezcntr(work(1,1),nxp,nyp)

   return
end subroutine par_est_time
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine par_node_paths()
   use cyclic_mod
   use grid_dims , only : maxgrds                         & ! intent(in)
                        , nxpmax                          & ! intent(in)
                        , nypmax                          ! ! intent(in)
   use mem_grid  , only : ngrids                          & ! intent(in)
                        , nnxp                            & ! intent(in)
                        , nnyp                            & ! intent(in)
                        , nxtnest                         & ! intent(in)
                        , ipm                             & ! intent(in)
                        , jpm                             & ! intent(in)
                        , ibnd                            & ! intent(in)
                        , jbnd                            ! ! intent(in)
   use rpara     , only : nmachs                          & ! intent(in)
                        , machnum                         & ! intent(in)
                        , ibcflg                          & ! intent(in)
                        , nxbeg                           & ! intent(in)
                        , nxend                           & ! intent(in)
                        , nybeg                           & ! intent(in)
                        , nyend                           & ! intent(in)
                        , ixb                             & ! intent(in)
                        , ixe                             & ! intent(in)
                        , iyb                             & ! intent(in)
                        , iye                             & ! intent(in)
                        , ibounds                         & ! intent(out)
                        , ipaths    => inode_paths_master & ! intent(out)
                        , igetpaths => iget_paths_master  ! ! intent(out)
                        
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(3) :: is0t
   integer, dimension(3) :: is0u
   integer, dimension(3) :: is0v
   integer, dimension(3) :: js0t
   integer, dimension(3) :: js0u
   integer, dimension(3) :: js0v
   integer               :: ngr
   integer               :: ncm
   integer               :: isend_type
   integer               :: idn
   integer               :: isn
   integer               :: i
   integer               :: j
   integer               :: iwid
   integer               :: info
   integer               :: nnn
   integer               :: indt
   integer               :: indu
   integer               :: indv
   integer               :: nxp
   integer               :: nyp
   integer               :: id
   integer               :: jd
   integer               :: nijst
   integer               :: nijsu
   integer               :: nijsv
   integer               :: iselft
   integer               :: iselfu
   integer               :: iselfv
   integer               :: mijs
   integer               :: is
   integer               :: js
   integer               :: ig
   integer               :: ih
   logical               :: neighbours
   !---------------------------------------------------------------------------------------!



   !------ If using cyclic boundary conditions, allocate ipaths_cyc array. ----------------!
   if (ibnd == 4 .or. jbnd == 4) then
      call ipaths_cyc_alloc(nnxp(1),nnyp(1),ibnd,jbnd)
   end if


   !---------------------------------------------------------------------------------------!
   !     Initialise ipaths array with zeroes.  This is how these arrays work:              !
   !                                                                                       !
   ! IPATHS:                                                                               !
   !  1st dimension : horizontal domain of the exchange between nodes.                     !
   !                  1-4: edges of the exchange (1=West, 2=East, 3=South, 4=North).       !
   !                  5:   MPI number of the node to where the information should go.      !
   !  2nd dimension : exchange type                                                        !
   !                  1 - thermodynamics                                                   !
   !                  2 - zonal wind                                                       !
   !                  3 - meridional wind                                                  !
   !                  4 - perturbation of the Exner function                               !
   !  3rd dimension : grid number                                                          !
   !  4th dimension : source node                                                          !
   !  5th dimension : destination node                                                     !
   !---------------------------------------------------------------------------------------!
   ipaths   (:,:,:,:,:) = 0
   igetpaths(  :,:,:,:) = 0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Grid loop.                                                                        !
   !---------------------------------------------------------------------------------------!
   gridloop:do ngr=1,ngrids
      
      !----- Aliases. ---------------------------------------------------------------------!
      nxp = nnxp(ngr)
      nyp = nnyp(ngr)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   This is the receptor (destination) node.                                         !
      !------------------------------------------------------------------------------------!
      destloop: do idn=1,nmachs

         !------ Initialise the ibounds array with the values that match the full domain. -!
         ibounds(ngr,1,idn) = nxp
         ibounds(ngr,2,idn) = 1
         ibounds(ngr,3,idn) = nyp
         ibounds(ngr,4,idn) = 1
         ibounds(ngr,5,idn) = nxp
         ibounds(ngr,6,idn) = 1
         ibounds(ngr,7,idn) = nyp
         ibounds(ngr,8,idn) = 1
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !   This is the source (send) node.                                               !
         !---------------------------------------------------------------------------------!
         srcloop: do isn = 1,nmachs
            !------------------------------------------------------------------------------!
            !     Check whether the nodes share the edge.  A source node shares the edge   !
            ! with the destination node when they are not the same node, and when          !
            ! the edge of the sub-domain of the destination node overlaps the edge of the  !
            ! computing region of the source point.                                        !
            !------------------------------------------------------------------------------!
            neighbours = (ixe(isn,ngr) >= nxbeg(idn,ngr) .and.                             &
                          iye(isn,ngr) >= nybeg(idn,ngr) .and.                             &
                          ixb(isn,ngr) <= nxend(idn,ngr) .and.                             &
                          iyb(isn,ngr) <= nyend(idn,ngr)       )
            !------------------------------------------------------------------------------!

            if (isn /= idn .and. neighbours) then

               !----- Long timestep overlap regions. --------------------------------------!
               igetpaths(1,ngr,isn,idn) = machnum(isn)

               !---------------------------------------------------------------------------!
               !     Save the edges for the thermodynamic lateral boundary condition.      !
               !---------------------------------------------------------------------------!
               ipaths(1,1,ngr,idn,isn)=max(ixb(isn,ngr),nxbeg(idn,ngr))
               ipaths(2,1,ngr,idn,isn)=min(ixe(isn,ngr),nxend(idn,ngr))
               ipaths(3,1,ngr,idn,isn)=max(iyb(isn,ngr),nybeg(idn,ngr))
               ipaths(4,1,ngr,idn,isn)=min(iye(isn,ngr),nyend(idn,ngr))
               ipaths(5,1,ngr,idn,isn)=machnum(idn)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     In case the boundaries are at the absolute edge of the domain, expand !
               ! ipaths so these edges are also included.                                  !
               !---------------------------------------------------------------------------!
               if (ipaths(1,1,ngr,idn,isn) ==     2) ipaths(1,1,ngr,idn,isn) = 1
               if (ipaths(2,1,ngr,idn,isn) == nxp-1) ipaths(2,1,ngr,idn,isn) = nxp
               if (ipaths(3,1,ngr,idn,isn) ==     2) ipaths(3,1,ngr,idn,isn) = 1
               if (ipaths(4,1,ngr,idn,isn) == nyp-1) ipaths(4,1,ngr,idn,isn) = nyp
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Small timestep overlap regions for zonal wind.                        !
               !---------------------------------------------------------------------------!
               if (ixe(isn,ngr) == ixb(idn,ngr)-1 .and.                                    &
                   iyb(isn,ngr) <= iye(idn,ngr)   .and.                                    &
                   iye(isn,ngr) >= iyb(idn,ngr)         ) then

                  igetpaths(2,ngr,isn,idn) = machnum(isn)

                  ipaths(1,2,ngr,idn,isn)  = ixe(isn,ngr)
                  ipaths(2,2,ngr,idn,isn)  = ixe(isn,ngr)
                  ipaths(3,2,ngr,idn,isn)  = max(iyb(isn,ngr),iyb(idn,ngr))
                  ipaths(4,2,ngr,idn,isn)  = min(iye(isn,ngr),iye(idn,ngr))
                  ipaths(5,2,ngr,idn,isn)  = machnum(idn)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Small timestep overlap regions for meridional wind                    !
               !---------------------------------------------------------------------------!
               if (iye(isn,ngr) == iyb(idn,ngr)-1 .and.                                    &
                   ixb(isn,ngr) <= ixe(idn,ngr)   .and.                                    &
                   ixe(isn,ngr) >= ixb(idn,ngr)         ) then

                  igetpaths(3,ngr,isn,idn) = machnum(isn)

                  ipaths(1,3,ngr,idn,isn)  = max(ixb(isn,ngr),ixb(idn,ngr))
                  ipaths(2,3,ngr,idn,isn)  = min(ixe(isn,ngr),ixe(idn,ngr))
                  ipaths(3,3,ngr,idn,isn)  = iye(isn,ngr)
                  ipaths(4,3,ngr,idn,isn)  = iye(isn,ngr)
                  ipaths(5,3,ngr,idn,isn)  = machnum(idn)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Small timestep overlap regions for the perturbation of the Exner      !
               ! function.                                                                 !
               !---------------------------------------------------------------------------!
               if (ixb(isn,ngr)  == ixe(idn,ngr)+1 .and.                                   &
                   iyb(isn,ngr)  <= iye(idn,ngr)   .and.                                   &
                   iye(isn,ngr)  >= iyb(idn,ngr)         ) then

                  igetpaths(4,ngr,isn,idn) = machnum(isn)

                  ipaths(1,4,ngr,idn,isn)  = ixb(isn,ngr)
                  ipaths(2,4,ngr,idn,isn)  = ixb(isn,ngr)
                  ipaths(3,4,ngr,idn,isn)  = max(iyb(isn,ngr),iyb(idn,ngr))
                  ipaths(4,4,ngr,idn,isn)  = min(iye(isn,ngr),iye(idn,ngr))
                  ipaths(5,4,ngr,idn,isn)  = machnum(idn)

               elseif (iyb(isn,ngr) == iye(idn,ngr)+1 .and.                                &
                       ixb(isn,ngr) <= ixe(idn,ngr)   .and.                                &
                       ixe(isn,ngr) >= ixb(idn,ngr)         ) then

                  igetpaths(4,ngr,isn,idn) = machnum(isn)

                  ipaths(1,4,ngr,idn,isn)  = max(ixb(isn,ngr),ixb(idn,ngr))
                  ipaths(2,4,ngr,idn,isn)  = min(ixe(isn,ngr),ixe(idn,ngr))
                  ipaths(3,4,ngr,idn,isn)  = iyb(isn,ngr)
                  ipaths(4,4,ngr,idn,isn)  = iyb(isn,ngr)
                  ipaths(5,4,ngr,idn,isn)  = machnum(idn)
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Check whether there should be communication between nested grids.       !
            ! There is a finer grid when nxtnest(ngr) is not zero.                         !
            !------------------------------------------------------------------------------!
            nnn = nxtnest(ngr)
            if (nnn /= 0) then 

               !---------------------------------------------------------------------------!
               !     The following flag will control the communication from coarse to fine !
               ! grid.                                                                     !
               !---------------------------------------------------------------------------!
               if (ixe(isn,nnn) >= ipm(ixb(idn,ngr)-1,ngr)-2 .and.                         &
                   iye(isn,nnn) >= jpm(iyb(idn,ngr)-1,ngr)-2 .and.                         &
                   ixb(isn,nnn) <= ipm(ixe(idn,ngr)+1,ngr)+1 .and.                         &
                   iyb(isn,nnn) <= jpm(iye(idn,ngr)+1,ngr)+1       ) then

                  igetpaths(5,ngr,isn,idn) = machnum(isn)

                  ipaths(1,5,ngr,idn,isn)  = max(ixb(isn,nnn),ipm(ixb(idn,ngr)-1,ngr)-2)
                  ipaths(2,5,ngr,idn,isn)  = min(ixe(isn,nnn),ipm(ixe(idn,ngr)+1,ngr)+1)
                  ipaths(3,5,ngr,idn,isn)  = max(iyb(isn,nnn),jpm(iyb(idn,ngr)-1,ngr)-2)
                  ipaths(4,5,ngr,idn,isn)  = min(iye(isn,nnn),jpm(iye(idn,ngr)+1,ngr)+1)
                  ipaths(5,5,ngr,idn,isn)  = machnum(idn)
               end if
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               ! Fine grid to coarse grid averaging communication.                         !
               !                                                                           !
               ! rtimh - micphys: in the following lines: change ixb, ixe, iyb, and iye to !
               !                  nxbeg, nxend, nybeg, and nyend on coarse grid            !
               !                  destination nodes in order to keep internal CG           !
               !                  boundaries up to date, since they have already been      !
               !                  updated by their CG neighbouring nodes??                 !
               !---------------------------------------------------------------------------!
               if (nxbeg(idn,nnn) <= ipm(ixe(isn,ngr),ngr) .and.                           &
                   nybeg(idn,nnn) <= jpm(iye(isn,ngr),ngr) .and.                           &
                   nxend(idn,nnn) >= ipm(ixb(isn,ngr),ngr) .and.                           &
                   nyend(idn,nnn) >= jpm(iyb(isn,ngr),ngr)       ) then

                  igetpaths(6,ngr,isn,idn) = machnum(isn)

                  ipaths(1,6,ngr,idn,isn)  = max(ipm(ixb(isn,ngr),ngr),nxbeg(idn,nnn))
                  ipaths(2,6,ngr,idn,isn)  = min(ipm(ixe(isn,ngr),ngr),nxend(idn,nnn))
                  ipaths(3,6,ngr,idn,isn)  = max(jpm(iyb(isn,ngr),ngr),nybeg(idn,nnn))
                  ipaths(4,6,ngr,idn,isn)  = min(jpm(iye(isn,ngr),ngr),nyend(idn,nnn))
                  ipaths(5,6,ngr,idn,isn)  = machnum(idn)

                  !------ This is used for binary reproducibility. ------------------------!
                  if (ngr >= 2) then
                     ncm = ngr - 1
                     ibounds(ncm,5,idn) = min(ibounds(ncm,5,idn),ipaths(1,6,ngr,idn,isn))
                     ibounds(ncm,6,idn) = max(ibounds(ncm,6,idn),ipaths(2,6,ngr,idn,isn))
                     ibounds(ncm,7,idn) = min(ibounds(ncm,7,idn),ipaths(3,6,ngr,idn,isn))
                     ibounds(ncm,8,idn) = max(ibounds(ncm,8,idn),ipaths(4,6,ngr,idn,isn))
                  end if
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     A second index value of 7 of the ipaths array is used to determine !
                  ! the loop limits in fdbackp for averaging the fine mesh over the over-  !
                  ! lap between the coarse mesh node and fine mesh node, rather than       !
                  ! always over the full fine mesh node.  It is not used for actually      !
                  ! sending stuff.  The ipaths(*,6,*,*,*) part of the array is still used  !
                  ! for sending the block of averaged coarse mesh points from the fine     !
                  ! mesh node to the coarse mesh node.                                     !
                  !------------------------------------------------------------------------!
                  forwiloop: do i = ixb(isn,ngr),ixe(isn,ngr)
                     if (ipm(i,ngr) == ipaths(1,6,ngr,idn,isn)) then
                        ipaths(1,7,ngr,idn,isn) = i
                        exit forwiloop
                     end if
                  end do forwiloop

                  backiloop: do i = ixe(isn,ngr),ixb(isn,ngr),-1
                     if (ipm(i,ngr) == ipaths(2,6,ngr,idn,isn)) then
                        ipaths(2,7,ngr,idn,isn) = i
                        exit backiloop
                     end if
                  end do backiloop

                  forwjloop: do j = iyb(isn,ngr),iye(isn,ngr)
                     if (jpm(j,ngr) == ipaths(3,6,ngr,idn,isn)) then
                        ipaths(3,7,ngr,idn,isn) = j
                        exit forwjloop
                     end if
                  end do forwjloop

                  backjloop: do j = iye(isn,ngr),iyb(isn,ngr),-1
                     if (jpm(j,ngr) == ipaths(4,6,ngr,idn,isn)) then
                        ipaths(4,7,ngr,idn,isn) = j
                        exit backjloop
                     end if
                  end do backjloop
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !----- Update ibounds, which will be used for binary reproducibility. ---------!
            if (ipaths(1,7,ngr,idn,isn) /= 0) then
               ibounds(ngr,1,idn) = min(ibounds(ngr,1,idn),ipaths(1,7,ngr,idn,isn))
               ibounds(ngr,2,idn) = max(ibounds(ngr,2,idn),ipaths(2,7,ngr,idn,isn))
               ibounds(ngr,3,idn) = min(ibounds(ngr,3,idn),ipaths(3,7,ngr,idn,isn))
               ibounds(ngr,4,idn) = max(ibounds(ngr,4,idn),ipaths(4,7,ngr,idn,isn))
            end if
            !------------------------------------------------------------------------------!
         end do srcloop
      end do destloop
   end do gridloop
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Cyclic boundary conditions (grid 1 only).                                        !
   !---------------------------------------------------------------------------------------!
   indt = 0
   indu = 0
   indv = 0
   nxp = nnxp(1)
   nyp = nnyp(1)


   !----- Loop over all destination ij columns. -------------------------------------------!
   jcycloop: do jd = 1,nyp
      icycloop: do id = 1,nxp

         if (id <= 2 .or. id >= nxp-1 .or. jd <= 2 .or. jd >= nyp-1) then
            !------------------------------------------------------------------------------!
            !     Determine 1, 2, or 3 source ij columns for this destination ij column.   !
            !------------------------------------------------------------------------------!
            nijst = 0
            nijsu = 0
            nijsv = 0
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Cyclic in the longitude.                                                 !
            !------------------------------------------------------------------------------!
            if (ibnd == 4) then
               !----- Thermodynamic point. ------------------------------------------------!
               if (id <= 2 .and. jd >= 2 .and. jd <= nyp-1) then
                  call fill_cyc(nijst,is0t,js0t,id+nxp-3,jd)
               elseif (id >= nxp-1 .and. jd >= 2 .and. jd <= nyp-1) then
                  call fill_cyc(nijst,is0t,js0t,id-nxp+3,jd)
               end if
               !---------------------------------------------------------------------------!



               !----- Zonal wind. ---------------------------------------------------------!
               if (id == 1 .and. jd >= 2 .and. jd <= nyp-1) then
                  call fill_cyc(nijsu,is0u,js0u,id+nxp-3,jd)
               elseif (id == nxp-1 .and. jd >= 2 .and. jd <= nyp-1) then
                  call fill_cyc(nijsu,is0u,js0u,id-nxp+3,jd)
               end if
               !---------------------------------------------------------------------------!



               !----- Meridional wind. ----------------------------------------------------!
               if (id <= 2 .and. jd >= 2 .and. jd < nyp-1) then
                  call fill_cyc(nijsv,is0v,js0v,id+nxp-3,jd)
               elseif (id >= nxp-1 .and. jd >= 2 .and. jd < nyp-1) then
                  call fill_cyc(nijsv,is0v,js0v,id-nxp+3,jd)
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Cyclic in the latitude.                                                  !
            !------------------------------------------------------------------------------!
            if (jbnd == 4) then
               !----- Thermodynamic. ------------------------------------------------------!
               if (jd <= 2 .and. id >= 2 .and. id <= nxp-1) then
                  call fill_cyc(nijst,is0t,js0t,id,jd+nyp-3)
               elseif (jd >= nyp-1 .and. id >= 2 .and. id <= nxp-1) then
                  call fill_cyc(nijst,is0t,js0t,id,jd-nyp+3)
               end if
               !---------------------------------------------------------------------------!



               !----- Zonal wind. ---------------------------------------------------------!
               if (jd <= 2 .and. id >= 2 .and. id < nxp-1) then
                  call fill_cyc(nijsu,is0u,js0u,id,jd+nyp-3)
               elseif (jd >= nyp-1 .and. id >= 2 .and. id < nxp-1) then
                  call fill_cyc(nijsu,is0u,js0u,id,jd-nyp+3)
               end if
               !---------------------------------------------------------------------------!



               !----- Meridional wind. ----------------------------------------------------!
               if (jd == 1 .and. id >= 2 .and. id <= nxp-1) then
                  call fill_cyc(nijsv,is0v,js0v,id,jd+nyp-3)
               elseif (jd == nyp-1 .and. id >= 2 .and. id <= nxp-1) then
                  call fill_cyc(nijsv,is0v,js0v,id,jd-nyp+3)
               end if
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Cyclic in the longitude and latitude, check corners.                     !
            !------------------------------------------------------------------------------!
            if (ibnd == 4 .and. jbnd == 4) then
               !----- Double corners (2 x 2 points each corner). --------------------------!



               !----- Thermodynamic. ------------------------------------------------------!
               if (id <= 2 .and. jd <= 2) then
                  call fill_cyc(nijst,is0t,js0t,id+nxp-3,jd+nyp-3)
               elseif (id >= nxp-1 .and. jd <= 2) then
                  call fill_cyc(nijst,is0t,js0t,id-nxp+3,jd+nyp-3)
               elseif (id <= 2 .and. jd >= nyp-1) then
                  call fill_cyc(nijst,is0t,js0t,id+nxp-3,jd-nyp+3)
               elseif (id >= nxp-1 .and. jd >= nyp-1) then
                  call fill_cyc(nijst,is0t,js0t,id-nxp+3,jd-nyp+3)
               end if
               !---------------------------------------------------------------------------!



               !----- Zonal wind. ---------------------------------------------------------!
               if (id == 1 .and. jd <= 2) then
                  call fill_cyc(nijsu,is0u,js0u,id+nxp-3,jd+nyp-3)
               elseif (id == nxp-1 .and. jd <= 2) then
                  call fill_cyc(nijsu,is0u,js0u,id-nxp+3,jd+nyp-3)
               elseif (id == 1 .and. jd >= nyp-1) then
                  call fill_cyc(nijsu,is0u,js0u,id+nxp-3,jd-nyp+3)
               elseif (id == nxp-1 .and. jd >= nyp-1) then
                  call fill_cyc(nijsu,is0u,js0u,id-nxp+3,jd-nyp+3)
               end if
               !---------------------------------------------------------------------------!



               !----- Meridional wind. ----------------------------------------------------!
               if (id <= 2 .and. jd == 1) then
                  call fill_cyc(nijsv,is0v,js0v,id+nxp-3,jd+nyp-3)
               elseif (id >= nxp-1 .and. jd == 1) then
                  call fill_cyc(nijsv,is0v,js0v,id-nxp+3,jd+nyp-3)
               elseif (id <= 2 .and. jd == nyp-1) then
                  call fill_cyc(nijsv,is0v,js0v,id+nxp-3,jd-nyp+3)
               elseif (id >= nxp-1 .and. jd == nyp-1) then
                  call fill_cyc(nijsv,is0v,js0v,id-nxp+3,jd-nyp+3)
               end if

            elseif (ibnd == 4) then
               !----- One corner (single point each corner). ------------------------------!



               !----- Thermodynamic, and zonal and meridional winds. ----------------------!
               if (id == 1 .and. (jd == 1 .or. jd == nyp)) then
                  call fill_cyc(nijst,is0t,js0t,id+nxp-3,jd)
                  call fill_cyc(nijsu,is0u,js0u,id+nxp-3,jd)
                  call fill_cyc(nijsv,is0v,js0v,id+nxp-3,jd)

               elseif (id == nxp .and. (jd == 1 .or. jd == nyp)) then
                  call fill_cyc(nijst,is0t,js0t,id-nxp+3,jd)
                  call fill_cyc(nijsu,is0u,js0u,id-nxp+3,jd)
                  call fill_cyc(nijsv,is0v,js0v,id-nxp+3,jd)

               end if

            elseif (jbnd == 4) then
               !----- One corner (single point each corner). ------------------------------!



               !----- Thermodynamic, and zonal and meridional winds. ----------------------!
               if ((id == 1 .or. id == nxp) .and. jd == 1) then
                  call fill_cyc(nijst,is0t,js0t,id,jd+nyp-3)
                  call fill_cyc(nijsu,is0u,js0u,id,jd+nyp-3)
                  call fill_cyc(nijsv,is0v,js0v,id,jd+nyp-3)
               elseif ((id == 1 .or. id == nxp) .and. jd == nyp) then
                  call fill_cyc(nijst,is0t,js0t,id,jd-nyp+3)
                  call fill_cyc(nijsu,is0u,js0u,id,jd-nyp+3)
                  call fill_cyc(nijsv,is0v,js0v,id,jd-nyp+3)
               end if
            end if


            !------------------------------------------------------------------------------!
            !     Define "self" message passing flags.                                     !
            !------------------------------------------------------------------------------!
            if (id == 1 .or. id == nxp .or. jd == 1 .or. jd == nyp) then
               iselft = 0
               iselfu = 0
               iselfv = 0

            elseif (id == nxp - 1) then
               iselft = 1
               iselfu = 0
               iselfv = 1
               if (jd == nyp - 1) iselfv = 0

            elseif (jd == nyp - 1) then
               iselft = 1
               iselfu = 1
               iselfv = 0

            else
               iselft = 1
               iselfu = 1
               iselfv = 1

            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Loop over all destination nodes and check if current ij destination      !
            ! column is in it.                                                             !
            !------------------------------------------------------------------------------!
            do idn = 1,nmachs
               if (id >= nxbeg(idn,1) .and. id <= nxend(idn,1) .and.                       &
                   jd >= nybeg(idn,1) .and. jd <= nyend(idn,1)       ) then

                  !----- Loop over up to 3 source ij columns. -----------------------------!
                  do mijs = 1,nijst
                     is = is0t(mijs)
                     js = js0t(mijs)

                     !---------------------------------------------------------------------!
                     !      Loop over all source nodes and check if current ij source      !
                     ! column is in it.                                                    !
                     !---------------------------------------------------------------------!
                     do isn = 1,nmachs
                        if (is >= ixb(isn,1) .and. is <= ixe(isn,1) .and.                  &
                            js >= iyb(isn,1) .and. js <= iye(isn,1)       ) then

                           indt                = indt + 1
                           ipathst_cyc(1,indt) = isn        ! source node #
                           ipathst_cyc(2,indt) = is         ! source node i
                           ipathst_cyc(3,indt) = js         ! source node j
                           ipathst_cyc(4,indt) = idn        ! destination node #
                           ipathst_cyc(5,indt) = id         ! destination node i
                           ipathst_cyc(6,indt) = jd         ! destination node j
                           ipathst_cyc(7,indt) = iselft     ! self count
                           ipathst_cyc(8,indt) = nijst      ! other point count
                        end if
                     end do
                  end do
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  do mijs = 1,nijsu
                     is = is0u(mijs)
                     js = js0u(mijs)


                     !---------------------------------------------------------------------!
                     !      Loop over all source nodes and check if current ij source      !
                     ! column is in it.                                                    !
                     !---------------------------------------------------------------------!
                     do isn = 1,nmachs
                        if (is >= ixb(isn,1) .and. is <= ixe(isn,1) .and.                  &
                            js >= iyb(isn,1) .and. js <= iye(isn,1)       ) then

                           indu                = indu + 1
                           ipathsu_cyc(1,indu) = isn        ! source node #
                           ipathsu_cyc(2,indu) = is         ! source node i
                           ipathsu_cyc(3,indu) = js         ! source node j
                           ipathsu_cyc(4,indu) = idn        ! destination node #
                           ipathsu_cyc(5,indu) = id         ! destination node i
                           ipathsu_cyc(6,indu) = jd         ! destination node j
                           ipathsu_cyc(7,indu) = iselfu     ! self count
                           ipathsu_cyc(8,indu) = nijsu      ! other point count
                        end if
                     end do
                  end do
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  do mijs = 1,nijsv
                     is = is0v(mijs)
                     js = js0v(mijs)

                     !---------------------------------------------------------------------!
                     !      Loop over all source nodes and check if current ij source      !
                     ! column is in it.                                                    !
                     !---------------------------------------------------------------------!
                    do isn = 1,nmachs
                        if (is >= ixb(isn,1) .and. is <= ixe(isn,1) .and.                  &
                            js >= iyb(isn,1) .and. js <= iye(isn,1)       ) then

                           indv                = indv + 1
                           ipathsv_cyc(1,indv) = isn        ! source node #
                           ipathsv_cyc(2,indv) = is         ! source node i
                           ipathsv_cyc(3,indv) = js         ! source node j
                           ipathsv_cyc(4,indv) = idn        ! destination node #
                           ipathsv_cyc(5,indv) = id         ! destination node i
                           ipathsv_cyc(6,indv) = jd         ! destination node j
                           ipathsv_cyc(7,indv) = iselfv     ! self count
                           ipathsv_cyc(8,indv) = nijsv      ! other point count
                        end if
                     end do
                  end do

               end if
            end do

         end if
      end do icycloop
   end do jcycloop

   return
end subroutine par_node_paths
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_cyc(nijs,is0,js0,is,js)
  implicit none
  integer :: nijs,is,js
  integer, dimension(*) :: is0,js0

  nijs = nijs + 1
  is0(nijs) = is
  js0(nijs) = js

  return
end subroutine fill_cyc

!******************************************************************************

subroutine decomp_plot(ngr)

  use mem_grid
  use rpara

  implicit none
  integer :: ngr
  character(len=2) :: cnum
  integer :: nm,nb,nn,ixxb,ixxe,iyyb,iyye,ng
  real :: pxb,pxe,pyb,pye,px,py

  call set(.1, .9,.1,.9,  &
       xtn(1,ngr),xtn(nnxp(ngr),ngr)  &
       ,ytn(1,ngr),ytn(nnyp(ngr),ngr),1)
  call perim(nnxp(ngr)-1,1,nnyp(ngr)-1,1)

  do nm=1,nmachs

     nb=2
     if(ngr.eq.ngrids) nb=1

     do nn=1,nb
        if(nn.eq.1) ng=ngr
        if(nn.eq.2) ng=ngr+1

        ixxb=ixb(nm,ng)-1
        ixxe=ixe(nm,ng)
        iyyb=iyb(nm,ng)-1
        iyye=iye(nm,ng)
        pxb=xtn(ixxb,ng)
        pxe=xtn(ixxe,ng)
        pyb=ytn(iyyb,ng)
        pye=ytn(iyye,ng)

        call line(pxb,pyb,pxe,pyb)
        call line(pxe,pyb,pxe,pye)
        call line(pxe,pye,pxb,pye)
        call line(pxb,pye,pxb,pyb)

        px=.5*(pxb+pxe)
        py=.5*(pyb+pye)
        write(cnum,'(i2)') nm
        call pwritx(px,py,cnum,2,15,0,0.)
     enddo

  enddo

  call frame
  return
end subroutine decomp_plot

!*******************************************************************************

subroutine PAR_decomp_input (maxmach, maxgrds, nmachs, ngr, nnxp, nnyp, ixb, ixe, iyb, iye, failed)

  ! par_decomp_input: 
  !  reads domain decomposition from input file

  use domain_decomp, only: &
       domain_fname         ! intent(in)

  implicit none
  integer, intent(in)  :: maxmach
  integer, intent(in)  :: maxgrds
  integer, intent(in)  :: nmachs
  integer, intent(in)  :: ngr
  integer, intent(in)  :: nnxp(maxgrds)
  integer, intent(in)  :: nnyp(maxgrds)
  integer, intent(out) :: ixb(maxmach,maxgrds)
  integer, intent(out) :: ixe(maxmach,maxgrds)
  integer, intent(out) :: iyb(maxmach,maxgrds)
  integer, intent(out) :: iye(maxmach,maxgrds)
  logical, intent(out) :: failed

  integer, parameter :: un=48
  logical :: ex
  integer :: err
  integer :: grid
  integer :: node
  integer :: xbeg
  integer :: xend
  integer :: ybeg
  integer :: yend
  integer :: occur
  integer :: igr
  integer :: imachs
  integer :: iy
  integer :: ix
  integer, parameter :: UndefIndex=-234
  !character(len=*), parameter :: fname="Domains"
  character(len=*), parameter :: h="**(par_decomp_input)**"
  character(len=10) :: c0, c1, c2
  logical, parameter :: dumpLocal=.FALSE.

  ! fail is the default
  failed = .true.

  ! Checking if domain_fname is valid
  if (len_trim(domain_fname)==0) then
     if (dumpLocal) then
        write(*,"(a)") h//" Explicit domain decomposition file invalid!"
        write(*,"(a)") h//" Original domain decomposition mode activated"
     end if
     return
  endif

  if (dumpLocal) then
     write(*,"(a)") h//" fetching file "//domain_fname
  end if

  ! Checking file existence
  inquire(file=domain_fname, exist=ex)
  if (.not. ex) then
     write(*,"(a)") h//" Explicit domain decomposition file "//domain_fname// &
          " not found"
     write(*,"(a)") h//" Original domain decomposition mode activated"
     return
  end if

  ! initialization
  ixb(1:nmachs, 1:ngr) = UndefIndex
  ixe(1:nmachs, 1:ngr) = UndefIndex
  iyb(1:nmachs, 1:ngr) = UndefIndex
  iye(1:nmachs, 1:ngr) = UndefIndex

  ! store file contents at desired places
  open(un, file=domain_fname, status="old", action="read")
  do

     ! read until eof
     read(un, *, iostat=err) grid, node, xbeg, xend, ybeg, yend
     if (err /= 0) then
        exit
     end if

     ! Checking data consistency
     if (grid < 1 .or. grid > ngr) then
        write(c0,"(i10)") grid
        write(c1,"(i10)") ngr
        write(*,"(a)") h//" unexpected grid number on input file "//domain_fname
        write(*,"(a)") h//" grid range [1:"//trim(adjustl(c1))//&
             &"]; unexpected grid number="//trim(adjustl(c0))
        write(*,"(a)") h//" Explicit domain decomposition file invalid!"
        write(*,"(a)") h//" Original domain decomposition mode activated"
        return
     end if
     if (node < 1 .or. node > nmachs) then
        write(c0,"(i10)") node
        write(c1,"(i10)") nmachs
        write(*,"(a)") h//" unexpected node number on input file "//domain_fname 
        write(*,"(a)") h//" node range [1:"//trim(adjustl(c1))//&
             &"]; unexpected node number="//trim(adjustl(c0))
        write(*,"(a)") h//" Explicit domain decomposition file invalid!"
        write(*,"(a)") h//" Original domain decomposition mode activated"
        return
     end if

     ! store

     ixb(node,grid) = xbeg
     ixe(node,grid) = xend
     iyb(node,grid) = ybeg
     iye(node,grid) = yend
  end do
  close(un)
  
  ! all data stored?
  ! it suffices to verify ixb (since ixe,iyb,iye are at the same read)
  do igr = 1, ngr
     do imachs = 1, nmachs
        if (ixb(imachs,igr) == UndefIndex) then
           write(c0,"(i10)") imachs
           write(c1,"(i10)") igr
           write(*,"(a)") h//" missing data for machine "&
                &//trim(adjustl(c0))//" and grid "//trim(adjustl(c1))
           write(*,"(a)") h//" on input file "//domain_fname
           write(*,"(a)") h//" Explicit domain decomposition file invalid!"
           write(*,"(a)") h//" Original domain decomposition mode activated"
           return
        end if
     end do
  end do

  ! input data partitions the domain?
  do igr = 1, ngr
     do iy = 2, nnyp(igr)-1
        do ix = 2, nnxp(igr)-1
           occur = count(&
                (ixb(1:nmachs,igr) <= ix) .and. &
                (ixe(1:nmachs,igr) >= ix) .and. &
                (iyb(1:nmachs,igr) <= iy) .and. &
                (iye(1:nmachs,igr) >= iy)          )
           if (occur /= 1) then
              write(*,"(a)") h//" domain decomposition is not a partition on file "//domain_fname
              write(*,"(a,3i10)") h//" point igr, ix, iy=",igr,ix,iy
              write(*,"(a,i10)") h//" number of occurences=",occur
              write(*,"(a)") h//" Explicit domain decomposition file invalid!"
              write(*,"(a)") h//" Original domain decomposition mode activated"
              return
           end if
        end do
     end do
  end do

  ! accepted domain decomposition by external file
  failed = .false.
  if (dumpLocal) then
     write(*,"(a)") h//" domain decomposition imposed by file "//domain_fname
  end if
end subroutine PAR_decomp_input
