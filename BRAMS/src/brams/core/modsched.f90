!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)

  ! modsched purpose: 
  !  (1) schedules nested time step iterations of multiple grids for synchronous 
  !      time advance;
  !  (2) schedules data exchange among grids to keep fields consistent among grids
  !
  ! modsched assumptions: 
  ! Grids are organized in a tree with:
  !  (1) exactly one node for each grid and 
  !  (2) all grids at any nonempty subtree are fully nested at the root grid
  ! 
  ! Grid nesting geography is restricted by:
  !  (a) There is a single outermost grid that contains all nested grids;
  !  (b) Sibling grids are disjoint (including boundary cells);
  !  (c) If grid b is nested into grid a, then a fixed integer number
  !      of gird b cells partitions a cell of grid a 
  !  (d) If grid b is nested into grid a, then the inner part
  !      of grid b (grid b without boundary cells) fully partitions a 
  !      set of inner cells of grid a.
  !
  ! Grid nested timestep is restricted by:
  !  (a) A nested grid timestep is an integer fraction of the parent's grid timestep.
  !
  ! modsched reasoning: 
  !   A single timestep of any grid requires multiple timesteps of nested grids 
  !   (since nested grids have a smaller timestep than coarser grids). Once a grid
  !   runs a timestep, all nested grids should run a set of timesteps to synchronous
  !   advance.
  !   To maintain fields consistent among grids, nested grids should receive
  !   boundary conditions from the parent grid (since the parent -- coarser -- grid is
  !   advanced in time) before running its set of own timesteps. After synchronizing
  !   in time with the parent grid, it should feedback the parent with all detailed 
  !   fields to be integrated back into the parent grid (feedback fields).
  !
  ! modsched algorithm, giving the tree and the timestep ratio of each grid 
  ! to the parent's grid:
  !   Noting to do on an empty tree;
  !   Visit the tree root as many times as the timestep ratio (for time synchronization
  !   of the root). At each visit, schedule a timestep, send boundary conditions 
  !   to all direct sons (to keep fields consistent) and recursivelly run the algorithm
  !   at all sub-trees rooted at sons. Finally, if the root has a parent grid, send feedback
  !   fields to keep parent's field consistent.
  !
  ! modsched data structure:
  !   nsubs:       is the total number of nested timesteps required to advance the full
  !                tree a single timestep;
  !   isched:      what to do at nested timestep i (first index, 1 <= i <= nsubs);
  !   isched(i,1): grid to run this nested timestep 
  !   isched(i,2): number of grids to send boundary conditions
  !   isched(i,3): nested timestep counter for this nested timestep grid 
  !                during a single timestep of the parent's grid
  !   isched(i,4): number of grids to send feedback fields 
  !   isched(i,5): nested timestep counter for this nested timestep grid 
  !                during a single timestep of grid 1 (that is, among all
  !                nested timesteps of a single outermost grid timestep)
  !
  ! modsched software assumptions: (not verified by modshched:)
  !   ngrids > 0                       (at least grid 1 exists)
  !   nxtnest(1)=0                     (grid 1 is the single outermost grid)
  !   1 <= nxtnest(2:ngrids) <= ngrids (tree is self contained and is not a forest)
  !   nndtrat(1)=1                     (grid 1 is timestep ratio reference)
  !   nndtrat(2:ngrids) >= 2           (nested grids should have smaller dt);
  !   maxsched >= ngrids               (space allocation)
  !   maxschent >= 5                   (space allocation)

  implicit none
  integer, intent(in)  :: maxsched                   ! maximum number of grids in the tree of nested grids
  integer, intent(in)  :: maxschent                  ! number of features per nested timestep stored at isched
  integer, intent(in)  :: ngrids                     ! number of grids in the tree of nested grids
  integer, intent(in)  :: nxtnest(ngrids)            ! parent grid: tree representation
  integer, intent(in)  :: nndtrat(ngrids)            ! delta t ratio: parent grid/this grid
  integer, intent(out) :: isched(maxsched,maxschent) ! scheduling table (see above)
  integer, intent(out) :: nsubs                      ! total number of nested timesteps for a dtlong

  integer :: nTS(ngrids)      ! number of nested timesteps at each grid

  ! initialize counters

  isched(:,:) = 0
  nsubs = 0
  nTS(:) = 0

  ! visit tree in pre-order to schedule timesteps
  ! (assumes node 1 is the root of a single tree)

  call operTS(1, isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs, nTS)
end subroutine modsched

! *************************************************************************

recursive subroutine operTS(root, isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs, nTS)

  ! operTS: 
  !   visit the tree of grids with root "root" in preorder;
  !   requires (1 <= root <= ngrids);

  implicit none
  integer, intent(in)    :: root                       ! grid id of the sub-tree root
  integer, intent(in)    :: maxsched                   ! maximum number of grids in the tree of nested grids
  integer, intent(in)    :: maxschent                  ! number of features per timestep stored at isched
  integer, intent(in)    :: ngrids                     ! number of grids in the tree of nested grids
  integer, intent(in)    :: nxtnest(ngrids)            ! parent grid 
  integer, intent(in)    :: nndtrat(ngrids)            ! delta t ratio: parent grid/this grid
  integer, intent(inout) :: isched(maxsched,maxschent) ! scheduling table
  integer, intent(inout) :: nsubs                      ! total number of timesteps for a dtlong
  integer, intent(inout) :: nTS(ngrids)                ! timesteps of each grid

  integer :: iter  ! nested timestep count for the current grid
  integer :: node  ! auxiliar grid count
  character(len=10) :: c0, c1
  character(len=*), parameter :: h="**(operTS)**"

  ! should be a non-empty tree

  if (root <= 0 .or. root > ngrids) then
     write(c0,"(i10)") root
     write(c0,"(i10)") ngrids
     write(*,"(a)") h//" INTERNAL ERROR: invoked with root ="//&
          &trim(adjustl(c0))//" outside bounds [1:"//&
          &trim(adjustl(c1))//"]"
     stop
  end if

  ! run all nested timesteps to synchronize

  do iter = 1, nndtrat(root)
     
     ! check bounds
     ! schedule execution of this grid nested timestep
     
     nsubs = nsubs + 1
     if (nsubs > maxsched) then
        write(c0,"(i10)") maxsched
        write(*,"(a)") h//" INTERNAL ERROR: increment exausted maxsched ="//&
             &trim(adjustl(c0))
        stop
     end if
     isched(nsubs,1) = root
     
     ! nested timestep counters relative to next coarser grid and to most coarser grid
     
     isched(nsubs,3) = iter
     nTS(root) = nTS(root)+1
     isched(nsubs,5) = nTS(root)
     
     ! # grids to receive boundary conditions
     ! (these are the direct sons)
     
     isched(nsubs,2) = count(nxtnest(1:ngrids) == root)
     
     ! schedule execution of trees rooted at each son
     
     do node = 1, ngrids
        if (nxtnest(node) == root) then
           call operTS(node, isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs, nTS)
        end if
     end do
  end do
  
  ! at the adequate nested timestep iteration (after all finer grids are synchronized)
  ! increase count of sending feedback fields to the father, if any.
  
  if (nxtnest(root) /= 0) then
     isched(nsubs,4) = isched(nsubs,4) + 1
  end if
end subroutine operTS

! *************************************************************************

subroutine dump_modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
  implicit none
  integer, intent(in) :: maxsched                   ! maximum number of grids in the tree of nested grids
  integer, intent(in) :: maxschent                  ! number of features per timestep stored at isched
  integer, intent(in) :: ngrids                     ! number of grids in the tree of nested grids
  integer, intent(in) :: nxtnest(ngrids)            ! parent grid 
  integer, intent(in) :: nndtrat(ngrids)            ! delta t ratio: parent grid/this grid
  integer, intent(in) :: isched(maxsched,maxschent) ! scheduling table (see below)
  integer, intent(in) :: nsubs

  integer :: is                   ! nested timestep count
  integer :: grid                 ! grid at this nested timestep
  integer :: ifg                  ! some finner grid
  integer :: icg                  ! some coarser grid
  integer :: ng                   ! auxiliar grid count
  integer :: bcgrid(ngrids)       ! grids to receive bc
  integer :: dtratouter(ngrids)   ! # grid timesteps to advance dtlong
  character(len=10) :: c0, c1, c2, c3, c4, c5, c6
  character(len=*), parameter :: h="**(dump_modsched)**"

  ! grid timesteps to advance dtlong

  dtratouter = 1
  do grid = 1, ngrids
     icg = grid
     do
        dtratouter(grid) = nndtrat(icg)*dtratouter(grid)
        icg = nxtnest(icg)
        if (icg == 0) exit
     end do
  end do

  ! banner

  write(*,"(a)")
  write(c0,"(i10)") nsubs
  write(*,"(a)") " === Timestep Schedule requires "//trim(adjustl(c0))//&
       &" nested timesteps: ===="

  ! for each timestep:

  do is = 1, nsubs
     grid = isched(is,1)

     ! dump grid to execute timestep and timestep counters

     write(c0,"(i10)") dtratouter(grid)
     write(c1,"(i10)") isched(is,1)
     write(c2,"(i10)") isched(is,3)
     write(c3,"(i10)") isched(is,5)
     write(c4,"(i10)") nndtrat(grid)
     write(c5,"(i10)") nxtnest(grid)
     write(c6,"(i10)") is
     write(*,"(a)", advance="NO") &
          &" nts "//trim(adjustl(c6))//&
          &": Grid "//trim(adjustl(c1))
     if (dtratouter(grid) == 1) then
        write(*,"(a)", advance="NO") " advances dtlong"
     else if (dtratouter(grid) == isched(is,5)) then
        write(*,"(a)", advance="NO") " reaches dtlong"
     else if (isched(is,3) == nndtrat(grid)) then
        write(*,"(a)", advance="NO") " reaches "//&
             &trim(adjustl(c3))//"/"//trim(adjustl(c0))//" of dtlong"//&
             &" (one grid "//trim(adjustl(c5))//" dt)"
     else
        write(*,"(a)", advance="NO") " reaches "//&
             &trim(adjustl(c3))//"/"//trim(adjustl(c0))//" of dtlong"//&
             &" ("//trim(adjustl(c2))//"/"//trim(adjustl(c4))//&
             &" of grid "//trim(adjustl(c5))//" dt)"
     end if

     ! find all direct sons

     ng=0
     do ifg = 1, ngrids
        if (nxtnest(ifg) == grid) then
           ng = ng + 1
           bcgrid(ng) = ifg
        end if
     end do

     ! dump grid(s) to receive boundary conditions

     if (isched(is,2) /= ng) then
        write(c0,"(i10)") is
        write(c1,"(i10)") isched(is,2)
        write(c2,"(i10)") ng
        write(*, "(/,a)") h//"**ERROR**: isched("//trim(adjustl(c0))//&
             &") = "//trim(adjustl(c1))//" and ng = "//trim(adjustl(c2))//&
             &" disagree"
        stop
     else if (isched(is,2) == 1) then
        write(c2,"(i10)") bcgrid(ng)
        write(*,"(a)", advance="NO") "; changes boundary of grid "//trim(adjustl(c2))
     else if (isched(is,2) > 1) then
        write(*,"(a)", advance="NO") "; changes boundary of grids "
        do ifg = 1, ng - 1
           write(c2,"(i10)") bcgrid(ifg)
           write(*,"(a)", advance="NO") trim(adjustl(c2))//", "
        end do
        write(c2,"(i10)") bcgrid(ng)
        write(*,"(a)", advance="NO") trim(adjustl(c2))
     end if

     ! dump grid(s) to receive feedback fields

     if (isched(is,4) /= 0) then
        ifg = isched(is,1)
        icg = nxtnest(ifg)
        write(c4,"(i10)") icg
        write(*,"(a)", advance="NO") "; feedbacks grid "//trim(adjustl(c4))
        do ng = 2, isched(is,4)
           ifg = icg
           icg = nxtnest(ifg)
           write(c4,"(i10)") icg
           write(*,"(a)", advance="NO") " that feedbacks grid "//trim(adjustl(c4))
        end do
     end if
     write(*,"(a)")
  end do
  write(*,"(a)")
end subroutine dump_modsched

!     *****************************************************************

subroutine cfl(n1,n2,n3,i0,j0,mynum)

  use mem_basic, only: &
       basic_g

  use mem_grid, only: &
       ngrid,         &
       grid_g

  implicit none

  integer :: n1,n2,n3,i0,j0,mynum

  call cfll(n1,n2,n3,i0,j0,mynum  &
       ,basic_g(ngrid)%up   (1,1,1)  ,basic_g(ngrid)%vp   (1,1,1)  &
       ,basic_g(ngrid)%wp   (1,1,1)  ,grid_g(ngrid)%rtgt    (1,1)  &
       ,grid_g(ngrid)%f13t    (1,1)  ,grid_g(ngrid)%f23t    (1,1)  &
       ,grid_g(ngrid)%dxt     (1,1)  ,grid_g(ngrid)%dyt     (1,1)  )
  return
end subroutine cfl

! **********************************************************************

subroutine cfll(n1,n2,n3,i0,j0,mynum,up,vp,wp,rtgt,f13t,f23t,dxt,dyt)

  use mem_grid, only: &
       nnxp,          &
       nnyp,          &
       cflxy,         &
       cflz,          &
       jdim,          &
       ngrids,        &
       ngrid,         &
       nxtnest,       &
       ipm,           &
       jpm,           &
       dtlt,          &
       ht,            &
       dzt

  use mem_scratch, only: &
       vctr1,            &
       vctr2,            &
       vctr3

  implicit none

  integer :: n1,n2,n3,i0,j0,mynum
  real, dimension(n1,n2,n3) :: up,vp,wp
  real, dimension(n2,n3)    :: rtgt,f13t,f23t,dxt,dyt

  integer :: i,j,k,ifm,icm,innest,nprints
  real :: c1x,c1y,c1z,cflnumh,cflnumv

  !     This routine returns flags the model to bring itself down when the CFL
  !     linear stability criteria on advection is exceeded.
  !     (Actually check on 90% of CFL)

5 format('cflx,ngrid,k,i,j,mynum = ',f5.1,5i5)
6 format('cfly,ngrid,k,i,j,mynum = ',f5.1,5i5)
7 format('cflz,ngrid,k,i,j,mynum = ',f5.1,5i5)

  nprints = 0
  cflnumh = .90
  cflnumv = .90
  cflxy(ngrid) = 0.
  cflz(ngrid) = 0.

  ! Let's try a new thing... if we have a grid point that is on a 
  !   coarse grid, but it is under a nested grid, we will ignore it
  !   under the assumption that the fine grid values will overwrite
  !   it, hence not allowing it to go numerically unstable.

  jloop: do j = 1+jdim,n3-jdim
     iloop: do i = 2,n2-1

        ! See if this is under a fine grid horizontally... ignore vertical for now
        innest=0
        if (ngrids > ngrid) then
           do ifm=ngrid+1,ngrids
              icm=nxtnest(ifm)
              if(icm == ngrid .and. &
                   i+i0 >= ipm(1,ifm) .and. i+i0 <= ipm(nnxp(ifm),ifm) .and. &
                   j+j0 >= jpm(1,ifm) .and. j+j0 <= jpm(nnyp(ifm),ifm) ) then
                 innest=1
                 exit
              endif
           enddo
        endif

        if(innest == 1) then
           !print '(a,5i4)', 'cfl under grid-' &
           !         ,ngrid,ifm,i+i0,j+j0,mynum
           cycle iloop
        endif

        kloop: do k = 2,n1-1

           vctr1(k) = .5*(up(k,i,j)+up(k,i-1,j))*dtlt*dxt(i,j)
           vctr2(k) = .5*(vp(k,i,j)+vp(k,i,j-jdim))*dtlt*dyt(i,j)
           vctr3(k) = ((wp(k,i,j)+wp(k-1,i,j))  &
                +(up(k,i,j)+up(k,i-1,j))*f13t(i,j)*ht(k)*rtgt(i,j)  &
                +(vp(k,i,j)+vp(k,i,j-jdim))*f23t(i,j)*ht(k)*rtgt(i,j)  &
                )*.5*dtlt*dzt(k)
        enddo kloop

        do k = 2,n1-1
           c1x = abs(vctr1(k))
           c1y = abs(vctr2(k))
           c1z = abs(vctr3(k))

           if (nprints .le. 10) then
              if (c1x .gt. cflnumh) then
                 nprints = nprints + 1
                 print 5, c1x,ngrid,k,i+i0,j+j0,mynum
                 print*,up(k,i,j),up(k,i-1,j),dtlt,dxt(i,j),vctr1(k)
              endif
              if (c1y .gt. cflnumh) then
                 nprints = nprints + 1
                 print 6, c1y,ngrid,k,i+i0,j+j0,mynum
              endif
              if (c1z .gt. cflnumv) then
                 nprints = nprints + 1
                 print 7, c1z,ngrid,k,i+i0,j+j0,mynum
              endif
           endif

           if (c1x .gt. cflxy(ngrid)) cflxy(ngrid) = c1x
           if (c1y .gt. cflxy(ngrid)) cflxy(ngrid) = c1y
           if (c1z .gt. cflz(ngrid)) cflz(ngrid) = c1z
        enddo
     enddo iloop
  enddo jloop

  return
end subroutine cfll

! ***************************************************************************

subroutine dump_dtset(nndtflg)
  use mem_grid, only: &
       ngrids,        &
       nndtrat,       &
       nnacoust,      &
       dtlongn

  implicit none
  integer, intent(in) :: nndtflg
  integer :: ifm

  write(*,"(a)")
  if (nndtflg == 0) then
     write(*,"(a)") " === Grids Delta T ===="
  else
     write(*,"(a)") " === Changed Grids Delta T ===="
  end if
  write(*,"(a)") " Grid, delta t, fraction of next coarser grid dt, acoustic steps per delta t" 
  do ifm = 1, ngrids
     write(*,"(i5,f9.2,15x,i4,24x,i4)") ifm, dtlongn(ifm), nndtrat(ifm), nnacoust(ifm)
  end do
end subroutine dump_dtset
