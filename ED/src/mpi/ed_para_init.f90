!------------------------------------------------------------------------------------------!
!    Subroutines based on the RAMS node decomposition. The main difference between the     !
! original code and this one is that when we split the domain we need to consider whether  !
! the polygon will fall on land or water. The water ones will be removed, so this should   !
! be taken into account for the standalone version.                                        !
!------------------------------------------------------------------------------------------!
subroutine ed_node_decomp(init,standalone,masterworks)

  use grid_coms,only: ngrids,nnxp,nnyp
  use ed_node_coms, only: mmxp,mmyp,mia,miz,mja,mjz,mi0,mj0,mibcon
  use ed_para_coms
  use mem_sites,only : n_ed_region
  use ed_work_vars,only : work_e,ed_alloc_work,ed_nullify_work
  use soil_coms, only: isoilflg

  implicit none

  integer , intent(in) :: init
  logical , intent(in) :: standalone,masterworks
  integer :: ngr,nsiz &
       ,ntotmachs

  !    This is a logical flag to test wheter the master should also do some
  ! processing

  if (masterworks) then
    ntotmachs=nmachs+1
  else
    ntotmachs=nmachs
  end if
    
  allocate(work_e(ngrids))

  !      Decompose all grids into subdomains
  
  do ngr = 1,ngrids
    ! SOI grids always have one point only. Since the structure will be sent
    ! I am filling the structures. It is just 4 extra numbers that will reach
    ! the other side anyway

     mmxp(ngr) = nnxp(ngr)
     mmyp(ngr) = nnyp(ngr)

     call ed_nullify_work(work_e(ngr))
     call ed_alloc_work(work_e(ngr),nnxp(ngr),nnyp(ngr))
  enddo

  call get_grid

  do ngr = 1,ngrids!n_ed_region

     ! Obtain estimates of the fraction of computational time (work) required
     ! for each column in the region of the domain.
     
     call get_work(ngr,mmxp(ngr),mmyp(ngr))
     
     ! Decompose the grid taking into account the work numbers.
     
     nsiz = nnxp(ngr)+nnyp(ngr)
     
     ! Here we decompose the domain into different working nodes, which may 
     ! include the master node.
     call ed_PAR_decomp(nnxp(ngr),nnyp(ngr),nsiz,ntotmachs,work_e(ngr)%work(1,1)  &
          ,ixb(1,ngr),ixe(1,ngr),iyb(1,ngr),iye(1,ngr))
     
  enddo

  do ngr=n_ed_region+1,ngrids
     call ed_newgrid(ngr)
     work_e(ngr)%work(1,1)=1.
     work_e(ngr)%land(1,1)=.true.
  end do
  
  ! Compute various bounds for the subdomains
  if (standalone) then
    call ed_PAR_decomp_bounds(n_ed_region,ngrids,nnxp,nnyp,0,0,masterworks)
  else
    call ed_PAR_decomp_bounds(n_ed_region,ngrids,nnxp,nnyp,1,1,masterworks)
  end if
  
  ! If your master node was set to do some work, then it's the time to fill the master node 
  ! with the specific node information.
  if (masterworks) then 
     do ngr=1,ngrids
        
        mmxp(ngr)=nxend(ntotmachs,ngr)-nxbeg(ntotmachs,ngr)+1
        mmyp(ngr)=nyend(ntotmachs,ngr)-nybeg(ntotmachs,ngr)+1
        mia(ngr)=nxbegc(ntotmachs,ngr)
        miz(ngr)=nxendc(ntotmachs,ngr)
        mja(ngr)=nybegc(ntotmachs,ngr)
        mjz(ngr)=nyendc(ntotmachs,ngr)
        mi0(ngr)=ixoff(ntotmachs,ngr)
        mj0(ngr)=iyoff(ntotmachs,ngr)
        mibcon(ngr)=ibcflg(ntotmachs,ngr)
     end do
  end if

  return
end subroutine ed_node_decomp

!==========================================================================================!
!==========================================================================================!

subroutine get_grid

  use mem_sites, only: grid_type,grid_res,n_ed_region &
                      ,soi_lat,soi_lon,ed_reg_lonmin,ed_reg_latmin

  use grid_coms, only: ngrids,nnxp,nnyp,nstratx,nstraty
  use ed_work_vars, only: work_e
  implicit none

  integer :: ifm,i,j

  select case (grid_type)
  case (0)          ! lat-lon type grid
     do ifm=1,n_ed_region
        do i=1,nnxp(ifm)
           do j=1,nnyp(ifm)
              work_e(ifm)%glon(i,j) = ed_reg_lonmin(ifm) + (float(i) - 0.5) * grid_res/real(nstratx(ifm))
              work_e(ifm)%glat(i,j) = ed_reg_latmin(ifm) + (float(j) - 0.5) * grid_res/real(nstraty(ifm))
           end do
        end do
     end do
  case (1) ! polar-stereo type grid
     call ed_gridset(1)
     do ifm=1,n_ed_region
        call ed_newgrid(ifm)
        call ed_polarst(nnxp(ifm),nnyp(ifm),work_e(ifm)%glat(1,1),work_e(ifm)%glon(1,1))
     end do
     
  case default
     !  Eventually we'll have polygons, but not yet.
     call fatal_error('Invalid grid_type in ED_grid_setup.' &
          ,'get_grid','ed_para_init.f90')
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

  use ed_work_vars,only : work_e
  use soil_coms, only: veg_database,soil_database,isoilflg,nslcon
  use mem_sites,only:n_soi

  implicit none
  integer, intent(in) :: ifm
  integer :: npoly
  integer, intent(in) :: nxp,nyp

  real, allocatable, dimension(:,:) :: lat_list
  real, allocatable, dimension(:,:) :: lon_list
  integer, allocatable, dimension(:) :: leaf_class_list
  integer, allocatable, dimension(:) :: ntext_soil_list
  integer, allocatable,dimension(:) :: ipcent_land
  integer :: datsoil,ipy,i,j
  integer :: jboff,jtoff,iloff,iroff
  integer,parameter :: min_land_pcent = 10
  real,   parameter :: soi_edge_deg = 0.05   ! 100th of a degree, about 5.5 km at the equator.

  npoly = nxp*nyp

  allocate(lat_list(3,npoly))
  allocate(lon_list(3,npoly))
  allocate(leaf_class_list(npoly))
  allocate(ipcent_land(npoly))
  
  ! Fill lat/lon lists

  ! j index is the North-South index, and it is inverted, ie larger index larger latitude

  if (n_soi.gt.0 .and. ifm.le.n_soi) then

     ipy = 0
     do i=1,nxp
        do j = 1,nyp
           ipy = ipy + 1
           
           lat_list(1,ipy) = work_e(ifm)%glat(i,j)
           lon_list(1,ipy) = work_e(ifm)%glon(i,j)
           
           ! Top latitude
           lat_list(2,ipy) = work_e(ifm)%glat(i,j) + soi_edge_deg
           
           ! Bottom latitude
           lat_list(3,ipy) = work_e(ifm)%glat(i,j) - soi_edge_deg
           
           ! Left longitude
           lon_list(2,ipy) = work_e(ifm)%glon(i,j) - soi_edge_deg
           
           ! Right longitude
           lon_list(3,ipy) = work_e(ifm)%glon(i,j) + soi_edge_deg
           
           if (lon_list(1,ipy) >=  180.) lon_list(1,ipy) = lon_list(1,ipy) - 360.
           if (lon_list(1,ipy) <= -180.) lon_list(1,ipy) = lon_list(1,ipy) + 360.
           if (lon_list(2,ipy) >=  180.) lon_list(2,ipy) = lon_list(2,ipy) - 360.
           if (lon_list(2,ipy) <= -180.) lon_list(2,ipy) = lon_list(2,ipy) + 360.
           if (lon_list(3,ipy) >=  180.) lon_list(3,ipy) = lon_list(3,ipy) - 360.
           if (lon_list(3,ipy) <= -180.) lon_list(3,ipy) = lon_list(3,ipy) + 360.

        enddo
     enddo

  else
     
     ipy = 0
     do i=1,nxp
        do j = 1,nyp
           ipy = ipy + 1
           
           iloff=1
           iroff=1
           jtoff=1
           jboff=1
           
           lat_list(1,ipy) = work_e(ifm)%glat(i,j)
           lon_list(1,ipy) = work_e(ifm)%glon(i,j)
           
           ! Top latitude
           if(j==nyp)jtoff=-1
           lat_list(2,ipy) = work_e(ifm)%glat(i,j) + real(jtoff)*0.5*(work_e(ifm)%glat(i,j+jtoff)-work_e(ifm)%glat(i,j))
           
           ! Bottom latitude
           if(j==1)jboff=-1
           lat_list(3,ipy) = work_e(ifm)%glat(i,j) + real(jboff)*0.5*(work_e(ifm)%glat(i,j-jboff)-work_e(ifm)%glat(i,j))
           
           ! Left longitude
           if(i==1)iloff=-1
           lon_list(2,ipy) = work_e(ifm)%glon(i,j) + real(iloff)*0.5*(work_e(ifm)%glon(i-iloff,j)-work_e(ifm)%glon(i,j))
           
           ! Right longitude
           if(i==nxp)iroff=-1
           lon_list(3,ipy) = work_e(ifm)%glon(i,j) + real(iroff)*0.5*(work_e(ifm)%glon(i+iroff,j)-work_e(ifm)%glon(i,j))
           
           if (lon_list(1,ipy) >=  180.) lon_list(1,ipy) = lon_list(1,ipy) - 360.
           if (lon_list(1,ipy) <= -180.) lon_list(1,ipy) = lon_list(1,ipy) + 360.
           if (lon_list(2,ipy) >=  180.) lon_list(2,ipy) = lon_list(2,ipy) - 360.
           if (lon_list(2,ipy) <= -180.) lon_list(2,ipy) = lon_list(2,ipy) + 360.
           if (lon_list(3,ipy) >=  180.) lon_list(3,ipy) = lon_list(3,ipy) - 360.
           if (lon_list(3,ipy) <= -180.) lon_list(3,ipy) = lon_list(3,ipy) + 360.
           
        enddo
     enddo

  endif
     

  ! Generate the land/sea mask

  write(unit=*,fmt=*) ' => Generating the land/sea mask.'

  call leaf_database(trim(veg_database), npoly, 'leaf_class', lat_list,  &
       lon_list, ipcent_land)

  if (isoilflg(ifm) == 1) then
     allocate(ntext_soil_list(npoly))
     call leaf_database(trim(soil_database), npoly, 'soil_text', lat_list,  &
          lon_list, ntext_soil_list)
  end if
 
  ! Re-map the land cover classes

  ipy = 0
  do i=1,nxp
     do j = 1,nyp
        ipy = ipy + 1

        work_e(ifm)%land(i,j) = ipcent_land(ipy) > min_land_pcent

        if (work_e(ifm)%land(i,j)) then
           work_e(ifm)%work(i,j)      = 1.0
           work_e(ifm)%landfrac(i,j)  = real(ipcent_land(ipy))/100.0

           if (isoilflg(ifm) == 1) then
              datsoil = ntext_soil_list(ipy)

              ! This is to prevent datsoil to be zero when the polygon was assumed land
              if (datsoil == 0) datsoil=nslcon
              work_e(ifm)%ntext(i,j) = datsoil
           else  !! set from ED2IN
              work_e(ifm)%ntext(i,j) = nslcon
           end if
        else
           !----- Making this grid point 100% water ---------------------------------------!
           work_e(ifm)%landfrac(i,j)  = 0.
           work_e(ifm)%work(i,j)      = epsilon(0.0)
           work_e(ifm)%ntext(i,j)     = 0
        end if
     end do
  end do

  ! PRINT OUT THE ARRAYS !
!  do j = nyp,1,-1
!     print*,(work_e(ifm)%glat(i,j),i=1,nxp)
!  enddo

!  do j = nyp,1,-1
!     print*,(work_e(ifm)%glon(i,j),i=1,nxp)
!  enddo

!  do j = nyp,1,-1
!     print*,(work_e(ifm)%landfrac(i,j),i=1,nxp)
!  enddo
 
  
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
subroutine ed_PAR_decomp(nxp,nyp,nsiz,nodes,work,ixb,ixe,iyb,iye)
  implicit none
  integer :: nxp,nyp,nsiz,nodes
  real :: work(nxp,nyp),workrow(nsiz),workload(nsiz),workblock(nsiz),workcol(nsiz)
  real :: anblocks(nsiz),ajrows(nsiz),ajrow(nsiz)
  integer :: ixb(*),ixe(*),iyb(*),iye(*)
  integer :: nblocks(nsiz),jrows(nsiz),jrow(nsiz)
  ! default relspeed = 1.0 for nodes of uniform speed.
  real :: relspeed(256)
  data relspeed/256*1./

  integer :: inode,i,j,islab,nslabs,min_blocks,nbigslabs,iblock &
       ,jnode,knode
  real :: anodes,aslabs,totspeed,workdom,workaccum,worksofar &
       ,slabspeed,workslab
      

  ! This routine decomposes grid domains of size (nnxp,nnyp) into a number,
  ! specified by nodes, of rectangular subdomains.  The convention is followed
  ! that any internal boundaries (between subdomains) that are parallel to
  ! the x-axis run continuously across the full domain, while boundaries
  ! parallel to the y-axis may or may not run the full distance across the
  ! domain.  For convenience, regions of the domain bounded by adjacent
  ! east-west internal boundaries are termed "slabs", while smaller divisions
  ! within each slab are termed "blocks".  Each block is required to have
  ! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
  ! with the given input parameters, the subroutine stops.


  ! Estimate the number of slabs to be used (aslabs), and compute a final
  ! nearest integer value (nslabs) which is limited to allowable values.
  ! Zero out array for accumulating number of columns for each node.

  if (nodes == 1) then
    ixb(1) = 1
    ixe(1) = nxp
    iyb(1) = 1
    iye(1) = nyp
    return
  end if


  anodes=float(nodes)
  aslabs = sqrt(anodes * float(nyp) / float(nxp))
  nslabs = min(nodes,max(1,nint(aslabs)))

  totspeed = 0.
  do inode = 1,nodes
     ixe(inode) = 0
     totspeed = totspeed + relspeed(inode)
  enddo

  !          print*, 'totspeed',totspeed

  ! Compute total work load over each row and over entire domain.

  workdom = 0.
  do j = 1,nyp
     workrow(j) = 0.
     do i = 1,nxp
        workrow(j) = workrow(j) + work(i,j)
     enddo
     workdom = workdom + workrow(j)

     !          print*, 'j,workdom,workrow(j)',j,workdom,workrow(j)

  enddo
!  workrow(2) = workrow(2) + workrow(1)
!  workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)

  ! Determine number of blocks and the average workload for each slab.

  min_blocks = nodes / nslabs
  nbigslabs = nodes - min_blocks * nslabs
  inode = 0
  do islab = 1,nslabs
     workload(islab) = 0.
     nblocks(islab) = min_blocks
     if (islab .le. nbigslabs) nblocks(islab) = min_blocks + 1
     do iblock = 1,nblocks(islab)
        inode = inode + 1
        workload(islab) = workload(islab)  &
             + workdom * relspeed(inode) / totspeed

        !           print*, 'islab,iblock,workload(islab),workdom,inode'
        !           print*,  islab,iblock,workload(islab),workdom,inode

     enddo
  enddo

  ! Assign all j-rows to their respective slabs in a way that balances the work
  ! load among slabs according to their respective numbers of nodes (blocks).
  ! The array jrows counts the number of rows in each slab, and the array
  ! jrow is the index of the southernmost row in each slab.

  do islab = 1,nslabs
     jrows(islab) = 0
  enddo

  workaccum = 0.
  worksofar = 0.
  islab = 0

  do j = 1,nyp
     workaccum = workaccum + workrow(j)
     if (workaccum - .5 * workrow(j) .gt. worksofar .and.  &
          islab .lt. nslabs) then
        islab = islab + 1
        jrow(islab) = j
        worksofar = worksofar + workload(islab)
     endif
     jrows(islab) = jrows(islab) + 1
  enddo

  inode = 0
  jnode = 0
  knode = 0
  do islab = 1,nslabs

     ! Compute the total work load for each slab and for each i-column in the
     ! slab.

     slabspeed = 0.
     workslab = 0.
     do i = 1,nxp
        workcol(i) = 0.
        do j = jrow(islab),jrow(islab)+jrows(islab)-1
           workcol(i) = workcol(i) + work(i,j)
        enddo
        workslab = workslab + workcol(i)
     enddo
!     workcol(2) = workcol(2) + workcol(1)
!     workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)

     ! Determine average workload for each block.

     do iblock = 1,nblocks(islab)
        jnode = jnode + 1
        slabspeed = slabspeed + relspeed(jnode)

        !           print*, 'r1:iblock,jnode,slabspeed,relspeed(jnode)'
        !           print*,     iblock,jnode,slabspeed,relspeed(jnode)

     enddo
     do iblock = 1,nblocks(islab)
        knode = knode + 1
        workblock(iblock) = workslab  &
             * relspeed(knode) / slabspeed

        !       print*, 'islab,iblock,workblock,workslab,relspeed,slabspeed'
        !       print*, islab,iblock,workblock(iblock),workslab,relspeed(knode)
        !     +       ,slabspeed
        !       print*, 'knode',knode

     enddo

     ! Assign the i-columns of each slab to their respective blocks in a way that
     ! balances the work load among the blocks.  The array ncols counts the number
     ! of i-columns on each node, and the array ncol is the index of the
     ! westernmost i-column on each node.

     workaccum = 0.
     worksofar = 0.

     iblock = 0
     do i = 1,nxp
        workaccum = workaccum + workcol(i)

        !        print*, 'islab',islab
        !        print*, 'i,workaccum,workcol(i),worksofar,iblock,nblocks'
        !        print*, i,workaccum,workcol(i),worksofar,iblock,nblocks(islab)

        if (workaccum - .5 * workcol(i) .gt. worksofar .and.  &
             iblock .lt. nblocks(islab)) then
           iblock = iblock + 1

           !ccccccc defining node variables here ccccccccccccccccccccccccccccccccccc
           inode = inode + 1
           iyb(inode) = jrow(islab)
           ixb(inode) = i
           iye(inode) = iyb(inode) + jrows(islab) - 1

           !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           worksofar = worksofar + workblock(iblock)
        endif
        ixe(inode) = ixe(inode) + 1
     enddo
  enddo

  !ccccccc defining node variable here ccccccccccccccccccccccccccccccccccc
  do jnode = 1,nodes
     ixe(jnode) = ixb(jnode) + ixe(jnode) - 1

     !           print*, 'jnode,ixb,ixe',jnode,ixb(jnode),ixe(jnode)
     !     +        ,(ixe(jnode)-ixb(jnode)+1),(iye(jnode)-iyb(jnode)+1)
     !     +        ,(ixe(jnode)-ixb(jnode)+1)*(iye(jnode)-iyb(jnode)+1)

  enddo
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Check to make sure that each subdomain has at least 2 interior 
  ! rows and columns.

  do jnode = 1,nodes
     if ((iye(jnode) - iyb(jnode)+1)*(ixe(jnode) - ixb(jnode)+1) .lt. 1) then
        print*, 'grid:',nxp,nyp,'  subdomain too small on node ',jnode
        print*, '(ixb,ixe,iyb,iye) = '  &
             ,ixb(jnode),ixe(jnode),iyb(jnode),iye(jnode)
        call fatal_error('Region is too small','ed_PAR_decomp','ed_para_init.f90')
     endif
  enddo

![MLO - Saving the integers into the real arrays
  do i=1,nsiz
    anblocks(i)=float(nblocks(i))
    ajrows(i)=float(jrows(i))
    ajrow(i)=float(jrow(i))
  end do


  return
end subroutine ed_PAR_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_PAR_decomp_bounds(n_ed_region,ngrids,nnxp,nnyp,nbndx,nbndy,masterworks)

   use ed_para_coms

   implicit none
   integer, intent(in) :: n_ed_region,ngrids
   logical, intent(in) :: masterworks

   integer :: nnxp(*),nnyp(*),nbndx,nbndy
   integer :: ng,nm,nx,ny,ntotmachs

   !        Compute various subdomain boundary numbers for the nodes
   !             nxbeg,nybeg,nxend,nyend - portions of full domain that nodes will have
   !                                     - includes overlap region
   !             ixoff,iyoff  - subdomain offsets relative to full domain
   !             nxbegc,nxendc,nybegc,nyendc - subdomain "compute" points,
   !                         or normal thermodynamic tendency points (2-nx, 2-ny for
   !                          non-parallel run
   !             ibcflag - flag denoting if real boundary is on subdomain
   !                       bit 1=west, bit 2=east, bit 4=south, bit 8=north

   if (masterworks) then
      ntotmachs=nmachs+1
   else 
      ntotmachs=nmachs
   end if


   do nm=1,ntotmachs
      do ng=1,n_ed_region
         nx=nnxp(ng)
         ny=nnyp(ng)
         ibcflg(nm,ng)=0
         !----- West Boundary
         if(ixb(nm,ng) == (nbndx+1)) then
            nxbeg(nm,ng)=1
            ixoff(nm,ng)=0
            nxbegc(nm,ng)=nbndx+1
            ibcflg(nm,ng)=ibcflg(nm,ng)+1
         else
            nxbeg(nm,ng)=ixb(nm,ng)-nbndx
            ixoff(nm,ng)=nxbeg(nm,ng)-1
            nxbegc(nm,ng)=nbndx+1
         endif

         !----- East Boundary
         if(ixe(nm,ng) == nx-nbndx) then
            nxend(nm,ng)=nx
            ibcflg(nm,ng)=ibcflg(nm,ng)+2
            nxendc(nm,ng)=(ixe(nm,ng)-ixb(nm,ng))+nxbegc(nm,ng)
         else
            nxend(nm,ng)=ixe(nm,ng)+nbndx
            nxendc(nm,ng)=(ixe(nm,ng)-ixb(nm,ng))+nxbegc(nm,ng)
         endif

         !----- South Boundary
         if(iyb(nm,ng) == (nbndy+1)) then
            nybeg(nm,ng)=1
            iyoff(nm,ng)=0
            nybegc(nm,ng)=nbndy+1
            ibcflg(nm,ng)=ibcflg(nm,ng)+4
         else
            nybeg(nm,ng)=iyb(nm,ng)-nbndy
            iyoff(nm,ng)=nybeg(nm,ng)-1
            nybegc(nm,ng)=nbndy+1
         endif

         !----- North Boundary
         if(iye(nm,ng) == (ny-nbndy)) then
            nyend(nm,ng)=ny
            ibcflg(nm,ng)=ibcflg(nm,ng)+8
            nyendc(nm,ng)=(iye(nm,ng)-iyb(nm,ng))+nybegc(nm,ng)
         else
            nyend(nm,ng)=iye(nm,ng)+nbndy
            nyendc(nm,ng)=(iye(nm,ng)-iyb(nm,ng))+nybegc(nm,ng)
         endif

         npxy(nm,ng)=(nxend(nm,ng)-nxbeg(nm,ng)+1)  &
              *(nyend(nm,ng)-nybeg(nm,ng)+1)
      end do
      ! These are the SOI runs, everything should be set to 1 (grid with a single polygon).
      do ng=n_ed_region+1,ngrids
         nxbeg  (nm,ng) = 1
         nxbegc (nm,ng) = 1
         nybeg  (nm,ng) = 1
         nybegc (nm,ng) = 1
         nxend  (nm,ng) = 1
         nxendc (nm,ng) = 1
         nyend  (nm,ng) = 1
         nyendc (nm,ng) = 1
         ixoff  (nm,ng) = 0
         iyoff  (nm,ng) = 0
         ibcflg (nm,ng) = 15
         npxy   (nm,ng) = 1
      end do
   end do

   return
end subroutine ed_PAR_decomp_bounds
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_dump_Domain_Decomposition(masterworks)
  use grid_coms, only: &
       ngrids

  use ed_para_coms, only: &
              nmachs,     &
              ixb,        &
              ixe,        &
              iyb,        &
              iye,        &
              ixoff,      &
              iyoff
  use ed_state_vars, only: gdpy,py_off

  implicit none
  logical, intent(in) :: masterworks
  integer :: ngr
  integer :: jnode
  
  
  write(*,'(a)') '================================= Domain decomposition ==================================='
  write(*,'(10(a8,1x))') adjustr('grid'),adjustr('node'),adjustr('x-beg'),adjustr('x-end') &
                       ,adjustr('y-beg'),adjustr('y-end'),adjustr('x-off'),adjustr('y-off') &
                       ,adjustr('npolys'),adjustr('offset')
  do ngr = 1, ngrids
     do jnode = 1,nmachs
        write(unit=*,fmt='(10(i8,1x))') ngr,jnode,ixb(jnode,ngr),ixe(jnode,ngr),iyb(jnode,ngr)  &
                                      ,iye(jnode,ngr),ixoff(jnode,ngr),iyoff(jnode,ngr)         &
                                      ,gdpy(jnode,ngr),py_off(jnode,ngr)
     enddo
     if (masterworks) then
        jnode=nmachs+1
        write(unit=*,fmt='(i8,1x,a8,1x,8(i8,1x))') ngr,adjustr('master'),ixb(jnode,ngr)      &
                                      ,ixe(jnode,ngr),iyb(jnode,ngr),iye(jnode,ngr)          &
                                      ,ixoff(jnode,ngr),iyoff(jnode,ngr)                     &
                                      ,gdpy(jnode,ngr),py_off(jnode,ngr)
     end if
  end do
  write(*,'(a)') '=========================================================================================='
  write(*,'(a)')
end subroutine ed_dump_Domain_Decomposition
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_PAR_est_time(nxp,nyp,work,cput,init,standalone)
  implicit none
  logical, intent(in) :: standalone
  integer :: nxp,nyp,init
  real :: work(nxp,nyp),cput(nxp,nyp)

  integer :: i,j
  real :: bfact

  ! Sample routine to fill work elements with values proportional to the time
  ! required to perform model operations.


  if(init==1) then
     do j = 2,nyp-1
        do i = 2,nxp-1
           work(i,j) = 1.
        enddo
     enddo
  else
     do j = 2,nyp-1
        do i = 2,nxp-1
           work(i,j) = cput(i,j)
        enddo
     enddo
  endif


  if (standalone) then
     bfact=1.0
  else
     ! Fill real boundaries with .2 of interior points, but just for the case in which
     ! we are running coupled runs. Standalone version should be all the same.
     bfact=.2
  end if
  do j = 1,nyp
     work(1,j) = bfact * work(2,j)
     work(nxp,j) = bfact * work(nxp-1,j)
  enddo

  do i = 1,nxp
     work(i,1) = bfact * work(i,2)
     work(i,nyp) = bfact * work(i,nyp-1)
  enddo

  return
end subroutine ed_PAR_est_time
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine dump_gridwork(ifm,xmax,ymax,x0,y0,glon,glat,work,text,land)
   use ed_node_coms, only: mynum
   implicit none
   integer, intent(in)                       :: ifm,xmax,ymax,x0,y0
   real   , intent(in), dimension(xmax,ymax) :: glon,glat,work
   integer, intent(in), dimension(xmax,ymax) :: text
   logical, intent(in), dimension(xmax,ymax) :: land
   integer           :: unit,x,y,xa,ya,xz,yz
   character(len=14) :: fmt_real,fmt_int
   character(len=12) :: fmt_log
   character(len=33) :: dumpname

   unit=60+mynum
   xa=x0+1
   ya=y0+1
   xz=x0+xmax
   yz=y0+ymax
   write(fmt_real,fmt='(a,i3.3,a)') '(',xmax,'(f8.3,1x))'
   write(fmt_log ,fmt='(a,i3.3,a)') '(',xmax,'(l1,1x))'
   write(fmt_int ,fmt='(a,i3.3,a)') '(',xmax,'(i2.2,1x))'
   
   write(dumpname,fmt='(a,i3.3,a,i3.3,a)') './grid_dump_node-',mynum,'_grid-',ifm,'.txt'
   
   open(unit=unit,file=dumpname,status='replace',action='write')
   write(unit=unit,fmt='(6(a,1x,i4,1x))') 'xmax=',xmax,'xa=',xa,'xz=',xz,'ymax=',ymax,'ya=',ya,'yz=',yz
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') 'Longitude:'
   do y=ymax,1,-1
         write(unit=unit,fmt=fmt_real) (glon(x,y),x=1,xmax)
   end do
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') 'Latitude:'
   do y=ymax,1,-1
         write(unit=unit,fmt=fmt_real) (glat(x,y),x=1,xmax)
   end do
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') 'Work:'
   do y=ymax,1,-1
         write(unit=unit,fmt=fmt_real) (work(x,y),x=1,xmax)
   end do
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') 'Soil_texture:'
   do y=ymax,1,-1
         write(unit=unit,fmt=fmt_int) (text(x,y),x=1,xmax)
   end do
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') 'Land:'
   do y=ymax,1,-1
         write(unit=unit,fmt=fmt_log) (land(x,y),x=1,xmax)
   end do
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   write(unit=unit,fmt='(a)') '=================================================================================='
   close(unit=unit,status='keep')
   return
end subroutine dump_gridwork
