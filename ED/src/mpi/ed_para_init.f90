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
     
     call ed_parvec_work(ngr,mmxp(ngr),mmyp(ngr),work_e(ngr)%land)
     
  enddo

  do ngr=n_ed_region+1,ngrids
     call ed_newgrid(ngr)
     work_e(ngr)%work(1,1)=1.
     work_e(ngr)%land(1,1)=.true.

     call ed_parvec_work(ngr,mmxp(ngr),mmyp(ngr),work_e(ngr)%land)

  end do


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
        call ed_polarst(nnxp(ifm),nnyp(ifm),work_e(ifm)%glat,work_e(ifm)%glon)
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
  use mem_sites,only:n_soi,grid_res,grid_type

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
  integer,parameter :: min_land_pcent = 25
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
           if(grid_type == 0)then
              lat_list(2,ipy) = work_e(ifm)%glat(i,j) + 0.5 * grid_res
              lat_list(3,ipy) = work_e(ifm)%glat(i,j) - 0.5 * grid_res
              lon_list(2,ipy) = work_e(ifm)%glon(i,j) - 0.5 * grid_res
              lon_list(3,ipy) = work_e(ifm)%glon(i,j) + 0.5 * grid_res
           elseif(grid_type == 1)then
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
           endif

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

        work_e(ifm)%land(i,j)      = ipcent_land(ipy) > min_land_pcent

        if (work_e(ifm)%land(i,j)) then
           work_e(ifm)%work(i,j)      = 1.0
           work_e(ifm)%landfrac(i,j)  = real(ipcent_land(ipy))/100.0

           select case (isoilflg(ifm))
           case (1)  !! set from data base or LEAF-3
              datsoil = ntext_soil_list(ipy)

              ! This is to prevent datsoil to be zero when the polygon was assumed land
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

subroutine ed_parvec_work(ifm,nxp,nyp,land)

  use ed_work_vars,  only: work_e
  
  implicit none

  integer,intent(in) :: nxp,nyp
  logical,intent(in),dimension(nxp,nyp) :: land
  integer :: npolygons
  integer :: poly
  integer,intent(in) :: ifm
  integer :: i,j
  
  ! Compute total work load over each row and over entire domain.
  
  npolygons = 0
  do j = 1,nyp
     do i = 1,nxp
        if(land(i,j)) then
           npolygons = npolygons + 1
        endif
     enddo
  enddo
  
  ! Allocate the polygon vectors
  
  allocate(work_e(ifm)%vec_glon(npolygons))
  allocate(work_e(ifm)%vec_glat(npolygons))
  allocate(work_e(ifm)%vec_landfrac(npolygons))
  allocate(work_e(ifm)%vec_ntext(npolygons))
  allocate(work_e(ifm)%vec_xid(npolygons))
  allocate(work_e(ifm)%vec_yid(npolygons))

  poly = 0
  do j = 1,nyp
     do i = 1,nxp
        
        if(land(i,j)) then
           poly = poly + 1
           
           work_e(ifm)%vec_glon(poly) = work_e(ifm)%glon(i,j)
           work_e(ifm)%vec_glat(poly) = work_e(ifm)%glat(i,j)
           work_e(ifm)%vec_landfrac(poly) = work_e(ifm)%landfrac(i,j)
           work_e(ifm)%vec_ntext(poly) = work_e(ifm)%ntext(i,j)
           work_e(ifm)%vec_xid(poly) = i
           work_e(ifm)%vec_yid(poly) = j

        endif
     end do
  end do
  
  return
end subroutine ed_parvec_work


