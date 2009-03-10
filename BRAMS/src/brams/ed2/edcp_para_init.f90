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
subroutine edcp_get_work(ifm,mxp,myp,iwest,ieast,jsouth,jnorth)

  use ed_work_vars,only : work_e
  use soil_coms, only: veg_database,soil_database,nslcon
  use io_params, only: b_isoilflg => isoilflg & ! intent(in)
                     , b_ivegtflg => ivegtflg ! ! intent(in)
  use mem_sites,only:n_soi

  implicit none
  integer, intent(in) :: ifm,iwest,ieast,jsouth,jnorth
  integer :: npoly
  integer, intent(in) :: mxp,myp

  real, allocatable, dimension(:,:) :: lat_list
  real, allocatable, dimension(:,:) :: lon_list
  integer, allocatable, dimension(:) :: leaf_class_list
  integer, allocatable, dimension(:) :: ntext_soil_list
  integer, allocatable,dimension(:) :: ipcent_land
  integer :: datsoil,ipy,i,j
  integer :: jboff,jtoff,iloff,iroff
  integer,parameter :: min_land_pcent = 25
  real,   parameter :: soi_edge_deg = 0.05   ! 100th of a degree, about 5.5 km at the equator.

  npoly = mxp*myp

  allocate(lat_list(3,npoly))
  allocate(lon_list(3,npoly))
  allocate(leaf_class_list(npoly))
  allocate(ipcent_land(npoly))
  
  ! Fill lat/lon lists

  ! j index is the North-South index, and it is inverted, ie larger index larger latitude

  if (n_soi.gt.0 .and. ifm.le.n_soi) then

     ipy = 0
     do i=1,mxp
        do j = 1,myp
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
     do i=1,mxp
        do j = 1,myp
           ipy = ipy + 1
           
           iloff=1
           iroff=1
           jtoff=1
           jboff=1
           
           lat_list(1,ipy) = work_e(ifm)%glat(i,j)
           lon_list(1,ipy) = work_e(ifm)%glon(i,j)
           
           ! Top latitude
           if(j==myp)jtoff=-1
           lat_list(2,ipy) = work_e(ifm)%glat(i,j) + real(jtoff)*0.5*(work_e(ifm)%glat(i,j+jtoff)-work_e(ifm)%glat(i,j))
           
           ! Bottom latitude
           if(j==1)jboff=-1
           lat_list(3,ipy) = work_e(ifm)%glat(i,j) + real(jboff)*0.5*(work_e(ifm)%glat(i,j-jboff)-work_e(ifm)%glat(i,j))
           
           ! Left longitude
           if(i==1)iloff=-1
           lon_list(2,ipy) = work_e(ifm)%glon(i,j) + real(iloff)*0.5*(work_e(ifm)%glon(i-iloff,j)-work_e(ifm)%glon(i,j))
           
           ! Right longitude
           if(i==mxp)iroff=-1
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
  write(unit=*,fmt='(3(a,1x,i5,1x))') 'GRID:',ifm,'IVEGTFLG',b_ivegtflg(ifm)               &
                                     ,'ISOILFLG:',b_isoilflg(ifm)
  select case (b_ivegtflg(ifm))
  case (0,1,2)
     call leaf3init_overwrite(iwest,ieast,jsouth,jnorth,npoly,ifm,'leaf',ipcent_land)
  case (3)
     call leaf_database(trim(veg_database), npoly, 'leaf_class', lat_list,  &
                        lon_list, ipcent_land)
  end select

  select case (b_isoilflg(ifm))
  !----------------------------------------------------------------------------------------!
  !    This allows us to use the vegetation and soil type defined for LEAF in coupled      !
  ! runs. Note that to use the ED dataset in coupled runs we need to use isoilflg =3.      !
  !----------------------------------------------------------------------------------------!
  case (0,1,2)
     allocate(ntext_soil_list(npoly))
     call leaf3init_overwrite(iwest,ieast,jsouth,jnorth,npoly,ifm,'soil',ntext_soil_list)
  case (3)
     allocate(ntext_soil_list(npoly))
     call leaf_database(trim(soil_database), npoly, 'soil_text', lat_list,  &
          lon_list, ntext_soil_list)
  end select
 
  ! Re-map the land cover classes

  ipy = 0
  do i=1,mxp
     do j = 1,myp
        ipy = ipy + 1

        work_e(ifm)%land(i,j)      = ipcent_land(ipy) > min_land_pcent

        if (work_e(ifm)%land(i,j)) then
           work_e(ifm)%work(i,j)      = 1.0
           work_e(ifm)%landfrac(i,j)  = real(ipcent_land(ipy))/100.0

           if (b_isoilflg(ifm) == 0) then !! set from ED2IN/RAMSIN
              work_e(ifm)%ntext(i,j) = nslcon
           else  !! set from data base or LEAF-3
              datsoil = ntext_soil_list(ipy)

              ! This is to prevent datsoil to be zero when the polygon was assumed land
              if (datsoil == 0) datsoil=nslcon
              work_e(ifm)%ntext(i,j) = datsoil
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
  !  do j = myp,1,-1
  !     print*,(work_e(ifm)%glat(i,j),i=1,mxp)
  !  enddo
  
  !  do j = myp,1,-1
  !     print*,(work_e(ifm)%glon(i,j),i=1,mxp)
  !  enddo
  
  !  do j = myp,1,-1
  !     print*,(work_e(ifm)%landfrac(i,j),i=1,mxp)
  !  enddo
  
  
  deallocate(lat_list)
  deallocate(lon_list)
  deallocate(leaf_class_list)
  deallocate(ipcent_land)
  if (allocated(ntext_soil_list)) deallocate (ntext_soil_list)

  return
end subroutine edcp_get_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine will copy LEAF-3 initial condition to ED, so we can use other LEAF     !
! databases to decide whether a polygon is inland or offshore, and aslo other soil texture !
! dataset.                                                                                 !
!------------------------------------------------------------------------------------------!
subroutine leaf3init_overwrite(iwest,ieast,jsouth,jnorth,nlandsea,ifm,varname,varout)
   use mem_leaf, only: leaf_g
   use mem_grid, only: nzg
   use node_mod, only: mynum
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)  :: iwest,ieast,jsouth,jnorth,nlandsea,ifm
   character(len=4)            , intent(in)  :: varname
   integer, dimension(nlandsea), intent(out) :: varout
   !----- Local variables -----------------------------------------------------------------!
   integer                                   :: i,j,ipy
   !---------------------------------------------------------------------------------------!

   select case (varname)
   case ('leaf')
      ipy=0
      do i=iwest,ieast
         do j=jsouth,jnorth
            ipy=ipy+1
            varout(ipy) = nint(100.* leaf_g(ifm)%patch_area(i,j,2))
         end do
      end do
   case ('soil')
      ipy=0
      do i=iwest,ieast
         do j=jsouth,jnorth
            ipy=ipy+1
            varout(ipy) = nint(leaf_g(ifm)%soil_text(nzg,i,j,2))
         end do
      end do
   case default
      call abort_run('Invalid key: '//varname,'leaf3init_overwrite.f90','edcp_met_init.f90')
   end select
   return
end subroutine leaf3init_overwrite
!==========================================================================================!
!==========================================================================================!
