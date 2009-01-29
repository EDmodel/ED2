
subroutine set_polygon_coordinates_ar()
  
  use grid_coms, only: ngrids,nzg
  use ed_work_vars, only: work_e
  use soil_coms, only: isoilflg
  use ed_node_coms, only: mxp, myp,mynum,iwest,jsouth
  use ed_state_vars, only: edgrid_g,gdpy
  
   implicit none 
   integer                :: ifm,ipy,x,y,npoly
   
   gloop: do ifm=1,ngrids

      npoly=gdpy(mynum,ifm)

      do ipy=1,npoly
         
         edgrid_g(ifm)%lon(ipy) = work_e(ifm)%vec_glon(ipy)
         edgrid_g(ifm)%lat(ipy) = work_e(ifm)%vec_glat(ipy)
         edgrid_g(ifm)%ntext_soil(1:nzg,ipy) = work_e(ifm)%vec_ntext(ipy)
         
         print*,mynum,edgrid_g(ifm)%lat(ipy),edgrid_g(ifm)%ntext_soil(nzg,ipy)

         ! This should be meaningless right now
         !         edgrid_g(ifm)%xatm(ipy) = work_e(ifm)%xatm(x,y) + iwest(ifm)  - 2
         !         edgrid_g(ifm)%yatm(ipy) = work_e(ifm)%yatm(x,y) + jsouth(ifm) - 2
         
         
      end do
   end do gloop

   
   return
 end subroutine set_polygon_coordinates_ar
 
 !==========================================================================================!
 
 subroutine soil_depth_fill_ar(cgrid,igr)
   !---------------------------------------------------------------------------------------!
   !    This subroutine fills the lsl variables based on the soil_depth file. In case         !
   ! isoildepthflg was zero, then the layer_index matrix was filled with zeroes, so we do not !
   ! need to worry about this here.                                                           !
   !------------------------------------------------------------------------------------------!
   
   use soil_coms, only: soildepth_db, slz, isoildepthflg,layer_index
   use grid_coms,only:ngrids, nzg
   use ed_state_vars,only: edtype
  
  implicit none

  type(edtype),target  :: cgrid
  integer, intent(in) :: igr
  integer :: ilat_bin
  integer :: ilon_bin
  integer :: ipy
    
  do ipy = 1,cgrid%npolygons
     ilat_bin = min(180,int(90.0 - cgrid%lat(ipy)) + 1)
     ilon_bin = int(180.0 + cgrid%lon(ipy)) + 1

     !------------------------------------------------------------------------------------------!
     !    Require at least 2 layers. This requirement was taken in consideration when           !
     ! layer_index was filled at the first initialization, so it is safe to just copy.          !
     !------------------------------------------------------------------------------------------!
     cgrid%lsl(ipy) =layer_index(ilat_bin,ilon_bin) 

  enddo
!----- layer_index is unecessary beyond this point. Deallocating it... ------=-------------!
  !if (igr == ngrids) deallocate(layer_index)

  return
end subroutine soil_depth_fill_ar

!==========================================================================================!

subroutine load_ecosystem_state

  use phenology_coms, only: iphen_scheme
  use misc_coms, only: ied_init_mode
  use phenology_startup, only: phenology_init
  use ed_node_coms, only: mynum,nmachs,nnodetot,mchnum,machs,master_num,sendnum,recvnum
  use grid_coms,only : ngrids
  use ed_state_vars,only : edgrid_g

  implicit none
  include 'mpif.h'
  
  integer,       dimension(MPI_STATUS_SIZE) :: status
  integer                                   :: ierr
  integer :: igr
  integer :: ping 

  ping = 741776
  

  ! ----------------------------------------------------
  ! Several Procedures requiring ASCII reads follow
  ! If this is a parallel run, then the nodes must queue.
  ! Since we will access sequential format files, each node needs to wait its turn
  ! to access it... MPI_File commands won't work with ascii files, so that's a 
  ! bottleneck here. If the run is serial mynum=nnodetot, so I don't need to 
  ! wait.
  ! ----------------------------------------------------


  if (mynum .eq. 1)print'(/,a)','    Doing sequential initialization over nodes'
  !----------------------------------------------------
  ! STEP 0: Find lowest soil layer for each site 
  !         (derived from soil depth)!
  !----------------------------------------------------
  
  do igr=1,ngrids
     call ed_newgrid(igr)
     call soil_depth_fill_ar(edgrid_g(igr),igr)
  end do


  ! ----------------------------------------------------
  ! STEP 1: Read in Site files and initialize 
  !         hydrologic adjacencies
  ! ----------------------------------------------------

  if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,734,MPI_COMM_WORLD,status,ierr)
  
  do igr = 1,ngrids
     call read_site_file_array(edgrid_g(igr))
  enddo
  
  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,734,MPI_COMM_WORLD,ierr)
!  if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! ----------------------------------------------------
  ! STEP 2: Do ascii type restart initialization of site
  !         patch and cohort biophysical states
  ! ----------------------------------------------------
  
  if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,799,MPI_COMM_WORLD,status,ierr)
  
  
  select case (ied_init_mode)
  case(0)
     
     ! Initialize everything with near-bare ground
     if (mynum /= 1) print'(/,a)','    Doing bare ground initialization'
     do igr=1,ngrids
          call bare_ground_init(edgrid_g(igr))
     end do
     
  case(-1,1,2,3)
     
     ! Initialize with ED1-type restart info
     
     write(*,'(a,i3.3)')'    Initializing from ED restart file. Node: ',mynum
     
     call read_ed1_history_file_array
     
  end select

  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,799,MPI_COMM_WORLD,ierr)
!  if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  

  ! ----------------------------------------------------
  ! STEP 3: Initialize phenology parameters and thermal sums
  ! ----------------------------------------------------
  
  if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,735,MPI_COMM_WORLD,status,ierr)
  
  write(*,'(a,i3.3)')'    Initializing phenology. Node: ',mynum
  call phenology_init
  
  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,735,MPI_COMM_WORLD,ierr)
!  if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  
  ! ----------------------------------------------------
  ! STEP 4: Initialize anthropogenic disturbance
  ! ----------------------------------------------------
  
  if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,736,MPI_COMM_WORLD,status,ierr)
  

  write(*,'(a,i3.3)')'    Initializing anthropogenic disturbance forcing. Node: ',mynum

  call landuse_init_array

  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,736,MPI_COMM_WORLD,ierr)
!  if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  

  ! Initialize the lapse rates that transfer meteorologic forcing variables
  ! from the atmospheric reference height to a site specific quanity.



  return
end subroutine load_ecosystem_state

!=============================================================================!

subroutine sfcdata_ed()
  use grid_coms, only: nzg
  use soil_coms, only: ed_nstyp, slz, dslz, dslzo2, dslzi, &
                       dslzidt, slzt, dslzt, dslzti, dslztidt,  &
                       fhydraul, slcons1, slden,  emisg, &
 ! Using the table defined in soil_coms instead of redefining locally.
                       soil
                       
  use misc_coms, only: dtlsm
  implicit none

  integer :: k
  integer :: nnn

  real :: refdepth


  ! Soil vertical grid spacing arrays (some with timestep info)

  slz(nzg+1) = 0.

  do k = 1,nzg
     dslz   (k) = slz(k+1) - slz(k)
     dslzo2 (k) = .5 * dslz(k)
     dslzi  (k) = 1. / dslz(k)
     dslzidt(k) = dslzi(k) * dtlsm
     slzt   (k) = .5 * (slz(k) + slz(k+1))
  enddo

  do k = 2,nzg
     dslzt   (k) = slzt(k) - slzt(k-1)
     dslzti  (k) = 1. / dslzt(k)
     dslztidt(k) = dslzti(k) * dtlsm
  enddo

  ! Soil constants

  refdepth = -2.0

  do nnn = 1,ed_nstyp
     fhydraul(nnn) = log (soil(nnn)%slcons / soil(nnn)%slcons0) / refdepth

     do k = 1,nzg
        slcons1(k,nnn) = soil(nnn)%slcons     ! ORIGINAL form - const with depth
  !     slcons1(k,nnn) = soilparms(5,nnn)  &  ! TOPMODEL form - large at surface
  !        * exp(slz(k) * fhydraul(nnn))      !    and exp decrease with depth
     enddo

     slden    (nnn) =  soil(nnn)%slden    
     emisg (nnn) = .98
  enddo

  return
end subroutine sfcdata_ed

