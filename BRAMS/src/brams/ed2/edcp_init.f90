!==========================================================================================!
!==========================================================================================!
!    This subroutine controls the initialisation at the head node.                         !
!------------------------------------------------------------------------------------------!
subroutine master_ed_init(iparallel)

   use rpara         , only : mainnum            & ! intent(in)
                            , nmachs             ! ! intent(out)
   use node_mod      , only : mynum              & ! intent(out)
                            , machs              & ! intent(out)
                            , mchnum             ! ! intent(out)
   use soil_coms     , only : layer_index        & ! intent(out)
                            , nlon_lyr           & ! intent(out)
                            , nlat_lyr           ! ! intent(out)
   use mem_leaf      , only : isfcl              ! ! intent(in)
   use mem_grid      , only : ngrids             ! ! intent(in)
   use ed_state_vars , only : gdpy               & ! intent(inout)
                            , py_off             & ! intent(inout)
                            , allocate_edglobals & ! subroutine
                            , allocate_edtype    & ! subroutine
                            , edgrid_g           ! ! intent(inout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: iparallel
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ifm
   integer             :: ierr
   !---------------------------------------------------------------------------------------!
   include 'mpif.h'

   if(isfcl /= 5) return

   if (iparallel == 1) then
      !----- Initialize the work arrays. --------------------------------------------------!
      call init_master_work(iparallel)

      !----- Read the soil depth database. ------------------------------------------------!
      call read_soil_depth()

      !----- Send soil depths to the nodes. -----------------------------------------------!
      call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Just to wait until the matrix is allocated
      call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   else
      !----- Setting up a serial run. -----------------------------------------------------!
      mynum  = 1
      nmachs = 1
      call copy_in_bramsmpi(mainnum,mchnum,mynum,nmachs,machs,iparallel)

      !----- Initialize the work arrays. --------------------------------------------------!
      call init_master_work(iparallel)

      !----- Read the soil depth database. ------------------------------------------------!
      call read_soil_depth()

      !----- Allocate the polygons on edgrid. ---------------------------------------------!
      write (unit=*,fmt='(a,i5,a)') ' + Polygon array allocation, node ',mynum,';'
      call allocate_edglobals(ngrids)
      do ifm=1,ngrids
         call ed_newgrid(ifm)
         call allocate_edtype(edgrid_g(ifm),gdpy(mynum,ifm))
      end do

      write (unit=*,fmt='(a,i5,a)') ' + Memory successfully allocated on none ',mynum,';'
      call onenode()
      call ed_coup_driver()

   end if
   return
end subroutine master_ed_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine node_ed_init
  
  
  ! BRAMS MODULES
  use node_mod,only: master_num,mchnum,mynum,nmachs,machs
  use mem_grid,only: ngrids
  use rpara,only: mainnum
  use mem_leaf,only:isfcl
  ! ED MODULES
  use soil_coms,only: layer_index,nlon_lyr,nlat_lyr
  use ed_state_vars,only:  &
       gdpy,   &
       py_off, &
       allocate_edglobals, &
       allocate_edtype,    &
       edgrid_g

  implicit none
  include 'mpif.h'
  integer :: ifm,ierr
  
  ! First we must transfer over the parallel information
  ! from BRAMS to ED2
  
  if(isfcl /= 5) return

  call copy_in_bramsmpi(master_num,mchnum,mynum,nmachs,machs,1)
  

  ! Calculate the polygon list on the current node
  call init_node_work

  ! Recieve the lowest soil layer index array
  if (allocated(layer_index)) deallocate(layer_index)
  allocate(layer_index(nlat_lyr,nlon_lyr))
  call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Safe to receive the data.
  call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
  ! Allocate the polygons on edgrid
  write (unit=*,fmt='(a,i5,a)') ' + Polygon array allocation, node ',mynum,';'
  call allocate_edglobals(ngrids)
  do ifm=1,ngrids
     call ed_newgrid(ifm)
     call allocate_edtype(edgrid_g(ifm),gdpy(mynum,ifm))
  end do
  
  write (unit=*,fmt='(a,i5,a)') ' + Memory successfully allocated on none ',mynum,';'

  call ed_coup_driver()


  return
end subroutine node_ed_init

! ===========================================================================================

subroutine copy_in_bramsmpi(master_num_b,mchnum_b,mynum_b,nmachs_b,machs_b,ipara)

  use ed_node_coms,only: master_num,mchnum,mynum,nmachs,machs,nnodetot,recvnum,sendnum
  use ed_para_coms,only: iparallel,mainnum
  use ed_max_dims,only:maxgrds
  use grid_coms, only : ngrids,nnxp,nnyp

  implicit none
  include 'mpif.h'
  integer :: ierr
  integer :: master_num_b
  integer :: mchnum_b
  integer :: mynum_b
  integer :: nmachs_b
  integer :: igr,ngrids_b
  integer,dimension(maxgrds) :: nnxp_b,nnyp_b
  integer,dimension(nmachs_b) :: machs_b
  integer, intent(in) :: ipara

  ! Saving parallel status (now it may not be parallel)
  iparallel = ipara
  
  master_num = master_num_b
  mchnum     = mchnum_b
  mynum      = mynum_b
  nmachs     = nmachs_b
  machs(1:nmachs)      = machs_b(1:nmachs)
  nnodetot   = nmachs_b

  recvnum = mynum-1
  sendnum = mynum+1

  return
end subroutine copy_in_bramsmpi

!==========================================================================================!

subroutine init_node_work
  
  
  use ed_work_vars       , only : work_e                 & ! intent(out)
                                , work_v                 & ! intent(out)
                                , ed_alloc_work          & ! subroutine
                                , ed_nullify_work        ! ! subroutine
  
  use grid_coms, only : ngrids

  use ed_node_coms, only: mynum,nmachs,nnodetot,mchnum,machs,master_num,sendnum,recvnum,mmxp,mmyp, &
       iwest,jsouth

  use rpara,only : mainnum
  
  use ed_state_vars,only:  gdpy,   &
       py_off, &
       allocate_edglobals, &
       allocate_edtype,    &
       edgrid_g

  implicit none
  include 'mpif.h'
  integer :: ifm,nm,nm2,ierr,i,j
  integer :: ipy
  integer,       dimension(MPI_STATUS_SIZE) :: status

  allocate(work_e(ngrids),work_v(ngrids))
  do ifm = 1,ngrids
          
     do nm=1,nmachs 

        call MPI_Bcast(nm2,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
        
        if(nm2/=nm) then
           print*,"THE BROADCAST NODE DOES NOT MATCH THE LOOP NODE"
           stop
        endif

        call MPI_Bcast(gdpy(nm,ifm),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(py_off(nm,ifm),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)


        if(nm==mynum) then
           
           call MPI_Recv(mmxp(ifm),1,MPI_INTEGER,mainnum,61,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(mmyp(ifm),1,MPI_INTEGER,mainnum,62,MPI_COMM_WORLD,status,ierr)

           call MPI_Recv(iwest(ifm),1,MPI_INTEGER,mainnum,63,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(jsouth(ifm),1,MPI_INTEGER,mainnum,64,MPI_COMM_WORLD,status,ierr)

           
           call ed_nullify_work(work_e(ifm))
           call ed_alloc_work(work_e(ifm),mmxp(ifm),mmyp(ifm))

           call MPI_Recv(work_e(ifm)%glat,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,65,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%glon,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,66,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%work,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,67,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%land,mmxp(ifm)*mmyp(ifm),MPI_LOGICAL,mainnum,68,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%landfrac,mmxp(ifm)*mmyp(ifm),MPI_REAL,mainnum,69,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%ntext,mmxp(ifm)*mmyp(ifm),MPI_INTEGER,mainnum,70,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%xatm,mmxp(ifm)*mmyp(ifm),MPI_INTEGER,mainnum,71,MPI_COMM_WORLD,status,ierr)
           call MPI_Recv(work_e(ifm)%yatm,mmxp(ifm)*mmyp(ifm),MPI_INTEGER,mainnum,72,MPI_COMM_WORLD,status,ierr)

        endif

     enddo


     ! Fill the work vectors - these will be used in the ed2 initialization procedures
     ! to populate the first polygons
     ! -------------------------------------------------------------------------------
     call edcp_parvec_work(ifm,mmxp(ifm),mmyp(ifm))


  enddo
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
        
end subroutine init_node_work


!==========================================================================================!

subroutine init_master_work(ipara)
  
  use mem_grid,only: grid_g,ngrids,nnxp,nnyp,jdim
  use rpara,only: ixb,ixe,iyb,iye,nmachs,mainnum,machnum
  use ed_work_vars,only : &
       work_e,                 & ! intent(out)
       work_v,                 & ! intent(out)
       ed_alloc_work,          & ! subroutine
       ed_nullify_work,        &         ! ! subroutine
       ed_dealloc_work

  use ed_state_vars,only:  &
       gdpy,   &
       py_off, &
       allocate_edglobals, &
       allocate_edtype,    &
       edgrid_g
  use ed_node_coms, only: mmxp,mmyp

  implicit none
  include 'mpif.h'
  integer, intent(in) :: ipara
  integer :: ierr
  integer :: nm,ifm
  integer :: offset,npolys
  integer :: xmax,ymax
  integer :: i,j,il,jl
  integer :: iwest,ieast,jsouth,jnorth
  integer :: ipy

  allocate(work_e(ngrids))
  if (ipara /= 1) allocate(work_v(ngrids))
  
  gridloop: do ifm = 1,ngrids
     call newgrid(ifm)
     npolys=0
     offset=0
     machloop: do nm=1,nmachs 
        
        if (ipara == 1) then
           iwest  = ixb(nm,ifm)
           ieast  = ixe(nm,ifm)
           jsouth = iyb(nm,ifm)
           jnorth = iye(nm,ifm)
        else
           iwest  = 2
           ieast  = nnxp(ifm)-1
           jsouth = 1+jdim
           jnorth = nnyp(ifm)-jdim
        end if

        xmax = ieast  - iwest  + 1
        ymax = jnorth - jsouth + 1

        call ed_nullify_work(work_e(ifm))
        call ed_alloc_work(work_e(ifm),xmax,ymax)

        jl=0
        jloop: do j=jsouth,jnorth
           jl=jl+1
           il=0
           iloop: do i=iwest,ieast
              il=il+1
              work_e(ifm)%glon(il,jl) = grid_g(ifm)%glon(i,j)
              work_e(ifm)%glat(il,jl) = grid_g(ifm)%glat(i,j)
              
              ! The following will give you the tile (local) coordinates
              work_e(ifm)%xatm(il,jl)  = i - iwest  + 2 ! Remember that all tiles have a 
              work_e(ifm)%yatm(il,jl)  = j - jsouth + 2 ! buffer cell so add one extra 

           end do iloop
        end do jloop
        
        call edcp_get_work(ifm,xmax,ymax,iwest,ieast,jsouth,jnorth)

        offset=offset+npolys
        npolys=count(work_e(ifm)%land)
           
        gdpy(nm,ifm)   = npolys
        py_off(nm,ifm) = offset

        if (ipara == 1) then
           call MPI_Bcast(nm,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
           call MPI_Bcast(npolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
           call MPI_Bcast(offset,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

           ! Send the work grids to the nodes

           call MPI_Send(xmax,1,MPI_INTEGER,machnum(nm),61,MPI_COMM_WORLD,ierr)
           call MPI_Send(ymax,1,MPI_INTEGER,machnum(nm),62,MPI_COMM_WORLD,ierr)
           
           call MPI_Send(iwest,1,MPI_INTEGER,machnum(nm),63,MPI_COMM_WORLD,ierr)
           call MPI_Send(jsouth,1,MPI_INTEGER,machnum(nm),64,MPI_COMM_WORLD,ierr)

           call MPI_Send(work_e(ifm)%glat,xmax*ymax,MPI_REAL,machnum(nm),65,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%glon,xmax*ymax,MPI_REAL,machnum(nm),66,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%work,xmax*ymax,MPI_REAL,machnum(nm),67,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%land,xmax*ymax,MPI_LOGICAL,machnum(nm),68,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%landfrac,xmax*ymax,MPI_REAL,machnum(nm),69,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%ntext,xmax*ymax,MPI_INTEGER,machnum(nm),70,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%xatm,xmax*ymax,MPI_INTEGER,machnum(nm),71,MPI_COMM_WORLD,ierr)
           call MPI_Send(work_e(ifm)%yatm,xmax*ymax,MPI_INTEGER,machnum(nm),72,MPI_COMM_WORLD,ierr)
    
           call ed_dealloc_work(work_e(ifm))
        else
           mmxp(ifm) = xmax
           mmyp(ifm) = ymax

           ! Fill the work vectors - these will be used in the ed2 initialization procedures
           ! to populate the first polygons
           ! -------------------------------------------------------------------------------
           call edcp_parvec_work(ifm,mmxp(ifm),mmyp(ifm))
     
        end if
     
     end do machloop
      
  end do gridloop
  if (ipara == 1) then
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     deallocate(work_e)
  end if

  return
end subroutine init_master_work
!==========================================================================================!
!==========================================================================================!
