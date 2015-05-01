!==========================================================================================!
!==========================================================================================!
!    This subroutine controls the initialisation at the head node.                         !
!------------------------------------------------------------------------------------------!
subroutine master_ed_init(iparallel)
   !----- BRAMS modules. ------------------------------------------------------------------!
   use rpara         , only : mainnum            & ! intent(in)
                            , nmachs             ! ! intent(out)
   use node_mod      , only : mynum              & ! intent(out)
                            , machs              & ! intent(out)
                            , mchnum             ! ! intent(out)
   use mem_leaf      , only : isfcl              ! ! intent(in)
   use mem_grid      , only : ngrids             ! ! intent(in)
   !----- ED modules. ---------------------------------------------------------------------!
   use soil_coms     , only : layer_index        & ! intent(out)
                            , nlon_lyr           & ! intent(out)
                            , nlat_lyr           ! ! intent(out)
   use ed_state_vars , only : gdpy               & ! intent(inout)
                            , py_off             & ! intent(inout)
                            , allocate_edglobals & ! subroutine
                            , allocate_edtype    & ! subroutine
                            , edgrid_g           ! ! intent(inout)
   implicit none
   !------Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: iparallel
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ifm
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !----- No need to do anything here if this is not a coupled ED2-BRAMS run. -------------!
   if (isfcl /= 5) return


   if (iparallel == 1) then
      !----- Initialize the work arrays. --------------------------------------------------!
      call init_master_work(iparallel)

      !----- Read the soil depth database. ------------------------------------------------!
      call read_soil_depth()

      !----- Send soil depths to the nodes. -----------------------------------------------!
#if defined(RAMS_MPI)
      call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Just to wait until the matrix is allocated
      call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
#endif

   else
      !----- Set up a serial run. ---------------------------------------------------------!
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
!    This subroutine controls the initialisation at the slave node.                        !
!------------------------------------------------------------------------------------------!
subroutine node_ed_init
   !----- BRAMS modules. ------------------------------------------------------------------!
   use node_mod     , only : master_num         & ! intent(in)
                           , mchnum             & ! intent(in)
                           , mynum              & ! intent(in)
                           , nmachs             & ! intent(in)
                           , machs              ! ! intent(in)
   use mem_grid     , only : ngrids             ! ! intent(in)
   use rpara        , only : mainnum            ! ! intent(in)
   use mem_leaf     , only : isfcl              ! ! intent(in)
   !----- ED modules. ---------------------------------------------------------------------!
   use soil_coms    , only : layer_index        & ! intent(inout)
                           , nlon_lyr           & ! intent(in)
                           , nlat_lyr           ! ! intent(in)
   use ed_state_vars, only : gdpy               & ! intent(in)
                           , py_off             & ! intent(in)
                           , allocate_edglobals & ! sub-routine
                           , allocate_edtype    & ! sub-routine
                           , edgrid_g           ! ! structure
   implicit none
   !------Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ifm
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !----- Nothing to be done here if this is not an ED2-BRAMS run. ------------------------!
   if(isfcl /= 5) return

   !----- First we must transfer over the parallel information from BRAMS to ED2. ---------!
   call copy_in_bramsmpi(master_num,mchnum,mynum,nmachs,machs,1)
   

   !----- Calculate the polygon list on the current node. ---------------------------------!
   call init_node_work()

   !----- Receive the lowest soil layer index array. --------------------------------------!
   if (allocated(layer_index)) deallocate(layer_index)
   allocate(layer_index(nlat_lyr,nlon_lyr))
#if defined(RAMS_MPI)
   call MPI_Barrier(MPI_COMM_WORLD,ierr) ! Safe to receive the data.
   call MPI_Bcast(layer_index,nlat_lyr*nlon_lyr,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
#endif
   
   !----- Allocate the polygons on edgrid. ------------------------------------------------!
   write (unit=*,fmt='(a,i5,a)') ' + Polygon array allocation, node ',mynum,';'
   call allocate_edglobals(ngrids)
   do ifm=1,ngrids
      call ed_newgrid(ifm)
      call allocate_edtype(edgrid_g(ifm),gdpy(mynum,ifm))
   end do
   
   write (unit=*,fmt='(a,i5,a)') ' + Memory successfully allocated on none ',mynum,';'

   !----- Call the main initialisation driver. --------------------------------------------!
   call ed_coup_driver()

   return
end subroutine node_ed_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will copy some of the BRAMS parallel information to the ED          !
! structures.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine copy_in_bramsmpi(master_num_b,mchnum_b,mynum_b,nmachs_b,machs_b,ipara)
   use ed_node_coms, only : master_num & ! intent(out)
                          , mchnum     & ! intent(out)
                          , mynum      & ! intent(out)
                          , nmachs     & ! intent(out)
                          , machs      & ! intent(out)
                          , nnodetot   & ! intent(out)
                          , recvnum    & ! intent(out)
                          , sendnum    ! ! intent(out)
   use ed_para_coms, only : iparallel  & ! intent(out)
                          , mainnum    ! ! intent(out)
   use ed_max_dims , only : maxgrds    ! ! intent(out)
   use grid_coms   , only : ngrids     & ! intent(out)
                          , nnxp       & ! intent(out)
                          , nnyp       ! ! intent(out)

   implicit none
   !----- Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in) :: master_num_b
   integer                     , intent(in) :: mchnum_b
   integer                     , intent(in) :: mynum_b
   integer                     , intent(in) :: nmachs_b
   integer, dimension(nmachs_b), intent(in) :: machs_b
   integer                     , intent(in) :: ipara
   !---------------------------------------------------------------------------------------!

   !----- Save the parallel status (i.e. whether this is a serial or a parallel run). -----!
   iparallel = ipara

   !----- Copy the BRAMS information to the ED structures. --------------------------------!
   master_num      = master_num_b
   mchnum          = mchnum_b
   mynum           = mynum_b
   nmachs          = nmachs_b
   machs(1:nmachs) = machs_b(1:nmachs)
   nnodetot        = nmachs_b

   !---------------------------------------------------------------------------------------!
   !     Define the number of the nodes that will send information about when it is safe   !
   ! for this node to start the initialisation (recvnum), and which node the current node  !
   ! should inform that the initialisation is done locally (sendnum).                      !
   !---------------------------------------------------------------------------------------!
   recvnum = mynum-1
   sendnum = mynum+1

   return
end subroutine copy_in_bramsmpi
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine controls the initialisation of the "work" structures (which will    !
! tell which grid points should call ED), and send this information to the nodes.          !
!------------------------------------------------------------------------------------------!
subroutine init_master_work(ipara)
   !----- BRAMS modules. ------------------------------------------------------------------!
   use mem_grid     , only : grid_g              & ! intent(in)
                           , ngrids              & ! intent(in)
                           , nnxp                & ! intent(in)
                           , nnyp                & ! intent(in)
                           , jdim                ! ! intent(in)
   use rpara        , only : nxbeg               & ! intent(in)
                           , nxend               & ! intent(in)
                           , nybeg               & ! intent(in)
                           , nyend               & ! intent(in)
                           , nxbegc              & ! intent(in)
                           , nxendc              & ! intent(in)
                           , nybegc              & ! intent(in)
                           , nyendc              & ! intent(in)
                           , ixoff               & ! intent(in)
                           , iyb                 & ! intent(in)
                           , iye                 & ! intent(in)
                           , iyoff               & ! intent(in)
                           , nmachs              & ! intent(in)
                           , mainnum             & ! intent(in)
                           , machnum             ! ! intent(in)
   !----- ED modules. ---------------------------------------------------------------------!
   use ed_work_vars , only : work_e              & ! intent(out)
                           , work_v              & ! intent(out)
                           , ed_alloc_work       & ! sub-routine
                           , ed_nullify_work     & ! sub-routine
                           , ed_dealloc_work     ! ! sub-routine
   use ed_state_vars, only : gdpy                & ! intent(out)
                           , py_off              & ! intent(out)
                           , allocate_edglobals  & ! sub-routine
                           , allocate_edtype     ! ! sub-routine
   use mem_polygons , only : maxsite             ! ! intent(in)
   implicit none
   !----- Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ipara
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ierr
   integer             :: nm
   integer             :: ifm
   integer             :: offset
   integer             :: npolys
   integer             :: xmax
   integer             :: ymax
   integer             :: i
   integer             :: j
   integer             :: il
   integer             :: jl
   integer             :: iwest
   integer             :: ieast
   integer             :: iskip
   integer             :: jsouth
   integer             :: jnorth
   integer             :: jskip
   integer             :: ipy
   !---------------------------------------------------------------------------------------!


   !----- Allocate the work structures (the vector one only if it is a serial run). -------!
   allocate(work_e(ngrids))
   if (ipara /= 1) allocate(work_v(ngrids))

   gridloop: do ifm = 1,ngrids
      call newgrid(ifm)
      npolys=0
      offset=0
      machloop: do nm=1,nmachs 
         !---------------------------------------------------------------------------------!
         !     Determine the sub-domain edges for each machine, and send this information  !
         ! to them.  In case this is a serial run, the sub-domain is the full domain.      !
         ! Now the work sctructures will be allocated exactly like the other structures,   !
         ! to make sure that the grid cell actually matches the polygon.                   !
         !---------------------------------------------------------------------------------!
         if (ipara == 1) then
            iwest  = nxbegc  (nm,ifm)
            ieast  = nxendc  (nm,ifm)
            iskip  = ixoff   (nm,ifm)
            jsouth = nybegc  (nm,ifm)
            jnorth = nyendc  (nm,ifm)
            jskip  = iyoff   (nm,ifm)
            xmax   = nxend   (nm,ifm) - nxbeg(nm,ifm) + 1
            ymax   = nyend   (nm,ifm) - nybeg(nm,ifm) + 1
         else
            iwest  = 2
            ieast  = nnxp(ifm)-1
            iskip  = 0
            jsouth = 1+jdim
            jnorth = nnyp(ifm)-jdim
            jskip  = 0
            xmax   = nnxp(ifm)
            ymax   = nnyp(ifm)
         end if
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a)')       '    ==== ED node splitting to sub-domains ====    '
         write(unit=*,fmt='(a,1x,i6)') ' - GRID =',ifm
         write(unit=*,fmt='(a,1x,i6)') ' - NXP  =',nnxp(ifm)
         write(unit=*,fmt='(a,1x,i6)') ' - NYP  =',nnyp(ifm)
         write(unit=*,fmt='(a,1x,i6)') ' - NODE =',nm
         write(unit=*,fmt='(a,1x,i6)') ' - MXP  =',xmax
         write(unit=*,fmt='(a,1x,i6)') ' - MYP  =',ymax
         write(unit=*,fmt='(a,1x,i6)') ' - IA   =',iwest
         write(unit=*,fmt='(a,1x,i6)') ' - IZ   =',ieast
         write(unit=*,fmt='(a,1x,i6)') ' - I0   =',iskip
         write(unit=*,fmt='(a,1x,i6)') ' - JA   =',jsouth
         write(unit=*,fmt='(a,1x,i6)') ' - JZ   =',jnorth
         write(unit=*,fmt='(a,1x,i6)') ' - J0   =',jskip

         !----- Allocate the arrays of the work structure. --------------------------------!
         write(unit=*,fmt='(a)')       ' - Nullify and allocate work_e.'
         call ed_nullify_work(work_e(ifm))
         call ed_alloc_work(work_e(ifm),xmax,ymax,maxsite)

         write(unit=*,fmt='(a)')       ' - Assign the grid co-ordinates.'
         jloop: do jl=1,ymax
            j = jl + jskip
            iloop: do il=1,xmax
               i = il + iskip
               work_e(ifm)%glon(il,jl) = grid_g(ifm)%glon(i,j)
               work_e(ifm)%glat(il,jl) = grid_g(ifm)%glat(i,j)
               
               !---------------------------------------------------------------------------!
               !     The following variables will contain the node's local coordinates     !
               ! (the co-ordinates of the node sub-domain.                                 !
               !---------------------------------------------------------------------------!
               work_e(ifm)%xatm(il,jl)  = il
               work_e(ifm)%yatm(il,jl)  = jl
               
               !---------------------------------------------------------------------------!
               !     Initialise the remaining variables.  Default is to not have a polygon !
               ! because the edges should never be solved by ED-2, instead the fluxes and  !
               ! state should be grabbed from the neighbour node (or copied from the       !
               ! next-to-the-edge polygon in case of the true edge).                       !
               !---------------------------------------------------------------------------!
               work_e(ifm)%work    (il,jl)   = 0.
               work_e(ifm)%land    (il,jl)   = .false.
               work_e(ifm)%landfrac(il,jl)   = 0.
               work_e(ifm)%nscol   (il,jl)   = 0
               work_e(ifm)%ntext   (:,il,jl) = 0
               work_e(ifm)%soilfrac(:,il,jl) = 0.
               !---------------------------------------------------------------------------!
            end do iloop
         end do jloop

         write(unit=*,fmt='(a)')       ' - Get the workload.'
         call edcp_get_work(ifm,nnxp(ifm),nnyp(ifm),nm                                     &
                           ,xmax,ymax,iwest,ieast,iskip,jsouth,jnorth,jskip)

         write(unit=*,fmt='(a)')       ' - Define actual number of polygons and offset.'
         offset = offset + npolys
         npolys = count(work_e(ifm)%land)
            
         gdpy(nm,ifm)   = npolys
         py_off(nm,ifm) = offset
         write(unit=*,fmt='(a,1x,i6)') '   + NPOLYS   =',npolys
         write(unit=*,fmt='(a,1x,i6)') '   + OFF-SET  =',offset

         !---------------------------------------------------------------------------------!
         !     If this is a parallel run, send the work grids and information to the       !
         ! nodes.                                                                          !
         !---------------------------------------------------------------------------------!
         if (ipara == 1) then
#if defined(RAMS_MPI)
            write(unit=*,fmt='(a)')       ' - Send the ED globals to all nodes.'
            call MPI_Bcast(nm,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(npolys,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(offset,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

            !----- Send the sub-domain dimension. -----------------------------------------!
            write(unit=*,fmt='(a)')       ' - Send the lon/lat information to the node.'
            call MPI_Send(xmax   ,1,MPI_INTEGER,machnum(nm),57,MPI_COMM_WORLD,ierr)
            call MPI_Send(ymax   ,1,MPI_INTEGER,machnum(nm),58,MPI_COMM_WORLD,ierr)
            call MPI_Send(iwest  ,1,MPI_INTEGER,machnum(nm),59,MPI_COMM_WORLD,ierr)
            call MPI_Send(ieast  ,1,MPI_INTEGER,machnum(nm),60,MPI_COMM_WORLD,ierr)
            call MPI_Send(iskip  ,1,MPI_INTEGER,machnum(nm),61,MPI_COMM_WORLD,ierr)
            call MPI_Send(jsouth ,1,MPI_INTEGER,machnum(nm),62,MPI_COMM_WORLD,ierr)
            call MPI_Send(jnorth ,1,MPI_INTEGER,machnum(nm),63,MPI_COMM_WORLD,ierr)
            call MPI_Send(jskip  ,1,MPI_INTEGER,machnum(nm),64,MPI_COMM_WORLD,ierr)

            !----- Send the matrices. -----------------------------------------------------!
            write(unit=*,fmt='(a)')       ' - Send the matrices to the node.'
            call MPI_Send(work_e(ifm)%glat    ,        xmax*ymax,MPI_REAL   ,machnum(nm)   &
                         ,65,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%glon    ,        xmax*ymax,MPI_REAL   ,machnum(nm)   &
                         ,66,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%work    ,        xmax*ymax,MPI_REAL   ,machnum(nm)   &
                         ,67,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%land    ,        xmax*ymax,MPI_LOGICAL,machnum(nm)   &
                         ,68,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%landfrac,        xmax*ymax,MPI_REAL   ,machnum(nm)   &
                         ,69,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%xatm    ,        xmax*ymax,MPI_INTEGER,machnum(nm)   &
                         ,70,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%yatm    ,        xmax*ymax,MPI_INTEGER,machnum(nm)   &
                         ,71,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%nscol   ,        xmax*ymax,MPI_INTEGER,machnum(nm)   &
                         ,72,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%ntext   ,maxsite*xmax*ymax,MPI_INTEGER,machnum(nm)   &
                         ,73,MPI_COMM_WORLD,ierr)
            call MPI_Send(work_e(ifm)%soilfrac,maxsite*xmax*ymax,MPI_REAL   ,machnum(nm)   &
                         ,74,MPI_COMM_WORLD,ierr)

            !------------------------------------------------------------------------------!
            !     De-allocate the internal matrices, so we can re-allocate them for the    !
            ! sub-domain that goes to the next machine (or free memory once we are done).  !
            !------------------------------------------------------------------------------!
            write(unit=*,fmt='(a)')       ' - De-allocate work_e on head node.'
            call ed_dealloc_work(work_e(ifm))
#endif
         else

            !------------------------------------------------------------------------------!
            !      Fill the work vectors - these will be used in the ed2 initialisation    !
            ! procedures to populate the first polygons.                                   !
            !------------------------------------------------------------------------------!
            write(unit=*,fmt='(a)')       ' - Send the grid information to myself.'
            call edcp_nodebounds_self(ifm,xmax,ymax,iwest,ieast,iskip,jsouth,jnorth,jskip)
            write(unit=*,fmt='(a)')       ' - Copy workload variables to the vector work.'
            call edcp_parvec_work(ifm,xmax,ymax,iwest,ieast,iskip,jsouth,jnorth,jskip)
         end if

         write(unit=*,fmt='(a)')       ' - Successfully sent information to the node.'
         write(unit=*,fmt='(a)')       '--------------------------------------------------'

      end do machloop
   end do gridloop
   if (ipara == 1) then
      !----- Once all slaves have received the information, it is safe to exit. -----------!
#if defined(RAMS_MPI)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      deallocate(work_e)
   end if

   return
end subroutine init_master_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine controls the initialisation of the "work" structures at the slave   !
! nodes.  The work structures tell which grid points should call ED.  The information will !
! be collected from the message sent by the master node.                                   !
!------------------------------------------------------------------------------------------!
subroutine init_node_work()
   use ed_work_vars , only : work_e                 & ! intent(out)
                           , work_v                 & ! intent(out)
                           , ed_alloc_work          & ! subroutine
                           , ed_nullify_work        ! ! subroutine
   use grid_coms    , only : ngrids                 ! ! intent(in)
   use ed_node_coms , only : mynum                  & ! intent(in)
                           , nmachs                 & ! intent(in)
                           , machs                  & ! intent(in)
                           , mmxp                   & ! intent(out)
                           , mmyp                   & ! intent(out)
                           , iwest                  & ! intent(out)
                           , ieast                  & ! intent(out)
                           , iskip                  & ! intent(out)
                           , jsouth                 & ! intent(out)
                           , jnorth                 & ! intent(out)
                           , jskip                  ! ! intent(out)
   use rpara        , only : mainnum                ! ! intent(in)
   use ed_state_vars, only : gdpy                   & ! intent(inout)
                           , py_off                 & ! intent(inout)
                           , allocate_edglobals     & ! sub-routine
                           , allocate_edtype        ! ! sub-routine
   use mem_polygons , only : maxsite                ! ! sub-routine
   implicit none
   !----- Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: ifm
   integer                             :: nm
   integer                             :: nm2
   integer                             :: ierr
   integer                             :: i
   integer                             :: j
   integer                             :: ipy
#if defined(RAMS_MPI)
   integer, dimension(MPI_STATUS_SIZE) :: status
#endif
   !---------------------------------------------------------------------------------------!

   !----- Allocate the work structures. ---------------------------------------------------!
   allocate(work_e(ngrids),work_v(ngrids))


   !----- Loop over all grids and all machines. -------------------------------------------!
   do ifm = 1,ngrids
#if defined(RAMS_MPI)
      do nm=1,nmachs 

         !----- Sanity check, make sure that the IDs are not messed up. -------------------!
         call MPI_Bcast(nm2,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
         if (nm2 /= nm) then
            call abort_run('The broadcast node doesn''t match the loop node.'              &
                          ,'init_node_work','edcp_init.f90')
         end if

         !----- Get and send the polygon size and off-set. --------------------------------!
         call MPI_Bcast(gdpy(nm,ifm)  ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(py_off(nm,ifm),1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

         !----- If this is the information this node should get, save it. -----------------!
         if (nm == mynum) then

            !----- Get the sub-domain dimensions. -----------------------------------------!
            call MPI_Recv(mmxp(ifm)  ,1,MPI_INTEGER,mainnum,57,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(mmyp(ifm)  ,1,MPI_INTEGER,mainnum,58,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(iwest(ifm) ,1,MPI_INTEGER,mainnum,59,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(ieast(ifm) ,1,MPI_INTEGER,mainnum,60,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(iskip(ifm) ,1,MPI_INTEGER,mainnum,61,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(jsouth(ifm),1,MPI_INTEGER,mainnum,62,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(jnorth(ifm),1,MPI_INTEGER,mainnum,63,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(jskip(ifm) ,1,MPI_INTEGER,mainnum,64,MPI_COMM_WORLD,status,ierr)
            write(unit=*,fmt='(a)')       '   + Lon/lat information received, thanks!'

            !----- Allocate the work structure before receiving the sub-domain. -----------!
            call ed_nullify_work(work_e(ifm))
            call ed_alloc_work(work_e(ifm),mmxp(ifm),mmyp(ifm),maxsite)

            call MPI_Recv(work_e(ifm)%glat    ,        mmxp(ifm)*mmyp(ifm),MPI_REAL        &
                         ,mainnum,65,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%glon    ,        mmxp(ifm)*mmyp(ifm),MPI_REAL        &
                         ,mainnum,66,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%work    ,        mmxp(ifm)*mmyp(ifm),MPI_REAL        &
                         ,mainnum,67,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%land    ,        mmxp(ifm)*mmyp(ifm),MPI_LOGICAL     &
                         ,mainnum,68,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%landfrac,        mmxp(ifm)*mmyp(ifm),MPI_REAL        &
                         ,mainnum,69,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%xatm    ,        mmxp(ifm)*mmyp(ifm),MPI_INTEGER     &
                         ,mainnum,70,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%yatm    ,        mmxp(ifm)*mmyp(ifm),MPI_INTEGER     &
                         ,mainnum,71,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%nscol   ,        mmxp(ifm)*mmyp(ifm),MPI_INTEGER     &
                         ,mainnum,72,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%ntext   ,maxsite*mmxp(ifm)*mmyp(ifm),MPI_INTEGER     &
                         ,mainnum,73,MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(work_e(ifm)%soilfrac,maxsite*mmxp(ifm)*mmyp(ifm),MPI_REAL        &
                         ,mainnum,74,MPI_COMM_WORLD,status,ierr)
            write(unit=*,fmt='(a)')       ' - Matrices arrived here, thanks!'
         end if
      end do
#endif



      !------------------------------------------------------------------------------------!
      !      Fill the work vectors - these will be used in the ed2 initialisation proce-   !
      ! dures to populate the first polygons.                                              !
      !------------------------------------------------------------------------------------!
      write(unit=*,fmt='(a,1x,i6)')  ' - Copy work_e to work_v on node:',mynum
      call edcp_parvec_work(ifm,mmxp(ifm),mmyp(ifm),iwest(ifm),ieast(ifm),iskip(ifm)       &
                           ,jsouth(ifm),jnorth(ifm),jskip(ifm))
   end do

   !----- Wait until all nodes have the information to move on. ---------------------------!
#if defined(RAMS_MPI)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   return
end subroutine init_node_work
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine copies the full domain boundary information to the ed_node_coms     !
! structure in case this is a serial run.                                                  !
!------------------------------------------------------------------------------------------!
subroutine edcp_nodebounds_self(ifm,xmax,ymax,ia,iz,i0,ja,jz,j0)
   use ed_node_coms, only : mmxp   & ! intent(out)
                          , mmyp   & ! intent(out)
                          , iwest  & ! intent(out)
                          , ieast  & ! intent(out)
                          , iskip  & ! intent(out)
                          , jsouth & ! intent(out)
                          , jnorth & ! intent(out)
                          , jskip  ! ! intent(out)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: ifm
   integer, intent(in) :: xmax
   integer, intent(in) :: ymax
   integer, intent(in) :: ia
   integer, intent(in) :: iz
   integer, intent(in) :: i0
   integer, intent(in) :: ja
   integer, intent(in) :: jz
   integer, intent(in) :: j0
   !---------------------------------------------------------------------------------------!


   !----- Copy the variables... -----------------------------------------------------------!
   mmxp  (ifm) = xmax
   mmyp  (ifm) = ymax
   iwest (ifm) = ia
   ieast (ifm) = iz
   iskip (ifm) = i0
   jsouth(ifm) = ja
   jnorth(ifm) = jz
   jskip (ifm) = j0
   !---------------------------------------------------------------------------------------!

   return
end subroutine edcp_nodebounds_self
!==========================================================================================!
!==========================================================================================!
