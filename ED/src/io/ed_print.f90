!==========================================================================================!
!==========================================================================================!
!      This subroutine prints patch, cohort, polygon or site level data, upscales that     !
! data to the site level and stores the data in its spatial array coordinate.  The data is !
! then printed to the screen, based on a specified window of data.  Note, this may be      !
! printing windows on various nodes, or this may be called from a master process.  Be      !
! conscious of this; as it will dictate what part of the domain you are accessing vari-    !
! ables from, and whether or not the variable of interest is stored in memory at that      !
! time.  For instance, many variables are stored only on the slave nodes, and need not be  !
! passed back to the master.  Likewise, many slave node data will accumulate after each    !
! lsm call, until they are passed back to the master, where they are normalized.  These    !
! variables will be immediately zeroed on the slaves after being sent to the master.       !
! Don't forget to adjust the number precision on the format string at the end.             !
!------------------------------------------------------------------------------------------!
subroutine print_fields(ifm,cgrid)
 
   use ed_node_coms , only : mynum         & ! intent(in)
                           , nnodetot      & ! intent(in)
                           , sendnum       & ! intent(in)
                           , recvnum       & ! intent(in)
                           , master_num    & ! intent(in)
                           , machs         ! ! intent(in)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use ed_misc_coms , only : printvars     & ! intent(in)
                           , ipmax         & ! intent(in)
                           , ipmin         & ! intent(in)
                           , pfmtstr       & ! intent(in)
                           , iprintpolys   ! ! intent(in)
   use ed_var_tables, only : vt_info       & ! intent(in)
                           , num_var       ! ! intent(in)
   use ed_max_dims  , only : str_len_short ! ! intent(in)
   implicit none
   !----- Standard common blocks. ---------------------------------------------------------!
   include 'mpif.h'
   !----- Arguments. ----------------------------------------------------------------------!
   integer                            , intent(in) :: ifm
   type(edtype)                       , target     :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(MPI_STATUS_SIZE)             :: status
   real   , dimension(:)              , pointer    :: pvar_l
   real   , dimension(:)              , pointer    :: pvar_g
   character(len=30)                               :: fmtstr
   character(len=str_len_short)                    :: pvar_name
   integer                                         :: nv
   integer                                         :: np
   integer                                         :: i
   integer                                         :: ip
   integer                                         :: ping
   integer                                         :: npolys
   integer                                         :: g_idmin
   integer                                         :: g_idmax
   integer                                         :: l_idmin
   integer                                         :: l_idmax
   integer                                         :: ierr
   integer                                         :: node_idmin
   integer                                         :: node_idmax
   integer                                         :: mast_idmin
   integer                                         :: mast_idmax
   integer                                         :: g_id
   integer                                         :: g_ln
   integer                                         :: nm
   integer                                         :: ncols
   integer                                         :: row
   integer                                         :: maxrows
   integer                                         :: col
   logical                                         :: pvartrue
   logical                                         :: ptr_recv
   logical                                         :: ptr_send
   !----- Local constants. ----------------------------------------------------------------!
   integer                            , parameter  :: maxcols  = 10
   real                               , parameter  :: undef    = -99.9
   !---------------------------------------------------------------------------------------!

   
   
   
   !----- Ping is just any flag to be used by MPI. ----------------------------------------!
   ping   = 8675309

   !----- Find the number of polygons to be plotted. --------------------------------------!
   npolys = ipmax - ipmin + 1
   

   !----- Adjust the format string according to the chosen variables. ---------------------!

   !----- Check the window size. ----------------------------------------------------------!
   if (ipmax > cgrid%npolygons_global) then
      write(unit=*,fmt='(a)')       '====================================================='
      write(unit=*,fmt='(a,1x,i5)') ' IPMAX                     = ',ipmax
      write(unit=*,fmt='(a,1x,i5)') ' Global number of polygons = ',cgrid%npolygons_global
      write(unit=*,fmt='(a)')       ' '
      write(unit=*,fmt='(a)')       ' You have specified a print index greater than the'
      write(unit=*,fmt='(a)')       ' total number of polygons.  Please reduce IPMAX...'
      write(unit=*,fmt='(a)')       '====================================================='
      call fatal_error('Too many polygons to be written on standard output.'               &
                      ,'print_fields','ed_print.f90')
   end if


   !---------------------------------------------------------------------------------------!
   !     Allocate the print (always) and scratch vector (if this is the "master" node or a ! 
   ! serial run.                                                                           !
   !---------------------------------------------------------------------------------------!
   allocate(pvar_l(npolys))
   if (mynum == nnodetot .or. nnodetot == 1) allocate(pvar_g(npolys))
   !---------------------------------------------------------------------------------------!

   
   !----- Loop through the printvar entries from the namelist. ----------------------------!
   ip = 0
   count_pvars: do

      ip = ip+1
      pvar_name = printvars(ip)
      if (len_trim(pvar_name) == str_len_short) then
         exit count_pvars
      endif

      if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

      pvartrue = .false.
      do nv = 1,num_var(ifm)

         if (trim(vt_info(nv,ifm)%name) .eq. trim(pvar_name)) then
            pvartrue = .true.
            
            !------------------------------------------------------------------------------!
            !     If this is root, then collect the sends, keep reading to find out what   !
            ! it is receiving.                                                             !
            !------------------------------------------------------------------------------!
            if (nnodetot > 1) then
               
               if (mynum == nnodetot) then
                  
                  pvar_g = undef
                  !----- Loop through the variable table to match the print variable. -----!
                  write(unit=*,fmt='( a)') ' '
                  write(unit=*,fmt='(3a)') ' =========== ',trim(pvar_name),' =========== '
                  write(unit=*,fmt='( a)') ' '

                  do nm = 1,nnodetot-1

                     call MPI_Recv(ptr_recv,1,MPI_LOGICAL,machs(nm),120,MPI_COMM_WORLD     &
                                  ,status,ierr)

                     if (ptr_recv) then
                        call MPI_Recv(mast_idmin,1,MPI_INTEGER,machs(nm),121               &
                                     ,MPI_COMM_WORLD,status,ierr)

                        call MPI_Recv(mast_idmax,1,MPI_INTEGER,machs(nm),122               &
                                     ,MPI_COMM_WORLD,status,ierr)

                        call MPI_Recv(pvar_g(mast_idmin:mast_idmax)                        &
                                     ,mast_idmax-mast_idmin+1,MPI_REAL                     &
                                     ,machs(nm),123,MPI_COMM_WORLD,status,ierr)
                     end if

                  end do
               end if
            else
               
               !----- Loop through the variable table to match the print variable. --------!
               write(unit=*,fmt='( a)') ' '
               write(unit=*,fmt='(3a)') ' =========== ',trim(pvar_name),' =========== '
               write(unit=*,fmt='( a)') ' '
               pvar_g = undef
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     The namelist print entry has been matched with the var_table entry.  Now !
            ! lets cycle through our machines and determine if those machines hold data    !
            ! that should be printed.  If so, then send that data to the master (last      !
            ! machine).                                                                    !
            !------------------------------------------------------------------------------!
            !----- This is scratch space that all machines will use. ----------------------!
            pvar_l = undef
            !----- Set the blocking recieve to allow ordering, start with machine 1. ------!
            if (mynum /= 1) then
               call MPI_Recv(ping,1,MPI_INTEGER,recvnum,93,MPI_COMM_WORLD,status,ierr)
            end if

            !------------------------------------------------------------------------------!
            !     Cycle through this node's pointers for the current variable.  If the     !
            ! index falls within the printable range. Save that pointer to a local array.  !
            ! Once all the pointers have been cycled, send the local array to the master   !
            ! to populate a global array and print.                                        !
            !------------------------------------------------------------------------------!
            ptr_send   = .false.
            node_idmin = -1
            node_idmax = -1
            do np = 1,vt_info(nv,ifm)%nptrs
               g_id = vt_info(nv,ifm)%vt_vector(np)%globid+1
               g_ln = vt_info(nv,ifm)%vt_vector(np)%varlen
               !---------------------------------------------------------------------------!
               !      Determine if any of this segment falls within the range desired for  !
               ! output, globid+1 is the global index.                                     !
               !---------------------------------------------------------------------------!
               if (g_id <= ipmax .and. g_id+g_ln-1 >= ipmin ) then

                  !---- OK, this segment is good, set the send flag. ----------------------!
                  ptr_send = .true.
                  !------------------------------------------------------------------------!
                  !     These are the indices of the data in the current segment to use    !
                  ! and the indices in the global array they will be sent to.              !
                  !------------------------------------------------------------------------!
                  if (g_id >= ipmin) then
                     l_idmin = 1
                     g_idmin = g_id - ipmin + 1
                  else
                     l_idmin = ipmin - g_id + 1
                     g_idmin = 1
                  end if

                  if (g_id+g_ln-1 < ipmax) then
                     l_idmax = g_ln
                     g_idmax = g_id + g_ln - ipmin
                  else
                     l_idmax = ipmax - g_id + 1
                     g_idmax = ipmax - ipmin + 1
                  end if
                  
                  !------------------------------------------------------------------------!
                  !     These should be the same size, if not stop.                        !
                  !------------------------------------------------------------------------!
                  if (l_idmax - l_idmin /= g_idmax - g_idmin ) then
                     write(unit=*,fmt='(a)'      ) '--------------------------------------'
                     write(unit=*,fmt='(a,1x,i5)') ' L_IDMIN  = ',l_idmin
                     write(unit=*,fmt='(a,1x,i5)') ' L_IDMAX  = ',l_idmax
                     write(unit=*,fmt='(a,1x,i5)') ' G_IDMIN  = ',g_idmin
                     write(unit=*,fmt='(a,1x,i5)') ' G_IDMAX  = ',g_idmax
                     write(unit=*,fmt='(a,1x,i5)') ' L_IDSIZE = ',l_idmax - l_idmin + 1
                     write(unit=*,fmt='(a,1x,i5)') ' G_IDSIZE = ',g_idmax - g_idmin + 1
                     write(unit=*,fmt='(a)'      ) '--------------------------------------'
                     call fatal_error('L_IDSIZE and G_IDSIZE don''t match!'                &
                                     ,'print_fields','ed_print.f90')
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Shift the global dataset so that it is applied to the first index. !
                  ! Note that only real arrays can be used this time.                      !
                  !------------------------------------------------------------------------!
                  select case (vt_info(nv,ifm)%dtype)
                  case ('R')
                     call fillvar_l(pvar_l,vt_info(nv,ifm)%vt_vector(np)%var_rp,npolys     &
                                   ,g_ln,g_idmin,l_idmin,l_idmax)
                  case default
                     write(unit=*,fmt='(a)')      '---------------------------------------'
                     write(unit=*,fmt='(a,1x,a)') ' VNAME: ',vt_info(nv,ifm)%name
                     write(unit=*,fmt='(a,1x,a)') ' DTYPE: ',vt_info(nv,ifm)%dtype
                     write(unit=*,fmt='(a)')      '---------------------------------------'
                     call fatal_error('Only variables of DTYPE=R can be printed...'        &
                                     ,'print_fields','ed_print.f90')
                  end select
                  !------------------------------------------------------------------------!

                  !----- Determine the minimum and maximum indices that will be sent. -----!
                  if (g_idmin < node_idmin .or. node_idmin == -1) node_idmin = g_idmin
                  if (g_idmax > node_idmax .or. node_idmax == -1) node_idmax = g_idmax
                  !------------------------------------------------------------------------!
               end if
            end do

            !------------------------------------------------------------------------------!
            !     Decide what to do based on whether this is a serial or parallel run.     !
            !------------------------------------------------------------------------------!
            if (nnodetot > 1) then

               !------ Parallel run, decide whether this is a master or a slave node. -----!
               if (mynum /= nnodetot) then
                  !------------------------------------------------------------------------!
                  !      Slave node.  The local array for this machine has been created.   !
                  ! Send it off to the master.                                             !
                  !------------------------------------------------------------------------!

                  call MPI_Send(ptr_send,1,MPI_LOGICAL,machs(nnodetot),120,MPI_COMM_WORLD  &
                               ,ierr)

                  if (ptr_send) then
                     call MPI_Send(node_idmin,1,MPI_INTEGER,machs(nnodetot),121            &
                                  ,MPI_COMM_WORLD,ierr)

                     call MPI_Send(node_idmax,1,MPI_INTEGER,machs(nnodetot),122            &
                                  ,MPI_COMM_WORLD,ierr)

                     call MPI_Send(pvar_l(node_idmin:node_idmax),node_idmax-node_idmin+1   &
                                  ,MPI_REAL,machs(nnodetot),123,MPI_COMM_WORLD,ierr)
                  end if
                  

                  !------------------------------------------------------------------------!
                  !     When this node is finished, send the blocking MPI_Send to the next !
                  ! machine.                                                               !
                  !------------------------------------------------------------------------!
                  call MPI_Send(ping,1,MPI_INTEGER,sendnum,93,MPI_COMM_WORLD,ierr)

               else
                  !----- Master node, just copy the local array to the global. ------------!
                  if (ptr_send) then
                     pvar_g(node_idmin:node_idmax) = pvar_l(node_idmin:node_idmax)
                  end if
               end if
            else
               !----- Serial run, simply copy the entire local array to the global. -------!
               pvar_g = pvar_l
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      The data over the sought range of indices have been collected.  If this !
            ! is the only machine or the root machine, then print it to standard output.   !
            !------------------------------------------------------------------------------!
            if (mynum == nnodetot .or. nnodetot == 1) then

               !----- Print out a maximum of 10 variables per row... ----------------------!
               maxrows = ceiling(real(npolys)/real(maxcols))

               do row = 1,maxrows
                  ncols = min( maxcols,npolys-((row-1)*maxcols)   )
                  col   = ( (row-1)*maxcols)+1
                  
                  write(fmtstr,'(i3)') ncols
                  fmtstr = '('// trim(fmtstr)  // '(2x,' // trim(pfmtstr(ip)) // '))'
                  write(unit=*,fmt=trim(fmtstr)) (pvar_g(i),i=col,col+ncols-1)
               end do
               write (unit=*,fmt=*) ''
               write (unit=*,fmt=*) ''
            end if
         end if
      end do

      !----- Check to see if we matched the variable. -------------------------------------!
      if (.not.pvartrue) then
         write (unit=*,fmt='(a)')      '--------------------------------------------------'
         write (unit=*,fmt='(a,1x,a)') ' SOUGHT VARIABLE : ',trim(pvar_name)
         write (unit=*,fmt='(a)')      '  '
         write (unit=*,fmt='(a)')      ' The diagnostic variable doesn''t match any of the'
         write (unit=*,fmt='(a)')      ' var_table variables.  Check your namelist '
         write (unit=*,fmt='(a)')      ' entries, and the variable registry and/or remove'
         write (unit=*,fmt='(a)')      ' this diagnostic variable.'
         call fatal_error('Bad variable name.','print_fields','ed_print.f90')
      end if
   end do count_pvars
   !---------------------------------------------------------------------------------------!



   !----- Don't proceed until everything is written out. ----------------------------------!
   if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
   !---------------------------------------------------------------------------------------!

   return
end subroutine print_fields
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine simply copy the variable that is stored in the main variable struc- !
! ture to the scratch array that will be exchanged.                                        !
!------------------------------------------------------------------------------------------!
subroutine fillvar_l(pvar_l,vt_ptr,npts_out,npts_in,out1,in1,in2)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, dimension(npts_out), intent(inout) :: pvar_l
   real, dimension(npts_in) , intent(in)    :: vt_ptr
   integer                  , intent(in)    :: npts_in
   integer                  , intent(in)    :: npts_out
   integer                  , intent(in)    :: out1
   integer                  , intent(in)    :: in1
   integer                  , intent(in)    :: in2
   !----- Local variables. ----------------------------------------------------------------!
   integer :: i,j
   !---------------------------------------------------------------------------------------!


   j = out1
   do i = in1,in2     
      pvar_l(j) = vt_ptr(i)
      j = j + 1
   end do

   return
end subroutine fillvar_l
!==========================================================================================!
!==========================================================================================!
