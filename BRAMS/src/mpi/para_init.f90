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
!      This sub-routine is the main driver for the domain decomposition.                   !
!------------------------------------------------------------------------------------------!
subroutine node_decomp(init)

   use mem_grid
   use var_tables
   use mem_scratch
   use mem_basic
   use rpara
   use mem_cuparm, only : nclouds
   use grid_dims , only : ndim_types
   use mem_aerad , only : nwave

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical                       , intent(in) :: init
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(ndim_types)             :: npvar
   integer                                    :: ngr
   integer                                    :: idn
   integer                                    :: isn
   integer                                    :: ng
   integer                                    :: mp_nzp
   integer                                    :: numbuff
   integer                                    :: icm
   integer                                    :: nsiz
   integer                                    :: jnode
   integer                                    :: ncols
   integer                                    :: ifm
   integer                                    :: nestvar
   integer                                    :: nc
   integer                                    :: nf
   integer                                    :: nv
   integer                                    :: num_lbc_buff
   integer                                    :: num_nest_buff
   integer                                    :: num_six_buff
   integer                                    :: num_feed_buff
   integer                                    :: itype
   integer                                    :: i1
   integer                                    :: i2
   integer                                    :: j1
   integer                                    :: j2
   integer                                    :: ixy
   integer                                    :: ixyz
   integer                                    :: fdzp
   integer                                    :: fdep
   integer                                    :: memf
   integer                                    :: idim
   logical                                    :: failed
   !---------------------------------------------------------------------------------------!



   !----- Try to read domain decomposition data from input file. --------------------------!
   call PAR_decomp_input(maxmach,maxgrds,nmachs,ngrids,nnxp,nnyp,ixb,ixe,iyb,iye,failed)
   !---------------------------------------------------------------------------------------!



   !----- Decompose all grids into subdomains. --------------------------------------------!
   do ngr = 1,ngrids
      !----- Decompose only if domain decomposition couldn't be specified by input file. --!
      if (failed) then
         !---------------------------------------------------------------------------------!
         !      Obtain estimates of the fraction of computational time (work) required for !
         ! each column in the region of the domain.                                        !
         !---------------------------------------------------------------------------------!
         call par_est_time(nnxp(ngr),nnyp(ngr),scratch%scr1,basic_g(ngr)%cputime,init)

         !----- Decompose the grid taking into account the work numbers. ------------------!
         nsiz = nnxp(ngr)+nnyp(ngr)
         call par_decomp(nnxp(ngr),nnyp(ngr),nsiz,nmachs,scratch%scr1                      &
                        ,ixb(:,ngr),ixe(:,ngr),iyb(:,ngr),iye(:,ngr))
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !----- Compute various bounds for the sub-domains. -------------------------------------!
   call par_decomp_bounds(1,1)
   !---------------------------------------------------------------------------------------!



   !----- Determine node sending paths and numbers of receiving nodes. --------------------!
   call par_node_paths()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute  send and receive buffer sizes. These will be maximum of long timestep,   !
   ! turbulence, nest boundaries, and nest feedback.  Small timestep will use same buffers !
   ! as they are always smaller.                                                           !
   !---------------------------------------------------------------------------------------!
   !----- Initialise the arrays with the size for Lateral Boundary Condition buffers. -----!
   do idn=1,nmachs
      do isn=1,nmachs
         lbc_buffs(1,idn,isn)=0
         lbc_buffs(2,idn,isn)=0
      end do
   end do
   !----- Find the maximum vertical size. -------------------------------------------------!
   mp_nzp = maxval(nnzp(1:ngrids))
   !---------------------------------------------------------------------------------------!
   !     Grid loop.                                                                        !
   !---------------------------------------------------------------------------------------!
   gridloop: do ng=1,ngrids

      !----- Find number of nested variables to be communicated. --------------------------!
      icm = ng
      ifm = ng
      if (ng /= 1) icm = nxtnest(ifm)
      nestvar = 4
      do nf=1,num_scalar(ifm)
         do nc=1,num_scalar(icm)
            if (scalar_tab(nf,ifm)%name==scalar_tab(nc,icm)%name) nestvar=nestvar+1
         end do
      end do

      !---- Find number of LBC variables to be communicated. ------------------------------!
      npvar(:) = 0
      do nv = 1,num_var(ng)
         if (vtab_r(nv,ng)%impt1 == 1 ) then
            idim        = vtab_r(nv,ng)%idim_type
            npvar(idim) = npvar(idim) + 1
         end if
      end do


      sourcemach1: do isn=1,nmachs
         destmach1: do idn=1,nmachs
            num_lbc_buff  = 0
            num_nest_buff = 0
            num_six_buff  = 0
            num_feed_buff = 0

            itype=1
            i1=inode_paths_master(1,itype,ng,idn,isn)
            i2=inode_paths_master(2,itype,ng,idn,isn)
            j1=inode_paths_master(3,itype,ng,idn,isn)
            j2=inode_paths_master(4,itype,ng,idn,isn)

            !----- If i1 is not 0, this means that the two nodes are direct neighbours. ---!
            if (i1 /= 0) then

               ixy = (i2-i1+1)*(j2-j1+1)

               !------ Add the total number of points to be sent. -------------------------!
               num_lbc_buff = 0
               do idim =2,ndim_types
                  call ze_dims(ng,idim,.true.,fdzp,fdep)
                  num_lbc_buff = num_lbc_buff + ixy * fdzp * fdep * npvar(idim)
               end do
               
               num_lbc_buff = num_lbc_buff + 2 *(sum(npvar(:)) + 100)
            end if
            !------------------------------------------------------------------------------!



            !----- Find the size of the nesting exchange. ---------------------------------!
            itype = 5
            i1 = inode_paths_master(1,itype,ng,idn,isn)
            i2 = inode_paths_master(2,itype,ng,idn,isn)
            j1 = inode_paths_master(3,itype,ng,idn,isn)
            j2 = inode_paths_master(4,itype,ng,idn,isn)
            if (i1 /= 0) then
               ixyz          = (i2-i1+1) * (j2-j1+1) * mp_nzp
               num_nest_buff = ixyz * nestvar + 2*(nestvar+100)
            end if
            !------------------------------------------------------------------------------!



            !----- Find the size of the nesting exchange. ---------------------------------!
            itype = 6
            i1 = inode_paths_master(1,itype,ng,idn,isn)
            i2 = inode_paths_master(2,itype,ng,idn,isn)
            j1 = inode_paths_master(3,itype,ng,idn,isn)
            j2 = inode_paths_master(4,itype,ng,idn,isn)
            if (i1 /= 0) then
               ixyz          = (i2-i1+1) * (j2-j1+1) * mp_nzp
               num_six_buff  = ixyz * nestvar + 2*(nestvar+100)
            end if
            !------------------------------------------------------------------------------!



            !----- Find the size of the nesting exchange. ---------------------------------!
            itype = 7
            i1 = inode_paths_master(1,itype,ng,idn,isn)
            i2 = inode_paths_master(2,itype,ng,idn,isn)
            j1 = inode_paths_master(3,itype,ng,idn,isn)
            j2 = inode_paths_master(4,itype,ng,idn,isn)
            if (i1 /= 0) then
               ixyz          = (i2-i1+1) * (j2-j1+1) * mp_nzp
               num_feed_buff = ixyz * nestvar + 2*(nestvar+100)
            end if
            !------------------------------------------------------------------------------!



            !------ Update the buffer size with the highest value so far. -----------------!
            lbc_buffs(1,idn,isn) = max( lbc_buffs(1,idn,isn)                               &
                                      , num_lbc_buff                                       &
                                      , num_nest_buff                                      &
                                      , num_six_buff                                       &
                                      , num_feed_buff)
            !------------------------------------------------------------------------------!
         end do destmach1
      end do sourcemach1



      sourcemach2: do isn=1,nmachs
         destmach2: do idn=1,nmachs
            num_lbc_buff  = 0
            num_nest_buff = 0
            num_six_buff  = 0
            num_feed_buff = 0

            itype = 1
            i1 = inode_paths_master(1,itype,ng,isn,idn)
            i2 = inode_paths_master(2,itype,ng,isn,idn)
            j1 = inode_paths_master(3,itype,ng,isn,idn)
            j2 = inode_paths_master(4,itype,ng,isn,idn)
            if (i1 /= 0) then
               ixy = (i2-i1+1)*(j2-j1+1)

               num_lbc_buff = 0
               do idim =2,ndim_types
                  call ze_dims(ng,idim,.true.,fdzp,fdep)
                  num_lbc_buff = num_lbc_buff + ixy * fdzp * fdep * npvar(idim)
               end do
               num_lbc_buff = num_lbc_buff + 2*(sum(npvar(:)) + 100)
            end if

            itype = 5
            i1 = inode_paths_master(1,itype,ng,isn,idn)
            i2 = inode_paths_master(2,itype,ng,isn,idn)
            j1 = inode_paths_master(3,itype,ng,isn,idn)
            j2 = inode_paths_master(4,itype,ng,isn,idn)
            if(i1 /= 0) then
               ixyz          = (i2-i1+1) * (j2-j1+1) * mp_nzp
               num_nest_buff = ixyz * nestvar + 2*(nestvar+100)
            end if

            itype = 6
            i1 = inode_paths_master(1,itype,ng,isn,idn)
            i2 = inode_paths_master(2,itype,ng,isn,idn)
            j1 = inode_paths_master(3,itype,ng,isn,idn)
            j2 = inode_paths_master(4,itype,ng,isn,idn)
            if(i1 /= 0) then
               ixyz          = (i2-i1+1) * (j2-j1+1) * mp_nzp
               num_six_buff  = ixyz * nestvar + 2*(nestvar+100)
            end if

            itype = 7
            i1 = inode_paths_master(1,itype,ng,isn,idn)
            i2 = inode_paths_master(2,itype,ng,isn,idn)
            j1 = inode_paths_master(3,itype,ng,isn,idn)
            j2 = inode_paths_master(4,itype,ng,isn,idn)
            if(i1 /= 0) then
               ixyz          = (i2-i1+1) * (j2-j1+1) * mp_nzp
               num_feed_buff = ixyz * nestvar + 2*(nestvar+100)
            end if

            lbc_buffs(2,idn,isn) = max( lbc_buffs(2,idn,isn)                               &
                                      , num_lbc_buff                                       &
                                      , num_nest_buff                                      &
                                      , num_six_buff                                       &
                                      , num_feed_buff)
         end do destmach2
      end do sourcemach2
   end do gridloop

   !------ Check nest boundary receive buffer size. ---------------------------------------!
   itype=5
   do idn=1,nmachs
      newbuff_nest1(idn)=1
      nbuff_nest1(idn)=0

      do ng=1,ngrids
         numbuff = 0
         icm     = nxtnest(ng)
         do isn=1,nmachs
            i1      = inode_paths_master(1,itype,ng,idn,isn)
            i2      = inode_paths_master(2,itype,ng,idn,isn)
            j1      = inode_paths_master(3,itype,ng,idn,isn)
            j2      = inode_paths_master(4,itype,ng,idn,isn)
            memf    = (i2-i1+1) * (j2-j1+1) * (mp_nzp) * nestvar
            numbuff = numbuff + memf
         end do
         nbuff_nest1(idn) = max(nbuff_nest1(idn),numbuff)
      end do
   end do


   return
end subroutine node_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine defines the domain in the case in which only one node is used.       !
!------------------------------------------------------------------------------------------!
subroutine onenode()
   use mem_grid
   use node_mod

   implicit none

   mmzp = nnzp
   mmxp = nnxp
   mmyp = nnyp

   mxp = nnxp(ngrid)
   myp = nnyp(ngrid)
   mzp = nnzp(ngrid)

   ia   = 2
   iz   = mxp-1
   ia_1 = max(ia-1,1)
   ia_2 = max(ia-2,1)
   ia_3 = max(ia-3,1)
   ia1  = ia+1
   ia2  = ia+2
   ia3  = ia+3
   iz_1 = iz-1
   iz_2 = iz-2
   iz_3 = iz-3
   iz1  = min(iz+1,mxp)
   iz2  = min(iz+2,mxp)
   iz3  = min(iz+3,mxp)
   izu  = iz-1

   if (nnyp(ngrid) > 1) then
      ja  = 2
      jz  = myp-1
      ja_1= max(ja-1,1)
      ja_2= max(ja-2,1)
      ja_3= max(ja-3,1)
      ja1 = ja+1
      ja2 = ja+2
      ja3 = ja+3
      jz_1= jz-1
      jz_2= jz-2
      jz_3= jz-3
      jz1 = min(jz+1,myp)
      jz2 = min(jz+2,myp)
      jz3 = min(jz+3,myp)
      jzv = jz-1
   else
      ja   = 1
      jz   = 1
      ja_1 = 1
      ja_2 = 1
      ja_3 = 1
      ja1  = 1
      ja2  = 1
      ja3  = 1
      jz_1 = 1
      jz_2 = 1
      jz_3 = 1
      jz1  = 1
      jz2  = 1
      jz3  = 1
      jzv  = 1
   end if

   i0    = 0
   j0    = 0
   ibcon = 1+2+4+8
   ipara = 0

   mi0(1:ngrids) = i0
   mj0(1:ngrids) = j0
   mia(1:ngrids) = ia
   miz(1:ngrids) = mmxp(1:ngrids)-1
   mja(1:ngrids) = ja
   mjz(1:ngrids) = mmyp(1:ngrids)-1



   return
end subroutine onenode
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     Compute various subdomain boundary numbers for the nodes                             !
!  - NXBEG, NYBEG, NXEND, NYEND     - portions of full domain that nodes will have         !
!                                     (includes overlap region and absolute edges)         !
!                                                                                          !
!  - IXOFF, IYOFF                   - subdomain offsets relative to full domain            !
!                                                                                          !
!  - NXBEGC, NXENDC, NYBEGC, NYENDC - subdomain "compute" points, or normal thermodynamic  !
!                                     tendency points (2-nx, 2-ny for non-parallel run)    !
!                                                                                          !
!  - IBCFLAG                        - flag denoting if real boundary is on subdomain       !
!                                     bit 1=west, bit 2=east, bit 3=south, bit 4=north     !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine par_decomp_bounds(nbndx,nbndy)
   use mem_grid, only : ngrids & ! intent(in)
                      , nnxp   & ! intent(in)
                      , nnyp   ! ! intent(in)
   use rpara

   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                   , intent(in) :: nbndx
   integer                   , intent(in) :: nbndy
   !------ Local variables. ---------------------------------------------------------------!
   integer                                :: ng
   integer                                :: nm
   integer                                :: nx
   integer                                :: ny
   !---------------------------------------------------------------------------------------!


   do ng=1,ngrids
      !----- Copy current grid information to a scalar alias. -----------------------------!
      nx=nnxp(ng)
      ny=nnyp(ng)
      !------------------------------------------------------------------------------------!


      do nm=1,nmachs

         !----- Initialise the boundary flag with zero (inner domain is default). ---------!
         ibcflg(nm,ng) = 0

         !---------------------------------------------------------------------------------!
         !     Check western side.                                                         !
         !---------------------------------------------------------------------------------!
         if (ixb(nm,ng) == 2) then
            !----- Sub-domain is in the western edge. -------------------------------------!
            nxbeg(nm,ng)  = 1
            ixoff(nm,ng)  = 0
            nxbegc(nm,ng) = 2
            ibcflg(nm,ng) = ibcflg(nm,ng) + 1
         else
            !----- Western edge of the sub-domain is not in the absolute edge. ------------!
            nxbeg(nm,ng)  = ixb(nm,ng) - nbndx
            ixoff(nm,ng)  = nxbeg(nm,ng) - 1
            nxbegc(nm,ng) = nbndx + 1
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check eastern side.                                                         !
         !---------------------------------------------------------------------------------!
         if(ixe(nm,ng) == nx-1) then
            !----- Sub-domain is in the eastern edge. -------------------------------------!
            nxend(nm,ng)  = nx
            nxendc(nm,ng) = nxbegc(nm,ng) + ixe(nm,ng) - ixb(nm,ng)
            ibcflg(nm,ng) = ibcflg(nm,ng) + 2
         else
            !----- Eastern edge of the sub-domain is not in the absolute edge. ------------!
            nxend(nm,ng)  = ixe(nm,ng) + nbndx
            nxendc(nm,ng) = nxbegc(nm,ng) + ixe(nm,ng) - ixb(nm,ng)
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check southern side.                                                        !
         !---------------------------------------------------------------------------------!
         if (iyb(nm,ng) == 2) then
            !----- Sub-domain is in the southern edge. ------------------------------------!
            nybeg(nm,ng)  = 1
            iyoff(nm,ng)  = 0
            nybegc(nm,ng) = 2
            ibcflg(nm,ng) = ibcflg(nm,ng) + 4
         else
            !----- Southern edge of the sub-domain is not in the absolute edge. -----------!
            nybeg(nm,ng)  = iyb(nm,ng) - nbndy
            iyoff(nm,ng)  = nybeg(nm,ng) - 1
            nybegc(nm,ng) = nbndy + 1
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check northern side.                                                        !
         !---------------------------------------------------------------------------------!
         if(iye(nm,ng) == ny-1) then
            nyend(nm,ng)  = ny
            ibcflg(nm,ng) = ibcflg(nm,ng) + 8
            nyendc(nm,ng) = nybegc(nm,ng) + iye(nm,ng) - iyb(nm,ng)
         else
            nyend(nm,ng)  = iye(nm,ng)+nbndy
            nyendc(nm,ng) = nybegc(nm,ng) + iye(nm,ng) - iyb(nm,ng)
         end if
      end do

      !----- Find the total number of horizontal points. ----------------------------------!
      do nm=1,nmachs
         npxy(nm,ng)= (nxend(nm,ng)-nxbeg(nm,ng)+1) * (nyend(nm,ng)-nybeg(nm,ng)+1)
      end do
      !------------------------------------------------------------------------------------!
   end do

   return
end subroutine par_decomp_bounds
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine dump_domain_decomposition()
   use mem_grid, only : ngrids & ! intent(in)
                      , grid_g ! ! intent(in)
   use rpara   , only : nmachs & ! intent(in)
                      , ixb    & ! intent(in)
                      , ixe    & ! intent(in)
                      , iyb    & ! intent(in)
                      , iye    ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: ng
   integer :: nm
   integer :: npts
   integer :: i
   integer :: j
   !---------------------------------------------------------------------------------------!
   write(unit=*,fmt='(a)')       '================= Domain decomposition ================='
   write(unit=*,fmt='(7(a,1x))') '  GRID','  NODE',' X_BEG',' X_END',' Y_BEG',' Y_END'     &
                                ,'  NPTS'
   do ng = 1, ngrids
      do nm = 1,nmachs
         !----- Define the output variable with the node number. --------------------------!
         do i=ixb(nm,ng),ixe(nm,ng)
            do j=iyb(nm,ng),iye(nm,ng)
               grid_g(ng)%fmynum(i,j) = real(nm)
            end do
         end do

         !----- Find the total number of points. ------------------------------------------!
         npts = (ixe(nm,ng)-ixb(nm,ng) + 1) * (iye(nm,ng)-iyb(nm,ng) + 1)

         !----- Print the table on screen. ------------------------------------------------!
         write (unit=*,fmt='(7(i6,1x))') ng,nm,ixb(nm,ng),ixe(nm,ng),iyb(nm,ng),iye(nm,ng) &
                                        ,npts
      end do
      write(unit=*,fmt='(a)') ' '
   end do
   write(unit=*,fmt='(a)')       '========================================================'
   write(unit=*,fmt='(a)')       ' '

   return
end subroutine dump_domain_decomposition
!==========================================================================================!
!==========================================================================================!


