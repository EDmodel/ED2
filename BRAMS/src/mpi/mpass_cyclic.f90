!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################

subroutine ipaths_cyc_alloc(nxp,nyp,ibnd,jbnd)
  use cyclic_mod
  implicit none
  integer :: nxp,nyp,ibnd,jbnd

  npts_cyc = 8 * maxmach
  if (ibnd == 4) npts_cyc = npts_cyc + nyp * 4
  if (jbnd == 4) npts_cyc = npts_cyc + nxp * 4
  if (ibnd == 4 .and. jbnd == 4) npts_cyc = npts_cyc + 16

  allocate (ipathst_cyc(8,npts_cyc))
  allocate (ipathsu_cyc(8,npts_cyc))
  allocate (ipathsv_cyc(8,npts_cyc))

  call izero(8*npts_cyc,ipathst_cyc)
  call izero(8*npts_cyc,ipathsu_cyc)
  call izero(8*npts_cyc,ipathsv_cyc)


  return
end subroutine ipaths_cyc_alloc

!*****************************************************************************

subroutine node_cycinit(nzp,nxp,nyp,npvar,nmachs,ibnd,jbnd,mynum)
  use cyclic_mod
  use mem_grid, only : npatch, nzg, nzs
  use mem_cuparm, only : nclouds
  use mem_aerad , only : nwave
  use grid_dims, only : ndim_types
  implicit none
  include 'mpif.h'
  integer :: nmachs,mynum,icypts,nzp,nxp,nyp,icycpts,iadd,mdn,msn,ndn,nsn  &
     ,ibnd,jbnd,nm,maxijrecv_cyc
  integer, dimension(ndim_types) :: npvar

  integer, save, allocatable :: ijcount(:),jdn(:),jsn(:)

  integer :: idn,isn,ii,iijj
  integer :: ierr
  integer :: float_size
  integer :: hugedim
  integer :: npvar_sum
  

  ! This subroutine is called by each node process.  It uses information from
  ! the master cyclic parallel array [ipaths_cyc] to construct integer scalars 
  ! and arrays that the node needs for parallel data sends and receives for 
  ! cyclic boundary conditions

  npvar_sum = sum(npvar)

  allocate (ndn_cyc(nmachs),nsn_cyc(nmachs),msn_cyc(nmachs),mdn_cyc(nmachs)  &
     ,nijsendt_cyc(nmachs),nijsendu_cyc(nmachs),nijsendv_cyc(nmachs)  &
     ,nijrecvt_cyc(nmachs),nijrecvu_cyc(nmachs),nijrecvv_cyc(nmachs)  &
     ,ijcount(nmachs)  &
     ,jdn(nmachs),jsn(nmachs),isend_req_cyc(6,nmachs),irecv_req_cyc(6,nmachs))

  do nm = 1,nmachs
     nijsendt_cyc(nm) = 0.
     nijsendu_cyc(nm) = 0.
     nijsendv_cyc(nm) = 0.
     nijrecvt_cyc(nm) = 0.
     nijrecvu_cyc(nm) = 0.
     nijrecvv_cyc(nm) = 0.
     isend_req_cyc(1:6,nm)=-999
     irecv_req_cyc(1:6,nm)=-999
  enddo

  do icycpts = 1,npts_cyc

     msn = ipathst_cyc(1,icycpts)
     mdn = ipathst_cyc(4,icycpts)

     if (msn .eq. mynum .and. mdn .gt. 0) then
        nijsendt_cyc(mdn) = nijsendt_cyc(mdn) + 1
     endif
     
     if (mdn .eq. mynum .and. msn .gt. 0) then
        nijrecvt_cyc(msn) = nijrecvt_cyc(msn) + 1
     endif
  enddo

  do icycpts = 1,npts_cyc
     msn = ipathsu_cyc(1,icycpts)
     mdn = ipathsu_cyc(4,icycpts)

     if (msn .eq. mynum .and. mdn .gt. 0) then
        nijsendu_cyc(mdn) = nijsendu_cyc(mdn) + 1
     endif
     if (mdn .eq. mynum .and. msn .gt. 0) then
        nijrecvu_cyc(msn) = nijrecvu_cyc(msn) + 1
     endif
  enddo

  do icycpts = 1,npts_cyc
     msn = ipathsv_cyc(1,icycpts)
     mdn = ipathsv_cyc(4,icycpts)

     if (msn .eq. mynum .and. mdn .gt. 0) then
        nijsendv_cyc(mdn) = nijsendv_cyc(mdn) + 1
     endif
     if (mdn .eq. mynum .and. msn .gt. 0) then
        nijrecvv_cyc(msn) = nijrecvv_cyc(msn) + 1
     endif

  enddo

  nsn = 0
  ndn = 0
  maxijsendt_cyc = 0
  maxijsendu_cyc = 0
  maxijsendv_cyc = 0
  maxijrecvt_cyc = 0
  maxijrecvu_cyc = 0
  maxijrecvv_cyc = 0

  do nm = 1,nmachs
     if (nijsendt_cyc(nm) .gt. 0) then
        ndn = ndn + 1
        ndn_cyc(nm) = ndn
        mdn_cyc(ndn) = nm
        nijsendt_cyc(ndn) = nijsendt_cyc(nm)
        nijsendu_cyc(ndn) = nijsendu_cyc(nm)
        nijsendv_cyc(ndn) = nijsendv_cyc(nm)
        maxijsendt_cyc = max(maxijsendt_cyc,nijsendt_cyc(ndn))
        maxijsendu_cyc = max(maxijsendu_cyc,nijsendu_cyc(ndn))
        maxijsendv_cyc = max(maxijsendv_cyc,nijsendv_cyc(ndn))
     endif

     if (nijrecvt_cyc(nm) .gt. 0) then
        nsn = nsn + 1
        nsn_cyc(nm) = nsn
        msn_cyc(nsn) = nm
        nijrecvt_cyc(nsn) = nijrecvt_cyc(nm)
        nijrecvu_cyc(nsn) = nijrecvu_cyc(nm)
        nijrecvv_cyc(nsn) = nijrecvv_cyc(nm)
        maxijrecvt_cyc = max(maxijrecvt_cyc,nijrecvt_cyc(nsn))
        maxijrecvu_cyc = max(maxijrecvu_cyc,nijrecvu_cyc(nsn))
        maxijrecvv_cyc = max(maxijrecvv_cyc,nijrecvv_cyc(nsn))
     endif
  enddo

  ndns_cyc = ndn
  nsns_cyc = nsn

  maxijrecv_cyc = max(maxijrecvt_cyc,maxijrecvu_cyc,maxijrecvv_cyc)

  allocate (ijsendt_cyc(6,maxijsendt_cyc,ndns_cyc) &
           ,ijsendu_cyc(6,maxijsendu_cyc,ndns_cyc) &
           ,ijsendv_cyc(6,maxijsendv_cyc,ndns_cyc) &
           ,ijrecv_cyc(6,maxijrecv_cyc))

  do ndn = 1,ndns_cyc
     ijcount(ndn) = 0
  enddo

  do icycpts = 1,npts_cyc
     msn = ipathst_cyc(1,icycpts)
     mdn = ipathst_cyc(4,icycpts)
     if (msn .eq. mynum .and. mdn .gt. 0) then
        ndn = ndn_cyc(mdn)
        ijcount(ndn) = ijcount(ndn) + 1
        ijsendt_cyc(1,ijcount(ndn),ndn) = ipathst_cyc(2,icycpts)
        ijsendt_cyc(2,ijcount(ndn),ndn) = ipathst_cyc(3,icycpts)
        ijsendt_cyc(3,ijcount(ndn),ndn) = ipathst_cyc(5,icycpts)
        ijsendt_cyc(4,ijcount(ndn),ndn) = ipathst_cyc(6,icycpts)
        ijsendt_cyc(5,ijcount(ndn),ndn) = ipathst_cyc(7,icycpts)
        ijsendt_cyc(6,ijcount(ndn),ndn) = ipathst_cyc(8,icycpts)
     endif
     
  enddo

  do ndn = 1,ndns_cyc
     ijcount(ndn) = 0
  enddo

  do icycpts = 1,npts_cyc
     msn = ipathsu_cyc(1,icycpts)
     mdn = ipathsu_cyc(4,icycpts)
     if (msn == mynum .and. mdn > 0) then
        ndn = ndn_cyc(mdn)
        ijcount(ndn) = ijcount(ndn) + 1
        ijsendu_cyc(1,ijcount(ndn),ndn) = ipathsu_cyc(2,icycpts)
        ijsendu_cyc(2,ijcount(ndn),ndn) = ipathsu_cyc(3,icycpts)
        ijsendu_cyc(3,ijcount(ndn),ndn) = ipathsu_cyc(5,icycpts)
        ijsendu_cyc(4,ijcount(ndn),ndn) = ipathsu_cyc(6,icycpts)
        ijsendu_cyc(5,ijcount(ndn),ndn) = ipathsu_cyc(7,icycpts)
        ijsendu_cyc(6,ijcount(ndn),ndn) = ipathsu_cyc(8,icycpts)
     endif
  enddo

  do ndn = 1,ndns_cyc
     ijcount(ndn) = 0
  enddo

  do icycpts = 1,npts_cyc
     msn = ipathsv_cyc(1,icycpts)
     mdn = ipathsv_cyc(4,icycpts)
     if (msn == mynum .and. mdn > 0) then
        ndn = ndn_cyc(mdn)
        ijcount(ndn) = ijcount(ndn) + 1
        ijsendv_cyc(1,ijcount(ndn),ndn) = ipathsv_cyc(2,icycpts)
        ijsendv_cyc(2,ijcount(ndn),ndn) = ipathsv_cyc(3,icycpts)
        ijsendv_cyc(3,ijcount(ndn),ndn) = ipathsv_cyc(5,icycpts)
        ijsendv_cyc(4,ijcount(ndn),ndn) = ipathsv_cyc(6,icycpts)
        ijsendv_cyc(5,ijcount(ndn),ndn) = ipathsv_cyc(7,icycpts)
        ijsendv_cyc(6,ijcount(ndn),ndn) = ipathsv_cyc(8,icycpts)
     endif
  enddo

  ! allocate cyclic buffers
  call MPI_Pack_size(1,MPI_REAL,MPI_COMM_WORLD,float_size,ierr)

  !---- Define the maximum size for cyclic buffers.
  hugedim = max(nzp,nzg,nzs)*max(npatch,nwave,nclouds)

  nbuffsend_cyc = 10+max(maxijsendt_cyc * (max(2,npvar_sum) * hugedim  + 6) + 1  &
                     ,maxijsendu_cyc * (hugedim + 6) + 1                 &
                    + maxijsendv_cyc * (hugedim + 6) + 1                 ) 
  nbuffsend_cyc = nbuffsend_cyc * float_size

  nbuffrecv_cyc = 10+max(maxijrecvt_cyc * (max(2,npvar_sum) * hugedim + 6) + 1  &
                     ,maxijrecvu_cyc * (hugedim + 6) + 1                 &
                    + maxijrecvv_cyc * (hugedim + 6) + 1                 ) 
  nbuffrecv_cyc = nbuffrecv_cyc * float_size

  nbuffrecv_cyc = max(nbuffrecv_cyc, nbuffsend_cyc)
  nbuffsend_cyc = nbuffrecv_cyc

  do nm = 1,maxmach
     allocate (node_buffs_cyc(nm)%buffsend_cyc(nbuffsend_cyc))
     allocate (node_buffs_cyc(nm)%buffrecv_cyc(nbuffrecv_cyc))
  enddo

  if (ibnd .eq. 4 .or. jbnd .eq. 4) then
     allocate(lstart_cyc(npvar_sum+4,nxp,nyp))
  endif

  return
end subroutine node_cycinit
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine node_sendcyclic(isflag)
   use mem_grid
   use var_tables
   use cyclic_mod
   use mem_basic
   use mem_scratch
   use node_mod
   use grid_dims
   use mem_cuparm, only: nclouds
   use mem_aerad , only : nwave    ! ! intent(in)

   implicit none
   !----- Included variables --------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: isflag
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   integer             :: ipos
   integer             :: nsn
   integer             :: icycpts
   integer             :: msn
   integer             :: mdn
   integer             :: isn
   integer             :: jsn
   integer             :: ind
   integer             :: nmp
   integer             :: ndn
   integer             :: ibob
   integer             :: ii
   integer             :: nm
   integer             :: iijj
   integer             :: nv
   integer             :: ijr
   integer             :: recv_int
   integer             :: mpiid
   integer             :: fdzp
   integer             :: fdep
   !---------------------------------------------------------------------------------------!

   if (ibnd /= 4 .and. jbnd /= 4) return

   !----- First, before we send anything, let's post the receives -------------------------!

   do nsn = 1,nsns_cyc
      msn = msn_cyc(nsn)
      !----- Making a unique ID for this package ------------------------------------------!
      mpiid = 1000000000+maxmach*10*(machs(msn)-1)+10*(machs(mynum)-1)+ isflag
      call MPI_Irecv(node_buffs_cyc(nsn)%buffrecv_cyc, nbuffrecv_cyc,MPI_PACKED            &
                    ,machs(msn),mpiid,MPI_COMM_WORLD,irecv_req_cyc(isflag,msn),ierr)
   end do

   !----- Now we can actually go on to sending the stuff ----------------------------------!

   do ndn = 1,ndns_cyc

      mdn = mdn_cyc(ndn)

      ipos = 1
      select case (isflag)
      case (1,4,6)
         call MPI_Pack(nijsendt_cyc(ndn),1,MPI_INTEGER,node_buffs_cyc(ndn)%buffsend_cyc(1) &
                      ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(ijsendt_cyc(:,:,ndn),6*nijsendt_cyc(ndn),MPI_INTEGER                &
                      ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos,MPI_COMM_WORLD  &
                      ,ierr)

         do icycpts = 1,nijsendt_cyc(ndn)
            isn = ijsendt_cyc(1,icycpts,ndn) - mi0(1)
            jsn = ijsendt_cyc(2,icycpts,ndn) - mj0(1)

            select case (isflag)
            case (1)
               do nv = 1,num_var(1)
                  if ( vtab_r(nv,1)%impt1 == 1) then
                     !----- Find the variable dimensions in Z and E axes. -----------------!
                     call ze_dims(1,vtab_r(nv,1)%idim_type,.false.,fdzp,fdep)

                     !----- Copy the variable to a scratch array. -------------------------!
                     call mk_cyc_buff(fdzp,mmxp(1),mmyp(1),fdep,vtab_r(nv,1)%var_p         &
                                   ,scratch%scr1,isn,jsn)

                     !----- Pack the variable into the buffer. ----------------------------!
                     call MPI_Pack(scratch%scr1,fdzp*fdep,MPI_REAL                         & 
                                  ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos     &
                                  ,MPI_COMM_WORLD,ierr)
                  end if
               end do

            case (4)
               call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%pp,scratch%scr1       &
                             ,isn,jsn)
               call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL                                 &
                            ,node_buffs_cyc(ndn)%buffsend_cyc                              &
                            ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)

            case (6)
               call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%wp,scratch%scr1       &
                             ,isn,jsn)
               call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL                                 &
                            ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos           &
                            ,MPI_COMM_WORLD,ierr)
               call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%pp,scratch%scr1       &
                             ,isn,jsn)
               call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL                                 &
                            ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos           &
                            ,MPI_COMM_WORLD,ierr)
            end select
         end do



      case (2)
         call MPI_Pack(nijsendu_cyc(ndn),1,MPI_INTEGER,node_buffs_cyc(ndn)%buffsend_cyc    &
                      ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(ijsendu_cyc(:,:,ndn),6*nijsendu_cyc(ndn),MPI_INTEGER                &
                      ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos,MPI_COMM_WORLD  &
                      ,ierr)

         do icycpts = 1,nijsendu_cyc(ndn)
            isn = ijsendu_cyc(1,icycpts,ndn) - mi0(1)
            jsn = ijsendu_cyc(2,icycpts,ndn) - mj0(1)

            call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%up,scratch%scr1,isn,jsn)
            call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL,node_buffs_cyc(ndn)%buffsend_cyc   &
                         ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         end do

      case (3)
         call MPI_Pack(nijsendv_cyc(ndn),1,MPI_INTEGER,node_buffs_cyc(ndn)%buffsend_cyc    &
                      ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(ijsendv_cyc(:,:,ndn),6*nijsendv_cyc(ndn),MPI_INTEGER                &
                      ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos,MPI_COMM_WORLD  &
                      ,ierr)

         do icycpts = 1,nijsendv_cyc(ndn)
            isn = ijsendv_cyc(1,icycpts,ndn) - mi0(1)
            jsn = ijsendv_cyc(2,icycpts,ndn) - mj0(1)

            call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%vp,scratch%scr1,isn,jsn)
            call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL,node_buffs_cyc(ndn)%buffsend_cyc   &
                         ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         end do

      case (5)
         call MPI_Pack(nijsendu_cyc(ndn),1,MPI_INTEGER,node_buffs_cyc(ndn)%buffsend_cyc    &
                      ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(ijsendu_cyc(:,:,ndn),6*nijsendu_cyc(ndn),MPI_INTEGER                &
                      ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos,MPI_COMM_WORLD  &
                      ,ierr)

         do icycpts = 1,nijsendu_cyc(ndn)
            isn = ijsendu_cyc(1,icycpts,ndn) - mi0(1)
            jsn = ijsendu_cyc(2,icycpts,ndn) - mj0(1)

            call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%up,scratch%scr1,isn,jsn)
            call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL,node_buffs_cyc(ndn)%buffsend_cyc   &
                         ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         end do

         call MPI_Pack(nijsendv_cyc(ndn),1,MPI_INTEGER,node_buffs_cyc(ndn)%buffsend_cyc    &
                      ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(ijsendv_cyc(:,:,ndn),6*nijsendv_cyc(ndn),MPI_INTEGER                &
                      ,node_buffs_cyc(ndn)%buffsend_cyc,nbuffsend_cyc,ipos,MPI_COMM_WORLD  &
                      ,ierr)

         do icycpts = 1,nijsendv_cyc(ndn)
            isn = ijsendv_cyc(1,icycpts,ndn) - mi0(1)
            jsn = ijsendv_cyc(2,icycpts,ndn) - mj0(1)

            call mk_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1,basic_g(1)%vp,scratch%scr1,isn,jsn)
            call MPI_Pack(scratch%scr1,mmzp(1),MPI_REAL,node_buffs_cyc(ndn)%buffsend_cyc   &
                         ,nbuffsend_cyc,ipos,MPI_COMM_WORLD,ierr)
         end do
      end select

      mpiid = 1000000000+maxmach*10*(machs(mynum)-1)+10*(machs(mdn)-1)+ isflag
      call MPI_Isend(node_buffs_cyc(ndn)%buffsend_cyc,ipos    &
                    ,MPI_PACKED,mdn,mpiid ,MPI_COMM_WORLD, isend_req_cyc(isflag,mdn),ierr)
   end do

   return
end subroutine node_sendcyclic
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine node_getcyclic(isflag)
   use mem_grid
   use var_tables
   use cyclic_mod
   use mem_basic
   use mem_scratch
   use node_mod
   use grid_dims
   use mem_cuparm, only : nclouds
   use mem_aerad , only : nwave    ! ! intent(in)

   implicit none
   !----- Included variables --------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer                            , intent(in) :: isflag
   !----- Local variables -----------------------------------------------------------------!
   integer, dimension(MPI_STATUS_SIZE)             :: status
   integer                                         :: ndn
   integer                                         :: nsn
   integer                                         :: mdn
   integer                                         :: msn
   integer                                         :: numcols
   integer                                         :: icol
   integer                                         :: mtp
   integer                                         :: mtc
   integer                                         :: ind
   integer                                         :: nmp
   integer                                         :: mijrecv
   integer                                         :: ijr
   integer                                         :: istart
   integer                                         :: jstart
   integer                                         :: ivar
   integer                                         :: idn
   integer                                         :: jdn
   integer                                         :: iijj
   integer                                         :: nv
   integer                                         :: kg
   integer                                         :: nid
   integer                                         :: ierr
   integer                                         :: ipos
   integer                                         :: fdzp
   integer                                         :: fdep
   !---------------------------------------------------------------------------------------!


   if (ibnd /= 4 .and. jbnd /= 4) return

   !----- First, let's make sure our sends are all finished and de-allocated. -------------!
   do ndn = 1,ndns_cyc
      mdn = mdn_cyc(ndn)
      call MPI_Wait(isend_req_cyc(isflag,mdn),status,ierr)
   end do

   !----- Now, let's wait on our receives. ------------------------------------------------!


   !-----  Initialize start flags. --------------------------------------------------------!
   lstart_cyc(:,:,:) = 0

   mtp = nnzp(1)
   do nsn = 1,nsns_cyc
      msn = msn_cyc(nsn)
         
      call MPI_Wait(irecv_req_cyc(isflag,msn),status,ierr)

      !----- We got all our stuff.  Now unpack it into appropriate space. -----------------!
      ipos = 1
      select case (isflag)
      case (1,4,6)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos,mijrecv,1     &
                        ,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos,ijrecv_cyc    &
                        ,6*mijrecv,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         do ijr = 1,mijrecv

            select case (isflag)
            case (1)
               ivar = 4
               do nv = 1,num_var(1)
                  if ( vtab_r(nv,1)%impt1 == 1) then
                     ivar = ivar + 1
                     
                     !----- Find the variable dimensions in Z and E axes. -----------------!
                     call ze_dims(1,vtab_r(nv,1)%idim_type,.false.,fdzp,fdep)

                     !----- Copy the variable to a scratch array. -------------------------!
                     call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc        &
                                    ,ipos,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)

                     !----- Copy the variable from the scratch array to the pointer. ------!
                     call ex_cyc_buff(fdzp,mmxp(1),mmyp(1),fdep,vtab_r(nv,ngrid)%var_p     &
                                     ,scratch%scr2 ,ijr,mi0(1),mj0(1),ivar,mynum)
                  end if
               end do
            case (4)
               call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos         &
                              ,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
               call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                  &
                               ,basic_g(1)%pp,scratch%scr2 ,ijr,mi0(1),mj0(1),4,mynum)
            case (6)
               call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos         &
                              ,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
               call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                  &
                               ,basic_g(1)%wp,scratch%scr2 ,ijr,mi0(1),mj0(1),3,mynum)

               call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos,        &
                        scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
               call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                  &
                               ,basic_g(1)%pp,scratch%scr2 ,ijr,mi0(1),mj0(1),4,mynum)

            end select
         end do

      case (2)

         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,mijrecv,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,ijrecv_cyc,6*mijrecv,MPI_INTEGER,MPI_COMM_WORLD,ierr)

         do ijr = 1,mijrecv
            call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos            &
                           ,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                     &
                            ,basic_g(1)%up,scratch%scr2 ,ijr,mi0(1),mj0(1),1,mynum)
         end do

      case (3)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,mijrecv,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,ijrecv_cyc,6*mijrecv,MPI_INTEGER,MPI_COMM_WORLD,ierr)

         do ijr = 1,mijrecv
            call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos            &
                           ,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                     &
                            ,basic_g(1)%vp,scratch%scr2 ,ijr,mi0(1),mj0(1),2,mynum)
         end do
         
      case (5)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,mijrecv,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,ijrecv_cyc,6*mijrecv,MPI_INTEGER,MPI_COMM_WORLD,ierr)

         do ijr = 1,mijrecv
            call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos            &
                           ,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                     &
                            ,basic_g(1)%up,scratch%scr2 ,ijr,mi0(1),mj0(1),1,mynum)
         end do

         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,mijrecv,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos               &
                        ,ijrecv_cyc,6*mijrecv,MPI_INTEGER,MPI_COMM_WORLD,ierr)

         do ijr = 1,mijrecv
            call MPI_Unpack(node_buffs_cyc(nsn)%buffrecv_cyc,nbuffrecv_cyc,ipos            &
                           ,scratch%scr2 ,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_cyc_buff(mmzp(1),mmxp(1),mmyp(1),1                                     &
                            ,basic_g(1)%vp,scratch%scr2 ,ijr,mi0(1),mj0(1),2,mynum)
         end do
      end select

   end do

   return
end subroutine node_getcyclic
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will create the buffer for the cyclic boundary.                     !
!------------------------------------------------------------------------------------------!
subroutine mk_cyc_buff(nz,nx,ny,ne,mydata,buff,x,y)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: ne
   real, dimension(nz,nx,ny,ne), intent(in)    :: mydata
   real, dimension(*)          , intent(inout) :: buff
   integer                     , intent(in)    :: x
   integer                     , intent(in)    :: y
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: z
   integer                                     :: e
   integer                                     :: n
   !---------------------------------------------------------------------------------------!

   n = 0
   do e = 1,ne
      do z = 1,nz
         n       = n + 1
         buff(n) = mydata(z,x,y,e)
      end do
   end do

   return
end subroutine mk_cyc_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ex_cyc_buff(mz,mx,my,me,mydata,buff,ijr,i0,j0,ivar,mynum)
   use cyclic_mod
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                      , intent(in)    :: mz
   integer                      , intent(in)    :: mx
   integer                      , intent(in)    :: my
   integer                      , intent(in)    :: me
   integer                      , intent(in)    :: ijr
   integer                      , intent(in)    :: i0
   integer                      , intent(in)    :: j0
   integer                      , intent(in)    :: ivar
   integer                      , intent(in)    :: mynum
   real, dimension(mz,mx,my,me) , intent(inout) :: mydata
   real, dimension(*)           , intent(in)    :: buff
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: xdn
   integer                                      :: ydn
   integer                                      :: z
   integer                                      :: iijj
   integer                                      :: i
   integer                                      :: j
   integer                                      :: xsn
   integer                                      :: ysn
   integer                                      :: e
   integer                                      :: m
   real                                         :: sum
   real                                         :: self
   real                                         :: other
   !---------------------------------------------------------------------------------------!

   xsn   = ijrecv_cyc(1,ijr)
   ysn   = ijrecv_cyc(2,ijr)
   xdn   = ijrecv_cyc(3,ijr) - i0
   ydn   = ijrecv_cyc(4,ijr) - j0
   sum   = real(ijrecv_cyc(5,ijr) + ijrecv_cyc(6,ijr))
   self  = real(ijrecv_cyc(5,ijr)) / sum
   other = 1./ sum

   if (lstart_cyc(ivar,xdn,ydn) == 0) then
      lstart_cyc(ivar,xdn,ydn) = 1
      do e = 1,me
         do z = 1,mz
            mydata(z,xdn,ydn,e) = mydata(z,xdn,ydn,e) * self
         end do
      end do
   end if
   m = 0
   do e = 1, me
      do z = 1,mz
         m = m + 1
         mydata(z,xdn,ydn,e) = mydata(z,xdn,ydn,e) + buff(m) * other
      end do
   end do

   return
end subroutine ex_cyc_buff
!==========================================================================================!
!==========================================================================================!

