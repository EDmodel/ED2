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

subroutine masterput_oda(master_num)

  use grid_dims
  use mem_oda
  use rpara

  implicit none

  !   +------------------------------------------------------------------
  !   ! This routine gives obs data to the nodes for ODA. All obs are sent 
  !   !    to all nodes.
  !   +------------------------------------------------------------------
  include 'interface.h'
  include 'mpif.h'
  integer :: ns
  integer :: master_num
  integer :: ierr

  ! Namelist info
  call MPI_Bcast(wt_oda_grid,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(frqoda,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(todabeg,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(todaend,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tnudoda,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_uv,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_th,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_pi,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_rt,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_sfc_til,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_sfc_tel,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_upa_til,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_upa_tel,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_sfce,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_sfc0,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_upae,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_upa0,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_zfact,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_hgt,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

  ! Sizes for arrays
  call MPI_Bcast  (num_oda_sfc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (nsfcfiles,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (maxtimes_sfc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (num_oda_upa,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (nupafiles,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (maxtimes_upa,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)


  ! Surface obs

  do ns=1,num_oda_sfc
     call MPI_Bcast (oda_sfc_info(ns)%id,8,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_sfc_info(ns)%intid,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_sfc_info(ns)%ntimes,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_sfc_info(ns)%iactive(1),maxodagrids,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xista(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xjsta(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xlat,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xlon,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xsta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%ysta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%stopo,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

     !Changed to double RGK 5-29-07
     call MPI_Bcast(oda_sfc_obs(ns)%time(1) ,maxtimes_sfc,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%temp(1) ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%dewpt(1),maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%us(1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%vs(1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%u (1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%v (1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%ps(1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  enddo

  ! Upper air obs


  do ns=1,num_oda_upa
     call MPI_Bcast (oda_upa_info(ns)%id,8,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_upa_info(ns)%intid,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_upa_info(ns)%ntimes,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_upa_info(ns)%iactive(1),maxodagrids,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xista(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xjsta(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xlat,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xlon,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xsta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%ysta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%stopo,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

     !Change to double RGK 5-27-07
     call MPI_Bcast(oda_upa_obs(ns)%time(1),maxtimes_upa,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%lp  (1),maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%lz  (1),maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%theta(1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%rv   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%us   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%vs   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%zz   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%u    (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%v    (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%pi   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%zgeo (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  enddo

  return
end subroutine masterput_oda


!     ****************************************************************

subroutine nodeget_oda()

  use node_mod
  use mem_oda

  implicit none

  include 'interface.h'
  include 'mpif.h'
  integer :: ierr
  integer :: ns

  ! Namelist info
  call MPI_Bcast(wt_oda_grid,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(frqoda,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(todabeg,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(todaend,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tnudoda,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_uv,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_th,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_pi,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(wt_oda_rt,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_sfc_til,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_sfc_tel,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_upa_til,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(oda_upa_tel,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_sfce,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_sfc0,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_upae,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_upa0,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_zfact,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(roda_hgt,maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

  ! Sizes for arrays
  call MPI_Bcast  (num_oda_sfc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (nsfcfiles,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (maxtimes_sfc,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (num_oda_upa,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (nupafiles,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast  (maxtimes_upa,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

  ! Allocate main obs arrays

  call oda_obs_alloc()

  ! Surface obs

  do ns=1,num_oda_sfc
     call MPI_Bcast (oda_sfc_info(ns)%id,8,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_sfc_info(ns)%intid,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_sfc_info(ns)%ntimes,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_sfc_info(ns)%iactive(1),maxodagrids,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xista(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xjsta(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xlat,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xlon,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%xsta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%ysta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_info(ns)%stopo,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     
     call MPI_Bcast(oda_sfc_obs(ns)%time(1) ,maxtimes_sfc,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%temp(1) ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%dewpt(1),maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%us(1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%vs(1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%u (1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%v (1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_sfc_obs(ns)%ps(1)   ,maxtimes_sfc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  enddo


  ! Upper air obs

  do ns=1,num_oda_upa
     call MPI_Bcast (oda_upa_info(ns)%id,8,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_upa_info(ns)%intid,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_upa_info(ns)%ntimes,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast  (oda_upa_info(ns)%iactive(1),maxodagrids,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xista(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xjsta(1),maxodagrids,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xlat,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xlon,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%xsta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%ysta,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_info(ns)%stopo,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

     call MPI_Bcast(oda_upa_obs(ns)%time(1),maxtimes_upa,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%lp  (1),maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%lz  (1),maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%theta(1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%rv   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%us   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%vs   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%zz   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%u    (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%v    (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%pi   (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(oda_upa_obs(ns)%zgeo (1,1),maxupalevs*maxtimes_upa,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  enddo

  return
end subroutine nodeget_oda
