module dtset

  use io_params, only : &
       maxgrds           ! INTENT(IN)

  implicit none

  real    :: ssodx(maxgrds)
  integer :: idelx(maxgrds)

contains

  subroutine dtset_new(mynum, nndtflg, dxtmax_local)

    use mem_grid, only : &
         ideltat,        & ! INTENT(IN)
         ngrids,         & ! INTENT(IN)
         deltaxn,        & ! INTENT(IN)
         nnxp,           & ! INTENT(IN)
         nnyp,           & ! INTENT(IN)
         nnzp,           & ! INTENT(IN)
         grid_g,         & ! INTENT(IN)
         sspct,          & ! INTENT(INOUT) ! Computed locally
         nxtnest,        & ! INTENT(IN)
         dtlongn,        & ! INTENT(INOUT) ! Computed locally 
         dtlong,         & ! INTENT(IN)    ! Passed to slaves
         nndtrat,        & ! INTENT(INOUT) ! Initial value passed to slaves
         ! Computed when this subroutine is updated.
         ! Retirar da com.:Mest->escravo
         nnacoust,       & ! INTENT(INOUT) ! Computed locally
         nacoust,        & ! INTENT(IN)
         cflxy,          & ! INTENT(IN)    ! Computed at cfl
         cflz,           & ! INTENT(IN)    ! To be computed at modsched local
         iflag,          & ! INTENT(INOUT) ! Defined locally
         timmax,         & ! INTENT(IN)
         time              ! INTENT(IN)
    
    ! "cflxy" and "cflz" are locally computed, and must be passed to slaves
    ! to determine the largest value.
    
    use rconstants, only : &
         cpdry,            & ! INTENT(IN)
         cvdry,            & ! INTENT(IN)
         rdry                ! INTENT(IN)
    
    use ref_sounding, only : &
         th01dn,             & ! INTENT(IN)
         pi01dn                ! INTENT(IN)
    ! Receive after the update/call of varf_read at RAMS_OUTPUT
    
    use io_params, only : &
         maxgrds,         & ! INTENT(IN)
         nzpmax,          & ! INTENT(IN)
         frqanl             ! INTENT(IN)
    use therm_lib, only : &
         extheta2temp       ! function

    implicit none
    
    ! Arguments:
    real(kind=8),external :: dmin2
    integer, intent(in)  :: mynum
    integer, intent(out) :: nndtflg
    real, intent(in)     :: dxtmax_local(maxgrds)
    
    ! Local variables:
    integer, parameter  :: ndx=37,ndt=42
    integer, dimension(maxgrds) :: idelt,nndtrat1,nnacoust1
    real, dimension(maxgrds) :: sscourn,dtlongn1
    real, dimension(nzpmax) :: vctr1
    real :: cflnumh, cflnumv, delx(ndx), delt(ndt)

    integer :: iabsdt,ifm,id,n2,n3,k,nn2,nn3,icm,ntf,ii
    real :: ssmax,tmax,dxtmax,sscnmax,sspct0,cflxyz,timeleft
    real(kind=8) :: dtt,dft
    
    real :: dxta, dxtb, dxtc, dxtd

    delx(:) = (/ &
         200000.,150000.,100000.,80000.,70000.,60000.,40000.,  &
         030000., 20000., 10000., 6000., 4000., 3000., 2000.,  &
         001000.,   800.,   600.,  500.,  400.,  300.,  200.,  &
         000150.,   100.,    80.,   60.,   50.,   40.,   30.,  &
         000020.,    10.,     8.,    6.,    5.,    4.,    3.,  &
         000002.,     1. /)

    delt(:) = (/ &
         300.,  240.,   180.,  150.,  120.,   90.,   60.,  &
         050.,   40.,    30.,   20.,   15.,   12.,   10.,  &
         006.,    5.,     4.,    3.,   2.5,   2.0,   1.5,  &
         01.2,    1.0,     .8,    .6,   .5,    .4,    .3,  &
         00.2,     .1,     .08,   .06,  .05,   .04,   .03,  &
         00.02,    .01,    .008,  .006, .005,  .004,  .003/)

    idelx(1) = 0
    
    iabsdt = abs(ideltat)
    
    ! On the first call to this subroutine, initialize idelx, ssodx, dtlongn,
    ! nnacoust, and if required, nndtrat.
    
    if ( idelx(1)==0 ) then
       
       cflnumh = .90
       cflnumv = .90
       
       do ifm = 1,ngrids
          
          do id = ndx,1,-1
             if (delx(id)<=deltaxn(ifm))  idelx(ifm) = id
          enddo
          
          n2 = nnxp(ifm)
          n3 = nnyp(ifm)
          do k = 1,nnzp(ifm)
             vctr1(k) = extheta2temp(pi01dn(k,1),th01dn(k,1))
          enddo
          tmax = maxval(vctr1(1:nnzp(ifm)))
          ssmax = sqrt(cpdry / cvdry * rdry * tmax)
          
          nn2 = nnxp(ifm)
          nn3 = nnyp(ifm)

          if (mynum==0) then
             ! Master process
             dxtmax = max(grid_g(ifm)%dxt(1,1), &
                  grid_g(ifm)%dxt(nn2,1),       &
                  grid_g(ifm)%dxt(nn2,nn3),     &
                  grid_g(ifm)%dxt(1,nn3)        )
          else
             ! Slave (node) process
             dxtmax = dxtmax_local(ifm)
          endif
          
          ssodx(ifm) = ssmax * dxtmax
          
       enddo
       
       if ( ideltat==0 ) then
          
          sspct = 1.
          do ifm = 1,ngrids
             icm = nxtnest(ifm)
             if ( icm==0 ) then
                dtlongn(ifm) = dtlong
             else
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
             nnacoust(ifm) = nacoust
             sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
             sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          enddo
          
       else
          
          sscnmax = 0.
          do ifm = 1,ngrids
             icm = nxtnest(ifm)
             dtlongn(ifm) = delt(idelx(ifm)+iabsdt-1)
             
             ! For coarse grid(s) adjust dtlongn so that it is an integer
             ! divisor of FRQANL.  For nested grids, compute nndtrat(ifm) as the
             ! first integer greater than or equal to the timestep ratio between
             ! a nested grid's parent and the nested grid. Then adjust 
             ! dtlongn(ifm) for the nested grid to be the parent grid timestep 
             ! divided by nndtrat(ifm).
             
             if ( icm==0 ) then
                ntf = nint(frqanl / dtlongn(1))
                dtlongn(ifm) = frqanl / ntf
             else
                nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
             
             ! Compute sst courant numbers (sstcourn(ifm)) for long timestep
             ! dtlongn.
             
             sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
             if ( sscourn(ifm)>sscnmax)  sscnmax = sscourn(ifm)
             
          enddo
          
          ! Define trial sspct0 in terms of sscnmax using nonlinear formula 
          ! intended to increase nnacoust as sspct decreases, but limit sspct0
          ! to a minimum of .2.
          
          sspct0 = min(1., (.95/sscnmax)**.5)
          
          if ( sspct0<.2 ) then
             call abort_run ('Sound speed percent is forced to be too low'                 &
                            ,'dtset_new','local_proc.f90')
          end if
          
          sspct = 1.
          do ifm = 1,ngrids
             nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
             sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          end do
          
       endif
       
    endif
    
    
    ! check Courant numbers
    
    nndtflg = 0
    
    if ( ideltat>=0 ) then
       
       do ifm = 1,ngrids
          cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
          if ( cflxyz>1. ) then
             iflag = 1
             print*, 'Model will stop because CFL limit exceeded (LOCAL).'
             print*, 'Considering (ideltat>=0): time, cflxyz =', &
                  time, cflxyz
             call flush(6)
             call abort_run('CFL limit exceeded...','dtset_new','local_proc.f90')
          endif
       enddo
       
    else
       
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          cflxyz = max(cflxy(ifm)/cflnumh, cflz(ifm)/cflnumv)
          
          nndtrat1(ifm) = nndtrat(ifm)
          dtlongn1(ifm) = dtlongn(ifm)
          nnacoust1(ifm) = nnacoust(ifm)
          
          do id = ndt,1,-1
             if ( delt(id)*cflxyz<=dtlongn(ifm) ) idelt(ifm) = id
          enddo
          
          if ( idelt(ifm)>idelx(ifm)+4 ) then
             print*, 'Adjustable timestep forced to become too small'
             print*, 'on grid ',ifm,' idelx= ', idelx(ifm),' idelt =',idelt(ifm)
             print*, 'Model will stop'
             iflag = 1
          else
             ii = max(idelx(ifm)+iabsdt-1,idelt(ifm))
             dtlongn(ifm) = delt(ii)
             
             ! For the coarse grid(s), adjust dtlongn(1) to not jump past an
             ! analysis write time or the end of the model run time.
             
             if ( icm==0 ) then

                !Use double precision intrinsics, bc timmax and time are RGK 5-29-07

                dtt = timmax-time
                dft = dble(frqanl) - dmod(time,dble(frqanl))
                timeleft = real( dmin2(dtt,dft))


                if ( dtlongn(1)>.95 * timeleft) then
                   dtlongn(ifm) = timeleft
                endif
             else
                nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
          endif
          
          ! Compute sst courant numbers (sstcourn(ifm))for long timestep dtlongn.
          
          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)
       enddo
       
       ! Define trial sspct0 in terms of sscnmax using nonlinear formula
       ! intended to increase nnacoust as sspct decreases, but limit sspct0 to 
       ! a minimum of .2.
       
       sspct0 = min(1., (.95/sscnmax) ** .5)
       if ( sspct0<.2 ) then
          print*, 'Sound speed percent is forced to be too low'
          stop 'low_sspct0'
       endif
       
       sspct = 1.
       do ifm = 1,ngrids
          nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          
          ! If there are any updates to dtlongn(ifm), print out new values.
          
          if (abs(dtlongn(ifm) - dtlongn1(ifm)) > 1.e-3) then
             write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),dtlongn(ifm)
122          format('Timestep update: ngrid, nndtrat, nnacoust,'  &
                  ,' dtlongn = ',3i3,f10.3)
          endif
          
          ! If there are any updates to nndtrat(ifm) or others, set 
          ! nndtflg = 1 to flag new call to subroutine modsched and send new 
          ! stuff to nodes.
          
          if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
               nnacoust(ifm) /= nnacoust1(ifm) .or. &
               dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1
          
       enddo
       
    endif
    
    return
  end subroutine dtset_new

end module dtset

!--------------------------------------------------------------------

subroutine master_putdxt(master_num)

  use mem_grid, only : &
       nnxp,           & ! INTENT(IN)
       nnyp,           & ! INTENT(IN)
       ngrids,         & ! INTENT(IN)
       grid_g            ! INTENT(IN)
  use io_params, only : &
       maxgrds           ! INTENT(IN)
  use rpara, only : &
       nmachs,      &    ! INTENT(IN)
       machnum           ! INTENT(IN)

  implicit none

  ! Local Variables:
  integer :: ifm, nn2, nn3
  real    :: dxtmax(maxgrds)
  integer :: master_num,ierr
  include 'mpif.h'

  ! Calculating DXTMAX
  do ifm = 1,ngrids
     nn2 = nnxp(ifm)
     nn3 = nnyp(ifm)
     dxtmax(ifm) = max(grid_g(ifm)%dxt(1,1), grid_g(ifm)%dxt(nn2,1), &
          grid_g(ifm)%dxt(nn2,nn3),grid_g(ifm)%dxt(1,nn3))
  enddo
  do ifm = ngrids+1,maxgrds
    dxtmax(ifm) = 0.
  end do

  call MPI_Bcast(dxtmax,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  return
end subroutine master_putdxt

!--------------------------------------------------------------------

subroutine master_putcflmax(master_num)

  use io_params, only : &
       maxgrds           ! INTENT(IN)
  use mem_grid, only : &
       cflxy,          & ! INTENT(IN)
       cflz              ! INTENT(IN)
  use rpara, only : &
       nmachs,      &    ! INTENT(IN)
       machnum,     &    ! INTENT(IN)
       mainnum           ! intent(in)
       
  include 'mpif.h'
  integer :: master_num,ierr

  call MPI_Bcast(cflxy,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cflz,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  return
end subroutine master_putcflmax

!--------------------------------------------------------------------

subroutine node_getdxt(dxtmax_local)

  use io_params, only : &
       maxgrds           ! INTENT(IN)
  use node_mod, only: master_num
  implicit none

  include 'mpif.h'
  real, intent(out) :: dxtmax_local(maxgrds)
  integer :: ierr
  
  call MPI_Bcast(dxtmax_local,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  return

end subroutine node_getdxt

!--------------------------------------------------------------------

subroutine node_getcflmax()

  use io_params, only : &
       maxgrds           ! INTENT(IN)
  use mem_grid, only : &
       cflxy,          & ! INTENT(OUT)
       cflz              ! INTENT(OUT)
  use node_mod, only : master_num

  implicit none

  include 'mpif.h'
  integer :: ierr


  call MPI_Bcast(cflxy,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cflz,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

  return
end subroutine node_getcflmax
