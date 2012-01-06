!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine opspec1

  ! this routine checks the option specifications in the $model_grids
  !   namelist for consistency, and overrides settings of icloud,
  !   irain, ipris, isnow, iaggr, igraup, and ihail, setting them
  !   all to zero if level is less than 3.

  use mem_grid
  use therm_lib, only: level
  use micphys

  implicit none

  integer :: ierr
  integer :: ifm
  integer :: icm
  integer :: ng
  integer :: ifaterr
  integer :: lev4bins
  integer :: lev5bins
  integer :: nhemgrds
  logical :: twod
  character(len=10) :: c0, c1, c2
  character(len=*), parameter :: h="**(opspec1)**"
 
  ifaterr=0

  ! check if number of grids is within bounds;
  ! if not, stop analyzing input data since it may be corrupted

  if (ngrids < 1 .or. ngrids > maxgrds) then
     write(unit=*,fmt='(a,1x,i6)') 'Ngrids=',ngrids
     write(unit=*,fmt='(a,1x,i6)') 'Maxgrds=',maxgrds
     call opspec_mess('Number of grids cannot be larger than the maximum allowed.' &
                   ,'opspec1')
     ifaterr = ifaterr + 1
  endif

  ! flag 2d simulation

  twod = all(nnyp(1:ngrids) == 1)

  ! check if number of grid points is within bounds

  do ngrid=1,ngrids
     if (nnxp(ngrid) < 4 .or. nnxp(ngrid) > nxpmax) then
        write(unit=*,fmt='(a,1x,i6)') 'ngrid  =',ngrid
        write(unit=*,fmt='(a,1x,i6)') 'nnxp   =',nnxp(ngrid)
        write(unit=*,fmt='(a,1x,i6)') 'nxpmin =',4
        write(unit=*,fmt='(a,1x,i6)') 'nxpmax =',nxpmax
        call opspec_mess('Number of x points is outside allowed interval.','opspec1')
        ifaterr = ifaterr + 1
     end if
     if (.not. twod .and. (nnyp(ngrid) < 4 .or. nnyp(ngrid) > nypmax)) then
        write(unit=*,fmt='(a,1x,i6)') 'ngrid  =',ngrid
        write(unit=*,fmt='(a,1x,i6)') 'nnyp   =',nnyp(ngrid)
        write(unit=*,fmt='(a,1x,i6)') 'nypmin =',4
        write(unit=*,fmt='(a,1x,i6)') 'nypmax =',nypmax
        call opspec_mess('Number of y points is outside allowed interval.','opspec1')
        ifaterr = ifaterr + 1
     end if
     if (nnzp(ngrid) < 11 .or. nnzp(ngrid) > nzpmax) then
        write(unit=*,fmt='(a,1x,i6)') 'ngrid  =',ngrid
        write(unit=*,fmt='(a,1x,i6)') 'nnzp   =',nnzp(ngrid)
        write(unit=*,fmt='(a,1x,i6)') 'nzpmin =',11
        write(unit=*,fmt='(a,1x,i6)') 'nzpmax =',nzpmax
        call opspec_mess('Number of z points is outside allowed interval.','opspec1')
        ifaterr = ifaterr + 1
     end if
  enddo

  ! check if number of soil levels is below maximum
  ! **(jp)**:  what is the minumum number of soil levels? i imposed 0.

  if (nzg < 1 .or. nzg > nzgmax) then
        write(unit=*,fmt='(a,1x,i6)') 'nzg    =',nzg        
        write(unit=*,fmt='(a,1x,i6)') 'nzgmin =',1
        write(unit=*,fmt='(a,1x,i6)') 'nzgmax =',nzgmax
        call opspec_mess('Number of soil layers is outside allowed interval.','opspec1')
     ifaterr = ifaterr + 1
  end if

  ! atmospheric levels must exceed soil levels

  do ngrid=1,ngrids
     if (nnzp(ngrid) <= nzg+nzs) then
        write(unit=*,fmt='(a,1x,i6)') 'ngrid  =',ngrid      
        write(unit=*,fmt='(a,1x,i6)') 'nnzp   =',nnzp(ngrid)
        write(unit=*,fmt='(a,1x,i6)') 'nzg    =',nzg   
        write(unit=*,fmt='(a,1x,i6)') 'nzs    =',nzs   
        call opspec_mess('You must have more atmospheric vertical levels than nzg+nzs.'&
                        ,'opspec1')
        ifaterr = ifaterr + 1
     endif
  enddo

  ! consistency of nstratx, nstraty

  do ifm=1,ngrids
     if (nstratx(ifm) < 1) then
        write(unit=*,fmt='(a,1x,i6)') 'ngrid  =',ifm      
        write(unit=*,fmt='(a,1x,i6)') 'nstratx=',nstratx(ifm)
        call opspec_mess('Nesting ratio in x direction (nstratx) must be at least 1.'&
                        ,'opspec1')
        ifaterr = ifaterr + 1
     endif
     if (nstraty(ifm) < 1) then
        write(unit=*,fmt='(a,1x,i6)') 'ngrid  =',ifm      
        write(unit=*,fmt='(a,1x,i6)') 'nstraty=',nstraty(ifm)
        call opspec_mess('Nesting ratio in y direction (nstraty) must be at least 1.'&
                        ,'opspec1')
        ifaterr = ifaterr + 1
     endif
  end do

  if (twod .and. any(nstraty(1:ngrids) /= 1)) then
     write(*,"(a)") h//" for 2d simulations, all nstraty should be 1; at least one is not"
     ifaterr = ifaterr + 1
  end if

  do ifm=1,ngrids
     icm = nxtnest(ifm)
     if (icm >= 1 .and. nnyp(ifm) == 1 .and.  &
          (ninest(ifm) < 3 .or. njnest(ifm) < 1)) then
        print*, ' fatal - nested 2d grid must have ninest > 2 '  &
             ,'and njnest = 1 in namelist.'
        ifaterr=ifaterr+1
     endif
  enddo

  ! allowable values of centlat, centlon, polelat, polelon
  !   (severity - f)

  if(polelat < -90..or.polelat > 90.) then
     print*,' fatal - polelat outside of legal bounds.'
     ifaterr=ifaterr+1
  endif

  if(polelon < -180..or.polelon > 180.) then
     print*,' fatal - polelon outside of legal bounds.'
     ifaterr=ifaterr+1
  endif

  do ng=1,ngrids
     if(centlat(ng) < -90..or.centlat(ng) > 90.) then
        print*,' fatal - centlat outside of legal bounds.'
        ifaterr=ifaterr+1
     endif

     if(centlon(ng) < -180..or.centlon(ng) > 180.) then
        print*,' fatal - centlon outside of legal bounds.'
        ifaterr=ifaterr+1
     endif
  enddo

  ! check nxtnest values for validity and whether this is a global simulation
  !   (severity - f)

  if (nxtnest(1) /= 0) then
     print*, ' fatal - grid # 1 must have its parent mesh'  &
          ,' designated as grid # 0 (nxtnest(1) = 0)'
     ifaterr=ifaterr+1
  endif

  nhemgrds = 0
  do ifm = 1,ngrids
     icm = nxtnest(ifm)

     if (icm >= ifm) then
        print 1, ifm
1       format (' fatal - nest #',i3,' has specified parent'  &
             ,' mesh of equal or higher number')
        ifaterr=ifaterr+1
     endif

     if (icm < 0) then
        print 2, ifm
2       format (' fatal - nest #',i3,' has specified parent'  &
             ,' mesh of number less than 0')
        ifaterr=ifaterr+1
     endif

     if (icm == 0) then
        nhemgrds = nhemgrds + 1
        nhemgrd2 = ifm
     endif

  enddo

  if (nhemgrds > 2) then
     print*, ' fatal - more than two grids have grid # 0 specified'  &
          ,' as their parent'
     ifaterr=ifaterr+1
  endif

  ! if this is a global simulation, print message that deltax and deltaz will be
  ! redefined by nnxp(1) and nnyp(1).  check that nnxp, nnyp, and nnzp of the
  ! two top grids are identical, and that nnxp and nnyp are equal to each
  ! other.  check to make sure that cyclic lateral boundary conditions are
  ! not specified.

  !print*, 'nhemgrd2,nhemgrds',nhemgrd2,nhemgrds
  !print*, 'nnxp(1),nnxp(nhemgrd2)',nnxp(1),nnxp(nhemgrd2)
  !print*, 'nnyp(1),nnyp(nhemgrd2)',nnyp(1),nnyp(nhemgrd2)
  !print*, 'nnzp(1),nnzp(nhemgrd2)',nnzp(1),nnzp(nhemgrd2)

  if (nhemgrds == 2) then
     if (nnxp(1) /= nnxp(nhemgrd2) .or.  &
         nnyp(1) /= nnyp(nhemgrd2) .or.  &
         nnzp(1) /= nnzp(nhemgrd2)) then
        print*, ' fatal - for a global simulation, nnxp, nnyp, and nnzp'
        print*, ' must be identical between both hemispheric grids'
        ifaterr=ifaterr+1
     endif

     if (nnxp(1) /= nnyp(1)) then
        print*, ' '
        print*, 'fatal - for a global simulation, nnxp must equal'
        print*, ' nnyp in both hemispheric grids'
        ifaterr = ifaterr + 1
     endif

     if (ibnd == 4 .or. jbnd == 4) then
        print*, ' '
        print*, 'fatal - for a global simulation, ibnd and jbnd'
        print*, ' must not be set to 4'
        ifaterr = ifaterr + 1
     endif

     print*, ' '
     print*, 'because two values of nxtnest are set to 0, this'
     print*, ' is configured as a global simulation.'
     print*, 'consequently, ihtran will automatically be set'
     print*, ' to 1, and deltax and deltay will be redefined'
     print*, ' in terms of nnxp(1) and nnyp(1), ignoring'
     print*, ' the values specified in the namelist.'

  endif

  ! check to make sure that top grids have nsttop and nstbot set to 1
  !   (severity - f)

  if (nnsttop(1) /= 1 .or. nnsttop(nhemgrd2) /= 1) then
     print*,'nsttop not set to 1 for a hemispheric grid'
     ifaterr=ifaterr+1
  endif

  if (nnstbot(1) /= 1 .or. nnstbot(nhemgrd2) /= 1) then
     print*,'nstbot not set to 1 for a hemispheric grid'
     ifaterr=ifaterr+1
  endif

  ! if level is less than 3, set microphysics parameters to zero.
  ! if level equals 3, check for values of microphysics parameters
  ! that are out of bounds.  if level is equal to 4, set microphysics
  ! parameters other than icloud to zero.

  if (level == 2) then

     icloud = 0
     irain = 0
     ipris = 0
     isnow = 0
     iaggr = 0
     igraup = 0
     ihail = 0

  elseif (level == 3) then

     if (icloud < 0 .or. icloud > 7) then
        print*,'fatal - icloud out of range'
        ifaterr = ifaterr + 1
     endif
     if (irain < 0 .or. irain > 5) then
        print*,'fatal - irain out of range'
        ifaterr = ifaterr + 1
     endif
     if (ipris <  0 .or. ipris > 7) then
        print*,'fatal - ipris out of range'
        ifaterr = ifaterr + 1
     endif
     if (isnow < 0 .or. isnow > 5) then
        print*,'fatal - isnow out of range'
        ifaterr = ifaterr + 1
     endif
     if (iaggr < 0 .or. iaggr > 5) then
        print*,'fatal - iaggr out of range'
        ifaterr = ifaterr + 1
     endif
     if (igraup < 0 .or. igraup > 5) then
        print*,'fatal - igraup out of range'
        ifaterr = ifaterr + 1
     endif
     if (ihail < 0 .or. ihail > 5) then
        print*,'fatal - ihail out of range'
        ifaterr = ifaterr + 1
     endif
     !----- Moved from mic_init to here ---------------------------------------------------!
     if(mkcoltab < 0.or.mkcoltab > 1)then
        print*, 'mkcoltab set to ',mkcoltab, 'which is out of bounds'
        ifaterr = ifaterr + 1
     endif
     !-------------------------------------------------------------------------------------!
     

  elseif (level == 4) then

     ipris = 0
     isnow = 0
     iaggr = 0
     igraup = 0
     ihail = 0

  endif

  ! if level is 4, make sure that naddsc is large enough for number
  !   of bins specified in icloud.

  if (level == 4) then
     if (irain == 0) then
        lev4bins = 2 * icloud + 1
        if (naddsc .lt. lev4bins) then
           print*, 'fatal - naddsc is not large enough for icloud'
           print*, 'value with level = 4.'
           print*, 'naddsc must be at least ',lev4bins
           ifaterr = ifaterr + 1
        endif
     else
        lev4bins = 2 * icloud + 1 + 67
        if (naddsc .lt. lev4bins) then
           print*, 'fatal - naddsc is not large enough for icloud'
           print*, 'value with level = 4 and irain = 1.'
           print*, 'naddsc must be at least ',lev4bins
           ifaterr = ifaterr + 1
        endif
     endif
  endif

  ! if level is 5, make sure that naddsc is large enough for number
  !   of bins specified in icloud.

  if (level == 5) then
     if (irain == 0) then
        lev5bins = 2 * (icloud + ipris + iaggr + igraup) + 5
        if (naddsc .lt. lev5bins) then
           print*, 'fatal - naddsc is not large enough for icloud,'
           print*, 'ipris, iaggr, and igraup values with level = 5.'
           print*, 'naddsc must be at least ',lev5bins
           ifaterr = ifaterr + 1
        endif
     else
        lev5bins = 2 * (icloud + ipris + iaggr + igraup) + 5 + 67
        if (naddsc .lt. lev5bins) then
           print*, 'fatal - naddsc is not large enough for icloud,'
           print*, 'ipris, iaggr, and igraup values with level = 5'
           print*, 'and irain = 1.'
           print*, 'naddsc must be at least ',lev5bins
           ifaterr = ifaterr + 1
        endif
     endif
  endif

  ! stop the run if there are any fatal errors.  list how many
  !   warning and informative errors.

  if (ifaterr > 0) then
     write(unit=*,fmt='(a)'      ) '----------------------- OPSPEC1 -----------------------'
     write(unit=*,fmt='(a,1x,i5)') ' fatal     errors - ',ifaterr
     write(unit=*,fmt='(a,1x,i5)') ' ------------------------------------------------------'
     call abort_run('Fatal errors at namelist.','opspec1','opspec.f90')
  end if
end subroutine opspec1

! ************************************************************************

subroutine opspec2

  ! check that fine mesh is a valid subset of its coarser mesh.
  !   and that top and bottom boundary flags are set correctly.
  !   (severity - f)

  use mem_grid
  use mem_varinit
  use mem_radiate, only: ISWRTYP, ILWRTYP ! Intent(in)
  use therm_lib, only: level

  implicit none

  integer :: icm,ifm,ifaterr,ncx,ncy,nfx,nfxp,nfy,nfyp
  integer :: ng,nesta,nfz,kc
  character(len=*), parameter :: h="**(opspec2)**"
  logical :: ex

  ifaterr=0

  do ifm=1,ngrids

     icm = nxtnest(ifm)
     if (icm .ge. 1) then

        ncx=(nnxp(ifm)-2)/nstratx(ifm)
        ncy=(nnyp(ifm)-2)/nstraty(ifm)

        if ((nnyp(ifm).eq.1.and.nnyp(icm).ne.1).or.  &
             (nnyp(ifm).ne.1.and.nnyp(icm).eq.1)) then
           print*,' fatal - grids must be either all 3-d or all 2-d'
           ifaterr=ifaterr+1
        endif

        if (ninest(ifm).lt.3) then
           print 11, ifm
11         format (' fatal - nest #',i3,' too close to western'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (njnest(ifm).lt.3.and.nnyp(ifm).gt.1) then
           print 12, ifm
12         format (' fatal - nest #',i3,' too close to southern'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (nknest(ifm).lt.3.and.nnstbot(ifm).eq.0) then
           print 13, ifm
13         format (' fatal - nest #',i3,' too close to lower'  &
                ,' boundary of coarser mesh or nnstbot incorrect')
           ifaterr=ifaterr+1
        endif

        if (nknest(ifm).ne.1.and.nnstbot(ifm).eq.1) then
           print 14, ifm
14         format (' fatal - nest #',i3,' not to lower boundary of'  &
                ,' coarser mesh or nnstbot flag set incorrectly')
           ifaterr=ifaterr+1
        endif

        if (ninest(ifm)+ncx.gt.nnxp(icm)-3) then
           print 15, ifm
15         format (' fatal - nest #',i3,' too close to eastern'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (njnest(ifm)+ncy.gt.nnyp(icm)-3.and.nnyp(ifm).gt.1) then
           print 16, ifm
16         format (' fatal - nest #',i3,' too close to northern'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (ncx.lt.2.or.(ncy.lt.2.and.nnyp(ifm).ne.1)) then
           print 17, ifm
17         format (' fatal - nest #',i3,' dimensioned too small'  &
                ,' in at least one horizontal direction')
           ifaterr=ifaterr+1
        endif

        nfx=ncx*nstratx(ifm)+2
        if (nnxp(ifm).ne.nfx) then
           nfxp=nfx+nstratx(ifm)
           print 18, ifm,nfxp,nfx
18         format (' fatal - nest #',i3,' nnxp incompatible with'  &
                ,' nstratx:  may increase nnxp to',i5  &
                ,' or decrease it to',i5)
           ifaterr=ifaterr+1
        endif

        nfy=ncy*nstraty(ifm)+2
        if (nnyp(ifm).ne.nfy.and.nnyp(ifm).gt.1) then
           nfyp=nfy+nstraty(ifm)
           print 19, ifm,nfyp,nfy
19         format (' fatal - nest #',i3,' nnyp incompatible with'  &
                ,' nstraty:  may increase nnyp to',i5  &
                ,' or decrease it to',i5)
           ifaterr=ifaterr+1
        endif

        if (nnstbot(ifm).eq.1.and.nnstbot(icm).ne.1) then
           print 20, ifm
20         format (' fatal - nest #',i3,' nnstbot flag'  &
                ,' incompatible with coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (nnsttop(ifm).eq.1.and.nnsttop(icm).ne.1) then
           print 21, ifm
21         format (' fatal - nest #',i3,' nnsttop flag incompatible'  &
                ,' with coarser mesh')
           ifaterr=ifaterr+1
        endif

     endif
  enddo

  nesta=abs(nestz1)
  if(nestz1.ne.0.and.nesta.le.ngrids)then
     nfz=nnzp(nesta)-2
     kc=nknest(nesta)
1002 continue
     kc=kc+1
     nfz=nfz-nrz(kc,nesta)
     if(nfz.lt.0)then
        print 195,nesta
195     format(' fatal - vertically nested grid #',i3,  &
             ' has illegal number of levels for given nstratz values')
        ifaterr=ifaterr+1
     endif
     if(nfz.gt.0)go to 1002
     if(nfz.eq.0)then
        if(kc.gt.nnzp(nxtnest(nesta))-3.and.nnsttop(nesta).ne.1)then
           print 22, nesta
22         format (' fatal - nest #',i3,' too high'  &
                ,' or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif

        if(kc.ne.nnzp(nxtnest(nesta))-1.and.nnsttop(nesta).eq.1)then
           print 23, nesta
23         format (' fatal - nest #',i3,' not to upper boundary of'  &
                ,' coarser mesh or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif
     endif
  endif

  do ifm = 1,ngrids
     icm = nxtnest(ifm)
     if (ifm .ne. nesta .and. icm .ge. 1)then
        if(nnzp(ifm).gt.nnzp(icm)-nknest(ifm)-1.and.nnsttop(ifm).eq.0)then
           print 24, ifm
24         format (' fatal - nest #',i3,' nnzp incompatible with'  &
                ,' parent grid or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif

        if(nnzp(ifm).ne.nnzp(icm)-nknest(ifm)+1.and.nnsttop(ifm).eq.1)then
           print 25, ifm
25         format (' fatal - nest #',i3,' not to upper boundary of'  &
                ,' coarser mesh or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif
     endif
  enddo
    

  ! this need to be done here since varfiles are filled before opspec3
  if (vwaittot.lt.vwait1) then
     print*,'total wait time must be <= individual varfile wait'
     print*,'      resetting vwaittot to ',vwait1
     vwaittot=vwait1
  endif

  ! stop the run if there are any fatal errors.  list how many
  !   warning and informative errors.


  if (ifaterr > 0) then
     print*,' -----------opspec2--------------------------'
     print*,' fatal     errors - ',ifaterr
     print*,' -----------------------------------------------'
     call abort_run('Fatal errors at namelist.','opspec2','opspec.f90')
  end if
end subroutine opspec2

! ********************************************************************

subroutine opspec3

  use mem_varinit
  use mem_grid
  use micphys
  use io_params
  use mem_radiate
  use mem_cuparm
  use mem_turb
  use mem_leaf
  use leaf_coms, only :     &
          ubmin,            & ! intent(in)
          ugbmin,           & ! intent(in)
          ustmin,           & ! intent(in)
          lc_gamm => gamm,  & ! intent(in)
          lc_gamh => gamh,  & ! intent(in)
          tprandtl,         & ! intent(in)
          ribmax,           & ! intent(in)
          leaf_maxwhc,      & ! intent(in)
          min_patch_area,   & ! intent(in)
          nstyp,            & ! intent(in)
          nvtyp             ! ! intent(in)
          
  use therm_lib , only:  &
          level          ! ! intent(in)
  use grell_coms, only:  &
          closure_type,  & ! intent(in)
          cap_maxs,      & ! intent(in)
          cld2prec,      & ! intent(in)
          maxclouds,     & ! intent(in)
          depth_min,     & ! intent(in)
          maxens_lsf,    & ! intent(in)
          maxens_eff,    & ! intent(in)
          maxens_dyn,    & ! intent(in)
          maxens_cap,    & ! intent(in)
          iupmethod,     & ! intent(in)
          radius,        & ! intent(in)
          zkbmax,        & ! intent(in)
          max_heat,      & ! intent(in)
          zcutdown,      & ! intent(in)
          z_detr         ! ! intent(in)

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM ! INTENT(IN)
  use mem_emiss, only: ichemi, isource ! INTENT(IN)

  ! CAT
  use catt_start, only: CATT ! INTENT(IN)

  ![MLO - mass check and exner function check
  use mem_mass, only : iexev, imassflx
  
  use mem_basic, only : ico2, co2con

  use mem_mnt_advec, only : iadvec
  implicit none

  integer :: ip,k,ifaterr,iwarerr,infoerr,ng,ngr,nc,nzz
  character(len=*), parameter :: h="**(opspec3)**"
  logical :: grell_on
  integer :: grella, grellz

  ifaterr=0
  iwarerr=0
  infoerr=0

  ! TEB_SPM
  !##########################################################################
  !EDF - Adition tho check if isource is activated if ichemi is
  !##########################################################################  
  if (TEB_SPM==1) then
     if (ichemi==1) then
        if (isource==0) then
           print*, 'FATAL - The SPM can not be activated without sources.'
           print*, '        ISOURCE must be equal to 1.'
           IFATERR = IFATERR + 1
        endif
     endif
  endif
  !##########################################################################

  ! CATT
  if (CATT==1) then
     ! Consistency in CATT
     ! Checking the tracers
     if (naddsc < 4) then
        print*, 'FATAL - If using CATT, the variable NADDSC must be >= 4.'
        IFATERR = IFATERR + 1
     endif
![MLO - make sure that CATT would work. Maybe these aren't strong requirements.
     do ng=1,ngrids
       if (ndeepest(ng) == 1 .or. ndeepest(ng) == 3) then
         print *, 'FATAL - You cannot run Kuo or Old Grell deep cumulus parameterization closure with CATT.'
         print *, 'Change ndeepest to 0 (off) or 2 (Grell).'
         IFATERR=IFATERR+1
       end if
       if (nshallowest(ng) == 1 .or. nshallowest(ng) == 3) then
         print *, 'FATAL - You cannot run Souza or old Grell shallow cumulus parametrization with CATT'
         print *, 'Change nshallowest to 0 (off) or 2 (Grell).'
         IFATERR=IFATERR+1
       end if
       if (nclouds > 2) then
         print *, 'FATAL - CATT expects up to two clouds only.'
         print *, 'Make your NCLOUDS to 1 or 2.'
         IFATERR=IFATERR+1
       end if
     end do
   end if
   
   !---------------------------------------------------------------------------------------!
   !    Checking whether the combination of akmin,akmax,hgtmin, hgtmax makes sense...      !
   !---------------------------------------------------------------------------------------!
   do ng=1,ngrids
      if (akmin(ng) <= 0.) then
         print *, 'FATAL - AKMIN must be positive'
         print *, 'Please change your setup for grid ',ng,'...'
         IFATERR=IFATERR+1
      end if

      if (akmax(ng) > 0.) then
         if (akmax(ng) < akmin(ng)) then
            print *, 'FATAL - If AKMAX is positive, then it must be greater than AKMIN...'
            print *, 'Please change your setup for grid ',ng,'...'
            IFATERR=IFATERR+1
         end if 
         if (hgtmax(ng) <= 0.) then
            print *, 'FATAL - If AKMAX is positive, HGTMAX must be positive as well...'
            print *, 'Please change your setup for grid ',ng,'...'
            IFATERR=IFATERR+1
         end if
         if (hgtmin(ng) < 0.) then
            print *, 'FATAL - If AKMAX is positive, HGTMIN must be non-negative...'
            print *, 'Please change your setup for grid ',ng,'...'
            IFATERR=IFATERR+1
         end if
         if (hgtmin(ng) >= hgtmax(ng)) then
            print *, 'FATAL - If AKMAX is positive, HGTMIN must less than HGTMAX...'
            print *, 'Please change your setup for grid ',ng,'...'
            IFATERR=IFATERR+1
         end if
      elseif (akmax(ng) < 0.) then
         print *, 'FATAL - AKMAX can''t be negative...'
         print *, 'Please change your setup for grid ',ng,'...'
         IFATERR=IFATERR+1
      end if
   end do

   !---------------------------------------------------------------------------------------!
   !   Determining whether Grell will be called or not, so the Grell-related               !
   ! configuration tests will be performed only if Grell is to be called.                  !
   !---------------------------------------------------------------------------------------!
   grell_on = any(nnqparm(1:ngrids) > 0)                                                   &
              .and. (nclouds > 2 .or. any(ndeepest(1:ngrids) == 2)                         &
                    .or. any(nshallowest(1:ngrids) == 2))

   ! Making sure that there aren't more clouds than the maximum
   if (nclouds > maxclouds) then
      print *, 'FATAL - Too many clouds, reduce nclouds'
      print *, 'Please change your setup for grid ',ng,'...'
      IFATERR=IFATERR+1
   end if
   ! Making sure that there aren't more clouds than the maximum
   do ng=1,ngrids
      if (nclouds < 1 .and. nnqparm(ng) > 0) then
         print *, 'FATAL - You need at least 1 cloud'
         print *, 'Please change your setup for grid ',ng,'...'
         IFATERR=IFATERR+1
      end if
      if (nnqparm(ng) > 0 .and. (ndeepest(ng) < 0 .or. ndeepest(ng) > 3)) then
         print *, 'FATAL - Ndeepest is out of range. Valid options are 0,1,2, or 3.'
         print *, 'Please change your setup for grid ',ng,'...'
         IFATERR=IFATERR+1
      end if 
      if (nnqparm(ng) > 0 .and. (nshallowest(ng) < 0 .or. nshallowest(ng) > 3)) then
         print *, 'FATAL - Nshallowest is out of range. Valid options are 0,1,2, or 3.'
         print *, 'Please change your setup for grid ',ng,'...'
         IFATERR=IFATERR+1
      end if 
   end do
   !  Blocking Grell convection without TKE. In the future this can be relaxed by 
   !  imposing iupmethod to be 1 in case the user wants to run idiffk=2 or 3. 
   do ng=1,ngrids
       if (grell_on) then
          select case (idiffk(ng))
          case(2,3)
            print *, 'FATAL - Grell cumulus requires turbulence with TKE (1,4,5,6,7,8)'
            print *, 'Please change your setup for grid ',ng,', currently set to '         &
                     ,idiffk(ng),'...'
            IFATERR=IFATERR+1
          end select
       end if
   end do
  if (any(nnqparm > 0) .and. confrq <= 0.) then
     print *, 'FATAL - CONFRQ must be positive when cumulus is on... '
     print *, '        Your CONFRQ= ',confrq
     IFATERR = IFATERR + 1
  end if

  if (grell_on) then
     do nc=1,nclouds-1
        if (radius(nc) < radius (nc+1)) then
           print *, 'FATAL - Cloud radii must be in decreasing sequence'
           print *, '        Deeper= ',radius(nc),' and Shallower=',radius(nc+1)
           IFATERR=IFATERR+1
        end if
     end do

     select case (closure_type)
     case ('en','nc','gr','lo','mc','kf','as')
       continue
     case default
       print *, 'FATAL - Invalid closure_type for Grell''s convection.'
       print *, 'Yours is currently set to ',closure_type
       IFATERR=IFATERR+1
     end select

     if (maxens_lsf <= 0) then
       print *, 'FATAL - maxens_lsf must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',maxens_lsf,'...'
       IFATERR=IFATERR+1
     end if
     if (maxens_eff <= 0) then
       print *, 'FATAL - maxens_eff must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',maxens_eff,'...'
       IFATERR=IFATERR+1
     end if
     if (maxens_cap <= 0) then
       print *, 'FATAL - maxens_cap must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',maxens_cap,'...'
       IFATERR=IFATERR+1
     end if


     if (iupmethod < 1 .or. iupmethod > 5) then
         print *, 'FATAL - If Cumulus parameterization is used, iupmethod must be between 1 and 4.'
         print *, 'Yours is currently set to ',iupmethod
         IFATERR=IFATERR+1
     elseif(iupmethod == 5 .and. cap_maxs > 0.) then
         print *, 'FATAL - updraft method 5 cannot be applied when cap_maxs > 0.'
         IFATERR=IFATERR+1
     end if

     if (cld2prec < 0. .or. cld2prec > 1.0) then
       print *, 'FATAL - cld2prec must be between 0 and 1.'
       print *, 'Yours is currently set to ',cld2prec
       IFATERR=IFATERR+1
     end if
     if (zkbmax <= 0.) then
       print *, 'FATAL - zkbmax must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',zkbmax
       IFATERR=IFATERR+1
     end if
     if (zcutdown <= 0.) then
       print *, 'FATAL - zcutdown must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',zcutdown
       IFATERR=IFATERR+1
     end if
     if (z_detr <= 0.) then
       print *, 'FATAL - z_detr must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',z_detr
       IFATERR=IFATERR+1
     end if
     if (max_heat <= 0.) then
       print *, 'FATAL - max_heat must be positive when Cuparm is activated.'
       print *, 'Yours is currently set to ',max_heat
       IFATERR=IFATERR+1
     end if
     !----------------------------------------------------------------------------------!
     !     CAP_MAXS must be always non-zero. In addition, it is allowed to be negative  !
     ! (fractional area method) only when when sigma-w is computed.                     !
     !----------------------------------------------------------------------------------!
     if (cap_maxs == 0.) then
        print *, 'FATAL - cap_maxs must be non-zero when cuparm is activated.'
        ifaterr = ifaterr + 1
     elseif (cap_maxs < 0) then
        do ng=1,ngrids
           if (idiffk(ng) /= 1 .and. idiffk(ng) /= 7 .and. idiffk(ng) /= 8) then
              print *, 'FATAL - cap_maxs(nc) can''t be < 0. if sigma-w is unavailable.'
              print *, 'Yours is currently set to ',cap_maxs
              print *, 'And your turbulence is set to ',idiffk(ng),' for grid ',ng
              ifaterr = ifaterr + 1
           elseif(ndeepest(ng)    == 3 .or. &
                  (nclouds > 1 .and. nshallowest(ng) == 3)) then
              print *, 'FATAL - cap_maxs can''t be < 0. for old Grell.'
              print *, 'Yours is currently set to ',cap_maxs
              print *, 'And your turbulence is set to ',idiffk(ng),' for grid ',ng
           end if
        end do
     end if

     do nc=1,nclouds
        if (depth_min(nc) <= 0.) then
          print *, 'FATAL - depth_min(nc) must be positive.'
          print *, 'Your is currently set to ',depth_min(nc),' for type ',nc
          IFATERR=IFATERR+1
        end if
     end do
  end if


  ! check that moisture is turned on if radiation is used.
  !   (severity - f)

  if(ilwrtyp+iswrtyp.gt.0.and.level.eq.0)then
     print*,' fatal  - radiation scheme must be run with moisture.'
     ifaterr=ifaterr+1
  endif

  ! microphysics flags and parameter settings

  if((irain.ge.2.and.irain.le.4.and.rparm.le.0.)  &
       .or.(icloud.ge.2.and.icloud.le.5.and.cparm.le.0.)  &
       .or.(ipris.ge.2.and.ipris.le.4.and.pparm.le.0.)  &
       .or.(isnow.ge.2.and.isnow.le.4.and.sparm.le.0.)  &
       .or.(igraup.ge.2.and.igraup.le.4.and.gparm.le.0.)  &
       .or.(iaggr.ge.2.and.iaggr.le.4.and.aparm.le.0.)  &
       .or.(ihail.ge.2.and.ihail.le.4.and.hparm.le.0.)) then
     print 26,ng,rparm,pparm,sparm,gparm,aparm,hparm
26   format (' fatal - microphysics - xparm must be positive'  &
          ,' if micro flags are set to 2, 3, or 4,'  &
          ,' or up to 5 for icloud. ',i3,5f10.7)
     ifaterr=ifaterr+1
  endif

  ! convective parameterization flags and parameter settings

  do ng=1,ngrids
     if(nnqparm(ng).gt.0.and.level.eq.0) then
        print 27
27      format (' fatal - level must be at least'  &
             ,' 1 for the cumulus parameterization')
        ifaterr=ifaterr+1
     endif
  enddo

  ! moving grids and topography interpolation

  do ng=1,ngrids
     if((abs(gridu(ng)) .gt. 1.e-20 .or. abs(gridv(ng)) .gt. 1.e-20)  &
          .and. itoptflg(ng) .ne. 0) then
        print 28
28      format (' fatal - nested grid topography must be interpolated'  &
             ,' from parent grid if nested grid moves')
        ifaterr=ifaterr+1
     endif
  enddo

  ! check horizontal and vertical grid spacings.

  if(dzmax.lt.deltaz)then
     print*,' warning - deltaz is being reduced by a low value',  &
          ' of dzmax.'
     iwarerr=iwarerr+1
  endif

  if(dzrat.gt.1.2)then
     print*,' warning - large vertical stretch ratios sacrifice',  &
          ' second order accuracy in the vertical differencing.'
     iwarerr=iwarerr+1
  endif

  ! check numerical schemes.

  if((sspct.lt.0.2.or.sspct.gt.1.0).and.sspct.ne.0.0)then
     print*,' warning - sspct should normally range from 0.2 to 1.0'
     iwarerr=iwarerr+1
  endif

  if(nfpt.gt.nnzp(1))then
     print*,' fatal - nfpt must be less than nnzp(1).'
     ifaterr=ifaterr+1
  endif

  if(iadvl.ne.2.and.iadvl.ne.4)then
     print*,' fatal - iadvl must be 2 or 4'
     ifaterr=ifaterr+1
  endif

  if(iadvf.ne.2.and.iadvf.ne.6)then
     print*,' fatal - iadvf must be 2 or 6'
     ifaterr=ifaterr+1
  endif

  ! check turbulence parameterization.

  do ngr=1,ngrids
     if(idiffk(ngr) < 1.or.idiffk(ngr) > 8)then
        print*,' fatal - idiffk must be between 1 and 8.'
        ifaterr=ifaterr+1
     endif
  enddo

  if(ibruvais < 1 .or. ibruvais > 3)then
     print*,' fatal - ibruvais must be either 1, 2, or 3. Yours is set to ',ibruvais,'...'
     ifaterr=ifaterr+1
  end if

  if(ibotflx < 0 .or. ibotflx > 1)then
     print*,' fatal - ibotflx must be either 0 or 1. Yours is set to ',ibotflx,'...'
     ifaterr=ifaterr+1
  end if

  ! check that diffusion flags are compatible if using ihorgrad=1
  if(ihorgrad.eq.2)then
     select case (idiffk(ngr))
     case (3:6)
        print*,' fatal - can''t use ihorgrad=2 if idiffk 3, 4, 5 or 6'
        ifaterr=ifaterr+1
     end select
  endif
  
  
  if (ubmin < 1.e-4 .or. ubmin > 1.0) then
     print *, 'FATAL - UBMIN must be between 0.0001 and 1.0'
     ifaterr = ifaterr + 1
  end if

  if (ustmin < 1.e-4 .or. ustmin > 1.0) then
     print *, 'FATAL - USTMIN must be between 0.0001 and 1.0'
     ifaterr = ifaterr + 1
  elseif (ustmin > ubmin) then
     print *, 'FATAL - USTMIN must be less than UBMIN.  Check your namelist'
  end if

  if (ugbmin < ustmin .or. ugbmin > ubmin) then
     print *, 'FATAL - UGBMIN must be between USTMIN and UBMIN'
     ifaterr = ifaterr + 1
  end if


  
  if (lc_gamm < 0.1 .or. lc_gamm > 100.0) then
     print *, 'FATAL - GAMM must be between 0.1 and 100.'
     ifaterr = ifaterr + 1
  end if
  
  if (lc_gamh < 0.1 .or. lc_gamh > 100.0) then
     print *, 'FATAL - GAMH must be between 0.1 and 100.'
     ifaterr = ifaterr + 1
  end if
  
  if (tprandtl < 0.01 .or. tprandtl > 100.0) then
     print *, 'FATAL - TPRANDTL must be between 0.01 and 100.'
     ifaterr = ifaterr + 1
  end if
  
  if (ribmax < 0.01 .or. ribmax > 20.0) then
     print *, 'FATAL - RIBMAX must be greater than 0.01 and less than 20.'
     print *, ribmax
     ifaterr = ifaterr + 1
  end if
  
  if (leaf_maxwhc < 0.0 .or. leaf_maxwhc > 10.0) then
     print *, 'FATAL - LEAF_MAXWHC must be greater than 0.0 and less than 10.'
     print *, leaf_maxwhc
     ifaterr = ifaterr + 1
  end if

  if (isoilbc < 0 .or. isoilbc > 4) then
     write (unit=*,fmt='(a,1x,i4,a)')                                                      &
       'Invalid ISOILBC, it must be between 0 and 4. Yours is set to',isoilbc,'...'
     ifaterr = ifaterr +1
  end if
 
  if (ipercol < 0 .or. ipercol > 2) then
     write (unit=*,fmt='(a,1x,i4,a)')                                                      &
       'Invalid IPERCOL, it must be between 0 and 2. Yours is set to',ipercol,'...'
     ifaterr = ifaterr +1
  end if
   
  if (runoff_time < 0.0) then
     write (unit=*,fmt='(a,1x,es14.7,a)')                                                  &
           'Invalid RUNOFF_TIME, it can''t be negative. Yours is set to',runoff_time,'...'
     ifaterr = ifaterr +1
  end if

  if (iadvec < 1 .or. iadvec > 2) then 
     print *, 'FATAL - IADVEC must be either 1 or 2.'
     ifaterr=ifaterr+1
  else if (iadvec == 2 .and. if_adap /= 0) then
     print *, 'FATAL - Monotonic advection is only allowed with sigma-z coordinates...'
     print *, '        Either set iadvec to 1 or if_adap to 0!'
     ifaterr=ifaterr+1
  end if

![MLO - Some extra checks for mass and Medvidy's fix on Exner tendency
! Complete Exner tendency and vertical coordinate.
  if (iexev < 1 .or. iexev > 2) then 
    print *, 'FATAL - IEXEV must be either 1 or 2.'
    ifaterr=ifaterr+1 
  elseif (iexev == 2 .and. if_adap /= 0) then 
    print *, 'FATAL - IEXEV cannot be set to 2 with adaptive coordinate'
    ifaterr=ifaterr+1 
  end if

  ! Just adding a warning message that cumulus parameterization feedback will be ignored 
  ! because the user didn't set iswrtyp or ilwrtyp to 3 (Harrington);
  if (iswrtyp == 4 .and. icumfdbk /= 0) then
    print *, '------------------------------------------------------------------------'
    print *, ' FATAL - CARMA shortwave radiation won''t work with partial cumulus     '
    print *, '         scheme. The ice mixing ratio in cumulus clouds can be much     '
    print *, '         larger than the resolved clouds, and bands 32-35 don''t        '
    print *, '         give reasonable results.  This may be just a bug and it will   '
    print *, '         be eventually revisited and addressed, but for the time being  '
    print *, '         either use Harrington or turn off the cumulus feedback.        '
    print *, '------------------------------------------------------------------------'
    ifaterr = ifaterr + 1
  elseif (iswrtyp /= 3 .and. icumfdbk /=0) then
    print *, '------------------------------------------------------------------------'
    print *, 'INFO - Shortwave radiation won''t have cumulus parameterization feedback'
    print *, '       You should use Harrington scheme to have this effect.            '
    print *, '------------------------------------------------------------------------'
    print *, ' '
    infoerr = infoerr + 1
  end if
  if (ilwrtyp /= 3 .and. ilwrtyp /= 4 .and. icumfdbk /=0) then
    print *, '------------------------------------------------------------------------'
    print *, 'INFO - Longwave radiation won''t have cumulus parameterization feedback '
    print *, '       You should use Harrington scheme or CARMA to have this effect.'
    print *, '------------------------------------------------------------------------'
    print *, ' '
    infoerr = infoerr + 1
  end if
  do ng=1,ngrids
    if (nnqparm(ng) == 0 .and. icumfdbk == 1 .and. & 
        (iswrtyp == 3 .or. ilwrtyp == 3)) then
       print *, '-------------------------------------------------------------------------'
       print *, 'INFO - Cumulus parameterization will have no effect on Harrington scheme '
       print *, '       on grid',ng,' because the cumulus parameterizations are off.      '
       print *, '-------------------------------------------------------------------------'
       print *, ' '
       infoerr = infoerr + 1
    end if 
    if (nnqparm(ng) == 0 .and. icumfdbk == 1 .and. & 
        (ilwrtyp == 4)) then
       print *, '-------------------------------------------------------------------------'
       print *, 'INFO - Cumulus parameterization will have no effect on CARMA             '
       print *, '       on grid',ng,' because the cumulus parameterizations are off.      '
       print *, '-------------------------------------------------------------------------'
       print *, ' '
       infoerr = infoerr + 1
    end if 
  end do

  select case (isfcl)
  case (0:2,5)
     continue
  case (3)
     print*,' fatal - SSiB is not available in this version, use LEAF or ED instead...'
     ifaterr=ifaterr+1
  case (4)
     print*,' fatal - GEMTM is not available in this version, use LEAF or ED instead...'
     ifaterr=ifaterr+1
  case default
     print*,' fatal - isfcl must be either 1, 2, or 5. Yours is set to ',isfcl,'...'
     ifaterr=ifaterr+1
  end select
  
  ! Check whether the surface layer exchange scheme the user chose is okay. 
  if (isfcl /= 0 .and. (istar < 1 .or. istar > 4)) then
     print *, 'fatal - ISTAR must be between 1 and 4, and yours is set to ',istar,'...'
     ifaterr=ifaterr+1
  end if
  
  ! Check whether the ground vapour method scheme the user chose is fine or not. 
  if (isfcl /= 0 .and. (igrndvap < 0 .or. igrndvap > 5)) then
     print *, 'fatal - IGRNDVAP must be between 0 and 5, and yours is set to '            &
            ,igrndvap,'...'
     ifaterr=ifaterr+1
  end if

  ! Check whether the time step makes sense for LEAF or ED.
  if (isfcl /= 0 .and. dtleaf <= 0.) then
     print *, 'fatal - DTLEAF must be positive, and yours is set to ',dtleaf,'...'
     ifaterr=ifaterr+1
  elseif (isfcl /= 0 .and. dtleaf > 30.) then
     print *, 'fatal - DTLEAF must be less than 30 sec, and yours is set to ',dtleaf,'...'
     ifaterr=ifaterr+1
  end if

  ! check whether the soil model will be run and make sure that the
  !   number of soil levels are correct.(severity - f,i )

  if(isfcl.eq.0.and.nzg.gt.1)then
     print*,' info  - more soil levels specified than needed.'
     infoerr=infoerr+1
  endif

  if (isfcl == 0 .and. npatch /= 2) then
     print*, ' fatal  - when isfcl = 0, npatch must be 2.'
     ifaterr = ifaterr + 1
  else if (npatch < 2) then
     print*, ' fatal - npatch must be at least 2.'
     ifaterr = ifaterr + 1
  endif

  if(isfcl.gt.0.and.nzg.le.2)then
     print*,  &
          ' fatal  - at least 2 soil levels are needed for soil'  &
          ,' model.'
     ifaterr=ifaterr+1
  endif

  !----------------------------------------------------------------------------------------!
  !     Check whether the number of patches make sense.                                    !
  !----------------------------------------------------------------------------------------!
  select case (isfcl)
  case (1,2)
     if (npatch > nvtyp+1) then
        print *, ' fatal - Running LEAF-3, and NPATCH is greater than the number of '//    &
                          'vegetation types.'
        ifaterr=ifaterr+1
     end if
  case (5)
     if (npatch > nstyp+1) then
        print *, ' fatal - Running ED-2, and NPATCH is greater than the number of '//      &
                          'vegetation types.'
        ifaterr=ifaterr+1
     end if
     if (nvegpat /= npatch - 1) then
        print *, ' fatal - Running ED-2, NVEGPAT must be equal to NPATCH-1'
        ifaterr=ifaterr+1
     end if
  end select
  !----------------------------------------------------------------------------------------!



  !----------------------------------------------------------------------------------------!
  !     Check whether the minimum area is reasonable.                                      !
  !----------------------------------------------------------------------------------------!
  if (min_patch_area < 0.0001 .or. min_patch_area > 0.1) then
     print *, ' fatal - MIN_PATCH_AREA must be between 0.0001 and 0.1'
     ifaterr = ifaterr + 1
  end if
  !----------------------------------------------------------------------------------------!


  !----- Check CO2 settings. --------------------------------------------------------------!
  select case (ico2)
  case (0,1)
     !----- Value is okay, now check whether the initial co2con makes sense... ------------!
     if (co2con(1) < 0.) then
        print *, 'FATAL - When ICO2 < 2, CO2CON(1) cannot be negative and yours is '       &
               ,co2con(1),'...'
        ifaterr = ifaterr + 1
     elseif (co2con(1) > 1.e6) then
        print *, 'FATAL - When ICO2 < 2, CO2CON(1) cannot be more than 1.E6 and yours is ' &
               ,co2con(1),'...'
        ifaterr = ifaterr + 1
     end if
  case (2)
     !----- Value is okay, now check whether the initial co2con makes sense... ------------!
     nzz = maxval(nnzp(1:ngrids))
     if (any(co2con(1:nzz) < 0.)) then
        print *, 'FATAL - When ICO2=2, CO2CON must be positive for all levels...'
        print *, 'Your values: ',co2con(1:nzz)
        ifaterr = ifaterr + 1
     elseif (any(co2con(1:nzz) > 1.e6)) then
        print *, 'FATAL - When ICO2=2, CO2CON must be less than 1 million for all levels...'
        print *, 'Your values: ',co2con(1:nzz)
        ifaterr = ifaterr + 1
     end if
  case (3)
        print *, 'FATAL - Assimilation of CO2 through RALPH2 files is not available yet.'
        ifaterr = ifaterr + 1
  case default 
        print *, 'FATAL - Invalid ICO2... Yours is set to ',ico2,'...'
        ifaterr = ifaterr + 1
  end select

  do k=1,nzg
     if (slz(k) > -.001) then
        write (unit=*,fmt='(a,200(es12.5,1x))') 'SLZ = ',slz(1:nzg)
        write (unit=*,fmt='(a,1x,i5,1x,a,a)') 'FATAL - soil level',k,'is greater than '    &
                                             ,'-0.001 m.  Make it deeper...'
        ifaterr=ifaterr+1
     endif
  enddo


  do k=1,nzg-1
     if (slz(k)-slz(k+1) > .001) then
        write (unit=*,fmt='(a)') ' FATAL - Soil layers must be defined from bottom to top.'
        write (unit=*,fmt='(a,1x,i5,1x,a,1x,es12.5)') '   * Layer=',k  ,' SLZ =',slz(k  )
        write (unit=*,fmt='(a,1x,i5,1x,a,1x,es12.5)') '   * Layer=',k+1,' SLZ =',slz(k+1)
        ifaterr=ifaterr+1
     end if
  end do


  do k=1,nzg
     if (slmstr(k) < -1. .or. slmstr(k) > 2.) then
        write (unit=*,fmt='(a,200(es12.5,1x))') 'SLMSTR = ',slmstr(1:nzg)
        write (unit=*,fmt='(a,1x,i5,1x,a)') ' FATAL - Soil moisture index in layer',k      &
                                           ,'is outside the allowed range [-1.0;2.0].'
        ifaterr=ifaterr+1
     end if
  end do

  ! if the soil model will be run with no radiation, make a suggestion
  !   that the radiation be turned on. (severity - f )

  do ngr=1,ngrids
     if(isfcl.gt.0.and.ilwrtyp+iswrtyp.eq.0)then
        print*,' fatal  - radiation scheme must be run with soil',  &
             ' model.'
        ifaterr=ifaterr+1
     endif
  enddo


  ! make sure that if nudging, nudging time scales are greater than
  ! the model coarse grid timestep, and that rayleigh friction nudging
  ! is not done with variable initialization.

  if (initial .eq. 1) then
     if (nfpt .gt. 0 .and. distim .gt. 0. .and.  &
          distim .lt. dtlongn(1)) then
        print*, 'rayleigh friction nudging is being done'
        print*, 'and distim is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif
  endif

  if (initial .eq. 2) then

     if (nfpt .gt. 0 .and. distim .gt. 0.) then
        print*, 'rayleigh friction nudging may not be used when'
        print*, 'variable initialization is used.'
        ifaterr=ifaterr+1
     endif

     if (nudlat .ge. 1 .and. tnudlat .gt. 0. .and.  &
          tnudlat .lt. dtlongn(1)) then
        print*, 'lateral boundary nudging is being done'
        print*, 'and tnudlat is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif

     if (tnudcent .gt. 0. .and. tnudcent .lt. dtlongn(1)) then
        print*, 'center nudging is being done'
        print*, 'and tnudcent is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif
     if (tnudtop .gt. 0. .and. tnudtop .lt. dtlongn(1)) then
        print*, 'top boundary nudging is being done'
        print*, 'and tnudtop is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif

  endif


  !     check the averaging and analysis frequencies for consistency.

  if (abs(avgtim).gt.0.0.and.frqmean.le.0.0.and.frqboth.le.0.) then
     print*,'have frqmean=0 & frqboth=0 even though avgtim=',avgtim
     print*,'respecifying avgtim=0.'
     avgtim=0.
     iwarerr=iwarerr+1
  endif
  if (frqlite.gt.0.) then
     !   if ( nl3d(1)+nl2d(1)+nl3ds(1).eq.0)then
     !      print*,'have no lite variables even though frqlite=',frqlite
     !      print*,'respecify in vtables or set frqlite=0 in namelist'
     !      ifaterr=ifaterr+1
     !   endif
  endif
  if (frqmean.gt.0.0.and.abs(avgtim).gt.0.) then
     if ( abs(avgtim).gt.frqmean ) then
        print*,'avgtim must be <= frqmean'
        ifaterr=ifaterr+1
     endif
     !   if ( nm3d(1)+nm2d(1)+nm3ds(1).eq.0)then
     !      print*,'have no mean variables even though frqmean=',frqmean
     !      print*,'respecify in vtables or set frqmean=0 in namelist'
     !      ifaterr=ifaterr+1
     !   endif
  endif
  if (frqboth.gt.0.0.and.abs(avgtim).gt.0.) then
     if ( abs(avgtim).gt.frqboth ) then
        print*,'avgtim must be <= frqboth'
        ifaterr=ifaterr+1
     endif
     !   if ( nb3d(1)+nb2d(1)+nb3ds(1).eq.0)then
     !      print*,'have no both variables even though frqboth=',frqboth
     !      print*,'respecify in vtables or set frqboth=0 in namelist'
     !      ifaterr=ifaterr+1
     !   endif
  endif
  if (frqmean.gt.0.0.and.frqboth.gt.0.0.and.abs(avgtim).gt.0.) then
     if ( (frqmean.gt.frqboth.and.mod(frqmean,frqboth).ne.0.).or.  &
          (frqmean.lt.frqboth.and.mod(frqboth,frqmean).ne.0.) ) then
        print*,'frqmean must be a multiple of frqboth or vice versa'
        ifaterr=ifaterr+1
     endif
  endif

  ! check printout parameters.

  do ngr=1,ngrids
     do ip=1,nplt
        if(ixsctn(ip).eq.1.and.isbval(ip).gt.nnyp(ngr))then
           print 1,ip,ngr
1          format (' fatal - isbval(',i2,') is out of bounds in'  &
                ,' y-direction for grid number ',i2,'.')
           ifaterr=ifaterr+1
        elseif(ixsctn(ip).eq.2.and.isbval(ip).gt.nnxp(ngr))then
           print 2,ip,ngr
2          format (' fatal - isbval(',i2,') is out of bounds in'  &
                ,' x-direction for grid number ',i2,'.')
           ifaterr=ifaterr+1
        elseif(ixsctn(ip).eq.3.and.isbval(ip).gt.nnzp(ngr))then
           print 3,ip,ngr
3          format (' fatal - isbval(',i2,') is out of bounds in'  &
                ,' z-direction for grid number ',i2,'.')
           ifaterr=ifaterr+1
        endif
     enddo
  enddo

  ! stop the run if there are any fatal errors.  list how many
  !   warning and informative errors.


  if (ifaterr > 0) then
     print*,' -----------opspec3--------------------------'
     print*,' fatal     errors - ',ifaterr
     print*,' warning   errors - ',iwarerr
     print*,' inform  messages - ',infoerr
     print*,' -----------------------------------------------'
     call abort_run('Fatal errors at namelist.','opspec3','opspec.f90')
  end if
end subroutine opspec3
