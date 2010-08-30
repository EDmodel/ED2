!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module io_params

  use grid_dims

  integer, parameter :: maxlite=80

  character(len=32)  :: lite_vars(maxlite)

  character(len=80)  :: afilin

  character(len=256) :: hfilout, afilout, pastfn, hfilin
  ! Modif. by ALF

  character(len=20)  :: xlite,ylite,zlite
  integer :: ipastin
  !----------------------------------------------------------------------------
  integer :: ioutput,iinput,iopunt,kwrite,ihistdel,iclobber,nlite_vars
  real    :: frqhis,frqanl,avgtim,frqlite,frqmean,frqboth  

  ! Increase timstr to a 8-byte and include namelist variables
  real(kind=8) :: timstr
  integer :: iyearh,imonthh,idateh,itimeh


  !----------------------------------------------------------------------------

  integer, dimension(maxgrds) :: itoptflg, isstflg, ivegtflg, isoilflg  &
       ,ndviflg, nofilflg, itopsflg, iz0flg
  !TEB
  integer, dimension(maxgrds) :: ifusflg

  real                        :: z0fact
  integer                     :: ntopsmth,izflat
  real, dimension(maxgrds)    :: z0max,toptenh,toptwvl
  !----------------------------------------------------------------------------
  character(len=256), dimension(maxgrds) :: itoptfn,isstfn,ivegtfn,isoilfn  &
       ,ndvifn
  ! TEB
  character(len=256), dimension(maxgrds) :: ifusfn

  ! Modif.by ALF
  !----------------------------------------------------------------------------

  integer, parameter                               :: maxsstfiles=2000
  character(len=256)                               :: sstfpfx, sfcfiles, &
       topfiles
  ! TEB
  character(len=256)                               :: fusfiles
  ! Modif.by ALF
  character(len=256),dimension(maxsstfiles,maxgrds) :: fnames_sst
  ! Modif.by ALF
  character(len=14),dimension(maxsstfiles,maxgrds) :: itotdate_sst
  ! Created by ALF
  character(len=14),dimension(maxgrds) :: lastdate_sst
  !----------------------------------------------------------------------------
  integer                             :: iupdsst,isstcyclic,isstcycdata
  integer,dimension(maxgrds)          :: nsstfiles,isstflp,isstflf
  real(kind=8),dimension(maxgrds)             :: ssttime1,ssttime2
  !----------------------------------------------------------------------------
  integer, parameter                                     :: maxndvifiles=2000
  character(len=256)                                     :: ndvifpfx
  ! Modif.by ALF
  character(len=256), dimension(maxndvifiles,maxgrds)     :: fnames_ndvi
  ! Modif.by ALF
  character(len=14), dimension(maxndvifiles,maxgrds)     :: itotdate_ndvi
  !----------------------------------------------------------------------------
  integer                             :: iupdndvi,indvicyclic,indvicycdata
  integer,dimension(maxgrds)          :: nndvifiles,indviflp,indviflf
  real(kind=8),dimension(maxgrds)             :: ndvitime1,ndvitime2

  !----------------------------------------------------------------------------
  character(len=8), dimension(50)  :: plfmt,pltit
  character(len=16), dimension(50) :: iplfld
  !----------------------------------------------------------------------------
  integer                :: nplt,initfld
  integer, dimension(50) :: ixsctn,iplvect,isbval,iaa,iab,joa,job,naavg,noavg
  real                   :: frqprt
  real, dimension(50)    :: plconlo,plconhi,plconin


end module io_params
