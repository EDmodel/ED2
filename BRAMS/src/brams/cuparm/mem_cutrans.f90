module mem_cutrans

  !Function: Share and Allocate variables to CUPARM_GRELL routines
  !Author: Luiz Flavio 	Date: 09-25-2003 	CPTEC/INPE
  !

  use mem_grell_param

  implicit none

  logical :: cutrans_alloc=.false.
  integer :: iruncon = 0

  !5d dependece  (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
  !ngrids_cp = numero de grades onde a parametrizacao de cumulus e' usada
  real,  allocatable :: zcup5d(:,:,:,:,:),  & ! z level
       pcup5d(:,:,:,:,:),  & !p level
       prup5d(:,:,:,:,:),  &
       clwup5d(:,:,:,:,:), &
       tup5d(:,:,:,:,:)

contains

  subroutine alloc_cutrans()

    use mem_grell_param, only: mgmxp,mgmyp,maxiens,ngrids_cp

    implicit none

    if(cutrans_alloc) then
       print *,'ERROR: cutrans already allocated'
       print *,'Routine: alloc_cutrans File: var_cutrans.f90'
       print *,'Dir: .../shared/tools/brams20/src/rams/5.04/modules'
       stop
    end if

    allocate(zcup5d(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp),       &
         pcup5d(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp), &
         prup5d(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp), &
         clwup5d(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp), &
         tup5d(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp))

    cutrans_alloc=.true.
  end subroutine alloc_cutrans

  subroutine zero_cutrans()

    implicit none

    zcup5d=0.
    pcup5d=0.
    prup5d=0.
    clwup5d=0.
    tup5d=0.

  end subroutine zero_cutrans

end module mem_cutrans
