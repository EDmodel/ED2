MODULE mem_sib

  !itb_usp---------------S I B 2 ----------------------------------------------

  !itb_usp.....here are the arrays for using SiB in RAMS...


  INTEGER, PARAMETER:: nstyp_sib = 12, nvtyp_sib = 13             &
       ,nzg_sib = 7

  REAL, DIMENSION(nstyp_sib) :: bee_sib,phsat_sib,satco_sib       &
       ,poros_sib,slope_sib,wopt_sib,skew_sib,respsat_sib

  REAL, DIMENSION(nvtyp_sib,2,2) :: tran_sib, ref_sib
  REAL, DIMENSION(nvtyp_sib,2)   :: soref_sib

  REAL, DIMENSION(nvtyp_sib) :: z2_sib,z1_sib,fvcover_sib,chil_sib  &
       ,sodep_sib,rootd_sib,phc_sib,vmax0_sib,effcon_sib    &
       ,gslope_sib,gsmin_sib,atheta_sib,btheta_sib,trda_sib &
       ,trdm_sib,trop_sib,respcp_sib,slti_sib,shti_sib      &
       ,hlti_sib,hhti_sib

  REAL,DIMENSION(50) :: laig_sib,fvcg_sib

  REAL,DIMENSION(nvtyp_sib,50,50) :: a_zo_sib,a_zp_sib,a_rbc_sib    &
       ,a_rdc_sib

  REAL,DIMENSION(nvtyp_sib) :: zc_sib,zlw_sib,zlen_sib,ltmax_sib    &
       ,stem_sib,nd98_sib,nd02_sib,srmax_sib,srmin_sib

  !itb_usp.................................................................

  TYPE sib_brams

     ! Variables to be dimensioned by (nxp, nyp)
     REAL, POINTER, DIMENSION(:, :) :: pco2ap, pco2m, rst,    &
          capac1, capac2, snow1, snow2, fss, fws, assimn, respg,   &
          rstfac1,rstfac2,rstfac3,rstfac4,ect,eci,egi,egs,hc,hg,   &
          w1,w2,w3,ww1,ww2,ww3,exo,ta,tc,tg,td1,td2,td3,td4,td5,   &
          td6,ra,rb,rc,rd,roff,zlt,green,apar,nee,cu,ct,ventmf,    &
          pco2c,pco2i,pco2s,ea,sha,em,rha,radvbc,radvdc,radnbc,    &
          radndc,dlwbot,cp,rho,psy,cupr,lspr


  END TYPE sib_brams

  TYPE (sib_brams), ALLOCATABLE :: sib_brams_g(:), sib_bramsm_g(:)

CONTAINS

  SUBROUTINE alloc_sib_brams(sib, n2, n3)

    IMPLICIT NONE
    TYPE (sib_brams) :: sib
    INTEGER, INTENT(in) :: n2, n3

    ! Allocate arrays based on options (if necessary)

    ALLOCATE (sib%pco2ap(n2, n3))
    ALLOCATE (sib%pco2m(n2, n3))
    ALLOCATE (sib%rst(n2, n3))
    ALLOCATE (sib%capac1(n2, n3))
    ALLOCATE (sib%capac2(n2, n3))
    ALLOCATE (sib%snow1(n2, n3))
    ALLOCATE (sib%snow2(n2, n3))
    ALLOCATE (sib%fss(n2, n3))
    ALLOCATE (sib%fws(n2, n3))
    ALLOCATE (sib%assimn(n2, n3))
    ALLOCATE (sib%respg(n2, n3))
    ALLOCATE (sib%rstfac1(n2, n3))
    ALLOCATE (sib%rstfac2(n2, n3))
    ALLOCATE (sib%rstfac3(n2, n3))
    ALLOCATE (sib%rstfac4(n2, n3))
    ALLOCATE (sib%ect(n2, n3))
    ALLOCATE (sib%eci(n2, n3))
    ALLOCATE (sib%egi(n2, n3))
    ALLOCATE (sib%egs(n2, n3))
    ALLOCATE (sib%hc(n2, n3))
    ALLOCATE (sib%hg(n2, n3))
    ALLOCATE (sib%w1(n2, n3))
    ALLOCATE (sib%w2(n2, n3))
    ALLOCATE (sib%w3(n2, n3))
    ALLOCATE (sib%ww1(n2, n3))
    ALLOCATE (sib%ww2(n2, n3))
    ALLOCATE (sib%ww3(n2, n3))
    ALLOCATE (sib%exo(n2, n3))
    ALLOCATE (sib%ta(n2, n3))
    ALLOCATE (sib%tc(n2, n3))
    ALLOCATE (sib%tg(n2, n3))
    ALLOCATE (sib%td1(n2, n3))
    ALLOCATE (sib%td2(n2, n3))
    ALLOCATE (sib%td3(n2, n3))
    ALLOCATE (sib%td4(n2, n3))
    ALLOCATE (sib%td5(n2, n3))
    ALLOCATE (sib%td6(n2, n3))
    ALLOCATE (sib%ra(n2, n3))
    ALLOCATE (sib%rb(n2, n3))
    ALLOCATE (sib%rc(n2, n3))
    ALLOCATE (sib%rd(n2, n3))
    ALLOCATE (sib%roff(n2, n3))
    ALLOCATE (sib%zlt(n2, n3))
    ALLOCATE (sib%green(n2, n3))
    ALLOCATE (sib%apar(n2, n3))
    ALLOCATE (sib%nee(n2, n3))
    ALLOCATE (sib%cu(n2, n3))
    ALLOCATE (sib%ct(n2, n3))
    ALLOCATE (sib%ventmf(n2, n3))
    ALLOCATE (sib%pco2c(n2, n3))
    ALLOCATE (sib%pco2i(n2, n3))
    ALLOCATE (sib%pco2s(n2, n3))
    ALLOCATE (sib%ea(n2, n3))
    ALLOCATE (sib%sha(n2, n3))
    ALLOCATE (sib%em(n2, n3))
    ALLOCATE (sib%rha(n2, n3))
    ALLOCATE (sib%radvbc(n2, n3))
    ALLOCATE (sib%radvdc(n2, n3))
    ALLOCATE (sib%radnbc(n2, n3))
    ALLOCATE (sib%radndc(n2, n3))
    ALLOCATE (sib%dlwbot(n2, n3))
    ALLOCATE (sib%cp(n2, n3))
    ALLOCATE (sib%rho(n2, n3))
    ALLOCATE (sib%psy(n2, n3))
    ALLOCATE (sib%cupr(n2, n3))
    ALLOCATE (sib%lspr(n2, n3))


    RETURN
  END SUBROUTINE alloc_sib_brams

  SUBROUTINE nullify_sib_brams(sib)

    IMPLICIT NONE
    TYPE (sib_brams) :: sib

    IF (ASSOCIATED(sib%pco2ap))  NULLIFY (sib%pco2ap)
    IF (ASSOCIATED(sib%pco2m))   NULLIFY (sib%pco2m)
    IF (ASSOCIATED(sib%rst))     NULLIFY (sib%rst)
    IF (ASSOCIATED(sib%capac1))  NULLIFY (sib%capac1)
    IF (ASSOCIATED(sib%capac2))  NULLIFY (sib%capac2)
    IF (ASSOCIATED(sib%snow1))   NULLIFY (sib%snow1)
    IF (ASSOCIATED(sib%snow2))   NULLIFY (sib%snow2)
    IF (ASSOCIATED(sib%fss))   NULLIFY (sib%fss)
    IF (ASSOCIATED(sib%fws))   NULLIFY (sib%fws)
    IF (ASSOCIATED(sib%assimn))   NULLIFY (sib%assimn)
    IF (ASSOCIATED(sib%respg))   NULLIFY (sib%respg)
    IF (ASSOCIATED(sib%rstfac1))   NULLIFY (sib%rstfac1)
    IF (ASSOCIATED(sib%rstfac2))   NULLIFY (sib%rstfac2)
    IF (ASSOCIATED(sib%rstfac3))   NULLIFY (sib%rstfac3)
    IF (ASSOCIATED(sib%rstfac4))   NULLIFY (sib%rstfac4)
    IF (ASSOCIATED(sib%ect))   NULLIFY (sib%ect)
    IF (ASSOCIATED(sib%eci))   NULLIFY (sib%eci)
    IF (ASSOCIATED(sib%egi))   NULLIFY (sib%egi)
    IF (ASSOCIATED(sib%egs))   NULLIFY (sib%egs)
    IF (ASSOCIATED(sib%hc))   NULLIFY (sib%hc)
    IF (ASSOCIATED(sib%hg))   NULLIFY (sib%hg)
    IF (ASSOCIATED(sib%w1))   NULLIFY (sib%w1)
    IF (ASSOCIATED(sib%w2))   NULLIFY (sib%w2)
    IF (ASSOCIATED(sib%w3))   NULLIFY (sib%w3)
    IF (ASSOCIATED(sib%ww1))   NULLIFY (sib%ww1)
    IF (ASSOCIATED(sib%ww2))   NULLIFY (sib%ww2)
    IF (ASSOCIATED(sib%ww3))   NULLIFY (sib%ww3)
    IF (ASSOCIATED(sib%exo))   NULLIFY (sib%exo)
    IF (ASSOCIATED(sib%ta))   NULLIFY (sib%ta)
    IF (ASSOCIATED(sib%tc))   NULLIFY (sib%tc)
    IF (ASSOCIATED(sib%tg))   NULLIFY (sib%tg)
    IF (ASSOCIATED(sib%td1))   NULLIFY (sib%td1)
    IF (ASSOCIATED(sib%td2))   NULLIFY (sib%td2)
    IF (ASSOCIATED(sib%td3))   NULLIFY (sib%td3)
    IF (ASSOCIATED(sib%td4))   NULLIFY (sib%td4)
    IF (ASSOCIATED(sib%td5))   NULLIFY (sib%td5)
    IF (ASSOCIATED(sib%td6))   NULLIFY (sib%td6)
    IF (ASSOCIATED(sib%ra))   NULLIFY (sib%ra)
    IF (ASSOCIATED(sib%rb))   NULLIFY (sib%rb)
    IF (ASSOCIATED(sib%rc))   NULLIFY (sib%rc)
    IF (ASSOCIATED(sib%rd))   NULLIFY (sib%rd)
    IF (ASSOCIATED(sib%roff))   NULLIFY (sib%roff)
    IF (ASSOCIATED(sib%zlt))   NULLIFY (sib%zlt)
    IF (ASSOCIATED(sib%green))   NULLIFY (sib%green)
    IF (ASSOCIATED(sib%apar))   NULLIFY (sib%apar)
    IF (ASSOCIATED(sib%nee))   NULLIFY (sib%nee)
    IF (ASSOCIATED(sib%cu))   NULLIFY (sib%cu)
    IF (ASSOCIATED(sib%ct))   NULLIFY (sib%ct)
    IF (ASSOCIATED(sib%ventmf))   NULLIFY (sib%ventmf)
    IF (ASSOCIATED(sib%pco2c))   NULLIFY (sib%pco2c)
    IF (ASSOCIATED(sib%pco2i))   NULLIFY (sib%pco2i)
    IF (ASSOCIATED(sib%pco2s))   NULLIFY (sib%pco2s)
    IF (ASSOCIATED(sib%ea))   NULLIFY (sib%ea)
    IF (ASSOCIATED(sib%sha))   NULLIFY (sib%sha)
    IF (ASSOCIATED(sib%em))   NULLIFY (sib%em)
    IF (ASSOCIATED(sib%rha))   NULLIFY (sib%rha)
    IF (ASSOCIATED(sib%radvbc))   NULLIFY (sib%radvbc)
    IF (ASSOCIATED(sib%radvdc))   NULLIFY (sib%radvdc)
    IF (ASSOCIATED(sib%radnbc))   NULLIFY (sib%radnbc)
    IF (ASSOCIATED(sib%radndc))   NULLIFY (sib%radndc)
    IF (ASSOCIATED(sib%dlwbot))   NULLIFY (sib%dlwbot)
    IF (ASSOCIATED(sib%cp))   NULLIFY (sib%cp)
    IF (ASSOCIATED(sib%rho))   NULLIFY (sib%rho)
    IF (ASSOCIATED(sib%psy))   NULLIFY (sib%psy)
    IF (ASSOCIATED(sib%cupr))   NULLIFY (sib%cupr)
    IF (ASSOCIATED(sib%lspr))   NULLIFY (sib%lspr)


    RETURN
  END SUBROUTINE nullify_sib_brams

  SUBROUTINE dealloc_sib_brams(sib)

    IMPLICIT NONE
    TYPE (sib_brams) :: sib

    IF (ASSOCIATED(sib%pco2ap))  DEALLOCATE (sib%pco2ap)
    IF (ASSOCIATED(sib%pco2m))   DEALLOCATE (sib%pco2m)
    IF (ASSOCIATED(sib%rst))     DEALLOCATE (sib%rst)
    IF (ASSOCIATED(sib%capac1))  DEALLOCATE (sib%capac1)
    IF (ASSOCIATED(sib%capac2))  DEALLOCATE (sib%capac2)
    IF (ASSOCIATED(sib%snow1))   DEALLOCATE (sib%snow1)
    IF (ASSOCIATED(sib%snow2))   DEALLOCATE (sib%snow2)
    IF (ASSOCIATED(sib%fss))   DEALLOCATE (sib%fss)
    IF (ASSOCIATED(sib%fws))   DEALLOCATE (sib%fws)
    IF (ASSOCIATED(sib%assimn))   DEALLOCATE (sib%assimn)
    IF (ASSOCIATED(sib%respg))   DEALLOCATE (sib%respg)
    IF (ASSOCIATED(sib%rstfac1))   DEALLOCATE (sib%rstfac1)
    IF (ASSOCIATED(sib%rstfac2))   DEALLOCATE (sib%rstfac2)
    IF (ASSOCIATED(sib%rstfac3))   DEALLOCATE (sib%rstfac3)
    IF (ASSOCIATED(sib%rstfac4))   DEALLOCATE (sib%rstfac4)
    IF (ASSOCIATED(sib%ect))   DEALLOCATE (sib%ect)
    IF (ASSOCIATED(sib%eci))   DEALLOCATE (sib%eci)
    IF (ASSOCIATED(sib%egi))   DEALLOCATE (sib%egi)
    IF (ASSOCIATED(sib%egs))   DEALLOCATE (sib%egs)
    IF (ASSOCIATED(sib%hc))   DEALLOCATE (sib%hc)
    IF (ASSOCIATED(sib%hg))   DEALLOCATE (sib%hg)
    IF (ASSOCIATED(sib%w1))   DEALLOCATE (sib%w1)
    IF (ASSOCIATED(sib%w2))   DEALLOCATE (sib%w2)
    IF (ASSOCIATED(sib%w3))   DEALLOCATE (sib%w3)
    IF (ASSOCIATED(sib%ww1))   DEALLOCATE (sib%ww1)
    IF (ASSOCIATED(sib%ww2))   DEALLOCATE (sib%ww2)
    IF (ASSOCIATED(sib%ww3))   DEALLOCATE (sib%ww3)
    IF (ASSOCIATED(sib%exo))   DEALLOCATE (sib%exo)
    IF (ASSOCIATED(sib%ta))   DEALLOCATE (sib%ta)
    IF (ASSOCIATED(sib%tc))   DEALLOCATE (sib%tc)
    IF (ASSOCIATED(sib%tg))   DEALLOCATE (sib%tg)
    IF (ASSOCIATED(sib%td1))   DEALLOCATE (sib%td1)
    IF (ASSOCIATED(sib%td2))   DEALLOCATE (sib%td2)
    IF (ASSOCIATED(sib%td3))   DEALLOCATE (sib%td3)
    IF (ASSOCIATED(sib%td4))   DEALLOCATE (sib%td4)
    IF (ASSOCIATED(sib%td5))   DEALLOCATE (sib%td5)
    IF (ASSOCIATED(sib%td6))   DEALLOCATE (sib%td6)
    IF (ASSOCIATED(sib%ra))   DEALLOCATE (sib%ra)
    IF (ASSOCIATED(sib%rb))   DEALLOCATE (sib%rb)
    IF (ASSOCIATED(sib%rc))   DEALLOCATE (sib%rc)
    IF (ASSOCIATED(sib%rd))   DEALLOCATE (sib%rd)
    IF (ASSOCIATED(sib%roff))   DEALLOCATE (sib%roff)
    IF (ASSOCIATED(sib%zlt))   DEALLOCATE (sib%zlt)
    IF (ASSOCIATED(sib%green))   DEALLOCATE (sib%green)
    IF (ASSOCIATED(sib%apar))   DEALLOCATE (sib%apar)
    IF (ASSOCIATED(sib%nee))   DEALLOCATE (sib%nee)
    IF (ASSOCIATED(sib%cu))   DEALLOCATE (sib%cu)
    IF (ASSOCIATED(sib%ct))   DEALLOCATE (sib%ct)
    IF (ASSOCIATED(sib%ventmf))   DEALLOCATE (sib%ventmf)
    IF (ASSOCIATED(sib%pco2c))   DEALLOCATE (sib%pco2c)
    IF (ASSOCIATED(sib%pco2i))   DEALLOCATE (sib%pco2i)
    IF (ASSOCIATED(sib%pco2s))   DEALLOCATE (sib%pco2s)
    IF (ASSOCIATED(sib%ea))   DEALLOCATE (sib%ea)
    IF (ASSOCIATED(sib%sha))   DEALLOCATE (sib%sha)
    IF (ASSOCIATED(sib%em))   DEALLOCATE (sib%em)
    IF (ASSOCIATED(sib%rha))   DEALLOCATE (sib%rha)
    IF (ASSOCIATED(sib%radvbc))   DEALLOCATE (sib%radvbc)
    IF (ASSOCIATED(sib%radvdc))   DEALLOCATE (sib%radvdc)
    IF (ASSOCIATED(sib%radnbc))   DEALLOCATE (sib%radnbc)
    IF (ASSOCIATED(sib%radndc))   DEALLOCATE (sib%radndc)
    IF (ASSOCIATED(sib%dlwbot))   DEALLOCATE (sib%dlwbot)
    IF (ASSOCIATED(sib%cp))   DEALLOCATE (sib%cp)
    IF (ASSOCIATED(sib%rho))   DEALLOCATE (sib%rho)
    IF (ASSOCIATED(sib%psy))   DEALLOCATE (sib%psy)
    IF (ASSOCIATED(sib%cupr))   DEALLOCATE (sib%cupr)
    IF (ASSOCIATED(sib%lspr))   DEALLOCATE (sib%lspr)


    RETURN
  END SUBROUTINE dealloc_sib_brams

  SUBROUTINE filltab_sib_brams(sib, sibm, imean, n2, n3, ng)

    USE var_tables

    IMPLICIT NONE
    TYPE (sib_brams) :: sib, sibm
    INTEGER, INTENT(in) :: imean, n2, n3, ng
    INTEGER :: npts
    ! REAL, POINTER :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3

    IF (ASSOCIATED(sib%pco2ap))  &
         CALL vtables2 (sib%pco2ap(1,1), sibm%pco2ap(1,1) &
         , ng, npts, imean                                &
         , 'pco2ap :2:hist:anal:mpti:mpt3')
    IF (ASSOCIATED(sib%pco2m))  &
         CALL vtables2 (sib%pco2m(1,1), sibm%pco2m(1,1) &
         , ng, npts, imean                                &
         , 'pco2m :2:hist:anal:mpti:mpt3')
    IF (ASSOCIATED(sib%rst))  &
         CALL vtables2 (sib%rst(1,1), sibm%rst(1,1) &
         , ng, npts, imean                                &
         , 'rst :2:hist:anal:mpti:mpt3')
    IF (ASSOCIATED(sib%capac1))  &
         CALL vtables2 (sib%capac1(1,1), sibm%capac1(1,1) &
         , ng, npts, imean                                &
         , 'capac1 :2:hist:anal:mpti:mpt3')
    IF (ASSOCIATED(sib%capac2))  &
         CALL vtables2 (sib%capac2(1,1), sibm%capac2(1,1) &
         , ng, npts, imean                                &
         , 'capac2 :2:hist:anal:mpti:mpt3')
    IF (ASSOCIATED(sib%snow1))  &
         CALL vtables2 (sib%snow1(1,1), sibm%snow1(1,1) &
         , ng, npts, imean                                &
         , 'snow1 :2:hist:anal:mpti:mpt3')
    IF (ASSOCIATED(sib%snow2))  &
         CALL vtables2 (sib%snow2(1,1), sibm%snow2(1,1) &
         , ng, npts, imean                                &
         , 'snow2 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%fss))  &
         CALL vtables2 (sib%fss(1,1), sibm%fss(1,1) &
         , ng, npts, imean                                &
         , 'fss :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%fws))  &
         CALL vtables2 (sib%fws(1,1), sibm%fws(1,1) &
         , ng, npts, imean                                &
         , 'fws :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%assimn))  &
         CALL vtables2 (sib%assimn(1,1), sibm%assimn(1,1) &
         , ng, npts, imean                                &
         , 'assimn :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%respg))  &
         CALL vtables2 (sib%respg(1,1), sibm%respg(1,1) &
         , ng, npts, imean                                &
         , 'respg :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rstfac1))  &
         CALL vtables2 (sib%rstfac1(1,1), sibm%rstfac1(1,1) &
         , ng, npts, imean                                &
         , 'rstfac1 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rstfac2))  &
         CALL vtables2 (sib%rstfac2(1,1), sibm%rstfac2(1,1) &
         , ng, npts, imean                                &
         , 'rstfac2 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rstfac3))  &
         CALL vtables2 (sib%rstfac3(1,1), sibm%rstfac3(1,1) &
         , ng, npts, imean                                &
         , 'rstfac3 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rstfac4))  &
         CALL vtables2 (sib%rstfac4(1,1), sibm%rstfac4(1,1) &
         , ng, npts, imean                                &
         , 'rstfac4 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ect))  &
         CALL vtables2 (sib%ect(1,1), sibm%ect(1,1) &
         , ng, npts, imean                                &
         , 'ect :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%eci))  &
         CALL vtables2 (sib%eci(1,1), sibm%eci(1,1) &
         , ng, npts, imean                                &
         , 'eci :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%egi))  &
         CALL vtables2 (sib%egi(1,1), sibm%egi(1,1) &
         , ng, npts, imean                                &
         , 'egi :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%egs))  &
         CALL vtables2 (sib%egs(1,1), sibm%egs(1,1) &
         , ng, npts, imean                                &
         , 'egs :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%hc))  &
         CALL vtables2 (sib%hc(1,1), sibm%hc(1,1) &
         , ng, npts, imean                                &
         , 'hc :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%hg))  &
         CALL vtables2 (sib%hg(1,1), sibm%hg(1,1) &
         , ng, npts, imean                                &
         , 'hg :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%w1))  &
         CALL vtables2 (sib%w1(1,1), sibm%w1(1,1) &
         , ng, npts, imean                                &
         , 'w1 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%w2))  &
         CALL vtables2 (sib%w2(1,1), sibm%w2(1,1) &
         , ng, npts, imean                                &
         , 'w2 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%w3))  &
         CALL vtables2 (sib%w3(1,1), sibm%w3(1,1) &
         , ng, npts, imean                                &
         , 'w3 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ww1))  &
         CALL vtables2 (sib%ww1(1,1), sibm%ww1(1,1) &
         , ng, npts, imean                                &
         , 'ww1 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ww2))  &
         CALL vtables2 (sib%ww2(1,1), sibm%ww2(1,1) &
         , ng, npts, imean                                &
         , 'ww2 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ww3))  &
         CALL vtables2 (sib%ww3(1,1), sibm%ww3(1,1) &
         , ng, npts, imean                                &
         , 'ww3 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%exo))  &
         CALL vtables2 (sib%exo(1,1), sibm%exo(1,1) &
         , ng, npts, imean                                &
         , 'exo :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ta))  &
         CALL vtables2 (sib%ta(1,1), sibm%ta(1,1) &
         , ng, npts, imean                                &
         , 'ta :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%tc))  &
         CALL vtables2 (sib%tc(1,1), sibm%tc(1,1) &
         , ng, npts, imean                                &
         , 'tc :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%tg))  &
         CALL vtables2 (sib%tg(1,1), sibm%tg(1,1) &
         , ng, npts, imean                                &
         , 'tg :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%td1))  &
         CALL vtables2 (sib%td1(1,1), sibm%td1(1,1) &
         , ng, npts, imean                                &
         , 'td1 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%td2))  &
         CALL vtables2 (sib%td2(1,1), sibm%td2(1,1) &
         , ng, npts, imean                                &
         , 'td2 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%td3))  &
         CALL vtables2 (sib%td3(1,1), sibm%td3(1,1) &
         , ng, npts, imean                                &
         , 'td3 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%td4))  &
         CALL vtables2 (sib%td4(1,1), sibm%td4(1,1) &
         , ng, npts, imean                                &
         , 'td4 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%td5))  &
         CALL vtables2 (sib%td5(1,1), sibm%td5(1,1) &
         , ng, npts, imean                                &
         , 'td5 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%td6))  &
         CALL vtables2 (sib%td6(1,1), sibm%td6(1,1) &
         , ng, npts, imean                                &
         , 'td6 :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ra))  &
         CALL vtables2 (sib%ra(1,1), sibm%ra(1,1) &
         , ng, npts, imean                                &
         , 'ra :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rb))  &
         CALL vtables2 (sib%rb(1,1), sibm%rb(1,1) &
         , ng, npts, imean                                &
         , 'rb :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rc))  &
         CALL vtables2 (sib%rc(1,1), sibm%rc(1,1) &
         , ng, npts, imean                                &
         , 'rc :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rd))  &
         CALL vtables2 (sib%rd(1,1), sibm%rd(1,1) &
         , ng, npts, imean                                &
         , 'rd :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%roff))  &
         CALL vtables2 (sib%roff(1,1), sibm%roff(1,1) &
         , ng, npts, imean                                &
         , 'roff :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%zlt))  &
         CALL vtables2 (sib%zlt(1,1), sibm%zlt(1,1) &
         , ng, npts, imean                                &
         , 'zlt :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%green))  &
         CALL vtables2 (sib%green(1,1), sibm%green(1,1) &
         , ng, npts, imean                                &
         , 'green :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%apar))  &
         CALL vtables2 (sib%apar(1,1), sibm%apar(1,1) &
         , ng, npts, imean                                &
         , 'apar :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%nee))  &
         CALL vtables2 (sib%nee(1,1), sibm%nee(1,1) &
         , ng, npts, imean                                &
         , 'nee :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%cu))  &
         CALL vtables2 (sib%cu(1,1), sibm%cu(1,1) &
         , ng, npts, imean                                &
         , 'cu :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%ct))  &
         CALL vtables2 (sib%ct(1,1), sibm%ct(1,1) &
         , ng, npts, imean                                &
         , 'ct :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%ventmf))  &
         CALL vtables2 (sib%ventmf(1,1), sibm%ventmf(1,1) &
         , ng, npts, imean                                &
         , 'ventmf :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%pco2c))  &
         CALL vtables2 (sib%pco2c(1,1), sibm%pco2c(1,1) &
         , ng, npts, imean                                &
         , 'pco2c :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%pco2i))  &
         CALL vtables2 (sib%pco2i(1,1), sibm%pco2i(1,1) &
         , ng, npts, imean                                &
         , 'pco2i :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%pco2s))  &
         CALL vtables2 (sib%pco2s(1,1), sibm%pco2s(1,1) &
         , ng, npts, imean                                &
         , 'pco2s :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%ea))  &
         CALL vtables2 (sib%ea(1,1), sibm%ea(1,1) &
         , ng, npts, imean                                &
         , 'ea :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%sha))  &
         CALL vtables2 (sib%sha(1,1), sibm%sha(1,1) &
         , ng, npts, imean                                &
         , 'sha :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%em))  &
         CALL vtables2 (sib%em(1,1), sibm%em(1,1) &
         , ng, npts, imean                                &
         , 'em :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rha))  &
         CALL vtables2 (sib%rha(1,1), sibm%rha(1,1) &
         , ng, npts, imean                                &
         , 'rha :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%radvbc))  &
         CALL vtables2 (sib%radvbc(1,1), sibm%radvbc(1,1) &
         , ng, npts, imean                                &
         , 'radvbc :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%radvdc))  &
         CALL vtables2 (sib%radvdc(1,1), sibm%radvdc(1,1) &
         , ng, npts, imean                                &
         , 'radvdc :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%radnbc))  &
         CALL vtables2 (sib%radnbc(1,1), sibm%radnbc(1,1) &
         , ng, npts, imean                                &
         , 'radnbc :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%radndc))  &
         CALL vtables2 (sib%radndc(1,1), sibm%radndc(1,1) &
         , ng, npts, imean                                &
         , 'radndc :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%dlwbot))  &
         CALL vtables2 (sib%dlwbot(1,1), sibm%dlwbot(1,1) &
         , ng, npts, imean                                &
         , 'dlwbot :2:hist:anal:mpti:mpt3')


    IF (ASSOCIATED(sib%cp))  &
         CALL vtables2 (sib%cp(1,1), sibm%cp(1,1) &
         , ng, npts, imean                                &
         , 'cp :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%rho))  &
         CALL vtables2 (sib%rho(1,1), sibm%rho(1,1) &
         , ng, npts, imean                                &
         , 'rho :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%psy))  &
         CALL vtables2 (sib%psy(1,1), sibm%psy(1,1) &
         , ng, npts, imean                                &
         , 'psy :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%cupr))  &
         CALL vtables2 (sib%cupr(1,1), sibm%cupr(1,1) &
         , ng, npts, imean                                &
         , 'cupr :2:hist:anal:mpti:mpt3')

    IF (ASSOCIATED(sib%lspr))  &
         CALL vtables2 (sib%lspr(1,1), sibm%lspr(1,1) &
         , ng, npts, imean                                &
         , 'lspr :2:hist:anal:mpti:mpt3')


    RETURN
  END SUBROUTINE filltab_sib_brams

  SUBROUTINE zero_sib_brams(sib, n2, n3)

    IMPLICIT NONE
    TYPE (sib_brams) :: sib
    INTEGER, INTENT(in) :: n2, n3

    INTEGER :: i, j

    DO j=1,n3
       DO i=1,n2
          sib%pco2ap(i, j) = 0.
          sib%pco2m(i, j)  = 0.
          sib%rst(i, j)    = 0.
          sib%capac1(i, j) = 0.
          sib%capac2(i, j) = 0.
          sib%snow1(i, j)  = 0.
          sib%snow2(i, j)  = 0.
          sib%fss(i, j)    = 0.
          sib%fws(i, j)    = 0.
          sib%assimn(i, j) = 0.
          sib%respg(i, j)  = 0.
          sib%rstfac1(i, j) = 0.
          sib%rstfac2(i, j) = 0.
          sib%rstfac3(i, j) = 0.
          sib%rstfac4(i, j) = 0.
          sib%ect(i, j)      = 0.
          sib%eci(i, j)      = 0.
          sib%egi(i, j)      = 0.
          sib%egs(i, j)      = 0.
          sib%hc(i, j)       = 0.
          sib%hg(i, j)       = 0.
          sib%w1(i,j)          = 0.   
          sib%w2(i,j)          = 0.   
          sib%w3(i,j)         = 0.   
          sib%ww1(i,j)        = 0.   
          sib%ww2(i,j)        = 0.   
          sib%ww3(i,j)        = 0.    
          sib%exo(i,j)        = 0.   
          sib%ta(i,j)         = 0.   
          sib%tc(i,j)         = 0.   
          sib%td1(i,j)        = 0.  
          sib%td2(i,j)        = 0.   
          sib%td3(i,j)        = 0.  
          sib%td4(i,j)        = 0.  
          sib%td5(i,j)        = 0.   
          sib%td6(i,j)        = 0.  
          sib%ra(i,j)         = 0.  
          sib%rb(i,j)         = 0.  
          sib%rc(i,j)         = 0.  
          sib%rd(i,j)        = 0.   
          sib%roff(i,j)       = 0.   
          sib%zlt(i,j)        = 0.  
          sib%green(i,j)       = 0. 
          sib%apar(i,j)         = 0.
          sib%nee(i,j)          = 0. 
          sib%cu(i,j)        = 0.   
          sib%ct(i,j)        = 0.   
          sib%ventmf(i,j)       = 0. 
          sib%pco2c(i,j)        = 0. 
          sib%pco2i(i,j)        = 0. 
          sib%pco2s(i,j)        = 0.  
          sib%ea(i,j)           = 0.
          sib%sha(i,j)          = 0.
          sib%em(i,j)           = 0.
          sib%rha(i,j)          = 0.
          sib%radvbc(i,j)       = 0.
          sib%radvdc(i,j)       = 0.
          sib%radnbc(i,j)       = 0.
          sib%radndc(i,j)       = 0.
          sib%dlwbot(i,j)       = 0.
          sib%cp(i,j)           = 0.
          sib%rho(i,j)          = 0.
          sib%psy(i,j)          = 0.
          sib%cupr(i,j)         = 0.
          sib%lspr(i,j)         = 0.
          sib%roff(i,j)         = 0.
      
       ENDDO   
    ENDDO

  END SUBROUTINE zero_sib_brams

END MODULE MEM_SIB
