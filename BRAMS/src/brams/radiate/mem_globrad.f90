!==========================================================================================!
!==========================================================================================!
!     Module mem_globrad.                                                                  !
!     This module contains some parameters and flags used by CARMA radiation model.        !
!------------------------------------------------------------------------------------------!
module mem_globrad  

   use mem_aerad    , only : nsol    & ! intent(in) 
                           , nwave   & ! intent(in)
                           , ngroup  ! ! intent(in)
   use rconstants   , only : gcgs    & ! intent(in)
                           , day_sec ! ! intent(in)
   implicit none

   !---------------------------------------------------------------------------------------!
   !     Define the dimensions used by the radiation model.                                !
   !---------------------------------------------------------------------------------------!
   integer :: nvert                            ! Maximum number of layers.
   integer :: nrad                             ! Maximum number of aerosol radius bins
   integer :: nlayer                           ! Maximum number of layer boundaries
   integer :: ndbl                             ! 2 * nlayer
   integer :: nradver                          ! nrad * nvert
   integer :: nradlay                          ! nrad * nlayer
   integer, parameter :: nsolp  = 83           ! Number of solar probability intervals.
   integer, parameter :: nirp   = 71           ! Number of infra-red probability intervals.
   integer, parameter :: ntotal = nsolp + nirp ! Total number of probability intervals. 
   integer, parameter :: ngauss = 3            ! Number of Gauss quadrature points.
 
   integer, parameter :: nlow = 12500          ! Used in the Planck function calculation.
   integer, parameter :: nhigh = 35000         ! Used in the Planck function calculation.
   integer            :: ncount                ! nhigh - nlow 
 
   !---------------------------------------------------------------------------------------!
   !       Define values of flag used for reading/writing Mie coefficients.                !
   !---------------------------------------------------------------------------------------!
   logical, parameter :: l_read  = .false.
   logical, parameter :: l_write = .true.
 
   !----- Constant parameters that might be specified by an external model. ---------------!
   real, dimension(:), allocatable :: o3mix

   !---------------------------------------------------------------------------------------!
   !     Ozone mass mixing ratio (g/g) at pressure levels defined by PJ vector.  Current   !
   ! values are taken from the US standard atmosphere mid-latitude profile.                !
   !---------------------------------------------------------------------------------------!
   real   , dimension(6)                 :: o3mixp
   real                                  :: o3c
   real                                  :: vrat
   real                                  :: ptop
   real   , dimension(:)   , allocatable :: rmin
   real   , dimension(:,:) , allocatable :: r
   logical, dimension(:)   , allocatable :: is_grp_ice
   real   , dimension(:,:) , allocatable :: core_rad
   real   , dimension(:,:) , allocatable :: coreup_rad

   real   , dimension(  6)               :: pj
   real   , dimension(4,6)               :: ako3
   real   , dimension(6,6)               :: akco2
   real   , dimension(54,6)              :: akh2o
 
   !----- Time-dependent variables that might be specified by an external model. ----------!
   real                                  :: u0ext
   real                                  :: albedo_sfc
   real                                  :: emisir
   real   , dimension(:)   , allocatable :: p
   real   , dimension(:)   , allocatable :: t
   real   , dimension(:)   , allocatable :: q
 
   !---------------------------------------------------------------------------------------!
   !     Output variables, calculated by the radiation model, that might be used by the    !
   ! external model.                                                                       !
   !---------------------------------------------------------------------------------------!
   real   , dimension(:)   , allocatable :: heati
   real   , dimension(:)   , allocatable :: heats
   real   , dimension(:)   , allocatable :: heat
   real                                  :: solnet
   real                                  :: xirdown
   real                                  :: xirup
 
   !----- Variables initialised in setuprad. ----------------------------------------------!
   integer                               :: jdble
   integer                               :: jn
   real   , dimension(:)   , allocatable :: sflx
   real   , dimension(:)   , allocatable :: wvln
   real   , dimension(:)   , allocatable :: emis
   integer, dimension(:)   , allocatable :: ltemp
   real   , dimension(:)   , allocatable :: sol
   real   , dimension(:)   , allocatable :: tauray
   real   , dimension(:,:) , allocatable :: gol
   real   , dimension(:,:) , allocatable :: wol
   integer, dimension(:,:) , allocatable :: nprobi

   !----- Roundoff error precision (aka. machine epsilon). --------------------------------!
   real        , parameter :: roundoff  = epsilon(1.)
   real(kind=8), parameter :: roundoff8 = dble(roundoff)

   !----- Bergstrom's water vapour continuum fix. -----------------------------------------!
   real, dimension(nirp) :: contnm
  
   !---------------------------------------------------------------------------------------!
   !     CORERAD  - Radius of core of aerosol particles.                                   !
   !     COREREAL - Real part of refractive index of cores.                                !
   !     COREIMAG - Imaginary part of refractive index of cores.                           !
   !---------------------------------------------------------------------------------------!
   real, parameter :: corereal=2.0
   real, parameter :: coreimag=1.0
  
   !---------------------------------------------------------------------------------------!
   !     Gauss angles and Gauss weights for Gaussian integration moments (use first moment !
   ! values) n=3.                                                                          !
   !---------------------------------------------------------------------------------------!
   real, dimension(ngauss)  :: gangle
   real, dimension(ngauss)  :: gratio
   real, dimension(ngauss)  :: gweight
   !---------------------------------------------------------------------------------------!
    
   real, dimension( ntotal) :: weight

   !----- Real and imaginary parts of refractive index for liquid water. ------------------!
   real, dimension(2,nwave) :: treal
   real, dimension(2,nwave) :: ttmag

   !---------------------------------------------------------------------------------------!
   !     Nprob is the number of probability intervals in each wavelength interval.  Notice !
   ! that wave bins 11 and 12 are swapped(this is done for historical reasons).  There-    !
   ! fore, cross sections, weights, refractive indicies etc. for bins 11 and 12 must be    !
   ! swapped as well.                                                                      !
   !---------------------------------------------------------------------------------------!
   integer, dimension(ntotal) :: nprob
 
   real, allocatable :: aco2(:)
   real, allocatable :: ah2o(:)
   real, allocatable :: ao2(:)
   real, allocatable :: ao3(:)
   real, allocatable :: planck(:,:)
   real, allocatable :: xsecta(:,:)
   real, allocatable :: rup(:,:)
   real, allocatable :: qscat(:,:,:)
   integer           :: iblackbody_above
   real, allocatable :: qbrqs(:,:,:)
   real              :: t_above
   real, allocatable :: rdqext(:,:,:)

   real, allocatable :: co2(:)
   real, allocatable :: rdh2o(:)
   real, allocatable :: o2(:)
   real, allocatable :: o3(:)
   real, allocatable :: pbar(:)
   real, allocatable :: tt(:)

   real :: tgrnd
   real :: u0
   integer :: isl
   integer :: ir
   integer :: irs
   real, parameter :: cpcon =   1.006
   real, parameter :: fdegday = 1.0E-4 * gcgs * day_sec/cpcon

   !
   !      DATA (PSO2(I),I=1,77)  /  10*0.0, 4*0.75, 63*0.0      /
   !
   real :: pso2(ntotal)
                                            
   real :: pso3(ntotal)
    
   !      DATA (PSH2O(I),I=1,77) /
   !     1      14*0.0, 3*0.54, 3*0.54, 0.0, 4*0.54, 0.0, 4*0.52,
   !     2      4*0.44, 3*0.00, 4*0.62, 5*0.0, 20*0.60,   4*0.60,
   !     3      7*0.0                     /
 
   real :: psh2o(ntotal)
  
   !      DATA (PSCO2(I),I=1,77) /
   !     1       41*0.0, 4*0.82, 0.0, 20*0.88, 5*0.0, 6*0.93    /
   real :: psco2(ntotal)
   !
   !     WAVE REFERS TO THE WAVE LENGTHS OF THE CENTERS OF THE SOLAR FLUX
   !     BINS AND THE EDGES OF THE INFRARED BINS.
   !     FOR THE CURRENT MODEL SETUP, WAVE BINS 11 AND 12 ARE REVERSED
   !     IN ORDER, THAT IS 12 PRECEEDS 11. THEREFORE, CROSS SECTIONS,
   !     WEIGHTS, REFRACTIVE INDICIES ETC FOR BINS 11 AND 12 IN DATA
   !     STATEMENTS MUST BE REVERSED ALSO. THIS IS DONE FOR HISTORICAL
   !     REASONS.
   !
   !      DATA WAVE / 0.256, 0.280, 0.296,0.319,0.335,0.365,0.420,0.482,
   !     1     0.598, 0.690, 0.762, 0.719, 0.813, 0.862, 0.926, 1.005,
   !     2     1.111, 1.333, 1.562, 1.770, 2.051, 2.210, 2.584, 3.284,
   !     3     3.809, 4.292,
   !     4     4.546, 4.878, 5.128, 5.405, 5.714, 6.061, 6.452, 6.897,
   !     5     7.407, 8.333, 9.009, 10.309, 12.500, 13.889, 16.667,
   !     6     20.000, 26.316, 35.714, 62.50     /
   !
   real :: wave(nwave+1)
   !
   !     SOLAR FLUXES ( W/M**2)
   !      DATA SOLFX  /  4.1712E0, 2.5074E0, 1.2024E1, 1.7296E1, 1.2299E1,
   !     1        5.6975E1, 1.0551E2, 1.3250E2, 2.7804E2, 2.8636E1,
   !     2        5.9268E1, 5.0747E1, 5.7410E1, 4.3280E1, 7.4598E1,
   !     3        5.2732E1, 8.6900E1, 1.2126E2, 2.5731E1, 6.0107E1,
   !     4        1.8400E1, 9.5952E0, 3.5664E1, 1.2764E1, 4.0354E0,
   !     5        4.5364E0             /
   !
   !!     SOLAR FLUXES ( W/M**2)
   real :: solfx(nsol)
  
   logical :: lmie = .false.
  
   !
   !     HERE ARE THE CROSS SECTIONS FOR VARIOUS GASES. THERE ARE NWAVE
   !     OF THESE. A IS CROSS SECTION, W IS WEIGHT, PS IS PRESSURE SCALE.
   !
   !     ***********************
   !     *  DATA FOR SOLAR     *
   !     ***********************
   !
   !     UNITS ARE (CM**2)/GM
   !srf - inclusao da letra X no nome para eliminar bug
   !      DATA (XAH2O(I),I=1,77)  /      14*0.0, 0.0000, 0.1965, 9.2460,
   !
   real ::  xah2o(nsolp)
   !
   !      UNITS ARE (CM AMAGAT)
   !       DATA (XACO2(I),I=1,77)  /
   !     1       34*0.0,  0.0,    0.0035, 0.1849, 4*0.,   0.0,
  
   real ::  xaco2(nsolp)
   !
   !     UNITS ARE (CM AMAGAT)
   !      DATA (XAO2(I),I=1,77)   /
   !     1      10*0.0,  0.00, 0.00, 0.0001, 0.0022, 63*0.0  /
 
   real :: xao2(nsolp)
 
   !
   !     UNITS ARE (CM AMAGAT)
   !      DATA (XAO3(I),I=1,77)   /
   !     1      260.0, 100.9, 11.93, 0.7370, 0.0872, 0.0, 0.0,
   !     2      0.0, 0.1180, 0.0, 67*0.0   /
 
   real :: xao3(nsolp)
  
   real :: ta(ntotal)
   real :: tb(ntotal)
   real :: wa(ntotal)
   real :: wb(ntotal)
   real :: ga(ntotal)
   real :: gb(ntotal)
   real :: tia(ntotal)
   real :: tib(ntotal)
   real :: wia(ntotal)
   real :: wib(ntotal)
   real :: gia(ntotal)
   real :: gib(ntotal)
   real :: alpha(ntotal)
   real :: gama(ntotal)
   !kml2
   !real :: caseE(9,nwave+1)
   !real :: caseW(9,nwave+1)
   !real :: caseG(9,nwave+1)
   real :: caseE(9,nwave)
   real :: caseW(9,nwave)
   real :: caseG(9,nwave)

   ! DEFINED IN 'OPPROP'
   real :: wot
   real :: got
   real, allocatable :: ptempg(:)
   real, allocatable :: ptempt(:)
   !REAL, allocatable :: g0(:,:)
   !REAL, allocatable :: opd(:,:)
   !REAL, allocatable :: ptemp(:,:)
   !REAL, allocatable :: taul(:,:)

   real, allocatable :: tauh2o(:,:)
   real, allocatable :: taus(:,:)
   real, allocatable :: taua(:,:)
   real, allocatable :: g01(:,:)
   real, allocatable :: ug0(:,:)
   real, allocatable :: utaul(:,:)
   !REAL, allocatable :: w0(:,:)
   real, allocatable :: uw0(:,:)
   !REAL, allocatable :: uopd(:,:)

   ! DEFINED IN 'TWOSTR'
   real, allocatable :: u1s(:)
   real, allocatable :: u1i(:)
   real, allocatable :: acon(:,:)
   !REAL, allocatable :: ak(:,:)
   real, allocatable :: bcon(:,:)
   !REAL, allocatable :: b1(:,:)
   !REAL, allocatable :: b2(:,:)
   !REAL, allocatable :: ee1(:,:)
   !REAL, allocatable :: em1(:,:) 
   !REAL, allocatable :: em2(:,:)
   !REAL, allocatable :: el1(:,:)
   !REAL, allocatable :: el2(:,:)
   !REAL, allocatable :: gami(:,:)
   !REAL, allocatable :: af(:,:)
   !REAL, allocatable :: bf(:,:)
   !REAL, allocatable :: ef(:,:)
 
   ! DEFINED IN 'ADD'
   real, allocatable :: sfcs(:)
   !REAL, allocatable :: b3(:,:)
   !REAL, allocatable :: ck1(:,:)
   !REAL, allocatable :: ck2(:,:)
   !REAL, allocatable :: cp(:,:)
   !REAL, allocatable :: cpb(:,:)

   real, allocatable :: cm(:,:)
   !REAL, allocatable :: cmb(:,:)
   real, allocatable :: direct(:,:)
   real, allocatable :: ee3(:,:)
   real, allocatable :: el3(:,:)
   !REAL, allocatable :: fnet(:,:)
   !REAL, allocatable :: tmi(:,:)
   real, allocatable :: as(:,:)
   real, allocatable :: df(:,:)
   real, allocatable :: ds(:,:)
   real, allocatable :: xk(:,:)

   ! DEFINED IN 'NEWFLUX1'
   real, allocatable :: weit(:)
   !REAL, allocatable :: direc(:,:)
   !REAL, allocatable :: directu(:,:)
   !REAL, allocatable :: slope(:,:)
   !REAL, allocatable :: dintent(:,:,:)
   !REAL, allocatable :: uintent(:,:,:)
   !REAL, allocatable :: tmid(:,:)
   !REAL, allocatable :: tmiu(:,:)
 
   ! printed in 'radout' (defined in 'radtran')
   real :: tslu
   real :: tsld
   real :: alb_tot
   real :: tiru
   real, allocatable :: firu(:)
   real, allocatable :: firn(:)
   real, allocatable :: fslu(:)
   real, allocatable :: fsld(:)
   real, allocatable :: fsln(:)
   real, allocatable :: alb_toa(:)
   real, allocatable :: fupbs(:)
   real, allocatable :: fdownbs(:)
   real, allocatable :: fnetbs(:)
   real, allocatable :: fupbi(:)
   real, allocatable :: fdownbi(:)
   real, allocatable :: fnetbi(:)
   real, allocatable :: qrad(:,:,:)

   real :: alb_tomi
   real :: alb_toai

   character(LEN=256) :: raddatfn
   integer :: rdatfnum=22  !0
   logical :: rad_data_not_read

   contains

   subroutine initial_definitions_globrad()

     use mem_aerad, only: &
          nz_rad,         &  !INTENT(IN)
          nbin,           &  !INTENT(IN)
          nir                !INTENT(IN)

     implicit none

     ! Defining Variables:
     nvert   = nz_rad
     nrad    = nbin 
     nlayer  = nvert+1
     ndbl    = 2*nlayer
     nradver = nrad*nvert
     nradlay = nrad*nlayer
     ncount = nhigh - nlow
     
     rdatfnum = 22
     rad_data_not_read = .true.

     ! Allocating arrays:
     allocate(o3mix(nlayer))
     allocate(rmin(ngroup))
     allocate(r(nrad,ngroup))
     allocate(is_grp_ice(ngroup))
     allocate(core_rad(nrad,ngroup))
     allocate(coreup_rad(nrad,ngroup))
     allocate(p(nvert))
     allocate(t(nvert))
     allocate(q(nvert))
     allocate(heati(nlayer))
     allocate(heats(nlayer))
     allocate(heat(nlayer))
     allocate(sflx(nsol))
     allocate(wvln(nsol))
     allocate(emis(ntotal))
     !ALLOCATE(rsfx(ntotal)
     allocate(ltemp(ntotal))
     allocate(sol(ntotal))
     allocate(tauray(ntotal))
     !ALLOCATE(gcld(  ntotal,nlayer)))
     allocate(gol(ntotal,nlayer))
     allocate(wol(ntotal,nlayer))
     allocate(nprobi(nwave,2))
     allocate(aco2(ntotal))
     allocate(ah2o(ntotal))
     allocate(ao2(ntotal))
     allocate(ao3(ntotal))
     !ALLOCATE(paco2(ntotal,nlayer))
     !ALLOCATE(pah2o(ntotal,nlayer))
     !ALLOCATE(pao2(ntotal,nlayer))
     !ALLOCATE(pao3(ntotal,nlayer))
     allocate(planck(nir+1,ncount))
     !ALLOCATE(taugas(ntotal,nlayer))
     allocate(xsecta(nrad,ngroup))
     allocate(rup(nrad,ngroup))
     allocate(qscat(nrad,ngroup,nwave))
     allocate(qbrqs(nrad,ngroup,nwave))
     allocate(rdqext(nrad,ngroup,nwave))
     allocate(co2(nlayer))
     allocate(rdh2o(nlayer))
     allocate(o2(nlayer))
     allocate(o3(nlayer))
     !allocate(caer(nrad,nlayer,ngroup))
     !allocate(press(nlayer))
     allocate(pbar(nlayer))
     !allocate(dpg(nlayer))
     allocate(tt(nlayer))
     !allocate(y3(ntotal,ngauss,nlayer))
     allocate(ptempg(ntotal))
     allocate(ptempt(ntotal))
     !ALLOCATE(g0(ntotal,nlayer))
     !ALLOCATE(opd(ntotal,nlayer))
     !ALLOCATE(ptemp(ntotal,nlayer))
     !ALLOCATE(taul(ntotal,nlayer))
     allocate(tauh2o(ntotal,nlayer))
     allocate(taus(nwave,nlayer))
     allocate(taua(nwave,nlayer))
     allocate(g01(nwave,nlayer))
     allocate(ug0(ntotal,nlayer))
     allocate(utaul(ntotal,nlayer))
     !allocate(w0(ntotal,nlayer))
     allocate(uw0(ntotal,nlayer))
     !allocate(uopd(ntotal,nlayer))
     allocate(u1s( ntotal))
     allocate(u1i( ntotal))
     allocate(acon(ntotal,nlayer))
     !ALLOCATE(ak(ntotal,nlayer)))
     allocate(bcon(ntotal,nlayer))
     !ALLOCATE(b1(  ntotal,nlayer))
     !ALLOCATE(b2(  ntotal,nlayer))
     !ALLOCATE(ee1( ntotal,nlayer))
     !ALLOCATE(em1(ntotal,nlayer) )
     !ALLOCATE(em2(ntotal,nlayer))
     !ALLOCATE(el1( ntotal,nlayer))
     !ALLOCATE(el2(ntotal,nlayer))
     !ALLOCATE(gami(ntotal,nlayer))
     !ALLOCATE(af(ntotal,ndbl))
     !ALLOCATE(bf(ntotal,ndbl))
     !ALLOCATE(ef(ntotal,ndbl))
     allocate(sfcs(ntotal))
     !ALLOCATE(b3(  ntotal,nlayer))
     !ALLOCATE(ck1(   ntotal,nlayer))
     !ALLOCATE(ck2( ntotal,nlayer))
     !ALLOCATE(cp(    ntotal,nlayer))
     !ALLOCATE(cpb( ntotal,nlayer))
     allocate(cm(ntotal,nlayer))
     !allocate(cmb( ntotal,nlayer))
     allocate(direct(ntotal,nlayer))
     allocate(ee3(ntotal,nlayer))
     allocate(el3(ntotal,nlayer))
     !allocate(fnet(ntotal,nlayer))
     !allocate(tmi(ntotal,nlayer))
     allocate(as(ntotal,ndbl))
     allocate(df(ntotal,ndbl))
     allocate(ds(ntotal,ndbl))
     allocate(xk(ntotal,ndbl))
     allocate(weit(ntotal))
     !ALLOCATE(direc(ntotal,nlayer))
     !ALLOCATE(directu(ntotal,nlayer))
     !ALLOCATE(slope(ntotal,nlayer))
     !ALLOCATE(dintent(ntotal,ngauss,nlayer))
     !ALLOCATE(uintent(ntotal,ngauss,nlayer))
     !ALLOCATE(tmid(ntotal,nlayer))
     !ALLOCATE(tmiu(ntotal,nlayer))
     allocate(firu(nir))
     allocate(firn(nir))
     allocate(fslu(nsol))
     allocate(fsld(nsol))
     allocate(fsln(nsol))
     allocate(alb_toa(nsol))
     allocate(fupbs(nlayer))
     allocate(fdownbs(nlayer))
     allocate(fnetbs(nlayer))
     allocate(fupbi(nlayer))
     allocate(fdownbi(nlayer))
     allocate(fnetbi(nlayer))
     allocate(qrad(ngroup,nlayer,nrad))

   end subroutine initial_definitions_globrad

   ! ***************************************************************************

   subroutine final_definitions_globrad()

     use mem_aerad, only: &
          nz_rad              !INTENT(IN)

     implicit none

     ! Deallocating arrays:
     deallocate(o3mix)
     deallocate(rmin)
     deallocate(r)
     deallocate(is_grp_ice)
     deallocate(core_rad)
     deallocate(coreup_rad)
     deallocate(p)
     deallocate(t)
     deallocate(q)
     deallocate(heati)
     deallocate(heats)
     deallocate(heat)
     deallocate(sflx)
     deallocate(wvln)
     deallocate(emis)
     !deALLOCATE(rsfx
     deallocate(ltemp)
     deallocate(sol)
     deallocate(tauray)
     !deALLOCATE(gcld)
     deallocate(gol)
     deallocate(wol)
     deallocate(nprobi)
     deallocate(aco2)
     deallocate(ah2o)
     deallocate(ao2)
     deallocate(ao3)
     !deALLOCATE(paco2)
     !deALLOCATE(pah2o)
     !deALLOCATE(pao2)
     !deALLOCATE(pao3)
     deallocate(planck)
     !deALLOCATE(taugas)
     deallocate(xsecta)
     deallocate(rup)
     deallocate(qscat)
     deallocate(qbrqs)
     deallocate(rdqext)
     deallocate(co2)
     deallocate(rdh2o)
     deallocate(o2)
     deallocate(o3)
     !deallocate(caer)
     !deallocate(press)
     deallocate(pbar)
     !deallocate(dpg)
     deallocate(tt)
     !deallocate(y3)
     deallocate(ptempg)
     deallocate(ptempt)
     !deALLOCATE(g0)
     !deALLOCATE(opd)
     !deALLOCATE(ptemp)
     !deALLOCATE(taul)
     deallocate(tauh2o)
     deallocate(taus)
     deallocate(taua)
     deallocate(g01)
     deallocate(ug0)
     deallocate(utaul)
     !deallocate(w0)
     deallocate(uw0)
     !deallocate(uopd)
     deallocate(u1s)
     deallocate(u1i)
     deallocate(acon)
     !deALLOCATE(ak)
     deallocate(bcon)
     !deALLOCATE(b1)
     !deALLOCATE(b2)
     !deALLOCATE(ee1)
     !deALLOCATE(em1)
     !deALLOCATE(em2)
     !deALLOCATE(el1)
     !deALLOCATE(el2)
     !deALLOCATE(gami)
     !deALLOCATE(af)
     !deALLOCATE(bf)
     !deALLOCATE(ef)
     deallocate(sfcs)
     !deALLOCATE(b3)
     !deALLOCATE(ck1)
     !deALLOCATE(ck2)
     !deALLOCATE(cp)
     !deALLOCATE(cpb)
     deallocate(cm)
     !deallocate(cmb)
     deallocate(direct)
     deallocate(ee3)
     deallocate(el3)
     !deallocate(fnet)
     !deallocate(tmi)
     deallocate(as)
     deallocate(df)
     deallocate(ds)
     deallocate(xk)
     deallocate(weit)
     !deALLOCATE(direc)
     !deALLOCATE(directu)
     !deALLOCATE(slope)
     !deALLOCATE(dintent)
     !deALLOCATE(uintent)
     !deALLOCATE(tmid)
     !deALLOCATE(tmiu)
     deallocate(firu)
     deallocate(firn)
     deallocate(fslu)
     deallocate(fsld)
     deallocate(fsln)
     deallocate(alb_toa)
     deallocate(fupbs)
     deallocate(fdownbs)
     deallocate(fnetbs)
     deallocate(fupbi)
     deallocate(fdownbi)
     deallocate(fnetbi)
     deallocate(qrad)

   end subroutine final_definitions_globrad

   ! **************************************************************************

   subroutine master_read_carma_data()

     use mem_radiate, only: ISWRTYP, ILWRTYP ! Intent(in)

     implicit none

     ! Local Variables
     integer, parameter :: input_unit=22

     namelist /rad/ &
          o3mixp,pj,ako3,akco2,akh2o,contnm,gangle, &
          gratio,gweight,weight,treal,ttmag,nprob,  &
          pso2,pso3,psh2o,psco2,wave,solfx,xah2o, &
          xaco2,xao2,xao3,ta,tb,wa,wb,ga,gb,tia,tib, &
          wia,wib,gia,gib,alpha,gama,caseE,caseW,caseG

     ! Check if CARMA Radiation is selected
     if (ISWRTYP/=4 .and. ILWRTYP/=4) return
     
     
     write (unit=*,fmt='(a)') 'Reading CARMA settings...'
     open (unit=input_unit,file=raddatfn,status='old')
     read (unit=input_unit,nml=rad)
     close(unit=input_unit)
     
     rad_data_not_read=.false.

     !   DERIVED PARAMETERS
     !
     jdble   =   2*nlayer
     jn      =   jdble-1
     

   end subroutine master_read_carma_data
  
end module mem_globrad 
