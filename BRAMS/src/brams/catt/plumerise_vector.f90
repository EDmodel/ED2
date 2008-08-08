!------------------------------------------------------------------------------
! Plume rise model for vegetation fires (CPTEC/INPE 2005-2006)
! Refs.:
! Freitas, S. R., K. M. Longo, R. Chatfield, D. Latham, M. A. F. Silva Dias, M. O. Andreae, 
! E. Prins, J. C. Santos, R. Gielow and J. A. Carvalho Jr.: Including the sub-grid scale 
! plume rise of vegetation fires in low resolution atmospheric transport models. 
!  Atmospheric Chemistry and Physics and discussion, 6, 2006.
!-
! Freitas, S. R.; Longo, K. M.; M. Andreae. Impact of including the plume rise of vegetation 
! fires in numerical simulations of associated atmospheric pollutants. Geophys. Res. Lett., 
! 33, L17808, doi:10.1029/2006GL026608, 2006. 
!------------------------------------------------------------------------------


!#####################################################
!Plumerise Subroutines
!Modified by Luiz Flavio Rodrigues - 20/02/2006
!Compiling instructions:
!   including inline and vectorization (allows OpenMP)
!
!   Ex: 1) Cluster:  ifort  -static -ip -O3 -arch SSE2 -openmp -o .......
!       2) SX6    :  sxf90 -Popenmp -pi exp=esat_l,Melt,convert,glaciate,sublimate,evaporate -o .....
!       3) pgf90            -Minline
!
!      The variable "block_size" must be corretly defined:
!      If using cluster (512 Mb cache): block_size=32
!                       (1Gb    cache): block_size=64
!      If using NEC SX6               : block_size=440
!
!#######################################################

MODULE plume_utils

  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: indexi,indexj
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ijindex
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: bindex
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: block_end
  INTEGER                            :: nob,maxblock_size
  REAL                               :: prfrq             ! Defined in RAMSIN
  
  CONTAINS

    SUBROUTINE AllocIndex(block_size,ia,ja,iz,jz,IsAlloc,plume,m1,m2,m3,k_CO_smold)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: block_size,ia,iz,ja,jz,IsAlloc
      INTEGER,INTENT(IN) :: m1,m2,m3,k_CO_smold
      REAL,INTENT(IN),DIMENSION(m1,m2,m3) :: plume
      INTEGER :: i1,j1,cc,inb
      INTEGER :: resto,i,ij
      INTEGER :: iveg_ag

!!$      ! Consistency check - ALF
!!$      do j1 = ja, jz
!!$         do i1 = ia, iz
!!$            do iveg_ag= 1, 3
!!$               if (plume(K_CO_SMOLD+iveg_ag,i1,j1)>0. .and. &
!!$                    plume(iveg_ag,i1,j1)<0.) then
!!$                  print *, "PLUMERISE: Inconsistent Data: ", &
!!$                       "i, j, iveg_ag, area, fire=", &
!!$                       i1, j1, iveg_ag, plume(iveg_ag,i1,j1), &
!!$                       plume(K_CO_SMOLD+iveg_ag,i1,j1)
!!$                  print *, "PLUMERISE: Putting Zero on inconsistent data"
!!$                  call flush(6)
!!$                  ! Putting Zero on inconsistent data
!!$                  plume(iveg_ag,i1,j1)            = 0.
!!$                  plume(K_CO_SMOLD+iveg_ag,i1,j1) = 0.
!!$               endif
!!$            enddo
!!$         enddo
!!$      enddo
!!$      !
  
      IF(IsAlloc==1) THEN  
        ij=0;cc=0
        DO j1=ja,jz
          DO i1=ia,iz
            cc=cc+1
            IF(plume(k_CO_smold+1,i1,j1)+plume(k_CO_smold+2,i1,j1)+ &
                                plume(k_CO_smold+3,i1,j1) < 1.e-6 ) CYCLE
                ij=ij+1
          END DO
        END DO

        IF(block_size>=ij) THEN
          nob=1
          ALLOCATE(block_end(1))
          block_end(1)=ij
        ELSE
          nob=ij/block_size
          resto=ij-(nob*block_size)
          IF(resto==0) THEN
             ALLOCATE(block_end(nob))
             block_end=block_size
          ELSE
             nob=nob+1
             ALLOCATE(block_end(nob))
             block_end=block_size
             block_end(nob)=resto
          END IF
       END IF
        maxblock_size=maxval(block_end)

        ! Alocating only the necessary number of collumns
        ALLOCATE(indexi(maxblock_size,nob))
        ALLOCATE(indexj(maxblock_size,nob))
        ALLOCATE(ijindex(ia:iz,ja:jz))
        ALLOCATE(bindex(ia:iz,ja:jz))
        ijindex=0
        ij=0
        inb=1
        DO j1=ja,jz
          DO i1=ia,iz
            IF(plume(k_CO_smold+1,i1,j1)+plume(k_CO_smold+2,i1,j1)+ &
                                plume(k_CO_smold+3,i1,j1) < 1.e-6 ) CYCLE
                ij=ij+1
            IF(ij>block_end(inb)) THEN
              ij=1
              inb=inb+1
            END IF            
                indexi(ij,inb)=i1
                indexj(ij,inb)=j1
            ijindex(i1,j1)=ij
            bindex(i1,j1)=inb
          END DO
        END DO
      ELSE
        DEALLOCATE(indexi)
        DEALLOCATE(indexj)
        DEALLOCATE(ijindex)
        DEALLOCATE(bindex)
      END IF

    END SUBROUTINE AllocIndex
  
END MODULE plume_utils

!*****************************************************************************

!subroutine plumerise_driver(mzp,mxp,myp,ia,iz,ja,jz,k_CO_smold,k_PM25_smold,block_size)
subroutine plumerise_driver(mzp,mxp,myp,ia,iz,ja,jz,k_CO_smold,k_PM25_smold)!,block_size)
  
  use mem_basic ,only: basic_g
  use mem_grid  ,only: grid_g,ngrid,time,dtlt,itimea,zt_rams=>zt,zm_rams=>zm
  use extras    ,only: extra3d !plume array
  use node_mod  ,only: ibcon,mynum
  use mem_scalar,only: scalar_g
  USE plume_utils, ONLY: AllocIndex,nob,maxblock_size,block_end, &
       indexi,indexj,ijindex,bindex,prfrq
!  use machine_arq, only: machine ! INTENT(IN)
  
  IMPLICIT NONE

  !  INTEGER,INTENT(IN) :: block_size
  !!INTEGER ,PARAMETER :: block_size=450
  INTEGER            :: block_size
  
  INTEGER,INTENT(IN) :: mzp,mxp,myp,ia,iz,ja,jz,k_CO_smold,k_PM25_smold
  INTEGER,PARAMETER :: nkp = 200
  INTEGER,PARAMETER :: ntime = 200
  INTEGER,PARAMETER :: mdur = 53        ! duration of burn, minutes
  REAL   ,PARAMETER :: dz=100.          ! set constant grid spacing of plume grid model(meters)
  INTEGER,PARAMETER :: iveg=3
  INTEGER,PARAMETER :: im=2
  INTEGER,PARAMETER :: maxtime =(mdur+2)*60  ! model time, in seconds
  REAL   ,PARAMETER :: tdur = mdur * 60.     !- number of seconds in the burn
  !!REAL,PARAMETER :: pr_time=3600.
  INTEGER :: ijend,ibeg,iend
  REAL   ,DIMENSION(nkp)                    :: dzm,dzt,zm,zt
  REAL   ,DIMENSION(im,iveg,ntime)          :: heating
  REAL                                      :: zsurf

  REAL                                      :: passo,count

  TYPE plume_type
     !3D Real
     REAL,POINTER,DIMENSION(:,:) :: rv_b  
     REAL,POINTER,DIMENSION(:,:) :: theta_b
     REAL,POINTER,DIMENSION(:,:) :: pp_b  
     REAL,POINTER,DIMENSION(:,:) :: pi0_b  
     REAL,POINTER,DIMENSION(:,:) :: dn0_b  
     REAL,POINTER,DIMENSION(:,:) :: plume_b
     REAL,POINTER,DIMENSION(:,:) :: srcCO_b  
     REAL,POINTER,DIMENSION(:,:) :: srcPM25_b
     !2D Real
     REAL,POINTER,DIMENSION(:) :: rtgt_b  
     !2D Integer
     INTEGER,POINTER,DIMENSION(:) :: lpw_b    
  END TYPE plume_type
  
  TYPE(plume_type), ALLOCATABLE,DIMENSION(:) :: Tplume
  
  REAL, DIMENSION (:),allocatable :: XRAY,YRAY
  REAL,DIMENSION(10) :: zlev
  INTEGER :: i,j,ib,i1,j1,k,nt,tn,np,ij

  ! Checking frequency calls
  
  if(mod(time+dtlt+.001,max(tiny(dble(1.)),dble(prfrq))).le.dtlt.or.time.lt.dble(.01)) then

     ! Defining block size
    ! if (machine==1) then
        ! If using SX-6
    !    block_size = 450
    ! else
        ! Using a generic IA32
        block_size =  32
   !  endif
     !

     zsurf=0.0  

     ! define vertical grid of plume model  
     CALL Set_Grid(dz,nkp,zt,zsurf,zm,dzm,dzt)

     CALL Get_Fire_Properties(mdur,ntime,heating,iveg,im)

     CALL AllocIndex(block_size,ia,ja,iz,jz,1,extra3d(6,ngrid)%d3, &
          mzp,mxp,myp,k_CO_smold)

     ALLOCATE(Tplume(nob))
  
!$OMP PARALLEL DO  
     DO i=1,nob
        !3D Variables
        ALLOCATE(Tplume(i)%theta_b(1:maxblock_size,mzp))
        ALLOCATE(Tplume(i)%rv_b(1:maxblock_size,mzp))
        ALLOCATE(Tplume(i)%pp_b(1:maxblock_size,mzp))  
        ALLOCATE(Tplume(i)%pi0_b(1:maxblock_size,mzp))  
        ALLOCATE(Tplume(i)%dn0_b(1:maxblock_size,mzp))
        ALLOCATE(Tplume(i)%plume_b(1:maxblock_size,mzp))
        ALLOCATE(Tplume(i)%srcCO_b(1:maxblock_size,mzp))  
        ALLOCATE(Tplume(i)%srcPM25_b(1:maxblock_size,mzp))
        !2D Variables
        ALLOCATE(Tplume(i)%rtgt_b(1:maxblock_size))  
        ALLOCATE(Tplume(i)%lpw_b(1:maxblock_size))
        DO ij=1,block_end(i)
           DO k=1,mzp
              Tplume(i)%theta_b(ij,k)   = &
                   basic_g(ngrid)%theta(k,indexi(ij,i),indexj(ij,i))
              Tplume(i)%rv_b(ij,k)      = &
                   basic_g(ngrid)%rv(k,indexi(ij,i),indexj(ij,i))
              Tplume(i)%pp_b(ij,k)      = &
                   basic_g(ngrid)%pp(k,indexi(ij,i),indexj(ij,i))
              Tplume(i)%pi0_b(ij,k)     = &
                   basic_g(ngrid)%pi0(k,indexi(ij,i),indexj(ij,i))
              Tplume(i)%dn0_b(ij,k)     = &
                   basic_g(ngrid)%dn0(k,indexi(ij,i),indexj(ij,i))    
              Tplume(i)%plume_b(ij,k)   = &
                   extra3d(6,ngrid)%d3(k,indexi(ij,i),indexj(ij,i))
              Tplume(i)%srcCO_b(ij,k)   = &
                   scalar_g(1,ngrid)%srcsc(k,indexi(ij,i),indexj(ij,i))
              Tplume(i)%srcPM25_b(ij,k) = &
                   scalar_g(3,ngrid)%srcsc(k,indexi(ij,i),indexj(ij,i))    
           END DO
           Tplume(i)%rtgt_b(ij)=grid_g(ngrid)%rtgt(indexi(ij,i),indexj(ij,i))
           Tplume(i)%lpw_b(ij)=nint(grid_g(ngrid)%flpw(indexi(ij,i),indexj(ij,i)))
        END DO
     END DO
      
     !timming begin
     !!srf  CALL CPU_TIME(t1)

!$OMP PARALLEL DO  
     DO ib=1,nob
        !print*,'------------- ib=',ib,nob,maxblock_size
        CALL Plumerise(mzp,mxp,myp,ibcon,mynum  &
             ,Tplume(ib)%theta_b(1,1)   &
             ,Tplume(ib)%pp_b(1,1)      &
             ,Tplume(ib)%pi0_b(1,1)     &
             ,Tplume(ib)%dn0_b(1,1)     &
             ,Tplume(ib)%rv_b(1,1)      &
             ,Tplume(ib)%rtgt_b(1)    &
             ,Tplume(ib)%lpw_b(1)     &
             ,zt_rams(1) &
             ,zm_rams(1) &
             ,Tplume(ib)%plume_b(1,1)   & ! plume array
             ,Tplume(ib)%srcCO_b(1,1)   & ! CO source
             ,Tplume(ib)%srcPM25_b(1,1) & ! PM25 source
             ,k_CO_smold,k_PM25_smold,&
             1,block_end(ib), &
             nkp,ntime,mdur,iveg,im, &
             dzm,dzt,zm,zt, &
             heating,tdur,maxtime,dz,&
             maxblock_size,ib)
     END DO
!$OMP END PARALLEL DO  

     DO i=1,nob
        DO ij=1,block_end(i)
           DO k=1,mzp
              scalar_g(1,ngrid)%srcsc(k,indexi(ij,i),indexj(ij,i)) = &
                   Tplume(i)%srcCO_b(ij,k)
              scalar_g(3,ngrid)%srcsc(k,indexi(ij,i),indexj(ij,i)) = &
                   Tplume(i)%srcPM25_b(ij,k)  
           END DO
        END DO
     END DO
  
!$OMP PARALLEL DO  
     DO i=1,nob
        !3D Variables
        DEALLOCATE(Tplume(i)%theta_b)
        DEALLOCATE(Tplume(i)%rv_b)
        DEALLOCATE(Tplume(i)%pp_b)  
        DEALLOCATE(Tplume(i)%pi0_b)  
        DEALLOCATE(Tplume(i)%dn0_b)
        DEALLOCATE(Tplume(i)%plume_b)
        DEALLOCATE(Tplume(i)%srcCO_b)  
        DEALLOCATE(Tplume(i)%srcPM25_b)
        !2D Variables
        DEALLOCATE(Tplume(i)%rtgt_b)  
        DEALLOCATE(Tplume(i)%lpw_b)
     END DO
     DEALLOCATE(Tplume)
     DEALLOCATE(indexi,indexj)
     DEALLOCATE(ijindex)
     DEALLOCATE(bindex)
     DEALLOCATE(block_end)

  endif
     
end subroutine plumerise_driver

!-------------------------------------------------------------------------

SUBROUTINE plumerise(m1,m2,m3,ibcon,mynum   &
                    ,theta_2d &
                    ,pp_2d  &
                    ,pi0_2d &
                    ,dn0_2d &
                    ,rv_2d &
                    ,rtgt_1d &
                    ,lpw_1d &
                    ,zt_rams &
                    ,zm_rams &
                    ,plume_2d &
                    ,srcCO_2d &
                    ,srcPM25_2d &
                    ,k_CO_smold,k_PM25_smold &
                    ,ijbeg,ijend, &
                    nkp,ntime,mdur,iveg,im &
                    ,dzm,dzt,zm,zt &
                    ,heating,tdur,maxtime,dz, &
                    maxblock_size,ib)

  
  USE rconstants
  USE plume_utils, ONLY: ijindex,bindex,indexi,indexj
  IMPLICIT NONE
  
  INTEGER,PARAMETER :: moist = 10                 ! fuel moisture, %. average fuel moisture,percent dry
  REAL   ,PARAMETER :: fmoist=moist/100.         !- fuel moisture fractio
  REAL   ,PARAMETER :: fcu =1.e+6                !=> mg [gas/part] /kg [ar]
  REAL   ,PARAMETER :: bload = 10.               ! total loading, kg/m**2
  REAL   ,PARAMETER :: heat = 19.3e6             !joules/kg - floresta em alta floresta (mt)
  REAL   ,PARAMETER :: bfract = 1.                  !- combustion factor
  REAL   ,PARAMETER :: effload = bload * bfract  !- patchy burning
  
  INTEGER,INTENT(IN)                        :: m1,m2,m3
  INTEGER,INTENT(IN)                        :: ibcon,mynum,maxtime
  INTEGER,INTENT(IN)                        :: k_CO_smold,k_PM25_smold,maxblock_size
  INTEGER,INTENT(IN)                        :: ijbeg,ijend,nkp,ntime,mdur,iveg,im,ib
  
  INTEGER,INTENT(IN)   ,DIMENSION(1:maxblock_size)    :: lpw_1d
  REAL   ,INTENT(IN)   ,DIMENSION(1:maxblock_size,m1) :: theta_2d
  REAL   ,INTENT(IN)   ,DIMENSION(1:maxblock_size,m1) :: pp_2d,pi0_2d,dn0_2d,rv_2d,plume_2d
  REAL   ,INTENT(IN)   ,DIMENSION(1:maxblock_size)    :: rtgt_1d
  REAL   ,INTENT(IN)   ,DIMENSION(m1)       :: zt_rams,zm_rams
  REAL   ,INTENT(IN)   ,DIMENSION(nkp)      :: dzm,dzt,zm,zt
  REAL   ,INTENT(IN)                        :: tdur,dz

  REAL   ,INTENT(INOUT),DIMENSION(1:maxblock_size,m1) :: srcCO_2d,srcPM25_2d
  REAL   ,INTENT(INOUT),DIMENSION(im,iveg,ntime) :: heating
  
  
  INTEGER,DIMENSION(ijbeg:ijend)                  :: kmt
  INTEGER,DIMENSION(ijbeg:ijend,iveg)             :: k1,k2
  REAL   ,DIMENSION(ijbeg:ijend,im,iveg)          :: ztopmax
  REAL   ,DIMENSION(ijbeg:ijend,nkp)              :: thtcon,picon,rvcon,zcon,zzcon,qvenv,pe,te
  REAL   ,DIMENSION(ijbeg:ijend,iveg)             :: rsurf
  REAL   ,DIMENSION(nkp)                    :: w,t,qv,sc,vth,vti
  REAL   ,DIMENSION(nkp)                    :: qpas,qtotal
  REAL   ,DIMENSION(nkp)                    :: tt,qvt,qct,qht,qit,sct,visc
  REAL   ,DIMENSION(nkp)                    :: thee,rhe,sce ! environment at plume grid
  REAL   ,DIMENSION(nkp)                    :: ucon,vcon,tmpcon,dncon,prcon,scon ! environment at RAMS  grid
  REAL                                      :: advw,advt,advv,advc,advh,advi
  REAL                                      :: vhrel,virel,zbase,lbase,area
  REAL                                      :: rainbucket,dz_flam,rhodzi
  INTEGER                                   :: i,j,k,iveg_ag,imm,ixx,n
  INTEGER                                   :: ncall = 0
  
  ! Unit conversion factors
  !!fcu=1.          !=> kg [gas/part] /kg [ar]
  !!fcu =1.e+12   !=> ng [gas/part] /kg [ar]
  !----------------------------------------------------------------------
  !               index to array "plume(k,i,j)"
  ! k
  ! 1   => area media (m^2) dos focos  em biomas floresta dentro do gribox i,j
  ! 2   => area media (m^2) dos focos  em biomas savana dentro do gribox i,j
  ! 3   => area media (m^2) dos focos  em biomas pastagem dentro do gribox i,j
  ! 4   => desvio padrao da area media (m^2) dos focos : floresta
  ! 5   => desvio padrao da area media (m^2) dos focos : savana
  ! 6   => desvio padrao da area media (m^2) dos focos : pastagem
  ! 7 a 9 =>  sem uso
  !10(=k_CO_smold) => parte da emissao total de CO correspondente a fase smoldering
  !11, 12 e 13 => este array guarda a relacao entre
  !               qCO( flaming, floresta) e a quantidade total emitida
  !               na fase smoldering, isto e;
  !               qCO( flaming, floresta) =  plume(11,i,j)*plume(10,i,j)
  !               qCO( flaming, savana  ) =  plume(12,i,j)*plume(10,i,j)
  !               qCO( flaming, pastagem) =  plume(13,i,j)*plume(10,i,j)
  !20(=k_PM25_smold),21,22 e 23 o mesmo para PM25              
  !
  !24-n1 =>  sem uso
  !----------------------------------------------------------------------
  INTEGER :: ij
    
  kmt=0
  ztopmax=0.0
  thtcon=0.0;picon=0.0;rvcon=0.0;zcon=0.0;zzcon=0.0;qvenv=0.0;pe=0.0;te=0.0
  rsurf=0.0;w=0.0;t=0.0;qv=0.0
  sc=0.0;vth=0.0;vti=0.0;qpas=0.0;qtotal=0.0;tt=0.0;qvt=0.0;qct=0.0;qht=0.0
  qit=0.0;sct=0.0;visc=0.0;vcon=0.0
  advw=0.0;advt=0.0;advv=0.0;advc=0.0;advh=0.0;advi=0.0
  zbase=0.0;lbase=0.0;rainbucket=0.0
    
  DO k=1,m1
    DO ij=ijbeg,ijend
      thtcon(ij,k)=theta_2d(ij,k)          ! pot temperature
      picon(ij,k)=(pp_2d(ij,k)+pi0_2d(ij,k)) ! exner function
      rvcon(ij,k)=rv_2d(ij,k)                   ! water vapor mixing ratio
      zcon(ij,k)=zt_rams(k) *rtgt_1d(ij)   ! termod-point height
      zzcon(ij,k)=zm_rams(k) *rtgt_1d(ij)    ! W-point height
    END DO
  END DO
  
  CALL Get_Env_Condition(lpw_1d,m1-1,kmt,thtcon,picon,rvcon,zcon,qvenv, &
                         pe,te,zt,nkp,ijbeg,ijend)

  DO iveg_ag=1,iveg
    DO ij=ijbeg,ijend
      IF( plume_2d(ij,k_CO_smold + iveg_ag) < 1.e-6 ) CYCLE
      !!rsurf(ij,iveg_ag) = sqrt(max(plume_2d(ij,iveg_ag),0)/ 3.14159)
      rsurf(ij,iveg_ag) = sqrt(plume_2d(ij,iveg_ag) / 3.14159)
    END DO
  END DO
            
  CALL Makeplume(kmt,ztopmax,zm,dzm,zt,dz,maxtime,qvenv,pe,te,ijbeg,ijend,nkp,rsurf,tdur, &
                       fmoist,heating,ntime,iveg,im,plume_2d,dzt,m1,k_co_smold,ib,maxblock_size)                  
          
  DO iveg_ag=1,iveg
    DO ij=ijbeg,ijend
      !- verifica a existencia de emissao flaming para um bioma especifico
      IF( plume_2d(ij,k_co_smold + iveg_ag) < 1.e-6 ) CYCLE
      !- define o dominio vertical onde a emissao flaming ira ser colocada
      CALL Set_Flam_Vert(ztopmax(ij,:,iveg_ag),k1(ij,iveg_ag),k2(ij,iveg_ag),&
                                 m1,m2,m3,ij,zzcon,nkp,ijbeg,ijend)
    END DO
  END DO        

  !- zera o termo fonte associado `as emissoes com plumerise (k>2)      
  srcCO_2d(:,3:m1)   = 0.
  srcPM25_2d(:,3:m1) = 0.
                          
  DO iveg_ag=1,iveg
    DO k=1,nkp
      DO ij=ijbeg,ijend
        IF( plume_2d(ij,k_co_smold + iveg_ag) < 1.e-6 ) CYCLE
        IF(k<k1(ij,iveg_ag) .OR. k>k2(ij,iveg_ag)) CYCLE
        dz_flam=zzcon(ij,k2(ij,iveg_ag))-zzcon(ij,k1(ij,iveg_ag)-1)      
            rhodzi= 1./(dn0_2d(ij,k) * dz_flam)
            srcco_2d(ij,k)=   srcco_2d(ij,k) + plume_2d(ij,k_co_smold+iveg_ag)*&
                                             plume_2d(ij,k_co_smold)*&
                                            rhodzi*fcu
            srcpm25_2d(ij,k)= srcpm25_2d(ij,k) + plume_2d(ij,k_pm25_smold+iveg_ag)*&
                                             plume_2d(ij,k_pm25_smold)*&
                                            rhodzi*fcu
      END DO
    END DO
  END DO

END SUBROUTINE Plumerise

!-------------------------------------------------------------------------

SUBROUTINE Get_Env_Condition(k1,k2,kmt,thtcon,picon,rvcon,zcon,qvenv, &
                             pe,te,zt,nkp,ijbeg,ijend)

  USE rconstants
  
  IMPLICIT NONE
  INTEGER,INTENT(IN)                       :: k2,nkp,ijend,ijbeg
  INTEGER,INTENT(IN) ,DIMENSION(ijbeg:ijend)     :: k1
  REAL   ,INTENT(IN) ,DIMENSION(ijbeg:ijend,nkp) :: thtcon,picon,rvcon,zcon
  REAL   ,INTENT(IN) ,DIMENSION(nkp)       :: zt
  INTEGER,INTENT(OUT),DIMENSION(ijbeg:ijend)     :: kmt
  REAL   ,INTENT(OUT),DIMENSION(ijbeg:ijend,nkp) :: qvenv,pe,te
 
 
  INTEGER :: k,kcon,klcl,nk,nkmid,i,ij
  REAL :: themax,tlll,plll,rlll,zlll,dzdd,dzlll,tlcl,plcl,dzlcl,dummy

  REAL,DIMENSION(ijbeg:ijend,nkp) :: the,pke,thve,dne
  INTEGER :: znz(ijbeg:ijend)
 
  znz=0
  
  DO k=nkp,1,-1
    DO ij=ijbeg,ijend
      IF(zt(k)<zcon(ij,k2) .AND. znz(ij)==0) THEN
        kmt(ij)=k
        znz(ij)=1
        CYCLE
      END IF
    END DO
  END DO
 
  DO ij=ijbeg,ijend
    IF(znz(ij)==0) STOP ' envir stop 12'
  END DO
 
  DO ij=ijbeg,ijend
   nk=k2-k1(ij)+1
   CALL htint(nk,thtcon(ij,:),zcon(ij,k1(ij):),kmt(ij),the(ij,:)  ,zt)
   CALL htint(nk,rvcon (ij,:),zcon(ij,k1(ij):),kmt(ij),qvenv(ij,:),zt)
  END DO
 
  DO k=1,nkp
    DO ij=ijbeg,ijend
      IF(k>kmt(ij)) CYCLE
      qvenv(ij,k)=max(qvenv(ij,k),1e-8)
    END DO
  END DO  

  pke(:,1)=picon(:,1)

  DO k=1,nkp
    DO ij=ijbeg,ijend
      IF(k>kmt(ij)) CYCLE
      thve(ij,k)=the(ij,k)*(1.+.61*qvenv(ij,k)) ! virtual pot temperature
    END DO
  END DO

  DO k=2,nkp
    DO ij=ijbeg,ijend
      IF(k>kmt(ij)) CYCLE
      pke(ij,k)=pke(ij,k-1)-g*2.*(zt(k)-zt(k-1))  & ! exner function
                 /(thve(ij,k)+thve(ij,k-1))
    END DO
  END DO

  DO k=1,nkp
    DO ij=ijbeg,ijend
      IF(k>kmt(ij)) CYCLE
      te(ij,k)  = the(ij,k)*pke(ij,k)/cp         ! temperature (K)
      pe(ij,k)  = (pke(ij,k)/cp)**cpor*p00    ! pressure (Pa)
      dne(ij,k)= pe(ij,k)/(rgas*te(ij,k)*(1.+.61*qvenv(ij,k))) !  dry air density (kg/m3)
      !  print*,'ENV=',qvenv(k)*1000., te(k)-273.15,zt(k)
    END DO
  ENDDO

  !--------- converte press de Pa para kPa para uso modelo de plumerise
  DO k=1,nkp
    DO ij=ijbeg,ijend
      IF(k>kmt(ij)) CYCLE
      pe(ij,k) = pe(ij,k)*1.e-3
    END DO
   END DO

  ! themax=0.
  ! nkmid=20
  ! kcon=2
  ! do k=2,nkmid
  !   if(thee(k).gt.themax)then
  !        themax=thee(k)
  !        kcon=k
  !   endif
  ! enddo
  !
  ! kcon=max(2,kcon)
  ! tlll=(te(kcon)+te(kcon+1)+te(kcon-1))/3.
  ! plll=pe(kcon)*1.e+3 ! converte de volta para Pa
  ! rlll=(qvenv(kcon)+qvenv(kcon+1)+qvenv(kcon-1))/3.
  ! zlll=zt(kcon)
  !
  ! call lcl(tlll,plll,rlll,tlcl,plcl,dzlcl)
  !
  ! !  FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL
  !
  ! dzlll=1e20
  ! do k=1,kmt
  !   dzdd=abs(zt(k)-(zlll+dzlcl))
  !   if(dzdd.lt.dzlll)then
  !        dzlll=dzdd
  !        klcl=k
  !   endif
  ! enddo
  ! zbase=zt(klcl)
  !
  ! print*,' ---------- LCL --------------'
  ! print*,kmt,klcl,zt(klcl)

END SUBROUTINE Get_Env_Condition

!-------------------------------------------------------------------------

SUBROUTINE Set_Grid(dz,nkp,zt,zsurf,zm,dzm,dzt)

  IMPLICIT NONE

  INTEGER,INTENT(IN)                 :: nkp
  REAL   ,INTENT(IN)                 :: zsurf,dz
  REAL   ,INTENT(OUT),DIMENSION(nkp) :: dzm,dzt,zm,zt
  INTEGER                            :: k,mzp

  mzp=nkp
  zt(1) = zsurf
  zm(1) = zsurf
  zt(2) = zt(1) + 0.5*dz
  zm(2) = zm(1) + dz
  DO k=3,mzp
   zt(k) = zt(k-1) + dz ! thermo and water levels
   zm(k) = zm(k-1) + dz ! dynamical levels        
  ENDDO
  !print*,zsurf
  !Print*,zt(:)
  DO k = 1,mzp-1
     dzm(k) = 1. / (zt(k+1) - zt(k))
  ENDDO
  dzm(mzp)=dzm(mzp-1)

  DO k = 2,mzp
    dzt(k) = 1. / (zm(k) - zm(k-1))
  ENDDO
  dzt(1) = dzt(2) * dzt(2) / dzt(3)
  
  !   dzm(1) = 0.5/dz
  !   dzm(2:mzp) = 1./dz
END SUBROUTINE Set_Grid

!-------------------------------------------------------------------------

SUBROUTINE Set_Flam_Vert(ztopmax,k1,k2,m1,m2,m3,ij,zzcon,nkp,ijbeg,ijend)

  USE plume_utils, ONLY: ijindex

  implicit none
  INTEGER,INTENT(IN)                       :: m1,m2,m3,ij,nkp,ijbeg,ijend
  REAL   ,INTENT(IN) ,DIMENSION(ijbeg:ijend,nkp) :: zzcon
  REAL   ,INTENT(IN) ,DIMENSION(2)         :: ztopmax
  INTEGER,INTENT(OUT)                      :: k1,k2

  INTEGER                 :: imm,k
  INTEGER, DIMENSION(2)  :: k_lim

  
  DO imm=1,2
    DO k=1,nkp-1
!    do k=1,m1-1
      IF(zzcon(ij,k) > ztopmax(imm) ) EXIT
    END DO
    k_lim(imm) = k
  END DO          
  k1=max(3,k_lim(1))
  k2=max(3,k_lim(2))
  
  IF(k2 < k1) THEN
    k2=k1
  END IF
    
END SUBROUTINE Set_Flam_Vert

!-------------------------------------------------------------------------

SUBROUTINE Get_Fire_Properties(mdur,ntime,heating,iveg,im)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                           :: iveg,im,ntime
  INTEGER,INTENT(IN)                           :: mdur
  REAL   ,INTENT(OUT),DIMENSION(im,iveg,ntime) :: heating
  INTEGER                 ::  i,icount,imm,iveg_ag
  REAL                    ::   hinc ,heat_fluxw
  REAL,DIMENSION(2,3) :: heat_flux
                    
  !Heat Flux
  !IGBP Land Cover:Floresta Cerrado/woody/Savana pastagem/lavoura
  !reference:      igbp 1a4 igbp 5a9             igbp 10 a 17
  heat_flux(1,:)=(/30.0    , 4.4                , 3.3/)          !Min
  heat_flux(2,:)=(/80.0    ,23.0                , 3.3/)          !Max

  DO iveg_ag=1,iveg
    DO imm=1,im
      !fluxo de calor para o bioma
      heat_fluxw = heat_flux(imm,iveg_ag) * 1000. ! converte para W/m^2                    
      ! calculate the energy flux and water content at lboundary.
      ! fills heating() on a minute basis. could ask for a file at this po
      ! in the program. whatever is input has to be adjusted to a one
      ! minute timescale.
      !
      DO i = 1, ntime         !- make sure of energy release
        heating(imm,iveg_ag,i) = 0.0001  !- avoid possible divide by 0
      END DO  
      !                                    
      !     spread the burning evenly over the interval
      !     except for the first few minutes for stability
      icount = 1  
      !
      IF(mdur > ntime) STOP 'increase time duration (ntime) in min - see file "plumerise_mod.f90"'

      DO WHILE (icount.le.mdur)                            
        !  HEATING (ICOUNT) = HEAT * EFFLOAD / TDUR  ! W/m**2
        !  HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2
        heating (imm,iveg_ag,icount) = heat_fluxw  * 0.55     ! w/m**2 (0.55 converte para energia convectiva)
        icount = icount + 1  
      END DO  
      !     ramp for 5 minutes
      hinc = heating (imm,iveg_ag,1) / 4.  
      heating (imm,iveg_ag,1) = 0.1  
      heating (imm,iveg_ag,2) = hinc  
      heating (imm,iveg_ag,3) = 2. * hinc  
      heating (imm,iveg_ag,4) = 3. * hinc  
      !
    END DO
  END DO

END SUBROUTINE Get_Fire_Properties

!
SUBROUTINE Makeplume(kmt,ztopmax,zm,dzm,zt,dz, &
                     maxtime,qvenv,pe,te,ijbeg,ijend,nkp,rsurf,tdur, &
                     fmoist,heating,ntime,iveg,im,plume_2d,dzt,&
                     m1,k_co_smold,ib,maxblock_size)
  !
  ! *********************************************************************
  !
  !    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
  !    Alan Weinstein, JAS V.27 pp 246-255. (W),
  !    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
  !    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
  !    Originally developed by: Don Latham (USFS)
  !
  !
  ! ************************ VARIABLE ID ********************************
  !
  !        DT=COMPUTING TIME INCREMENT (SEC)
  !        DZ=VERTICAL INCREMENT (M)
  !        LBASE=LEVEL ,CLOUD BASE
  !
  !        CONSTANTS:
  !          G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
  !          R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
  !          CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
  !          HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
  !          HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
  !          HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
  !          EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
  !          DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
  !          TFREEZE = FREEZING TEMPERATURE (K)
  !
  !
  !        PARCEL VALUES:
  !          TEMP = TEMPERATURE (K)
  !          TXS = TEMPERATURE EXCESS (K)
  !          QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
  !          QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
  !          QC = WATER CONTENT (G/G DRY AIR)
  !          QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
  !          QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
  !          RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
  !          ES = SATURATION VAPOR PRESSURE (kPa)
  !
  !        ENVIRONMENT VALUES:
  !          TE = TEMPERATURE (K)
  !          PE = PRESSURE (kPa)
  !          QVENV = WATER VAPOR (G/G)
  !          RHE = RELATIVE HUMIDITY FRACTION (e/esat)
  !          DNE = dry air density (kg/m^3)
  !
  !        HEAT VALUES:
  !          HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
  !          MDUR = DURATION OF BURN, MINUTES
  !
  !          VVEL = VERTICAL VELOCITY (M/S)
  !          RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
  !          RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
  !          ALPHA = ENTRAINMENT CONSTANT
  !          MAXTIME = TERMINATION TIME (MIN)
  !
  !
  !**********************************************************************
  USE plume_utils, ONLY: ijindex,indexj,indexi
  
  IMPLICIT NONE
  ! ******************* SOME CONSTANTS **********************************
  !
  !         XNO=10.0E06 median volume diameter raindrop (K table 4)
  !         VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)
  !
  REAL   ,PARAMETER :: vc=5.107387
  REAL   ,PARAMETER :: g=9.80796,r=287.04,cp=1004.,eps=0.622,tmelt =273.3
  REAL   ,PARAMETER :: heatsubl=2.834e6,heatfus=3.34e5,heatcond=2.501e6
  REAL   ,PARAMETER :: tfreeze=269.3
  REAL   ,PARAMETER :: alpha = 0.1     !- entrainment constant
  REAL   ,PARAMETER :: tstpf = 2.0     !- timestep factor
  REAL   ,PARAMETER :: viscosity = 500.!- viscosity constant (original value: 0.001)
  INTEGER,PARAMETER :: deltaK  = 20
  INTEGER,PARAMETER :: nrectotal=150  

  INTEGER,INTENT(IN)                                 :: ntime,iveg,im
  INTEGER,INTENT(IN)                                 :: m1,k_co_smold,ib,maxblock_size
  INTEGER,INTENT(IN)                                 :: ijbeg,ijend,nkp,maxtime
  REAL   ,INTENT(IN),DIMENSION(im,iveg,ntime)    :: heating
  REAL   ,INTENT(IN),DIMENSION(1:maxblock_size,m1)   :: plume_2d
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)  :: qvenv,pe,te
  REAL   ,INTENT(IN),DIMENSION(nkp)                 :: zm,dzm,zt,dzt
  REAL   ,INTENT(IN)                                 :: dz,tdur,fmoist
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,iveg) :: rsurf
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend)      ::  kmt

  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,im,iveg) :: ztopmax

  REAL      ,DIMENSION(ijbeg:ijend,ntime) :: ztop_
  REAL      ,DIMENSION(ijbeg:ijend,nkp)   :: tt,qvt,qct,qit,sct,wc
  REAL      ,DIMENSION(ijbeg:ijend,nkp)   :: qht,wt,rho,visc,cvi
  REAL      ,DIMENSION(ijbeg:ijend,nkp)          :: vth,vti,qsat,txs,qv,radius,est
  REAL      ,DIMENSION(ijbeg:ijend,nkp)          :: qh,qi,qc,vvel,temp
  REAL      ,DIMENSION(ijbeg:ijend)       :: wbar,dqsdz
  REAL      ,DIMENSION(ijbeg:ijend)       :: time,dt,wmax
  REAL      ,DIMENSION(ijbeg:ijend)       :: ztop,dt_save

  INTEGER   ,DIMENSION(ijbeg:ijend)          :: nm1,mintime,kkmax,kk
  LOGICAL   ,DIMENSION(ijbeg:ijend)          :: isdone
  LOGICAL :: isready
  
  INTEGER                   :: k,l,kci  !,ilastprint,izprint
  INTEGER                   :: i_micro,n_sub_step,ij,iveg_ag,imm
  REAL                     :: rmaxtime,es,adiabat
  
  REAL,EXTERNAL :: esat_l

  !CALL ftrace_region_begin('plume')
  
  rmaxtime = float(maxtime)
  n_sub_step=3
    
  DO iveg_ag=1,iveg
    DO imm=1,im  
    
      IF(iveg_ag == 3 .and. imm == 2) CYCLE
    
      DO k=1,nkp
        DO ij=ijbeg,ijend
          cvi=0.0
          qv=0.0
        END DO
      END DO
      
      wbar=0.0
      dqsdz=0.0
      ztop=0.0
      time=0.0
      !*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
      mintime = 1  
      ztopmax(:,imm,iveg_ag) = 0.
      ztop    = 0.
        time = 0.  
          dt = 1.
       wmax = 1.
      kkmax = 10
      kci=k_co_smold + iveg_ag
      l = 1 ! L initialization
      isdone=.false.  
      DO ij=ijbeg,ijend
         IF(plume_2d(ij,kci) < 1.e-6) THEN
           isdone(ij)=.true.
           CYCLE
         END IF
      END DO
      isready=.false.              

      !CALL ftrace_region_begin('initial')
      !initialization
      CALL Initial(kmt,imm,iveg_ag,iveg,ijbeg,ijend,nkp,rsurf,qvenv,pe,te,txs,vvel, &
             temp,wc,wt,qv,vth,vti,qh,qi,qc,est,qsat,rho,radius,visc,&
             viscosity,alpha,zt,dt,mintime,time,tdur,fmoist,heating,&
             ntime,l,wbar,dqsdz,cvi,plume_2d,m1,im,isdone,ib,maxblock_size)
      !CALL ftrace_region_end('initial')
              
     ! ******************* model evolution ******************************
     DO WHILE(.not. Isready)!beginning of time loop
      
       !CALL ftrace_region_begin('begin_loop')  
        DO ij=ijbeg,ijend
          IF(isdone(ij)) CYCLE
  
          IF(time(ij)>rmaxtime) THEN
            isdone(ij)=.true.
            CYCLE
          END IF
          !-- set model top integration
          nm1(ij) = min(kmt(ij), kkmax(ij) + deltak)
          !-- set timestep
          dt(ij) = min(5.,(zm(2)-zm(1)) / (tstpf * wmax(ij)))
          !-- elapsed time, sec
          time(ij) = time(ij)+dt(ij)
          !-- elapsed time, minutes                                        
          mintime(ij) = 1 + int (time(ij)) / 60        
          wmax(ij) = 1.  !no zeroes allowed.
  
        END DO
       !CALL ftrace_region_end('begin_loop')  

       !CALL ftrace_region_begin('Tend0_Plumerise')  
        !-- zerout all model tendencies
        CALL Tend0_Plumerise(nm1,wt,tt,qvt,qct,qht,qit,sct,nkp,ijbeg,ijend,isdone)
       !CALL ftrace_region_end('Tend0_Plumerise')


        l=1
        !CALL ftrace_region_begin('Lbound')
        CALL Lbound(imm,iveg_ag,qh,qi,qc,rsurf,plume_2d,iveg,ijbeg,ijend,alpha,&
                qvenv,pe,te,wc,vvel,temp,vth,vti,txs,visc,rho,qv,viscosity,dt, &
                mintime,time,tdur,fmoist,heating,ntime,l,wbar,dqsdz,cvi,est, &
                qsat,nkp,im,isdone,ib,maxblock_size,m1)
        !CALL ftrace_region_end('Lbound')

        !-- dynamics for the level k>1
        !-- W advection
        !   call vel_advectc_plumerise(NM1,wc(ij,:),wt(ij,:),DNE,DZM)
        !CALL ftrace_region_begin('Vel_Advectc_Plumerise')
        CALL Vel_Advectc_Plumerise(nm1,wc,wt,rho,dzm,nkp,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Vel_Advectc_Plumerise')

        !-- scalars advection 1
        !CALL ftrace_region_begin('Scl_Advectc_Plumerise')
        CALL Scl_Advectc_Plumerise('SC',nm1,nkp,dt,wc,vvel,rho,dzm,dzt,zm,zt, &
                                   temp,tt,qv,qvt,qc,qct,qh,qht,qi,qit, &
                                   ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Scl_Advectc_Plumerise')
  
        !-- scalars advection 2
        !call scl_advectc_plumerise2('SC',NM1)
  
        !-- scalars entrainment, adiabatic
        !-- scalars entrainment, adiabatic
        !CALL ftrace_region_begin('Scl_Misc')  
        CALL Scl_Misc(nm1,nkp,ijbeg,ijend,qvenv,te,vvel,temp,qv,qc,qh,qi, &
                      tt,qvt,qct,qht,qit,radius,alpha,adiabat,wbar,isdone)
        !CALL ftrace_region_end('Scl_Misc')  

        
        !-- gravity wave damping using Rayleigh friction layer fot T
        !CALL ftrace_region_begin('Damp_Grav_Wave1')  
        CALL Damp_Grav_Wave(1,nm1,deltak,dt,zt,zm,vvel,temp,tt,qv,qh,qi,&
                                qc,te,pe,qvenv,nkp,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Damp_Grav_Wave1')  
        
        !CALL ftrace_region_begin('microphysics')  
        !-- microphysics
        DO ij=ijbeg,ijend
          IF(isdone(ij)) CYCLE
          dt_save(ij)=dt(ij)
          dt(ij)=dt(ij)/float(n_sub_step)
        END DO
        
        
        DO i_micro=1,n_sub_step
          !-- sedim ?
        !CALL ftrace_region_begin('micro_Fallpart')  
          CALL Fallpart(nm1,nkp,qvt,qct,qht,qit,rho,vvel,qh,qi,zm, &
                          vth,vti,cvi,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('micro_Fallpart')  
          !-- microphysics
        !CALL ftrace_region_begin('micro_Fase1')  
          DO l=2,nkp
            DO ij=ijbeg,ijend                                                  
              IF(isdone(ij) .OR. l>nm1(ij)-1) CYCLE
                wbar(ij)    = 0.5*(vvel(ij,l)+vvel(ij,l-1))
                es          = 0.1*esat_l(temp(ij,l))               !blob saturation vapor pressure, em kpa
                qsat(ij,l) = (eps * es) / (pe(ij,l) - es)  !blob saturation lwc(ij,:) g/g dry air
                est (ij,l) = es  
                rho (ij,l) = 3483.8 * pe(ij,l) / temp(ij,l) ! air parcel density , g/m**3
              !srf18jun2005
              !          IF (vvel(L) .ge. 0.) DQSDZ = (QSAT(L  ) - QSAT(L-1)) / (ZT(L  ) -ZT(L-1))
              !          IF (vvel(L) .lt. 0.) DQSDZ = (QSAT(L+1) - QSAT(L  )) / (ZT(L+1) -ZT(L  ))
            END DO
        !CALL ftrace_region_end('micro_Fase1')  
        !CALL ftrace_region_begin('micro_Fase2')  
            DO ij=ijbeg,ijend                                                  
              IF(isdone(ij) .OR. l>nm1(ij)-1) CYCLE
                IF (vvel(ij,l) >= 0.) then
                  dqsdz(ij) = (qsat(ij,l+1) - qsat(ij,l-1)) / (zt(l+1 )-zt(l-1))
                ELSE
                  dqsdz(ij) = (qsat(ij,l+1) - qsat(ij,l-1)) / (zt(l+1) -zt(l-1))
                END IF
              END DO
            !CALL ftrace_region_end('micro_Fase2')  
        
            !CALL ftrace_region_begin('micro_Waterbal')                
            CALL Waterbal(nkp,l,qc,qh,qi,qv,wbar,dqsdz,dt,temp,rho,est,cvi,qsat,ijbeg,ijend,isdone,nm1)  
            !CALL ftrace_region_end('micro_Waterbal')  
          END DO
        END DO
  
        DO ij=ijbeg,ijend
          IF(isdone(ij)) CYCLE
          dt(ij)=dt_save(ij)
        END DO
        !CALL ftrace_region_end('microphysics')  

        !-- W-viscosity for stability
        !CALL ftrace_region_begin('Visc_W')  
        CALL Visc_W(nm1,deltak,kmt,nkp,vvel,temp,qv,qh,qc,qi,zm,zt,visc, &
                     wt,tt,qvt,qct,qht,qit,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Visc_W')  

        !-- update scalars
        !CALL ftrace_region_begin('Update_Plumerise1')  
        CALL Update_Plumerise(nm1,'S',nkp,dt,wt,tt,qvt,qct,qht,qit,&
                            vvel,temp,qv,qh,qc,qi,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Update_Plumerise1')  
      

        !CALL ftrace_region_begin('Hadvance_Plumerise1')  
        CALL Hadvance_Plumerise(1,nm1,dt,wc,wt,vvel,mintime,nkp,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Hadvance_Plumerise1')  
          
        !-- Buoyancy
        !CALL ftrace_region_begin('Buoyancy_Plumerise')  
        CALL Buoyancy_Plumerise(nm1,temp,te,qv,qvenv, &
                                qh,qi,qc,wt,nkp,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Buoyancy_Plumerise')  
  
        
        !-- Entrainment
        !CALL ftrace_region_begin('Entrainment')  
        CALL Entrainment(nm1,vvel,wt,radius,alpha,nkp,wbar,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Entrainment')  


        !-- update vvel
        !CALL ftrace_region_begin('Update_Plumerise2')  
        CALL Update_Plumerise(nm1,'W',nkp,dt,wt,tt,qvt,qct,qht,qit,&
                            vvel,temp,qv,qh,qc,qi,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Update_Plumerise2')  
    
        !CALL ftrace_region_begin('Hadvance_Plumerise2')  
        CALL Hadvance_Plumerise(2,nm1,dt,wc,wt,vvel,mintime,nkp,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Hadvance_Plumerise2')  
          
        !-- misc
        !CALL ftrace_region_begin('misc')  
        DO k=2,nkp
          DO ij=ijbeg,ijend
            IF(isdone(ij) .OR. k>nm1(ij)) CYCLE
            !    pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
            es        = 0.1*esat_l(temp(ij,k)) !blob saturation vapor pressure, em kPa
            !    rotina do plumegen calcula em kPa
            !    es        = esat_pr (t(k))  !blob saturation vapor pressure, em kPa
            qsat(ij,k) = (eps * es) / (pe(ij,k) - es)  !blob saturation lwc(ij,:) g/g dry air
            est (ij,k) = es  
            txs (ij,k) = temp(ij,k) - te(ij,k)
            rho (ij,k) = 3483.8 * pe(ij,k) / temp (ij,k) ! air parcel density , g/m**3
                                                                ! no pressure diff with radius
            IF((abs(wc(ij,k)))>wmax(ij)) wmax(ij) = abs(wc(ij,k)) ! keep wmax largest w
          END DO
        END DO  
        !CALL ftrace_region_end('misc')  
    
        ! Gravity wave damping using Rayleigh friction layer for vvel
        !CALL ftrace_region_begin('Damp_Grav_Wave2')  
        CALL Damp_Grav_Wave(2,nm1,deltak,dt,zt,zm,vvel,temp,tt,qv,qh,qi,qc,&      
                             te,pe,qvenv,nkp,ijbeg,ijend,isdone)
        !CALL ftrace_region_end('Damp_Grav_Wave2')  
        
        !CALL ftrace_region_begin('findtop')  
        DO ij=ijbeg,ijend
          IF(isdone(ij)) CYCLE
          !-- try to find the plume top (above surface height)
          kk(ij) = 1
          DO WHILE (vvel(ij,kk(ij)) > 1.)  
            kk(ij) = kk(ij) + 1  
            ztop(ij) =  zm(kk(ij))
          END DO  
        END DO
        
        DO ij=ijbeg,ijend
          IF(isdone(ij)) CYCLE
            ztop_(ij,mintime(ij)) = ztop(ij)
            ztopmax(ij,imm,iveg_ag) = max(ztop(ij), ztopmax(ij,imm,iveg_ag))
            kkmax(ij)   = max (kk(ij), kkmax(ij))
        END DO
        !CALL ftrace_region_end('findtop')  
      
        !
        !srf-27082005
        ! if the solution is going to a stationary phase, exit
        !CALL ftrace_region_begin('endloop')  
        DO ij=ijbeg,ijend
          IF(isdone(ij)) CYCLE
          IF(mintime(ij) > 10) THEN
            IF( abs(ztop_(ij,mintime(ij))-ztop_(ij,mintime(ij)-10)) < dz ) isdone(ij)=.true.
          END IF
        END DO

        !Criterio de parada: Checa se todas as colunas foram calculadas
        isready=.true.
        DO ij=ijbeg,ijend
          IF( plume_2d(ij,kci) < 1.e-6) CYCLE
          IF(.not. isdone(ij)) isready=.false.
        END DO
        !CALL ftrace_region_end('endloop')  

      END DO !next timestep
    END DO !Next imm
  END DO   !Next iveg_ag

  DO ij=ijbeg,ijend
    IF( plume_2d(ij,k_co_smold + iveg_ag) < 1.e-6 ) CYCLE
    !-- so ha um valor para eflux de pastagem=> ztopmax(2)=ztopmax(1)
    ztopmax(ij,2,3)=ztopmax(ij,1,3)
  END DO
  
  !CALL ftrace_region_end('plume')

END SUBROUTINE Makeplume

!-------------------------------------------------------------------------------

SUBROUTINE Initial(kmt,imm,iveg_ag,iveg,ijbeg,ijend,nkp,rsurf,qvenv,pe,te,txs,vvel, &
                       temp,wc,wt,qv,vth,vti,qh,qi,qc,est,qsat,rho,radius,visc,&
                       viscosity,alpha,zt,dt,mintime,time,tdur,fmoist,heating, &
                     ntime,l,wbar,dqsdz,cvi,plume_2d,m1,im,isdone,ib,maxblock_size)  
        
  USE plume_utils, ONLY: ijindex
  
  IMPLICIT NONE
  REAL, PARAMETER                           :: tfreeze = 269.3

  INTEGER,INTENT(IN)                                 :: imm,iveg_ag,m1,ib,maxblock_size
  INTEGER,INTENT(IN)                                 :: iveg,ijbeg,ijend,nkp
  INTEGER,DIMENSION(ijbeg:ijend),INTENT(IN)    :: kmt,mintime
  INTEGER,INTENT(IN)                                 :: ntime,l,im
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,iveg)     :: rsurf
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)      :: qvenv,pe,te
  REAL   ,INTENT(IN),DIMENSION(nkp)               :: zt
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)::   cvi
  REAL   ,INTENT(IN)                                 :: viscosity,alpha,tdur
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend)    :: time,dt
  REAL   ,INTENT(IN)                                 :: fmoist
  REAL   ,DIMENSION(ijbeg:ijend),INTENT(IN)    :: wbar,dqsdz
  LOGICAL,DIMENSION(ijbeg:ijend),INTENT(IN)    :: isdone
  REAL   ,INTENT(IN) ,DIMENSION(1:maxblock_size,m1)      :: plume_2d
  REAL   ,INTENT(IN) ,DIMENSION(im,iveg,ntime) :: heating

  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: wc,wt,vth,vti,rho,visc,qsat,txs,qv
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: radius,est,qh,qi,qc,vvel,temp
    
  INTEGER :: n,ij
  INTEGER :: isub,k,n1,n2,n3,lbuoy,itmp,isubm1
  REAL    :: xn1,xi
  REAL   ,DIMENSION(ijbeg:ijend) :: es
  
  REAL,EXTERNAL :: Esat_l
  REAL,EXTERNAL :: Esat_Pr
  !
    
  !n=kmt(ij)
  !ij=ijindex(i1,j1)
  
  DO k = 1, nkp                            
    DO ij=ijbeg,ijend
       IF(isdone(ij)) CYCLE
       IF(K>kmt(ij)) CYCLE
       txs (ij,k) = 0.0  
       vvel (ij,k) = 0.0              
       temp(ij,k) = te(ij,k)       !blob set to environment                    
       wc(ij,k) = 0.0
       wt(ij,k) = 0.0
       qv(ij,k) = qvenv(ij,k)   !blob set to environment                
       vth(ij,k) = 0.          !initial rain velocity = 0                              
       vti(ij,k) = 0.          !initial ice  velocity = 0                              
       qh(ij,k) = 0.          !no rain                              
       qi(ij,k) = 0.          !no ice                              
       qc(ij,k) = 0.          !no cloud drops                      
    END DO
  END DO  

  DO k = 1, nkp                            
    DO ij=ijbeg,ijend
       IF(isdone(ij)) CYCLE
       IF(K>kmt(ij)) CYCLE
       !  pe esta em kpa  - esat do rams esta em mbar = 100 pa = 0.1 kpa
       !  rotina do plumegen calcula em kpa
       es(ij)       = 0.1*Esat_l(temp(ij,k)) !blob saturation vapor pressure, em kpa
    END DO
  END DO  

  ! initialize temperature structure,to the end of equal spaced sounding,
  DO k = 1, nkp                            
    DO ij=ijbeg,ijend
       IF(isdone(ij)) CYCLE
       IF(K>kmt(ij)) CYCLE
       !  es       = esat_pr (t(k))  !blob saturation vapor pressure, em kpa
       est  (ij,k) = es(ij)  
       qsat (ij,k) = (.622 * es(ij)) / (pe(ij,k) - es(ij)) !saturation lwc g/g
       rho  (ij,k) = 3483.8 * pe(ij,k) / temp(ij,k)         !dry air density g/m**3    
    END DO
  END DO  

  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
      ! Initialize the entrainment radius, Turner-style plume
      radius(ij,1) = rsurf(ij,iveg_ag)
      !  Initialize the viscosity
      visc(ij,1) = viscosity
  END DO
  
  DO k=2,nkp
    DO ij=ijbeg,ijend
       IF(isdone(ij)) CYCLE
       IF(K>kmt(ij)) CYCLE
       radius(ij,k) = radius(ij,k-1)+(6./5.)*alpha*(zt(k)-zt(k-1))
    END DO
  END DO
    
   DO k=2,nkp
    DO ij=ijbeg,ijend
       IF(isdone(ij)) CYCLE
       IF(K>kmt(ij)) CYCLE
       visc(ij,k) = viscosity!max(1.e-3,visc(k-1) - 1.* VISCOSITY/float(nkp))
     END DO
   ENDDO
   !--   Initialize gas/concentration
  !DO k =10,20
  !   SC(k) = 20.
  !ENDDO
  !stop 333

  CALL Lbound(imm,iveg_ag,qh,qi,qc,rsurf,plume_2d,iveg,ijbeg,ijend,alpha,&
             qvenv,pe,te,wc,vvel,temp,vth,vti,txs,visc,rho,qv,viscosity,dt, &
             mintime,time,tdur,fmoist,heating,ntime,l,wbar,dqsdz,cvi,est, &
             qsat,nkp,im,isdone,ib,maxblock_size,m1)
                
END SUBROUTINE Initial

!-------------------------------------------------------------------------------

SUBROUTINE Lbound(imm,iveg_ag,qh,qi,qc,rsurf,plume_2d,iveg,ijbeg,ijend,alpha,&
                  qvenv,pe,te,wc,vvel,temp,vth,vti,txs,visc,rho,qv,viscosity,dt, &
                  mintime,time,tdur,fmoist,heating,ntime,l,wbar,dqsdz,cvi, &
                  est,qsat,nkp,im,isdone,ib,maxblock_size,m1)  
  !
  ! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
  !
  ! source of equations: J.S. Turner Buoyancy Effects in Fluids
  !                         Cambridge U.P. 1973 p.172,
  !                         G.A. Briggs Plume Rise, USAtomic Energy Commissio
  !                         TID-25075, 1969, P.28
  !
  ! fundamentally a point source below ground. at surface, this produces
  ! a velocity w(1) and temperature T(1) which vary with time. There is
  ! also a water load which will first saturate, then remainder go into
  ! QC(1).
  ! EFLUX = energy flux at ground,watt/m**2 for the last DT
  !
  
  USE plume_utils, ONLY: indexi,indexj
  
  implicit none
  REAL, PARAMETER :: g=9.80796,r=287.04,cp=1004.6,eps=0.622,tmelt=273.3
  REAL, PARAMETER :: tfreeze=269.3,pi=3.14159,e1=1./3.,e2=5./3.
  REAL, PARAMETER :: heat = 19.3e6    !joules/kg - floresta em alta floresta (mt)

  INTEGER,INTENT(IN) :: imm,iveg_ag,ijbeg,ijend,iveg
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend) :: mintime
  INTEGER,INTENT(IN) :: ntime,l,nkp,im,ib,maxblock_size,m1

  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp)           :: cvi
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,iveg)    :: rsurf
  REAL   ,INTENT(IN)   ,DIMENSION(1:maxblock_size,m1)    :: plume_2d
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp)     :: qvenv,pe,te
  REAL   ,INTENT(IN)   ,DIMENSION(im,iveg,ntime) :: heating
  REAL   ,INTENT(IN)                             :: alpha,viscosity
  REAL   ,INTENT(IN)                             :: fmoist,tdur
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend)   :: wbar,dqsdz,time,dt
  LOGICAL,INTENT(IN)   ,DIMENSION(ijbeg:ijend)   :: isdone  
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)  :: est,qh,qi,qc,vvel,temp
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)  :: visc,rho,qsat,wc,vth,vti,txs,qv

  REAL,EXTERNAL :: Esat_l
  REAL,EXTERNAL :: Esat_Pr
  
  REAL,DIMENSION(ijbeg:ijend)   :: es
  INTEGER,DIMENSION(ijbeg:ijend)   :: nm1_l
  
  REAL    :: eflux, water,  pres, c1,  c2, f, zv,  denscor, xwater
  INTEGER :: ij
  !
  nm1_l=1
  !
  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    qh(ij,1) = qh(ij,2)   !soak up hydrometeors
    qi(ij,1) = qi(ij,2)              
    qc(ij,1) = 0.            !no cloud here
    !
    !CALL Burn(eflux, water,imm,iveg_ag)  
    IF (time(ij)>tdur) THEN    !is the burn over?  
       eflux = 0.000001         !prevent a potential divide by zero
       water = 0.  
    ELSE  
    !                                                      
       eflux = heating(imm,iveg_ag,mintime(ij))                            ! watts/m**2                                                  
    !  water = eflux * (dt / heat) * (0.5 + fmoist)        ! kg/m**2
      water = eflux * (dt(ij)/heat) * (0.5 + fmoist) /0.55 ! kg/m**2
      water = water * 1000.                                ! g/m**2
    ENDIF  
    !
    !  calculate parameters at boundary from a virtual buoyancy point source
    !
    !   IF(ij==1 .AND. time<=65.0 .AND. time>=45.0)  WRITE (81,FMT='(3(I3.3),1X,5(F18.8,1X))') iveg_ag,imm,mintime,time,tdur,eflux,&
    !                                                                                  water,heating(imm,iveg_ag,mintime)
    pres = pe(ij,1) * 1000.   !need pressure in n/m**2
    c1 = 5. / (6. * alpha)  !alpha is entrainment constant
    c2 = 0.9 * alpha  
    f = eflux / (pres * cp * pi)  
    f = g * r * f * plume_2d(ij,iveg_ag)  !buoyancy flux
    zv = c1 * rsurf(ij,iveg_ag)  !virtual boundary height

    !WRITE(99,FMT='(4(I3.3,1X),6(E18.8,1X))') ij,indexi(ij,ib),indexj(ij,ib),iveg_ag,c1,c2,f,e1,zv,rsurf(ij,iveg_ag)
    vvel(ij,1) = c1 * ((c2 * f) **e1) / zv**e1  !boundary velocity
    denscor = c1 * f / g / (c2 * f) **e1 / zv**e2   !density correction
    temp(ij,1) = te(ij,1) / (1. - denscor)    !temperature of virtual plume at zsurf
    !
    wc(ij,1) = vvel(ij,1)
  
    !SC(1) = SCE(1)+F/1000.*dt  ! gas/particle (g/g)
  
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        match dw/dz,dt/dz at the boundary. F is conserved.
    !
    !wbar = w (1) * (1. - 1. / (6. * zv) )  
    !advw = wbar * w (1) / (3. * zv)  
    !advt = wbar * (5. / (3. * zv) ) * (denscor / (1. - denscor) )  
    !advc = 0.  
    !advh = 0.  
    !advi = 0.  
    !adiabat = - wbar * g / cp  
    vth(ij,1) = - 4.  
    vti(ij,1) = - 3.  
    txs(ij,1) = temp(ij,1) - te(ij,1)  
    visc(ij,1) = viscosity  
    rho(ij,1) = 3483.8 * pe(ij,1) / temp(ij,1)        !air density at level 1, g/m**3
    xwater = water/ (vvel(ij,1) * dt(ij) * rho(ij,1) )   !firewater mixing ratio
    qv(ij,1) = xwater + qvenv(ij,1)  !plus what's already there
  END DO
  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    ! pe esta em kpa  - esat do rams esta em mbar = 100 pa = 0.1 kpa
    es(ij)           = 0.1*Esat_l(temp(ij,1)) !blob saturation vapor pressure, em kpa
  END DO
  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    !rotina do plumegen ja calcula em kpa
    !es           = esat_pr (t(1))  !blob saturation vapor pressure, em kpa
    est(ij,1)  = es(ij)                                  
    qsat(ij,1) = (eps * es(ij)) / (pe(ij,1) - es(ij))   !blob saturation lwc g/g dry air
    IF (qv(ij,1) .gt. qsat (ij,1) ) THEN  
       qc(ij,1) = qv(ij,1) - qsat(ij,1) + qc(ij,1)  !remainder goes into cloud drops
       qv(ij,1) = qsat(ij,1)  
    ENDIF  
  END DO
  
  CALL Waterbal(nkp,l,qc,qh,qi,qv,wbar,dqsdz,dt,temp,rho,est,cvi,qsat,ijbeg,ijend,isdone,nm1_l)

END SUBROUTINE Lbound

!-------------------------------------------------------------------------------

SUBROUTINE Waterbal(nkp,l,qc,qh,qi,qv,wbar,dqsdz,dt,temp,rho,est,cvi,qsat,ijbeg,ijend,isdone,nm1)  

  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: nkp,l,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend):: wbar,dqsdz,dt
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend):: isdone
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend):: nm1

  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)  :: qsat,rho,cvi,est
  
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)        :: qv,qh,qi,qc,temp
  LOGICAL,DIMENSION(ijbeg:ijend):: nottodo
  INTEGER :: ij
 
  DO ij=ijbeg,ijend
    IF(isdone(ij) .OR. l>nm1(ij)-1) THEN
      nottodo(ij)=.true.
    ELSE
      nottodo(ij)=.false.
    END IF
  END DO
 
  DO ij=ijbeg,ijend
    IF(nottodo(ij)) CYCLE
    IF (qc(ij,l) .le.1.0e-10) qc(ij,l) = 0.  !defeat underflow problem
    IF (qh(ij,l) .le.1.0e-10) qh(ij,l) = 0.  
    IF (qi(ij,l) .le.1.0e-10) qi(ij,l) = 0.  
  END DO
  
  !CALL ftrace_region_begin('Evaporate')
  !vapor to cloud,cloud to vapor  
  CALL Evaporate(l,nkp,qsat,qv,wbar,dqsdz,dt,qh,qi,qc,temp,rho,est,cvi,ijbeg,ijend,nottodo)    
  !CALL ftrace_region_end('Evaporate')

  !CALL ftrace_region_begin('Sublimate')
  !vapor to ice                              
  CALL Sublimate(l,nkp,qv,qi,temp,dt,qsat,rho,est,ijbeg,ijend,nottodo)    
  !CALL ftrace_region_end('Sublimate')
  
  !CALL ftrace_region_begin('Glaciate')
  !rain to ice                            
  CALL Glaciate(l,nkp,qh,qi,temp,qsat,dt,qv,ijbeg,ijend,nottodo)    
  !CALL ftrace_region_end('Glaciate')
  
  !CALL ftrace_region_begin('Melt')
  !ice to rain                          
  CALL Melt(l,nkp,rho,cvi,dt,qh,qi,temp,ijbeg,ijend,nottodo)
  !CALL ftrace_region_end('Melt')

  !CALL ftrace_region_begin('Convert')
  !(auto)conversion and accretion
  CALL Convert(l,nkp,qc,qh,rho,temp,dt,ijbeg,ijend,nottodo)
 ! CALL ftrace_region_end('Convert')
  
END SUBROUTINE Waterbal

!-------------------------------------------------------------------------------

SUBROUTINE Evaporate(l,nkp,qsat,qv,wbar,dqsdz,dt,qh,qi,qc,temp,rho,est,cvi,ijbeg,ijend,nottodo)  
  !
  !- evaporates cloud,rain and ice to saturation
  !
  IMPLICIT NONE
  !
  !        XNO=10.0E06
  !        HERC = 1.93*1.E-6*XN035        !evaporation constant
  !
  REAL,PARAMETER :: herc = 5.44e-4, cp = 1.004, heatcond = 2.5e3  
  REAL,PARAMETER :: heatsubl = 2834., tmelt = 273., tfreeze = 269.3
  REAL,PARAMETER :: frc = heatcond / cp, src = heatsubl / cp
  
  INTEGER,INTENT(IN) :: l,nkp,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend):: wbar,dqsdz,dt
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend):: nottodo
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)  :: qsat,rho,cvi,est

  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: qv,qh,qi,qc,temp

  REAL :: evrate, quant, dividend, divisor, devidt
  REAL   ,DIMENSION(ijbeg:ijend) :: sd,evhdt,evidt,evap
  INTEGER,DIMENSION(ijbeg:ijend) :: ijcount
  LOGICAL,DIMENSION(ijbeg:ijend) :: test1,test2,test3
  INTEGER :: ij,icount,i,ic
  !
  !
  test1=.false.
  test2=.false.
  test3=.false.
  
  icount=ijbeg-1
  DO ij=ijbeg,ijend
    IF(nottodo(ij)) CYCLE
    sd(ij) = qsat (ij,l) - qv (ij,l)  !vapor deficit
    IF(sd(ij)==0.0) CYCLE
    icount=icount+1
    ijcount(icount)=ij
  END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    evhdt(ij) = 0.  
    evidt(ij) = 0.  
    evrate = abs (wbar(ij) * dqsdz(ij))   !evaporation rate (Kessler 8.32)
    evap(ij) = evrate * dt(ij)            !what we can get in DT
  END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    IF(sd(ij).le.0.0 .AND. evap(ij) >= abs(sd(ij))) test1(ij)=.true.
    IF(sd(ij).le.0.0 .AND. evap(ij) <  abs(sd(ij))) test2(ij)=.true.
  END DO
  
  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (test1(ij)) THEN  
      qc (ij,l) = qc  (ij,l) - sd(ij)  !deficit,remember?
      qv (ij,l) = qsat(ij,l)       !set the vapor to saturation  
      temp(ij,l) = temp(ij,l) - sd(ij) * frc  !heat gained through condensation
    END IF
  END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (test2(ij)) THEN  
      qc (ij,l) = qc (ij,l) + evap(ij)         !get what we can in DT
      qv (ij,l) = qv (ij,l) - evap(ij)         !remove it from the vapor
      temp(ij,l) = temp(ij,l) + evap(ij) * frc   !get some heat
    END IF
  END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (sd(ij)>0.0) THEN    !SD is positive, need some water
    !
    ! not saturated. saturate if possible. use everything in order
    ! cloud, rain, ice. SD is positive
       IF (evap(ij).le.qc (ij,l) ) THEN        !enough cloud to last DT  
            IF (sd(ij).le.evap(ij)) THEN          !enough time to saturate
               qc (ij,l) = qc (ij,l) - sd(ij)          !remove cloud                                        
               qv (ij,l) = qsat (ij,l)                !saturate
               temp(ij,l) = temp(ij,l) - sd(ij) * frc   !cool the parcel                                    
               CYCLE  !done
    !
                                              
            ELSE   !not enough time
                                              
               sd(ij) = sd(ij)-evap(ij)                !use what there is
               qv (ij,l) = qv (ij,l) + evap(ij)          !add vapor
               temp(ij,l) = temp(ij,l) - evap(ij) * frc !lose heat
               qc (ij,l) = qc (ij,l) - evap(ij)          !lose cloud
                                          !go on to rain.                                      
            ENDIF    
       ELSE                   !not enough cloud to last DT
            IF (sd(ij).le.qc (ij,l) ) THEN   !but there is enough to sat
               qv (ij,l) = qsat (ij,l)  !use it
               qc (ij,l) = qc (ij,l) - sd(ij)  
               temp(ij,l) = temp(ij,l) - sd(ij) * frc  
               CYCLE  
            ELSE            !not enough to sat
               sd(ij) = sd(ij)-qc (ij,l)  
               qv (ij,l) = qv (ij,l) + qc (ij,l)  
               temp(ij,l) = temp(ij,l) - qc (ij,l) * frc    
               qc (ij,l) = 0.0  !all gone
            ENDIF       !on to rain                            
       ENDIF              !finished with cloud
   END IF
 END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (sd(ij)>0.0) THEN    !SD is positive, need some water
    !
    !  but still not saturated, so try to use some rain
    !  this is tricky, because we only have time DT to evaporate. if there
    !  is enough rain, we can evaporate it for dt. ice can also sublimate
    !  at the same time. there is a compromise here.....use rain first, then
    !  ice. saturation may not be possible in one DT time.
    !  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
    !  sd is still positive or we wouldn't be here.
    
       IF (qh (ij,l) .le.1.e-10) CYCLE                              
    
    !srf-25082005
    !  QUANT = ( QC (L)  + QV (L) - QSAT (L) ) * RHO (L)   !g/m**3
       quant = ( qsat (ij,l)- qc (ij,l) - qv (ij,l)   ) * rho (ij,l)   !g/m**3
    !
       evhdt(ij) = (dt(ij) * herc * (quant) * (qh (ij,l) * rho (ij,l) ) **.65) / rho (ij,l)
    !                  rain evaporation in time DT

    END IF
  END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (sd(ij)>0.0) THEN    !SD is positive, need some water
       IF (qh (ij,l) .le.1.e-10) CYCLE                              
       IF (evhdt(ij).le.qh (ij,l) ) THEN    !enough rain to last DT
            IF (sd(ij).le.evhdt(ij)) THEN             !enough time to saturate              
               qh (ij,l) = qh (ij,l) - sd(ij)            !remove rain      
               qv (ij,l) = qsat (ij,l)                    !saturate              
               temp(ij,l) = temp(ij,l) - sd(ij) * frc            !cool the parcel                      
               CYCLE                            !done
            ELSE                               !not enough time
               sd(ij) = sd(ij)-evhdt(ij)               !use what there is
               qv (ij,l) = qv (ij,l) + evhdt(ij)         !add vapor
               temp(ij,l) = temp(ij,l) - evhdt(ij) * frc     !lose heat
               qh (ij,l) = qh (ij,l) - evhdt(ij)         !lose rain
            ENDIF                               !go on to ice.
       ELSE  !not enough rain to last DT
            IF (sd(ij).le.qh (ij,l) ) THEN      !but there is enough to sat
               qv (ij,l) = qsat (ij,l)                      !use it
               qh (ij,l) = qh (ij,l) - sd(ij)  
               temp(ij,l) = temp(ij,l) - sd(ij) * frc  
               CYCLE  
            ELSE                              !not enough to sat
               sd(ij) = sd(ij)-qh (ij,l)  
               qv (ij,l) = qv (ij,l) + qh (ij,l)  
               temp(ij,l) = temp(ij,l) - qh (ij,l) * frc    
               qh (ij,l) = 0.0                    !all gone
            ENDIF                             !on to ice
       ENDIF                                    !finished with rain
    !
    !
    !  now for ice
    !  equation from (OT); correction factors for units applied
    !
   END IF
 END DO

  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (sd(ij)>0.0) THEN    !SD is positive, need some water
       IF (qi (ij,l) .le.1.e-10) CYCLE         !no ice there
       dividend = ( (1.e6 / rho (ij,l) ) **0.475) * (sd(ij) / qsat (ij,l) &
                  - 1) * (qi (ij,l) **0.525) * 1.13
       divisor = 7.e5 + 4.1e6 / (10. * est (ij,l) )  
       devidt = - cvi(ij,l) * dividend / divisor   !rate of change
       evidt(ij) = devidt * dt(ij)                        !what we could get
   END IF
 END DO


  DO i=ijbeg,icount
    ij=ijcount(i)
    IF (sd(ij)>0.0) THEN    !SD is positive, need some water
       IF (qi (ij,l) .le.1.e-10) CYCLE         !no ice there
    !
    ! logic here is identical to rain. could get fancy and make subroutine
    ! but duplication of code is easier. God bless the screen editor.
    !
       IF (evidt(ij).le.qi (ij,l) ) THEN      !enough ice to last DT
            IF (sd(ij).le.evidt(ij)) THEN               !enough time to saturate
               qi (ij,l) = qi (ij,l) - sd(ij)               !remove ice
               qv (ij,l) = qsat (ij,l)                     !saturate
               temp(ij,l) = temp(ij,l) - sd(ij) * src           !cool the parcel
               CYCLE                              !done
            ELSE                                !not enough time
               sd(ij) = sd(ij)-evidt(ij)                !use what there is
               qv (ij,l) = qv (ij,l) + evidt(ij)        !add vapor
                temp(ij,l) =  temp(ij,l) - evidt(ij) * src      !lose heat
               qi (ij,l) = qi (ij,l) - evidt(ij)        !lose ice
            ENDIF                               !go on,unsatisfied
       ELSE                                      !not enough ice to last DT
            IF (sd(ij).le.qi (ij,l) ) THEN      !but there is enough to sat
               qv (ij,l) = qsat (ij,l)                      !use it
               qi (ij,l) = qi   (ij,l) - sd(ij)  
                temp(ij,l) =  temp(ij,l) - sd(ij) * src  
               CYCLE  
            ELSE                                 !not enough to sat
               sd(ij) = sd(ij)-qi (ij,l)  
               qv (ij,l) = qv (ij,l) + qi (ij,l)  
               temp(ij,l) = temp(ij,l) - qi (ij,l) * src        
               qi (ij,l) = 0.0                       !all gone
            ENDIF                                !on to better things
                                                 !finished with ice
       ENDIF  
    ENDIF                                       !finished with the SD decision
  END DO

END SUBROUTINE Evaporate

!-----------------------------------

SUBROUTINE Sublimate(l,nkp,qv,qi,temp,dt,qsat,rho,est,ijbeg,ijend,nottodo)  
  !
  ! ********************* VAPOR TO ICE (USE EQUATION OT22)***************
  !
  REAL,PARAMETER :: eps = 0.622, heatfus = 334., heatsubl = 2834., cp = 1.004
  REAL,PARAMETER :: src = heatsubl / cp, frc = heatfus / cp, tmelt = 273.3
  REAL,PARAMETER :: tfreeze = 269.3

  INTEGER,INTENT(IN)                       :: l,nkp,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend)            :: dt
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend)            :: nottodo
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)        :: qsat,rho,est
  
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)     :: qv,qi,temp

  REAL ::dtsubh,  dividend,divisor
  REAL,DIMENSION(ijbeg:ijend) :: subl
  INTEGER,DIMENSION(ijbeg:ijend) :: ijcount
  INTEGER :: ij,count,i
  !                                      !
  !selection rules
  
  count=ijbeg-1
  DO ij=ijbeg,ijend
    IF(nottodo(ij)) CYCLE  
    IF (temp(ij,l)  .gt. tfreeze  ) CYCLE  
    IF (qv (ij,l) .le. qsat (ij,l) ) CYCLE  
    count=count+1
    ijcount(count)=ij
  END DO
  !
  DO i=ijbeg,count
    ij=ijcount(i)
    !
    !          from (OT); correction factors for units applied
    !
    dividend = ( (1.e6 / rho (ij,l) ) **0.475) * (qv (ij,l) / qsat (ij,l) &
                - 1) * (qi (ij,l) **0.525) * 1.13
    divisor = 7.e5 + 4.1e6 / (10. * est (ij,l) )  
    !
                                            
    dtsubh = abs (dividend / divisor)        !sublimation rate
    subl(ij) = dtsubh * dt(ij)                  !and amount possible
  END DO
  
  DO i=ijbeg,count
    ij=ijcount(i)
    !
    !          again check the possibilities
    !
    IF (subl(ij).lt.qv (ij,l) ) THEN  
      qv (ij,l) = qv (ij,l) - subl(ij)             !lose vapor
      qi (ij,l) = qi (ij,l) + subl(ij)          !gain ice
      temp(ij,l) = temp(ij,l) + subl(ij) * src         !energy change, warms air
    ELSE  
      qi (ij,l) = qv (ij,l)                       !use what there is
      temp(ij,l) = temp(ij,l) + qv (ij,l) * src      !warm the air
      qv (ij,l) = 0.0  
    END IF  
  END DO

END SUBROUTINE Sublimate

!------------------------------------------------------

SUBROUTINE Glaciate(l,nkp,qh,qi,temp,qsat,dt,qv,ijbeg,ijend,nottodo)  
  !
  ! *********************** CONVERSION OF RAIN TO ICE *******************
  !     uses equation OT 16, simplest. correction from W not applied, but
  !     vapor pressure differences are supplied.  
  !
  REAL,PARAMETER :: heatfus = 334., cp = 1.004, eps = 0.622, heatsubl = 2834.
  REAL,PARAMETER :: frc = heatfus / cp, frs = heatsubl / cp, tfreeze =  269.3
  REAL,PARAMETER :: glconst = 0.025   !glaciation time constant, 1/sec
  
  INTEGER,INTENT(IN) :: l,nkp,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend)            :: dt
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend)            :: nottodo
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)        :: qsat,qv

  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)        :: qh,qi,temp
  REAL,DIMENSION(ijbeg:ijend) :: dfrzh !rate of mass gain in ice
  INTEGER,DIMENSION(ijbeg:ijend) :: ijcount
  INTEGER :: ij,count,i
  !                                      !
  !selection rules
  
  count=ijbeg-1
  DO ij=ijbeg,ijend
    IF(nottodo(ij)) CYCLE  
    IF (qh (ij,l) .le. 0.    ) CYCLE
    IF (qv (ij,l) .lt. qsat (ij,l) ) CYCLE                                            
    IF (temp(ij,l) .gt. tfreeze  ) CYCLE  
    count=count+1
    ijcount(count)=ij
  END DO
  !
  !
  DO i=ijbeg,count
    ij=ijcount(i)
    dfrzh(ij) = dt(ij) * glconst * qh (ij,l)        ! from OT(16)
  END DO

  DO i=ijbeg,count
    ij=ijcount(i)
    !
    IF (dfrzh(ij).lt.qh (ij,l) ) THEN  
      qi (ij,l) = qi (ij,l) + dfrzh(ij)  
      qh (ij,l) = qh (ij,l) - dfrzh(ij)  
      temp(ij,l) = temp(ij,l) + frc * dfrzh(ij)  !warms air
    ELSE  
      qi (ij,l) = qi (ij,l) + qh (ij,l)  
      temp(ij,l) = temp(ij,l) + frc * qh (ij,l)  
      qh (ij,l) = 0.0  
    END IF  
  END DO

END SUBROUTINE Glaciate

!------------------------------------------------------

SUBROUTINE Melt(l,nkp,rho,cvi,dt,qh,qi,temp,ijbeg,ijend,nottodo)  
  IMPLICIT NONE
  ! makes water out of ice
                                                
  REAL,PARAMETER :: frc = 332.27, tmelt = 273., f0 = 0.75   !ice velocity factor
  INTEGER,INTENT(IN) :: l,nkp,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend)            :: dt
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend)            :: nottodo
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp)        :: rho,cvi
  
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)  ::  qh,qi,temp
  REAL,DIMENSION(ijbeg:ijend) :: dtmelt
  INTEGER,DIMENSION(ijbeg:ijend) :: ijcount
  INTEGER :: ij,count,i
  !                                      !
  !selection rules
  
  count=ijbeg-1
  DO ij=ijbeg,ijend
    IF(nottodo(ij)) CYCLE  
    IF(qi (ij,l) .le. 0.0) CYCLE
    IF(temp(ij,l)  .lt. tmelt) CYCLE
    count=count+1
    ijcount(count)=ij
  END DO
    
  DO i=ijbeg,count
    ij=ijcount(i)
    dtmelt(ij) = dt(ij) * (2.27 / rho (ij,l) ) * cvi(ij,l) * (temp(ij,l) - tmelt) * ( (rho(ij,l)  &
           * qi (ij,l) * 1.e-6) **0.525) * (f0** ( - 0.42) )
  END DO
  !
  !     check the possibilities  
  !
  DO i=ijbeg,count
    ij=ijcount(i)
    IF (dtmelt(ij).lt.qi (ij,l) ) THEN  
      !
      qh (ij,l) = qh (ij,l) + dtmelt(ij)  
      qi (ij,l) = qi (ij,l) - dtmelt(ij)  
      temp(ij,l) = temp(ij,l) - frc * dtmelt(ij)     !cools air
    ELSE  
      qh (ij,l) = qh (ij,l) + qi (ij,l)   !get all there is to get
      temp(ij,l) = temp(ij,l) - frc * qi (ij,l)  
      qi (ij,l) = 0.0  
    END IF  
  END DO

END SUBROUTINE Melt

!------------------------------------------------------

SUBROUTINE Convert(l,nkp,qc,qh,rho,temp,dt,ijbeg,ijend,nottodo)  
  !- accretion and autoconversion
  IMPLICIT NONE
  !
  REAL   , PARAMETER :: ak1 = 0.001    !conversion rate constant
  REAL   , PARAMETER :: ak2 = 0.0052   !collection (accretion) rate
  REAL   , PARAMETER :: th  = 0.5   !kessler threshold
  INTEGER, PARAMETER :: iconv = 1        !- Kessler conversion (=0)
  !REAL  , PARAMETER :: anbase =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime)
  REAL   , PARAMETER :: anbase =100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental)
  !REAL  , PARAMETER :: bdisp = 0.366     !Berry--size dispersion (maritime)
  REAL   , PARAMETER :: bdisp = 0.146    !Berry--size dispersion (continental)
  REAL   , PARAMETER :: tfreeze = 269.3  !ice formation temperature  
  
  INTEGER,INTENT(IN)                   :: l,nkp,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend)        :: dt
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend)        :: nottodo
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: rho
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: temp
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: qh,qc

  REAL ::   accrete, con, q, h, bc1,   bc2
  REAL,DIMENSION(ijbeg:ijend) :: total
  INTEGER,DIMENSION(ijbeg:ijend) :: ijcount
  INTEGER :: ij,count,i
  !                                      !
  !selection rules
  
  count=ijbeg-1
  DO ij=ijbeg,ijend
    IF(nottodo(ij)) CYCLE
    IF (temp(ij,l)<=tfreeze) CYCLE  !process not allowed above ice
    IF (qc(ij,l)==0.  ) CYCLE  
    count=count+1
    ijcount(count)=ij
  END DO
  
  DO i=ijbeg,count
    ij=ijcount(i)
    accrete = 0.  
    con = 0.  
    q = rho (ij,l) * qc (ij,l)  
    h = rho (ij,l) * qh (ij,l)  
    
    !          selection rules
    IF (qh (ij,l) .gt. 0.  ) accrete = ak2 * q * (h**.875)  !accretion, Kessler
    !
    !IF (iconv.ne.0) THEN   !select Berry or Kessler
      !old   bc1 = 120.  
      !old   bc2 = .0266 * anbase * 60.  
      !old   con = bdisp * q * q * q / (bc1 * q * bdisp + bc2)      
    con = q*q*q*bdisp/(60.*(5.*q*bdisp+0.0366*anbase))
    !ELSE  
      !   con = ak1 * (q - th)   !kessler autoconversion rate
      !   if (con.lt.0.0) con = 0.0   !havent reached threshold
    !  con = max(0.,ak1 * (q - th)) ! versao otimizada
    !END IF  
    
    total(ij) = (con + accrete) * dt(ij) / rho (ij,l)  
  END DO

  DO i=ijbeg,count
    ij=ijcount(i)
    
    IF (total(ij).lt.qc (ij,l) ) THEN  
      qc (ij,l) = qc (ij,l) - total(ij)  
      qh (ij,l) = qh (ij,l) + total(ij)    !no phase change involved
    ELSE  
      qh (ij,l) = qh (ij,l) + qc (ij,l)    !uses all there is
      qc (ij,l) = 0.0  
    ENDIF  
  
  END DO
  
END SUBROUTINE Convert

!-------------------------------------------------------------------------------

SUBROUTINE Tend0_Plumerise(nm1,wt,tt,qvt,qct,qht,qit,sct,nkp,ijbeg,ijend,isdone)
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)  :: nkp,ijbeg,ijend
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend)  :: nm1
  LOGICAL,DIMENSION(ijbeg:ijend) :: isdone
  REAL,DIMENSION(ijbeg:ijend,nkp),INTENT(INOUT) ::  wt,tt,qvt,qct,qht,qit,sct
  INTEGER :: nma,ij
  
  DO nma=1,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. nma>nm1(ij)) CYCLE
      wt(ij,nma)   = 0.
      tt(ij,nma)   = 0.
      qvt(ij,nma)  = 0.
      qct(ij,nma)  = 0.
      qht(ij,nma)  = 0.
      qit(ij,nma)  = 0.
      !sct(1:nm1)  = 0.
    END DO
  END DO
  
END SUBROUTINE Tend0_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE Vel_Advectc_Plumerise(nm1,wc,wt,rho,dzm,nkp,ijbeg,ijend,isdone)
  IMPLICIT NONE

  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend)  :: nm1
  INTEGER,               INTENT(IN)    ::nkp,ijbeg,ijend
  REAL   ,DIMENSION(ijbeg:ijend,nkp),INTENT(IN)    ::  wc,rho
  REAL   ,DIMENSION(nkp),INTENT(IN)    ::  dzm
  LOGICAL,DIMENSION(ijbeg:ijend),INTENT(IN) :: isdone
  REAL   ,DIMENSION(ijbeg:ijend,nkp),INTENT(INOUT) ::  wt

  REAL,DIMENSION(ijbeg:ijend,nkp) :: dn0,flxw ! var local

  REAL :: c1z
  INTEGER :: k,ij
  
  !dzm(:)= 1./dz
  DO k=1,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)) CYCLE
      dn0(ij,k)=rho(ij,k)*1.e-3 ! converte de cgs para mks
    END DO
  END DO
  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    flxw(ij,1) = wc(ij,1) * dn0(ij,1)
  END DO    
  DO k=2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
     flxw(ij,k) = wc(ij,k) * .5 * (dn0(ij,k) + dn0(ij,k+1))
    END DO
  END DO

  ! Compute advection contribution to W tendency
  c1z = .5
  DO k = 2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)-2) CYCLE      
      wt(ij,k) = wt(ij,k)  &
      + c1z * dzm(k) / (dn0(ij,k) + dn0(ij,k+1)) *     (   &
      (flxw(ij,k) + flxw(ij,k-1))  * (wc(ij,k) + wc(ij,k-1))   &
      - (flxw(ij,k) + flxw(ij,k+1))  * (wc(ij,k) + wc(ij,k+1))   &
      + (flxw(ij,k+1) - flxw(ij,k-1)) * 2.* wc(ij,k)       )
    END DO
  END DO

END SUBROUTINE Vel_Advectc_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE Scl_Advectc_Plumerise(varn,nm1,nkp,dt,wc,vvel,rho,dzm,dzt,zm,zt, &
                                 temp,tt,qv,qvt,qc,qct,qh,qht,qi,qit,ijbeg,ijend,isdone)

  IMPLICIT NONE
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend) :: nm1
  INTEGER,INTENT(IN)                    :: nkp,ijbeg,ijend
  CHARACTER(LEN=*)   ,INTENT(IN)    :: varn
  REAL,DIMENSION(ijbeg:ijend),INTENT(IN)    :: dt
  LOGICAL,DIMENSION(ijbeg:ijend),INTENT(IN)    :: isdone
  REAL,DIMENSION(ijbeg:ijend,nkp),INTENT(IN)    :: wc,rho,vvel
  REAL,DIMENSION(nkp),INTENT(IN)    :: dzm,dzt,zm,zt

  REAL,DIMENSION(ijbeg:ijend,nkp),INTENT(INOUT) :: qc,qh,qi,temp
  REAL,DIMENSION(ijbeg:ijend,nkp),INTENT(INOUT) :: tt,qvt,qct,qht,qit,qv
  
  REAL,DIMENSION(ijbeg:ijend)  :: dtlto2
  INTEGER           :: k,ij
  REAL,DIMENSION(ijbeg:ijend,nkp) ::  vctr1,vctr2,vt3dc,vt3df,vt3dk,vt3dg,scr1
  
  DO k=1,nkp
    DO ij=ijbeg,ijend
      vctr1(ij,k)=0.0
      vctr2(ij,k)=0.0
      vt3dc(ij,k)=0.0
      vt3df(ij,k)=0.0
      vt3dk(ij,k)=0.0
      vt3dg(ij,k)=0.0
      scr1(ij,k)=0.0
    END DO
   END DO
  
   DO ij=ijbeg,ijend
     IF(isdone(ij)) CYCLE
     !  wp => w
     !- Advect  scalars
     dtlto2(ij)   = .5 * dt(ij)
   END DO
  
   DO ij=ijbeg,ijend
     IF(isdone(ij)) CYCLE
     !  vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * dne(1)
     vt3dc(ij,1) =      (vvel(ij,1) + wc(ij,1)) * dtlto2(ij) * rho(ij,1)*1.e-3 !converte de CGS p/ MKS
     vt3df(ij,1) = .5 * (vvel(ij,1) + wc(ij,1)) * dtlto2(ij) * dzm(1)
   END DO
  
   DO k = 2,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE    
       !     vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (dne(k) + dne(k+1))
       vt3dc(ij,k) =  (vvel(ij,k) + wc(ij,k)) * dtlto2(ij) *.5 * (rho(ij,k) + rho(ij,k+1))*1.e-3
       vt3df(ij,k) =  (vvel(ij,k) + wc(ij,k)) * dtlto2(ij) *.5 *  dzm(k)
     END DO
   END DO
 
   !-srf-24082005
   !  do k = 1,nm1-1
   DO k = 1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE        
       vctr1(ij,k) = (zt(k+1) - zm(k)) * dzm(k)
       vctr2(ij,k) = (zm(k)   - zt(k)) * dzm(k)
       !   vt3dk(k) = dzt(k) / dne(k)
       vt3dk(ij,k) = dzt(k) /(rho(ij,k)*1.e-3)
     END DO
   END DO

   !      scalarp => scalar_tab(n,ngrid)%var_p
   !      scalart => scalar_tab(n,ngrid)%var_t

   !- temp advection tendency (TT)
    
   DO k=1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE        
       scr1(ij,k)=temp(ij,k)
     END DO
   END DO
  
   CALL Fa_Zc_Plumerise(nm1,temp,scr1,vt3dc,vt3df  &
                            ,vt3dg,vt3dk,vctr1,vctr2,nkp,ijbeg,ijend,isdone)
   CALL Advtndc_Plumerise(nm1,temp,scr1,tt,dt,nkp,ijbeg,ijend,isdone)

  !- water vapor advection tendency (QVT)
   DO k=1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE        
       scr1(ij,k)=qv(ij,k)
     END DO
   END DO
      
  CALL Fa_Zc_Plumerise(nm1,qv,scr1,vt3dc,vt3df  &
                            ,vt3dg,vt3dk,vctr1,vctr2,nkp,ijbeg,ijend,isdone)
  CALL Advtndc_Plumerise(nm1,qv,scr1,qvt,dt,nkp,ijbeg,ijend,isdone)
  !- liquid advection tendency (QCT)
   DO k=1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE        
       scr1(ij,k)=qc(ij,k)
     END DO
   END DO
      
  CALL Fa_Zc_Plumerise(nm1,qc,scr1,vt3dc,vt3df  &
                            ,vt3dg,vt3dk,vctr1,vctr2,nkp,ijbeg,ijend,isdone)
  CALL Advtndc_Plumerise(nm1,qc,scr1,qct,dt,nkp,ijbeg,ijend,isdone)

  !- ice advection tendency (QIT)
   DO k=1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE        
       scr1(ij,k)=qi(ij,k)
     END DO
   END DO
      
  CALL Fa_Zc_Plumerise(nm1,qi,scr1,vt3dc,vt3df  &
                            ,vt3dg,vt3dk,vctr1,vctr2,nkp,ijbeg,ijend,isdone)
  CALL Advtndc_Plumerise(nm1,qi,scr1,qit,dt,nkp,ijbeg,ijend,isdone)

  !- hail/rain advection tendency (QHT)
   DO k=1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. K>nm1(ij)) CYCLE        
       scr1(ij,k)=qh(ij,k)
     END DO
   END DO
      
  CALL Fa_Zc_Plumerise(nm1,qh,scr1,vt3dc,vt3df  &
                        ,vt3dg,vt3dk,vctr1,vctr2,nkp,ijbeg,ijend,isdone)
  CALL Advtndc_Plumerise(nm1,qh,scr1,qht,dt,nkp,ijbeg,ijend,isdone)


END SUBROUTINE Scl_Advectc_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE Fa_Zc_Plumerise(nm1,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,nkp,ijbeg,ijend,isdone)

  IMPLICIT NONE
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend) :: isdone
  INTEGER,INTENT(IN)                   :: nkp,ijbeg,ijend
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: vctr1,vctr2,vt3dc,vt3df,vt3dk
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: scp

  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: scr1,vt3dg

  INTEGER             :: k,ij
  REAL                :: dfact
  
  dfact = .5

  ! Compute scalar flux VT3DG
   DO k = 1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
       vt3dg(ij,k) = vt3dc(ij,k)*(vctr1(ij,k)*scr1(ij,k)+vctr2(ij,k)*scr1(ij,k+1)+&
                 vt3df(ij,k)*(scr1(ij,k)-scr1(ij,k+1)))
     END DO
   END DO
      
  ! Modify fluxes to retain positive-definiteness on scalar quantities.
  !    If a flux will remove 1/2 quantity during a timestep,
  !    reduce to first order flux. This will remain positive-definite
  !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
  !    both fluxes are evacuating the box.
  DO k = 1,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE    
         IF (vt3dc(ij,k) .gt. 0.) THEN
           IF (vt3dg(ij,k)*vt3dk(ij,k) .gt. dfact*scr1(ij,k)) THEN
              vt3dg(ij,k) = vt3dc(ij,k)*scr1(ij,k)
           END IF
         ELSEIF (vt3dc(ij,k) .lt. 0.) THEN
           IF (-vt3dg(ij,k)*vt3dk(ij,k+1) .gt. dfact*scr1(ij,k+1)) THEN
             vt3dg(ij,k) = vt3dc(ij,k)*scr1(ij,k+1)
           END IF
         END IF
     END DO
  END DO

  ! Compute flux divergence
  DO k = 2,nkp
     DO ij=ijbeg,ijend
       IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE    
       scr1(ij,k) = scr1(ij,k)+vt3dk(ij,k)*(vt3dg(ij,k-1)-vt3dg(ij,k)+ &
              scp(ij,k)*(vt3dc(ij,k)-vt3dc(ij,k-1)))
     END DO
  END DO

END SUBROUTINE Fa_Zc_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE Advtndc_Plumerise(nm1,scp,sca,sct,dtl,nkp,ijbeg,ijend,isdone)

  IMPLICIT NONE

  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend) :: isdone
  INTEGER,INTENT(IN)                    :: nkp,ijbeg,ijend
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend)     :: dtl
  REAL,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: scp
  REAL,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: sca
  REAL,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: sct

  INTEGER              :: k,ij
  REAL,DIMENSION(ijbeg:ijend) :: dtli

  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    dtli(ij) = 1. / dtl(ij)
  END DO
  
  DO k = 2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
      sct(ij,k) = sct(ij,k) + (sca(ij,k)-scp(ij,k)) * dtli(ij)
    END DO
  END DO

END SUBROUTINE Advtndc_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE scl_misc(nm1,nkp,ijbeg,ijend,qvenv,te,vvel,temp,qv,qc,qh,qi, &
                    tt,qvt,qct,qht,qit,radius,alpha,adiabat,wbar,isdone)

  USE plume_utils, ONLY: ijindex

  IMPLICIT NONE
  REAL, PARAMETER    :: g = 9.81, cp=1004.

  INTEGER,INTENT(IN) ,DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN) ,DIMENSION(ijbeg:ijend) :: isdone
  INTEGER,INTENT(IN)                      :: nkp,ijbeg,ijend
  REAL   ,INTENT(IN)                      :: alpha
  REAL   ,INTENT(IN),DIMENSION(ijbeg:ijend,nkp) :: qvenv,te,qv,radius,qc,qh,qi,vvel,temp

  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp)  :: tt,qvt,qct,qht,qit
  REAL   ,INTENT(OUT)                     :: adiabat
  REAL   ,INTENT(OUT),DIMENSION(ijbeg:ijend) :: wbar

  INTEGER            :: k,ij
  REAL              :: dmdtm

  DO k=2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
      wbar(ij)    = 0.5*(vvel(ij,k)+vvel(ij,k-1))  
      !-- dry adiabat
      adiabat = - wbar(ij) * g / cp
      !-- entrainment    
      dmdtm = 2. * alpha * abs (wbar(ij)) / radius (ij,k)  != (1/m)dm/dt
      !-- tendency temperature = adv + adiab + entrainment
      tt(ij,k) = tt(ij,k) + adiabat - dmdtm * (temp(ij,k)-te(ij,k))
      !-- tendency water vapor = adv  + entrainment
      qvt(ij,k) = qvt(ij,k)- dmdtm * (qv(ij,k) - qvenv(ij,k))
      qct(ij,k) = qct(ij,k)- dmdtm * (qc(ij,k))
      qht(ij,k) = qht(ij,k)- dmdtm * (qh(ij,k))
      qit(ij,k) = qit(ij,k)- dmdtm * (qi(ij,k))
      !-- tendency gas/particle = adv  + entrainment
      !      sct(k) = sct(k)         - dmdtm * ( sc (k) -   sce (k) )
    END DO
  END DO

END SUBROUTINE Scl_Misc


!-------------------------------------------------------------------------------

SUBROUTINE Damp_Grav_Wave(ifrom,nm1,deltak,dt,zt,zm,vvel,temp,tt,qv,qh,qi,qc,&
                          te,pe,qvenv,nkp,ijbeg,ijend,isdone)

  IMPLICIT NONE
  REAL,PARAMETER :: distim = 60.
  INTEGER,INTENT(IN)                    :: ifrom,deltak,nkp,ijbeg,ijend
  INTEGER,INTENT(IN)   , DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN)   , DIMENSION(ijbeg:ijend) :: isdone
  REAL   ,INTENT(IN)   , DIMENSION(ijbeg:ijend)     :: dt
  REAL   ,INTENT(IN)   , DIMENSION(ijbeg:ijend,nkp) :: qh,qi,qc
  REAL   ,INTENT(IN)   , DIMENSION(nkp) :: zt,zm
  REAL   ,INTENT(IN)   , DIMENSION(ijbeg:ijend,nkp) :: te,pe,qvenv,qv
  REAL   ,INTENT(INOUT), DIMENSION(ijbeg:ijend,nkp) :: tt,vvel,temp

  REAL,DIMENSION(ijbeg:ijend,nkp) :: dummy
  REAL,DIMENSION(ijbeg:ijend) :: zmkf,ztop,c1,c2
  INTEGER               :: k,nfpt,kf,ij
  
  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    kf = nm1(ij) - int(deltak/2)
    zmkf(ij) = zm(kf) !old: float(kf )*dz
    ztop(ij) = zm(nm1(ij))
    c1(ij) = 1. / (distim * (ztop(ij) - zmkf(ij)))
    c2(ij) = dt(ij) * c1(ij)
  END DO

  IF(ifrom==1) THEN !Friction
    DO ij=ijbeg,ijend
      IF(isdone(ij)) CYCLE
      DO k = nm1(ij),2,-1
        IF (zt(k) .le. zmkf(ij)) CYCLE
        tt(ij,k) = tt(ij,k)   + c1(ij) * (zt(k) - zmkf(ij))*(te(ij,k) - temp(ij,k))
      END DO
    END DO
  END IF
  
  IF(ifrom==2) THEN !Friction
    DO k=1,nkp
      DO ij=ijbeg,ijend
        IF(isdone(ij)) CYCLE
        dummy(ij,k) = 0.
      END DO
    END DO  
    DO ij=ijbeg,ijend
      IF(isdone(ij)) CYCLE
      DO k = nm1(ij),2,-1
        IF (zt(k) .le. zmkf(ij)) CYCLE
        vvel(ij,k) =  vvel(ij,k) + c2(ij) * (zt(k) - zmkf(ij))*(dummy(ij,k) - vvel(ij,k))
      END DO
    END DO
 END IF  

END SUBROUTINE Damp_Grav_Wave


SUBROUTINE Fallpart(nm1,nkp,qvt,qct,qht,qit,rho,vvel,qh,qi,zm, &
                    vth,vti,cvi,ijbeg,ijend,isdone)

  IMPLICIT NONE

  REAL, PARAMETER :: vconst = 5.107387, eps = 0.622, f0 = 0.75  
  REAL, PARAMETER :: g = 9.81, cp = 1004.

  INTEGER,INTENT(IN) :: nkp,ijbeg,ijend
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend) :: isdone
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: qvt,qct,rho,qh,qi,vvel
  REAL   ,INTENT(IN)   ,DIMENSION(nkp) :: zm
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: qht,qit
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: vth,vti,cvi

  INTEGER          :: k,ij
  REAL            :: vtc,dfhz,dfiz,dz1,vhrel,virel

  !srf==================================
  !   verificar se o gradiente esta correto
  !  
  !srf==================================
  !
  !     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
  !     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
  !     for ice, see (OT18), use F0=0.75 per argument there. rho*q
  !     values are in g/m**3, velocities in m/s
  !
  DO k=2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
      vtc = vconst * rho (ij,k) **.125   ! median volume fallspeed (ktable4)
      !  hydrometeor assembly velocity calculations (k table4)
      !  vth(k)=-vtc*qh(k)**.125  !median volume fallspeed, water            
      vth (ij,k) = - 4.            !small variation with qh
      vhrel = vvel(ij,k) + vth (ij,k)  !relative to surrounding cloud
      !  rain ventilation coefficient for evaporation
      !cvh(k) = 1.6 + 0.57e-3 * (abs (vhrel) ) **1.5  
      !  vti(k)=-vtc*f0*qi(k)**.125    !median volume fallspeed,ice            
      vti (ij,k) = - 3.                !small variation with qi
      virel = vvel(ij,k) + vti (ij,k)       !relative to surrounding cloud
      !  ice ventilation coefficient for sublimation
      cvi(ij,k) = 1.6 + 0.57e-3 * (abs (virel) ) **1.5 / f0  
      IF (vhrel.ge.0.0) THEN  
        dfhz=qh(ij,k)*(rho(ij,k  )*vth(ij,k  )-rho(ij,k-1)*vth(ij,k-1))/rho(ij,k-1)
      ELSE  
       dfhz=qh(ij,k)*(rho(ij,k+1)*vth(ij,k+1)-rho(ij,k  )*vth(ij,k  ))/rho(ij,k)
      END IF  
      IF (virel.ge.0.0) THEN  
        dfiz=qi(ij,k)*(rho(ij,k  )*vti(ij,k  )-rho(ij,k-1)*vti(ij,k-1))/rho(ij,k-1)
      ELSE  
        dfiz=qi(ij,k)*(rho(ij,k+1)*vti(ij,k+1)-rho(ij,k  )*vti(ij,k  ))/rho(ij,k)
      END IF
      dz1=zm(k)-zm(k-1)
      qht(ij,k) = qht(ij,k) - dfhz / dz1 !hydrometeors don't
      qit(ij,k) = qit(ij,k) - dfiz / dz1  !nor does ice? hail, what about
    END DO
  END DO

END SUBROUTINE Fallpart

!-------------------------------------------------------------------------------

SUBROUTINE Visc_W(nm1,deltak,kmt,nkp,vvel,temp,qv,qh,qc,qi,zm,zt,visc, &
                  wt,tt,qvt,qct,qht,qit,ijbeg,ijend,isdone)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nkp,deltak,ijbeg,ijend
  INTEGER,INTENT(IN),DIMENSION(ijbeg:ijend) :: kmt,nm1
  LOGICAL,INTENT(IN),DIMENSION(ijbeg:ijend) :: isdone
  REAL   ,INTENT(IN)   ,DIMENSION(nkp) :: zm,zt
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: visc,qv,qh,qc,qi,vvel,temp
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: wt,qct,qht,qit
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: tt,qvt
  
  INTEGER        :: k,ij
  INTEGER,DIMENSION(ijbeg:ijend) :: m2
  REAL          :: dz1t,dz1m,dz2t,dz2m,d2wdz,d2tdz,d2qvdz,d2qhdz
  REAL          :: d2qcdz ,d2qidz ,d2scdz

  !srf--- 17/08/2005
  !m2=min(m1+deltak,kmt)
  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    m2(ij)=min(nm1(ij),kmt(ij))
  END DO
  
  DO k=2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>m2(ij)-1) CYCLE
      dz1t   = 0.5*(zt(k+1)-zt(k-1))
      dz2t   = visc (ij,k) / (dz1t * dz1t)  
      dz1m   = 0.5*(zm(k+1)-zm(k-1))
      dz2m   = visc (ij,k) / (dz1m * dz1m)  
      d2wdz  = (vvel(ij,k + 1) - 2 * vvel(ij,k) + vvel(ij,k - 1) ) * dz2m  
      d2tdz  = (temp(ij,k + 1) - 2 * temp(ij,k) + temp(ij,k - 1) ) * dz2t  
      d2qvdz = (qv (ij,k + 1) - 2 * qv (ij,k) + qv (ij,k - 1) ) * dz2t  
      d2qhdz = (qh (ij,k + 1) - 2 * qh (ij,k) + qh (ij,k - 1) ) * dz2t
      d2qcdz = (qc (ij,k + 1) - 2 * qc (ij,k) + qc (ij,k - 1) ) * dz2t  
      d2qidz = (qi (ij,k + 1) - 2 * qi (ij,k) + qi (ij,k - 1) ) * dz2t  
      !d2scdz = (sc (k + 1) - 2 * sc (k) + sc (k - 1) ) * dz2t
 
      wt(ij,k) =   wt(ij,k) + d2wdz
      tt(ij,k) =   tt(ij,k) + d2tdz                          
      qvt(ij,k) =  qvt(ij,k) + d2qvdz
      qct(ij,k) =  qct(ij,k) + d2qcdz
      qht(ij,k) =  qht(ij,k) + d2qhdz
      qit(ij,k) =  qit(ij,k) + d2qidz    
      !sct(k) =  sct(k) + d2scdz
    END DO  
  END DO
  
END SUBROUTINE Visc_W

!-------------------------------------------------------------------------------

SUBROUTINE Update_Plumerise(nm1,varn,nkp,dt,wt,tt,qvt,qct,qht,qit,&
                            vvel,temp,qv,qh,qc,qi,ijbeg,ijend,isdone)
  IMPLICIT NONE

  INTEGER  ,INTENT(IN)                   :: nkp,ijbeg,ijend
  INTEGER  ,INTENT(IN)    ,DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL  ,INTENT(IN)    ,DIMENSION(ijbeg:ijend) :: isdone
  REAL     ,INTENT(IN)   ,DIMENSION(ijbeg:ijend)     :: dt
  CHARACTER,INTENT(IN)                   :: varn
  REAL     ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: wt,qct,qht,qit
  REAL     ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: tt,qvt
  REAL     ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: qv,qh,qc,qi,vvel,temp

  INTEGER :: k,ij
 
  IF(varn == 'W') THEN
    DO k=2,nkp
      DO ij=ijbeg,ijend
        IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
        vvel(ij,k) =  vvel(ij,k) +  wt(ij,k) * dt(ij)  
      END DO
    END DO
    RETURN
  ELSE
    DO k=2,nkp
      DO ij=ijbeg,ijend
        IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
        temp(ij,k) =  temp(ij,k) +  tt(ij,k)  * dt(ij)  
        qv(ij,k) = qv(ij,k) + qvt(ij,k) * dt(ij)  
        qc(ij,k) = qc(ij,k) + qct(ij,k) * dt(ij) !cloud drops travel with air
        qh(ij,k) = qh(ij,k) + qht(ij,k) * dt(ij)  
        qi(ij,k) = qi(ij,k) + qit(ij,k) * dt(ij)
        ! sc(k) = sc(k) + sct(k) * dt
        !srf---18jun2005  
        qv(ij,k) = max(0., qv(ij,k))
        qc(ij,k) = max(0., qc(ij,k))
        qh(ij,k) = max(0., qh(ij,k))
        qi(ij,k) = max(0., qi(ij,k))
        !sc(k) = max(0., sc(k))
      END DO
     END DO
  END IF
  
END SUBROUTINE Update_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE Hadvance_Plumerise(iac,nm1,dt,wc,wt,vvel,mintime,nkp,ijbeg,ijend,isdone)

  IMPLICIT NONE
  INTEGER,INTENT(IN)                   :: iac,nkp,ijbeg,ijend
  INTEGER,INTENT(IN)   ,DIMENSION(ijbeg:ijend) :: nm1,mintime
  LOGICAL,INTENT(IN)   ,DIMENSION(ijbeg:ijend) :: isdone
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend)     :: dt
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: wc,wt,vvel
  
  REAL   ,DIMENSION(ijbeg:ijend,nkp) :: dummy
  REAL,DIMENSION(ijbeg:ijend) :: eps
  INTEGER                :: k,ij
  !     It is here that the Asselin filter is applied.  For the velocities
  !     and pressure, this must be done in two stages, the first when
  !     IAC=1 and the second when IAC=2.

  DO ij=ijbeg,ijend
    IF(isdone(ij)) CYCLE
    eps(ij) = .2
    IF(mintime(ij) == 1) eps(ij)=0.5
  END DO
  !     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.
  IF (iac .eq. 1) THEN
    DO k = 1,nkp
      DO ij=ijbeg,ijend
        IF(isdone(ij) .OR. k>nm1(ij)) CYCLE
        wc(ij,k) = wc(ij,k) + eps(ij) * (vvel(ij,k) - 2. * wc(ij,k))
      END DO
    END DO
    RETURN
  ELSEIF (iac .eq. 2) THEN
    DO k = 1,nkp
      DO ij=ijbeg,ijend
        IF(isdone(ij) .OR. k>nm1(ij)) CYCLE
        dummy(ij,k) = vvel(ij,k)
        vvel(ij,k) = wc(ij,k) + eps(ij) * dummy(ij,k)
      END DO
    END DO
  END IF

  DO k = 1,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)) CYCLE
      wc(ij,k) = dummy(ij,k)
    END DO
  END DO

END SUBROUTINE Hadvance_Plumerise


SUBROUTINE Buoyancy_Plumerise(nm1,temp,te,qv,qvenv,qh,qi,qc,wt,nkp,ijbeg,ijend,isdone)
  IMPLICIT NONE
  REAL, PARAMETER     :: mu = 0.15
  REAL, PARAMETER     :: g = 9.8, eps = 0.622, gama = 0.5 ! mass virtual coeff.
  ! compensa a falta do termo de aceleracao associado `as
  ! das pertubacoes nao-hidrostaticas no campo de pressao  
  REAL, PARAMETER     :: umgamai = 1./(1.+gama)
  !- new                 ! Siesbema et al, 2004
  !REAL, PARAMETER :: umgamai = 1./(1.-2.*mu)

  INTEGER,INTENT(IN)                   :: nkp,ijbeg,ijend
  INTEGER,INTENT(IN)   ,DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN)   ,DIMENSION(ijbeg:ijend) :: isdone
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: te,qvenv,qv,qh,qi,qc,temp
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: wt
 
  REAL   ,DIMENSION(ijbeg:ijend,nkp)               :: scr1

  REAL                :: tv,tve,qwtotl
  INTEGER :: k,ij
  
  DO k = 2,nkp
    DO ij=ijbeg,ijend
      IF (isdone(ij) .OR. k>nm1(ij)-1) CYCLE
      tv =   temp(ij,k) * (1. + (qv(ij,k)   /eps))/(1. + qv(ij,k)   )  !blob virtual temp.                                                  
      tve = te(ij,k) * (1. + (qvenv(ij,k)/eps))/(1. + qvenv(ij,k))  !and environment
      qwtotl = qh(ij,k) + qi(ij,k) + qc(ij,k)                         ! QWTOTL*G is drag
      !- orig
      !scr1(ij,k)= g*( umgamai*(  tv - tve) / tve   - qwtotl)
      scr1(ij,k)= g*  umgamai*( (tv - tve) / tve   - qwtotl)
    END DO
  END DO

  DO k = 2,nkp
    DO ij=ijbeg,ijend
      IF (isdone(ij) .OR. k>nm1(ij)-2) CYCLE
      wt(ij,k) = wt(ij,k)+0.5*(scr1(ij,k)+scr1(ij,k+1))
    END DO
  END DO

END SUBROUTINE Buoyancy_Plumerise

!-------------------------------------------------------------------------------

SUBROUTINE Entrainment(nm1,vvel,wt,radius,alpha,nkp,wbar,ijbeg,ijend,isdone)
  IMPLICIT NONE
  REAL,PARAMETER :: mu = 0.15 ,gama = 0.5 ! mass virtual coeff.
  ! compensa a falta do termo de aceleracao associado `as
  ! das pertubacoes nao-hidrostaticas no campo de pressao  
  REAL,PARAMETER :: umgamai = 1./(1.+gama)
  
  INTEGER,INTENT(IN)                   :: nkp,ijbeg,ijend
  INTEGER,INTENT(IN)  ,DIMENSION(ijbeg:ijend) :: nm1
  LOGICAL,INTENT(IN)  ,DIMENSION(ijbeg:ijend) :: isdone
  REAL   ,INTENT(IN)                   :: alpha
  REAL   ,INTENT(IN)   ,DIMENSION(ijbeg:ijend,nkp) :: radius,vvel
  REAL   ,INTENT(INOUT),DIMENSION(ijbeg:ijend,nkp) :: wt
  REAL   ,INTENT(OUT),DIMENSION(ijbeg:ijend)       :: wbar
  
  REAL    :: dmdtm,radius_bar
  INTEGER :: k,ij
  
  DO k=2,nkp
    DO ij=ijbeg,ijend
      IF(isdone(ij) .OR. k>nm1(ij)-1) CYCLE
      !-- for w: wbar is only w(k)
      !     wbar=0.5*(w(k)+w(k-1))          
      wbar(ij)=vvel(ij,k)          
      radius_bar = 0.5*(radius(ij,k) + radius(ij,k+1))
      ! orig
      !dmdtm =           2. * alpha * abs (wbar) / radius_bar  != (1/m)dm/dt
      dmdtm = umgamai * 2. * alpha * abs (wbar(ij)) / radius_bar  != (1/m)dm/dt
      ! dmdtm*w(l) entrainment,
      wt(ij,k) = wt(ij,k)  - dmdtm*abs (wbar(ij))
      !print*,'w-entr=',k,w(k),- dmdtm*abs (wbar)
    END DO
  END DO

END SUBROUTINE  Entrainment

!     ******************************************************************
REAL FUNCTION  Esat_l(temp)
  IMPLICIT NONE
  REAL, PARAMETER :: abz=273.15
  REAL,INTENT(IN) :: temp
  REAL :: tc
  !     esat(millibars),t(kelvin)
  tc=temp-abz
  esat_l=6.1078*exp((17.2693882*tc)/(tc+237.3))

END FUNCTION Esat_l

!-----------------------------------------------------------------------------
FUNCTION Esat_Pr(tem)  
  !
  ! ******* Vapor Pressure  A.L. Buck JAM V.20 p.1527. (1981) ***********
  !
  REAL, PARAMETER :: ci1 = 6.1115, ci2 = 22.542, ci3 = 273.48
  REAL, PARAMETER :: cw1 = 6.1121, cw2 = 18.729, cw3 = 257.87, cw4 = 227.3
  REAL, PARAMETER :: tmelt = 273.3

  REAL,INTENT(IN) :: tem

  REAL :: Esat_Pr,temc,esatm
  !
  !     formulae from Buck, A.L., JAM 20,1527-1532
  !     custom takes esat wrt water always. formula for h2o only
  !     good to -40C so:
  !  
  
  temc = tem - tmelt  
  IF (temc.gt. - 40.0) GOTO 230  
  esatm = ci1 * exp (ci2 * temc / (temc + ci3) )  !ice, millibars  
  esat_pr = esatm / 10.        !kpa                          

  RETURN  
  
  230 esatm = cw1 * exp ( ( (cw2 - (temc / cw4) ) * temc) / (temc + cw3))                          
  esat_pr = esatm / 10.        !kpa                          

END FUNCTION Esat_Pr
