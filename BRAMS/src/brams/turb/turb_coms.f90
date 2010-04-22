!==========================================================================================!
!==========================================================================================!
!    This module contains the parameters used by the Nakanishi and Niino's improved        !
! Mellor-Yamada parametrisation.  The "A", "B", and "C" constants are now provided in the  !
! namelist, so it can be tuned/optimised more easily.                                      !
!                                                                                          !
! References:                                                                              !
!                                                                                          !
! JANJIC, Z. I. Nonsingular implementation of the Mellor-Yamada level 2.5 scheme in the    !
!   NCEP meso model, office note # 437, National Centers for Environmental Prediction,     !
!   2001, 61 pp.                                                                           !
!                                                                                          !
! MELLOR, G. L.; YAMADA, T. Development of a turbulence closure model for geophysical      !
!   fluid models. Rev. Geophys. Space Phys., vol. 20, p. 851-875, 1982.                    !
!                                                                                          !
! NAKANISHI, M. Improvement of the Mellor-Yamada turbulence closure model based on         !
!   large-eddy simulation data. Boundary-Layer Meteor., v. 99, p. 349-378, 2001.           !
!                                                                                          !
! NAKANISHI, M.; NIINO, H. An improved Mellor-Yamada level-3 model with condens-           !
!   ation physics: its design and verification. Boundary-Layer Meteor., v.112,             !
!   p. 1-31, 2004.                                                                         !
!                                                                                          !
! NAKANISHI, M.; NIINO, H. An improved Mellor-Yamada level-3 model with condensation       !
!   physics: its numerical stability and application to a regional prediction of advection !
!   fog. Boundary-Layer Meteor., vol. 119, p. 397-407, 2006.                               !
!                                                                                          !
! HANNA, S. R. Application in air pollution modeling. In: NIEUSWSTADT, F. M. T.;           !
!   VAN DOP, H. Atmospheric turbulence and air pollution modelling. Dordrecht: D.          !
!   Reidel Publishing Company, 1982, chap. 7, p. 275-310.                                  !
!                                                                                          !
! HELFAND, H. M.; LABRAGA, J. C. Design of a non-singular level 2.5 second-order closure   !
!   model for the prediction of atmospheric turbulence. J. Atmos. Sci., vol. 45,           !
!   p. 113-132, 1988.                                                                      !
!                                                                                          !
! VOGEZELANG, D. H. P.; HOLTSLAG, A. M. Evaluation and model impacts of alternati-         !
!   ve boundary-layer height formulations. Boundary-Layer Meteor., v. 81, p. 245-          !
!   269, 1996.                                                                             !
!------------------------------------------------------------------------------------------!
module turb_coms

   !--- Nakanishi and Niino (2004) A, B, and C variables that are provided by the RAMSIN. -!
   real, dimension(2) :: nna
   real, dimension(2) :: nnb
   real, dimension(5) :: nnc
   
   !---------------------------------------------------------------------------------------!
   !      Local variables that will receive the namelist variables and used by Nakanishi-  !
   ! Niino clousure.                                                                       !
   !---------------------------------------------------------------------------------------!
   !----- "A" constants. ------------------------------------------------------------------!
   real :: nna1
   real :: nna2
   !----- "B" constants. ------------------------------------------------------------------!
   real :: nnb1
   real :: nnb2
   !----- "C" constants. ------------------------------------------------------------------!
   real :: nnc1
   real :: nnc2
   real :: nnc3
   real :: nnc4
   real :: nnc5
   !---------------------------------------------------------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   !      Local constants for Helfand and Labraga (1988) clousure.  This could eventually  !
   ! use the NNA, NNB, and NNC variables, but that requires removing several hard-coded    !
   ! constants in subroutine tkemy.                                                        !
   !---------------------------------------------------------------------------------------!
   !----- "A" constants. ------------------------------------------------------------------!
   real, parameter :: hla1 = 0.92
   real, parameter :: hla2 = 0.74
   !----- "B" constants. ------------------------------------------------------------------!
   real, parameter :: hlb1 = 16.6
   real, parameter :: hlb2 = 10.1
   !----- "C" constants. ------------------------------------------------------------------!
   real, parameter :: hlc1 = 0.08
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Constants to define the Turbulent length scale, and came for Nakanishi and Niino   !
   ! (2004), appendix A, equations.  These will be kept constant for the time being.       !
   !---------------------------------------------------------------------------------------!
   real, parameter :: nns1 =   3.70
   real, parameter :: nns2 =   2.70
   real, parameter :: nns3 = 100.00
   real, parameter :: nns4 =   0.20
   real, parameter :: nns5 =   0.23
   real, parameter :: nns6 =   5.00
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Turbulent parameters from Helfand and Labraga (Mellor-Yamada).  Another good      !
   ! reference to explain these variables is the HYPACT technical manual.                  !
   !---------------------------------------------------------------------------------------!
   real :: hle1
   real :: hle2
   real :: hle3
   real :: hle4
   real :: hle5
   real :: hlch
   real :: hlcm
   real :: hlrf1
   real :: hlrf2
   real :: hlrf3
   real :: hlrf4
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following parameters are derived from A, B, and C families.                   !
   !---------------------------------------------------------------------------------------!
   real :: nngama1
   real :: nngama2
   real :: nnf1
   real :: nnf2
   real :: nnrf1
   real :: nnrf2
   real :: nnri1
   real :: nnri2
   real :: nnri3
   real :: nnrfc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Some constants to compute E1, E2, etc. from Nakanishi and Niino (2004) equation   !
   ! 20.                                                                                   !
   !---------------------------------------------------------------------------------------!
   real :: nnce1a
   real :: nnce1b
   real :: nnce2
   real :: nnce3
   real :: nnce4
   real :: nncr1
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Part of the constants defined in the Janjic (2001) report, adapted to             !
   ! Nakanishi and Niino routines.                                                         !
   !---------------------------------------------------------------------------------------!
   real :: nno1
   real :: nno2
   real :: nno3
   real :: nno4
   real :: nno5
   real :: nno6
   real :: nno7
   real :: nno8
   real :: nnaeh
   real :: nnaem
   real :: nnreq
   real :: nnrsl
   real :: nnmacheps
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine assigns the derived parameters at the beginning.  This way we     !
   ! don't need to compute these constants over and over during the integration.           !
   !---------------------------------------------------------------------------------------!
   subroutine assign_turb_params()

      implicit none
      !----- Local variables. -------------------------------------------------------------!
      logical, parameter :: printinscreen = .false.
      !------------------------------------------------------------------------------------!

      !----- Transfer the "A", "B", and "C" vectors to the scalars. -----------------------!
      nna1 = nna(1)
      nna2 = nna(2)
      nnb1 = nnb(1)
      nnb2 = nnb(2)
      nnc1 = nnc(1)
      nnc2 = nnc(2)
      nnc3 = nnc(3)
      nnc4 = nnc(4)
      nnc5 = nnc(5)

      !----- Compute the derived constants. -----------------------------------------------!
      nngama1= 1./3.0-2.*nna1/nnb1
      nngama2=(2.*nna1*(3.-2.*nnc2)+nnb2*(1.-nnc3))/nnb1
      nnf1=nnb1*(nngama1-nnc1)+2.*nna1*(3.-2.*nnc2)+3.*nna2*(1.-nnc2)*(1.-nnc5)
      nnf2=nnb1*(nngama1+nngama2)-3.*nna1*(1.-nnc2)
      nnrf1=nnb1*(nngama1-nnc1)/nnf1
      nnrf2=nnb1*nngama1/nnf2
      nnri1=nna2*nnf2/(2.*nna1*nnf1)
      nnri2=nnrf1/(2.*nnri1)
      nnri3=(2.*nnrf2-nnrf1)/nnri1
      nnrfc=nngama1/(nngama1+nngama2)

      !------------------------------------------------------------------------------------!
      !      Parameters to compute E1,E2,E3,E4, and R1, from Nakanishi and Niino (2004)    ! 
      ! equation 20.                                                                       !
      !------------------------------------------------------------------------------------!
      nnce1a=6.*nna1*nna1
      nnce1b=9.*nna1*nna2*(1.-nnc2)
      nnce2 =3.*nna1*(4.*nna1+3.*nna2*(1.-nnc5))*(1.-nnc2)
      nnce3 =6.*nna1*nna2
      nnce4 =12.*nna1*nna2*(1.-nnc2)+3.*nna2*nnb2*(1.-nnc3)
      nncr1 =nna1*(1.-3*nnc1)

      !----- Turbulent parameters from Helfand and Labraga (1988).  Another good ----------!
      hle1  = hlb1 - 6.* hla1
      hle2  = 12. * hla1 + hlb1 + 3. * hlb2
      hle3  = hlb1 * (1. - 3. * hlc1) - 6. * hla1
      hle4  = hlb1 * (1. - 3. * hlc1) + 12. * hla1 + 9. * hla2
      hle5  = hlb1 + 3. * hla1 + 3. * hlb2
      hlch  = hla2 * hle2 / hlb1
      hlcm  = hla1 * hle4 / (hla2 * hle5)
      hlrf1 = 1.
      hlrf2 = hle1 / hle2
      hlrf3 = hle1 / hle5
      hlrf4 = hle3 / hle4
      !------------------------------------------------------------------------------------!

      !-------------------------------------------------------------------------------------!
      !     New parameters based on Janjic (2001) constraint in L/q for unstable cases      !
      ! They correspond to the parts of C and D                                             !
      !-------------------------------------------------------------------------------------!
      nno1  =  9. * nna1 * nna2 * nna2 * (12. * nna1 * (1.-nnc2) + 3. * nnb2 * (1.-nnc2))
      nno2  = 18. * nna1 * nna1 * nna2                                                     &
            * (nnb2 * (1.-nnc3) - 3. * nna2 * (1.-nnc2) * (1.-nnc5))
      nno3  =  6. * nna1 * nna1
      nno4  =  3. * nna2 * (7. * nna1 * (1.-nnc2) + nnb2 * (1.-nnc3))
      nno5  = 27. * nna1 * nna2 * nna2 * nnb2 * (1.-nnc2) * (1.-nnc3)
      nno6  = 54. * nna1 * nna1 * nna2 * nnb2 * nnc1 * (1.-nnc3)
      nno7  = 18. * nna1 * nna1 * nnc1
      nno8  =  9. * nna1 * nna2 * (1-nnc2) + 3. * nna2 * nnb2 * (1.-nnc3)
      nnaeh =  9. * nna1 * nna2 * nna2 * (3. * (1.-nnc2) * (4.*nna1+nnb2) + nna2 * nnb1)
      nnaem =  3. * nna1 * nna2 * nnb1                                                     &
            * ( 3. * nna2 * (1.-nnc2) * (1.-nnc5) + 3. * nnb2 * nnc1 * (1.-nnc3)           &
              + 6. * nna1 * nnc1 * (3.-2.*nnc2) - nnb2 * (1.-nnc3) )                       &
            + 18. * nna1 * nna1 * nna2 * (nnb2 * (1.-nnc3)                                 &
                                         - 3. * nna2 * (1.-nnc2) * (1.-nnc5))
      nnreq = - nnaeh / nnaem

      !------------------------------------------------------------------------------------!
      !  nnmacheps is the machine precision                                                !
      !------------------------------------------------------------------------------------!
      nnmacheps = epsilon(1.)
      nnrsl     = (1.+nnmacheps) * ((nno5 + nno6 * nnreq)/(3. * nno1 + 3. * nno2 * nnreq))


      !------------------------------------------------------------------------------------!
      !    Print the banner with the turbulence constants.                                 !
      !------------------------------------------------------------------------------------!
      if (printinscreen) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Mellor-Yamada / Nakanishi-Niino / Janjic constants '
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNA1      =',nna1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNA2      =',nna2
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNB1      =',nnb1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNB2      =',nnb2
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNC1      =',nnc1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNC2      =',nnc2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNC3      =',nnc3
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNC4      =',nnc4
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNC5      =',nnc5
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNS1      =',nns1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNS2      =',nns2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNS3      =',nns3
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNS4      =',nns4
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNS5      =',nns5
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNS6      =',nns6
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNGAMA1   =',nngama1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNGAMA2   =',nngama2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNF1      =',nnf1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNF2      =',nnf2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRF1     =',nnrf1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRF2     =',nnrf2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRFC     =',nnrfc
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRI1     =',nnri1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRI2     =',nnri2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRI3     =',nnri3
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNCE1A    =',nnce1a
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNCE1B    =',nnce1b
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNCE2     =',nnce2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNCE3     =',nnce3
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNCE4     =',nnce4
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNCR1     =',nncr1
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO1      =',nno1
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO2      =',nno2
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO3      =',nno3
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO4      =',nno4
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO5      =',nno5
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO6      =',nno6
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO7      =',nno7
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNO8      =',nno8
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNAEH     =',nnaeh
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNAEM     =',nnaem
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNREQ     =',nnreq
         write (unit=*,fmt='(1(a,1x,f12.5,1x))')  'NNRSL     =',nnrsl
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(1(a,1x,es12.5,1x))') 'NNMACHEPS =',nnmacheps
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
      end if

      return
   end subroutine assign_turb_params
   !=======================================================================================!
   !=======================================================================================!

end module turb_coms
!==========================================================================================!
!==========================================================================================!
