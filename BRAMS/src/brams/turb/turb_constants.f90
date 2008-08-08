!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module turb_constants

! Nakanishi and Niino (2004) set of constants
real, parameter :: nna1=1.18, nna2=0.665, nnb1=24.0, nnb2=15.0, nnc1=0.137
real, parameter :: nnc2=0.70, nnc3=0.323,  nnc4=0.0,  nnc5=0.2
! Constants to define the Turbulent length scale
real, parameter :: nnq1=3.7, nnq2=2.7, nnq3=100., nnq4=0.2, nnq5=0.23, nnq6=5. 
! Constants defined
real :: nngama1, nngama2, nnf1, nnf2, nnrf1, nnrf2, nnri1, nnri2, nnri3, nnrfc
! Some constants to compute E1, E2, etc. from Nakanishi and Niino (2004) equation 20
real :: nnce1a,nnce1b, nnce2,nnce3,nnce4,nncr1
! Part of the constants defined in the Janjic (2001) report, adapted to Nakanishi/Niino
real :: nno1, nno2, nno3, nno4, nno5, nno6, nno7, nno8, nnaeh, nnaem, nnreq,nnrsl, nnmacheps


Contains

   subroutine assign_const_nakanishi()

   implicit none
   logical :: printinscreen
   printinscreen=.false.

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
!------------------------------------------------------------------------------------------!
! Parameters to compute E1,E2,E3,E4, and R1, from Nakanishi and Niino (2004) equation 20.  !
!------------------------------------------------------------------------------------------!
   nnce1a=6.*nna1*nna1
   nnce1b=9.*nna1*nna2*(1.-nnc2)
   nnce2 =3.*nna1*(4.*nna1+3.*nna2*(1.-nnc5))*(1.-nnc2)
   nnce3 =6.*nna1*nna2
   nnce4 =12.*nna1*nna2*(1.-nnc2)+3.*nna2*nnb2*(1.-nnc3)
   nncr1 =nna1*(1.-3*nnc1)
!------------------------------------------------------------------------------------------!
! New parameters based on Janjic (2001) constraint in L/q for unstable cases               !
! They correspond to the parts of C and D                                                  !
!------------------------------------------------------------------------------------------!
   nno1=9.*nna1*nna2*nna2*(12.*nna1*(1.-nnc2)+3.*nnb2*(1.-nnc2))
   nno2=18.*nna1*nna1*nna2*(nnb2*(1.-nnc3)-3.*nna2*(1.-nnc2)*(1.-nnc5))
   nno3=6.*nna1*nna1
   nno4=3.*nna2*(7.*nna1*(1.-nnc2)+nnb2*(1.-nnc3))
   nno5=27.*nna1*nna2*nna2*nnb2*(1.-nnc2)*(1.-nnc3)
   nno6=54.*nna1*nna1*nna2*nnb2*nnc1*(1.-nnc3)
   nno7=18.*nna1*nna1*nnc1
   nno8=9.*nna1*nna2*(1-nnc2)+3.*nna2*nnb2*(1.-nnc3)
   nnaeh=9.*nna1*nna2*nna2*(3.*(1.-nnc2)*(4.*nna1+nnb2)+nna2*nnb1)
   nnaem= 3.*nna1*nna2*nnb1*(3.*nna2*(1.-nnc2)*(1.-nnc5)+3.*nnb2*nnc1*(1.-nnc3)+           &
          6.*nna1*nnc1*(3.-2.*nnc2)-nnb2*(1.-nnc3))+                                       &
         18.*nna1*nna1*nna2*(nnb2*(1.-nnc3)-3.*nna2*(1.-nnc2)*(1.-nnc5))
   nnreq=-nnaeh/nnaem
!------------------------------------------------------------------------------------------!
! nnmacheps is the machine precision                                                       !
!------------------------------------------------------------------------------------------!
   nnmacheps=epsilon(1.)
   nnrsl=(1.+nnmacheps)*((nno5+nno6*nnreq)/(3.*nno1+3.*nno2*nnreq))
   if (printinscreen) then
     write (unit=*,fmt='(a)') '--------------- Nakanishi-Niino / Janjic constants ---------------'
     write (unit=*,fmt='(4(a,1x,f12.5,1x))') 'nna1=   ',nna1,    'nna2=   ',nna2,   'nnb1=   ',nnb1,  'nnb2=   ',nnb2
     write (unit=*,fmt='(4(a,1x,f12.5,1x))') 'nnc1=   ',nnc1,    'nnc2=   ',nnc2,   'nnc3=   ',nnc3,  'nnc5=   ',nnc5
     write (unit=*,fmt='(3(a,1x,f12.5,1x))') 'nnq1=   ',nnq1,    'nnq2=   ',nnq2,   'nnq3=   ',nnq3
     write (unit=*,fmt='(3(a,1x,f12.5,1x))') 'nnq4=   ',nnq4,    'nnq5=   ',nnq5,   'nnq6=   ',nnq6
     write (unit=*,fmt='(4(a,1x,f12.5,1x))') 'nngama1=',nngama1, 'nngama2=',nngama2,'nnf1=   ',nnf1,  'nnf2=   ',nnf2
     write (unit=*,fmt='(3(a,1x,f12.5,1x))') 'nnrf1=  ',nnrf1,   'nnrf2=  ',nnrf2,  'nnrfc=  ',nnrfc
     write (unit=*,fmt='(3(a,1x,f12.5,1x))') 'nnri1=  ',nnri1,   'nnri2=  ',nnri2,  'nnri3=  ',nnri3
     write (unit=*,fmt='(3(a,1x,f12.5,1x))') 'nnce1a= ',nnce1a,  'nnce1b= ',nnce1b, 'nnce2=  ',nnce2
     write (unit=*,fmt='(3(a,1x,f12.5,1x))') 'nnce3=  ',nnce3,   'nnce4=  ',nnce4,  'nncr1=  ',nncr1
     write (unit=*,fmt='(4(a,1x,f12.5,1x))') 'nno1=   ',nno1,    'nno2=   ',nno2,   'nno3=   ',nno3,  'nno4=   ',nno4
     write (unit=*,fmt='(4(a,1x,f12.5,1x))') 'nno5=   ',nno5,    'nno6=   ',nno6,   'nno7=   ',nno7,  'nno8=   ',nno8
     write (unit=*,fmt='(4(a,1x,f12.5,1x))') 'nnaeh=  ',nnaeh,   'nnaem=  ',nnaem,  'nnreq=  ',nnreq, 'nnrsl=  ',nnrsl
     write (unit=*,fmt='(a,1x,es12.5)')   'nnmacheps=',nnmacheps
     write (unit=*,fmt='(a)') '------------------------------------------------------------------'
   end if
   
   end subroutine assign_const_nakanishi
End Module turb_constants
