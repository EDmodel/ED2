!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

Module ozone_const

  !       ozone.h:   especificação das constantes utilizadas para correção de k
  !                   para serem utilizados no módulo químico
  
  ! constantes referentes ao cálculo de correção de k em função de T, T=300K

  real, dimension(15) :: A, Ea, B

  data A /1.185E-10, 0.0, 2.04E-05,2.642E+03, 2.055e+02,0.0,3.229E+05, &
       3.068E+02, 4.991E+03, 2.055E+01, 6.1E+03,4.11E+03, 2.936E+04,&
       4.213E+03, 5.578E+02/


  data Ea /-1.05,  0.0,  0.0, 2.72, 4.91, 0.0, 0.0, 0.0, -0.54,&
       1.19, -1.0,-0.57,  0.0,  0.0, -1.57/

  data B  /-2.00, 0.0, -4.80, -1.00, -1.00,0.0,-1.00,-1.00,-1.00,-1.00,-1.00,&
       -1.00, -1.00,-1.00, -1.00/
    
  real, parameter ::        &
       Rcal   = 0.0019872,  &
       RJOULE = 8.314,      &
       M      = 1.0e+06,    &
       o2     = 2.09e+05,   &
       Tref   = 300.0,      &
       Mih2o  = 18.01
 
  ! unidades das constantes
  ! A=fator de frequência de choques (ppm-1min-1)
  ! R= constante dos gases (kcal/mol/K)
  ! B= coeficiente de correção adicionado a equação de 
  !    Arrhenius para correção da T (adimensional)
  ! Ea= energia de ativação (kca/mol)
  ! M=concentração de ar (ppm)
  ! o2=concentração do oxigênio constante (ppm)

end Module ozone_const
