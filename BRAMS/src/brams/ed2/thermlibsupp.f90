real function rhovsil(tc)
  
  implicit none
  
  real, intent(in) :: tc
  real, external :: rhovsl,rhovsi
  
  !     This function calculates the density of water vapor at saturation,
  !     over liquid or ice depending on temperature, as a function of
  !     Celsius temperature
  
  if (tc >= 0.) then
     rhovsil = rhovsl(tc)
  else
     rhovsil = rhovsi(tc)
  endif
  
  return
end function rhovsil
!============================================================================

real function rhovsi(tc)

  !     This function calculates the density of water vapor at saturation
  !     over ice as a function of Celsius temperature
  
  use consts_coms, only: t00,rvap
  implicit none
  real, intent(in) :: tc
  
  real, parameter :: c0 = .6114327e+03 ,c1 = .5027041e+02 ,c2 = .1875982e+01
  real, parameter :: c3 = .4158303e-01 ,c4 = .5992408e-03 ,c5 = .5743775e-05
  real, parameter :: c6 = .3566847e-07 ,c7 = .1306802e-09 ,c8 = .2152144e-12
  real :: esi,x
  
  x = max(-80.,tc)
  esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  rhovsi = esi / (rvap * (tc + t00))
  
  return
end function rhovsi
