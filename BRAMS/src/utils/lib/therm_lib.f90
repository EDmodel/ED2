!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

function rslf(p,t)
use rconstants, only: ep,t00
!     This function calculates the liquid saturation vapor mixing ratio as
!     a function of pressure and Kelvin temperature

implicit none
real esl,rslf,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01)
parameter (c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05)
parameter (c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13)

x=max(-80.,t-t00)

esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rslf=ep*esl/(p-esl)

return
end

!     ******************************************************************

real function rsif(p,t)

!     This function calculates the ice saturation vapor mixing ratio as a
!     function of pressure and Kelvin temperature
use rconstants, only: ep,t00
implicit none
real esi,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01)
parameter (c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05)
parameter (c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12)

x=max(-80.,t-t00)
esi=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rsif=ep*esi/(p-esi)

return
end

!     ******************************************************************

real function rslif(p,t)
use rconstants, only: t00
implicit none
real :: p,t
real, external :: rslf,rsif

!     This function calculates the saturation vapor mixing ratio, over
!     liquid or ice depending on temperature, as a function of pressure 
!     and Kelvin temperature

if (t >= t00) then
   rslif = rslf(p,t)
else
   rslif = rsif(p,t)
endif

return
end

!     ******************************************************************

real function eslf(t)

!     This function calculates the liquid saturation vapor pressure as a
!     function of Celcius temperature

implicit none
real :: x,t
real, parameter ::c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01
real, parameter ::c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05
real, parameter ::c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13

x=max(-80.,t)
eslf=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end

!     ******************************************************************

real function esif(t)

!     This function calculates the ice saturation vapor pressure as a
!     function of Celsius temperature

implicit none
real :: x,t
real, parameter ::c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01
real, parameter ::c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05
real, parameter ::c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12

x=max(-80.,t)
esif=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end

!     ******************************************************************

real function eslpf(t)

!     This function calculates the partial derivative of liquid saturation vapor
!     pressure with respect to temperature as a function of Celsius temperature

implicit none
real :: x,t
real, parameter ::d0= .4443216e+02,d1= .2861503e+01,d2= .7943347e-01
real, parameter ::d3= .1209650e-02,d4= .1036937e-04,d5= .4058663e-07
real, parameter ::d6=-.5805342e-10,d7=-.1159088e-11,d8=-.3189651e-14

x=max(-80.,t)
eslpf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
end

!     ******************************************************************

real function esipf(t)

!     This function calculates the partial derivative of ice saturation vapor
!     pressure with respect to temperature as a function of Celsius temperature

implicit none
real  :: x,t
real, parameter ::d0= .5036342e+02,d1= .3775758e+01,d2= .1269736e+00
real, parameter ::d3= .2503052e-02,d4= .3163761e-04,d5= .2623881e-06
real, parameter ::d6= .1392546e-08,d7= .4315126e-11,d8= .5961476e-14

x=max(-80.,t)
esipf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
end

!     ******************************************************************

!     This function calculates the partial derivative of liquid saturation vapor
!     mixing ratio with respect to temperature as a function of pressure and
!     Kelvin temperature

real function rslfp(p,t)
use rconstants, only: ep,t00

implicit none
real :: eslpf,x,t,p
real, parameter ::d0= .4443216e+02,d1= .2861503e+01,d2= .7943347e-01
real, parameter ::d3= .1209650e-02,d4= .1036937e-04,d5= .4058663e-07
real, parameter ::d6=-.5805342e-10,d7=-.1159088e-11,d8=-.3189651e-14

x=max(-80.,t-t00)
eslpf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
rslfp=ep*eslpf/(p-eslpf)

return
end

!     ******************************************************************

!     This function calculates the partial derivative of ice saturation vapor
!     mixing ratio with respect to temperature as a function of pressure and
!     Kelvin temperature

real function rsifp(p,t)
use rconstants, only: ep,t00
implicit none
real :: esipf,x,t,p
real, parameter ::d0= .5036342e+02,d1= .3775758e+01,d2= .1269736e+00
real, parameter ::d3= .2503052e-02,d4= .3163761e-04,d5= .2623881e-06
real, parameter ::d6= .1392546e-08,d7= .4315126e-11,d8= .5961476e-14

x=max(-80.,t-t00)
esipf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
rsifp=ep*esipf/(p-esipf)

return
end

!     ******************************************************************

subroutine mrsl(n1,p,t,rsl)

implicit none
integer :: n,n1
real :: rsl(n1),rslf,p(n1),t(n1)

do n=1,n1
   rsl(n)=rslf(p(n),t(n))
enddo

return
end

!     ******************************************************************

subroutine mrsi(n1,p,t,rsi)

implicit none
integer :: n,n1
real :: rsi(n1),rsif,p(n1),t(n1)

do n=1,n1
   rsi(n)=rsif(p(n),t(n))
enddo

return
end

!     ******************************************************************

subroutine thvtoth(nn,theta,rv,rtp)

implicit none
integer :: nn,k
real :: theta(nn),rv(nn),rtp(nn)

do k=1,nn
  theta(k)=theta(k)*(1.+rtp(k))/(1.+1.61*rv(k))
enddo

return
end

!     ******************************************************************

real function td(p,rs)
use rconstants, only: ep

implicit none
real :: rr,rs,es,esln,p

rr=rs+1e-8
es=p*rr/(ep+rr)
esln=log(es)
td=(35.86*esln-4947.2325)/(esln-23.6837)

return
end

!     ******************************************************************

real function rs(p,t)
use rconstants, only: ep,t00
implicit none
real :: p,t,es

es=610.78*exp(17.269*(t-t00)/(t-35.86))
rs=ep*es/(p-es)

return
end

!     ******************************************************************

subroutine thetae(p,t,rv,the)
use rconstants, only: ep,t00,rgas,alvl,cp,cpog,rocp,g
implicit none
real :: p,t,rv,the
real :: pit,tupo,ttd,dz,tupn,tmn
integer :: itter
real, external :: td
logical :: converged

pit=p
tupo=t
ttd=td(p,rv)
dz=cpog*(t-ttd)
if(dz > 0.) then
  converged=.false.
  itloop: do itter=1,50
     tupn=t-g/cp*dz
     tmn=(tupn+t)*.5*(1.+.61*rv)
     pit=p*exp(-g*dz/(rgas*tmn))
     if(abs(tupn-tupo).lt.0.001) then
       converged=.true.
       exit 
     end if
     ttd=td(pit,rv)
     tupo=tupn
     dz=dz+cpog*(tupn-ttd)
  end do itloop
  if (.not.converged) stop 'Theta_e diverged (therm_lib.f90, subroutine thetae)' 
end if
the=tupo*(1e5/pit)**rocp*exp(alvl*rv/(cp*tupo))

return
end subroutine thetae

!     ******************************************************************

real function  tw( rvp,thet,p)
use rconstants, only:rocp
implicit none
real :: rvp,thet,p
real :: press,rvap,piter,temper,x,aos
real, external :: os,tsa,tmr
integer :: id

!     abs is absolute value
!     all arguments and tw (kelvin)

press=p*1.e-2
rvap=rvp*1.e3
piter =  press
itloop: do id=1,10
   temper=thet*(piter*1.e-3)**rocp
   x  =  .02*( tmr(rvap,piter) - temper)
   if( abs(x).lt. 0.01  ) exit itloop
   piter = piter* ( 2.**(x)  )
end do itloop
temper=thet*(piter*1.e-3)**rocp

aos  =   os(temper,piter)
tw   =  tsa( aos,press)

return
end

!     ******************************************************************

real function virtt(t,rv)
   !---------------------------------------------------------------------------------------!
   ! This function computes the virtual temperature.                                       !
   ! Inputs: t  - temperature   [    K];                                                   !
   !         rv - mixing ratio  [kg/kg];                                                   !
   !---------------------------------------------------------------------------------------!
   use rconstants, only: ep
   implicit none
   real, intent(in) :: t
   real, intent(in) :: rv
   
   virtt = t * (rv + ep)/(ep*(rv+1.))
   return
end function virtt


!     ******************************************************************

real function   os(t,p)
use rconstants, only: rocp
implicit none
real :: t,p
real,external :: w

!     os and t (kelvin) , p (millibars )

os=t*((1000./p)**rocp)/(exp(-2.6518986*w(t,p)/t))

return
end

!     ******************************************************************

real function tsa(os,p)
use rconstants, only : rocp
implicit none
real :: os,p
real :: a,tq,d,x
real,external :: w
integer :: id

!     tsa and os(kelvin),p(millibars)
!     sign(a,b) rreplaces the algebretic sign of a with the sign of b

a  =  os
tq =  253.16
d  =  120
itloop: do id= 1,12
   d = d/2.
!     if the temperature difference,x, is small,exit this loop
   x=a*exp(-2.6518986*w(tq,p)/tq)-tq*((1000./p)**.286)
   if(abs(x).lt.0.01) exit itloop
   tq = tq + sign(d,x)
end do itloop
tsa=tq

return
end

!     ******************************************************************

real function  esat(t)
use rconstants, only: t00
implicit none
real :: t

!     esat(millibars),t(kelvin)

real :: tc

tc=t-t00
esat=6.1078*exp((17.2693882*tc)/(tc+237.3))

return
end

!     ******************************************************************

real function  w(t,p)
use rconstants, only: ep 
implicit none
real :: t,p,x
real,external :: esat

!     w(grams water vapor/kilogram dry air ), p(millibar )

if(t < 999.) then
  x  =  esat(t)
  w  =  ep*x/(p-x)
else
  w=0.0
end if
return
end

!     ******************************************************************

real function  tmr(w,p)
use rconstants, only: ep
implicit none
real :: w,p,x

!     tmr(kelvin),w(grams water vapor/kilogram dry air),p(millibar)
!     log10  15   log to the base  ten.

x =  log10(   w*p/(ep+ w)  )
tmr=10.**(.0498646455*x+2.4082965)-7.07475+38.9114*((10.**(  &
  .0915*x ) - 1.2035 )**2 )

return
end

!     ******************************************************************

subroutine the2t(the,p,th,t,r)
use rconstants, only : cp, alvl,rocp,p00i
implicit none
real :: the,p,th,t,r

real :: pi,to,tn
integer :: itter
real, external :: rs
logical :: converged


pi=(p*p00i)**rocp
to=the/exp(alvl*.012/(cp*295.))*pi
converged=.false.
itloop: do itter=1,50
   r=rs(p,to)
   th=the/exp(alvl*r/(cp*to))
   tn=th*pi
   if(abs(to-tn).lt.0.005) then
     converged=.true.
     exit itloop
   end if
   to=to+(tn-to)*.3
end do itloop
if (.not. converged) then 
  write(6,1) the,p,to,tn,th,r
  1 format(' stop in routine the2t '/' the,p,to,tn,th,r',6e15.6)
  stop 10
end if
t=tn

return
end

!     ******************************************************************

subroutine qtk(q,tempk,fracliq)
use rconstants, only: cliqi,cicei,allii,alli,t00
implicit none
real, intent(in)  :: q
real, intent(out) :: tempk, fracliq

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!       tempk    temperature [K]
!       fracliq  liquid fraction [dimensionless]
!     Local Constants:
!       4186     specific heat of liquid [J/(kg K)]
!       2093     specific heat of ice [J/(kg K)]
!       334000   latent heat of fusion [J/kg]
!       273.15   conversion from temp [C] to temp [K]

if (q .le. 0.) then
   fracliq = 0.
   tempk = q * cicei + t00
elseif (q .ge. alli) then
   fracliq = 1.
   tempk = q * cliqi + 193.36
else
   fracliq = q * allii
   tempk = t00
endif

return
end

!     ******************************************************************

subroutine qtc(q,tempc,fracliq)
use rconstants, only: cliqi,cicei,allii,alli,t00
implicit none
real :: q,tempc,fracliq
real,parameter :: r4186=1./4186.,r2093=1./2093.,r334000=1./334000.

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!        tempc    temperature [C]
!        fracliq  liquid fraction [dimensionless]
!     Local Constants:
!        4186     specific heat of liquid [J/(kg K)]
!        2093     specific heat of ice [J/(kg K)]
!        334000   latent heat of fusion [J/kg]
!        273.15   conversion from temp [C] to temp [K]

if (q .le. 0.) then
   fracliq = 0.
   tempc = q * cicei
elseif (q .ge. alli) then
   fracliq = 1.
   tempc = q * cliqi - 80.
else
   fracliq = q * allii
   tempc = 0.
endif

return
end

!     ******************************************************************

subroutine qwtk(qw,w,dryhcap,tempk,fracliq)
use rconstants, only: cliqi,cliq,cicei,cice,allii,alli,t00
implicit none
real :: qw,w,dryhcap,tempk,fracliq
real :: qwliq0
!     Inputs:
!        qw       internal energy [J/m^2] or [J/m^3]
!        w        mass [kg/m^2] or [kg/m^3]
!        dryhcap  heat capacity of nonwater part [J/(m^2 K)] or [J/(m^3 K)]
!     Outputs:
!        tempk    temperature [K]
!        fracliq  liquid fraction [dimensionless]
!     Local Constants:
!        4186     specific heat of liquid [J/(kg K)]
!        2093     specific heat of ice [J/(kg K)]
!        334000   latent heat of fusion [J/kg]
!        273.15   conversion from temp [C] to temp [K]

qwliq0 = w * alli
if (qw < 0.) then
   fracliq = 0.
   tempk = qw / (cice * w + dryhcap) + t00
elseif (qw >= qwliq0) then
   fracliq = 1.
   tempk = (qw - qwliq0) / (cliq * w + dryhcap) + t00
else
   fracliq = qw / qwliq0
   tempk = t00
endif
return
end

