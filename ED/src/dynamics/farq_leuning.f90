subroutine lphysiol_full(T_L,  &
     e_A,  &
     C_A,  &
     PAR,  &
     rb,  &
     adens,  &
     A_open,  &
     A_cl,  &
     rsw_open,  &
     rsw_cl,  &
     pft,  &
     prss,  &
     leaf_resp,  &
     green_leaf_factor,  &
     leaf_aging_factor, &
     old_st_data)

  use c34constants
  use pft_coms, only: D0, cuticular_cond, dark_respiration_factor,   &
       stomatal_slope, quantum_efficiency, photosyn_pathway, Vm0, Vm_low_temp
  use physiology_coms, only: istoma_scheme
  use therm_lib, only : rslif
  use consts_coms, only : t00,mmdov
  implicit none

  real, intent(in) :: T_L
  real, intent(in) :: e_A
  real, intent(in) :: C_A
  real, intent(in) :: PAR
  real, intent(in) :: rb
  real, intent(in) :: adens
  real, intent(out) :: A_open
  real, intent(out) :: A_cl
  real, intent(out) :: rsw_open
  real, intent(out) :: rsw_cl
  real, intent(in) :: prss
  real, intent(out) :: leaf_resp
  real, intent(in) :: green_leaf_factor
  real, intent(in) :: leaf_aging_factor
  integer, intent(in) :: pft
  type(stoma_data), intent(inout) :: old_st_data

  type(farqdata) :: gsdata
  type(metdat) :: met
  type(solution) :: sol
  integer :: recalc
  type(glim) :: apar
  integer :: ilimit

  real :: co2cp,arrhenius
  real :: gsw_update
  integer :: success

  ! load physiological parameters into structure
  gsdata%D0 = D0(pft)
  gsdata%b = cuticular_cond(pft)
  gsdata%gamma = dark_respiration_factor(pft)
  gsdata%m = stomatal_slope(pft)
  gsdata%alpha = quantum_efficiency(pft)
  ! Load met into structure
  met%tl = T_L
  met%ea = e_A
  met%ca = C_A
  met%par = PAR*(1.0e6)
  met%gbc = adens / (rb*4.06e-8)
  met%gbci = 1.0/met%gbc
  met%el = mmdov * rslif(prss, met%tl + t00)
  met%compp = co2cp(met%tl)
  met%gbw = 1.4*met%gbc
  met%eta = 1.0 + (met%el-met%ea)/gsdata%d0

  ! Set up how we look for the solution
  sol%eps = 3.0e-8
  sol%ninterval = 6

  ! Prepare derived terms for both exact and approximate solutions
  call prep_lphys_solution(photosyn_pathway(pft), Vm0(pft), met,   &
       Vm_low_temp(pft), leaf_aging_factor, green_leaf_factor, leaf_resp,   &
       gsdata, apar)

  ! Decide whether to do the exact solution or the approximation
  recalc = 1
  if(istoma_scheme == 1)then
     if(old_st_data%recalc == 0)recalc = 0
  endif

  if(recalc == 1)call exact_lphys_solution(photosyn_pathway(pft), met,  &
       apar, gsdata, sol, ilimit)

  if(istoma_scheme == 1 .and. recalc == 1)then
     call store_exact_lphys_solution(old_st_data, met, prss,   &
          leaf_aging_factor, green_leaf_factor, sol, ilimit, gsdata, apar,  &
          photosyn_pathway(pft), Vm0(pft), Vm_low_temp(pft))
  endif

  if(recalc == 1)then

     call fill_lphys_sol_exact(A_open, rsw_open, A_cl, rsw_cl, sol, adens)
     if(istoma_scheme == 1)old_st_data%recalc = 0

  else

     call fill_lphys_sol_approx(gsdata, met, apar, old_st_data, sol,   &
          A_cl, rsw_cl, adens, rsw_open, A_open, photosyn_pathway(pft), prss)

  endif

  return
end subroutine lphysiol_full

!================================================================

subroutine c3solver(met,apar,gsdata,sol,ilimit)
  use c34constants
  implicit none

  
  type(farqdata) :: gsdata
  type(metdat) :: met
  type(glim) :: apar
  type(solution) :: sol
  integer :: success_flag
  integer :: ilimit

  ! Solve par case
  ilimit = 1
  call setapar_c3(gsdata,met,apar,1)
  call solve_closed_case_c3(gsdata,met,apar,sol,1)
  ! Return if nighttime
  if(met%par < 1.0e-3)then
     call closed2open(sol,1)
     ilimit = -1
     return
  endif
  success_flag = 1
  call solve_open_case_c3(gsdata,met,apar,sol,1,success_flag)
  if(success_flag == 0)then
     call closed2open(sol,1)
     ilimit = -1
     return
  endif

  ! Solve the vm case
  call setapar_c3(gsdata,met,apar,2)
  call solve_closed_case_c3(gsdata,met,apar,sol,2)
  call solve_open_case_c3(gsdata,met,apar,sol,2,success_flag)
  if(success_flag == 0)then
     call closed2open(sol,1)
     ilimit = -1
     return
  endif

  if(sol%a(1,2) < sol%a(1,1))then
     sol%gsw(1,1) = sol%gsw(1,2)
     sol%es(1,1) = sol%es(1,2)
     sol%ci(1,1) = sol%ci(1,2)
     sol%cs(1,1) = sol%cs(1,2)
     sol%a(1,1) = sol%a(1,2)
  endif
  if(sol%a(2,2) < sol%a(2,1))then
     sol%gsw(2,1) = sol%gsw(2,2)
     sol%es(2,1) = sol%es(2,2)
     sol%ci(2,1) = sol%ci(2,2)
     sol%cs(2,1) = sol%cs(2,2)
     sol%a(2,1) = sol%a(2,2)
     ilimit = 2
  endif

  if(sol%cs(2,1) > 1.25e7)then
     success_flag = 0
     
!     print*,'Stomatal conductance dangerously close to upper limit.  Stopping.'
!     print*,met%ea,met%ca,met%rn,met%tl,met%par,met%gbc,met%gbw,met%ta  &
!          ,met%el,met%compp,met%eta

  endif

  return

end subroutine c3solver
!================================================================

subroutine c4solver(met,apar,gsdata,sol,ilimit)
  use c34constants
  implicit none
  

  type(farqdata) :: gsdata
  type(metdat) :: met
  type(glim) :: apar
  type(solution) :: sol
  integer :: success_flag
  integer :: ilimit

  if(apar%vm > gsdata%alpha*met%par)then
     ! Solve par case
     ilimit = 1
     call setapar_c4(gsdata,met,apar,1)
     call solve_closed_case_c4(gsdata,met,apar,sol,1)
     if(met%par < 1.0e-3)then
        call closed2open(sol,1)
        ilimit = -1
        return ! nighttime; return
     endif 
     success_flag = 1
     call solve_open_case_c4(gsdata,met,apar,sol,1,success_flag)
     if(success_flag == 0)then
        call closed2open(sol,1)
        ilimit = -1
        return
     endif
  else
     ! Solve the vm case
     ilimit = 2
     call setapar_c4(gsdata,met,apar,2)
     ! yes, these should be ones below.  
     call solve_closed_case_c4(gsdata,met,apar,sol,1)
     success_flag = 1
     call solve_open_case_c4(gsdata,met,apar,sol,1,success_flag)
     if(success_flag == 0)then
        call closed2open(sol,1)
        ilimit = -1
        return
     endif
  endif
  ! Solve the third case
  call setapar_c4(gsdata,met,apar,3)
  call solve_closed_case_c4(gsdata,met,apar,sol,2)
  call solve_open_case_c4(gsdata,met,apar,sol,2,success_flag)
  if(success_flag == 0)then
     call closed2open(sol,1)
     ilimit = -1
     return
  endif

  if(sol%a(1,2) < sol%a(1,1))then
     sol%gsw(1,1) = sol%gsw(1,2)
     sol%es(1,1) = sol%es(1,2)
     sol%ci(1,1) = sol%ci(1,2)
     sol%cs(1,1) = sol%cs(1,2)
     sol%a(1,1) = sol%a(1,2)
  endif
  if(sol%a(2,2) < sol%a(2,1))then
     sol%gsw(2,1) = sol%gsw(2,2)
     sol%es(2,1) = sol%es(2,2)
     sol%ci(2,1) = sol%ci(2,2)
     sol%cs(2,1) = sol%cs(2,2)
     sol%a(2,1) = sol%a(2,2)
     ilimit = 3
  endif

  if(sol%cs(2,1) > 1.25e7)then
     print*,'Stomatal conductance dangerously close to upper limit.  Stopping.'
     success_flag = 0
     print*,met%ea,met%ca,met%rn,met%tl,met%par,met%gbc,met%gbw,met%ta  &
          ,met%el,met%compp,met%eta
  endif

  return

end subroutine c4solver

!=====================================
subroutine setapar_c3(gsdata,met,apar,i)
  use c34constants
  implicit none
  

  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  integer :: i

  if(i == 1)then
     apar%rho = gsdata%alpha*met%par
     apar%sigma = -gsdata%alpha*met%par*met%compp
     apar%tau = 2.0*met%compp
  elseif(i == 2)then
     apar%rho = apar%vm
     apar%sigma = -apar%vm*met%compp
     apar%tau = apar%k1*(1.0+apar%k2)
  endif

  return
end subroutine setapar_c3

!=====================================
subroutine setapar_c4(gsdata,met,apar,i)
  use c34constants
  implicit none
  

  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  integer :: i

  if(i == 1)then
     apar%rho = 0.0
     apar%sigma = gsdata%alpha*met%par-gsdata%gamma*apar%vm
  elseif(i == 2)then
     apar%rho = 0.0
     apar%sigma = apar%vm*(1.0-gsdata%gamma)
  elseif(i == 3)then
     apar%rho = 18000.0*apar%vm
     apar%sigma = -gsdata%gamma*apar%vm
  endif

  return
end subroutine setapar_c4

!================================================
real function aflux_c3(apar,gsw,ci)
  use c34constants
  implicit none
  
  type(glim) :: apar
  real :: gsw,ci

  aflux_c3 = (apar%rho*ci+apar%sigma)/(ci+apar%tau)+apar%nu

  return
end function aflux_c3

!================================================
real function aflux_c4(apar,met,gsw)
  use c34constants
  implicit none
  
  type(glim) :: apar
  type(metdat) :: met
  real :: gsw

  aflux_c4 = apar%rho * gsw * (-gsw*apar%sigma   &
       + met%gbc*met%ca*(1.6*apar%rho+gsw))  &
       /((1.6*apar%rho+gsw)*(apar%rho*gsw+met%gbc*(1.6*apar%rho+gsw))) &
       +apar%sigma*gsw/(1.6*apar%rho+gsw)

  return
end function aflux_c4

!================================================
real function csc_c4(apar,met,gsw)
  use c34constants
  implicit none
  
  type(metdat) :: met
  type(glim) :: apar
  real :: gsw

  csc_c4 = (-gsw*apar%sigma   &
       + met%gbc*met%ca*(1.6*apar%rho+gsw))  &
       /(apar%rho*gsw+met%gbc*(1.6*apar%rho+gsw))


  return
end function csc_c4

!================================================
real function csc_c3(met,a)
  use c34constants
  implicit none
  
  type(metdat) :: met
  real :: a

  csc_c3 = met%ca-a*met%gbci


  return
end function csc_c3

!================================================
real function residual_c3(gsdata,met,apar,x)
  use c34constants
  implicit none
  

  type(farqdata) :: gsdata
  type(metdat) :: met
  type(glim) :: apar
  real :: x,ci,a,cs
  integer :: success
  real, external :: quad4ci
  real, external :: aflux_c3
  real, external :: csc_c3

  ci = quad4ci(gsdata,met,apar,x,success)
  if(success == 0)then
     residual_c3 = 9.9e9
     return
  endif
  a = aflux_c3(apar,x,ci)
  cs = csc_c3(met,a)

  residual_c3 = (cs-met%compp)*x**2   &
       + ((met%gbw*met%eta-gsdata%b)*(cs-met%compp)-gsdata%m*a)*x  &
       -gsdata%b*met%eta*met%gbw*(cs-met%compp)-gsdata%m*a*met%gbw

  return
end function residual_c3

!================================================
real function residual_c4(gsdata,met,apar,x)
  use c34constants
  implicit none
  

  type(farqdata) :: gsdata
  type(metdat) :: met
  type(glim) :: apar
  real :: a,cs,x
  real, external :: aflux_c4
  real, external :: csc_c4

  a = aflux_c4(apar,met,x)
  cs = csc_c4(apar,met,x)

  residual_c4 = (cs-met%compp)*x**2   &
       + ((met%gbw*met%eta-gsdata%b)*(cs-met%compp)-gsdata%m*a)*x  &
       -gsdata%b*met%eta*met%gbw*(cs-met%compp)-gsdata%m*a*met%gbw


  return
end function residual_c4

!===============================
subroutine zbrak_c3(gsdata,met,apar,x1,x2,n,xb1,xb2,nb)
  use c34constants
  implicit none
  

  real :: x1,x2
  integer :: n,nb,nbb,i
  real :: x,dx,fp,fc
  real, dimension(nb) :: xb1,xb2
  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  real, external :: residual_c3

  nbb = nb
  nb = 0
  x = x1
  dx = (x2-x1)/n
  fp = residual_c3(gsdata,met,apar,x)
  do i=1,n
     x = x + dx
     fc = residual_c3(gsdata,met,apar,x)
     if(fc*fp < 0.0)then
        nb = nb + 1
        xb1(nb) = x-dx
        xb2(nb) = x
     endif
     fp = fc
     if(nbb == nb)return
  enddo

  return
end subroutine zbrak_c3

!-----------------------------------
real function zbrent_c3(gsdata,met,apar,x1,x2,tol,success_flag)
  use c34constants
  implicit none
  

  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  integer, parameter :: itmax = 100
  real, parameter :: eps = 3.0e-8
  real :: x1,x2,tol,a,b,fa,fb,fc,c,d,e,tol1,s,q,p,r,xm
  integer :: iter,it,success_flag
  real, external :: residual_c3

  a = x1
  b = x2

  fa = residual_c3(gsdata,met,apar,a)
  fb = residual_c3(gsdata,met,apar,b)
  if(fb*fa > 0.0)then
     print*,'Root must be bracketed for ZBRENT_C3.'
     print*,fa,fb
     success_flag = 0
     zbrent_c3=0.0
     return
  endif
  fc=fb
  do iter = 1,itmax
     if(fb*fc > 0.0)then
        c=a
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc) < abs(fb))then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.0*eps*abs(b)+0.5*tol
     xm = 0.5*(c-b)
     if(abs(xm) <= tol1.or.fb == 0.0)then
        zbrent_c3=b
        return
     endif
     if(abs(e) >= tol1 .and. abs(fa) > abs(fb))then
        s=fb/fa
        if(a == c)then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
           q=(q-1.0)*(r-1.0)*(s-1.0)
        endif
        if(p > 0.0) q = -q
        p=abs(p)
        if(2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q)))then
           e=d
           d=p/q
        else
           d=xm
           e=d
        endif
     else
        d=xm
        e=d
     endif
     a=b
     fa=fb
     if(abs(d) > tol1)then
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb = residual_c3(gsdata,met,apar,b)
  enddo
  write (unit=*,fmt='(a)') 'ZBRENT_C3 exceeding maximum iterations.'
  zbrent_c3=b
  return
end function zbrent_c3

!===============================
subroutine zbrak_c4(gsdata,met,apar,x1,x2,n,xb1,xb2,nb)
  use c34constants
  implicit none
  

  real :: x1,x2
  integer :: n,nb,nbb,i
  real :: x,dx,fp,fc
  real, dimension(nb) :: xb1,xb2
  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  real, external :: residual_c4

  nbb = nb
  nb = 0
  x = x1
  dx = (x2-x1)/n
  fp = residual_c4(gsdata,met,apar,x)
  do i=1,n
     x = x + dx
     fc = residual_c4(gsdata,met,apar,x)
     if(fc*fp < 0.0)then
        nb = nb + 1
        xb1(nb) = x-dx
        xb2(nb) = x
     endif
     fp = fc
     if(nbb == nb)return
  enddo

  return
end subroutine zbrak_c4

!-----------------------------------
real function zbrent_c4(gsdata,met,apar,x1,x2,tol,success_flag)
  use c34constants
  implicit none
  

  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  integer, parameter :: itmax = 100
  real, parameter :: eps = 3.0e-8
  real :: x1,x2,tol,a,b,fa,fb,fc,c,d,e,tol1,s,q,p,r,xm
  integer :: iter,it,success_flag
  real, external :: residual_c4

  a = x1
  b = x2

  fa = residual_c4(gsdata,met,apar,a)
  fb = residual_c4(gsdata,met,apar,b)
  if(fb*fa > 0.0)then
     print*,'Root must be bracketed for ZBRENT_C4.'
     print*,fa,fb
     success_flag = 0
     zbrent_c4 = 0.0
     print*,met%ea,met%ca,met%rn,met%tl,met%par,met%gbc,met%gbw,met%ta  &
          ,met%el,met%compp,met%eta
     return
  endif
  fc=fb
  do iter = 1,itmax
     if(fb*fc > 0.0)then
        c=a
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc) < abs(fb))then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.0*eps*abs(b)+0.5*tol
     xm = 0.5*(c-b)
     if(abs(xm) <= tol1.or.fb == 0.0)then
        zbrent_c4=b
        return
     endif
     if(abs(e) >= tol1 .and. abs(fa) > abs(fb))then
        s=fb/fa
        if(a == c)then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
           q=(q-1.0)*(r-1.0)*(s-1.0)
        endif
        if(p > 0.0) q = -q
        p=abs(p)
        if(2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q)))then
           e=d
           d=p/q
        else
           d=xm
           e=d
        endif
     else
        d=xm
        e=d
     endif
     a=b
     fa=fb
     if(abs(d) > tol1)then
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb = residual_c4(gsdata,met,apar,b)
  enddo
  pause 'ZBRENT_C4 exceeding maximum iterations.'
  zbrent_c4=b
  return
end function zbrent_c4

!=====================================================
subroutine solve_closed_case_c3(gsdata,met,apar,sol,ilimit)
  use c34constants
  implicit none
  

  integer :: ilimit
  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  type(solution) :: sol
  real :: b,c,x1,x2,q
  real, dimension(2) :: ci,cs,a
  integer :: j
  

  ! Do not allow assimilation

  sol%gsw(1,ilimit) = gsdata%b
  sol%es(1,ilimit) = (met%ea*met%gbw+gsdata%b*met%el)/(gsdata%b+met%gbw)
  sol%a(1,ilimit) = -gsdata%gamma*apar%vm
  sol%cs(1,ilimit) = met%ca - sol%a(1,ilimit)/met%gbc
  sol%ci(1,ilimit) = sol%cs(1,ilimit) - sol%a(1,ilimit) * 1.6 / gsdata%b
  return

  ! Allow assimilation


  b=apar%tau-met%ca+(apar%rho+apar%nu)*(gsdata%b+1.6*met%gbc)  &
       /(gsdata%b*met%gbc)
  c=(apar%sigma+apar%nu*apar%tau)*(gsdata%b+1.6*met%gbc)/(gsdata%b*met%gbc)  &
       -apar%tau*met%ca
  q=-0.5*b*(1.0+sqrt(1.0-4.0*c/b**2))
  ci(1)=q
  ci(2)=c/q
  if(abs(ci(1)-met%ca) < abs(ci(2)-met%ca))then
     j=1
  else
     j=2
  endif

  sol%gsw(1,ilimit) = gsdata%b
  sol%es(1,ilimit) = (met%ea*met%gbw+gsdata%b*met%el)/(gsdata%b+met%gbw)
  sol%ci(1,ilimit) = ci(j)
  sol%cs(1,ilimit) = (gsdata%b*ci(j)+1.6*met%gbc*met%ca)/(gsdata%b+1.6*met%gbc)
  sol%a(1,ilimit) = met%gbc*(met%ca-sol%cs(1,ilimit))

  return
end subroutine solve_closed_case_c3
!=====================================================
subroutine solve_closed_case_c4(gsdata,met,apar,sol,ilimit)
  use c34constants
  implicit none
  

  integer :: ilimit
  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  type(solution) :: sol

  sol%gsw(1,ilimit) = gsdata%b
  sol%es(1,ilimit) = (met%ea*met%gbw+gsdata%b*met%el)/(gsdata%b+met%gbw)
  sol%ci(1,ilimit) = (gsdata%b*met%gbc*met%ca-apar%sigma  &
       *(gsdata%b+1.6*met%gbc)) &
       / (apar%rho*(gsdata%b+1.6*met%gbc)+gsdata%b*met%gbc)
  sol%cs(1,ilimit) = (gsdata%b*sol%ci(1,ilimit)+1.6*met%gbc*met%ca)  &
       /(gsdata%b+1.6*met%gbc)
  sol%a(1,ilimit) = apar%sigma + apar%rho * sol%ci(1,ilimit)

  return
end subroutine solve_closed_case_c4

!=========================================================

subroutine solve_open_case_c3(gsdata,met,apar,sol,ilimit,success_flag)
  use c34constants
  implicit none

  
  integer :: ilimit,success_flag
  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  type(solution) :: sol

  integer, parameter :: maxroots=5
  integer :: nroot,isol,success
  real :: gswmin,gswmax
  real, dimension(maxroots) :: xb1,xb2
  real, external :: aflux_c3
  real, external :: quad4ci
  real, external :: zbrent_c3
  real, external :: csc_c3

  gswmin = gsdata%b
  gswmax = 1.3e7
  nroot = maxroots

  call zbrak_c3(gsdata,met,apar,gswmin,gswmax,sol%ninterval,xb1,xb2,nroot)
  if(nroot == 0)then
     ! No open case solution.  Values revert to those from closed case.
!     call closed2open(sol,ilimit)
     success_flag = 0
  else
     ! We did find a solution
     do isol=1,nroot
        sol%gsw(2,ilimit) = zbrent_c3(gsdata,met,apar  &
             ,xb1(isol),xb2(isol),sol%eps,success_flag)
        if(success_flag == 0)return
        sol%ci(2,ilimit) = quad4ci(gsdata,met,apar,sol%gsw(2,ilimit),success)
        if(success /= 0)then
           sol%a(2,ilimit)= aflux_c3(apar,sol%gsw(2,ilimit),sol%ci(2,ilimit))
           sol%cs(2,ilimit) = csc_c3(met,sol%a(2,ilimit))
           sol%es(2,ilimit) = (met%ea*met%gbw+sol%gsw(2,ilimit)*met%el)  &
                /(sol%gsw(2,ilimit)+met%gbw)
        elseif(nroot == 1)then
 !          call closed2open(sol,ilimit)
           success_flag = 0
       endif

     enddo
  endif

  return
end subroutine solve_open_case_c3

!=====================================================

subroutine solve_open_case_c4(gsdata,met,apar,sol,ilimit,success_flag)
  use c34constants
  implicit none
  

  integer :: ilimit
  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  type(solution) :: sol

  integer, parameter :: maxroots=5
  integer :: nroot,isol
  real :: gswmin,gswmax
  real, dimension(maxroots) :: xb1,xb2
  integer :: success_flag
  real, external :: zbrent_c4

  gswmin = gsdata%b
  gswmax = 1.3e7
  nroot = maxroots

  call zbrak_c4(gsdata,met,apar,gswmin,gswmax,sol%ninterval,xb1,xb2,nroot)
  if(nroot == 0)then
     ! No open case solution.  Values revert to those from closed case.
!     call closed2open(sol,ilimit)
     success_flag = 0
  else
     ! We did find a solution
     do isol=1,nroot
        sol%gsw(2,ilimit) = zbrent_c4(gsdata,met,apar  &
             ,xb1(isol),xb2(isol),sol%eps,success_flag)
        if(success_flag == 0)return
        sol%cs(2,ilimit) = (-apar%sigma*sol%gsw(2,ilimit)  &
             +met%gbc*met%ca*(1.6*apar%rho+sol%gsw(2,ilimit)))  &
             /(apar%rho*sol%gsw(2,ilimit)  &
             +met%gbc*(1.6*apar%rho+sol%gsw(2,ilimit)))
        sol%a(2,ilimit) = met%gbc*(met%ca-sol%cs(2,ilimit))
        sol%ci(2,ilimit) = sol%cs(2,ilimit)  &
             -1.6*sol%a(2,ilimit)/sol%gsw(2,ilimit)
        sol%es(2,ilimit) = (met%ea*met%gbw+sol%gsw(2,ilimit)*met%el)  &
             /(sol%gsw(2,ilimit)+met%gbw)

     enddo
  endif

  return
end subroutine solve_open_case_c4

!=====================================================
subroutine closed2open(sol,ilimit)
  use c34constants
  implicit none
  

  type(solution) :: sol
  integer :: ilimit

  sol%gsw(2,ilimit) = sol%gsw(1,ilimit)
  sol%es(2,ilimit) = sol%es(1,ilimit)
  sol%ci(2,ilimit) = sol%ci(1,ilimit)
  sol%cs(2,ilimit) = sol%cs(1,ilimit)
  sol%a(2,ilimit) = sol%a(1,ilimit)

  return
end subroutine closed2open
!=====================================================
real function quad4ci(gsdata,met,apar,x,success)
  use c34constants
  implicit none
  

  type(glim) :: apar
  type(farqdata) :: gsdata
  type(metdat) :: met
  real :: b,c,x1,x2,q,x,ciout,q_1
  real, dimension(2) :: ci
  integer :: j,success
  integer, dimension(2) :: sol_flag  
  integer :: isol
  logical,external :: isnan_ext

  b = (apar%rho + apar%nu) * (x + 1.6 * met%gbc) / (x * met%gbc) -   &
       met%ca + apar%tau
  c = (apar%sigma + apar%nu * apar%tau) * (x + 1.6 * met%gbc) /   &
       (x * met%gbc) - apar%tau * met%ca
  if (b == 0.) then
     q = tiny(0.)
  else
     q = -0.5 * b * (1.0 + sqrt(1.0 - 4.0 * c / b**2))
  endif
  ci(1) = q
  ci(2) = c / q
  success = 1

  ! Test to see if the solutions are greater than zero.
  sol_flag(1:2) = 1
  do isol = 1,2
     if(ci(isol) <= 0.0 .or. ci(isol) /= ci(isol))sol_flag(isol) = 0
  enddo

  if(sol_flag(1) == 0 .and. sol_flag(2) /= 0)then
     quad4ci = ci(2)
  elseif(sol_flag(1) /= 0 .and. sol_flag(2) == 0)then
     quad4ci = ci(1)
  elseif(sol_flag(1) /= 0 .and. sol_flag(2) /= 0)then
     if(abs(ci(1)-met%ca) < abs(ci(2)-met%ca))then
        quad4ci = ci(1)
     else
        quad4ci = ci(2)
     endif
  else
     quad4ci = met%ca
     success = 0
  endif

  return
end function quad4ci

!=============================================

logical function isnan_ext(rv)
  implicit none
  real :: rv
  if (rv == rv) then
     isnan_ext=.false.
  else
     isnan_ext=.true.
  endif
  return
end function isnan_Ext

!==============================================

real function co2cp(T)
  implicit none
  real :: arrhenius,T
  co2cp = arrhenius(T,2.12e-5,5000.0)
  return
end function co2cp

!=================================================

real function arrhenius(T,c1,c2)
  use consts_coms, only: t00
  implicit none
  real :: T,c1,c2
  real(kind=8) :: arr8
  arr8 = dble(c1) * dexp( dble(c2)*(dble(1.)/dble(288.15)-dble(1.0)/dble(T+t00)))
  arrhenius = sngl(arr8)
  return
end function arrhenius

!=========================================
subroutine testsolution(gsdata,met,apar,x)

  use c34constants

  implicit none

  type(farqdata) :: gsdata
  type(metdat) :: met
  type(glim) :: apar
  real :: tl,cs,psi,ci,eta,vm,gamma,arrhenius,co2cp,a,aflux,csc
  real :: x,es

  eta = 1.0 + (met%el-met%ea)/gsdata%d0
  gamma = co2cp(met%tl)

  a = apar%rho * x * (-x*apar%sigma   &
       + met%gbc*met%ca*(1.6*apar%rho+x))  &
       /((1.6*apar%rho+x)*(apar%rho*x+met%gbc*(1.6*apar%rho+x))) &
       +apar%sigma*x/(1.6*apar%rho+x)
  cs = (-x*apar%sigma   &
       + met%gbc*met%ca*(1.6*apar%rho+x))  &
       /(apar%rho*x+met%gbc*(1.6*apar%rho+x))
  ci = (x*cs/1.6-apar%sigma)/(apar%rho+x/1.6)
  es = (met%ea*met%gbw+x*met%el)/(x+met%gbw)

  return
end subroutine testsolution

!===============================================================

subroutine prep_lphys_solution(photosyn_pathway, Vm0, met, Vm_low_temp,  &
     leaf_aging_factor, green_leaf_factor, leaf_resp, gsdata, apar)

  use c34constants

  implicit none

  

  integer, intent(in) :: photosyn_pathway
  real, intent(in) :: Vm0
  real, intent(in) :: Vm_low_temp
  type(metdat), intent(in) :: met
  real, intent(in) :: leaf_aging_factor
  real, intent(in) :: green_leaf_factor
  real, intent(out) :: leaf_resp
  type(farqdata), intent(in) :: gsdata
  type(glim), intent(inout) :: apar
  real(kind=8) :: vmdble
  real, external :: arrhenius

  if(photosyn_pathway == 3)then

     ! C3 parameters
     vmdble = dble(Vm0) * dble(arrhenius(met%tl, 1.0, 3000.0)) / &
              ( (dble(1.0) + dexp(dble(0.4)*dble(Vm_low_temp - met%tl) )) * &
                (dble(1.0) + dexp(dble(0.4)*dble(met%tl - 45.0 ))) )
     apar%vm = sngl(vmdble)

     ! Adjust Vm according to the aging factor.
     if(leaf_aging_factor > 0.01 .and. green_leaf_factor > 0.0001)then
        apar%vm = apar%vm * leaf_aging_factor / green_leaf_factor
     endif

     ! Compute leaf respiration and other constants.
     leaf_resp = apar%vm * gsdata%gamma
     apar%nu = -apar%vm * gsdata%gamma
     apar%k1 = arrhenius(met%tl, 1.5e-4, 6000.0)
     apar%k2 = arrhenius(met%tl, 0.836, -1400.0)

  else

     ! C4 parameters

     apar%vm = Vm0 * arrhenius(met%tl, 1.0, 3000.0) / (  &
          (1.0 + exp(0.4 * (5.0    - met%tl))) *  &
          (1.0 + exp(0.4 * (met%tl - 100.0 ))) ) 

     leaf_resp = apar%vm * gsdata%gamma

  endif

  return
end subroutine prep_lphys_solution

!===================================================================

subroutine exact_lphys_solution(photosyn_pathway, met, apar, gsdata, sol,   &
     ilimit)

  use c34constants

  implicit none

  

  integer, intent(in) :: photosyn_pathway
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  type(farqdata), intent(in) :: gsdata
  type(solution), intent(inout) :: sol
  integer, intent(out) :: ilimit


  if(photosyn_pathway == 3)then

     call c3solver(met,apar,gsdata,sol,ilimit)

  else

     call c4solver(met,apar,gsdata,sol,ilimit)

  endif

  return

end subroutine exact_lphys_solution

!========================================================================

subroutine store_exact_lphys_solution(old_st_data, met, prss,   &
     leaf_aging_factor, green_leaf_factor, sol, ilimit, gsdata, apar,  &
     photosyn_pathway, Vm0, Vm_low_temp)

  use c34constants
  use therm_lib, only: rslif
  use consts_coms, only: t00,mmdov
  implicit none

  

  integer, intent(in) :: photosyn_pathway
  type(stoma_data), intent(inout) :: old_st_data
  type(metdat), intent(inout) :: met
  real, intent(in) :: prss
  real, intent(in) :: leaf_aging_factor
  real, intent(in) :: green_leaf_factor
  type(solution), intent(in) :: sol
  integer, intent(in) :: ilimit
  type(farqdata), intent(in) :: gsdata
  type(glim), intent(inout) :: apar
  real, intent(in) :: Vm0
  real, intent(in) :: Vm_low_temp


  real, external :: co2cp
  real, external :: arrhenius
  real :: dprss
  real, external :: residual_c3
  real, external :: residual_c4

  ! Save old meteorological information
  old_st_data%T_L = met%tl
  old_st_data%e_a = met%ea
  old_st_data%par = met%par
  old_st_data%rb_factor = met%gbc
  old_st_data%prss = prss
  old_st_data%phenology_factor = leaf_aging_factor / green_leaf_factor
  old_st_data%gsw_open = sol%gsw(2,1)
  old_st_data%ilimit = ilimit
  
  if(ilimit == -1)then

     ! In this case, no open-stomata solution was found to exist.

     old_st_data%gsw_residual = 0.0
     old_st_data%t_l_residual = 0.0
     old_st_data%e_a_residual = 0.0
     old_st_data%par_residual = 0.0
     old_st_data%rb_residual = 0.0
     old_st_data%prss_residual = 0.0
     old_st_data%leaf_residual = 0.0

  else

     if(photosyn_pathway == 3)then

        ! Set parameters
        call setapar_c3(gsdata,met,apar,ilimit)

        ! stomatal conductance derivative
        old_st_data%gsw_residual = residual_c3(gsdata,met,apar,  &
             sol%gsw(2,1)*1.01) / (0.01*sol%gsw(2,1))
        
        ! Temperature derivative
        met%tl = met%tl + 0.1
        met%el = mmdov * rslif(prss,met%tl + t00)
        met%compp = co2cp(met%tl)
        apar%vm = Vm0 * arrhenius(met%tl,1.0,3000.0)  &
             /(1.0+exp(0.4*(Vm_low_temp-met%tl)))  &
             /(1.0+exp(0.4*(met%tl-45.0))) 
        if(leaf_aging_factor > 0.01)then
           apar%vm = apar%vm * leaf_aging_factor   &
                / green_leaf_factor
        endif
        apar%nu = -apar%vm * gsdata%gamma
        apar%k1 = arrhenius(met%tl,1.5e-4,6000.0)
        apar%k2 = arrhenius(met%tl,0.836,-1400.0)
        call setapar_c3(gsdata,met,apar,ilimit)
        old_st_data%T_L_residual = residual_c3(gsdata,met,apar,  &
             sol%gsw(2,1)) / (0.1 * old_st_data%gsw_residual)

        ! Reset parameters
        met%tl = met%tl - 0.1
        met%el = mmdov*rslif(prss,met%tl+t00)
        met%compp = co2cp(met%tl)
        apar%vm = Vm0 * arrhenius(met%tl,1.0,3000.0)  &
             /(1.0+exp(0.4*(Vm_low_temp-met%tl)))  &
             /(1.0+exp(0.4*(met%tl-45.0))) 
        if(leaf_aging_factor > 0.01)then
           apar%vm = apar%vm * leaf_aging_factor   &
                / green_leaf_factor
        endif
        apar%nu = -apar%vm * gsdata%gamma
        apar%k1 = arrhenius(met%tl,1.5e-4,6000.0)
        apar%k2 = arrhenius(met%tl,0.836,-1400.0)
        
        ! humidity derivative
        met%ea = met%ea * 0.99
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0
        call setapar_c3(gsdata,met,apar,ilimit)
        old_st_data%e_a_residual = residual_c3(gsdata,met,apar,sol%gsw(2,1)) &
             / (met%ea*(1.0-1.0/0.99)*old_st_data%gsw_residual)
        met%ea = met%ea / 0.99
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0

        ! PAR derivative
        met%par = met%par * 1.01
        call setapar_c3(gsdata,met,apar,ilimit)
        old_st_data%par_residual = residual_c3(gsdata,met,apar,sol%gsw(2,1)) &
             / (met%par*(1.0-1.0/1.01)*old_st_data%gsw_residual)
        met%par = met%par / 1.01

        ! aerodynamics resistance derivative
        met%gbc = met%gbc * 1.01
        call setapar_c3(gsdata,met,apar,ilimit)
        old_st_data%rb_residual = residual_c3(gsdata,met,apar,sol%gsw(2,1)) &
             / (met%gbc*(1.0-1.0/1.01)*old_st_data%gsw_residual)
        met%gbc = met%gbc / 1.01

        ! pressure derivative
        dprss = prss * 1.005
        met%el = mmdov*rslif(dprss,met%tl+t00)
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0
        call setapar_c3(gsdata,met,apar,ilimit)
        old_st_data%prss_residual = residual_c3(gsdata,met,apar,sol%gsw(2,1)) &
             / (dprss*(1.0-1.0/1.005)*old_st_data%gsw_residual)
        met%el = mmdov*rslif(prss,met%tl+t00)
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0

     else

        ! Now do the same for C4 plants
        
        call setapar_c4(gsdata,met,apar,ilimit)
        old_st_data%gsw_residual = residual_c4(gsdata,met,apar,  &
             sol%gsw(2,1)*1.01) / (0.01*sol%gsw(2,1))
        
        met%tl = met%tl + 0.1
        met%el = mmdov*rslif(prss,met%tl+t00)
        met%compp = co2cp(met%tl)
        apar%vm = Vm0 * arrhenius(met%tl,1.0,3000.0)  &
             /(1.0+exp(0.4*(5.0-met%tl)))/(1.0+exp(0.4*(met%tl-100.0))) 
        call setapar_c4(gsdata,met,apar,ilimit)
        old_st_data%T_L_residual =   &
             residual_c4(gsdata,met,apar,sol%gsw(2,1)) &
             / (0.1 * old_st_data%gsw_residual)
        met%tl = met%tl - 0.1
        met%el = mmdov*rslif(prss,met%tl+t00)
        met%compp = co2cp(met%tl)
        apar%vm = Vm0 * arrhenius(met%tl,1.0,3000.0)  &
             /(1.0+exp(0.4*(5.0-met%tl)))/(1.0+exp(0.4*(met%tl-100.0))) 
        
        met%ea = met%ea * 0.99
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0
        call setapar_c4(gsdata,met,apar,ilimit)
        old_st_data%e_a_residual =   &
             residual_c4(gsdata,met,apar,sol%gsw(2,1)) &
             / (met%ea*(1.0-1.0/0.99)*old_st_data%gsw_residual)
        met%ea = met%ea / 0.99
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0
        
        met%par = met%par * 1.01
        call setapar_c4(gsdata,met,apar,ilimit)
        old_st_data%par_residual =   &
             residual_c4(gsdata,met,apar,sol%gsw(2,1)) &
             / (met%par*(1.0-1.0/1.01)*old_st_data%gsw_residual)
        met%par = met%par / 1.01
        
        met%gbc = met%gbc * 1.01
        call setapar_c4(gsdata,met,apar,ilimit)
        old_st_data%rb_residual =   &
             residual_c4(gsdata,met,apar,sol%gsw(2,1)) &
             / (met%gbc*(1.0-1.0/1.01)*old_st_data%gsw_residual)
        met%gbc = met%gbc / 1.01
        
        dprss = prss * 1.005
        met%el = mmdov*rslif(dprss,met%tl+t00)
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0
        call setapar_c4(gsdata,met,apar,ilimit)
        old_st_data%prss_residual = residual_c4(gsdata,met,apar,sol%gsw(2,1)) &
             / (dprss*(1.0-1.0/1.005)*old_st_data%gsw_residual)
        met%el = mmdov*rslif(prss,met%tl+t00)
        met%eta = 1.0 + (met%el-met%ea)/gsdata%d0
        
     endif
  endif

  return
end subroutine store_exact_lphys_solution

!==================================================================

subroutine fill_lphys_sol_exact(A_open, rsw_open, A_cl, rsw_cl, sol, adens)

  use consts_coms, only : mmdry1000
  use c34constants

  implicit none

  

  real, intent(out) :: A_open
  real, intent(out) :: A_cl
  real, intent(out) :: rsw_open
  real, intent(out) :: rsw_cl
  type(solution), intent(in) :: sol
  real, intent(in) :: adens

  A_open = sol%a(2,1)
  rsw_open = 1.0e9 * adens / (mmdry1000 * sol%gsw(2,1))
!  rsw_open = (2.9e-8 * adens) / sol%gsw(2,1)
  A_cl = sol%a(1,1)
  rsw_cl = 1.0e9 * adens / (mmdry1000 * sol%gsw(1,1))
!  rsw_cl = (2.9e-8 * adens) / sol%gsw(1,1)

  return
end subroutine fill_lphys_sol_exact

!=======================================================================

subroutine fill_lphys_sol_approx(gsdata, met, apar, old_st_data, sol,   &
     A_cl, rsw_cl, adens, rsw_open, A_open, photosyn_pathway, prss)

  use c34constants
  use consts_coms, only : mmdry1000

  implicit none

  

  type(farqdata), intent(in) :: gsdata
  type(metdat), intent(in) :: met
  type(glim), intent(inout) :: apar
  type(stoma_data), intent(in) :: old_st_data
  type(solution), intent(inout) :: sol
  real, intent(out) :: A_cl
  real, intent(out) :: A_open
  real, intent(out) :: rsw_cl 
  real, intent(out) :: rsw_open
  real, intent(in) :: adens
  integer, intent(in) :: photosyn_pathway
  real, intent(in) :: prss

  integer :: success
  real :: gsw_update
  real :: ci_approx
  real, external :: aflux_c4
  real, external :: aflux_c3
  real, external :: quad4ci

  if(photosyn_pathway == 3)then
     call setapar_c3(gsdata,met,apar,max(1,old_st_data%ilimit))
     call solve_closed_case_c3(gsdata,met,apar,sol,1)
  else
     call setapar_c4(gsdata,met,apar,max(1,old_st_data%ilimit))
     call solve_closed_case_c4(gsdata,met,apar,sol,1)
  endif

  A_cl = sol%a(1,1)
  rsw_cl = 1.0e9 * adens / (mmdry1000 * sol%gsw(1,1))
!  rsw_cl = (2.9e-8 * adens) / sol%gsw(1,1)
  if(old_st_data%ilimit /= -1)then
     gsw_update = old_st_data%gsw_open - &
          old_st_data%t_l_residual * (met%tl - old_st_data%t_l) - &
          old_st_data%e_a_residual * (met%ea - old_st_data%e_a) - &
          old_st_data%par_residual * (met%par - old_st_data%par) - &
          old_st_data%rb_residual * (met%gbc - old_st_data%rb_factor) - &
          old_st_data%prss_residual * (prss - old_st_data%prss)
     
     if(photosyn_pathway == 3)then
        ci_approx = quad4ci(gsdata,met,apar,gsw_update,success)
        if(success == 1)then
           A_open = aflux_c3(apar,gsw_update,ci_approx)
           rsw_open = 1.0e9 * adens / (mmdry1000 * gsw_update)
!           rsw_open = (2.9e-8 * adens) / gsw_update
        else
           A_open = A_cl
           rsw_open = rsw_cl
        endif
     else
        A_open = aflux_c4(apar,met,gsw_update)
        rsw_open = 1.0e9 * adens / (mmdry1000 * gsw_update)
!        rsw_open = (2.9e-8 * adens) / gsw_update
     endif
  else
     A_open = A_cl
     rsw_open = rsw_cl
  endif

  return
end subroutine fill_lphys_sol_approx

