!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine fillscr(m1,m2,m3,n1,n2,n3,k1,k2,scr,var)

  ! fillscr:
  !  copy selected verticals of the "smaller" array var
  !  into "larger" array scr
  !
  !  copies var(k1:k2,1:n2,1:n3) into same portion of scr
  !  it assumes that scr is large enough to contain var - see declarations

  implicit none
  integer, intent(in   ) :: m1
  integer, intent(in   ) :: m2
  integer, intent(in   ) :: m3
  integer, intent(in   ) :: n1
  integer, intent(in   ) :: n2
  integer, intent(in   ) :: n3
  integer, intent(in   ) :: k1
  integer, intent(in   ) :: k2
  real,    intent(inout) :: scr(m1,m2,m3)   ! changed only at (k1:k2,1:n2,1:n3)
  real,    intent(in   ) :: var(n1,n2,n3)

  integer :: i,j,k

  do j = 1,n3
     do i = 1,n2
        do k = k1,k2
           scr(k,i,j) = var(k,i,j)
        end do
     end do
  end do
end subroutine fillscr






subroutine fillvar(m1,m2,m3,n1,n2,n3,k1,k2,scr,var)

  ! fillvar:
  !  copy selected verticals of the "larger" array scr
  !  into "smaller" array var
  !
  !  copies scr(k1:k2,1:n2,1:n3) into same portion of var
  !  it assumes that scr is at least as large as var
  !  at any dimension
  
  implicit none
  integer, intent(in   ) :: m1
  integer, intent(in   ) :: m2
  integer, intent(in   ) :: m3
  integer, intent(in   ) :: n1
  integer, intent(in   ) :: n2
  integer, intent(in   ) :: n3
  integer, intent(in   ) :: k1
  integer, intent(in   ) :: k2
  real,    intent(in   ) :: scr(m1,m2,m3)
  real,    intent(inout) :: var(n1,n2,n3)   ! changed only at (k1:k2,1:n2,1:n3)

  integer :: i,j,k

  do j = 1,n3
     do i = 1,n2
        do k = k1,k2
           var(k,i,j) = scr(k,i,j)
        end do
     end do
  end do
end subroutine fillvar






subroutine dnswt2(nzmx,nxmx,nymx,n1,n2,n3,var1,dn0x,vnam,idir)

  ! dnswt2:
  !  weights (or "unweights") selected portions of array var1 by density;
  !  weights are applied accordingly to field (indicated by vnam):
  !    vertical interpolation on "w" up to n1-1
  !    no interpolation otherwise, up to n1
  !  
  !  idir controls the direction: 
  !    idir=1 means weight (multiply)
  !    idir=2 means unweight (divide)
  !
  !  it assumes that var1, dimensioned (nzmx,nxmx,nymx) is "larger"
  !  than density, dimensioned (n1,n2,n3)
  !
  !  var1 is changed only at:
  !    (1:n1,1:n2,1:n3)   if vnam /= "w"
  !    (1:n1-1,1:n2,1:n3) if vnam == "w"

  implicit none
  integer,          intent(in   ) :: nzmx
  integer,          intent(in   ) :: nxmx
  integer,          intent(in   ) :: nymx
  integer,          intent(in   ) :: n1
  integer,          intent(in   ) :: n2
  integer,          intent(in   ) :: n3
  integer,          intent(in   ) :: idir
  real,             intent(in   ) :: dn0x(n1,n2,n3)
  character(len=*), intent(in   ) :: vnam
  real,             intent(inout) :: var1(nzmx,nxmx,nymx) ! portions changed (see above)

  integer :: i,j,k

  if (idir == 1) then

     if (vnam == 'w') then
        do j = 1,n3
           do i = 1,n2
              do k = 1,n1-1
                 var1(k,i,j) = var1(k,i,j)  &
                      * .5 * (dn0x(k,i,j) + dn0x(k+1,i,j))
              end do
           end do
        end do
     else
        do j = 1,n3
           do i = 1,n2
              do k = 1,n1
                 var1(k,i,j) = var1(k,i,j) * dn0x(k,i,j)
              end do
           end do
        end do
     end if

  elseif (idir == 2) then

     if (vnam == 'w') then
        do j = 1,n3
           do i = 1,n2
              do k = 1,n1-1
                 var1(k,i,j) = var1(k,i,j)  &
                      / (.5 * (dn0x(k,i,j) + dn0x(k+1,i,j)))
              end do
           end do
        end do
     else
        do j = 1,n3
           do i = 1,n2
              do k = 1,n1
                 var1(k,i,j) = var1(k,i,j) / dn0x(k,i,j)
              end do
           end do
        end do
     end if
  end if
end subroutine dnswt2
