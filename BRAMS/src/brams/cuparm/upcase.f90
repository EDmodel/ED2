subroutine upcase(s1,n)
  
  implicit none

  integer           :: n
  character (LEN=n) :: s1

  integer           :: i, c

  if (n < 1) then
     print *, "Conversion UpCase Error!"
     print *, "String lenght: ", n
     return
  endif
  do i=1,n
     c = iachar(s1(i:i))
     if (c >= 97) then
        c = c - 32
        s1(i:i) = achar(c)
     endif
  enddo

end subroutine upcase

subroutine VER_GRELL_CLOSURE(CLOSURE_TYPE)

  implicit none

  character (LEN=2) :: CLOSURE_TYPE

  if ((CLOSURE_TYPE /= 'EN').and.(CLOSURE_TYPE /= 'GR').and.   &
       (CLOSURE_TYPE /= 'LO').and.(CLOSURE_TYPE /= 'MC').and.  &
       (CLOSURE_TYPE /= 'SC').and.(CLOSURE_TYPE /= 'AS')) then
     print *, "****Grell Closure type ERROR!!"
     print *, "    CLOSURE_TYPE: ",CLOSURE_TYPE
     print *, "    Program will stop!!"
     stop
  endif
  
end subroutine VER_GRELL_CLOSURE
