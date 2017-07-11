subroutine alt_dia(idate1, imonth1, iyear1, INT_DIF_TIME,   &
     idate2, imonth2, iyear2)
  implicit none
  integer :: idate1, imonth1, iyear1, INT_DIF_TIME, idate2, imonth2, iyear2
  integer :: i, increm, DMES(12)

  data DMES/31,28,31,30,31,30,31,31,30,31,30,31/

  iyear2  = iyear1
  imonth2 = imonth1
  idate2  = idate1

  increm  = 1

  if (INT_DIF_TIME < 1) increm = increm*(-1)

  do i = 1, abs(INT_DIF_TIME)
     idate2 = idate2 + increm
     if (idate2 < 1) then
        imonth2 = imonth2 + increm
        if (imonth2 < 1) then
           imonth2 = 12
           iyear2 = iyear2 - 1
        endif
        idate2 = DMES(imonth2)
     elseif (idate2 > DMES(imonth2)) then
        imonth2 = imonth2 + increm
        if (imonth2 > 12) then
           imonth2 = 1
           iyear2 = iyear2 + 1
        endif
        idate2 = 1
     endif

  enddo

end subroutine alt_dia
