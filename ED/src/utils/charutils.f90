!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine deblank(str1,str2,nch)
implicit none
character(len=*) :: str1,str2
integer :: n,ln,nch

! strips blanks from a string and returns number of chars

str2=' '
ln=len(str1)
nch=0
do n=1,ln
   if(str1(n:n).ne.' ') then
      nch=nch+1
      str2(nch:nch)=str1(n:n)
   endif
enddo

return
end
!***************************************************************************

subroutine parse(str,tokens,ntok)
implicit none
integer :: ntok
character(len=*) :: str,tokens(*)
character(len=1) :: sep
integer, parameter :: ntokmax=100

integer :: n,nc,npt,nch,ntbeg,ntend

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

sep=' '
ntok=0
npt=1
nch=len_trim(str)
nc=1
do ntok=1,ntokmax
   do n=nc,nch
      if(str(n:n).ne.sep) then
         ntbeg=n
         goto 21
      endif
   enddo
   21 continue
   do n=ntbeg,nch
      if(str(n:n).eq.sep) then
         ntend=n-1
         goto 22
      endif
      if(n.eq.nch) then
         ntend=n
         goto 22
      endif
   enddo
   22 continue
   tokens(ntok)=str(ntbeg:ntend)
   nc=ntend+1
   if(nc.ge.nch) goto 25
enddo

25 continue


return
end

!***************************************************************************

subroutine tokenize1(str1,tokens,ntok,toksep)
implicit none
integer :: ntok
character(len=*) :: str1,tokens(*)
character(len=1), intent(in) :: toksep

character(len=256) :: str
integer :: nch,ist,npt,nc

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

call deblank(str1,str,nch)

ist=1
if(str(1:1).eq.toksep) ist=2
npt=ist
ntok=0
do nc=ist,nch
   if(str(nc:nc).eq.toksep.or.nc.eq.nch) then
      if(nc-npt.ge.1) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc-1)
         if(nc.eq.nch.and.str(nc:nc).ne.toksep) then
            tokens(ntok)=str(npt:nc)
            goto 10
         endif
         npt=nc+1
      endif
   endif
enddo
10 continue

return
end
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine tolower(word,dimword)
!------------------------------------------------------------------------------------------!
! Subroutine tolower                                                                       !
!                                                                                          !
!    This subroutine converts all common upper-case characters into lowercase.             !
!------------------------------------------------------------------------------------------!
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)                                 :: dimword
   character(len=*), dimension(dimword), intent(inout) :: word
   !----- Internal variables --------------------------------------------------------------!
   integer                                             :: wmax
   integer                                             :: w
   integer                                             :: d
   integer                                             :: inow
   !----- Local constant. -----------------------------------------------------------------!
   integer, parameter                                  :: iau   = iachar('A')
   integer, parameter                                  :: izu   = iachar('Z')
   integer, parameter                                  :: delta = iachar('e') - iachar('E')
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Outer loop, the words.                                                           !
   !---------------------------------------------------------------------------------------!
   wordloop: do d=1,dimword
      wmax=len_trim(word(d))

      !------------------------------------------------------------------------------------!
      !     Inner loop, word characters.                                                   !
      !------------------------------------------------------------------------------------!
      charloop: do w=1,wmax


         !----- Only letters between A-Z (no accented letters) will be converted. ---------!

         inow = iachar(word(d)(w:w))
         if (inow >= iau .and. inow <= izu) then
            word(d)(w:w) = achar(iachar(word(d)(w:w)) + delta)
         end if
      end do charloop
      !------------------------------------------------------------------------------------!
   end do wordloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine tolower
