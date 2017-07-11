!############################# Change Log ##################################
! 1.0.0.1
!
! 001002 MJB char_strip_var ##
!            Replaced index calls with f90 intrinsic len_trim. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

integer function lastchar(str)
character*(*) str

! returns last non-blank character position from a string

ln=len(str)
do n=ln,1,-1
   if(str(n:n).ne.' ') then
      lastchar=n
      return
   endif
enddo
lastchar=0

return
end

!***************************************************************************

integer function ifirstchar(str)
character*(*) str

! returns last non-blank character position from a string

ln=len(str)
do n=1,ln
   if(str(n:n).ne.' ') then
      ifirstchar=n
      return
   endif
enddo
ifirstchar=1

return
end

!***************************************************************************

subroutine deblank(str1,str2,nch)
character*(*) str1,str2

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

subroutine detab(str1,str2,nch)
character*(*) str1,str2
character*1  tab

tab=achar( 9)

! strips tabs from a string and returns number of chars

str2=' '
ln=lastchar(str1)
nch=0
do n=1,ln
   if(str1(n:n).ne.tab) then
      !print*,'no tab:',str1(n:n)
      nch=nch+1
      str2(nch:nch)=str1(n:n)
   else
      print*,'found one:',str1
      str2(nch+1:nch+6)='      '
      nch=nch+6
   endif
enddo

return
end

!***************************************************************************

integer function lastslash(str)
character*(*) str

! returns last slash character position from a string

ln=len(str)
do n=ln,1,-1
   if(str(n:n).eq.'/') then
      lastslash=n
      return
   endif
enddo
lastslash=0

return
end

!***************************************************************************

subroutine char_strip_var(line,var,line2)
character*(*) line,var,line2

! removes instances of a substring from a string

ncl=len(line)
do nn=1,ncl
   if(line(nn:nn).ne.' ') then
      nb=index(line(nn:),' ')
      var=line(nn:nn+nb-1)
      goto 25
   endif
enddo
25 continue
line2=line(nn+nb-1:)

return
end

!***************************************************************************

subroutine findln(text,ltext,order)
character*(*) text
integer ltext,order

! find first non-blank character if order=0, last non-blank if order=1

if(order.eq.1) then
   do i=len(text),1,-1
      if(text(i:i).ne.' ') then
         ltext=i
         goto 10
      endif
   enddo
   10 continue
else
   do i=1,len(text)
      if(text(i:i).ne.' ') then
         ltext=i
         goto 20
      endif
   enddo
   20 continue
endif

return
end

!***************************************************************************

subroutine parse(str,tokens,ntok)
character*(*) str,tokens(*)
character sep*1
data ntokmax/100/

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

sep=' '
ntok=0
npt=1
nch=lastchar(str)
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

!do nc=1,nch
!   if(str(nc:nc).eq.sep.or.nc.eq.nch)then
!      if(nc-npt.ge.1)then
!         ntok=ntok+1
!         tokens(ntok)=str(npt:nc-1)
!         if(nc.eq.nch.and.str(nc:nc).ne.sep)then
!            tokens(ntok)=str(npt:nc)
!            go to 10
!         endif
!      endif
!      ntok=ntok+1
!      tokens(ntok)=str(nc:nc)
!      npt=nc+1
!      go to 10
!   endif
!   10 continue
!enddo

return
end

!***************************************************************************

subroutine tokenize(str1,tokens,ntok,toksep,nsep)
use rpost_dims, only : str_len
character*(*) str1,tokens(*)
character(len=str_len) :: str
character*1 toksep(nsep)

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

ntok=0
npt=1
call deblank(str1,str,nch)
do nc=1,nch
   do ns=1,nsep
      if(str(nc:nc).eq.toksep(ns).or.nc.eq.nch) then
         if(nc-npt.ge.1)then
            ntok=ntok+1
            tokens(ntok)=str(npt:nc-1)
            if(nc.eq.nch.and.str(nc:nc).ne.toksep(ns)) then
               tokens(ntok)=str(npt:nc)
               goto 10
            endif
         endif
         ntok=ntok+1
         tokens(ntok)=str(nc:nc)
         npt=nc+1
         goto 10
      endif
   enddo
10      continue
enddo
return
end

!***************************************************************************

subroutine tokenize1(str1,tokens,ntok,toksep)
use rpost_dims, only : str_len
character*(*) str1,tokens(*)
character(len=str_len) :: str
character*1 toksep

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

!***************************************************************************

subroutine tokfind(toks,ntok,str,iff)
character*(*) toks(*),str

! looks for a number of tokens (substrings) within a string

do n=1,ntok
   !print*,'tokfind-',n,toks(n)(1:lastchar(toks(n)))  &
   !      ,'=====',str(1:lastchar(str))
   if(str(1:lastchar(str)).eq.toks(n)(1:lastchar(toks(n)))) then
      iff=1
      return
   endif
enddo
iff=0

return
end

!***************************************************************************

subroutine rams_intsort(ni,nums,cstr)
use rpost_dims, only : str_len
dimension nums(*)
character cstr(*)*(*)
character(len=str_len) :: cscr

! sort an array of character strings by an associated integer field

do n=1,ni
   mini=1000000
   do nm=n,ni
      if(nums(nm).lt.mini) then
         nmm=nm
         mini=nums(nm)
      endif
   enddo
   nscr=nums(n)
   nums(n)=nums(nmm)
   nums(nmm)=nscr
   cscr=cstr(n)
   cstr(n)=cstr(nmm)
   cstr(nmm)=cscr
enddo

return
end

!***************************************************************************

subroutine rams_fltsort(ni,xnums,cstr)
use rpost_dims, only : str_len
dimension xnums(*)
character cstr(*)*(*)
character(len=str_len) :: cscr

! sort an array of character strings by an associated float field

do n=1,ni
   xmini=1.e30
   do nm=n,ni
      if(xnums(nm).lt.xmini) then
         nmm=nm
         xmini=xnums(nm)
      endif
   enddo
   xnscr=xnums(n)
   xnums(n)=xnums(nmm)
   xnums(nmm)=xnscr
   cscr=cstr(n)
   cstr(n)=cstr(nmm)
   cstr(nmm)=cscr
enddo

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
!----- Arguments --------------------------------------------------------------------------!
  integer, intent(in)                                 :: dimword
  character(len=*), dimension(dimword), intent(inout) :: word
!----- Internal variables -----------------------------------------------------------------!
  integer                                       :: wmax,w,d
!------------------------------------------------------------------------------------------!
  do d=1,dimword
    wmax=len_trim(word(d))
    do w=1,wmax
      select case(word(d)(w:w))
      case('A') 
        word(d)(w:w)='a'
      case('B')
        word(d)(w:w)='b'
      case('C')
        word(d)(w:w)='c'
      case('D')
        word(d)(w:w)='d'
      case('E')
        word(d)(w:w)='e'
      case('F')
        word(d)(w:w)='f'
      case('G')
        word(d)(w:w)='g'
      case('H')
        word(d)(w:w)='h'
      case('I')
        word(d)(w:w)='i'
      case('J')
        word(d)(w:w)='j'
      case('K')
        word(d)(w:w)='k'
      case('L')
        word(d)(w:w)='l'
      case('M')
        word(d)(w:w)='m'
      case('N')
        word(d)(w:w)='n'
      case('O')
        word(d)(w:w)='o'
      case('P')
        word(d)(w:w)='p'
      case('Q')
        word(d)(w:w)='q'
      case('R')
        word(d)(w:w)='r'
      case('S')
        word(d)(w:w)='s'
      case('T')
        word(d)(w:w)='t'
      case('U')
        word(d)(w:w)='u'
      case('V')
        word(d)(w:w)='v'
      case('W')
        word(d)(w:w)='w'
      case('X')
        word(d)(w:w)='x'
      case('Y')
        word(d)(w:w)='y'
      case('Z')
        word(d)(w:w)='z'
      case('Á')
        word(d)(w:w)='á'
      case('É')
        word(d)(w:w)='é'
      case('Í')
        word(d)(w:w)='í'
      case('Ó')
        word(d)(w:w)='ó'
      case('Ú')
        word(d)(w:w)='ú'
      case('Ý')
        word(d)(w:w)='ý'
      case('À')
        word(d)(w:w)='à'
      case('È')
        word(d)(w:w)='è'
      case('Ì')
        word(d)(w:w)='ì'
      case('Ò')
        word(d)(w:w)='ò'
      case('Ù')
        word(d)(w:w)='ù'
      case('Â')
        word(d)(w:w)='â'
      case('Ê')
        word(d)(w:w)='ê'
      case('Î')
        word(d)(w:w)='î'
      case('Ô')
        word(d)(w:w)='ô'
      case('Û')
        word(d)(w:w)='û'
      case('Ä')
        word(d)(w:w)='ä'
      case('Ë')
        word(d)(w:w)='ë'
      case('Ï')
        word(d)(w:w)='ï'
      case('Ö')
        word(d)(w:w)='ö'
      case('Ü')
        word(d)(w:w)='ü'
      case('Ã')
        word(d)(w:w)='ã'
      case('Õ')
        word(d)(w:w)='õ'
      case('Ñ')
        word(d)(w:w)='ñ'
      case('Å')
        word(d)(w:w)='å'
      case('Ç')
        word(d)(w:w)='ç'
      end select
    end do
  end do
  return
end subroutine tolower
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine tolower                                                                       !
!                                                                                          !
!    This subroutine converts all common upper-case characters into lowercase.             !
!------------------------------------------------------------------------------------------!
subroutine tolower_sca(word)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*), intent(inout) :: word
   !----- Internal variables --------------------------------------------------------------!
   integer                         :: wmax
   integer                         :: w
   !---------------------------------------------------------------------------------------!

   wmax=len_trim(word)
   do w=1,wmax
     select case(word(w:w))
     case('A') 
       word(w:w)='a'
     case('B')
       word(w:w)='b'
     case('C')
       word(w:w)='c'
     case('D')
       word(w:w)='d'
     case('E')
       word(w:w)='e'
     case('F')
       word(w:w)='f'
     case('G')
       word(w:w)='g'
     case('H')
       word(w:w)='h'
     case('I')
       word(w:w)='i'
     case('J')
       word(w:w)='j'
     case('K')
       word(w:w)='k'
     case('L')
       word(w:w)='l'
     case('M')
       word(w:w)='m'
     case('N')
       word(w:w)='n'
     case('O')
       word(w:w)='o'
     case('P')
       word(w:w)='p'
     case('Q')
       word(w:w)='q'
     case('R')
       word(w:w)='r'
     case('S')
       word(w:w)='s'
     case('T')
       word(w:w)='t'
     case('U')
       word(w:w)='u'
     case('V')
       word(w:w)='v'
     case('W')
       word(w:w)='w'
     case('X')
       word(w:w)='x'
     case('Y')
       word(w:w)='y'
     case('Z')
       word(w:w)='z'
     case('Á')
       word(w:w)='á'
     case('É')
       word(w:w)='é'
     case('Í')
       word(w:w)='í'
     case('Ó')
       word(w:w)='ó'
     case('Ú')
       word(w:w)='ú'
     case('Ý')
       word(w:w)='ý'
     case('À')
       word(w:w)='à'
     case('È')
       word(w:w)='è'
     case('Ì')
       word(w:w)='ì'
     case('Ò')
       word(w:w)='ò'
     case('Ù')
       word(w:w)='ù'
     case('Â')
       word(w:w)='â'
     case('Ê')
       word(w:w)='ê'
     case('Î')
       word(w:w)='î'
     case('Ô')
       word(w:w)='ô'
     case('Û')
       word(w:w)='û'
     case('Ä')
       word(w:w)='ä'
     case('Ë')
       word(w:w)='ë'
     case('Ï')
       word(w:w)='ï'
     case('Ö')
       word(w:w)='ö'
     case('Ü')
       word(w:w)='ü'
     case('Ã')
       word(w:w)='ã'
     case('Õ')
       word(w:w)='õ'
     case('Ñ')
       word(w:w)='ñ'
     case('Å')
       word(w:w)='å'
     case('Ç')
       word(w:w)='ç'
     end select
   end do
   return
end subroutine tolower_sca
!==========================================================================================!
!==========================================================================================!
