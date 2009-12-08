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
character*(*) str1,tokens(*),str*256
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
character*(*) str1,tokens(*),str*256
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
dimension nums(*)
character cstr(*)*(*),cscr*200

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
dimension xnums(*)
character cstr(*)*(*),cscr*200

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
