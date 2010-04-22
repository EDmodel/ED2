!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
subroutine patch_array_size(npq,deltax  &
   ,ivegtflg,ivegtfn,isoilflg,isoilfn,ndviflg,ndvifn)
   use rconstants, only : spcon

   implicit none

   integer :: npq,ivegtflg,isoilflg,ndviflg
   character(len=*) :: ivegtfn,isoilfn,ndvifn
   real :: deltax

   integer :: iblksiz,no,isbeg,iwbeg
   real :: offlat,offlon,deltall,deltallo_min

   ! Find highest resolution dataset among all that are to be used

   deltallo_min = 1.e20
   if (ivegtflg == 1) then
      call read_header(ivegtfn,iblksiz,no,isbeg  &
         ,iwbeg,offlat,offlon,deltall,'veg')
      deltallo_min = min(deltallo_min,deltall)
   endif

   if (isoilflg == 1) then
      call read_header(isoilfn,iblksiz,no,isbeg  &
         ,iwbeg,offlat,offlon,deltall,'soil')
      deltallo_min = min(deltallo_min,deltall)
   endif

   if (ndviflg == 1) then
      call read_header(ndvifn,iblksiz,no,isbeg  &
         ,iwbeg,offlat,offlon,deltall,'ndvi')
      deltallo_min = min(deltallo_min,deltall)
   endif

   npq = min(10,max(1,nint(deltax / (deltallo_min * spcon))))

   return
end subroutine patch_array_size
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine patch_latlon(n2,n3,xt,yt,deltax,deltay,platn,plonn)

   use mem_mksfc
   implicit none

   integer :: n2,n3
   real :: xt(n2),yt(n3),platn,plonn

   integer :: jr,jp,ip,ir
   real :: yp,xp,deltax,deltay,deltaxp,deltayp

   ! Fill arrays with offset latitudes and longitudes of all p points

      deltaxp = deltax / float(npq)
      deltayp = deltay / float(npq)
      do jr = 1,n3
         do jp = 1,npq
            yp = yt(jr) + (float(jp) - .5 * float(npq+1)) * deltayp
            do ir = 1,n2
               do ip = 1,npq
                  xp = xt(ir) + (float(ip) - .5 * float(npq+1)) * deltaxp
                  call xy_ll(glatp(ip,jp,ir,jr),glonp(ip,jp,ir,jr)  &
                     ,platn,plonn,xp,yp)
               enddo
            enddo
         enddo
      enddo

   return
end subroutine patch_latlon
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Subroutine landuse_opr reads in one or more landuse data types and defines and fills  !
! subgrid patches from them.  Currently, landuse (vegetation) class, soil textural class,  !
! and ndvi value are implemented.  In this version of landuse_opr, landuse class must be   !
! used if any datasets are used.  Patch areas are determined from landuse class alone. If  !
! soil textural class data is used, the most dominant type occurring in each               !
! landuse-class-defined patch is assigned to that patch.  If ndvi data is used, the        !
! average value for each landuse-class-defined patch is assigned to that patch.            !
!------------------------------------------------------------------------------------------!
subroutine landuse_opqr(n2,n3,mzg,npat,nvegpat  &
   ,ivegtflg,ivegtfn,isoilflg,isoilfn,ndviflg,ndvifn,cndvifil   &
   ,iaction,platn,plonn     &
   ,soil_text,patch_area,leaf_class,veg_ndvif)

   use mem_mksfc
   use rconstants

   implicit none
   integer :: n2,n3,mzg,npat,nvegpat
   character(len=*) :: iaction,ivegtfn,isoilfn,ndvifn,cndvifil
   integer :: ivegtflg,isoilflg,ndviflg,ing_prt,nn,nmiss
   real :: platn,plonn
   real, dimension(mzg,n2,n3,npat) :: soil_text
   real, dimension(n2,n3,npat) :: patch_area,leaf_class,veg_ndvif
   integer, parameter :: maxmiss=1000
   character(len=80) :: fnmiss(maxmiss)
    
   integer :: i,j
   real :: checksum

   integer, parameter :: maxdatq=32,nsoil=12

   integer :: datq,lsp,ngrdpix(0:maxdatq,2),datq_pat  &
             ,datsoil,soil_count,nlev_soil      &
             ,jr,jp,ir,ip,iqv,ing,maxdq,jng,jng1,ngstor1    &
             ,ngstor2,npatpixs,nwat,ipat,idatq,isoil,jsoil,k      &
             ,no_veg ,iblksizo_veg ,isbego_veg ,iwbego_veg   &
             ,no_soil,iblksizo_soil,isbego_soil,iwbego_soil  &
             ,no_ndvi,iblksizo_ndvi,isbego_ndvi,iwbego_ndvi

   integer, dimension(0:maxdatq) :: sumpix
   integer, dimension(0:maxdatq,nsoil) :: datq_soil,soil_tab

   real :: deltallo_veg ,deltallo_soil ,deltallo_ndvi  &
          ,offlat_veg   ,offlat_soil   ,offlat_ndvi    &
          ,offlon_veg   ,offlon_soil   ,offlon_ndvi    &
          ,xp,yp,fracwat,plpp

   real, dimension(0:maxdatq) :: datq_ndvi,sumndvi


   nmiss=0

   !print*,n2,n3,mzg,npat,nvegpat,npq
   !print*,iaction,platn,plonn
   checksum=0

   if (iaction .eq. 'veg') then

      

      call read_header(ivegtfn,iblksizo_veg,no_veg,isbego_veg  &
         ,iwbego_veg,offlat_veg,offlon_veg,deltallo_veg,'veg')
         
      call fill_datp(n2,n3,no_veg,iblksizo_veg,isbego_veg,iwbego_veg  &
         ,platn,plonn,offlat_veg,offlon_veg,deltallo_veg,ivegtfn,iaction  &
         ,nmiss,fnmiss)

   ! 9/30/97:  Carry out the first translation of the input DATP values into a
   ! condensed set called DATQ_patch.  The range of DATQ_patch values represents
   ! the total variety of land surface conditions (patches) to be allowed for 
   ! the present simulation, and may be a broader class than the LEAF-2 
   ! vegetation classes for which all the vegetation physical parameters are 
   ! defined.  For example, two different DATQ_patch classes may be mapped to 
   ! the same LEAF-2 vegetation class, but be initialized with different soil 
   ! moistures or soil types, and therefore require different patches.

   ! Fill datq_patch (patch class) values from input landuse dataset.
   ! Currently, this data serves as the primary criterion for defining patches.

   !   print*,(datp(i,i,10,10),i=1,npq)


      do jr = 1,n3
         do ir = 1,n2

            do iqv = 0,maxdatq
               ngrdpix(iqv,1) = 0     ! Initialize counter for datq pixels
               ngrdpix(iqv,2) = iqv   ! Initialize array of consecutive datq values 
            enddo

            do jp = 1,npq
               do ip = 1,npq

   !write(6,202) ip,jp,ir,jr,datp(ip,jp,ir,jr)
   !202 format('pp2',4i5,f5.1)

                  call datp_datq(datp(ip,jp,ir,jr),datq_patch(ip,jp,ir,jr))
                  datq = datq_patch(ip,jp,ir,jr)
                  ngrdpix(datq,1) = ngrdpix(datq,1) + 1

   !write(6,203) ip,jp,ir,jr,datp(ip,jp,ir,jr),datq_patch(ip,jp,ir,jr)  &
   !    ,ngrdpix(datq,1)
   !203 format('pp3',4i5,f5.1,2i5)

               enddo
            enddo

   ! Sort values of ngrdpix by prevalence for non-water patches (datq .ge. 2)

            do ing = 2,maxdatq
               maxdq = -1
               do jng = ing,maxdatq
                  if (ngrdpix(jng,1) .gt. maxdq) then
                     jng1 = jng
                     maxdq = ngrdpix(jng,1)
                  endif
               enddo
               ngstor1 = ngrdpix(ing,1)
               ngstor2 = ngrdpix(ing,2)
               ngrdpix(ing,1) = ngrdpix(jng1,1)
               ngrdpix(ing,2) = ngrdpix(jng1,2)
               ngrdpix(jng1,1) = ngstor1
               ngrdpix(jng1,2) = ngstor2
               
            enddo

            !do ing_prt = 1,maxdatq
            !    write(6,204)ir,jr,ing_prt,ngrdpix(ing_prt,1),ngrdpix(ing_prt,2)
            !204 format('ppx',5i5)
            !enddo
            
            ! Fill patches numbered 2 through nvegpat+1 with the nvegpat most prevalent
            !    nonwater landuse types.  Count pixels for these patches for normalization
            !    of total grid cell area

            npatpixs = 1
            nwat = ngrdpix(0,1) + ngrdpix(1,1)

            if (nwat .lt. npq*npq) then

               npatpixs = 0
               do ipat = 2,nvegpat+1

                  datq = ngrdpix(ipat,2)
                  leaf_class(ir,jr,ipat) = float(datq)
                  npatpixs = npatpixs + ngrdpix(ipat,1)


                  !write(6,205) ir,jr,ipat,datq,leaf_class(ir,jr,ipat),npatpixs
                  !205 format('ppy',4i5,f10.2,i5)
               enddo
            endif

            fracwat = float(nwat) / float(npq * npq)
            plpp = (1. - fracwat) / float(npatpixs)
            patch_area(ir,jr,1) = fracwat

            if (ngrdpix(0,1) .ge. ngrdpix(1,1)) then
               leaf_class(ir,jr,1) = 0.
            else
               leaf_class(ir,jr,1) = 1.
            endif

            do ipat = 2,nvegpat+1
               patch_area(ir,jr,ipat) = plpp * float(ngrdpix(ipat,1))
               
               !write(6,207) ir,jr,ipat,patch_area(ir,jr,ipat),plpp  &
               !    ,ngrdpix(ipat,1),npatpixs,fracwat
               !207 format('pp4',3i5,2f7.3,2i5,f7.3)
               

            enddo

         enddo
      enddo

   elseif (iaction .eq. 'soil') then

      call read_header(isoilfn,iblksizo_soil,no_soil,isbego_soil  &
         ,iwbego_soil,offlat_soil,offlon_soil,deltallo_soil,'soil')
         
         
      call fill_datp(n2,n3,no_soil,iblksizo_soil,isbego_soil,iwbego_soil  &
         ,platn,plonn,offlat_soil,offlon_soil,deltallo_soil,isoilfn,iaction  &
         ,nmiss,fnmiss)

      do jr = 1,n3
         do ir = 1,n2
            if (patch_area(ir,jr,1) .le. .9999) then

               do idatq = 0,maxdatq
                  do isoil = 1,nsoil
                     datq_soil(idatq,isoil) = 0     ! Initialize counter for datq pixels
                  enddo
               enddo

               do jp = 1,npq
                  do ip = 1,npq

   ! Fill datq_soil values as secondary criterion.  This is for finding dominant
   ! soil class for each datq class.

                     


                     call datp_datsoil(datp(ip,jp,ir,jr),datsoil)
                     
                     datq_pat = datq_patch(ip,jp,ir,jr)
                     datq_soil(datq_pat,datsoil)  &
                        = datq_soil(datq_pat,datsoil) + 1
    

                     !write(6,440)  ip,jp,ir,jr,datp(ip,jp,ir,jr),datsoil,datq_pat  &
                     !   ,datq_soil(datq_pat,datsoil)
                     !440 format('action veg ',4i4,f6.1,i5,i5,i5)

                  enddo
               enddo

               do ipat = 2,nvegpat+1

                  if (patch_area(ir,jr,ipat) .ge. .0001) then   

                     datq_pat = nint(leaf_class(ir,jr,ipat))

   ! Find isoil value for which soil_tab(datq_pat,isoil) is a maximum

                     
                     soil_count = 0
                     do isoil = 1,nsoil
                        if (datq_soil(datq_pat,isoil) .gt. soil_count) then
                           soil_count = datq_soil(datq_pat,isoil)
                           jsoil = isoil
                        endif
                     enddo

                     checksum=checksum+float(jsoil)
                     
   ! For now, assume single level of input soil data (e.g., FAO) and
   ! fill all soil levels with this value.

                     do k = 1,mzg
                        if (jsoil /= 0) soil_text(k,ir,jr,ipat) = float(jsoil)
                     enddo

                  endif
               enddo

            endif
         enddo
      enddo



   elseif (iaction .eq. 'ndvi') then

      call read_header(ndvifn,iblksizo_ndvi,no_ndvi,isbego_ndvi  &
         ,iwbego_ndvi,offlat_ndvi,offlon_ndvi,deltallo_ndvi,'ndvi')
         
      call fill_datp(n2,n3,no_ndvi,iblksizo_ndvi,isbego_ndvi,iwbego_ndvi  &
         ,platn,plonn,offlat_ndvi,offlon_ndvi,deltallo_ndvi,cndvifil,iaction  &
         ,nmiss,fnmiss)

      do jr = 1,n3
         do ir = 1,n2
            if (patch_area(ir,jr,1) .le. .9999) then

               do idatq = 0,maxdatq
                  sumndvi(idatq) = 0.  ! initialize ndvi sum
                  sumpix(idatq) = 0    ! initialize ndvi pixel count
               enddo

               do jp = 1,npq
                  do ip = 1,npq

   ! Fill datq_ndvi values to compute average value for each datq class.



                     datq_pat = datq_patch(ip,jp,ir,jr)
                     sumndvi(datq_pat) = sumndvi(datq_pat) + datp(ip,jp,ir,jr)
                     sumpix(datq_pat) = sumpix(datq_pat) + 1

                  enddo
               enddo

               do ipat = 2,nvegpat+1

                     datq_pat = nint(leaf_class(ir,jr,ipat))
                     if (sumpix(datq_pat) == 0.) then
                        veg_ndvif(ir,jr,ipat) =  0.05
                     else
                        veg_ndvif(ir,jr,ipat) =  max(.05,sumndvi(datq_pat)/ sumpix(datq_pat))
                     end if
               enddo

            endif
         enddo
      enddo

   endif

   if(nmiss.gt.0) then
      print*,'-----------------------------------------------------'
      print*,'Input surface characteristic data file processing:',iaction
      print*,'-----------------------------------------------------'
      print*,'  Input data blocks not found (data assumed to be ocean or default):'
      do nn=1,nmiss
         print*,trim(trim(fnmiss(nn)))
      enddo
      print*,'-----------------------------------------------------'
   endif



   return
end subroutine landuse_opqr

!********************************************************************

subroutine read_header(ofn,iblksizo,no,isbego,iwbego,offlat,offlon,deltallo  &
   ,ifield)

implicit none
integer :: iblksizo,no,isbego,iwbego,lb
real :: offlat,offlon,deltallo
character :: ofn*(*),title*256,ifield*(*)

lb = len_trim(ofn)
if (lb .le. 0) then
   print*,'| ',ifield,' input data prefix incorrect !'
   print*,'|  file prefix:',ofn(1:lb)
   print*,'====================================================='
   stop 'landuse-file'
endif

title = ofn(1:lb)//'HEADER'
lb = len_trim(title)

call rams_f_open(29,title(1:lb),'FORMATTED','OLD','READ',0)
read(29,*) iblksizo,no,isbego,iwbego,offlat,offlon
close(29)
deltallo = float(iblksizo) / float(no-1)

return
end

!*************************************************************************

subroutine fill_datp(n2,n3,no,iblksizo,isbego,iwbego  &
   ,platn,plonn,offlat,offlon,deltallo,ofn,iaction,nmiss,fnmiss)

use mem_mksfc

implicit none
character(len=*) :: ofn,iaction,fnmiss(*)
integer :: nmiss

integer :: n2,n3,no,iblksizo,isbego,iwbego,isoc,iwoc,iofr  &
   ,isocpt,isocpo,iwocph,iwocpt,iwocpo,lb,io,jo  &
   ,ir,jr,ip,jp,ind,nc3,nc2,j3d,j2d,j1d,ind1,ind2,io1,jo1  &
   ,ifile_max,jfile_max,ifile,jfile,missing,ptab,ptab0,idatp,nn

real :: rio,rjo,rno,platn,plonn,offlat,offlon  &
       ,glatp1,glonp1,deltallo,wio2,wjo2,wio1,wjo1

character :: title1*3,title2*4,title3*80

logical l1,l2

include 'interface.h'

! Compute number of files in input dataset that span all latitudes and
! longitudes on earth.  Allocate nump and numpind arrays to this size and 
! initialize to 0.
! Allocate ptable array.

rno = float(no)
ifile_max = 360 / iblksizo
jfile_max = 180 / iblksizo

if (iaction == 'ndvi') then
   allocate (dato(no,no))
else
   allocate (cdato(no,no))
endif

allocate (nump    (ifile_max,jfile_max)  &
         ,numpind (ifile_max,jfile_max)  &
         ,numpind1(ifile_max,jfile_max)  &
         ,numpind2(ifile_max,jfile_max)  &
         ,ptable  (npq*npq*n2*n3))

do jfile = 1,jfile_max
   do ifile = 1,ifile_max
      nump(ifile,jfile) = 0
      numpind(ifile,jfile) = 0
   enddo
enddo

! Get file index (ifile,jfile) within full dataset and count number of p 
! points (nump) that occur in each file


do jr = 1,n3
   do jp = 1,npq
      do ir = 1,n2
         do ip = 1,npq

            glatp1 = max(-89.9999,min(89.9999,glatp(ip,jp,ir,jr) - offlat))
            glonp1 = glonp(ip,jp,ir,jr) - offlon
            if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
            if (glonp1 .le. -180.) glonp1 = glonp1 + 360.

            ifile = int((glonp1 - float(iwbego)) / float(iblksizo)) + 1
            jfile = int((glatp1 - float(isbego)) / float(iblksizo)) + 1

            nump(ifile,jfile) = nump(ifile,jfile) + 1

!write(6,202) ip,jp,ir,jr,ifile,jfile,nump(ifile,jfile)  &
!               ,glatp1,glonp1
!202 format('at2',7i5,2f10.4)
            
         enddo
      enddo
   enddo
enddo

! Set up array index values for ptable array

ind = 1
do jfile = 1,jfile_max
   do ifile = 1,ifile_max
      numpind1(ifile,jfile) = ind
      numpind2(ifile,jfile) = ind
      ind = ind + nump(ifile,jfile)
   enddo
enddo

! Fill ptable array

nc3 = n2 * npq * npq
nc2 = npq * npq

do jr = 1,n3
   j3d = (jr - 1) * nc3
   do ir = 1,n2
      j2d = (ir - 1) * nc2
      do jp = 1,npq
         j1d = (jp - 1) * npq
         do ip = 1,npq

            glatp1 = max(-89.9999,min(89.9999,glatp(ip,jp,ir,jr) - offlat))
            glonp1 = glonp(ip,jp,ir,jr) - offlon
            if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
            if (glonp1 .le. -180.) glonp1 = glonp1 + 360.

            ifile = int((glonp1 - float(iwbego)) / float(iblksizo)) + 1
            jfile = int((glatp1 - float(isbego)) / float(iblksizo)) + 1

            ind = numpind2(ifile,jfile)
            ptable(ind) = j3d + j2d + j1d + ip
            numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1
            
!write(6,203) ip,jp,ir,jr,ifile,jfile,ind,ptable(ind)  &
!    ,numpind2(ifile,jfile),glonp1,glatp1,float(isbego),glatp1 - float(isbego)
!203 format('dat44',9i5,4f14.8)
            
         enddo
      enddo
   enddo
enddo

! Read files and extract data

do jfile = 1,jfile_max
   do ifile = 1,ifile_max
   
      ind1 = numpind1(ifile,jfile)
      ind2 = numpind2(ifile,jfile)
   
   
      if (ind2 .gt. ind1) then
         isoc = (jfile - 1) * iblksizo + isbego
         iwoc = (ifile - 1) * iblksizo + iwbego

! Construct filename

         isocpt = abs(isoc) / 10
         isocpo = abs(isoc) - isocpt*10
         iwocph = abs(iwoc) / 100
         iwocpt = (abs(iwoc) - iwocph * 100) / 10
         iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10
         
         
         if (isoc .ge. 0) then
            write(title1,'(2i1,a1)') isocpt,isocpo,'N'
         else
            write(title1,'(2i1,a1)') isocpt,isocpo,'S'
         endif
         
         
         if (iwoc .ge. 0) then
            write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
         else
            write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
         endif

         lb = len_trim(ofn)
         title3 = ofn(1:lb)//title1//title2
         if (iaction == 'ndvi') title3=trim(title3)//'.hdf'
         lb = len_trim(title3)
         
         inquire(file=title3(1:lb),exist=l1,opened=l2)

! Read file or set missing flag to 1
         if (l1) then
            missing = 0
           print*, 'getting file ',title3(1:lb)
           
            if (iaction == 'ndvi') then
               call read_hdf(no,no,dato,title3(1:lb))
            else
               call rams_c_open(title3(1:lb)//char(0),'rb'//char(0))
               call rams_c_read_char(4,no*no,cdato(1,1))
               call rams_c_close()
            endif
         else
            do nn=1,nmiss
               if(trim(title3(1:lb)) == trim(fnmiss(nn)) ) goto 302
            enddo
            nmiss=nmiss+1
            fnmiss(nmiss)=title3(1:lb)
302         continue
            missing = 1
         endif


         do ind = ind1,ind2-1
         

            ptab = ptable(ind)         
            ptab0 = ptab - 1
            jr = ptab0 / nc3 + 1
            j3d = (jr - 1) * nc3
            ir = (ptab0 - j3d) / nc2 + 1
            j2d = (ir - 1) * nc2
            jp = (ptab0 - j3d - j2d) / npq + 1
            j1d = (jp - 1) * npq
            ip = ptab - j3d - j2d - j1d

            glatp1 = max(-89.9999,min(89.9999,glatp(ip,jp,ir,jr) - offlat))
            glonp1 = glonp(ip,jp,ir,jr) - offlon
            if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
            if (glonp1 .le. -180.) glonp1 = glonp1 + 360.

            rio = (glonp1 - float(iwoc)) / deltallo + 1.
            !if( abs(glonp1 - float(iwoc)) < 1.e-5 ) rio = 1.
            rjo = (glatp1 - float(isoc)) / deltallo + 1.
            !if( abs(glatp1 - float(isoc)) < 1.e-5 ) rjo = 1.
            
!write(6,205) ip,jp,ir,jr,rio,rjo,glatp1,glonp1,glatp(ip,jp,ir,jr)  &
!              ,glonp(ip,jp,ir,jr),float(isoc),glatp1 - float(isoc)
!205 format('ip,jp ',4i5,8f14.8)

            if (rio .lt. .9 .or. rio .gt. rno+.1 .or.  &
                rjo .lt. .9 .or. rjo .gt. rno+.1) then
                write (unit=*,fmt='(a)') '-----------------------------------------------------'
                write (unit=*,fmt='(a)') 'Error: rio or rjo out of range:'
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'rio=      ',rio
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'rjo=      ',rjo
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'glonp1=   ',glonp1
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'glatp1=   ',glatp1
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'offlon=   ',offlon
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'offlat=   ',offlat
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'deltallo= ',deltallo
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'glonp=    ',glonp(ip,jp,ir,jr)
                write (unit=*,fmt='(a,1x,es13.6,1x)') 'glatp=    ',glatp(ip,jp,ir,jr)
                write (unit=*,fmt='(a,1x,i13,1x)')    'isoc=     ',isoc
                write (unit=*,fmt='(a,1x,i13,1x)')    'iwoc=     ',iwoc
                write (unit=*,fmt='(a,1x,i13,1x)')    'ip=       ',ip
                write (unit=*,fmt='(a,1x,i13,1x)')    'jp=       ',jp
                write (unit=*,fmt='(a,1x,i13,1x)')    'ir=       ',ir
                write (unit=*,fmt='(a,1x,i13,1x)')    'jr=       ',jr
                stop 45
            endif

            if (missing .eq. 0) then

               if (iaction .eq. 'veg' .or. iaction .eq. 'soil') then
                  io = nint(rio)
                  jo = nint(rjo)
                  idatp = ichar(cdato(io,jo))

                  datp(ip,jp,ir,jr) = float(mod(idatp+256,256))
                  !if(datp(ip,jp,ir,jr)> 94) then
                  !   print*,'big:',ip,jp,ir,jr,io,jo,datp(ip,jp,ir,jr)
                  !        datp(ip,jp,ir,jr)=94
                  ! endif
                  
                  !write(6,207) ind,ptab,ip,jp,ir,jr,datp(ip,jp,ir,jr),io,jo,ichar(cdato(io,jo))  &
                  !       ,title3(1:lb)
                  !207 format('datp405',6i5,f7.2,3i5,a40)               
                  
               elseif (iaction .eq. 'ndvi') then
                  io1 = max(1,min(int(rio),no-1))
                  jo1 = max(1,min(int(rjo),no-1))
                  wio2 = rio - float(io1)
                  wjo2 = rjo - float(jo1)
                  wio1 = 1. - wio2
                  wjo1 = 1. - wjo2
                  
                  if(jo1==0) then
                     print*,'jp0:',ip,jp,ir,jr,io1,jo1,rio,rjo,rno
                     print*,'jp0:',rjo, glatp1 , float(isoc)  &
                        ,glatp1-float(isoc), deltallo
                     stop
                  endif
                  datp(ip,jp,ir,jr) =  &
                     wio1 * (wjo1 * dato(io1  ,jo1  )   &
                          +  wjo2 * dato(io1  ,jo1+1))  &
                   + wio2 * (wjo1 * dato(io1+1,jo1  )   &
                          +  wjo2 * dato(io1+1,jo1+1))
               endif

            else

               datp(ip,jp,ir,jr) = 0.

            endif

            !write(6,209) ip,jp,ir,jr,glonp(ip,jp,ir,jr),glatp(ip,jp,ir,jr)  &
            !   ,datp(ip,jp,ir,jr),rio,rjo,no,rno,iaction,io
            !209 format('fillp ',4i5,3f10.4,

         enddo
      endif
   enddo
enddo

if (iaction == 'ndvi') then
   deallocate (dato)
else
   deallocate (cdato)
endif

deallocate(nump,numpind,numpind1,numpind2,ptable)

return
end
!****************************************************************************

subroutine read_hdf(nxll,nyll,data1,file_name)

   implicit none

   integer :: nxll,nyll
   real, dimension(nxll,nyll) :: data1
   character :: file_name*(*)
#if USE_HDF4
   ! Parameter declaration.

   integer, parameter :: dfacc_read = 1 ,dfnt_float32 = 5
                       
   ! Function declaration

   integer :: sfstart, sfselect, sfrdata, sfendacc, sfend

   ! Variable declaration

   integer :: sd_id, sds_id, sds_index, status
   integer, dimension(2) :: start, edges, stride

   ! End of variable declaration

   ! Open the file and initialize the SD interface
   sd_id = sfstart(file_name, dfacc_read)

   !  Select the first data set.
   sds_index = 0
   sds_id = sfselect(sd_id, sds_index)

   ! Set elements of the array start to 0, elements of the
   ! array edges to SDS dimensions, and elements of the array stride
   ! to 1 to read the entire data.

   start(1) = 0
   start(2) = 0
   edges(1) = nxll
   edges(2) = nyll
   stride(1) = 1
   stride(2) = 1

   ! Read entire data into data array.  Note that sfrdata is used
   ! to read the numeric data.
   status = sfrdata(sds_id, start, stride, edges, data1)

   !-srf
   !print*,'data1=',data1
   !stop 32323
   !-srf



   ! Terminate access to the data set.
   status = sfendacc(sds_id)

   ! Terminate access to the SD interface and close the file.
   status = sfend(sd_id)

#else
   call abort_run('You must compile BRAMS with HDF4 libraries to read NDVI data'           &
                 ,'read_hdf','landuse_input.F90')
#endif

   return
end subroutine read_hdf
