!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine ndvi_read_dataheader(ifm)

use mem_mksfc
use io_params

implicit none
integer :: ifm

integer :: itime,nc
character(len=256) :: flnm,line,line2
character(len=1) :: dummy
logical :: there
integer, external :: lastslash

! Read header file for all ndvidata files (all times and locations).  The header
! file contains:
! first line:       geographic block size (degrees), file size, geographic
!                   starting point for the dataset (south pole), and offsets
! second line:      number of data times (NDTIM)
! next NDTIM lines: file prefix, year, month, day, and hour for one data time
! last 3 lines:     comments describing the first, second, and following lines
!                   above

! Construct header file name

flnm=trim(ndvifn(ifm))//'HEADER'

inquire(file=flnm,exist=there)
if (.not.there) then
   print*,'ndvifn data header file for grid ',ifm,' not there.'
   stop 'ndvi_read_fileheader-1'
endif

! Read this header file

call rams_f_open(25,flnm,'FORMATTED','OLD','READ',0)
rewind 25

! read number of data times in dataset

read(25,*) dummy
read(25,*) nvndvif(ifm)

if (nvndvif(ifm) <= 0) then
   print*, 'No ndvi input files found with specified prefix or incorrect header'
   close(25)
   stop 'ndvi_read_fileheader-2'
endif

! read prefix list and times 

do itime = 1,nvndvif(ifm)
   read(25,'(A80)') line
   call char_strip_var(line,flnm,line2)
   read (line2,*) iyearvn(itime,ifm),imonthvn(itime,ifm)  &
      ,idatevn(itime,ifm),ihourvn(itime,ifm)

   nc=lastslash(ndvifn(ifm))
   vndvifil(itime,ifm)=ndvifn(ifm)(1:nc)//trim(flnm)

enddo

close(25)

return
end

!******************************************************************************

subroutine ndvinest(ifm,ivtime)

use mem_mksfc
use mem_grid
use mem_leaf
use io_params

implicit none

integer :: ifm,icm,ivtime,mynum,i,j,k,ic,jc,ipat

icm = nxtnest(ifm)

! Initialize SEATP and SEATF in subroutine ndviinit

call ndviinit(nnxp(ifm),nnyp(ifm),npatch,ifm  &
   ,sfcfile_p(ifm)%veg_ndvif)

if (icm >= 1 .and. ndviflg(ifm) == 0) then

   ! Assign NDVIF from coarser grid

   do ipat = 2,npatch
      do k = 1,nzg
         do j = 1,nnyp(ifm)
            do i = 1,nnxp(ifm)
               ic = ipm(i,ifm)
               jc = jpm(j,ifm)

               sfcfile_p(ifm)%veg_ndvif(i,j,ipat)  =  &
                  sfcfile_p(icm)%veg_ndvif(ic,jc,ipat)

            enddo
         enddo
      enddo
   enddo
   
   nvndvif(ifm) = nvndvif(icm)
   iyearvn (1:nvndvif(ifm),ifm) = iyearvn (1:nvndvif(ifm),icm)
   imonthvn(1:nvndvif(ifm),ifm) = imonthvn(1:nvndvif(ifm),icm)
   idatevn (1:nvndvif(ifm),ifm) = idatevn (1:nvndvif(ifm),icm)
   ihourvn (1:nvndvif(ifm),ifm) = ihourvn (1:nvndvif(ifm),icm)

elseif (ndviflg(ifm) == 1) then

   ! Assign NDVIF from standard dataset:

   call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat  &
      ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm) &
      ,ndviflg(ifm),ndvifn(ifm),vndvifil(ivtime,ifm)  &
      ,'ndvi',platn(ifm),plonn(ifm)        &
      ,sfcfile_p(ifm)%soil_color  &
      ,sfcfile_p(ifm)%soil_text  &
      ,sfcfile_p(ifm)%patch_area   &
      ,sfcfile_p(ifm)%leaf_class   &
      ,sfcfile_p(ifm)%veg_ndvif)
else

   iyearvn (1,ifm) = iyeara ; imonthvn(1,ifm) = imontha
   idatevn (1,ifm) = idatea ; ihourvn (1,ifm) = ihoura

endif


! If desired, override current values of NDVIP and NDVIF in user-specified
! changes to subroutine ndviinit_user in the file ruser.f.

call ndviinit_user(nnxp(ifm),nnyp(ifm),npatch,ifm  &
   ,sfcfile_p(ifm)%veg_ndvif)

return
end

! ****************************************************************************

subroutine ndviinit(n2,n3,npat,ifm,veg_ndvif)
implicit none
integer :: n2,n3,npat,ifm,i,j,ipat
real, dimension(n2,n3,npat) :: veg_ndvif

! Fill the veg_ndvif array with a default value of .5.  This default is 
! used only when a standard RAMS ndvi dataset is not used and when no 
! overrides to ndvi values are defined in subroutine ndviinit_user in the
! file ruser.f.

do j = 1,n3
   do i = 1,n2
      veg_ndvif(i,j,1) = 0.
      veg_ndvif(i,j,2) = .5

      do ipat = 3,npat
         veg_ndvif(i,j,ipat) = veg_ndvif(i,j,2)
      enddo

   enddo
enddo
return
end

!****************************************************************************

subroutine ndvi_write(ifm,ivt)

use mem_mksfc
use mem_grid
use io_params

implicit none
integer :: ifm,ivt,ip

real :: glatr,glonr
character(len=256) :: flnm
character(len=2) :: cgrid
real(kind=8) :: zero

! Write ndvi data to ndvi file for one grid and one time

write(cgrid,'(a1,i1)') 'g',ifm


zero=0.0

call makefnam(flnm,ndvifpfx,zero,iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
      ,idatevn(ivt,ifm),ihourvn (ivt,ifm)*10000,'N',cgrid,'vfm')

print*,flnm

call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

call rams_f_open(25,flnm,'FORMATTED','REPLACE','WRITE',1)
rewind 25

write(25,99) 999999,2
99   format(2i8)

write(25,100) iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
             ,idatevn(ivt,ifm),ihourvn (ivt,ifm)
100  format(1x,i4.4,2(1x,i2.2),1x,i4.4)

write(25,101) nnxp(ifm),nnyp(ifm),npatch
101  format(4i5)

write(25,102) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
     ,glatr,glonr
102  format(6f16.5)

do ip = 1,npatch
   call vforec(25,sfcfile_p(ifm)%veg_ndvif(:,:,ip),nnxyp(ifm),24,scrx,'LIN')
enddo

close(25)

return
end
