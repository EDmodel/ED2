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
   character(len=256) :: h5name
   real :: deltax

   integer :: iblksiz,no,isbeg,iwbeg
   real :: offlat,offlon,deltall,deltallo_min

   ! Find highest resolution dataset among all that are to be used

   deltallo_min = 1.e20
   if (ivegtflg == 1) then
      call read_header(ivegtfn,iblksiz,no,isbeg  &
         ,iwbeg,offlat,offlon,deltall,'veg',h5name)
      deltallo_min = min(deltallo_min,deltall)
   endif

   if (isoilflg == 1) then

      call read_header(isoilfn,iblksiz,no,isbeg  &
         ,iwbeg,offlat,offlon,deltall,'soil',h5name)
      deltallo_min = min(deltallo_min,deltall)
   endif

   if (ndviflg == 1) then
      call read_header(ndvifn,iblksiz,no,isbeg  &
           ,iwbeg,offlat,offlon,deltall,'ndvi',h5name)
      deltallo_min = min(deltallo_min,deltall)
   endif

!!   if (ndviflg == 1) then
!!      call read_header(ndvifn,iblksiz,no,isbeg  &
!!         ,iwbeg,offlat,offlon,deltall,'ndvi')
!!      deltallo_min = min(deltallo_min,deltall)
!!   endif

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
subroutine landuse_opqr(n2,n3,mzg,npat,nvegpat,ivegtflg,ivegtfn,isoilflg,isoilfn           &
                       ,ndviflg,ndvifn,cndvifil,iaction,platn,plonn,soil_text,patch_area   &
                       ,leaf_class,veg_ndvif)

   use mem_mksfc
   use rconstants
   use mem_leaf , only : nslcon          & ! intent(in)
                       , isfcl           ! ! intent(in)
   use leaf_coms, only : min_patch_area  & ! intent(in)
                       , nstyp           & ! intent(in)
                       , nvtyp           & ! intent(in)
                       , nvtyp_teb       ! ! intent(in)
   use io_params, only : iuselai         ! ! intent(in)
   use grid_dims, only : str_len         ! ! intent(in)
   implicit none
   !----- Local parameters. ---------------------------------------------------------------!
   integer                                           , parameter     :: maxdatq =   32
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                           , intent(in)    :: n2
   integer                                           , intent(in)    :: n3
   integer                                           , intent(in)    :: mzg
   integer                                           , intent(in)    :: npat
   integer                                           , intent(in)    :: nvegpat
   character(len=*)                                  , intent(in)    :: iaction
   character(len=*)                                  , intent(in)    :: ivegtfn
   character(len=*)                                  , intent(in)    :: isoilfn
   character(len=*)                                  , intent(in)    :: ndvifn
   character(len=*)                                  , intent(in)    :: cndvifil
   integer                                           , intent(in)    :: ivegtflg
   integer                                           , intent(in)    :: isoilflg
   integer                                           , intent(in)    :: ndviflg
   real                                              , intent(in)    :: platn
   real                                              , intent(in)    :: plonn
   real, dimension(mzg,n2,n3,npat)                   , intent(inout) :: soil_text
   real, dimension(n2,n3,npat)                       , intent(inout) :: patch_area
   real, dimension(n2,n3,npat)                       , intent(inout) :: leaf_class
   real, dimension(n2,n3,npat)                       , intent(inout) :: veg_ndvif
   !----- Local variables. ----------------------------------------------------------------!
   character(len=str_len), dimension(maxndvidata)                    :: fnmiss
   character(len=str_len)                                            :: h5name
   integer               , dimension(0:maxdatq)                      :: sumpix
   integer               , dimension(0:maxdatq,2)                    :: ngrdpix
   integer               , dimension(0:maxdatq,nstyp)                :: datq_soil
   integer               , dimension(0:maxdatq,nstyp)                :: soil_tab
   real                  , dimension(0:maxdatq)                      :: datq_ndvi
   real                  , dimension(0:maxdatq)                      :: sumndvi
   integer                                                           :: ing_prt
   integer                                                           :: nn
   integer                                                           :: nmiss
   integer                                                           :: i
   integer                                                           :: j
   integer                                                           :: datq
   integer                                                           :: lsp
   integer                                                           :: datq_pat
   integer                                                           :: datsoil
   integer                                                           :: soil_count
   integer                                                           :: nlev_soil
   integer                                                           :: jr
   integer                                                           :: jp
   integer                                                           :: ir
   integer                                                           :: ip
   integer                                                           :: iqv
   integer                                                           :: ing
   integer                                                           :: maxdq
   integer                                                           :: jng
   integer                                                           :: jng1
   integer                                                           :: ngstor1
   integer                                                           :: ngstor2
   integer                                                           :: npatpixs
   integer                                                           :: nwat
   integer                                                           :: ipat
   integer                                                           :: idatq
   integer                                                           :: isoil
   integer                                                           :: jsoil
   integer                                                           :: k
   integer                                                           :: no_veg
   integer                                                           :: iblksizo_veg
   integer                                                           :: isbego_veg
   integer                                                           :: iwbego_veg
   integer                                                           :: no_soil
   integer                                                           :: iblksizo_soil
   integer                                                           :: isbego_soil
   integer                                                           :: iwbego_soil
   integer                                                           :: no_ndvi
   integer                                                           :: iblksizo_ndvi
   integer                                                           :: isbego_ndvi
   integer                                                           :: iwbego_ndvi
   real                                                              :: checksum
   real                                                              :: parea_tot
   real                                                              :: deltallo_veg
   real                                                              :: deltallo_soil
   real                                                              :: deltallo_ndvi
   real                                                              :: offlat_veg
   real                                                              :: offlat_soil
   real                                                              :: offlat_ndvi
   real                                                              :: offlon_veg
   real                                                              :: offlon_soil
   real                                                              :: offlon_ndvi
   real                                                              :: xp
   real                                                              :: yp
   real                                                              :: fracwat
   real                                                              :: plpp
   !---------------------------------------------------------------------------------------!



   !------ Initialise some counters. ------------------------------------------------------!
   nmiss    = 0
   checksum = 0
   !---------------------------------------------------------------------------------------!



   select case (trim(iaction))
   case ('veg')

      !------------------------------------------------------------------------------------!
      !     Read in the header and data files.                                             !
      !------------------------------------------------------------------------------------!
      call read_header(ivegtfn,iblksizo_veg,no_veg,isbego_veg,iwbego_veg,offlat_veg        &
                      ,offlon_veg,deltallo_veg,'veg',h5name)
         
      call fill_datp(n2,n3,no_veg,iblksizo_veg,isbego_veg,iwbego_veg,platn,plonn           &
                    ,offlat_veg,offlon_veg,deltallo_veg,ivegtfn,iaction,nmiss,fnmiss,h5name)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Carry out the first translation of the input DATP values into a condensed set  !
      ! called DATQ_patch.  The range of DATQ_patch values represents the total variety of !
      ! land surface conditions (patches) to be allowed for the present simulation, and    !
      ! may be a broader class than the LEAF-2 vegetation classes for which all the        !
      ! vegetation physical parameters are defined.  For example, two different DATQ_patch !
      ! classes may be mapped to the same LEAF-2 vegetation class, but be initialized with !
      ! different soil moistures or soil types, and therefore require different patches.   !
      !     Fill datq_patch (patch class) values from input landuse dataset.  Currently,   !
      ! this data serves as the primary criterion for defining patches if running LEAF,    !
      ! but not if running ED, in which case soil will overwrite the patch areas.          !
      !------------------------------------------------------------------------------------!
      jvegloop: do jr = 1,n3
         ivegloop: do ir = 1,n2

            do iqv = 0,maxdatq
               ngrdpix(iqv,1) = 0     ! Initialise counter for datq pixels
               ngrdpix(iqv,2) = iqv   ! Initialise array of consecutive datq values 
            enddo

            !------------------------------------------------------------------------------!
            !     Count the pixels for this grid point.                                    !
            !------------------------------------------------------------------------------!
            do jp = 1,npq
               do ip = 1,npq

                  call leaf_datp_datq(datp(ip,jp,ir,jr),datq_patch(ip,jp,ir,jr))
                  datq = datq_patch(ip,jp,ir,jr)
                  ngrdpix(datq,1) = ngrdpix(datq,1) + 1
               end do
            end do
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !    Sort values of ngrdpix by prevalence for non-water patches (datq >= 2).   !
            !------------------------------------------------------------------------------!
            do ing = 2,maxdatq
               maxdq = -1
               do jng = ing,maxdatq
                  if (ngrdpix(jng,1) .gt. maxdq) then
                     jng1 = jng
                     maxdq = ngrdpix(jng,1)
                  endif
               end do
               ngstor1 = ngrdpix(ing,1)
               ngstor2 = ngrdpix(ing,2)
               ngrdpix(ing,1) = ngrdpix(jng1,1)
               ngrdpix(ing,2) = ngrdpix(jng1,2)
               ngrdpix(jng1,1) = ngstor1
               ngrdpix(jng1,2) = ngstor2
            end do
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !    Fill patches numbered 2 through nvegpat+1 with the nvegpat most prevalent !
            ! nonwater landuse types.  Count pixels for these patches for normalization of !
            ! total grid cell area.                                                        !
            !------------------------------------------------------------------------------!
            npatpixs = 1
            nwat     = ngrdpix(0,1) + ngrdpix(1,1)
            if (nwat < npq*npq) then
               npatpixs = 0
               do ipat  = 2,nvegpat+1

                  datq                   = ngrdpix(ipat,2)
                  leaf_class(ir,jr,ipat) = float(datq)
                  npatpixs               = npatpixs + ngrdpix(ipat,1)
               end do
            end if

            fracwat             = float(nwat) / float(npq * npq)
            plpp                = (1. - fracwat) / float(npatpixs)
            patch_area(ir,jr,1) = max(min_patch_area,fracwat)

            if (ngrdpix(0,1) >= ngrdpix(1,1)) then
               leaf_class(ir,jr,1) = 0.
            else
               leaf_class(ir,jr,1) = 1.
            end if

            do ipat = 2,nvegpat+1
               patch_area(ir,jr,ipat) = plpp * float(ngrdpix(ipat,1))
               if (patch_area(ir,jr,ipat) < min_patch_area) patch_area(ir,jr,ipat) = 0.
            end do

            !----- Rescale the patch areas, eliminating those tiny patches. ---------------!
            parea_tot = 0.
            do ipat = 1,nvegpat+1
               parea_tot = parea_tot + patch_area(ir,jr,ipat)
            end do
            do ipat = 1,nvegpat+1
               patch_area(ir,jr,ipat) = patch_area(ir,jr,ipat) / parea_tot
            end do
         end do ivegloop
      end do jvegloop
      !------------------------------------------------------------------------------------!

   case ('soil')


      !------------------------------------------------------------------------------------!
      !     Read in the header and data files.                                             !
      !------------------------------------------------------------------------------------!
      call read_header(isoilfn,iblksizo_soil,no_soil,isbego_soil,iwbego_soil,offlat_soil   &
                      ,offlon_soil,deltallo_soil,'soil',h5name)

      call fill_datp(n2,n3,no_soil,iblksizo_soil,isbego_soil,iwbego_soil,platn,plonn       &
                    ,offlat_soil,offlon_soil,deltallo_soil,isoilfn,iaction,nmiss,fnmiss    &
                    ,h5name)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      jsoilloop: do jr = 1,n3
         isoilloop: do ir = 1,n2
            if (patch_area(ir,jr,1) < 1.0) then

               !---------------------------------------------------------------------------!
               !    Here we initialise both ngrdpix and datq_soil, because we will assign  !
               ! patches differently depending on the surface model (LEAF or ED).          !
               !---------------------------------------------------------------------------!
               do idatq = 0,maxdatq
                  ngrdpix(idatq,1) = 0     ! Initialise counter for datq pixels
                  ngrdpix(idatq,2) = idatq ! Initialise array of consecutive datq values 
                  do isoil = 1,nstyp
                     datq_soil(idatq,isoil) = 0  ! Initialise counter for datq pixels
                  end do
               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Count the pixels for this grid point.                                 !
               !---------------------------------------------------------------------------!
               do jp = 1,npq
                  do ip = 1,npq

                     !---------------------------------------------------------------------!
                     !      Fill datq_soil values as secondary criterion if running LEAF,  !
                     ! and fill ngrdpix as a primary criterion if running ED.              !
                     !---------------------------------------------------------------------!
                     datsoil  = nint(datp(ip,jp,ir,jr))

                     datq_pat = datq_patch(ip,jp,ir,jr)

                     datq_soil(datq_pat,datsoil) = datq_soil(datq_pat,datsoil) + 1
                     ngrdpix(datsoil,1)          = ngrdpix(datsoil,1) + 1
                  end do
               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Here we will do different things depending on whether we are running  !
               ! LEAF-3 or ED-2.                                                           !
               !---------------------------------------------------------------------------!
               select case (isfcl)
               case (1,2)
                  !------------------------------------------------------------------------!
                  !     LEAF-3.   In this case we select the commonest soil class for each !
                  ! vegetation type, as they may be correlated.                            !
                  !------------------------------------------------------------------------!
                  do ipat = 2,nvegpat+1

                     if (patch_area(ir,jr,ipat) >= min_patch_area) then

                        datq_pat = nint(leaf_class(ir,jr,ipat))

                        !------------------------------------------------------------------!
                        !      Find isoil value for which datq_soil(datq_pat,isoil) is a   !
                        ! maximum.                                                         !
                        !------------------------------------------------------------------!
                        soil_count = 0
                        do isoil = 1,nstyp
                           if (datq_soil(datq_pat,isoil) > soil_count) then
                              soil_count = datq_soil(datq_pat,isoil)
                              jsoil      = isoil
                           end if
                        end do

                        checksum = checksum + float(jsoil)
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !      For now, assume single level of input soil data (e.g., FAO) !
                        ! and fill all soil levels with this value.                        !
                        !------------------------------------------------------------------!
                        do k = 1,mzg
                           if (jsoil /= 0) then
                              soil_text(k,ir,jr,ipat) = float(jsoil)
                           else
                              soil_text(k,ir,jr,ipat) = float(nslcon)
                           end if
                        end do
                        !------------------------------------------------------------------!
                     end if
                  end do
                  !------------------------------------------------------------------------!

               case (5)
                  !------------------------------------------------------------------------!
                  !     ED-2 run.   In this case soil texture will be the primary          !
                  ! criterion to determine the patches.   We must organise the counter     !
                  ! array so it goes from the commonest soil class to the rarest one.  We  !
                  ! must also redefine the patch area (which will be the site area).       !
                  !------------------------------------------------------------------------!





                  !------------------------------------------------------------------------!
                  !    Sort values of ngrdpix by prevalence for non-water patches          !
                  ! (datq >= 1).                                                           !
                  !------------------------------------------------------------------------!
                  do ing = 1,maxdatq
                     maxdq = -1
                     do jng = ing,maxdatq
                        if (ngrdpix(jng,1) .gt. maxdq) then
                           jng1 = jng
                           maxdq = ngrdpix(jng,1)
                        endif
                     end do
                     ngstor1 = ngrdpix(ing,1)
                     ngstor2 = ngrdpix(ing,2)
                     ngrdpix(ing,1) = ngrdpix(jng1,1)
                     ngrdpix(ing,2) = ngrdpix(jng1,2)
                     ngrdpix(jng1,1) = ngstor1
                     ngrdpix(jng1,2) = ngstor2
                  end do
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !    Fill patches numbered 2 through nvegpat+1 with the nvegpat most     !
                  ! prevalent nonwater landuse types.  Count pixels for these patches for  !
                  ! normalization of total grid cell area.                                 !
                  !------------------------------------------------------------------------!
                  npatpixs = 1
                  nwat     = ngrdpix(0,1)
                  if (nwat < npq*npq) then
                     npatpixs = 0
                     do ipat  = 2,nvegpat+1

                        datq                   = ngrdpix(ipat-1,2)
                        !------------------------------------------------------------------!
                        !      For now, assume single level of input soil data (e.g., FAO) !
                        ! and fill all soil levels with this value.                        !
                        !------------------------------------------------------------------!
                        do k=1,mzg
                           soil_text(k,ir,jr,ipat) = float(datq)
                        end do
                        !------------------------------------------------------------------!

                        npatpixs               = npatpixs + ngrdpix(ipat-1,1)
                     end do
                     
                     !---------------------------------------------------------------------!
                     !     PLPP will scale the soil class so the sum of water + land       !
                     ! patches adds up to one.                                             !
                     !---------------------------------------------------------------------!
                     plpp = (1. - patch_area(ir,jr,1)) / real(npatpixs)
                     
                     !---------------------------------------------------------------------!
                     !     Find the area of each patch, and discard the leaf class         !
                     ! information by copying the second patch to the least common ones.   !
                     !---------------------------------------------------------------------!
                     do ipat=2,nvegpat+1
                        leaf_class(ir,jr,ipat) = leaf_class(ir,jr,2)
                        patch_area(ir,jr,ipat) = plpp * real(ngrdpix(ipat-1,1))
                        if (patch_area(ir,jr,ipat) < min_patch_area) then
                           patch_area(ir,jr,ipat) = 0.
                        end if
                     end do
                     !---------------------------------------------------------------------!

                  else
                     !---------------------------------------------------------------------!
                     !     The soil dataset indicated that all pixels were water, even     !
                     ! though the vegetation dataset had land.  If this is the case,       !
                     ! assign only one leaf patch/ed site with the default soil type and   !
                     ! set the others areas to zero.                                       !
                     !---------------------------------------------------------------------!
                     do k=1,mzg
                        soil_text(k,ir,jr,2) = nslcon
                     end do
                     patch_area(ir,jr,2) = 1. - patch_area(ir,jr,1)
                     do ipat = 3, nvegpat+1
                        do k=1, mzg
                           soil_text(k,ir,jr,ipat) = nslcon
                        end do
                        patch_area(ir,jr,ipat) = 0.
                     end do
                  end if
                  !------------------------------------------------------------------------!


                  !----- Rescale the patch areas, eliminating those tiny patches. ---------!
                  parea_tot = 0.
                  do ipat = 1,nvegpat+1
                     parea_tot = parea_tot + patch_area(ir,jr,ipat)
                  end do
                  do ipat = 1,nvegpat+1
                     patch_area(ir,jr,ipat) = patch_area(ir,jr,ipat) / parea_tot
                  end do
                  !------------------------------------------------------------------------!

               end select
            end if
         end do isoilloop
      end do jsoilloop



   case ('ndvi')

      !------------------------------------------------------------------------------------!
      !     Read in the header and NDVI data files.                                        !
      !------------------------------------------------------------------------------------!
      call read_header(ndvifn,iblksizo_ndvi,no_ndvi,isbego_ndvi,iwbego_ndvi,offlat_ndvi    &
                      ,offlon_ndvi,deltallo_ndvi,'ndvi',h5name)
      call fill_datp  (n2,n3,no_ndvi,iblksizo_ndvi,isbego_ndvi,iwbego_ndvi,platn,plonn     &
                      ,offlat_ndvi,offlon_ndvi,deltallo_ndvi,cndvifil,iaction,nmiss,fnmiss &
                      ,h5name)
      !------------------------------------------------------------------------------------!


      jndviloop: do jr = 1,n3
         indviloop: do ir = 1,n2
            if (patch_area(ir,jr,1) < 1.0) then

               do idatq = 0,maxdatq
                  sumndvi(idatq) = 0.   ! Initialise ndvi sum
                  sumpix(idatq)  = 0    ! Initialise ndvi pixel count
               end do

               do jp = 1,npq
                  do ip = 1,npq

                     !---------------------------------------------------------------------!
                     !     Fill datq_ndvi values to compute average value for each datq    !
                     ! class.                                                              !
                     !---------------------------------------------------------------------!
                     datq_pat          = datq_patch(ip,jp,ir,jr)
                     sumndvi(datq_pat) = sumndvi(datq_pat) + datp(ip,jp,ir,jr)
                     sumpix(datq_pat)  = sumpix(datq_pat) + 1
                  end do
               end do


               !---------------------------------------------------------------------------!
               !     Fill in the NDVI (or LAI) arrays using the average NDVI (or LAI) for  !
               ! each class.                                                               !
               !---------------------------------------------------------------------------!
               do ipat = 2,nvegpat+1
                  datq_pat = nint(leaf_class(ir,jr,ipat))

                  select case (iuselai)
                  case (1)
                     !----- LAI data. -----------------------------------------------------!
                     if (sumpix(datq_pat) == 0.) then
                        veg_ndvif(ir,jr,ipat) =  0.0
                     else
                        veg_ndvif(ir,jr,ipat) =  sumndvi(datq_pat)/ sumpix(datq_pat)
                     end if
                  case default
                     !----- NDVI data, make sure it's never zero (singularity). -----------!
                     if (sumpix(datq_pat) == 0.) then
                        veg_ndvif(ir,jr,ipat) =  0.05
                     else
                        veg_ndvif(ir,jr,ipat) =  max( 0.05                                 &
                                                    , sumndvi(datq_pat)/ sumpix(datq_pat))
                     end if
                  end select
               end do
            end if
         end do indviloop
      end do jndviloop

   end select

   if (nmiss > 0) then
      write (unit=*,fmt='(a)')      '-----------------------------------------------------'
      write (unit=*,fmt='(a,1x,a)') '  The following blocks for ',trim(iaction)            &
                                    ,' were missing.'
      do nn=1,nmiss
         write (unit=*,fmt='(a,1x,a)')  '  - ',trim(fnmiss(nn))
      end do
      write (unit=*,fmt='(a)')      '  Data were assumed ocean or default in these areas.'
      write (unit=*,fmt='(a)')      '-----------------------------------------------------'
      write (unit=*,fmt='(a)')      ' '
   end if
   return
end subroutine landuse_opqr
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_header(ofn,iblksizo,no,isbego,iwbego,offlat,offlon,deltallo  &
   ,ifield,h5name)

implicit none
integer :: iblksizo,no,isbego,iwbego,lb
real :: offlat,offlon,deltallo
character :: ofn*(*),title*256,ifield*(*)
character :: h5name*(*)


lb = len_trim(ofn)
if (lb .le. 0) then
   print*,'| ',ifield,' input data prefix incorrect !'
   print*,'|  file prefix:',ofn(1:lb)
   print*,'====================================================='
   stop 'landuse-file'
endif

title = ofn(1:lb)//'HEADER'

lb = len_trim(title)

if (trim(ifield).eq.'ndvi') then
   call rams_f_open(29,title(1:lb),'FORMATTED','OLD','READ',0)
   read(29,*,end=1) iblksizo,no,isbego,iwbego,offlat,offlon,h5name
1  continue
   close(29)
   deltallo = float(iblksizo) / float(no-1)

   print*, 'read_header1 ',ifield
   print*, 'read_header2 ',iblksizo,no,isbego,iwbego,offlat,offlon,deltallo
else
   
   call rams_f_open(29,title(1:lb),'FORMATTED','OLD','READ',0)
   read(29,*) iblksizo,no,isbego,iwbego,offlat,offlon
   close(29)
   deltallo = float(iblksizo) / float(no-1)
end if

return
end subroutine read_header

!*************************************************************************

subroutine fill_datp(n2,n3,no,iblksizo,isbego,iwbego  &
   ,platn,plonn,offlat,offlon,deltallo,ofn,iaction,nmiss,fnmiss,h5name)

use mem_mksfc
use io_params, only : iuselai

#if USE_HDF5
use hdf5_utils
#endif

implicit none
character(len=*) :: ofn,iaction,fnmiss(*)
character(len=*) :: h5name
integer :: nmiss

integer :: n2,n3,no,iblksizo,isbego,iwbego,isoc,iwoc,iofr  &
   ,isocpt,isocpo,iwocph,iwocpt,iwocpo,lb,io,jo  &
   ,ir,jr,ip,jp,ind,nc3,nc2,j3d,j2d,j1d,ind1,ind2,io1,jo1  &
   ,ifile_max,jfile_max,ifile,jfile,missing,ptab,ptab0,idatp,nn

real :: rio,rjo,rno,platn,plonn,offlat,offlon  &
       ,glatp1,glonp1,deltallo,wio2,wjo2,wio1,wjo1

character :: title1*3,title2*4,title3*256

logical l1,l2,h5
integer :: ndims,idims(4),ii,jj

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
   allocate (cdato(no,no),idato(no,no))
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
         if (iaction == 'ndvi') then
            select case (iuselai)
            case (0)
               h5 = .false.
               title3=trim(title3)//'.hdf'
            case (1)
               h5 = .true.
#if USE_HDF5
               title3=trim(title3)//'.h5'
#else
               call abort_run('You must compile BRAMS with HDF5 to use LAI data' &
                             ,'fill_datp','landuse_input.F90')
#endif
            end select
         else
            h5 = .false.
         end if

         lb = len_trim(title3)
         inquire(file=title3(1:lb),exist=l1,opened=l2)

#if USE_HDF5
         ! If file not found, then check for an hdf5 file.
         if (.not. l1) then
            if (iaction == 'ndvi') then
               title3=title3(1:lb-3)//'h5'
            else
               title3=trim(title3)//'.h5'
            end if

            inquire(file=trim(title3),exist=l1,opened=l2)
            h5 = l1
         end if
#endif
! Read file or set missing flag to 1
     
         if (l1) then
            missing = 0
           !print*, 'getting file ',title3(1:lb),ir,jr,ip,jp
           print*, 'getting file ',trim(title3)
           
            if (iaction == 'ndvi') then
               if (h5) then
#if USE_HDF5
                  call shdf5_open_f(title3,'R')
                  ndims=2 ; idims(1)=no ; idims(2)=no
                  call shdf5_irec_f(ndims,idims,trim(h5name),rvara=dato)
                  call shdf5_close_f()
#else
                  call abort_run('You must compile BRAMS with HDF5 to use HDF5 data' &
                                ,'fill_datp','landuse_input.F90')
#endif
               elseif (iuselai == 0) then
                  call read_hdf(no,no,dato,title3(1:lb))
               else
                  call abort_run('LAI data must be in HDF5 format...' &
                                ,'fill_datp','landuse_input.F90')
               end if

            else
#if USE_HDF5
               if (h5) then
                  call shdf5_open_f(title3,'R')
                  ndims=2 ; idims(1)=no ; idims(2)=no
                  call shdf5_irec_f(ndims,idims,trim(h5name),ivara=idato)
                  call shdf5_close_f()
               else
                  call rams_c_open(title3(1:lb)//char(0),'rb'//char(0))
                  call rams_c_read_char(4,no*no,cdato(1,1))
                  call rams_c_close()
                  do jj = 1,no
                     do ii = 1,no
                        idato(ii,jj)=ichar(cdato(ii,jj))
                     enddo
                  enddo
               endif
#else
               call rams_c_open(title3(1:lb)//char(0),'rb'//char(0))
               call rams_c_read_char(4,no*no,cdato(1,1))
               call rams_c_close()
               do jj = 1,no
                  do ii = 1,no
                     idato(ii,jj)=ichar(cdato(ii,jj))
                  end do
               end do
#endif
            end if
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
            if (glonp1 >=  180.) then
               glonp1 = glonp1 - 360.
            elseif (glonp1 <= -180.) then
               glonp1 = glonp1 + 360.
            end if
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
   deallocate (cdato,idato)
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
