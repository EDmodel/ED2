!==========================================================================================!
!==========================================================================================!
!    This subroutine will determine when should we copy the boundary conditions to the     !
! edges.  This should be done only when we are the full domain edge, not at the node edge, !
! because this is going to be used by routines like the turbulence.  In this case, the     !
! variable should not be copied here; instead, they should be send through MPI, and the    !
! variable table should include the mpt1 flag (lateral boundary condition).                !
!------------------------------------------------------------------------------------------!
subroutine leaf3_bcond(m2,m3,mzg,mzs,npat,ia,iz,ja,jz,jdim,ibcon,soil_water,sfcwater_mass  &
                      ,soil_energy,sfcwater_energy,soil_color,soil_text,psibar_10d         &
                      ,sfcwater_depth,ustar,tstar,rstar,cstar,zeta,ribulk,veg_albedo       &
                      ,veg_fracarea,veg_lai,veg_tai,veg_rough,veg_height,veg_displace      &
                      ,patch_area,patch_rough,patch_wetind,leaf_class,soil_rough           &
                      ,sfcwater_nlev,stom_condct,ground_rsat,ground_rvap,ground_temp       &
                      ,ground_fliq,veg_water,veg_hcap,veg_energy,can_prss,can_theiv        &
                      ,can_vpdef,can_theta,can_rvap,can_co2,hflxac,wflxac,qwflxac,eflxac   &
                      ,cflxac,hflxgc,wflxgc,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp   &
                      ,intercepted,qintercepted,wshed,qwshed,throughfall,qthroughfall      &
                      ,runoff,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp       &
                      ,veg_ndvip,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t       &
                      ,sflux_r,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)    :: m2
   integer                        , intent(in)    :: m3
   integer                        , intent(in)    :: mzg
   integer                        , intent(in)    :: mzs
   integer                        , intent(in)    :: npat
   integer                        , intent(in)    :: ia
   integer                        , intent(in)    :: iz
   integer                        , intent(in)    :: ja
   integer                        , intent(in)    :: jz
   integer                        , intent(in)    :: jdim
   integer                        , intent(in)    :: ibcon
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_water
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_energy
   real, dimension(    m2,m3,npat), intent(inout) :: soil_color
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_text
   real, dimension(    m2,m3,npat), intent(inout) :: psibar_10d
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_mass
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_energy
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_depth
   real, dimension(    m2,m3,npat), intent(inout) :: ustar
   real, dimension(    m2,m3,npat), intent(inout) :: tstar
   real, dimension(    m2,m3,npat), intent(inout) :: rstar
   real, dimension(    m2,m3,npat), intent(inout) :: cstar
   real, dimension(    m2,m3,npat), intent(inout) :: zeta
   real, dimension(    m2,m3,npat), intent(inout) :: ribulk
   real, dimension(    m2,m3,npat), intent(inout) :: veg_albedo
   real, dimension(    m2,m3,npat), intent(inout) :: veg_fracarea
   real, dimension(    m2,m3,npat), intent(inout) :: veg_lai
   real, dimension(    m2,m3,npat), intent(inout) :: veg_tai
   real, dimension(    m2,m3,npat), intent(inout) :: veg_rough
   real, dimension(    m2,m3,npat), intent(inout) :: veg_height
   real, dimension(    m2,m3,npat), intent(inout) :: veg_displace
   real, dimension(    m2,m3,npat), intent(inout) :: patch_area
   real, dimension(    m2,m3,npat), intent(inout) :: patch_rough
   real, dimension(    m2,m3,npat), intent(inout) :: patch_wetind
   real, dimension(    m2,m3,npat), intent(inout) :: leaf_class
   real, dimension(    m2,m3,npat), intent(inout) :: soil_rough
   real, dimension(    m2,m3,npat), intent(inout) :: sfcwater_nlev
   real, dimension(    m2,m3,npat), intent(inout) :: stom_condct
   real, dimension(    m2,m3,npat), intent(inout) :: ground_rsat
   real, dimension(    m2,m3,npat), intent(inout) :: ground_rvap
   real, dimension(    m2,m3,npat), intent(inout) :: ground_temp
   real, dimension(    m2,m3,npat), intent(inout) :: ground_fliq
   real, dimension(    m2,m3,npat), intent(inout) :: veg_water
   real, dimension(    m2,m3,npat), intent(inout) :: veg_energy
   real, dimension(    m2,m3,npat), intent(inout) :: veg_hcap
   real, dimension(    m2,m3,npat), intent(inout) :: can_prss
   real, dimension(    m2,m3,npat), intent(inout) :: can_theiv
   real, dimension(    m2,m3,npat), intent(inout) :: can_vpdef
   real, dimension(    m2,m3,npat), intent(inout) :: can_theta
   real, dimension(    m2,m3,npat), intent(inout) :: can_rvap
   real, dimension(    m2,m3,npat), intent(inout) :: can_co2
   real, dimension(    m2,m3,npat), intent(inout) :: hflxac
   real, dimension(    m2,m3,npat), intent(inout) :: wflxac
   real, dimension(    m2,m3,npat), intent(inout) :: qwflxac
   real, dimension(    m2,m3,npat), intent(inout) :: eflxac
   real, dimension(    m2,m3,npat), intent(inout) :: cflxac
   real, dimension(    m2,m3,npat), intent(inout) :: hflxgc
   real, dimension(    m2,m3,npat), intent(inout) :: wflxgc
   real, dimension(    m2,m3,npat), intent(inout) :: qwflxgc
   real, dimension(    m2,m3,npat), intent(inout) :: hflxvc
   real, dimension(    m2,m3,npat), intent(inout) :: wflxvc
   real, dimension(    m2,m3,npat), intent(inout) :: qwflxvc
   real, dimension(    m2,m3,npat), intent(inout) :: transp
   real, dimension(    m2,m3,npat), intent(inout) :: qtransp
   real, dimension(    m2,m3,npat), intent(inout) :: intercepted
   real, dimension(    m2,m3,npat), intent(inout) :: qintercepted
   real, dimension(    m2,m3,npat), intent(inout) :: wshed
   real, dimension(    m2,m3,npat), intent(inout) :: qwshed
   real, dimension(    m2,m3,npat), intent(inout) :: throughfall
   real, dimension(    m2,m3,npat), intent(inout) :: qthroughfall
   real, dimension(    m2,m3,npat), intent(inout) :: runoff
   real, dimension(    m2,m3,npat), intent(inout) :: qrunoff
   real, dimension(    m2,m3,npat), intent(inout) :: drainage
   real, dimension(    m2,m3,npat), intent(inout) :: qdrainage
   real, dimension(    m2,m3,npat), intent(inout) :: gpp
   real, dimension(    m2,m3,npat), intent(inout) :: plresp
   real, dimension(    m2,m3,npat), intent(inout) :: resphet
   real, dimension(    m2,m3,npat), intent(inout) :: growresp
   real, dimension(    m2,m3,npat), intent(inout) :: veg_ndvip
   real, dimension(    m2,m3,npat), intent(inout) :: veg_ndvic
   real, dimension(    m2,m3,npat), intent(inout) :: veg_ndvif
   real, dimension(    m2,m3     ), intent(inout) :: sflux_u
   real, dimension(    m2,m3     ), intent(inout) :: sflux_v
   real, dimension(    m2,m3     ), intent(inout) :: sflux_w
   real, dimension(    m2,m3     ), intent(inout) :: sflux_t
   real, dimension(    m2,m3     ), intent(inout) :: sflux_r
   real, dimension(    m2,m3     ), intent(inout) :: sflux_c
   real, dimension(    m2,m3     ), intent(inout) :: albedt
   real, dimension(    m2,m3     ), intent(inout) :: rlongup
   real, dimension(    m2,m3,npat), intent(inout) :: rshort_gnd
   real, dimension(    m2,m3,npat), intent(inout) :: rlong_gnd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: i
   integer                                        :: j
   !---------------------------------------------------------------------------------------!


   !----- Western Boundary ----------------------------------------------------------------!
   if (iand(ibcon,1) /= 0) then
      do j=ja,jz
         call leaf3_clone(m2,m3,mzg,mzs,npat,2,1,j,j                                       &
                         ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color  &
                         ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta &
                         ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough         &
                         ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind      &
                         ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat      &
                         ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap           &
                         ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap       &
                         ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc        &
                         ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted         &
                         ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff        &
                         ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip &
                         ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r      &
                         ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)
      end do
   end if
   !---------------------------------------------------------------------------------------!





   !----- Eastern Boundary ----------------------------------------------------------------!
   if (iand(ibcon,2) /= 0) then
      do j=ja,jz
         call leaf3_clone(m2,m3,mzg,mzs,npat,m2-1,m2,j,j                                   &
                         ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color  &
                         ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta &
                         ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough         &
                         ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind      &
                         ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat      &
                         ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap           &
                         ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap       &
                         ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc        &
                         ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted         &
                         ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff        &
                         ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip &
                         ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r      &
                         ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)
      end do
   end if
   !---------------------------------------------------------------------------------------!




   !----- Southern Boundary ---------------------------------------------------------------!
   if (jdim == 1 .and. iand(ibcon,4) /= 0) then
      do i = ia,iz
         call leaf3_clone(m2,m3,mzg,mzs,npat,i,i,2,1                                       &
                         ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color  &
                         ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta &
                         ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough         &
                         ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind      &
                         ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat      &
                         ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap           &
                         ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap       &
                         ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc        &
                         ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted         &
                         ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff        &
                         ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip &
                         ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r      &
                         ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)
      end do
   end if
   !---------------------------------------------------------------------------------------!




   !----- Northern Boundary ---------------------------------------------------------------!
   if (jdim == 1 .and. iand(ibcon,8) /= 0) then
      do i = ia,iz
         call leaf3_clone(m2,m3,mzg,mzs,npat,i,i,m3-1,m3                                   &
                         ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color  &
                         ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta &
                         ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough         &
                         ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind      &
                         ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat      &
                         ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap           &
                         ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap       &
                         ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc        &
                         ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted         &
                         ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff        &
                         ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip &
                         ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r      &
                         ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)
      end do
   end if
   !---------------------------------------------------------------------------------------!




   !----- Southwestern corner -------------------------------------------------------------!
   if (iand(ibcon,5) /= 0 .or. (iand(ibcon,1) /= 0 .and. jdim == 0)) then
      call leaf3_clone(m2,m3,mzg,mzs,npat,2,1,1+jdim,1                                     &
                      ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color     &
                      ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta    &
                      ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough            &
                      ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind         &
                      ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat         &
                      ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap              &
                      ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap          &
                      ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc           &
                      ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted            &
                      ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff           &
                      ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip    &
                      ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r         &
                      ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)                       
   end if
   !---------------------------------------------------------------------------------------!




   !----- Southeastern corner -------------------------------------------------------------!
   if (iand(ibcon,6) /= 0 .or. (iand(ibcon,2) /= 0 .and. jdim == 0)) then
      call leaf3_clone(m2,m3,mzg,mzs,npat,m2-1,m2,1+jdim,1                                 &
                      ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color     &
                      ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta    &
                      ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough            &
                      ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind         &
                      ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat         &
                      ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap              &
                      ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap          &
                      ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc           &
                      ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted            &
                      ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff           &
                      ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip    &
                      ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r         &
                      ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)                       
   end if
   !---------------------------------------------------------------------------------------!




   !----- Northwestern corner -------------------------------------------------------------!
   if (iand(ibcon,9) /= 0 .or. (iand(ibcon,1) /= 0 .and. jdim == 0)) then
      call leaf3_clone(m2,m3,mzg,mzs,npat,2,1,m3-jdim,m3                                   &
                      ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color     &
                      ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta    &
                      ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough            &
                      ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind         &
                      ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat         &
                      ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap              &
                      ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap          &
                      ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc           &
                      ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted            &
                      ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff           &
                      ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip    &
                      ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r         &
                      ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)                       
   end if
   !---------------------------------------------------------------------------------------!




   !----- Northeastern corner -------------------------------------------------------------!
   if (iand(ibcon,10) /= 0 .or. (iand(ibcon,2) /= 0 .and. jdim == 0)) then
      call leaf3_clone(m2,m3,mzg,mzs,npat,m2-1,m2,m3-jdim,m3                               &
                      ,soil_water,sfcwater_mass,soil_energy,sfcwater_energy,soil_color     &
                      ,soil_text,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta    &
                      ,ribulk,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough            &
                      ,veg_height,veg_displace,patch_area,patch_rough,patch_wetind         &
                      ,leaf_class,soil_rough,sfcwater_nlev,stom_condct,ground_rsat         &
                      ,ground_rvap,ground_temp,ground_fliq,veg_water,veg_hcap              &
                      ,veg_energy,can_prss,can_theiv,can_vpdef,can_theta,can_rvap          &
                      ,can_co2,hflxac,wflxac,qwflxac,eflxac,cflxac,hflxgc,wflxgc           &
                      ,qwflxgc,hflxvc,wflxvc,qwflxvc,transp,qtransp,intercepted            &
                      ,qintercepted,wshed,qwshed,throughfall,qthroughfall,runoff           &
                      ,qrunoff,drainage,qdrainage,gpp,plresp,resphet,growresp,veg_ndvip    &
                      ,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r         &
                      ,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)                       
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_bcond
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will actually copy the boundary conditions.                         !
!------------------------------------------------------------------------------------------!
subroutine leaf3_clone(m2,m3,mzg,mzs,npat,isrc,idest,jsrc,jdest,soil_water                 &
                      ,sfcwater_mass,soil_energy,sfcwater_energy,soil_color,soil_text      &
                      ,psibar_10d,sfcwater_depth,ustar,tstar,rstar,cstar,zeta,ribulk       &
                      ,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough,veg_height        &
                      ,veg_displace,patch_area,patch_rough,patch_wetind,leaf_class         &
                      ,soil_rough,sfcwater_nlev,stom_condct,ground_rsat,ground_rvap        &
                      ,ground_temp,ground_fliq,veg_water,veg_hcap,veg_energy,can_prss      &
                      ,can_theiv,can_vpdef,can_theta,can_rvap,can_co2,hflxac,wflxac        &
                      ,qwflxac,eflxac,cflxac,hflxgc,wflxgc,qwflxgc,hflxvc,wflxvc,qwflxvc   &
                      ,transp,qtransp,intercepted,qintercepted,wshed,qwshed,throughfall    &
                      ,qthroughfall,runoff,qrunoff,drainage,qdrainage,gpp,plresp,resphet   &
                      ,growresp,veg_ndvip,veg_ndvic,veg_ndvif,sflux_u,sflux_v,sflux_w      &
                      ,sflux_t,sflux_r,sflux_c,albedt,rlongup,rshort_gnd,rlong_gnd)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)    :: m2
   integer                        , intent(in)    :: m3
   integer                        , intent(in)    :: mzg
   integer                        , intent(in)    :: mzs
   integer                        , intent(in)    :: npat
   integer                        , intent(in)    :: isrc
   integer                        , intent(in)    :: idest
   integer                        , intent(in)    :: jsrc
   integer                        , intent(in)    :: jdest
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_water
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_energy
   real, dimension(    m2,m3,npat), intent(inout) :: soil_color
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_text
   real, dimension(    m2,m3,npat), intent(inout) :: psibar_10d
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_mass
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_energy
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_depth
   real, dimension(    m2,m3,npat), intent(inout) :: ustar
   real, dimension(    m2,m3,npat), intent(inout) :: tstar
   real, dimension(    m2,m3,npat), intent(inout) :: rstar
   real, dimension(    m2,m3,npat), intent(inout) :: cstar
   real, dimension(    m2,m3,npat), intent(inout) :: zeta
   real, dimension(    m2,m3,npat), intent(inout) :: ribulk
   real, dimension(    m2,m3,npat), intent(inout) :: veg_albedo
   real, dimension(    m2,m3,npat), intent(inout) :: veg_fracarea
   real, dimension(    m2,m3,npat), intent(inout) :: veg_lai
   real, dimension(    m2,m3,npat), intent(inout) :: veg_tai
   real, dimension(    m2,m3,npat), intent(inout) :: veg_rough
   real, dimension(    m2,m3,npat), intent(inout) :: veg_height
   real, dimension(    m2,m3,npat), intent(inout) :: veg_displace
   real, dimension(    m2,m3,npat), intent(inout) :: patch_area
   real, dimension(    m2,m3,npat), intent(inout) :: patch_rough
   real, dimension(    m2,m3,npat), intent(inout) :: patch_wetind
   real, dimension(    m2,m3,npat), intent(inout) :: leaf_class
   real, dimension(    m2,m3,npat), intent(inout) :: soil_rough
   real, dimension(    m2,m3,npat), intent(inout) :: sfcwater_nlev
   real, dimension(    m2,m3,npat), intent(inout) :: stom_condct
   real, dimension(    m2,m3,npat), intent(inout) :: ground_rsat
   real, dimension(    m2,m3,npat), intent(inout) :: ground_rvap
   real, dimension(    m2,m3,npat), intent(inout) :: ground_temp
   real, dimension(    m2,m3,npat), intent(inout) :: ground_fliq
   real, dimension(    m2,m3,npat), intent(inout) :: veg_water
   real, dimension(    m2,m3,npat), intent(inout) :: veg_energy
   real, dimension(    m2,m3,npat), intent(inout) :: veg_hcap
   real, dimension(    m2,m3,npat), intent(inout) :: can_prss
   real, dimension(    m2,m3,npat), intent(inout) :: can_theiv
   real, dimension(    m2,m3,npat), intent(inout) :: can_vpdef
   real, dimension(    m2,m3,npat), intent(inout) :: can_theta
   real, dimension(    m2,m3,npat), intent(inout) :: can_rvap
   real, dimension(    m2,m3,npat), intent(inout) :: can_co2
   real, dimension(    m2,m3,npat), intent(inout) :: hflxac
   real, dimension(    m2,m3,npat), intent(inout) :: wflxac
   real, dimension(    m2,m3,npat), intent(inout) :: qwflxac
   real, dimension(    m2,m3,npat), intent(inout) :: eflxac
   real, dimension(    m2,m3,npat), intent(inout) :: cflxac
   real, dimension(    m2,m3,npat), intent(inout) :: hflxgc
   real, dimension(    m2,m3,npat), intent(inout) :: wflxgc
   real, dimension(    m2,m3,npat), intent(inout) :: qwflxgc
   real, dimension(    m2,m3,npat), intent(inout) :: hflxvc
   real, dimension(    m2,m3,npat), intent(inout) :: wflxvc
   real, dimension(    m2,m3,npat), intent(inout) :: qwflxvc
   real, dimension(    m2,m3,npat), intent(inout) :: transp
   real, dimension(    m2,m3,npat), intent(inout) :: qtransp
   real, dimension(    m2,m3,npat), intent(inout) :: intercepted
   real, dimension(    m2,m3,npat), intent(inout) :: qintercepted
   real, dimension(    m2,m3,npat), intent(inout) :: wshed
   real, dimension(    m2,m3,npat), intent(inout) :: qwshed
   real, dimension(    m2,m3,npat), intent(inout) :: throughfall
   real, dimension(    m2,m3,npat), intent(inout) :: qthroughfall
   real, dimension(    m2,m3,npat), intent(inout) :: runoff
   real, dimension(    m2,m3,npat), intent(inout) :: qrunoff
   real, dimension(    m2,m3,npat), intent(inout) :: drainage
   real, dimension(    m2,m3,npat), intent(inout) :: qdrainage
   real, dimension(    m2,m3,npat), intent(inout) :: gpp
   real, dimension(    m2,m3,npat), intent(inout) :: plresp
   real, dimension(    m2,m3,npat), intent(inout) :: resphet
   real, dimension(    m2,m3,npat), intent(inout) :: growresp
   real, dimension(    m2,m3,npat), intent(inout) :: veg_ndvip
   real, dimension(    m2,m3,npat), intent(inout) :: veg_ndvic
   real, dimension(    m2,m3,npat), intent(inout) :: veg_ndvif
   real, dimension(    m2,m3     ), intent(inout) :: sflux_u
   real, dimension(    m2,m3     ), intent(inout) :: sflux_v
   real, dimension(    m2,m3     ), intent(inout) :: sflux_w
   real, dimension(    m2,m3     ), intent(inout) :: sflux_t
   real, dimension(    m2,m3     ), intent(inout) :: sflux_r
   real, dimension(    m2,m3     ), intent(inout) :: sflux_c
   real, dimension(    m2,m3     ), intent(inout) :: albedt
   real, dimension(    m2,m3     ), intent(inout) :: rlongup
   real, dimension(    m2,m3,npat), intent(inout) :: rshort_gnd
   real, dimension(    m2,m3,npat), intent(inout) :: rlong_gnd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: kzg
   integer                                        :: kzs
   integer                                        :: ipat
   !---------------------------------------------------------------------------------------!


   !----- First the variables that are not patch-dependent. -------------------------------!
   sflux_u(idest,jdest) = sflux_u(isrc,jsrc)
   sflux_v(idest,jdest) = sflux_v(isrc,jsrc)
   sflux_w(idest,jdest) = sflux_w(isrc,jsrc)
   sflux_t(idest,jdest) = sflux_t(isrc,jsrc)
   sflux_r(idest,jdest) = sflux_r(isrc,jsrc)
   sflux_c(idest,jdest) = sflux_c(isrc,jsrc)
   albedt (idest,jdest) = albedt (isrc,jsrc)
   rlongup(idest,jdest) = rlongup(isrc,jsrc)

   !----- Then the patch-dependent variables. ---------------------------------------------!
   patchloop: do ipat = 1,npat

      ustar          (idest,jdest,ipat) = ustar            (isrc,jsrc,ipat)
      tstar          (idest,jdest,ipat) = tstar            (isrc,jsrc,ipat)
      rstar          (idest,jdest,ipat) = rstar            (isrc,jsrc,ipat)
      cstar          (idest,jdest,ipat) = cstar            (isrc,jsrc,ipat)
      zeta           (idest,jdest,ipat) = zeta             (isrc,jsrc,ipat)
      ribulk         (idest,jdest,ipat) = ribulk           (isrc,jsrc,ipat)
      veg_fracarea   (idest,jdest,ipat) = veg_fracarea     (isrc,jsrc,ipat)
      veg_lai        (idest,jdest,ipat) = veg_lai          (isrc,jsrc,ipat)
      veg_tai        (idest,jdest,ipat) = veg_tai          (isrc,jsrc,ipat)
      veg_rough      (idest,jdest,ipat) = veg_rough        (isrc,jsrc,ipat)
      veg_height     (idest,jdest,ipat) = veg_height       (isrc,jsrc,ipat)
      veg_displace   (idest,jdest,ipat) = veg_displace     (isrc,jsrc,ipat)
      patch_area     (idest,jdest,ipat) = patch_area       (isrc,jsrc,ipat)
      patch_rough    (idest,jdest,ipat) = patch_rough      (isrc,jsrc,ipat)
      patch_wetind   (idest,jdest,ipat) = patch_wetind     (isrc,jsrc,ipat)
      leaf_class     (idest,jdest,ipat) = leaf_class       (isrc,jsrc,ipat)
      soil_rough     (idest,jdest,ipat) = soil_rough       (isrc,jsrc,ipat)
      sfcwater_nlev  (idest,jdest,ipat) = sfcwater_nlev    (isrc,jsrc,ipat)
      stom_condct    (idest,jdest,ipat) = stom_condct      (isrc,jsrc,ipat)
      ground_rsat    (idest,jdest,ipat) = ground_rsat      (isrc,jsrc,ipat)
      ground_rvap    (idest,jdest,ipat) = ground_rvap      (isrc,jsrc,ipat)
      ground_temp    (idest,jdest,ipat) = ground_temp      (isrc,jsrc,ipat)
      ground_fliq    (idest,jdest,ipat) = ground_fliq      (isrc,jsrc,ipat)
      veg_water      (idest,jdest,ipat) = veg_water        (isrc,jsrc,ipat)
      veg_hcap       (idest,jdest,ipat) = veg_hcap         (isrc,jsrc,ipat)
      veg_energy     (idest,jdest,ipat) = veg_energy       (isrc,jsrc,ipat)
      can_prss       (idest,jdest,ipat) = can_prss         (isrc,jsrc,ipat)
      can_theiv      (idest,jdest,ipat) = can_theiv        (isrc,jsrc,ipat)
      can_vpdef      (idest,jdest,ipat) = can_vpdef        (isrc,jsrc,ipat)
      can_theta      (idest,jdest,ipat) = can_theta        (isrc,jsrc,ipat)
      can_rvap       (idest,jdest,ipat) = can_rvap         (isrc,jsrc,ipat)
      can_co2        (idest,jdest,ipat) = can_co2          (isrc,jsrc,ipat)
      hflxac         (idest,jdest,ipat) = hflxac           (isrc,jsrc,ipat)
      wflxac         (idest,jdest,ipat) = wflxac           (isrc,jsrc,ipat)
      qwflxac        (idest,jdest,ipat) = qwflxac          (isrc,jsrc,ipat)
      eflxac         (idest,jdest,ipat) = eflxac           (isrc,jsrc,ipat)
      cflxac         (idest,jdest,ipat) = cflxac           (isrc,jsrc,ipat)
      hflxgc         (idest,jdest,ipat) = hflxgc           (isrc,jsrc,ipat)
      wflxgc         (idest,jdest,ipat) = wflxgc           (isrc,jsrc,ipat)
      qwflxgc        (idest,jdest,ipat) = qwflxgc          (isrc,jsrc,ipat)
      hflxvc         (idest,jdest,ipat) = hflxvc           (isrc,jsrc,ipat)
      wflxvc         (idest,jdest,ipat) = wflxvc           (isrc,jsrc,ipat)
      qwflxvc        (idest,jdest,ipat) = qwflxvc          (isrc,jsrc,ipat)
      transp         (idest,jdest,ipat) = transp           (isrc,jsrc,ipat)
      qtransp        (idest,jdest,ipat) = qtransp          (isrc,jsrc,ipat)
      intercepted    (idest,jdest,ipat) = intercepted      (isrc,jsrc,ipat)
      qintercepted   (idest,jdest,ipat) = qintercepted     (isrc,jsrc,ipat)
      wshed          (idest,jdest,ipat) = wshed            (isrc,jsrc,ipat)
      qwshed         (idest,jdest,ipat) = qwshed           (isrc,jsrc,ipat)
      throughfall    (idest,jdest,ipat) = throughfall      (isrc,jsrc,ipat)
      qthroughfall   (idest,jdest,ipat) = qthroughfall     (isrc,jsrc,ipat)
      runoff         (idest,jdest,ipat) = runoff           (isrc,jsrc,ipat)
      qrunoff        (idest,jdest,ipat) = qrunoff          (isrc,jsrc,ipat)
      drainage       (idest,jdest,ipat) = drainage         (isrc,jsrc,ipat)
      qdrainage      (idest,jdest,ipat) = qdrainage        (isrc,jsrc,ipat)
      gpp            (idest,jdest,ipat) = gpp              (isrc,jsrc,ipat)
      plresp         (idest,jdest,ipat) = plresp           (isrc,jsrc,ipat)
      resphet        (idest,jdest,ipat) = resphet          (isrc,jsrc,ipat)
      growresp       (idest,jdest,ipat) = growresp         (isrc,jsrc,ipat)
      veg_ndvip      (idest,jdest,ipat) = veg_ndvip        (isrc,jsrc,ipat)
      veg_ndvic      (idest,jdest,ipat) = veg_ndvic        (isrc,jsrc,ipat)
      veg_ndvif      (idest,jdest,ipat) = veg_ndvif        (isrc,jsrc,ipat)

      rshort_gnd     (idest,jdest,ipat) = rshort_gnd       (isrc,jsrc,ipat)
      rlong_gnd      (idest,jdest,ipat) = rlong_gnd        (isrc,jsrc,ipat)

      soil_color     (idest,jdest,ipat) = soil_color       (isrc,jsrc,ipat)
      psibar_10d     (idest,jdest,ipat) = psibar_10d       (isrc,jsrc,ipat)
      soilloop: do kzg = 1,mzg
         soil_water       (kzg,idest,jdest,ipat) = soil_water         (kzg,isrc,jsrc,ipat)
         soil_energy      (kzg,idest,jdest,ipat) = soil_energy        (kzg,isrc,jsrc,ipat)
         soil_text        (kzg,idest,jdest,ipat) = soil_text          (kzg,isrc,jsrc,ipat)
      end do soilloop

      sfcwloop: do kzs = 1,mzs
         sfcwater_mass    (kzs,idest,jdest,ipat) = sfcwater_mass      (kzs,isrc,jsrc,ipat)
         sfcwater_energy  (kzs,idest,jdest,ipat) = sfcwater_energy    (kzs,isrc,jsrc,ipat)
         sfcwater_depth   (kzs,idest,jdest,ipat) = sfcwater_depth     (kzs,isrc,jsrc,ipat)
      end do sfcwloop
   end do patchloop

   return
end subroutine leaf3_clone
!==========================================================================================!
!==========================================================================================!
