!==========================================================================================!
!==========================================================================================!
subroutine patch_interp_driver(icm,ifm)
   use mem_leaf
   use mem_basic
   use mem_scratch
   use mem_grid
!  use io_params
   implicit none
   integer, intent(in) :: icm, ifm

   call patch_interp(icm,ifm,nzg,nnxp(icm),nnyp(icm),npatch,nzg,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%soil_water,leaf_g(ifm)%soil_water,leaf_g(icm)%patch_area  &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db)
   call patch_interp(icm,ifm,nzg,nnxp(icm),nnyp(icm),npatch,nzg,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%soil_energy,leaf_g(ifm)%soil_energy                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area,scratch%vt3da           &
                    ,scratch%vt3db,scratch%vt2da,scratch%vt2db)
   call patch_interp(icm,ifm,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%sfcwater_mass,leaf_g(ifm)%sfcwater_mass                   &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area,scratch%vt3da           &
                    ,scratch%vt3db,scratch%vt2da,scratch%vt2db)
   call patch_interp(icm,ifm,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%sfcwater_energy,leaf_g(ifm)%sfcwater_energy               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%sfcwater_depth,leaf_g(ifm)%sfcwater_depth                 &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )

   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_fracarea,leaf_g(ifm)%veg_fracarea                     &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_agb,leaf_g(ifm)%veg_agb                               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_lai,leaf_g(ifm)%veg_lai                               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_tai,leaf_g(ifm)%veg_tai                               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_rough,leaf_g(ifm)%veg_rough                           &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_height,leaf_g(ifm)%veg_height                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_albedo,leaf_g(ifm)%veg_albedo                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%patch_rough,leaf_g(ifm)%patch_rough                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_fracarea,leaf_g(ifm)%veg_fracarea                     &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%patch_wetind,leaf_g(ifm)%patch_wetind                     &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%soil_rough,leaf_g(ifm)%soil_rough                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db)
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%sfcwater_nlev,leaf_g(ifm)%sfcwater_nlev                   &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%stom_resist,leaf_g(ifm)%stom_resist                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_rsat,leaf_g(ifm)%ground_rsat                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_rvap,leaf_g(ifm)%ground_rvap                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_temp,leaf_g(ifm)%ground_temp                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_fliq,leaf_g(ifm)%ground_fliq                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_water,leaf_g(ifm)%veg_water                           &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_hcap,leaf_g(ifm)%veg_hcap                             &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_energy,leaf_g(ifm)%veg_energy                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_prss,leaf_g(ifm)%can_prss,leaf_g(icm)%patch_area      &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_rvap,leaf_g(ifm)%can_rvap,leaf_g(icm)%patch_area      &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_co2,leaf_g(ifm)%can_co2,leaf_g(icm)%patch_area        &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_theiv,leaf_g(ifm)%can_theiv,leaf_g(icm)%patch_area    &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_theta,leaf_g(ifm)%can_theta,leaf_g(icm)%patch_area    &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_ndvic,leaf_g(ifm)%veg_ndvic,leaf_g(icm)%patch_area    &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%zeta,leaf_g(ifm)%zeta,leaf_g(icm)%patch_area              &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ribulk,leaf_g(ifm)%ribulk,leaf_g(icm)%patch_area          &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )

   return
end subroutine patch_interp_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This was broken into various subroutines to avoid compiler crashes...                !
!------------------------------------------------------------------------------------------!
subroutine coarse2fine_driver(icm,ifm)
   use mem_grid, only : nnxp   & ! intent(in)
                      , nnyp   & ! intent(in)
                      , nzg    & ! intent(in)
                      , nzs    & ! intent(in)
                      , npatch ! ! intent(in)
   use mem_leaf, only : leaf_g ! ! intent(inout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: icm
   integer, intent(in) :: ifm
   !---------------------------------------------------------------------------------------!

   call coarse2fine(ifm,nnxp(ifm),nnyp(ifm),icm,nnxp(icm),nnyp(icm),nzg,nzs,npatch         &
                                 ,leaf_g(ifm)%soil_water     , leaf_g(icm)%soil_water      &
                                 ,leaf_g(ifm)%soil_energy    , leaf_g(icm)%soil_energy     &
                                 ,leaf_g(ifm)%sfcwater_mass  , leaf_g(icm)%sfcwater_mass   &
                                 ,leaf_g(ifm)%sfcwater_energy, leaf_g(icm)%sfcwater_energy &
                                 ,leaf_g(ifm)%sfcwater_depth , leaf_g(icm)%sfcwater_depth  &
                                 ,leaf_g(ifm)%veg_fracarea   , leaf_g(icm)%veg_fracarea    &
                                 ,leaf_g(ifm)%veg_agb        , leaf_g(icm)%veg_agb         &
                                 ,leaf_g(ifm)%veg_lai        , leaf_g(icm)%veg_lai         &
                                 ,leaf_g(ifm)%veg_tai        , leaf_g(icm)%veg_tai         &
                                 ,leaf_g(ifm)%veg_rough      , leaf_g(icm)%veg_rough       &
                                 ,leaf_g(ifm)%veg_height     , leaf_g(icm)%veg_height      &
                                 ,leaf_g(ifm)%veg_albedo     , leaf_g(icm)%veg_albedo      &
                                 ,leaf_g(ifm)%patch_rough    , leaf_g(icm)%patch_rough     &
                                 ,leaf_g(ifm)%patch_wetind   , leaf_g(icm)%patch_wetind    &
                                 ,leaf_g(ifm)%soil_rough     , leaf_g(icm)%soil_rough      &
                                 ,leaf_g(ifm)%sfcwater_nlev  , leaf_g(icm)%sfcwater_nlev   &
                                 ,leaf_g(ifm)%stom_resist    , leaf_g(icm)%stom_resist     &
                                 ,leaf_g(ifm)%ground_rsat    , leaf_g(icm)%ground_rsat     &
                                 ,leaf_g(ifm)%ground_rvap    , leaf_g(icm)%ground_rvap     &
                                 ,leaf_g(ifm)%ground_temp    , leaf_g(icm)%ground_temp     &
                                 ,leaf_g(ifm)%ground_fliq    , leaf_g(icm)%ground_fliq     &
                                 ,leaf_g(ifm)%veg_water      , leaf_g(icm)%veg_water       &
                                 ,leaf_g(ifm)%veg_energy     , leaf_g(icm)%veg_energy      &
                                 ,leaf_g(ifm)%veg_hcap       , leaf_g(icm)%veg_hcap        &
                                 ,leaf_g(ifm)%can_rvap       , leaf_g(icm)%can_rvap        &
                                 ,leaf_g(ifm)%can_co2        , leaf_g(icm)%can_co2         &
                                 ,leaf_g(ifm)%can_theiv      , leaf_g(icm)%can_theiv       &
                                 ,leaf_g(ifm)%can_theta      , leaf_g(icm)%can_theta       &
                                 ,leaf_g(ifm)%can_prss       , leaf_g(icm)%can_prss        &
                                 ,leaf_g(ifm)%veg_ndvic      , leaf_g(icm)%veg_ndvic       &
                                 ,leaf_g(ifm)%sensible_gc    , leaf_g(icm)%sensible_gc     &
                                 ,leaf_g(ifm)%sensible_vc    , leaf_g(icm)%sensible_vc     &
                                 ,leaf_g(ifm)%evap_gc        , leaf_g(icm)%evap_gc         &
                                 ,leaf_g(ifm)%evap_vc        , leaf_g(icm)%evap_vc         &
                                 ,leaf_g(ifm)%transp         , leaf_g(icm)%transp          &
                                 ,leaf_g(ifm)%gpp            , leaf_g(icm)%gpp             &
                                 ,leaf_g(ifm)%plresp         , leaf_g(icm)%plresp          &
                                 ,leaf_g(ifm)%resphet        , leaf_g(icm)%resphet         )
   return
end subroutine coarse2fine_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine coarse2fine(ifm,mxpf,mypf,icm,mxpc,mypc,mzg,mzs,mpat                            &
                ,f_soil_water     , c_soil_water     ,f_soil_energy    , c_soil_energy     &
                ,f_sfcwater_mass  , c_sfcwater_mass  ,f_sfcwater_energy, c_sfcwater_energy &
                ,f_sfcwater_depth , c_sfcwater_depth ,f_veg_fracarea   , c_veg_fracarea    &
                ,f_veg_agb        , c_veg_agb        ,f_veg_lai        , c_veg_lai         &
                ,f_veg_tai        , c_veg_tai        ,f_veg_rough      , c_veg_rough       &
                ,f_veg_height     , c_veg_height     ,f_veg_albedo     , c_veg_albedo      &
                ,f_patch_rough    , c_patch_rough    ,f_patch_wetind   , c_patch_wetind    &
                ,f_soil_rough     , c_soil_rough     ,f_sfcwater_nlev  , c_sfcwater_nlev   &
                ,f_stom_resist    , c_stom_resist    ,f_ground_rsat    , c_ground_rsat     &
                ,f_ground_rvap    , c_ground_rvap    ,f_ground_temp    , c_ground_temp     &
                ,f_ground_fliq    , c_ground_fliq    ,f_veg_water      , c_veg_water       &
                ,f_veg_energy     , c_veg_energy     ,f_veg_hcap       , c_veg_hcap        &
                ,f_can_rvap       , c_can_rvap       ,f_can_co2        , c_can_co2         &
                ,f_can_theiv      , c_can_theiv      ,f_can_theta      , c_can_theta       &
                ,f_can_prss       , c_can_prss       ,f_veg_ndvic      , c_veg_ndvic       &
                ,f_sensible_gc    , c_sensible_gc    ,f_sensible_vc    , c_sensible_vc     &
                ,f_evap_gc        , c_evap_gc        ,f_evap_vc        , c_evap_vc         &
                ,f_transp         , c_transp         ,f_gpp            , c_gpp             &
                ,f_plresp         , c_plresp         ,f_resphet        , c_resphet         )
   use mem_grid, only : ipm & ! intent(in)
                      , jpm ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                            , intent(in)  :: ifm
   integer                            , intent(in)  :: mxpf
   integer                            , intent(in)  :: mypf
   integer                            , intent(in)  :: icm
   integer                            , intent(in)  :: mxpc
   integer                            , intent(in)  :: mypc
   integer                            , intent(in)  :: mzg
   integer                            , intent(in)  :: mzs
   integer                            , intent(in)  :: mpat
   real, dimension(mzg,mxpc,mypc,mpat), intent(in)  :: c_soil_water
   real, dimension(mzg,mxpc,mypc,mpat), intent(in)  :: c_soil_energy
   real, dimension(mzs,mxpc,mypc,mpat), intent(in)  :: c_sfcwater_mass
   real, dimension(mzs,mxpc,mypc,mpat), intent(in)  :: c_sfcwater_energy
   real, dimension(mzs,mxpc,mypc,mpat), intent(in)  :: c_sfcwater_depth
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_fracarea
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_agb
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_lai
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_tai
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_rough
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_height
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_albedo
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_patch_rough
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_patch_wetind
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_soil_rough
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_sfcwater_nlev
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_stom_resist
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_ground_rsat
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_ground_rvap
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_ground_temp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_ground_fliq
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_water
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_energy
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_hcap
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_rvap
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_co2
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_theiv
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_theta
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_prss
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_ndvic
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_sensible_gc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_sensible_vc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_evap_gc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_evap_vc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_transp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_gpp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_plresp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_resphet
   real, dimension(mzg,mxpf,mypf,mpat), intent(out) :: f_soil_water
   real, dimension(mzg,mxpf,mypf,mpat), intent(out) :: f_soil_energy
   real, dimension(mzs,mxpf,mypf,mpat), intent(out) :: f_sfcwater_mass
   real, dimension(mzs,mxpf,mypf,mpat), intent(out) :: f_sfcwater_energy
   real, dimension(mzs,mxpf,mypf,mpat), intent(out) :: f_sfcwater_depth
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_fracarea
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_agb
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_lai
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_tai
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_rough
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_height
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_albedo
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_patch_rough
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_patch_wetind
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_soil_rough
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_sfcwater_nlev
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_stom_resist
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_ground_rsat
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_ground_rvap
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_ground_temp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_ground_fliq
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_water
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_energy
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_hcap
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_rvap
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_co2
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_theiv
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_theta
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_prss
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_ndvic
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_sensible_gc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_sensible_vc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_evap_gc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_evap_vc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_transp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_gpp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_plresp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_resphet
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: ic
   integer                                          :: jc
   integer                                          :: k
   integer                                          :: ipat
   !---------------------------------------------------------------------------------------!


   patchloop: do ipat = 1,mpat
      latloop: do j = 1,mypf
         lonloop: do i = 1,mxpf
            ic = ipm(i,ifm)
            jc = jpm(j,ifm) 

            gndloop: do k = 1,mzg
               f_soil_water      (k,i,j,ipat) = c_soil_water      (k,ic,jc,ipat)
               f_soil_energy     (k,i,j,ipat) = c_soil_energy     (k,ic,jc,ipat)
            end do gndloop

            snowloop:do k = 1,mzs
               f_sfcwater_mass   (k,i,j,ipat) = c_sfcwater_mass   (k,ic,jc,ipat)
               f_sfcwater_energy (k,i,j,ipat) = c_sfcwater_energy (k,ic,jc,ipat)  
               f_sfcwater_depth  (k,i,j,ipat) = c_sfcwater_depth  (k,ic,jc,ipat)
            end do snowloop

            f_veg_fracarea         (i,j,ipat) = c_veg_fracarea    (ic,jc,ipat)
            f_veg_agb              (i,j,ipat) = c_veg_agb         (ic,jc,ipat)
            f_veg_lai              (i,j,ipat) = c_veg_lai         (ic,jc,ipat)
            f_veg_tai              (i,j,ipat) = c_veg_tai         (ic,jc,ipat)
            f_veg_rough            (i,j,ipat) = c_veg_rough       (ic,jc,ipat)
            f_veg_height           (i,j,ipat) = c_veg_height      (ic,jc,ipat)
            f_veg_albedo           (i,j,ipat) = c_veg_albedo      (ic,jc,ipat)
            f_patch_rough          (i,j,ipat) = c_patch_rough     (ic,jc,ipat)
            f_patch_wetind         (i,j,ipat) = c_patch_wetind    (ic,jc,ipat)
            f_soil_rough           (i,j,ipat) = c_soil_rough      (ic,jc,ipat)
            f_sfcwater_nlev        (i,j,ipat) = c_sfcwater_nlev   (ic,jc,ipat)
            f_stom_resist          (i,j,ipat) = c_stom_resist     (ic,jc,ipat) 
            f_ground_rsat          (i,j,ipat) = c_ground_rsat     (ic,jc,ipat)
            f_ground_rvap          (i,j,ipat) = c_ground_rvap     (ic,jc,ipat)
            f_ground_temp          (i,j,ipat) = c_ground_temp     (ic,jc,ipat)
            f_ground_fliq          (i,j,ipat) = c_ground_fliq     (ic,jc,ipat)
            f_veg_water            (i,j,ipat) = c_veg_water       (ic,jc,ipat)
            f_veg_energy           (i,j,ipat) = c_veg_energy      (ic,jc,ipat) 
            f_veg_hcap             (i,j,ipat) = c_veg_hcap        (ic,jc,ipat) 
            f_can_rvap             (i,j,ipat) = c_can_rvap        (ic,jc,ipat)
            f_can_co2              (i,j,ipat) = c_can_co2         (ic,jc,ipat)
            f_can_theiv            (i,j,ipat) = c_can_theiv       (ic,jc,ipat) 
            f_can_theta            (i,j,ipat) = c_can_theta       (ic,jc,ipat) 
            f_can_prss             (i,j,ipat) = c_can_prss        (ic,jc,ipat) 
            f_veg_ndvic            (i,j,ipat) = c_veg_ndvic       (ic,jc,ipat)
            f_sensible_gc          (i,j,ipat) = c_sensible_gc     (ic,jc,ipat)
            f_sensible_vc          (i,j,ipat) = c_sensible_vc     (ic,jc,ipat)
            f_evap_gc              (i,j,ipat) = c_evap_gc         (ic,jc,ipat)
            f_evap_vc              (i,j,ipat) = c_evap_vc         (ic,jc,ipat)
            f_transp               (i,j,ipat) = c_transp          (ic,jc,ipat)
            f_gpp                  (i,j,ipat) = c_gpp             (ic,jc,ipat)
            f_plresp               (i,j,ipat) = c_plresp          (ic,jc,ipat)
            f_resphet              (i,j,ipat) = c_resphet         (ic,jc,ipat)
         end do lonloop
      end do latloop
   end do patchloop

   return
end subroutine coarse2fine
!==========================================================================================!
!==========================================================================================!
