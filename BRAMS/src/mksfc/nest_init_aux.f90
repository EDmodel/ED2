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
                    ,leaf_g(icm)%veg_displace,leaf_g(ifm)%veg_displace                     &
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
                    ,leaf_g(icm)%stom_condct,leaf_g(ifm)%stom_condct                       &
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
                    ,leaf_g(icm)%can_vpdef,leaf_g(ifm)%can_vpdef,leaf_g(icm)%patch_area    &
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
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%psibar_10d,leaf_g(ifm)%psibar_10d,leaf_g(icm)%patch_area  &
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
                                 ,leaf_g(ifm)%veg_displace   , leaf_g(icm)%veg_displace    &
                                 ,leaf_g(ifm)%veg_albedo     , leaf_g(icm)%veg_albedo      &
                                 ,leaf_g(ifm)%patch_rough    , leaf_g(icm)%patch_rough     &
                                 ,leaf_g(ifm)%patch_wetind   , leaf_g(icm)%patch_wetind    &
                                 ,leaf_g(ifm)%soil_rough     , leaf_g(icm)%soil_rough      &
                                 ,leaf_g(ifm)%sfcwater_nlev  , leaf_g(icm)%sfcwater_nlev   &
                                 ,leaf_g(ifm)%stom_condct    , leaf_g(icm)%stom_condct     &
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
                                 ,leaf_g(ifm)%can_vpdef      , leaf_g(icm)%can_vpdef       &
                                 ,leaf_g(ifm)%can_theta      , leaf_g(icm)%can_theta       &
                                 ,leaf_g(ifm)%can_prss       , leaf_g(icm)%can_prss        &
                                 ,leaf_g(ifm)%veg_ndvic      , leaf_g(icm)%veg_ndvic       &
                                 ,leaf_g(ifm)%hflxac         , leaf_g(icm)%hflxac          &
                                 ,leaf_g(ifm)%wflxac         , leaf_g(icm)%wflxac          &
                                 ,leaf_g(ifm)%qwflxac        , leaf_g(icm)%qwflxac         &
                                 ,leaf_g(ifm)%eflxac         , leaf_g(icm)%eflxac          &
                                 ,leaf_g(ifm)%cflxac         , leaf_g(icm)%cflxac          &
                                 ,leaf_g(ifm)%hflxgc         , leaf_g(icm)%hflxgc          &
                                 ,leaf_g(ifm)%wflxgc         , leaf_g(icm)%wflxgc          &
                                 ,leaf_g(ifm)%qwflxgc        , leaf_g(icm)%qwflxgc         &
                                 ,leaf_g(ifm)%hflxvc         , leaf_g(icm)%hflxvc          &
                                 ,leaf_g(ifm)%wflxvc         , leaf_g(icm)%wflxvc          &
                                 ,leaf_g(ifm)%qwflxvc        , leaf_g(icm)%qwflxvc         &
                                 ,leaf_g(ifm)%transp         , leaf_g(icm)%transp          &
                                 ,leaf_g(ifm)%qtransp        , leaf_g(icm)%qtransp         &
                                 ,leaf_g(ifm)%intercepted    , leaf_g(icm)%intercepted     &
                                 ,leaf_g(ifm)%qintercepted   , leaf_g(icm)%qintercepted    &
                                 ,leaf_g(ifm)%wshed          , leaf_g(icm)%wshed           &
                                 ,leaf_g(ifm)%qwshed         , leaf_g(icm)%qwshed          &
                                 ,leaf_g(ifm)%throughfall    , leaf_g(icm)%throughfall     &
                                 ,leaf_g(ifm)%qthroughfall   , leaf_g(icm)%qthroughfall    &
                                 ,leaf_g(ifm)%runoff         , leaf_g(icm)%runoff          &
                                 ,leaf_g(ifm)%qrunoff        , leaf_g(icm)%qrunoff         &
                                 ,leaf_g(ifm)%drainage       , leaf_g(icm)%drainage        &
                                 ,leaf_g(ifm)%qdrainage      , leaf_g(icm)%qdrainage       &
                                 ,leaf_g(ifm)%psibar_10d     , leaf_g(icm)%psibar_10d      &
                                 ,leaf_g(ifm)%gpp            , leaf_g(icm)%gpp             &
                                 ,leaf_g(ifm)%plresp         , leaf_g(icm)%plresp          &
                                 ,leaf_g(ifm)%resphet        , leaf_g(icm)%resphet         &
                                 ,leaf_g(ifm)%growresp       , leaf_g(icm)%growresp        )
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
                ,f_veg_height     , c_veg_height     ,f_veg_displace   , c_veg_displace    &
                ,f_veg_albedo     , c_veg_albedo     ,f_patch_rough    , c_patch_rough     &
                ,f_patch_wetind   , c_patch_wetind   ,f_soil_rough     , c_soil_rough      &
                ,f_sfcwater_nlev  , c_sfcwater_nlev  ,f_stom_condct    , c_stom_condct     &
                ,f_ground_rsat    , c_ground_rsat    ,f_ground_rvap    , c_ground_rvap     &
                ,f_ground_temp    , c_ground_temp    ,f_ground_fliq    , c_ground_fliq     &
                ,f_veg_water      , c_veg_water      ,f_veg_energy     , c_veg_energy      &
                ,f_veg_hcap       , c_veg_hcap       ,f_can_rvap       , c_can_rvap        &
                ,f_can_co2        , c_can_co2        ,f_can_theiv      , c_can_theiv       &
                ,f_can_vpdef      , c_can_vpdef      ,f_can_theta      , c_can_theta       &
                ,f_can_prss       , c_can_prss       ,f_veg_ndvic      , c_veg_ndvic       &
                ,f_hflxac         , c_hflxac         ,f_wflxac         , c_wflxac          &
                ,f_qwflxac        , c_qwflxac        ,f_eflxac         , c_eflxac          &
                ,f_cflxac         , c_cflxac         ,f_hflxgc         , c_hflxgc          &
                ,f_wflxgc         , c_wflxgc         ,f_qwflxgc        , c_qwflxgc         &
                ,f_hflxvc         , c_hflxvc         ,f_wflxvc         , c_wflxvc          &
                ,f_qwflxvc        , c_qwflxvc        ,f_transp         , c_transp          &
                ,f_qtransp        , c_qtransp        ,f_intercepted    , c_intercepted     &
                ,f_qintercepted   , c_qintercepted   ,f_wshed          , c_wshed           &
                ,f_qwshed         , c_qwshed         ,f_throughfall    , c_throughfall     &
                ,f_qthroughfall   , c_qthroughfall   ,f_runoff         , c_runoff          &
                ,f_qrunoff        , c_qrunoff        ,f_drainage       , c_drainage        &
                ,f_qdrainage      , c_qdrainage      ,f_psibar_10d     , c_psibar_10d      &
                ,f_gpp            , c_gpp            ,f_plresp         , c_plresp          &
                ,f_resphet        , c_resphet        ,f_growresp       , c_growresp        )
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
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_displace
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_albedo
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_patch_rough
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_patch_wetind
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_soil_rough
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_sfcwater_nlev
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_stom_condct
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
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_vpdef
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_theta
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_can_prss
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_veg_ndvic
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_hflxac
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_wflxac
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qwflxac
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_eflxac
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_cflxac
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_hflxgc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_wflxgc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qwflxgc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_hflxvc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_wflxvc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qwflxvc
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_transp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qtransp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_intercepted
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qintercepted
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_wshed
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qwshed
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_throughfall
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qthroughfall
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_runoff
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qrunoff
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_drainage
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_qdrainage
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_psibar_10d
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_gpp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_plresp
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_resphet
   real, dimension(    mxpc,mypc,mpat), intent(in)  :: c_growresp
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
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_displace
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_albedo
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_patch_rough
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_patch_wetind
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_soil_rough
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_sfcwater_nlev
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_stom_condct
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
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_vpdef
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_theta
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_can_prss
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_veg_ndvic
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_hflxac
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_wflxac
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qwflxac
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_eflxac
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_cflxac
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_hflxgc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_wflxgc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qwflxgc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_hflxvc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_wflxvc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qwflxvc
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_transp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qtransp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_intercepted
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qintercepted
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_wshed
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qwshed
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_throughfall
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qthroughfall
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_runoff
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qrunoff
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_drainage
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_qdrainage
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_psibar_10d
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_gpp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_plresp
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_resphet
   real, dimension(    mxpf,mypf,mpat), intent(out) :: f_growresp
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
            f_veg_displace         (i,j,ipat) = c_veg_displace    (ic,jc,ipat)
            f_veg_albedo           (i,j,ipat) = c_veg_albedo      (ic,jc,ipat)
            f_patch_rough          (i,j,ipat) = c_patch_rough     (ic,jc,ipat)
            f_patch_wetind         (i,j,ipat) = c_patch_wetind    (ic,jc,ipat)
            f_soil_rough           (i,j,ipat) = c_soil_rough      (ic,jc,ipat)
            f_sfcwater_nlev        (i,j,ipat) = c_sfcwater_nlev   (ic,jc,ipat)
            f_stom_condct          (i,j,ipat) = c_stom_condct     (ic,jc,ipat) 
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
            f_can_vpdef            (i,j,ipat) = c_can_vpdef       (ic,jc,ipat) 
            f_can_theta            (i,j,ipat) = c_can_theta       (ic,jc,ipat) 
            f_can_prss             (i,j,ipat) = c_can_prss        (ic,jc,ipat) 
            f_veg_ndvic            (i,j,ipat) = c_veg_ndvic       (ic,jc,ipat)
            f_hflxac               (i,j,ipat) = c_hflxac          (ic,jc,ipat)
            f_wflxac               (i,j,ipat) = c_wflxac          (ic,jc,ipat)
            f_qwflxac              (i,j,ipat) = c_qwflxac         (ic,jc,ipat)
            f_eflxac               (i,j,ipat) = c_eflxac          (ic,jc,ipat)
            f_cflxac               (i,j,ipat) = c_cflxac          (ic,jc,ipat)
            f_hflxgc               (i,j,ipat) = c_hflxgc          (ic,jc,ipat)
            f_wflxgc               (i,j,ipat) = c_wflxgc          (ic,jc,ipat)
            f_qwflxgc              (i,j,ipat) = c_qwflxgc         (ic,jc,ipat)
            f_hflxvc               (i,j,ipat) = c_hflxvc          (ic,jc,ipat)
            f_wflxvc               (i,j,ipat) = c_wflxvc          (ic,jc,ipat)
            f_qwflxvc              (i,j,ipat) = c_qwflxvc         (ic,jc,ipat)
            f_transp               (i,j,ipat) = c_transp          (ic,jc,ipat)
            f_qtransp              (i,j,ipat) = c_qtransp         (ic,jc,ipat)
            f_intercepted          (i,j,ipat) = c_intercepted     (ic,jc,ipat)
            f_qintercepted         (i,j,ipat) = c_qintercepted    (ic,jc,ipat)
            f_wshed                (i,j,ipat) = c_wshed           (ic,jc,ipat)
            f_qwshed               (i,j,ipat) = c_qwshed          (ic,jc,ipat)
            f_throughfall          (i,j,ipat) = c_throughfall     (ic,jc,ipat)
            f_qthroughfall         (i,j,ipat) = c_qthroughfall    (ic,jc,ipat)
            f_runoff               (i,j,ipat) = c_runoff          (ic,jc,ipat)
            f_qrunoff              (i,j,ipat) = c_qrunoff         (ic,jc,ipat)
            f_drainage             (i,j,ipat) = c_drainage        (ic,jc,ipat)
            f_qdrainage            (i,j,ipat) = c_qdrainage       (ic,jc,ipat)
            f_psibar_10d           (i,j,ipat) = c_psibar_10d      (ic,jc,ipat)
            f_gpp                  (i,j,ipat) = c_gpp             (ic,jc,ipat)
            f_plresp               (i,j,ipat) = c_plresp          (ic,jc,ipat)
            f_resphet              (i,j,ipat) = c_resphet         (ic,jc,ipat)
            f_growresp             (i,j,ipat) = c_growresp        (ic,jc,ipat)
         end do lonloop
      end do latloop
   end do patchloop

   return
end subroutine coarse2fine
!==========================================================================================!
!==========================================================================================!
