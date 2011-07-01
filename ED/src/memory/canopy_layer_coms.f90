!==========================================================================================!
!==========================================================================================!
!     Module canopy_layer_coms.                                                            !
!                                                                                          !
!     This module contains some variables used by the different models we use to           !
! represent the canopy structure.  Currently the subroutines that use this model are those !
! related to radiation and roughness, and to some extent the rainfall interception.        !
!------------------------------------------------------------------------------------------!
module canopy_layer_coms


   !---------------------------------------------------------------------------------------!
   !     Crown model flag (set at the namelist):                                           !
   !     0. Original;                                                                      !
   !     1. Finite-crown mixing model (Dietze 2008).                                       !
   !     2. Finite-crown mixing model with horizontal competition.                         !
   !---------------------------------------------------------------------------------------!
   integer :: crown_mod
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Variables that define the number of layers in the canopy.                        !
   !   The height of the top of each layer is defined as:                                  !
   !   Ht(n) = Ht(0) * n^Eh                                                                !
   !   where H(0) is the top of the first layer (n = 1)                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4)    :: zztop0     ! Top of the first layer
   real(kind=8)    :: zztop08    ! Top of the first layer
   real(kind=4)    :: zztop0i    ! 1/zztop0
   real(kind=8)    :: zztop0i8   ! 1/zztop08
   real(kind=4)    :: ehgt       ! Eh - Exp. to define the rate of increase of delta Z
   real(kind=8)    :: ehgt8      ! Eh - Exp. to define the rate of increase of delta Z
   real(kind=4)    :: ehgti      ! 1/Eh
   real(kind=8)    :: ehgti8     ! 1/Eh
   integer         :: ncanlyr    ! Number of layers.
   integer         :: ncanlyrp1  ! Number of layers + 1
   integer         :: ncanlyrt2  ! Number of layers * 2
   !---------------------------------------------------------------------------------------!


   !----- Vertical levels (bottom, middle, and top) of each layer. ------------------------!
   real(kind=4), dimension(:)            , allocatable :: zztop    ! Top of each layer
   real(kind=4), dimension(:)            , allocatable :: zzmid    ! Middle of each layer
   real(kind=4), dimension(:)            , allocatable :: zzbot    ! Bottom of each layer
   real(kind=4), dimension(:)            , allocatable :: dzcan    ! Layer thickness
   real(kind=8), dimension(:)            , allocatable :: zztop8   ! Top of each layer
   real(kind=8), dimension(:)            , allocatable :: zzmid8   ! Middle of each layer
   real(kind=8), dimension(:)            , allocatable :: zzbot8   ! Bottom of each layer
   real(kind=8), dimension(:)            , allocatable :: dzcan8   ! Layer thickness
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Fraction of open canopy, used by both roughness and radiation.                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4), dimension(:)   , allocatable :: opencan
   real(kind=8), dimension(:)   , allocatable :: opencan8
   !---------------------------------------------------------------------------------------!


   !----- Variables used by the roughness scheme. -----------------------------------------!
   real(kind=4), dimension(:), allocatable :: lad          ! Leaf area density
   real(kind=4), dimension(:), allocatable :: cdrag        ! Drag coefficient
   real(kind=4), dimension(:), allocatable :: pshelter     ! Sheltering factor
   real(kind=4), dimension(:), allocatable :: cumldrag     ! Cumulative leaf drag area fctn.
   real(kind=4), dimension(:), allocatable :: windlyr      ! Wind profile
   real(kind=4), dimension(:), allocatable :: windext_full ! Full Wind extinction 
   real(kind=4), dimension(:), allocatable :: windext_half ! Half-layer wind extinction
   real(kind=8), dimension(:), allocatable :: lad8
   real(kind=8), dimension(:), allocatable :: cdrag8
   real(kind=8), dimension(:), allocatable :: pshelter8
   real(kind=8), dimension(:), allocatable :: cumldrag8
   real(kind=8), dimension(:), allocatable :: windlyr8
   real(kind=8), dimension(:), allocatable :: windext_full8
   real(kind=8), dimension(:), allocatable :: windext_half8
   !---------------------------------------------------------------------------------------!



   !----- Variables that are used by the radiation scheme. --------------------------------!
   integer     , dimension(:)   , allocatable :: indx
   logical     , dimension(:)   , allocatable :: populated
   real(kind=8), dimension(:,:) , allocatable :: matal
   real(kind=8), dimension(:,:) , allocatable :: mastermat
   real(kind=8), dimension(:,:) , allocatable :: masmatcp
   real(kind=8), dimension(:)   , allocatable :: layer_scatter
   real(kind=8), dimension(:)   , allocatable :: layer_backscatter
   real(kind=8), dimension(:)   , allocatable :: layer_clumping
   real(kind=8), dimension(:)   , allocatable :: expkl_top
   real(kind=8), dimension(:)   , allocatable :: expkl_bot
   real(kind=8), dimension(:)   , allocatable :: expamk_top
   real(kind=8), dimension(:)   , allocatable :: expamk_bot
   real(kind=8), dimension(:)   , allocatable :: expapk_top
   real(kind=8), dimension(:)   , allocatable :: expapk_bot
   real(kind=8), dimension(:)   , allocatable :: A_top
   real(kind=8), dimension(:)   , allocatable :: A_bot
   real(kind=8), dimension(:)   , allocatable :: B_top
   real(kind=8), dimension(:)   , allocatable :: B_bot
   real(kind=8), dimension(:)   , allocatable :: C_top
   real(kind=8), dimension(:)   , allocatable :: C_bot
   real(kind=8), dimension(:)   , allocatable :: F_top
   real(kind=8), dimension(:)   , allocatable :: F_bot
   real(kind=8), dimension(:)   , allocatable :: G_top
   real(kind=8), dimension(:)   , allocatable :: G_bot
   real(kind=8), dimension(:)   , allocatable :: H_top
   real(kind=8), dimension(:)   , allocatable :: H_bot
   real(kind=8), dimension(:)   , allocatable :: beam_bot
   real(kind=8), dimension(:)   , allocatable :: beam_mid
   real(kind=8), dimension(:)   , allocatable :: beam_bot_crown
   real(kind=8), dimension(:)   , allocatable :: upward_vis_beam
   real(kind=8), dimension(:)   , allocatable :: upward_vis_diffuse
   real(kind=8), dimension(:)   , allocatable :: upward_nir_beam
   real(kind=8), dimension(:)   , allocatable :: upward_nir_diffuse
   real(kind=8), dimension(:)   , allocatable :: downward_nir_beam
   real(kind=8), dimension(:)   , allocatable :: downward_nir_diffuse
   real(kind=8), dimension(:)   , allocatable :: downward_vis_beam
   real(kind=8), dimension(:)   , allocatable :: downward_vis_diffuse
   real(kind=8), dimension(:)   , allocatable :: mastervec_beam
   real(kind=8), dimension(:)   , allocatable :: masveccp_beam
   real(kind=8), dimension(:)   , allocatable :: mastervec_diffuse
   real(kind=8), dimension(:)   , allocatable :: masveccp_diffuse
   real(kind=4), dimension(:)   , allocatable :: PAR_beam_layer
   real(kind=4), dimension(:)   , allocatable :: PAR_diffuse_layer
   real(kind=4), dimension(:)   , allocatable :: SW_abs_beam_layer
   real(kind=4), dimension(:)   , allocatable :: SW_abs_diffuse_layer

   real(kind=8), dimension(:)   , allocatable :: dzcanpop
   real(kind=8), dimension(:)   , allocatable :: mastervec_surf
   real(kind=8), dimension(:)   , allocatable :: mastervec_incid
   real(kind=8), dimension(:)   , allocatable :: layer_emis
   real(kind=8), dimension(:)   , allocatable :: layer_temp
   real(kind=8), dimension(:)   , allocatable :: explai
   real(kind=8), dimension(:)   , allocatable :: exmlai
   real(kind=8), dimension(:)   , allocatable :: downward_lw_incid
   real(kind=8), dimension(:)   , allocatable :: downward_lw_surf
   real(kind=8), dimension(:)   , allocatable :: upward_lw_incid
   real(kind=8), dimension(:)   , allocatable :: upward_lw_surf
   real(kind=8), dimension(:)   , allocatable :: source_lw
   real(kind=8), dimension(:)   , allocatable :: forcing_lw
   real(kind=8), dimension(:)   , allocatable :: A_dw
   real(kind=8), dimension(:)   , allocatable :: B_dw
   real(kind=8), dimension(:)   , allocatable :: C_dw
   real(kind=8), dimension(:)   , allocatable :: D_dw
   real(kind=8), dimension(:)   , allocatable :: A_uw
   real(kind=8), dimension(:)   , allocatable :: B_uw
   real(kind=8), dimension(:)   , allocatable :: C_uw
   real(kind=8), dimension(:)   , allocatable :: D_uw
   real(kind=8), dimension(:)   , allocatable :: E_uw
   real(kind=8), dimension(:)   , allocatable :: F_uw
   real(kind=4), dimension(:)   , allocatable :: lw_v_surf_layer
   real(kind=4), dimension(:)   , allocatable :: lw_v_incid_layer


   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine allocates the scratch arrays.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_canopy_layer()
      implicit none
      !------------------------------------------------------------------------------------!

      allocate(zztop                (  ncanlyr)          )
      allocate(zzmid                (  ncanlyr)          )
      allocate(zzbot                (  ncanlyr)          )
      allocate(dzcan                (  ncanlyr)          )
      allocate(zztop8               (  ncanlyr)          )
      allocate(zzmid8               (  ncanlyr)          )
      allocate(zzbot8               (  ncanlyr)          )
      allocate(dzcan8               (  ncanlyr)          )

      allocate(opencan              (  ncanlyr)          )
      allocate(opencan8             (  ncanlyr)          )


      allocate(lad                  (  ncanlyr)          )
      allocate(cdrag                (  ncanlyr)          )
      allocate(pshelter             (  ncanlyr)          )
      allocate(cumldrag             (  ncanlyr)          )
      allocate(windlyr              (  ncanlyr)          )
      allocate(windext_full         (  ncanlyr)          )
      allocate(windext_half         (  ncanlyr)          )
      allocate(lad8                 (  ncanlyr)          )
      allocate(cdrag8               (  ncanlyr)          )
      allocate(pshelter8            (  ncanlyr)          )
      allocate(cumldrag8            (  ncanlyr)          )
      allocate(windlyr8             (  ncanlyr)          )
      allocate(windext_full8        (  ncanlyr)          )
      allocate(windext_half8        (  ncanlyr)          )

      allocate(indx                 (ncanlyrt2)          )
      allocate(populated            (  ncanlyr)          )
      allocate(dzcanpop             (  ncanlyr)          )
      allocate(layer_scatter        (  ncanlyr)          )
      allocate(layer_backscatter    (  ncanlyr)          )
      allocate(layer_clumping       (  ncanlyr)          )
      allocate(expkl_top            (  ncanlyr)          )
      allocate(expkl_bot            (  ncanlyr)          )
      allocate(expamk_top           (  ncanlyr)          )
      allocate(expamk_bot           (  ncanlyr)          )
      allocate(expapk_top           (  ncanlyr)          )
      allocate(expapk_bot           (  ncanlyr)          )
      allocate(A_top                (  ncanlyr)          )
      allocate(A_bot                (  ncanlyr)          )
      allocate(B_top                (  ncanlyr)          )
      allocate(B_bot                (  ncanlyr)          )
      allocate(C_top                (  ncanlyr)          )
      allocate(C_bot                (  ncanlyr)          )
      allocate(F_top                (  ncanlyr)          )
      allocate(F_bot                (  ncanlyr)          )
      allocate(G_top                (  ncanlyr)          )
      allocate(G_bot                (  ncanlyr)          )
      allocate(H_top                (  ncanlyr)          )
      allocate(H_bot                (  ncanlyr)          )
      allocate(beam_bot             (  ncanlyr)          )
      allocate(beam_mid             (  ncanlyr)          )
      allocate(beam_bot_crown       (  ncanlyr)          )
      allocate(upward_vis_beam      (ncanlyrp1)          )
      allocate(upward_vis_diffuse   (ncanlyrp1)          )
      allocate(upward_nir_beam      (ncanlyrp1)          )
      allocate(upward_nir_diffuse   (ncanlyrp1)          )
      allocate(downward_nir_beam    (ncanlyrp1)          )
      allocate(downward_nir_diffuse (ncanlyrp1)          )
      allocate(downward_vis_beam    (ncanlyrp1)          )
      allocate(downward_vis_diffuse (ncanlyrp1)          )
      allocate(mastervec_beam       (ncanlyrt2)          )
      allocate(masveccp_beam        (ncanlyrt2)          )
      allocate(mastervec_diffuse    (ncanlyrt2)          )
      allocate(masveccp_diffuse     (ncanlyrt2)          )
      allocate(PAR_beam_layer       (  ncanlyr)          )
      allocate(PAR_diffuse_layer    (  ncanlyr)          )
      allocate(SW_abs_beam_layer    (  ncanlyr)          )
      allocate(SW_abs_diffuse_layer (  ncanlyr)          )

      allocate(layer_emis           (  ncanlyr)          )
      allocate(layer_temp           (  ncanlyr)          )
      allocate(mastervec_surf       (ncanlyrt2)          )
      allocate(mastervec_incid      (ncanlyrt2)          )
      allocate(explai               (ncanlyrp1)          )
      allocate(exmlai               (ncanlyrp1)          )
      allocate(downward_lw_incid    (ncanlyrp1)          )
      allocate(downward_lw_surf     (ncanlyrp1)          )
      allocate(upward_lw_incid      (ncanlyrp1)          )
      allocate(upward_lw_surf       (ncanlyrp1)          )
      allocate(source_lw            (ncanlyrt2)          )
      allocate(forcing_lw           (ncanlyrt2)          )
      allocate(A_dw                 (ncanlyrt2)          )
      allocate(B_dw                 (ncanlyrt2)          )
      allocate(C_dw                 (ncanlyrt2)          )
      allocate(D_dw                 (ncanlyrt2)          )
      allocate(A_uw                 (ncanlyrt2)          )
      allocate(B_uw                 (ncanlyrt2)          )
      allocate(C_uw                 (ncanlyrt2)          )
      allocate(D_uw                 (ncanlyrt2)          )
      allocate(E_uw                 (ncanlyrt2)          )
      allocate(F_uw                 (ncanlyrt2)          )
      allocate(lw_v_surf_layer      (ncanlyrt2)          )
      allocate(lw_v_incid_layer     (ncanlyrt2)          )


      allocate(matal                (ncanlyrt2,        2))
      allocate(mastermat            (ncanlyrt2,        5))
      allocate(masmatcp             (ncanlyrt2,ncanlyrt2))

      return
   end subroutine alloc_canopy_layer
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine resets the canopy layer scratch variables.                       !
   !---------------------------------------------------------------------------------------!
   subroutine zero_canopy_layer(thiscall)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*), intent(in) :: thiscall
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Select which variables to flush to zero based on the call.                     !
      !------------------------------------------------------------------------------------!
      select case(trim(thiscall))
      case ('canopy_turbulence')
         opencan              (:)   = 0.
         lad                  (:)   = 0.
         cdrag                (:)   = 0.
         pshelter             (:)   = 0.
         cumldrag             (:)   = 0.
         windlyr              (:)   = 0.
         windext_full         (:)   = 0.
         windext_half         (:)   = 0.

      case ('canopy_turbulence8')
         opencan              (:)   = 0.d0
         lad8                 (:)   = 0.d0
         cdrag8               (:)   = 0.d0
         pshelter8            (:)   = 0.d0
         cumldrag8            (:)   = 0.d0
         windlyr8             (:)   = 0.d0
         windext_full8        (:)   = 0.d0
         windext_half8        (:)   = 0.d0

      case ('sw_twostream_layer')
         indx                 (:)   = 0
         populated            (:)   = .false.
         opencan8             (:)   = 0.d0
         dzcanpop             (:)   = 0.d0
         layer_scatter        (:)   = 0.d0
         layer_backscatter    (:)   = 0.d0
         layer_clumping       (:)   = 0.d0
         expkl_top            (:)   = 0.d0
         expkl_bot            (:)   = 0.d0
         expamk_top           (:)   = 0.d0
         expamk_bot           (:)   = 0.d0
         expapk_top           (:)   = 0.d0
         expapk_bot           (:)   = 0.d0
         A_top                (:)   = 0.d0
         A_bot                (:)   = 0.d0
         B_top                (:)   = 0.d0
         B_bot                (:)   = 0.d0
         C_top                (:)   = 0.d0
         C_bot                (:)   = 0.d0
         F_top                (:)   = 0.d0
         F_bot                (:)   = 0.d0
         G_top                (:)   = 0.d0
         G_bot                (:)   = 0.d0
         H_top                (:)   = 0.d0
         H_bot                (:)   = 0.d0
         beam_bot             (:)   = 0.d0
         beam_mid             (:)   = 0.d0
         beam_bot_crown       (:)   = 0.d0
         upward_vis_beam      (:)   = 0.d0
         upward_vis_diffuse   (:)   = 0.d0
         upward_nir_beam      (:)   = 0.d0
         upward_nir_diffuse   (:)   = 0.d0
         downward_nir_beam    (:)   = 0.d0
         downward_nir_diffuse (:)   = 0.d0
         downward_vis_beam    (:)   = 0.d0
         downward_vis_diffuse (:)   = 0.d0
         mastervec_beam       (:)   = 0.d0
         masveccp_beam        (:)   = 0.d0
         mastervec_diffuse    (:)   = 0.d0
         masveccp_diffuse     (:)   = 0.d0
         PAR_beam_layer       (:)   = 0.
         PAR_diffuse_layer    (:)   = 0.
         SW_abs_beam_layer    (:)   = 0.
         SW_abs_diffuse_layer (:)   = 0.
         matal                (:,:) = 0.d0
         mastermat            (:,:) = 0.d0
         masmatcp             (:,:) = 0.d0

      case ('lw_twostream_layer')
         indx                 (:)   = 0
         populated            (:)   = .false.
         layer_emis           (:)   = 0.d0
         layer_temp           (:)   = 0.d0
         mastervec_surf       (:)   = 0.d0
         mastervec_incid      (:)   = 0.d0
         explai               (:)   = 0.d0
         exmlai               (:)   = 0.d0
         downward_lw_incid    (:)   = 0.d0
         downward_lw_surf     (:)   = 0.d0
         upward_lw_incid      (:)   = 0.d0
         upward_lw_surf       (:)   = 0.d0
         source_lw            (:)   = 0.d0
         forcing_lw           (:)   = 0.d0
         A_dw                 (:)   = 0.d0
         B_dw                 (:)   = 0.d0
         C_dw                 (:)   = 0.d0
         D_dw                 (:)   = 0.d0
         A_uw                 (:)   = 0.d0
         B_uw                 (:)   = 0.d0
         C_uw                 (:)   = 0.d0
         D_uw                 (:)   = 0.d0
         E_uw                 (:)   = 0.d0
         F_uw                 (:)   = 0.d0
         lw_v_surf_layer      (:)   = 0.
         lw_v_incid_layer     (:)   = 0.
         matal                (:,:) = 0.d0
         mastermat            (:,:) = 0.d0
         masmatcp             (:,:) = 0.d0

      end select
      !------------------------------------------------------------------------------------!
      return
   end subroutine zero_canopy_layer
   !=======================================================================================!
   !=======================================================================================!
end module canopy_layer_coms
!==========================================================================================!
!==========================================================================================!
