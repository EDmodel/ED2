!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 08/31/08 Rearranging the cumulus structure, eliminating grell and shcu structures. !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
subroutine rams_mem_alloc(proc_type)

   use mem_all
   use node_mod
   use leaf_coms, only : alloc_leafcol

   use io_params, only : maxlite,lite_vars,nlite_vars,frqlite,frqanl

   !----- Needed for Grell Cumulus parameterization ---------------------------------------!
   use grell_coms, only:        & ! Scalar parameters:
           closure_type         & ! - intent(in)
          ,mgmzp                & ! - intent(out)
          ,maxens_lsf           & ! - intent(inout)
          ,maxens_eff           & ! - intent(inout)
          ,maxens_cap           & ! - intent(inout)
          ,maxens_dyn           & ! - intent(inout)
          ,define_grell_coms    ! ! - subroutine
            
   use mem_scratch_grell, only: & ! Grell's simple scratch variables:
           alloc_scratch_grell  & !  - subroutine
          ,zero_scratch_grell     !  - subroutine

   use mem_ensemble, only:         & ! Ensemble scratch variables:
           ensemble_e             & !  - intent(out)
          ,alloc_ensemble         & !  - subroutine
          ,nullify_ensemble       & !  - subroutine
          ,zero_ensemble            !  - subroutine

   ! Needed for Optmization in HTINT
   use mem_opt

   ! Needed for CATT
   use catt_start, only: CATT           ! intent(in)
   use mem_carma
   use mem_aerad, only: &
        nwave,          &           !INTENT(IN)
        initial_definitions_aerad !Subroutine
   use mem_globaer, only: &
        initial_definitions_globaer !Subroutine
   use mem_globrad, only: &
        initial_definitions_globrad !Subroutine
   use extras
   use mem_turb_scalar
   !25082004


   ! For specific optimization depending the type of machine
   use machine_arq, only: machine ! INTENT(IN)

   ! Global grid dimension definitions
   use mem_grid_dim_defs, only: define_grid_dim_pointer ! subroutine

   use mem_scalar

   ! TEB_SPM
   use teb_spm_start, only: TEB_SPM ! INTENT(IN)

   use mem_emiss, only: isource ! INTENT(IN)

   use mem_gaspart, only: &
        gaspart_g,        & ! INTENT(IN)
        gaspartm_g,       & ! INTENT(IN)
        gaspart_vars,     & ! Type
        nullify_gaspart,  & ! Subroutine
        alloc_gaspart,    & ! Subroutine
        filltab_gaspart     ! Subroutine

   use mem_teb_common, only: &
        tebc_g,              & ! INTENT(IN)
        tebcm_g,             & ! INTENT(IN)
        nullify_tebc,        & ! Subroutine
        alloc_tebc,          & ! Subroutine
        filltab_tebc           ! Subroutine

   use teb_vars_const, only: iteb ! INTENT(IN)

   use mem_teb, only: &
        teb_g,        & ! INTENT(IN)
        tebm_g,       & ! INTENT(IN)
        nullify_teb,  & ! Subroutine
        alloc_teb,    & ! Subroutine
        filltab_teb     ! Subroutine

   use mem_mass, only :  mass_g,        & !subroutine
                         massm_g,       & !subroutine
                         define_frqmassave, & ! subroutine 
                         nullify_mass,  & !subroutine
                         alloc_mass,    & !subroutine
                         zero_mass,     & !subroutine
                         filltab_mass     !subroutine

   ! Needed for Nakanishi and Niino turbulence closure
   use turb_constants, only: assign_const_nakanishi
   
   
   ! For old Grell 
   use mem_grell_param       ! scalar parameters
   use mem_scratch1_grell    ! scratch 1
   use mem_scratch2_grell    ! scratch 2
   use mem_scratch3_grell    ! scratch 3
   use mem_scratch2_grell_sh ! scratch 2
   use mem_scratch3_grell_sh ! scratch 3

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                        , intent(in) :: proc_type
   !----- Local Variables: ----------------------------------------------------------------!
   integer, pointer , dimension(:)             :: nmzp,nmxp,nmyp
   integer                                     :: ng,nv,imean,na,ne
   logical                                     :: Alloc_Old_Grell_Flag
   !----- Local variables because of TEB_SPM ----------------------------------------------!
   type(gaspart_vars), pointer :: gaspart_p
   !---------------------------------------------------------------------------------------!


   write (unit=*,fmt=*) '-----------------------------------------------------------------'
   write (unit=*,fmt=*) ' BRAMS memory allocation:'

   !----- First, depending on type of process, define grid point pointers correctly. ------!
   nullify(gaspart_p)

   select case (proc_type)
   case (0,1)
      !----- This is the call for either a single processor run or for the master process -!
      nmzp => nnzp
      nmxp => nnxp
      nmyp => nnyp
   case (2)
      !----- This is the call for a initial compute node process --------------------------!
      nmzp => mmzp
      nmxp => mmxp
      nmyp => mmyp
   case (3)
      !----- This is the call for a dynamic balance compute node process. -----------------!
      nmzp => mmzp
      nmxp => mmxp
      nmyp => mmyp
      call dealloc_all()
   end select


   !----- Call global grid dimension definitions. -----------------------------------------!
   call define_grid_dim_pointer(proc_type,ngrids,maxgrds,nnzp,nnxp,nnyp,mmzp,mmxp,mmyp)

   !----- If we are doing time-averaging for output, set flag... --------------------------!
   imean=0
   if (avgtim /= 0.) imean=1

   !----- Allocate universal variable tables. ---------------------------------------------!
   allocate (num_var(maxgrds))
   allocate (vtab_r(maxvars,maxgrds))

   num_var(1:maxgrds)=0
   nvgrids=ngrids

   !----- Allocate scalar table. ----------------------------------------------------------!
   allocate(num_scalar(maxgrds))
   allocate(scalar_tab(maxsclr,maxgrds))
   num_scalar(1:maxgrds)=0

   !----- Allocate Basic variables data type. ---------------------------------------------!
   write (unit=*,fmt=*) ' [+] Basic allocation on node ',mynum,'...'
   allocate(basic_g(ngrids),basicm_g(ngrids))
   do ng=1,ngrids
      call nullify_basic(basic_g(ng))
      call nullify_basic(basicm_g(ng))
      call alloc_basic(basic_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
      if (imean == 1) then
         call alloc_basic(basicm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
      elseif (imean == 0) then
         call alloc_basic(basicm_g(ng),1,1,1)
      endif
      call filltab_basic(basic_g(ng),basicm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !---------------------------------------------------------------------------------------!



   !----- Allocate Cuparm variables data type. --------------------------------------------!
   write (unit=*,fmt=*) ' [+] Cuparm allocation on node ',mynum,'...'
   !---------------------------------------------------------------------------------------!
   !    Before I proceed, check how many cloud types I need to allocate. This may redefine !
   ! some variables given by BRAMS just to make sure the flags are consistent.             !
   !---------------------------------------------------------------------------------------!
   call define_cumulus_dimensions(ngrids)
  
   allocate(cuparm_g(ngrids),cuparmm_g(ngrids))
   do ng=1,ngrids
      call nullify_cuparm(cuparm_g(ng))
      call nullify_cuparm(cuparmm_g(ng))
      call alloc_cuparm(cuparm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,co2_on)
      if (imean == 1) then  
         call alloc_cuparm(cuparmm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,co2_on)
      elseif (imean == 0) then
         call alloc_cuparm(cuparmm_g(ng),1,1,1,ng,co2_on)
      endif
      call initialize_cuparm(cuparm_g(ng))
      call initialize_cuparm(cuparmm_g(ng))
      call filltab_cuparm(cuparm_g(ng),cuparmm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !---------------------------------------------------------------------------------------!
   !     Now I check whether I need to allocate any of the Grell structures. If so,        !
   ! proceed with the allocation.                                                          !
   !---------------------------------------------------------------------------------------!
   if (nclouds >= 3 .or. any(ndeepest(1:ngrids) == 2) .or.                                 &
       any(nshallowest(1:ngrids) == 2)) then
      write (unit=*,fmt=*) ' [+] Grell allocation on node ',mynum,'...'
      ! Calculating the necessary space for scratch data
      call define_grell_coms(ngrids,nclouds,nmzp(1:ngrids),nnqparm(1:ngrids)               &
                            ,grell_1st(1:ngrids),grell_last(1:ngrids))
      ! Initializing Grell scratch
      call alloc_scratch_grell(mgmzp)
      call zero_scratch_grell(3)

      allocate(ensemble_e(nclouds))
      do ne=1, nclouds
         call nullify_ensemble(ensemble_e(ne))
         call alloc_ensemble(ensemble_e(ne),nclouds,mgmzp,maxens_dyn,maxens_lsf,maxens_eff &
                            ,maxens_cap)
         call zero_ensemble(ensemble_e(ne))
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Allocate Cuparm variables data type. --------------------------------------------!
   Alloc_Old_Grell_Flag = .false.
   grellplusloop: do ng=1, ngrids
      if (nnqparm(ng) == 1 .and. (ndeepest(ng) == 3 .or. nshallowest(ng) == 3)) then
         Alloc_Old_Grell_Flag = .true.
         exit grellplusloop
      end if
   end do grellplusloop
   if (Alloc_Old_Grell_Flag) then
      write (unit=*,fmt=*) ' [+] ''Old Grell'' scratch allocation on node ',mynum,'...'
      call define_memory(nmxp, nmyp, nmzp, ngrids, nnqparm,ndeepest,nshallowest)

      !----- Allocating data for scratch data ---------------------------------------------!
      call alloc_scratch2_grell
      call alloc_scratch3_grell
      call alloc_scratch2_grell_sh 
      call alloc_scratch3_grell_sh 
      allocate(sc1_grell_g(ngrids))
      do ng=1,ngrids
         if (ndeepest(ng) == 3 .or. nshallowest(ng) == 3) then
            call nullify_scratch1_grell(sc1_grell_g(ng))
            call alloc_scratch1_grell(sc1_grell_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),nclouds)
            call zero_scratch1_grell(sc1_grell_g(ng))
         end if
      end do
      call zero_scratch2_grell() !LFR
      call zero_scratch3_grell() !LFR
      call zero_scratch2_grell_sh() !LFR
      call zero_scratch3_grell_sh() !LFR

      Flag_Grell = 2
      !----- Translating closure_type to icoic --------------------------------------------!
      select case (closure_type)
      case ('en')
         icoic = 0
         icoic_sh = 0
      case ('gr')
         icoic = 1
         icoic_sh = 1
      case ('lo')
         icoic = 4
         icoic_sh = 8
      case ('mc')
         icoic = 7
         icoic_sh = 8
      case ('kf') 
         icoic = 10
         icoic_sh = 8
      case ('as') 
         icoic = 13
         icoic_sh = 4
      case default 
         icoic_sh = 1
      end select
   end if

   !----- Allocate Leaf type. -------------------------------------------------------------!
   write (unit=*,fmt=*) ' [+] Leaf allocation on node ',mynum,'...'
   allocate(leaf_g(ngrids),leafm_g(ngrids))
   do ng=1,ngrids
      call nullify_leaf(leaf_g(ng)) ; call nullify_leaf(leafm_g(ng))
      call alloc_leaf(leaf_g(ng),nmzp(ng),nmxp(ng),nmyp(ng)  &
           ,nzg,nzs,npatch,ng,teb_spm)
      if (imean == 1) then
         call alloc_leaf(leafm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),nzg,nzs,npatch,ng,teb_spm)
      elseif (imean == 0) then
         call alloc_leaf(leafm_g(ng),1,1,1,1,1,1,1,teb_spm)
      end if

      call filltab_leaf(leaf_g(ng),leafm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),nzg,nzs    &
                       ,npatch,ng)
   end do
   !---------------------------------------------------------------------------------------!
   ! Bob (1/10/2002) added the following line. Is this the right place for the long term?? !
   !---------------------------------------------------------------------------------------!
   call alloc_leafcol(nzg,nzs)
   !---------------------------------------------------------------------------------------!


   !----- Allocate Micro variables data type. ---------------------------------------------!
   write (unit=*,fmt=*) ' [+] Micro allocation on node ',mynum,'...'
   allocate(micro_g(ngrids),microm_g(ngrids))
   do ng=1,ngrids
      call nullify_micro(micro_g(ng))
      call nullify_micro(microm_g(ng))
      call alloc_micro(micro_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,isfcl)
      if (imean == 1) then
         call alloc_micro(microm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,isfcl)
      elseif (imean == 0) then
         call alloc_micro(microm_g(ng),1,1,1,ng,isfcl)
      end if
      call zero_micro(micro_g(ng))
      call zero_micro(microm_g(ng))
      call filltab_micro(micro_g(ng),microm_g(ng),imean  &
           ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !---------------------------------------------------------------------------------------!



   !----- Allocate radiate variables data type --------------------------------------------!
   write (unit=*,fmt=*) ' [+] Radiate allocation on node ',mynum,'...'
   allocate(radiate_g(ngrids),radiatem_g(ngrids))
   do ng=1,ngrids
      call nullify_radiate(radiate_g(ng))
      call nullify_radiate(radiatem_g(ng))
      call alloc_radiate(radiate_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
      if (imean == 1) then
         call alloc_radiate(radiatem_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
      elseif (imean == 0) then
         call alloc_radiate(radiatem_g(ng),1,1,1,ng)
      end if
      call zero_radiate(radiate_g(ng))
      call zero_radiate(radiatem_g(ng))

      call filltab_radiate(radiate_g(ng),radiatem_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng)   &
                          ,ng)
   end do
   !---------------------------------------------------------------------------------------!


   !----- Allocate turb variables data type -----------------------------------------------!
   write (unit=*,fmt=*) ' [+] Turb allocation on node ',mynum,'...'
   allocate(turb_g(ngrids),turbm_g(ngrids))
   do ng=1,ngrids
      call nullify_turb(turb_g(ng)) ; call nullify_turb(turbm_g(ng))
      call alloc_turb(turb_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,co2_on)
      if (imean == 1) then
         call alloc_turb(turbm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,co2_on)
      elseif (imean == 0) then
         call alloc_turb(turbm_g(ng),1,1,1,ng,co2_on)
      end if
      call zero_turb(turb_g(ng))
      call zero_turb(turbm_g(ng))

      call filltab_turb(turb_g(ng),turbm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !----- Compute some Nakanishi and Niino constants --------------------------------------!
   if (any(idiffk(1:ngrids) == 7)) call assign_const_nakanishi()
   !---------------------------------------------------------------------------------------!


   !---- Allocate varinit variables data type. These never need "mean" type. --------------!
   write (unit=*,fmt=*) ' [+] Varinit allocation on node ',mynum,'...'
   allocate(varinit_g(ngrids),varinitm_g(ngrids))
   do ng=1,ngrids
      call nullify_varinit(varinit_g(ng)) 
      call nullify_varinit(varinitm_g(ng))
      call alloc_varinit(varinit_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,co2_on)
      call alloc_varinit(varinitm_g(ng),1,1,1,ng,co2_on)
      call filltab_varinit(varinit_g(ng),varinitm_g(ng),0,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !---------------------------------------------------------------------------------------!


   !----- Allocate oda variables data type. These do not need "mean" type ever. -----------!
   write (unit=*,fmt=*) ' [+] Oda allocation on node ',mynum,'...'
   allocate(oda_g(ngrids),odam_g(ngrids))
   do ng=1,ngrids
      call nullify_oda(oda_g(ng))
      call nullify_oda(odam_g(ng))
      call alloc_oda(oda_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,proc_type)
      call alloc_oda(odam_g(ng),1,1,1,ng,proc_type)
      call filltab_oda(oda_g(ng),odam_g(ng),0,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !---------------------------------------------------------------------------------------!


   !----- Allocate grid variables data type. ----------------------------------------------!
   write (unit=*,fmt=*) ' [+] Grid allocation on node ',mynum,'...'
   allocate(grid_g(ngrids),gridm_g(ngrids))
   do ng=1,ngrids
      call nullify_grid(grid_g(ng))
      call nullify_grid(gridm_g(ng))
      call alloc_grid(grid_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng,if_adap)
      call alloc_grid(gridm_g(ng),1,1,1,ng,if_adap)
      call filltab_grid(grid_g(ng),gridm_g(ng),0,nmzp(ng),nmxp(ng),nmyp(ng),ng)
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Allocate any added Scalar types.                                                      !
   ! NOT ALLOWING DIFFERENT NUMBERS OF SCALARS ON DIFFERENT NESTS.                         !
   !   Allocate length 1 of these datatypes by default                                     !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt=*) ' [+] Scalar allocation on node ',mynum,'...'
   allocate(scalar_g(1,ngrids),scalarm_g(1,ngrids))

   if (naddsc > 0) then
      !----- Deallocate datatypes, then re-alloc to correct length ------------------------!
      deallocate(scalar_g,scalarm_g)
      allocate(scalar_g(naddsc,ngrids),scalarm_g(naddsc,ngrids))
      do ng=1,ngrids
         call nullify_scalar(scalar_g(:,ng),naddsc)
         call nullify_scalar(scalarm_g(:,ng),naddsc)
         call alloc_scalar(scalar_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),naddsc)
         if (imean == 1) then
            call alloc_scalar(scalarm_g(:,ng),nmzp(ng),nmxp(ng),nmyp(ng),naddsc)
         elseif (imean == 0) then
            call alloc_scalar(scalarm_g(:,ng),1,1,1,naddsc)
         end if
      end do
      ! For CATT
      do ng=1,ngrids
         do na=1,naddsc ! For CATT
            call filltab_scalar(scalar_g(na,ng),scalarm_g(na,ng),imean,nmzp(ng),nmxp(ng)   &
                               ,nmyp(ng),ng,na)
         end do
      end do
   else !----- Necessary for CARMA Radiation, assume naddsc==3 ----------------------------!
      !----- Deallocate datatypes, then re-alloc to correct length ------------------------!
      deallocate(scalar_g,scalarm_g)
      allocate(scalar_g(3,ngrids),scalarm_g(3,ngrids))
      do ng=1,ngrids
         call nullify_scalar(scalar_g(:,ng), 3)
         call nullify_scalar(scalarm_g(:,ng), 3)
         !----- Assume grid 1x1xz ---------------------------------------------------------!
         call alloc_scalar(scalar_g(:,ng),  1, 1, 1, 3)
         call alloc_scalar(scalarm_g(:,ng), 1, 1, 1, 3)
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !----- TEB_SPM -------------------------------------------------------------------------!
   if (TEB_SPM==1) then
      write (unit=*,fmt=*) ' [+] TEB allocation on node ',mynum,'...'
      if (isource == 1) then
         !----- Allocate  gaspart vars for emission. Defining pointers --------------------!
         allocate(gaspart_g(ngrids),gaspartm_g(ngrids))
         do ng=1,ngrids
            call nullify_gaspart(gaspart_g(ng))
            call nullify_gaspart(gaspartm_g(ng))
            call alloc_gaspart(gaspart_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
            if (imean == 1) then
               call alloc_gaspart(gaspartm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
            elseif (imean == 0) then
               call alloc_gaspart(gaspartm_g(ng),1,1,1,ng)
            end if
            call filltab_gaspart(gaspart_g(ng),gaspartm_g(ng),imean,nmzp(ng),nmxp(ng)      &
                                 ,nmyp(ng),ng)
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Allocate Tendency data type, filltab_tendency is responsible for filling the main !
   ! atmospheric model variables in the scalar table, so make sure to call any routines    !
   ! that define scalar variables first.                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Assuming same scalars on all grids ----------------------------------------------!
   call nullify_tend(naddsc)
   write (unit=*,fmt=*) ' [+] Tend allocation on node ',mynum,'...'
   call alloc_tend(nmzp,nmxp,nmyp,ngrids,naddsc,proc_type)
   do ng=1,ngrids
      if (TEB_SPM==1) then
         nullify(gaspart_p)
         gaspart_p => gaspart_g(ng)
      end if
      call filltab_tend(basic_g(ng),micro_g(ng),turb_g(ng),scalar_g(:,ng),gaspart_p,naddsc &
                       ,ng)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Allocate Scratch data type, This also fills the max's that are needed by nesting   !
   ! stuff.                                                                                !
   !---------------------------------------------------------------------------------------!
   call nullify_scratch()
   write (unit=*,fmt=*) ' [+] Scratch allocation on node ',mynum,'...'
   call alloc_scratch(ngrids,nmzp(1:ngrids),nmxp(1:ngrids),nmyp(1:ngrids)                  &
                            ,nnzp(1:ngrids),nnxp(1:ngrids),nnyp(1:ngrids)                  &
                            ,nzg,nzs,npatch,nclouds,maxnxp,maxnyp,maxnzp)
   !---------------------------------------------------------------------------------------!


   !----- For CATT - LFR ------------------------------------------------------------------!
   if (CATT == 1) then
      write (unit=*,fmt=*) ' [+] CATT allocation on node ',mynum,'...'
      allocate(turb_s(ngrids))
      do ng=1,ngrids
         call nullify_turb_s(turb_s(ng))
         call alloc_turb_s(turb_s(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
         call filltab_turb_s(turb_s(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- Reproducibility - Saulo Barros --------------------------------------------------!
   call nullify_scratch1()
   write (unit=*,fmt=*) ' [+] Scratch1 allocation on node ',mynum,'...'
   call alloc_scratch1(nodebounds,maxgrds,ngrids,mmzp,mynum)
   !----- For optmization - ALF -----------------------------------------------------------!
   call nullify_opt_scratch()
   if ((if_adap==0).and.(ihorgrad==2)) then
      write (unit=*,fmt=*) ' [+] Opt_scratch allocation on node ',mynum,'...'
      call alloc_opt_scratch(proc_type,ngrids,nnzp,nnxp,nnyp,1000,1000)
   end if
   !---------------------------------------------------------------------------------------!



   !----- Allocate nested boundary interpolation arrays. All grids will be allocated. -----!

   !---------------------------------------------------------------------------------------!
   !    Changed by Alvaro L.Fazenda. To correct a problem when running in a NEC SX-6.      !
   ! Master process needs allocation for nesting in a parallel run.                        !
   !---------------------------------------------------------------------------------------!
   if (proc_type == 0 .or. proc_type == 2 .or. proc_type == 1) then
      write (unit=*,fmt=*) ' [+] Nest allocation on node ',mynum,'...'
      do ng=1,ngrids
         if(nxtnest(ng) == 0 ) then
            call alloc_nestb(ng,1,1,1)
         else
            call alloc_nestb(ng,nnxp(ng),nnyp(ng),nnzp(ng))
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- Allocate mass variables data type -----------------------------------------------!
   write (unit=*,fmt=*) ' [+] Mass allocation on node ',mynum,'...'                        
   allocate(mass_g(ngrids), massm_g(ngrids))                                               
   do ng=1, ngrids                                                                         
      call nullify_mass(mass_g(ng))                                                        
      call nullify_mass(massm_g(ng))                                                       
      call alloc_mass(mass_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),nclouds,ng,nnqparm(ng)         &
                     ,idiffk(ng))                                                          
      if (imean == 1) then                                                                 
        call alloc_mass(massm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),nclouds,ng,nnqparm(ng)      &
                       ,idiffk(ng))                                                        
      else if (imean == 0) then                                                            
        call alloc_mass(massm_g(ng),1,1,1,1,ng,nnqparm(ng),idiffk(ng))                     
      end if                                                                               
      call zero_mass(mass_g(ng))                                                           
      call zero_mass(massm_g(ng))                                                          
      call filltab_mass(mass_g(ng),massm_g(ng),imean, &                                    
                         nmzp(ng),nmxp(ng),nmyp(ng),nclouds,ng)                            
   end do                                                                                  
   ! Checking the frequency I should use for averaging                                     
   call define_frqmassave(ngrids,frqlite,frqanl,idiffk(1:ngrids),maxlite,nlite_vars        &
                         ,lite_vars,mynum)
   !---------------------------------------------------------------------------------------!



   !----- CATT - Allocation for Transport -------------------------------------------------!
   if (CATT == 1) then
      write (unit=*,fmt=*) ' [+] CATT extras allocation on node ',mynum,'...'
      allocate(extra2d (na_extra2d,ngrids))
      allocate(extra3d (na_extra3d,ngrids))
      allocate(extra2dm(na_extra2d,ngrids))
      allocate(extra3dm(na_extra3d,ngrids))
      call nullify_extra2d(extra2d,na_extra2d,ngrids)
      call nullify_extra2d(extra2dm,na_extra2d,ngrids)
      call nullify_extra3d(extra3d,na_extra3d,ngrids)
      call nullify_extra3d(extra3dm,na_extra3d,ngrids)
      do ng=1,ngrids
         call alloc_extra2d(extra2d,nmxp(ng),nmyp(ng),na_extra2d,ng)
         call zero_extra2d(extra2d,na_extra2d,ng)
         if (imean == 1) then
            call alloc_extra2d(extra2dm,nmxp(ng),nmyp(ng),na_extra2d,ng)
            call zero_extra2d(extra2dm,na_extra2d,ng)
         else
            call alloc_extra2d(extra2dm,1,1,na_extra2d,ng)
            call zero_extra2d(extra2dm,na_extra2d,ng)
         end if
         call alloc_extra3d(extra3d,nmzp(ng),nmxp(ng),nmyp(ng),na_extra3d,ng)
         call zero_extra3d(extra3d,na_extra3d,ng)
         if (imean == 1) then
            call alloc_extra3d(extra3dm,  &
                 nmzp(ng),nmxp(ng),nmyp(ng),na_extra3d,ng)
            call zero_extra3d(extra3dm,na_extra3d,ng)
         else
            call alloc_extra3d(extra3dm,1,1,1,na_extra3d,ng)
            call zero_extra3d(extra3dm,na_extra3d,ng)
         end if
      end do
      do ng=1,ngrids
         do na=1,na_extra2d
            call filltab_extra2d(extra2d(na,ng),extra2dm(na,ng),imean, &
                 nmxp(ng),nmyp(ng),ng,na)
         end do
         do na=1,na_extra3d
            call filltab_extra3d(extra3d(na,ng),extra3dm(na,ng),imean, &
                 nmzp(ng),nmxp(ng),nmyp(ng),ng,na)
         end do
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- CATT - Carma Radiation ----------------------------------------------------------!
   if (ilwrtyp==4 .or. iswrtyp==4) then
       write (unit=*,fmt=*) ' [+] CARMA allocation on node ',mynum,'...'
      call initial_definitions_aerad()
      call initial_definitions_globrad()
      call initial_definitions_globaer()

      allocate(carma(ngrids))
      allocate(carma_m(ngrids))
      do ng=1,ngrids
         call nullify_carma(carma,ng)
         call alloc_carma(carma,ng,nmxp(ng),nmyp(ng),nwave)
         call zero_carma(carma,ng)
         call nullify_carma(carma_m,ng)
         if(imean == 1) then
            call alloc_carma(carma_m,ng,nmxp(ng),nmyp(ng),nwave)
            call zero_carma(carma_m,ng)
         else
            call alloc_carma(carma_m,ng,1,1,nwave)
            call zero_carma(carma_m,ng)
         end if

         call filltab_carma(carma(ng),carma_m(ng),ng,imean,nmxp(ng),nmyp(ng),nwave)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- Some extra TEB_SPM allocation ---------------------------------------------------!
   write (unit=*,fmt=*) ' [+] Tebc allocation on node ',mynum,'...'
   !---------------------------------------------------------------------------------------!
   !     We will always allocate teb_g and tebm_g types... In case TEB is turned off, we   !
   ! won't allocate the pointers belonging to this structure. This avoids the danger of    !
   ! sending a pointer that may or may not be associated as arguments in subroutines, that !
   ! sometimes causes problems.                                                            !
   !---------------------------------------------------------------------------------------!
   allocate(tebc_g(ngrids),tebcm_g(ngrids))
   allocate(teb_g(ngrids),tebm_g(ngrids))
   if (teb_spm == 1) then
      !----- Allocate common use variables (TEB_SPM,LEAF) ---------------------------------!
      do ng=1,ngrids
         call nullify_tebc(tebc_g(ng))
         call nullify_tebc(tebcm_g(ng))
         call alloc_tebc(tebc_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
         if (imean == 1) then
            call alloc_tebc(tebcm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
         elseif (imean == 0) then
            call alloc_tebc(tebcm_g(ng),1,1,1,ng)
         end if
         call filltab_tebc(tebc_g(ng),tebcm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng)
      end do

      !----- Allocate data for urban canopy parameterization ------------------------------!
      if (iteb == 1) then
         do ng=1,ngrids
            call nullify_teb(teb_g(ng))
            call nullify_teb(tebm_g(ng))
            call alloc_teb(teb_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
            if (imean == 1) then
               call alloc_teb(tebm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
            elseif (imean == 0) then
               call alloc_teb(tebm_g(ng),1,1,1,ng)
            end if
            call filltab_teb(teb_g(ng),tebm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng)
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!


   !----- Set "Lite" variable flags according to namelist input LITE_VARS. ----------------!
   if (proc_type == 0 .or. proc_type == 2 .or. proc_type == 1) then
      write (unit=*,fmt=*) ' [+] Lite allocation on node ',mynum,'...'
      call lite_varset(proc_type)
   end if


   !---------------------------------------------------------------------------------------!
   !     Set ALL variables in the vtab_r variable table to zero by default. These are var- !
   ! iables processed in the filltab_* routines with a call to vtables2. This does NOT in- !
   ! clude scratch arrays, tendencies, or mean arrays.                                     !
   !---------------------------------------------------------------------------------------!
   do ng = 1, ngrids
      do nv = 1,num_var(ng)
         call azero( vtab_r(nv,ng)%npts, vtab_r(nv,ng)%var_p)
      enddo
   enddo
   write (unit=*,fmt=*) '-----------------------------------------------------------------'

   return
end subroutine rams_mem_alloc
!==========================================================================================!
!==========================================================================================!
