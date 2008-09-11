module cup_output_vars

  implicit none

  logical          :: cup_output_vars_alloc = .false.
  real,allocatable :: xmb_ave(:), xmb_std(:)
  real,allocatable :: xmb_ske(:), xmb_cur(:)
  real,allocatable :: pr_ave(:),  pr_std(:)
  real,allocatable :: pr_ske(:),  pr_cur(:)
  real,allocatable :: x_ave(:,:), x_std(:,:)
  real,allocatable :: x_ske(:,:), x_cur(:,:), x_ave_cap(:,:)
!- srf: out2004
  real,allocatable ::  x_ave_cap1(:,:)  &
     		      ,x_ave_cap2(:,:)  &
     		      ,x_ave_cap3(:,:)  &
     		      ,x_ave_cap4(:,:)  &
     		      ,x_ave_cap5(:,:)  
!- srf ---fim----

contains

  subroutine alloc_cup_output_vars(mgmxp,maxens,maxens3)

    implicit none
    integer :: mgmxp, maxens, maxens3

    if(cup_output_vars_alloc) then
       print *,'ERROR: cut_output_vars already allocated'
       print *,'Routine: alloc_cut_output_vars File: cut_out_vars.f90'
       print *,'Dir: .../shared/tools/brams20/src/rams/5.04/modules'
       stop
    end if

    allocate(xmb_ave(mgmxp),xmb_std(mgmxp),xmb_ske(mgmxp),xmb_cur(mgmxp))
    allocate( pr_ave(mgmxp), pr_std(mgmxp), pr_ske(mgmxp), pr_cur(mgmxp))
    allocate(x_ave(mgmxp,maxens3),x_std(mgmxp,maxens3))
    allocate(x_ske(mgmxp,maxens3),x_cur(mgmxp,maxens3))
    allocate(x_ave_cap(mgmxp,maxens))

!- srf: out2004
    allocate(x_ave_cap1(mgmxp,maxens)&
     	    ,x_ave_cap2(mgmxp,maxens)&
     	    ,x_ave_cap3(mgmxp,maxens)&
     	    ,x_ave_cap4(mgmxp,maxens)&
     	    ,x_ave_cap5(mgmxp,maxens)&
	    )
!- srf ---fim----



    cup_output_vars_alloc=.true.

  end subroutine alloc_cup_output_vars

end module cup_output_vars
