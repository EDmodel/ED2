module optimiz_coms
  use ed_max_dims, only: str_len
  
  !------------------------------------------------------------------------------!
  !    This module has currently a single variable used nowhere in the code. I   !
  ! am assuming that when the optimization code is implemented, more variables   !
  ! will be added here. If it's not the case, please move it somewhere else      !
  ! (ed_misc_coms seems a good place), and remove this module.                   !
  !------------------------------------------------------------------------------!
  
  character(len=str_len) :: ioptinpt

end module optimiz_coms
