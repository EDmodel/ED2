/*!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
*/

#if defined(LINUX)

#define fegetenv fegetenv_
#define form_tmpname form_tmpname_
#define irsleep irsleep_
#define iralloc iralloc_
#define irfree irfree_
#define rams_c_open rams_c_open_
#define rams_c_close rams_c_close_
#define rams_c_pos rams_c_pos_
#define rams_c_tell rams_c_tell_
#define rams_c_write rams_c_write_
#define rams_c_read rams_c_read_
#define rams_c_read_char rams_c_read_char_
#define vfirecr vfirecr_
#define vforecr vforecr_
#define readdted1 readdted1_
#define rams_master rams_master_
#define rams_node rams_node_
#define walltime walltime_
#define par_init_fortran par_init_fortran_
#define par_init_put par_init_put_
#define par_send par_send_
#define par_bsend par_bsend_
#define par_put_int par_put_int_
#define par_put_float par_put_float_
#define par_put_char par_put_char_
#define par_put_double par_put_double_

#define par_send_noblock par_send_noblock_
#define par_send_noblock_c2f par_send_noblock_c2f_
#define par_get_noblock par_get_noblock_
#define par_get_noblock_c2f par_get_noblock_c2f_
#define par_assoc_buff par_assoc_buff_
#define par_wait par_wait_
#define par_wait_c2f par_wait_c2f_
#define par_get_new par_get_new_
#define par_get_int par_get_int_
#define par_get_float par_get_float_
#define par_get_char par_get_char_
#define par_get_double par_get_double_
#define par_init par_init_
#define par_enroll par_enroll_
#define par_exit par_exit_
#define par_error par_error_
#define par_pause par_pause_
#define par_ready par_ready_

#endif
