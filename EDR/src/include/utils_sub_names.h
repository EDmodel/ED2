/*!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
*/

#if defined(SUN) || defined(ALPHA) || defined(IBM) || defined (PC_LINUX1) || defined(CRAY)

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

#if defined(CRAY2)

#define fegetenv FEGETENV
#define form_tmpname FORM_TMPNAME
#define irsleep IRSLEEP
#define iralloc IRALLOC
#define irfree IRFREE
#define rams_c_open RAMS_C_OPEN
#define rams_c_close RAMS_C_CLOSE
#define rams_c_pos RAMS_C_POS
#define rams_c_tell RAMS_C_TELL
#define rams_c_read RAMS_C_READ
#define rams_c_read_char RAMS_C_READ_CHAR
#define rams_c_write RAMS_C_WRITE
#define vfirecr VFIRECR
#define vforecr VFORECR
#define readdted1 READDTED1_
#define rams_master RAMS_MASTER_
#define rams_node RAMS_NODE_
#define walltime WALLTIME_
#define par_init_fortran PAR_INIT_FORTRAN_
#define par_init_put PAR_INIT_PUT_
#define par_send PAR_SEND_
#define par_put_int PAR_PUT_INT_
#define par_put_float PAR_PUT_FLOAT_
#define par_put_char PAR_PUT_CHAR_

#define par_send_noblock PAR_SEND_NOBLOCK_

#define par_send_noblock_c2f PAR_SEND_NOBLOCK_C2F_

#define par_get_noblock PAR_GET_NOBLOCK_

#define par_get_noblock_c2f PAR_GET_NOBLOCK_C2F_

#define par_assoc_buff PAR_ASSOC_BUFF_

#define par_wait PAR_WAIT_

#define par_wait_c2f PAR_WAIT_C2F_

#define par_get_new PAR_GET_NEW_
#define par_get_int PAR_GET_INT_
#define par_get_float PAR_GET_FLOAT_
#define par_get_char PAR_GET_CHAR_
#define par_init PAR_INIT_
#define par_enroll PAR_ENROLL_
#define par_exit PAR_EXIT_
#define par_pause PAR_PAUSE_
#define par_ready PAR_READY_

#endif
