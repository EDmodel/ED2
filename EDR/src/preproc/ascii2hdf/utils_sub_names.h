/*!===============================================================================
! OLAM version 2.10  

! Copyright (C) 2002-2006; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================*/
#if defined(SUN) || defined(ALPHA) || defined(SGI) || defined (PC_LINUX1) || defined(NEC_SX)

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
#define par_put_int par_put_int_
#define par_put_float par_put_float_
#define par_put_char par_put_char_
#define par_send_noblock par_send_noblock_
#define par_get_noblock par_get_noblock_
#define par_assoc_buff par_assoc_buff_
#define par_wait par_wait_
#define par_get_new par_get_new_
#define par_get_int par_get_int_
#define par_get_float par_get_float_
#define par_get_char par_get_char_
#define par_init par_init_
#define par_enroll par_enroll_
#define par_exit par_exit_
#define par_pause par_pause_
#define par_ready par_ready_
#define fh5f_open          fh5f_open_
#define fh5f_create        fh5f_create_
#define fh5f_close         fh5f_close_
#define fh5d_open          fh5d_open_
#define fh5d_close         fh5d_close_
#define fh5s_get_ndims     fh5s_get_ndims_
#define fh5s_get_dims      fh5s_get_dims_
#define fh5_prepare_read   fh5_prepare_read_
#define fh5d_read          fh5d_read_
#define fh5_close_read     fh5_close_read_
#define fh5_prepare_write  fh5_prepare_write_
#define fh5_write          fh5_write_
#define fh5_close_write    fh5_close_write_

#endif

#if defined(CRAY)

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
#define readdted1 READDTED1
#define rams_master RAMS_MASTER
#define rams_node RAMS_NODE
#define walltime WALLTIME
#define par_init_fortran PAR_INIT_FORTRAN
#define par_init_put PAR_INIT_PUT
#define par_send PAR_SEND
#define par_put_int PAR_PUT_INT
#define par_put_float PAR_PUT_FLOAT
#define par_put_char PAR_PUT_CHAR
#define par_send_noblock PAR_SEND_NOBLOCK
#define par_get_noblock PAR_GET_NOBLOCK
#define par_assoc_buff PAR_ASSOC_BUFF
#define par_wait PAR_WAIT
#define par_get_new PAR_GET_NEW
#define par_get_int PAR_GET_INT
#define par_get_float PAR_GET_FLOAT
#define par_get_char PAR_GET_CHAR
#define par_init PAR_INIT
#define par_enroll PAR_ENROLL
#define par_exit PAR_EXIT
#define par_pause PAR_PAUSE
#define par_ready PAR_READY
#define fh5f_open          FH5F_OPEN
#define fh5f_create        FH5F_CREATE
#define fh5f_close         FH5F_CLOSE
#define fh5d_open          FH5D_OPEN
#define fh5d_close         FH5D_CLOSE
#define fh5s_get_ndims     FH5S_GET_NDIMS
#define fh5s_get_dims      FH5S_GET_DIMS
#define fh5_prepare_read   FH5_PREPARE_READ
#define fh5d_read          FH5D_READ
#define fh5_close_read     FH5_CLOSE_READ
#define fh5_prepare_write  FH5_PREPARE_WRITE
#define fh5_write          FH5_WRITE
#define fh5_close_write    FH5_CLOSE_WRITE

#endif

#if defined(PC_NT1)

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
#define readdted1 READDTED1
#define rams_master RAMS_MASTER
#define rams_node RAMS_NODE
#define walltime WALLTIME
#define par_init_fortran PAR_INIT_FORTRAN
#define par_init_put PAR_INIT_PUT
#define par_send PAR_SEND
#define par_put_int PAR_PUT_INT
#define par_put_float PAR_PUT_FLOAT
#define par_put_char PAR_PUT_CHAR
#define par_send_noblock PAR_SEND_NOBLOCK
#define par_get_noblock PAR_GET_NOBLOCK
#define par_assoc_buff PAR_ASSOC_BUFF
#define par_wait PAR_WAIT
#define par_get_new PAR_GET_NEW
#define par_get_int PAR_GET_INT
#define par_get_float PAR_GET_FLOAT
#define par_get_char PAR_GET_CHAR
#define par_init PAR_INIT
#define par_enroll PAR_ENROLL
#define par_exit PAR_EXIT
#define par_pause PAR_PAUSE
#define par_ready PAR_READY
#define fh5f_open_         FH5F_OPEN
#define fh5f_create_       FH5F_CREATE
#define fh5f_close_        FH5F_CLOSE
#define fh5d_open_         FH5D_OPEN
#define fh5d_close_        FH5D_CLOSE
#define fh5s_get_ndims_    FH5S_GET_NDIMS
#define fh5s_get_dims_     FH5S_GET_DIMS
#define fh5_prepare_read_  FH5_PREPARE_READ
#define fh5d_read_         FH5D_READ
#define fh5_close_read_    FH5_CLOSE_READ
#define fh5_prepare_write_ FH5_PREPARE_WRITE
#define fh5_write_         FH5_WRITE
#define fh5_close_write_   FH5_CLOSE_WRITE

#endif
