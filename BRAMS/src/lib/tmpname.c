/*!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
*/

#include "utils_sub_names.h"

/*************************************************************************/
  


int form_tmpname(char *filename,int len1)

{
  int i, ierr;
  extern int mkstemp(char*);

/*
 *   printf("C_tmpnam - %d %s \n",len1,filename);
 */
  ierr = mkstemp(filename);
/*
 *   printf("C_tmpnam - %d %s \n",len1,filename);
 */
  for (i=strlen(filename); i<len1; i++) {
    *(filename+i) = ' ';
  }
  return(ierr);
}
