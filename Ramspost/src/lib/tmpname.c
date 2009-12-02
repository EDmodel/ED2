/*############################ Change Log ##################################
! 1.0.0.1
!
! 000829 BT  many ##
!            Eliminated compiler warnings. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!#########################################################################*/

#include "utils_sub_names.h"

/*************************************************************************/
  
#ifdef STARDENT 

int form_tmpname(filename)
     struct  { char *string; int len; } *filename;
{
  extern void tmpnam(char*);

  tmpnam(filename);
  /*printf("C_tmpnam - %s \n",filename->string);*/
  return(0);
}

#else

int form_tmpname(char *filename,int len1)

{
  int i;
  extern void tmpnam(char*);

  tmpnam(filename);
  /*printf("C_tmpnam - %s \n",filename);*/
  for (i=strlen(filename); i<len1; i++) {
    *(filename+i) = ' ';
  }
  return(0);
}
#endif
