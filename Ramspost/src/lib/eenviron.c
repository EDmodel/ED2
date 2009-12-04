/*############################ Change Log ##################################
! 1.0.0.3
!
! 010809 CJT egetenv ##
!            Replaced index and strncasecmp since Windows doesn't have them. ##
! 000908 MJB egetenv ##
!            Changed warning messages. ##
! 000829 BT  many ##
!            Eliminated compiler warnings. ##
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!#########################################################################*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
char *egetenv(char *);

#ifndef PROGRAM
#include "utils_sub_names.h"
#endif

#define MAXLINELEN 100

#ifndef ENVNAME
#  define ENVNAME  "RAMS_PATHS"
#endif

#ifdef PROGRAM
main(argc,argv)
  int argc; char *argv[];
{
  extern char *egetenv(char *);
  if (argc < 2  || argc > 2)
  {
    fprintf(stderr,"Usage: %s VARNAME\n",argv[0]);
    exit(1);
  }
  printf("%s\n",egetenv(argv[1]));
  exit(0);
}
#endif

int bills_strncasecmp(const char *s1, const char *s2, int n) {

  char *c1,*c2;
  int  nn;

  for (nn=0; nn<n ; nn++) {
    c1=(char *)s1+nn;
    c2=(char *)s2+nn;
    if (tolower(*c1) != tolower(*c2)) return(1);
    if (*c1=='\0' || *c2=='\0') return(1);
  }
  return(0);
}

#ifdef SGI
char *getenv();
char *index();
#endif
#ifdef IBM
char *getenv();
#endif


char *egetenv(variable)
  char *variable;   /* variable to find value of */
{
  static char value[100];
  char var[40] , line[100], *env;
  FILE *pfile;
  int maxlen;
  extern char *getenv(char*);

  /* First check for an environmental variable */
  if ( (env=getenv(variable)) != NULL)
  {
    strncat(value,env,MAXLINELEN-1);
    return(value);
  }

  /* Now look in our Paths file for this variable */
  if ( (env=getenv(ENVNAME)) == NULL)
  {
    fprintf(stderr,"WARNING: environmental variable %s is not set\n",ENVNAME);
    return(NULL);
  }
    
  if ( (pfile=fopen(env,"r")) == NULL)
  {
    fprintf(stderr,"WARNING: could not open file %s\n",env);
    fprintf(stderr,"  Check value of environmental variable %s\n",ENVNAME);
    return(NULL);
  }

  while(fscanf(pfile,"%[^\n]\n",line) != EOF)  /* reads the whole line */
  {
    /* printf("line: %s\n",line);*/
    if (line[0] == '#') continue;                /* Skip comments */
    if (bills_strncasecmp(line,"export ",7)==0) continue;  /* Skip export statement */

    /* We have a path setting line so split the line */
    sscanf(line," %[^\040=] = %s",var,value); /*first reads to blank or = */

    /* Now check for our variable without trailing blanks */
    if (strlen(variable) > strlen(var))
      maxlen = strlen(variable);
    else
      maxlen = strlen(var);
    if (strncmp(var,variable,maxlen) != 0) continue;

    fclose(pfile);
    return(value);
  }
  fclose(pfile);
  return(NULL);
}

/* FORTRAN INTERFACE */
void fegetenv(var,env,ierr,vlen,elen)
  char *var, *env;
  int  *ierr,vlen,elen;
{
  char tvar[40],*val;
  #ifdef PC_NT1
  extern void exit(int);
  #endif
  
  /* fix up the variable string */
  strncpy(tvar,var,vlen);
  tvar[vlen]='\0';
  val = strchr(tvar,' ');
  if (val != NULL) *val = '\0';

  if ( (val = egetenv(tvar)) != NULL)
  {
    if ( (int)strlen(val) > elen-1)
    {
      fprintf(stderr,"Error: fgetenv: length of env variable too small\n");
      exit(1);
    }
    strcpy(env,val);

    /* now fix up the outgoing string */
    val = (char *)strchr(env,'\0');
    for ( ; val<env+elen; val++) *val=' ';
    *ierr = 0;
  }
  else
  *ierr = 1;
}
