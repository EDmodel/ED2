/*!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
*/

#include "utils_sub_names.h"
/*
     ***************************************************************
     *                                                             *
     * Copyright: Science Applications International Corporation,  *
     *            1710 Goodridge Dr., McLean, VA 22102. (C) 1994.  *
     *                                                             *
     *               >>---->>>> WARNING <<<<----<<                 *
     *  OMEGA has been provided to the  US Defense Nuclear Agency  *
     *  under  Contract  DNA001-92-C-0076  with a restriction for  *
     *  official government use only.   The dissemination of this  *
     *  software  to  non-governmental  entities  must conform to  *
     *  this restriction.  SAIC reserves all commercial rights to  *
     *  OMEGA.                                                     *
     *                                                             *
     ***************************************************************

      g_rdd0.c
      g_rdd0.c
      g_rdd0.c
 ---------------------------------------------------------------------
       <MODULE>
 
      This module contains the C routines that read dted database
      files.

      Most recent code developement done by Paul Boris and Mary Hall
 
 ---------------------------------------------------------------------
*/

/* Application to expand block of data from character to integer */

#include <ctype.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

/* #include "arg_pars.h" */

/*****************************************************************
        Include file that does argument parsing.
        To be used in the following way:
        #include "arg_parse.h"

        In the main routine you can check for any input argument
        provided the user has used a flag preceded by a "-".
        Ex.     -size <depth><width>
        required on input can be checked for in the program
        by doing the following:
        
        int depth,width;
        int n;
        
        TEST_ARG("-size",n);

        This will search the command line arguments for "-size"
        and return the integer "n" . "n" will contain the number
        of values given with the argument "-size".

        Upon return "n" will contain
                -1  if the command line did not have a -size
                 0  if the command line had the -size argument
                    without any other values 
                 1  if the command line contained -size followed
                    by one value (ie. depth)
                 2  if the command line contained -size followed
                    by two values (ie. depth and width)

        Once the value of "n" is known use GETARG_? to read in the correct
        number of integers, characters, floats, long integers or strings.
                 
        TEST_ARG("-size",n);
        GETARG_I(depth);
        GETARG_I(width);

        Will read in the depth and width in the above example.
        Or to do error checking on input do the following:

        TEST_ARG("-size",n);
        switch(n)
                {
                case 0:
                        error(1,"Must give depth and width"); break;
                case 1:
                        width = GETARG_I; width = depth;break;
                case 2:
                        width = GETARG_I; depth = GETARG_I;break;
                case -1:
                        depth = width = 0;break;
                default:
                        error(2,"Too many values");break;
                }
        
        In this example default values can be set by testing the "n"
        value returned from TEST_ARG.
*******************************************************************/
#ifndef _ARG_PARSE_H
#define _ARG_PARSE_H

/* #include <stdio.h> */ 

double atof();
                                
char **ARGV;
int  ARGC;
#define  GETARG_I (atoi(argv[ARGC++]))/*read integer from command line*/
#define  GETARG_F (atof(argv[ARGC++]))/*read float from command line*/
#define  GETARG_L (atol(argv[ARGC++]))/*read long int from command line*/
#define  GETARG_S (argv[ARGC]==NULL?NULL:strdup(argv[ARGC++])) /*read string from command line*/

#define  CMP_ARG(FLAG)  ((StringCompare(FLAG,ARGV[ARGC]) == 0) ? 1 : 0)
#define  IS_ARG(X)      (ARGV[(X)][0] == '-' && !isdigit(ARGV[X][1]))
#define  COUNT_ARG(X)   while((++ARGC < argc) && !IS_ARG(ARGC))X++;\
                                 if(X++ >= 0) ARGC -= X;
#define  TEST_ARG(FLAG,X)  ARGC = 0; ARGV = argv; X = -1;\
                           while(ARGC++ < argc)\
                                if(CMP_ARG(FLAG)){COUNT_ARG(X);break;}

#define  DOES_NOT_OCCUR   -1
#define  FLAG_ONLY         0

#endif

/*****************************************************************************
                          THE END OF EVERYTHING
*****************************************************************************/
/* #include "fixedstr.h" */


/********************* FixedString.h **************************/

#ifndef _FIXEDSTRING_H
#define _FIXEDSTRING_H
/* #include <string.h> */

#define StringCompare( s1,s2 )\
	( (s1) != NULL && (s2) != NULL ? strcmp( s1,s2 ) : -2 )
#define StringCompareN( s1,s2,n )\
	( (s1) != NULL && (s2) != NULL ? strncmp( s1,s2,n ) : -2 )

#define StringCopy( s1,s2 )\
	( (s1) != NULL && (s2) != NULL ? strcpy( s1,s2 ) : (s1) )
#define StringCopyN( s1,s2 )\
	( (s1) != NULL && (s2) != NULL ? strncpy( s1,s2,n ) : (s1) )

#define StringConcat( s1,s2 )\
	( (s1) != NULL && (s2) != NULL ? strcat( s1,s2 ) : (s1) )
#define StringConcatN( s1,s2,n )\
	( (s1) != NULL && (s2) != NULL ? strncat( s1,s2,n ) : (s1) )

#define StringLastOccur( s1,c )\
	( (s1) != NULL ? strrchr( (s1), (c) ) : NULL )
#define StringFirstOccur( s1,c )\
	( (s1) != NULL ? strchr( (s1), (c) ) : NULL )

#define StringLength( s1 )\
	( (s1) != NULL ? strlen(s1) : 0 )
#define StringDup( s1 )\
	( (s1) != NULL ? strdup( s1 ) : (s1) )

#endif	/* end of FixedString.h */


/* #include "util.h"  */

#ifndef _util_h
#define _util_h
#define PRIVATE static
#define EXIT_SUCCESS 0
#define EXIT_FAILU 1
#define BOOLEAN short int
#define LONG   long
#define TRUE 1
#define FALSE 0
#define MAX_PATH			128
#define min(a,b)  		( (a) < (b) ? (a) : (b) )
#define max(a,b)  		( (a) > (b) ? (a) : (b) )
#define nint(a)			( (int) ( (a) + 0.5 ) )

#ifndef BORLAND_C
#define Huge
#define huge
#define farcalloc calloc
#define farfree free
#define FAR
#define NEAR
#define far
#define near
#define SEEK_CUR    1
#define SEEK_END    2
#define SEEK_SET    0
#endif

#endif
#define MEAN(x,n,mean,n2) 	( ( n * mean + x ) / (n2) )

  BOOLEAN  input=0,output=0;
  char input_file[128],output_file[128];
  FILE *fpin,*fpout;


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float readdted1(int *iwres, float d[360000])
{
	FILE *fpin, *fpat;
        char fname[80];
        short i,j;
	short z[360000];
        int ik, jl;


   fpat = fopen("namefils.out","r");
   fscanf(fpat,"%s",fname);
/*    printf(" ***%s***\n ",fname);    */
   fclose(fpat);       

/*        fgets(fname, 8, fpa1); */
/*        fread(fname,sizeof(char),8,fpa1);  */

        
        jl = 600;
        ik = 600/ *iwres;

	if(( fpin = fopen(fname,"rb")) == NULL)
	{
		fprintf(stderr," Data file %s not found\n",input_file);
	        /*fprintf(stderr," Data file %s not found\n",fname);*/
                *iwres= -100000;
                return (-1.);
	}

	fread(z,sizeof(short),ik*jl,fpin);   

        for (j=0;j<600;j++) {
	   for (i=0;i<ik;i++) {
                  jl= j*ik + i; 
                  d[jl] = (float)z[jl];
		}
	    }
       	fclose(fpin);
        return (0);
} 
