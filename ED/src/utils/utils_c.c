/*!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
*/

#include "utils_sub_names.h"

#include <stdio.h>
#include <math.h>
#ifdef CRAY
#include <stdlib.h>
#endif
#ifdef IBM
#include <malloc/malloc.h>
#elif MAC_OS_X
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

/*#include <unistd.h>*/

#ifdef SGI
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#endif

/* Prototypes not needed except for C++
int vfscale(float *,int ,double *,double *);
void *malloc(int);
void *free(void *);
*/

/*********************************************************/

#define nlocrams 1000
void *iaddrams[nlocrams];

int iralloc(int *memtot,int *ia,int *ioff)
{
  extern void *iaddrams[nlocrams];
  int n, ifaddr, imaddr;
  void *iaddr;

/* Allocate the memory */

  iaddr = malloc((*memtot)*sizeof(float));

/* Compute the offset for Fortran */

#ifdef CRAY
  ifaddr = (int)ia;
  imaddr = (int)iaddr;
  *ioff=(imaddr-ifaddr);
#else
  ifaddr = (int )ia;
  imaddr = (int )iaddr;
  *ioff = (imaddr-ifaddr)/sizeof(float);
#endif

/* Find first empty location in address array */

  for(n=0; n < nlocrams ; n++) {
    if( iaddrams[n] == 0 ){
      iaddrams[n] = iaddr;
      break;
    }
  }

 /* printf("in ralloc- %i %i %i %i\n",n,*memtot,iaddr,*ioff);*/
  
/* Return this location to FORTRAN as an address "handle" */

  return(n);
}

/*********************************************************/

int irfree(int *nmem) 
{
  int ans=0;
  extern void *iaddrams[nlocrams];
  
  /* printf("in rfree- %i %i\n",*nmem,iaddrams[*nmem]); */
  free(iaddrams[*nmem]);                   
  iaddrams[*nmem]=0;
  return(ans);
  
}

/* ******************************************************* */

void irsleep(int *seconds)
{
   extern int sleep(int);

#if !defined (PC_NT1)
   sleep( *seconds );
#endif

   return;
}

/* ******************************************************* */

FILE *ramsfile;

#ifdef STARDENT 

int rams_c_open(filename,faccess)
     struct  { char *string; int len; } *filename,*faccess;
{
  extern FILE *ramsfile;
  
  /* printf(" C_open - %s %s \n",filename->string,faccess->string); */
  ramsfile=fopen(filename->string,faccess->string);
  return(0);
}
#else

int rams_c_open(char *filename,char *faccess)

{
  extern FILE *ramsfile;

  ramsfile=fopen(filename,faccess);
 /* perror("rams_c_open"); */
  return(0);
}
#endif

/*********************************************************/

int rams_c_close()
{
  extern FILE *ramsfile;
  int istat;

  istat=fclose(ramsfile);
  return(istat);
}

/*********************************************************/

int rams_c_pos(long int *fbyte)
{ 
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  return(retcode);
}

/*********************************************************/

void rams_c_tell(int *pos)
{ 
  extern FILE *ramsfile;

  *pos=ftell(ramsfile);
}

/*********************************************************/

int rams_c_read(int *fbyte,int *numbytes,int *a)
{
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  fread(a,1,*numbytes,ramsfile);
  return(retcode);
}

/*********************************************************/
int rams_c_read_char(int *fbyte,int *numbytes,int *a)
 
{
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  fread(a,1,*numbytes,ramsfile);
  return(retcode);
}

/*********************************************************/

int rams_c_write(int *fbyte,int *numbytes,int *a)
{
  int retcode;
  extern FILE *ramsfile;

  retcode=fseek(ramsfile,*fbyte,0);
  fwrite(a,1,*numbytes,ramsfile);
  return(retcode);
}

/**********************************************************************/
/*   C versions of vfirecr and vforecr written by Peter Olsson, 1993 */

#include <ctype.h>
#define BITFIELDLENGTH 6
#define SMALL_OFFSET 1.e-20
#define VFORMMASK 63 

void vfirecr(int *unit,float *a,int *n,char *type,char *b,int *irec)    
{
  extern FILE *ramsfile;
  int i, j, nn, nbits, nchs;
  float bias, fact, inverse_fact;
  unsigned vnum, char_count;

  fseek( ramsfile, *irec, 0);
  fread(b,1,80, ramsfile);
  sscanf(b,"%d %d %f %f",&nn, &nbits, &bias, &fact);
  inverse_fact = 1./fact;
  nchs=nbits/6;
  fread(b,1,*n * nchs,ramsfile);
  if( nn != *n )
    printf("Word count mismatch on vfirec record\n  Words on record - %d\n  Words expected  - %d\n ",nn,*n);
  for(i = 0, char_count = 0; i < *n; i++)
  {
    for(j = 0,vnum=0; j < nchs; j++, char_count++)
    {
      vnum = vnum << BITFIELDLENGTH;
      if( isdigit( b[char_count] ) ) 
               vnum = vnum | (unsigned)b[char_count] - 48;
      else if( isupper( b[char_count] ) ) 
               vnum = vnum |(unsigned)b[char_count] - 55;
      else 
               vnum = vnum | (unsigned)b[char_count] - 61;
    }
    a[i] = vnum*inverse_fact - bias;
  }
}

/*************************************************************************/

char vc[65] = 
     "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz{|";

void vforecr(int *unit,float *a,int *n,int *nbits,float *scr
            ,char *cscr,char *type,int *irec)
{
  extern FILE *ramsfile;
  extern char vc[];
  double amax, amin, bias, fact, divisor; 
  int i, j, char_count, nchs;
  float ftemp;
  extern int vfscale(float*, int, double*, double* );

  vfscale( a, *n, &amin, &amax);

  bias = -amin + SMALL_OFFSET;
//  printf("\t nbits %e bias %e amax %e so %e sum %e \n",(double)*nbits,bias,amax,SMALL_OFFSET,bias + amax + SMALL_OFFSET);
  /* MLO - Avoiding infinity */
  if ( bias + amax + SMALL_OFFSET != 0){
    divisor=bias+amax + SMALL_OFFSET;
  }else if (bias + amax != 0){
    divisor=bias+amax;
  }else{
    divisor=SMALL_OFFSET;
  }
  fact = (pow( 2.0, (double)*nbits)-1 ) / divisor;
  fprintf(ramsfile,"%8d%8d%20.10e%20.10e                        "
          ,*n,*nbits,bias,fact);

 /* assume for now that transformation is linear  */   

  for( i = 0; i < *n;  i++ )
    scr[i] = ( a[i] + bias ) * fact;

  nchs = *nbits / BITFIELDLENGTH;

  for( j = 0, char_count = 0 ; j < *n; j++ ) {
    ftemp=scr[j];
    for( i = nchs-1; i >= 0; i--, char_count++ )
      cscr[char_count]= 
          vc[( (unsigned)ftemp >> (i * BITFIELDLENGTH) ) & VFORMMASK];
  }

  fwrite(cscr , sizeof( char ), *n * nchs, ramsfile );
  *irec=ftell(ramsfile);
}

/*************************************************************************/

int vfscale(float *a,int n,double *min,double *max )
{ 
  int i;
  
  *min=  1.e20;
  *max = -*min;  

  for( i = 0; i < n; i++)
  {
    if( a[i] > *max ) *max = a[i];
    if( a[i] < *min ) *min = a[i];
  }
  return(0);
}

/************************************************************************/
#include <dirent.h>
#include <string.h>

void filelist_c_( int *inum, int *indices, char *prefix, char *chario, int dirlen, int charlen ){
  

  struct dirent **nameout;
  char filestr[200],filestr_p[200],filestr_a[200];
  char dir[200],tmpdir[200];
  char fpref0[80],fpref1[80],fpref2[80];
  char c1[1];
  int n,m,i;
  int good,val,num,lastlen,index,tfound;
  int nval;
  char *token,*found;

  char *delim = "/\0";
  char *delim2 = "*\0";
  
  /* First thing is to split the prefix into  a directory and a file prefix.  To do 
     this scan the string and detect the last  "/"                              */

  
  index=0;
  strcpy(tmpdir,"\0");

  /* Then we have an absolute path */
  if(strncmp(prefix,"/",1)==0){strcpy(tmpdir,"/\0");}

  
  token = strtok (prefix, delim);
  tfound=0;
  while (token !=  '\0') {
    tfound += 1;

    strcpy(dir,tmpdir);
    strcat(tmpdir,token);
    strcat(tmpdir,delim);
    strcpy(fpref0,token);

    // Fetch next token
    token = strtok('\0', delim);
  }

  /* Now we have the string parsed into the file prefix and the directory. 
     The next step is to break it into components and do a comparison with 
     the directory contents. */

  /* Find the star and break up file name. */
  tfound=0;

  strcpy(fpref1,"\0");
  strcpy(fpref2,"\0");

  if(strncmp(fpref0,"*",1)==0){
    
    /* Try the first token */
    strcpy(fpref1,"");
    tfound=1;
    
    /* Try the next token */
    token = strtok('\0',delim2);
    if (token != '\0'){
      tfound=2;
      strcpy(fpref2,token);
    }
    
    
  }
  else{


    /* Try the first token */
    token = strtok (fpref0, delim2);
    if (token != '\0'){
      tfound=1;
      strcpy(fpref1,token);
    }
    
    /* Try the next token */
    token = strtok('\0',delim2);
    if (token != '\0'){
      tfound=2;
      strcpy(fpref2,token);
    }
    
  }

  m=0;

  /*    Scan in the directory contents    */
  /****************************************/
  num = scandir(dir, &nameout, 0, alphasort);
  
  strcpy(filestr,nameout[1]->d_name);
  
  /*    Set the string vector to null     */
  strcpy(chario,"");
  
  /* Test if there are any entries, there should be at
     least two if the command was succesful         */
  if (num < 0){
    perror("scandir");
  }
  else if(num == 0) {
    
    /* Only one entry? I dont think its possible, but 
       deallocate it anyway...                       */
    free(nameout[0]);
  }
  else if(num == 1) {
    
    /* The directory was scanned, but there is nothing in it */
    /* Just the . and ..                                     */
    free(nameout[0]);
    free(nameout[1]);
  }
  else {
    
    /* So there is something in here besides the two base dirs */

    /* Set the first index of the string to 1 */
    indices[0]=1;

    /* n = 0,1 are the base dirs, skip them */
    for(n=2;n<num;n++){

      /* Copy the first file returned by scandir to filestr */
      strcpy(filestr,nameout[n]->d_name);
      
      /* Check for old files by looking for the ~ character */
      good=0;
      for(i=0;i<strlen(filestr);i++){
	val=strcmp("~",&filestr[i]);
	if(val==0){good=1;}
      }
      

      /* if the entry is not obsolete then compare it to the
	 file prefix strings fpref1 and fpref2 */
      if(good<1){
	
	
	/* Compare to prefix 1, fpref1 */
	val=0;
	if (tfound>0) {
	  val=-1;
	  val=strncmp(fpref1,filestr,strlen(fpref1));

	  

	  /* Now compare the end of the string */
	  
	  if (tfound>1 & val==0) {

	    val=-1;

	    /* Set a new variable filestr_p, to be the remaining file
	       string that has not been compared yet.              */
	    strcpy(filestr_p,"\0");
	    
	    strcat(filestr_p,&filestr[  strlen(filestr)-strlen(fpref2)   ]); 
	    strcat(filestr_p,"\0");


	    /* Search the new string for the second search prefix */
	    
	    val=strncmp(fpref2,filestr_p,strlen(fpref2));

	  }
	}
	
	
	if (val==0) {
	  
	  /* If val==0, then we matched all necessary
	     parts of the file prefix with the current
	     file. Add this strings name to the vector */

	  /* Add the vector and give it back its direcotry path*/
	  strcpy( filestr_a,"\0");
	  strcat( filestr_a,dir);
	  strcat( filestr_a,filestr);
	  strcat( filestr_a,"\0");
	  strcat( chario,filestr_a);

	  /* And give it a position index in the vector */
	  if(m>0){indices[m]=indices[m-1]+lastlen;}
	  lastlen=strlen(filestr_a);
	  
	  m++;
	}
	
      }

      free(nameout[n]);
	
    }
    //    m=m-1;
    
  }
  
  indices[m]=indices[m-1]+lastlen;
  

  /* Append a null character to the end of the string */
  strcat( chario,"\0");
  
  /* Return the number of entries */
  *inum=m;

  /* Release the scratch sting structure from memory */
  free(nameout);
    


}

/* This is for the omp thread/processor pinning check. */
#include <utmpx.h>
int sched_getcpu();
int findmycpu_ ()
{
    int cpu;
    cpu = sched_getcpu();
    return cpu;
}
