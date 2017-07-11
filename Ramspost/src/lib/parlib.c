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

#include <stdio.h>
#include "utils_sub_names.h"

#if defined (RAMS_MPI)
#include <mpi.h>
#endif

char *ibuff;
int ipos,nbuff;

/*==========================================================================*/

void par_init_put(buff,numbuff) 
     int *numbuff,buff;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos, nbuff;
#endif

#if defined (RAMS_MPI)
  ibuff=(char*)buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0 ;
#endif

}

/*==========================================================================*/

void par_send(mach,msgtype)
     int *mach, *msgtype;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos, nbuff;
#endif
  int ierr=0;
/*  int addr, nummach, mynum, msgid;*/ 

#if defined (RAMS_MPI)
/*  MPI_Request msgreq;*/ 
/*  printf(" par_send- %d %d %d %d \n", *mach, *msgtype, ibuff,ipos);*/ 
  ierr= MPI_Send( ibuff, ipos, MPI_PACKED, *mach, *msgtype, MPI_COMM_WORLD);
/*  printf(" par_sent- %d %d %d %d %d \n", *mach, *msgtype, msgreq,ipos,ierr);*/
#endif


  if(ierr < 0) 
    printf("Error in par_send - %d %d %d \n",*mach, *msgtype,ierr);

 /* printf(" par_send- %d %d %d %d \n", *mach, *msgtype, ierr,ipos);*/
}

/*==========================================================================*/

void par_put_int(iwords,numwords)
     int *iwords, *numwords;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos,nbuff;
#endif
  int ierr=0, isize=0;

#if defined (RAMS_MPI)
  ierr = MPI_Pack( iwords,*numwords,MPI_INT,ibuff,nbuff,&ipos
		  ,MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_put_int- %d %d %d \n",ierr,isize,ipos);
}

/*==========================================================================*/

void par_put_float(words,numwords)
     int *numwords;
     float *words;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos, nbuff;
#endif
  int ierr=0;

#if defined (RAMS_MPI)
  ierr = MPI_Pack( words,*numwords,MPI_FLOAT,ibuff,nbuff,&ipos
		  ,MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_put_float- %d \n",ierr);
}

/*==========================================================================*/

void par_put_chr(words,numbytes)
     int *numbytes;
     char *words;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos, nbuff;
#endif
  int ierr=0;

#if defined (RAMS_MPI)
  ierr = MPI_Pack( words,*numbytes,MPI_BYTE,ibuff,nbuff,&ipos
		  ,MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_put_float- %d \n",ierr);
}

/*==========================================================================*/

void par_send_noblock(mach,msgtype,msgtag)
     int *mach, *msgtype;
#if defined (RAMS_MPI)
MPI_Request *msgtag;
#else
int *msgtag;
#endif
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos, nbuff;
#endif
  int ierr=0;
/*  int addr, nummach, mynum, msgid, nbytes;*/

#if defined (RAMS_MPI)
/*  printf(" par_send- %d %d %d %d \n", *mach, *msgtype, ierr,ipos);*/
  ierr= MPI_Isend( ibuff, ipos, MPI_PACKED, *mach, *msgtype, MPI_COMM_WORLD
    ,msgtag);
/*  printf(" par_sent- %d %d %d %d \n", *mach, *msgtype, ierr,ipos);*/
#endif

  if(ierr < 0) 
    printf("Error in par_send_noblock - %d %d %d \n",*mach, *msgtype,ierr);

 /* printf(" par_send- %d %d %d %d \n", *mach, *msgtype, ierr,ipos);*/

}

/*==========================================================================*/

void par_get_noblock(buff,numbuff,mmtype,ihostnum,msgtag)
     int *numbuff,*mmtype,*ihostnum;
     void *buff;
#if defined (RAMS_MPI)
MPI_Request *msgtag;
#else
int *msgtag;
#endif
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos,nbuff;
  int ierr=0;
#endif

#if defined (RAMS_MPI)
/*  MPI_Status status;*/

  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;

  ierr=MPI_Irecv(ibuff,nbuff,MPI_PACKED,*ihostnum,*mmtype
	   ,MPI_COMM_WORLD,msgtag);
  if(ierr < 0) 
    printf("Error in par_get_noblock\n");
#endif

}

/*==========================================================================*/

void par_assoc_buff(buff,numbuff)
     int *numbuff;
     void *buff;
{
  extern char *ibuff;
  extern int ipos,nbuff;

  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;
}

/*==========================================================================*/

void par_wait(msgtag,ibytes,msgtype,ihostnum)
     int *ibytes,*msgtype,*ihostnum;
#if defined (RAMS_MPI)
MPI_Request *msgtag;
#else
int *msgtag;
#endif
{
#if defined (RAMS_MPI)
  MPI_Status status;
  int ierr=0;
#endif
/*  int mpl_source,nummach,mynum,msgid;*/

/*  printf("Node waiting- %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/

#if defined (RAMS_MPI)
  ierr=MPI_Wait(msgtag,&status);
  
  if(ierr < 0) 
    printf("Error in par_wait\n");
    
  MPI_Get_count(&status,MPI_PACKED,ibytes);
  *msgtype=status.MPI_TAG;
  *ihostnum=status.MPI_SOURCE;
#endif

/*printf("Node done wait- %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/

}

/*==========================================================================*/

void par_get_new(buff,numbuff,mmtype,ibytes,msgtype,ihostnum)
     int *numbuff,*mmtype,*ibytes,*msgtype,*ihostnum;
     void *buff;
{
  extern char *ibuff;
  extern int ipos,nbuff;
#if defined (RAMS_MPI)
  int ierr=0;
#endif
/*  int mpl_source,nummach,mynum,msgid;*/

#if defined (RAMS_MPI)
  MPI_Status status;
#endif

  ibuff=buff;
  nbuff=*numbuff*sizeof(float);
  ipos=0;

/* printf("Node waiting for- %d \n",*mmtype);*/

#if defined (RAMS_MPI)
  ierr=MPI_Recv(ibuff,nbuff,MPI_PACKED,MPI_ANY_SOURCE,*mmtype
	   ,MPI_COMM_WORLD,&status);
      
  if(ierr < 0) 
    printf("Error in par_get_new\n");
    
  MPI_Get_count(&status,MPI_PACKED,ibytes);
  *msgtype=status.MPI_TAG;
  *ihostnum=status.MPI_SOURCE;
#endif

/* printf("Node got- %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/

}

/*==========================================================================*/

void par_get_int(iwords,numwords)
     int *iwords,*numwords;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos,nbuff;
#endif
  int ierr=0;

/* printf("Unpack int- %d %d \n",*numwords,nbuff);*/
#if defined (RAMS_MPI)
  ierr = MPI_Unpack(ibuff,nbuff,&ipos,iwords,*numwords
	     ,MPI_INT,MPI_COMM_WORLD);
#endif

/* printf("Node got- %d %d %d %d\n",*msgtype,*ihostnum,*ibytes,mynum);*/
  if(ierr < 0) printf("Error in par_get_int-%d \n",ierr);
}

/*==========================================================================*/

void par_get_float(words,numwords)
     int *numwords;
     float *words;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos,nbuff;
#endif
  int ierr=0;

/* printf("Unpack flt- %d %d \n",*numwords,nbuff);*/
#if defined (RAMS_MPI)
  ierr= MPI_Unpack(ibuff,nbuff,&ipos,words,*numwords
		   ,MPI_FLOAT,MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_get_float-%d \n",ierr);
}

/*==========================================================================*/

void par_get_chr(words,numbytes)
     int *numbytes;
     char *words;
{
#if defined (RAMS_MPI)
  extern char *ibuff;
  extern int ipos,nbuff;
#endif
  int ierr=0;

/* printf("Unpack flt- %d %d \n",*numwords,nbuff);*/
#if defined (RAMS_MPI)
  ierr= MPI_Unpack(ibuff,nbuff,&ipos,words,*numbytes
		   ,MPI_BYTE,MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_get_float-%d \n",ierr);
}

/*==========================================================================*/

void par_init_fortran (argc,fargv,farglen,machnum,machsize)
     int *argc,*farglen; char *fargv;
     int *machnum, *machsize;
{
  int i,numarg,carglen;
  char *argvp[20];
#if defined (RAMS_MPI)
  char **argv;
#endif
  
     numarg=*argc;
     carglen=*farglen;
    printf("par init numargs: %d %s %d %d\n",numarg,fargv,carglen,*machnum);

  for (i = 0; i < numarg; i++) {
    argvp[i]=&(fargv[i*carglen]);
    printf("par init args: %d %s %s\n",i,"argvp[i]",argvp[i]);
    }

#if defined (RAMS_MPI)
    printf("par init RAMS_MPI defined \n");
  argv=&(argvp[0]);
  MPI_Init(&numarg, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,machnum);
  MPI_Comm_size(MPI_COMM_WORLD,machsize);
#endif

  printf("par_init: %d %d \n",*machnum,*machsize);
}

/*==========================================================================*/

void par_init(machnum,machsize)
     int *machnum, *machsize;
{
#if defined (RAMS_MPI)
     int argc; char **argv;

  printf("par_init0\n");
  MPI_Init(&argc, &argv);
  printf("par_init1\n");
  MPI_Comm_rank(MPI_COMM_WORLD,machnum);
  printf("par_init2\n");
  MPI_Comm_size(MPI_COMM_WORLD,machsize);
#endif

  printf("par_init: %d %d \n",*machnum,*machsize);
}

/*==========================================================================*/

void par_enroll(machnum)
     int *machnum;
{

}

/*==========================================================================*/

void par_exit()
{

#if defined (RAMS_MPI)
  MPI_Finalize();
#endif

  printf("MP exiting \n");

}

/*==========================================================================*/

void par_pause(machnum,ibarrier)
     int *machnum, *ibarrier;
{
  int ierr=0;

#if defined (RAMS_MPI)
  ierr = MPI_Barrier( MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_pause- %d %d %d \n"
          ,*machnum,*ibarrier,ierr);
}

/*==========================================================================*/

void par_ready(nmach,machnum,ibarrier)
     int *nmach,*machnum,*ibarrier;
{
  int ierr=0;

  printf("par_ready- %d %d \n",*ibarrier,*machnum);

#if defined (RAMS_MPI)
  MPI_Barrier( MPI_COMM_WORLD);
#endif

  if(ierr < 0) printf("Error in par_pause- %d %d %d \n"
          ,*machnum,*ibarrier,ierr);

}

