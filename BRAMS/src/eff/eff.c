#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#if defined(sun) || defined(__svr4__)
#  include <string.h>
#  include <dirent.h>
#  define  DIRENT dirent
#else
#  define _BSD
#  include <sys/dir.h>
#  undef _BSD
#  include <strings.h>
#  define  DIRENT direct
#endif
#include <memory.h>

#ifndef PROGNAME
#  define PROGNAME "eff"
#endif
#ifndef ENVNAME
#  define ENVNAME  "RAMS_ROOT"
#endif
#ifndef DEFDIRS
#  define DEFDIRS  {"bin","scripts","etc",NULL}
#endif


#define MAXLEN 128
char *eff(char *, char *);
void eff_(char *, char *,char *,int,int,int);

#ifdef PROGRAM
int main(argc,argv)
  int argc; char *argv[];
{
  extern char *egetenv(char *);
  char *result;
  
  if (argc < 2  || argc > 3)
  {
    fprintf(stderr,"Usage: %s file_to_find [project_directory]\n",PROGNAME);
    fprintf(stderr,"  Searches for file in the local directory first\n");
    fprintf(stderr,"  and then the bin, scripts and etc sub-directories of\n");
    fprintf(stderr,"  of the project and RAMS_ROOT directories\n");
    exit(1);
  }
  if (argc == 2) {
    result=eff(argv[1],NULL);
  } else {
    result=eff(argv[1],argv[2]);
  }
  if (result!=NULL) printf("%s\n",result);
  exit(0);
}
#endif

char *eff(filename,projdir)
  char *filename;   /* file to find */
  char *projdir;    /* optional project directory to search */
{
  static char *defdirs[]=DEFDIRS;
  char  dirs[20][MAXLEN];
  static char path[100], *env;
  char rams_root[MAXLEN];
  struct stat filestat;
  int idx, ldx, found;
  int haveramsroot, haveproject;
  DIR   *dir;
  struct DIRENT *direntry;
  char *ptr;
  
  /* First check for the RAMS_ROOT environmental variable */
  haveramsroot=1;
  if ( (env=getenv(ENVNAME)) != NULL)
  {
    ptr = env+strlen(env)-1;
    if (*ptr == '/') { *ptr = '\0';}  /* Remove ending dir. slash if there */
    strncpy(rams_root,env,(size_t)MAXLEN-1);
  } else {
    if ( (env=getenv("HOME")) != NULL)
    {
      strncpy(rams_root,env,(size_t)MAXLEN-1);
      strncat(rams_root,"/rams",(size_t)MAXLEN-1);
    } else {
      haveramsroot=0;
    }    

  }

  /* Next test the directories */
  haveproject=0;
  if (projdir != NULL) {
    ptr = projdir+strlen(projdir)-1;
    if (*ptr == '/') { *ptr = '\0';}  /* Remove ending dir. slash if there */
    if (stat(projdir,&filestat) == 0) {
      if (! S_ISDIR(filestat.st_mode)) {
         fprintf(stderr,"%s: Project directory, %s, is not a directory\n",PROGNAME,projdir);
         return(NULL);
       } else {
         haveproject=1;
       }
    } else {
       fprintf(stderr,"%s: Cannot read Project directory, %s\n",PROGNAME,projdir);
       return(NULL);
    }
  }

  if (haveramsroot) {
    if (stat(rams_root,&filestat) == 0) {
      if (! S_ISDIR(filestat.st_mode)) {
        haveramsroot=0;
      }
    } else {
      haveramsroot=0;
    }
  }

  if (!haveproject && !haveramsroot) {
    fprintf(stderr,"%s: No RAMS_ROOT or project directory found.\n",PROGNAME);
    return(NULL);
  }

  /* Next construct a list of directories to search */

     /* First the current directory for a development copy */
  idx=0;
  if ( (env=getenv("PWD")) != NULL)
  {
    strncpy(dirs[idx],env,(size_t)MAXLEN-1);
  } else {
    strncpy(dirs[idx],".",(size_t)MAXLEN-1);
  }
  idx++;
  
     /* Next the ~/bin directory */
  idx=0;
  if ( (env=getenv("HOME")) != NULL)
  {
    strncat(env,"/bin",(size_t)MAXLEN-1);
    strncpy(dirs[idx],env,(size_t)MAXLEN-1);
  }
  idx++;
  
     /* Next the local project directory */
     /* Next the local project directory */
  if (haveproject) {
    ldx=0;
    while (defdirs[ldx] != NULL) {
      strncpy(dirs[idx],projdir,(size_t)MAXLEN-1);
      strncat(dirs[idx],"/",(size_t)MAXLEN-1);
      strncat(dirs[idx],defdirs[ldx++],(size_t)MAXLEN-1);
      idx++;
    }
  }
     /* Finally the RAMS_ROOT directory */
  if (haveramsroot) {
    ldx=0;
    while (defdirs[ldx] != NULL) {
      strncpy(dirs[idx],rams_root,(size_t)MAXLEN-1);
      strncat(dirs[idx],"/",(size_t)MAXLEN-1);
      strncat(dirs[idx],defdirs[ldx++],(size_t)MAXLEN-1);
      idx++;
    }
  }
  dirs[idx][0] = 0;

  path[0]='\0';
  found=0;
  /* for (idx=0;dirs[idx][0] != 0;idx++) {printf("dirs[idx]=%s\n",dirs[idx]); }*/
  idx=0;
  while (found==0 && dirs[idx][0] != 0) {
    if ((dir = opendir(dirs[idx])) != NULL) {
      while (found==0 && (direntry = readdir(dir))) {
        /* printf("%s/%s\n",dirs[idx],direntry->d_name);*/
	if (strcmp(filename,direntry->d_name) == 0) {
          strcpy(path,dirs[idx]);
	  strcat(path,"/");
	  strcat(path,direntry->d_name);
          // Ignore if a directory
          if (stat(path,&filestat) != 0) continue;
          if (S_ISDIR(filestat.st_mode)) continue;
	  found=1;
	}
      }
      closedir(dir);
    }

    idx++;
  }

  if (strlen(path) > 0) return(path);
  else return(NULL);
}

/* FORTRAN interface */
void eff_(filename,projdir,result,siz1,siz2,siz3)
  char *filename;   /* file to find */
  char *projdir;    /* optional project directory to search */
  char *result;     /* String for result */
{

  char fn[128]="";
  char pj[128]="";
  char *res, *ptr;

  strncpy(fn,filename,128);
  ptr = index(fn,' ');
  *ptr = '\0';

  if (siz2 == 0) {
    res=eff(fn,NULL);
  } else {
    strncpy(pj,projdir, 128);
    ptr = index(pj,' ');
    *ptr = '\0';
    res=eff(fn,pj);
  }

  if (res != NULL) {
     strncpy(result, res,siz3);
  } else {
     strncpy(result, " ",siz3);

  }
}
