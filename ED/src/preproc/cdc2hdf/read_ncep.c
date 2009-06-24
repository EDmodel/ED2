/*  SUMMARY: 
      1. This programs reads the global NCEP reanalysis files (in netcdf 
      format) that you have downloaded.  
      2. It then prints out the data for a subregion over a specified time 
      period in ASCII. 
*/

/* VARIABLES OF INTEREST
   I hope these are the only things you will need to modify.  However, 
   search for GOTCHA below.
    LATMIN  --- minimum latitude of the subregion
    LATMAX  --- maximum latitude of the subregion
    LONMIN  --- minimum longitude of the subregion
    LONMAX  --- maximum longitude of the subregion
    MIN_YEAR --- first year you want to process
    MAX_YEAR --- last year you want to process
    INPUT_DIR --- location of your downloaded ncep input files
    OUTPUT_DIR --- directory where you want to put your output files 
*/



#define  _GNU_SOURCE		/* for NaN */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <netcdf.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

/* USER CONTROL */
/* ************ */
#define LATMIN  -60
#define LATMAX   30
#define LONMIN -120
#define LONMAX  -25
#define MIN_YEAR 2008
#define MAX_YEAR 2008
#define INPUT_DIR  "/n/data/moorcroft_lab/mlongo/NCEP_MET/netcdf"
#define OUTPUT_DIR "/n/data/moorcroft_lab/mlongo/NCEP_MET/ascii"
#define OUTPUT_PREF "SOUTHAM"
/* ************ */

#define NR_END 1
#define FREE_ARG char*
#define TIMESTEP 6
int IMAX(int a,int b);
int IMIN(int a,int b);
#define check_ncstat(stat, msg) _check_ncstat(stat, __FILE__, __LINE__, msg);

float rslf(float p, float t);
void _check_ncstat(int stat, const char *filename, unsigned line,
		   const char *msg);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
/**********************************************************/

struct ncvar {
  char *filename;	/* file data came from */
  int ncid;			/* file id */
  int varid;  			/* variable id */
  float minlat, maxlat;		/* latitude range */
  float latstep;
  float minlon, maxlon;		/* longitude range */
  size_t nlat, nlon;
  float lonstep;
  long int mintime, maxtime;
  size_t N;			/* number of timesteps */
  float offset, scale;		/* for decoding the variable */
  short nodata;			/* code for missing data */
};
struct ncvar * new_ncvar(const char *filename, const char *varname);
void _check_ncstat(int stat, const char *filename, unsigned line,const char *msg );
float *** ncgetvar(const struct ncvar *v);
void free_ncvar(struct ncvar *var);
void ncnearest(struct ncvar *v, float lat, float lon,
	       size_t *lati, size_t *loni);
void index2latlon(struct ncvar *v, size_t lati, size_t loni,
		  float *lat, float *lon);
void free_f3tensor(float ***t, const struct ncvar *v);

/****************************************************************/
int main(){
  int itime;
  FILE*outfile;
  char outname[256];
  FILE*routfile;
  char routname[256];
  struct ncvar *air_var = NULL;
  struct ncvar *pres_var = NULL;
  struct ncvar *rhum_var = NULL;
  struct ncvar *uwnd_var = NULL;
  struct ncvar *vwnd_var = NULL;
  struct ncvar *dlwrf_var = NULL;
  struct ncvar *nbdsf_var = NULL;
  struct ncvar *nddsf_var = NULL;
  struct ncvar *vbdsf_var = NULL;
  struct ncvar *vddsf_var = NULL;
  struct ncvar *prate_var = NULL;
  struct ncvar *tmpvar = NULL;
  float ***air = NULL;
  float ***pres = NULL;
  float ***rhum = NULL;
  float ***uwnd = NULL;
  float ***vwnd = NULL;
  float ***dlwrf = NULL;
  float ***nbdsf = NULL;
  float ***nddsf = NULL;
  float ***vbdsf = NULL;
  float ***vddsf = NULL;
  float ***prate = NULL;
  char infilename[256];
  int lati,loni;
  size_t latmini,lonmini,latmaxi,lonmaxi,rlati,rloni;
  float lat,lon;
  double time;
  unsigned int timei;
  float u;
  int beg_time[13];
  int imonth;
  int nlat,nlon;
  float lat_min,lon_min,lat_max,lon_max,rv;
  int idate;
  float stime;
  int offset_from_utc,id,ib,ih;
  int ntime,iyear;
  char tmpvar_name[256];
  char mname[3];
  char syscmd[256];
  float var_out[1464];
  int ndims;
  int idims[1];
  int i;
  char var_name[20];
  int hdferr;
  int timeo;

  /* PART I:  Find information about your grid */

  /* Read an arbitrary netcdf variable in order to find information about 
     the grid.  */

  /* GOTCHA:  this presumes that you have the file containing air temperature
     data.  If you do not, please change the next two lines to refer to a 
     variable that you do have.  CANNOT BE A VARIABLE WITH GAUSS IN THE 
     FILE NAME.  THAT IS WHY WE USE UWND HERE.  */

  sprintf(tmpvar_name,"%s/uwnd.sig995/uwnd.sig995.1948.nc",INPUT_DIR);
  tmpvar = new_ncvar(tmpvar_name,"uwnd");

  /* Find the grid cell indices nearest your specified LATMIN, LATMAX, LONMIN, 
     LONMAX */ 
  ncnearest(tmpvar,LATMIN,LONMIN,&latmini,&lonmini);
  ncnearest(tmpvar,LATMAX,LONMAX,&latmaxi,&lonmaxi);
  //  printf("%d %d\n",latmini,latmaxi);
  //printf("%d %d\n",lonmini,lonmaxi);

  /* Count the total number of grid cells */
  nlat = latmini - latmaxi + 1;
  nlon = lonmaxi - lonmini + 1;

  /* Find the coordinates corresponding the grid cell indices */
  index2latlon(tmpvar, latmini, lonmini, &lat_min, &lon_min);
  index2latlon(tmpvar, latmaxi, lonmaxi, &lat_max, &lon_max);

  /* PART II:  Read files and write the output */

  /* Make header file */
  sprintf(outname,"%s/%s_HEADER",OUTPUT_DIR,OUTPUT_PREF);
  outfile = fopen(outname,"w");
  fprintf(outfile,"%d %d\n",nlat,nlon);
  fprintf(outfile,"%f %f %f %f\n",sqrt(pow(tmpvar->latstep,2)),
	  sqrt(pow(tmpvar->lonstep,2)),lat_min,lon_min);
  fclose(outfile);

  /* Loop over years */
  for(iyear=MIN_YEAR;iyear<=MAX_YEAR;iyear++){

    /* Read in each variable */
    printf("trying air\n");
    sprintf(infilename,"%s/air.sig995/air.sig995.%d.nc",INPUT_DIR,iyear);
    air_var = new_ncvar(infilename,"air");
    air = ncgetvar(air_var);

    printf("trying pres\n");
    sprintf(infilename,"%s/pres.sfc/pres.sfc.%d.nc",INPUT_DIR,iyear);
    pres_var = new_ncvar(infilename,"pres");
    pres = ncgetvar(pres_var);
    
    printf("trying rhum\n");
    sprintf(infilename,"%s/rhum.sig995/rhum.sig995.%d.nc",INPUT_DIR,iyear);
    rhum_var = new_ncvar(infilename,"rhum");
    rhum = ncgetvar(rhum_var);
    
    printf("trying uwnd\n");
    sprintf(infilename,"%s/uwnd.sig995/uwnd.sig995.%d.nc",INPUT_DIR,iyear);
    uwnd_var = new_ncvar(infilename,"uwnd");
    uwnd = ncgetvar(uwnd_var);
    
    printf("trying vwnd\n");
    sprintf(infilename,"%s/vwnd.sig995/vwnd.sig995.%d.nc",INPUT_DIR,iyear);
    vwnd_var = new_ncvar(infilename,"vwnd");
    vwnd = ncgetvar(vwnd_var);
    
    printf("trying dlwrf\n");
    sprintf(infilename,"%s/dlwrf.sfc.gauss/dlwrf.sfc.gauss.%d.nc",INPUT_DIR,iyear);
    dlwrf_var = new_ncvar(infilename,"dlwrf");
    dlwrf = ncgetvar(dlwrf_var);
    
    printf("trying nbdsf\n");
    sprintf(infilename,"%s/nbdsf.sfc.gauss/nbdsf.sfc.gauss.%d.nc",INPUT_DIR,iyear);
    nbdsf_var = new_ncvar(infilename,"nbdsf");
    nbdsf = ncgetvar(nbdsf_var);

    printf("trying nddsf\n");
    sprintf(infilename,"%s/nddsf.sfc.gauss/nddsf.sfc.gauss.%d.nc",INPUT_DIR,iyear);
    nddsf_var = new_ncvar(infilename,"nddsf");
    nddsf = ncgetvar(nddsf_var);

    printf("trying vbdsf\n");
    sprintf(infilename,"%s/vbdsf.sfc.gauss/vbdsf.sfc.gauss.%d.nc",INPUT_DIR,iyear);
    vbdsf_var = new_ncvar(infilename,"vbdsf");
    vbdsf = ncgetvar(vbdsf_var);

    printf("trying vddsf\n");
    sprintf(infilename,"%s/vddsf.sfc.gauss/vddsf.sfc.gauss.%d.nc",INPUT_DIR,iyear);
    vddsf_var = new_ncvar(infilename,"vddsf");
    vddsf = ncgetvar(vddsf_var);

    printf("trying prate\n");
    sprintf(infilename,"%s/prate.sfc.gauss/prate.sfc.gauss.%d.nc",INPUT_DIR,iyear);
    prate_var = new_ncvar(infilename,"prate");
    prate = ncgetvar(prate_var);

    /*   The variable beg_time is used to determine what month you are in.
     */
    beg_time[0] = 0;
    beg_time[1] = 31*4;
    if(iyear%400 == 0 | (iyear%4 == 0 & iyear%100 != 0)){
      beg_time[2] = beg_time[1] + 29 * 4;
    }else{
      beg_time[2] = beg_time[1] + 28 * 4;
    }
    beg_time[3] = beg_time[2] + 31 * 4;
    beg_time[4] = beg_time[3] + 30 * 4;
    beg_time[5] = beg_time[4] + 31 * 4;
    beg_time[6] = beg_time[5] + 30 * 4;
    beg_time[7] = beg_time[6] + 31 * 4;
    beg_time[8] = beg_time[7] + 31 * 4;
    beg_time[9] = beg_time[8] + 30 * 4;
    beg_time[10] = beg_time[9] + 31 * 4;
    beg_time[11] = beg_time[10] + 30 * 4;
    beg_time[12] = beg_time[11] + 31 * 4;
    
    /* Loop over months */
    for(imonth=0;imonth<12;imonth++){
      if (imonth == 0){
         sprintf(mname,"JAN");
      }else if(imonth == 1){
         sprintf(mname,"FEB");
      }else if(imonth == 2){
         sprintf(mname,"MAR");
      }else if(imonth == 3){
         sprintf(mname,"APR");
      }else if(imonth == 4){
         sprintf(mname,"MAY");
      }else if(imonth == 5){
         sprintf(mname,"JUN");
      }else if(imonth == 6){
         sprintf(mname,"JUL");
      }else if(imonth == 7){
         sprintf(mname,"AUG");
      }else if(imonth == 8){
         sprintf(mname,"SEP");
      }else if(imonth == 9){
         sprintf(mname,"OCT");
      }else if(imonth == 10){
         sprintf(mname,"NOV");
      }else if(imonth == 11){
         sprintf(mname,"DEC");
      }
      
      /* Make output file name */
      sprintf(outname,"%s/%s_%.4d%s.dat",OUTPUT_DIR,OUTPUT_PREF,iyear,mname);
      printf("Trying new time: %s\n",outname);
      outfile = fopen(outname,"w");

      /* loop over sites */
      // Goes from north pole to south pole
      for(lati = latmaxi; lati <= latmini; lati++) { /* index decreases
							with lat */
	//	printf("lati: %d \n",lati);
	// Goes from east to west
	for(loni = lonmini; loni <= lonmaxi; loni++) {

	  index2latlon(tmpvar, lati, loni, &lat, &lon);
	  ncnearest(prate_var, lat ,lon, &rlati, &rloni);
	  
	  // Fill the output array
	  for(timei=beg_time[imonth];timei<beg_time[imonth+1];timei++){
	    fprintf(outfile,"%f %f %f %f %f %f %f %f %f %f %f\n", 
		    prate[timei][rlati][rloni],
		    dlwrf[timei][rlati][rloni], nbdsf[timei][rlati][rloni],
		    nddsf[timei][rlati][rloni], vbdsf[timei][rlati][rloni],
		    vddsf[timei][rlati][rloni], air[timei][rlati][rloni],
		    pres[timei][rlati][rloni], rhum[timei][rlati][rloni],
		    uwnd[timei][rlati][rloni], vwnd[timei][rlati][rloni]);
	  }
	}
      }

    }
    
    /* Free memory */

    free_f3tensor(air, air_var);
    free_f3tensor(pres, pres_var);
    free_f3tensor(rhum, rhum_var);
    free_f3tensor(uwnd, uwnd_var);
    free_f3tensor(vwnd, vwnd_var);
    free_f3tensor(dlwrf, dlwrf_var);
    free_f3tensor(nbdsf, nbdsf_var);
    free_f3tensor(nddsf, nddsf_var);
    free_f3tensor(vbdsf, vbdsf_var);
    free_f3tensor(vddsf, vddsf_var);
    free_f3tensor(prate, prate_var);

    free_ncvar(air_var);
    free_ncvar(pres_var);
    free_ncvar(rhum_var);
    free_ncvar(uwnd_var);
    free_ncvar(vwnd_var);
    free_ncvar(dlwrf_var);
    free_ncvar(nbdsf_var);
    free_ncvar(nddsf_var);
    free_ncvar(vbdsf_var);
    free_ncvar(vddsf_var);
    free_ncvar(prate_var);

  }

  free_ncvar(tmpvar);

  return 0;
}

/*******************************************************************/


struct ncvar * new_ncvar(const char *filename, const char *varname) {
  struct ncvar *v = NULL;
  int varid;			/* variable id */
  int dimid;			/* dimension id */
  float range[2];		/* min and max */
  int ncstat;
  
  v = calloc(1, sizeof(struct ncvar));
  
  /* open file */
  ncstat = nc_open(filename, NC_NOWRITE, &v->ncid);
  check_ncstat(ncstat, filename);

  /* get lat/lon/time sizes */
  ncstat = nc_inq_dimid(v->ncid, "lat", &dimid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_inq_dimlen(v->ncid, dimid, &v->nlat); 
  check_ncstat(ncstat, NULL);

  ncstat = nc_inq_dimid(v->ncid, "lon", &dimid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_inq_dimlen(v->ncid, dimid, &v->nlon); 
  check_ncstat(ncstat, NULL);

  ncstat = nc_inq_dimid(v->ncid, "time", &dimid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_inq_dimlen(v->ncid, dimid, &v->N); 
  check_ncstat(ncstat, NULL);

  /* get mins and maxes, calculate stepsizes */
  ncstat = nc_inq_varid (v->ncid, "lat", &varid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, varid, "actual_range", range);
  check_ncstat(ncstat, NULL);
  v->minlat = range[0]; v->maxlat = range[1];
  v->latstep = (v->maxlat - v->minlat) / (v->nlat -1.0);
  
  ncstat = nc_inq_varid (v->ncid, "lon", &varid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, varid, "actual_range", range);
  check_ncstat(ncstat, NULL);
  v->minlon = range[0]; v->maxlon = range[1];
  v->lonstep = (v->maxlon - v->minlon) / (v->nlon -1.0);
  
  ncstat = nc_inq_varid (v->ncid, "time", &varid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, varid, "actual_range", range);
  check_ncstat(ncstat, NULL);
  v->mintime = range[0]; v->maxtime = range[1];

  /* get varid */
  ncstat = nc_inq_varid (v->ncid, varname, &v->varid);
  check_ncstat(ncstat, NULL);

  /* get var unpacking info */
  ncstat = nc_get_att_float(v->ncid, v->varid, "scale_factor", &v->scale);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, v->varid, "add_offset", &v->offset); 
  check_ncstat(ncstat, NULL);

  /* get nodata val */
  ncstat = nc_get_att_short(v->ncid, v->varid, "missing_value", &v->nodata); 
  check_ncstat(ncstat, NULL);

  /* record filename for debugging */
  v->filename = strdup(filename);
  
  return v;
} /* new_ncvar */

void _check_ncstat(int stat, const char *filename, unsigned line,
		   const char *msg )
{
  if(stat == NC_NOERR) return; /* no error */
  fprintf(stderr, "%s:%d netCDF error: %s",
	  filename, line, nc_strerror(stat));
  if(msg != NULL) {
    fprintf(stderr, " %s", msg);
  }
  putchar('\n');
  exit(1);
}

float *** ncgetvar(const struct ncvar *v) {
  int ncstat;
  int i;
  size_t len;
  float ***arr = NULL;
  float *ptr = NULL;

  arr = f3tensor(0,  v->N-1, 0, v->nlat-1, 0, v->nlon-1);
  ptr = **arr;
  len = v->N * v->nlat * v->nlon;

  ncstat = nc_get_var_float(v->ncid, v->varid, ptr);
  check_ncstat(ncstat, NULL);
  
  /* scale the values */
  for(i = 0; i<len; i++) {
    if (ptr[i] == v->nodata) {	/* handle nodata */
      ptr[i] = NAN;
    } else {
      ptr[i] *= v->scale;
      ptr[i] += v->offset;
    }
  }
  return arr;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh){
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  t+= NR_END;
  t -= nrl;

  t[nrl] = (float **)malloc ((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  t[nrl] += NR_END;
  t[nrl] -= ncl; 

  t[nrl][ncl] = (float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  for(j=ncl+1;j<=nch;j++)t[nrl][j] = t[nrl][j-1] + ndep;
  for(i=nrl+1;i<=nrh;i++){
    t[i] = t[i-1]+ncol;
    t[i][ncl] = t[i-1][ncl] + ncol * ndep;
    for(j=ncl+1;j<=nch;j++)t[i][j] = t[i][j-1] + ndep;
  }
  return t;
}

void free_ncvar(struct ncvar *var) {
  if(var == NULL) return;
  nc_close(var->ncid);
  free(var->filename);
  free(var);
  var = NULL;
} /* free_ncvar */

void ncnearest(struct ncvar *v, float lat, float lon,
	       size_t *lati, size_t *loni)
{
  if(0)printf("lat %.1f lon %.1f > \n", lat, lon);  
  if (lon < 0) lon += 360;			/* convert to 0-360 range */
  *lati = rint( (lat - v->minlat) / v->latstep );
  *loni = rint( (lon - v->minlon) / v->lonstep );
  if( *lati >= v->nlat) *lati = v->nlat - 1;
  if( *loni >= v->nlon) *loni = v->nlon - 1;
  if(0) {
    /* debug */
    printf("lati %ld loni %ld\n", *lati, *loni);
  }
}

void index2latlon(struct ncvar *v, size_t lati, size_t loni,
		  float *lat, float *lon)
{
  *lat = v->minlat + (v->latstep * lati);
  *lon = v->minlon + (v->lonstep * loni);
  *lon += 0.5 * v->lonstep;

  if(*lon > 180) *lon -= 360;
  /* debug */
  if(0)printf("lat %.1f lon %.1f < lati %ld loni %ld\n", *lat, *lon, lati, loni);
}

float rslf(float p, float t){

  float esl,rslf,x,c0,c1,c2,c3,c4,c5,c6,c7,c8;
  float y;

  c0= .6105851e+03;
  c1= .4440316e+02;
  c2= .1430341e+01;
  c3= .2641412e-01;
  c4= .2995057e-03;
  c5= .2031998e-05;
  c6= .6936113e-08;
  c7= .2564861e-11;
  c8= -.3704404e-13;

  x = t - 273.16;
  if(x < -80.0)x = -80.0;

  esl = c0 + x * ( c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * (c6 + x * (c7 + x * c8)))))));
  y = .622 * esl / (p - esl);

  return(y);
}

int IMIN(int a,int b){
  int min_out;
  
  min_out = a;
  if(b < a)min_out = b;

  return(min_out);
}
int IMAX(int a,int b){
  int max_out;
  
  max_out = a;
  if(b > a)max_out = b;

  return(max_out);
}
void free_f3tensor(float ***t, const struct ncvar *v){
/* free a float f3tensor allocated by f3tensor() */

  long nrl, nrh, ncl, nch, ndl, ndh;

  nrl = 0;
  nrh = v->N-1;
  ncl = 0;
  nch = v->nlat-1;
  ndl = 0;
  ndh = v->nlon-1;

  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
