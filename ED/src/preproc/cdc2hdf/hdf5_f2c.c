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

#include "utils_sub_names.h"

#include "/n/Moorcroft_Lab/Users/dmm2/HDF5/hdf5-1.6.5/hdf5/include/hdf5.h"

hid_t fileid, dsetid, dspcid, mspcid, propid;

/** File routines *********************************************************/

void fh5f_open_(char *locfn, int *iaccess, int *hdferr)
{

unsigned flags;
hid_t access_id;
extern hid_t fileid;


access_id = H5P_DEFAULT;
if(*iaccess == 1) flags = H5F_ACC_RDONLY;
if(*iaccess == 2) flags = H5F_ACC_RDWR;

fileid = H5Fopen(locfn, flags, access_id);

//printf("fh5f_open_ - fileid: %d\n",fileid);

*hdferr = fileid;

return;
}

/******/

void fh5f_create_(char *locfn, int *iaccess, int *hdferr)
{

unsigned flags;
hid_t access_id,create_id;
extern hid_t fileid;


access_id = H5P_DEFAULT;
create_id = H5P_DEFAULT;
if(*iaccess == 1) flags = H5F_ACC_TRUNC;
if(*iaccess == 2) flags = H5F_ACC_EXCL ;

fileid = H5Fcreate(locfn, flags, create_id, access_id);

//printf("fh5f_open_ - fileid: %d\n",fileid);

*hdferr = fileid;

return;
}

/******/

void fh5f_close_(int *hdferr)
{

extern hid_t fileid;
herr_t herr;

herr = H5Fclose(fileid);
//printf("fh5f_close: %d\n",herr);

*hdferr = herr;

return;
}

/** Dataset routines ****************************************************/

void fh5d_open_(char *dname, int *hdferr)
{

extern hid_t fileid;
herr_t herr;

dsetid = H5Dopen(fileid, dname);
//printf("fh5d_open: %d\n",dsetid);

if (dsetid < 0) { *hdferr = dsetid; return;}


dspcid = H5Dget_space(dsetid);
//printf("fh5d_get_space: %d\n",dspcid);

*hdferr = dspcid;

return;
}

/******/

void fh5d_close_(int *hdferr)
{

extern hid_t dsetid;
herr_t herr;

herr = H5Dclose(dsetid);
//printf("fh5d_close: %d\n",herr);

return;
}

/******/

/** Dataset info routines ****************************************************/

void fh5s_get_ndims_(int *ndims)
{

extern hid_t dspcid;
int result;

result = H5Sget_simple_extent_ndims(dspcid);
/* printf("fh5d_get_ndims: %d\n",result); */

*ndims = result;

return;
}

/******/

void fh5s_get_dims_(int *dims)
{

extern hid_t dspcid;
hsize_t dimsc[7], maxdimsc[7];
int ndims,i;

ndims = H5Sget_simple_extent_dims(dspcid, dimsc, maxdimsc);
/* printf("fh5d_get_dims: %d %d %d %d\n",ndims,dimsc[0],dimsc[1],dimsc[2]);
 */
 
for (i = 0; i < ndims; i++) {  dims[i] = dimsc[i];   }
//printf("fh5d_get_dims: %d %d %d %d\n",ndims,dimsc[0],dimsc[1],dimsc[2]);
return;
}

/** Reading routines ****************************************************/

void fh5_prepare_read_(char *dname,int *ndims,int *dims,int *hdferr)
{

extern hid_t fileid, dsetid, dspcid, mspcid;

int i;
herr_t herr;

static hsize_t start[7]  = {0,0,0,0,0,0,0};
static hsize_t stride[7] = {1,1,1,1,1,1,1};
static hsize_t block[7]  = {1,1,1,1,1,1,1};

hsize_t dimsc[7]  = {1,1,1,1,1,1,1};
hsize_t maxdims[7]  = {1,1,1,1,1,1,1};
hsize_t count[7]  = {0,0,0,0,0,0,0};

//printf("fh5_prep - ndims: %d\n",*ndims);

*hdferr = 0;

for (i = 0; i < *ndims; i++) {  dimsc[i] = dims[i]; maxdims[i] = dims[i]; 
                                count[i] = dims[i]; }

dsetid = H5Dopen(fileid, dname);
//printf("fh5_prep - open: %d\n",dsetid);
if (dsetid < 0) { *hdferr = dsetid; return;}


dspcid = H5Dget_space(dsetid);
//printf("fh5d_get_space: %d\n",dspcid);
if (dspcid < 0) { *hdferr = dspcid; return;}


/* Select hyperslab in the dataset. */
//herr = H5Sselect_hyperslab(dspcid, H5S_SELECT_SET,start,count,stride,block);
//printf("fh5s_select_hyperslab - 1: %d\n",herr);
//if (herr < 0) { *hdferr = herr; return;}

/* Create memory dataspace for the dataset. */ 
//mspcid = H5Screate_simple(*ndims, dimsc, maxdims);
//printf("fh5s_create_simple: %d\n",mspcid);
//if (mspcid < 0) { *hdferr = mspcid; return;}

/* Select hyperslab in memory */
//herr = H5Sselect_hyperslab(mspcid, H5S_SELECT_SET,start,count,stride,block);
//printf("fh5s_select_hyperslab - 2: %d\n",herr);
//if (herr < 0) { *hdferr = herr; return;}


return;
}

/******/

void fh5d_read_(int *h5type, void *buf, int *hdferr)
{

extern hid_t dsetid, dspcid, mspcid;
herr_t herr;
hid_t memtype;

if (*h5type == 1) memtype=H5T_NATIVE_INT;
if (*h5type == 2) memtype=H5T_NATIVE_FLOAT;
if (*h5type == 3) memtype=H5T_NATIVE_CHAR;
if (*h5type == 4) memtype=H5T_NATIVE_DOUBLE;
if (*h5type == 5) memtype=H5T_NATIVE_HBOOL;


herr = H5Dread(dsetid,memtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf);
//printf("fh5d_read: %d\n",herr);
//exit(0);

*hdferr = herr;

return;
}

/******/

void fh5_close_read_(int *hdferr)
{

extern hid_t dsetid, dspcid, mspcid;
herr_t herr;

herr = H5Dclose(dsetid);
//printf("fh5d_close: %d\n",herr);

*hdferr = herr;

return;
}

/** Writing routines ****************************************************/

void fh5_prepare_write_(int *ndims,int *dims,int *hdferr)
{

extern hid_t fileid, dsetid, dspcid, mspcid, propid;
int i;
herr_t herr;

hsize_t dimsc[7]  = {1,1,1,1,1,1,1};
hsize_t maxdims[7]  = {1,1,1,1,1,1,1};
hsize_t chunk_size[7]  = {0,0,0,0,0,0,0};

for (i = 0; i < *ndims; i++) {  dimsc[i] = dims[i]; chunk_size[i] = dims[i];  
                                maxdims[i] = dims[i];}

/*  Create the data space for the dataset. */
mspcid = H5Screate_simple(*ndims, dimsc, maxdims);
//printf("fh5_prepw - create 1: %d\n",mspcid);

/* Create properties for gzip compression.*/
propid = H5Pcreate(H5P_DATASET_CREATE);
//printf("fh5_prepw - propid: %d\n",propid);

herr = H5Pset_chunk  (propid, *ndims, chunk_size);
herr = H5Pset_shuffle(propid);
herr = H5Pset_deflate(propid, 9);

*hdferr = herr;
//printf("fh5_prepw - compress: %d\n",mspcid);

return;
}

/******/

void fh5_write_(int *h5type, void *buf, char *dname, int *hdferr)
{

extern hid_t fileid, dsetid, dspcid, mspcid, propid;
herr_t herr;
hid_t memtype;

if (*h5type == 1) memtype=H5T_NATIVE_INT32;
if (*h5type == 2) memtype=H5T_NATIVE_FLOAT;
if (*h5type == 3) memtype=H5T_STRING;
if (*h5type == 4) memtype=H5T_NATIVE_DOUBLE;
if (*h5type == 5) memtype=H5T_NATIVE_HBOOL;
//printf("fh5_write - start: %d \n",memtype);

dsetid = H5Dcreate(fileid, dname, memtype, mspcid, propid);
//printf("fh5_write - create 1: %d\n",dsetid);

herr = H5Dwrite(dsetid, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
//printf("fh5_write - write 1: %d\n",herr);

*hdferr = herr;

return;
}
/******/

void fh5_close_write_(int *hdferr)
{

extern hid_t dsetid, dspcid, mspcid;
herr_t herr;

herr = H5Sclose(mspcid);
herr = H5Pclose(propid);
herr = H5Dclose(dsetid);
//printf("fh5_close_write: %d\n",herr);

*hdferr = herr;

return;
}


