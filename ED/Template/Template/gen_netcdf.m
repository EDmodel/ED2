firstFile ='thispoly-S-myyeara-mymontha-01-000000-g01.h5';
lastFile ='thispoly-S-myyearz-mymonthz-mydatez-000000-g01.h5';
path1 = 'paththere/thispoly/histo/';
path2 = 'paththere/nc_files/';
SaveName=['thispoly'];
Freq='HR';
VarList='VarList_DiurnalCheck';
CO=[1 1];

TMP=[path1 firstFile];
TMP2=[TMP(1:end-13) '010000-g01.h5'];
copyfile(TMP2, TMP)

cd /n/moorcroft_scratch/nlevine/data/mfiles/ED_genNCfile
Create_ncfile_fn (firstFile, lastFile, path1, path2, SaveName, Freq, VarList, CO)
