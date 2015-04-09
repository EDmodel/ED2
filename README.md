More information about ED can be found in the paper written by
[moorecroft et al.](http://flux.aos.wisc.edu/~adesai/documents/macrosys_papers-ankur/modeling/Moorcroft-EcolMono-EDmodel.pdf)

r956 aka 0c1bf644bd377bc0636c4f612b6f766f8e682599 from April 9th, is considered "somewhat stable", see report: https://github.com/EDmodel/ED2Documents/blob/master/EDTS/r956vr922rapid.pdf

This document describes the contents of the folders containing the 
ED code. You should have a directory called ED with the following 
structure: 
 
ED: This directory contains the ED source code (src) and the 
directory for compilation (build). For further instructions on how to 
compile and use the model, we strongly suggest accessing the ED 
Wiki website: https://github.com/EDmodel/ED2/wiki
 
RAPP: This directory contains the NCEP reanalysis pre-processor, that 
produces meteorological forcing in the ED-friendly format (HDF5) 
based on the NCEP/NCAR reanalysis (Kalnay et al 1996). The source 
code (src) and a build directory are included. The run directory 
contains the namelist and a shell script to help with the downloading 
process. A brief instruction can be found in the directory too.
