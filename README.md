TEST
<!-- NOTE: This page uses a combination of HTML and markdown to display properly in both GitHub and Doxygen. -->

## Important: Please see <a href="#doxygit"> Doxygen and Git Commits </a>

## Table of Contents
1. <a href="#overview"> Model Overview </a>
2. <a href="#contents"> Repository Contents </a> 
4. <a href="#implementation"> Implementation Notes </a> 
5. <a href="#info"> Further Info </a> 
   1. <a href="#doxygit"> Doxygen and Git Commits </a>
   1. <a href="#doxyinfo"> General Doxygen Info </a>
   
## <a name="overview"> Model Overview </a>
The Ecosystem Demography Biosphere Model (ED2) is an integrated terrestrial biosphere model incorporating hydrology, land-surface biophysics, vegetation dynamics, and soil carbon and nitrogen biogeochemistry (Medvigy et al., 2009). Like its predecessor, ED (Hurtt et al., 1998, Moorcroft et al., 2001), ED2 utilizes a set of size- and age-structured partial differential equations that track the changing structure and composition of the plant canopy. With the ED2 model, in contrast to conventional biosphere models in which ecosystems with climatological grid cells are represented in a highly aggregated manner, the state of the aboveground ecosystem is described by the density of trees of different sizes and how this varies across horizontal space for a series of plant functional types. This more detailed description of ecosystem composition and structure enables the ED2 model to make realistic projections of both the fast-timescale exchanges of carbon, water and energy between the land and atmosphere, and long-term vegetation dynamics incorporating effects of ecosystem heterogeneity, including disturbance history and recovery (Hurtt et al., 2012).

## <a name="contents"> Repository Contents </a> 
Copies of the ED2 repository should contain the following directories:
 - <b> BRAMS: </b> Contains the Brazilian Regional Atmospheric Model Somethingorother.
 - <b> Doc: </b> Contains the ED2 documentation.
 - <b> ED: </b> Contains the ED source code (src) and the directory for compilation (build). For further instructions on how to compile and use the model, we strongly suggest accessing the ED Wiki website: https://github.com/EDmodel/ED2/wiki
 - <b> EDR: </b> Contains the source code (src), build (build), and basic run files (run) for a stripped-down version of the ED2 models radiative transfer scheme.
 - <b> EDTS: </b> Contains the ED model test suite for evaluating the results of changes to the source code under a variety of run conditions.
 - <b> Ramspost: </b> The Regional Atmospheric Model's Post Processor 
 - <b> RAPP: </b> This directory contains the NCEP reanalysis pre-processor, that produces meteorological forcing in the ED-friendly format (HDF5) based on the NCEP/NCAR reanalysis (Kalnay et al 1996). The source code (src) and a build directory are included. The run directory contains the namelist and a shell script to help with the downloading process. A brief instruction can be found in the directory too.
 - <b> R-utils: </b> A collection of utilities for model pre- and post-processing. 

## <a name="implementation"> Implementation Notes </a> 
The primary data structure in ED, which can be found in ed_state_vars.f90, is a named, nested array of arrays. Each level of the heirarchy contains many fields of depth one, but the key large scale structure is as follows:
 - <b> grid: </b> The most coarse data in the model. Basically just a simulation book-keeping linking of polygons. 
 - <b> polygon: </b> A collection of sites sharing a meteorology.
 - <b> site: </b> A collection of patches sharing a common soil system and ground hydrology.
 - <b> patch: </b> A collection of cohorts sharing a disturbance history and age.
 - <b> cohort: </b> A collection of plants of identical PFT and height.
 
Note: height and age, being continuous variables, are "binned". "Identical" in the above refers to bin membership.

If you're not sure where to start in browsing the documentation, consider looking at ed_model.f90, which controls the actual simulation of ecosystem processes. ed_driver.f90 and edmain.f90 mostly do model initialization and coordination of things like mpi.

## <a name="info"> Further Info </a> 
This documentation includes clickable callgraphs and caller-graphs for each function in the code, except the routine "fatal_error" which is connected to just about everything.

More information about ED can be found in the paper written by
[Moorcroft et al.](http://flux.aos.wisc.edu/~adesai/documents/macrosys_papers-ankur/modeling/Moorcroft-EcolMono-EDmodel.pdf)

r956 aka 0c1bf644bd377bc0636c4f612b6f766f8e682599 from April 9th, is considered "somewhat stable", see report: https://github.com/EDmodel/ED2Documents/blob/master/EDTS/r956vr922rapid.pdf


### <a name="doxygit"> Doxygen and Git Commits: </a>
In order for further pull requests to the mainline to be accepted, modified subroutines will require the following doxygen tags:
 - Brief subroutine descriptions using the brief tag
 - Detailed subroutine descriptions using the details tag
 - Subroutine authorship statements using author tag
 - Inline subroutine argument descriptions using "!< comment"

Please do not commit changes to model code and documentation together if high numbers of documentation files have been modified. Instead ...
 - Submit a pull request to the mainline with a comment that docs needs regeneration or
 - Seperately pull request code changes and documentation updates

Finally, please note:
 - Seperating code and doc changes will make inter-branch comparison much easier when many files are modified. 
 - Regenerating documentation may produce many spurious line-end encoding differences which git will pick up on. "git add"-ing such files will often return them to the repo standard CRLF, and they will cease to be listed as modified.

### <a name="doxyinfo"> General Doxygen Info </a>
The following info may be helpful for familiarizing one's self with Doxygen, an auto-documentation program which utilizes a system of tags in source code. To tag a subroutine in the ED model as required above, try taking ed_model.f90 as a template and/or browsing the first few links below.
 - NASA Doxygen Quickstart: https://modelingguru.nasa.gov/docs/DOC-1811
 - Doxygen Keywords: https://www.stack.nl/~dimitri/doxygen/manual/commands.html
 - Doxygen Documentation: https://www.stack.nl/~dimitri/doxygen/manual/index.html
 
More General Info
 - http://www.msg.chem.iastate.edu/gamess/DoxygenRules.oct10.pdf
 - http://stackoverflow.com/questions/51667/best-tips-for-documenting-code-using-doxygen
 - https://www.rosettacommons.org/manuals/rosetta3.2_user_guide/doxygen_tips.html
 
Mainpage & custom file construction links:
 - http://stackoverflow.com/questions/10136201/mainpage-in-doxygen-documentation?rq=1
 - http://stackoverflow.com/questions/9502426/how-to-make-an-introduction-page-with-doxygen
 - http://stackoverflow.com/questions/13368350/use-the-readme-md-file-as-main-page-in-doxygen
 - http://stackoverflow.com/questions/3052036/how-to-include-custom-files-in-doxygen/6336368#6336368
 - http://stackoverflow.com/questions/2337307/doxygen-adding-a-custom-link-under-the-related-pages-section
 
Using grouping:
 - http://www.stack.nl/~dimitri/doxygen/manual/grouping.html

Configuration: 
 - http://www.doxygen.nl/config.html

Documentation generated by doxygen can be accessed (locally) on a machine with this repository under Doc/html/index.html
