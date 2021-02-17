
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
The Ecosystem Demography Biosphere Model (ED2) is an integrated terrestrial biosphere model incorporating hydrology, land-surface biophysics, vegetation dynamics, and soil carbon and nitrogen biogeochemistry (<a href="https://dx.doi.org/10.5194/gmd-12-4309-2019">Longo et al. 2019</a>;<a href="https://dx.doi.org/10.1029/2008JG000812">Medvigy et al., 2009</a>).  Like its predecessor, ED (<a href="https://dx.doi.org/10.1890/0012-9615(2001)071[0557:AMFSVD]2.0.CO;2">Moorcroft et al., 2001</a>), ED2 uses a set of size- and age-structured partial differential equations that track the changing structure and composition of the plant canopy. With the ED2 model, in contrast to conventional biosphere models in which ecosystems with climatological grid cells are represented in a highly aggregated manner, the state of the aboveground ecosystem is described by the density of trees of different sizes and how this varies across horizontal space for a series of plant functional types. This more detailed description of ecosystem composition and structure enables the ED2 model to make realistic projections of both the fast-timescale exchanges of carbon, water and energy between the land and atmosphere, and long-term vegetation dynamics incorporating effects of ecosystem heterogeneity, including disturbance history and recovery.  

## <a name="contents"> Repository Contents </a> 
Copies of the ED2 repository should contain the following directories:
 - <b> ED: </b> Contains the ED source code (src) and the directory for compilation (build). For further instructions on how to compile and use the model, we strongly suggest accessing the ED Wiki website: https://github.com/EDmodel/ED2/wiki
 - <b> EDR: </b> Contains the source code (src), build (build), and basic run files (run) for a stripped-down version of the ED2 models radiative transfer scheme.
 - <b> EDTS: </b> Contains the ED model test suite for evaluating the results of changes to the source code under a variety of run conditions.
 - <b> BRAMS: </b> Contains a version of the Brazilian Developments on the Regional Atmospheric Model System (<a href="https://dx.doi.org/10.5194/gmd-10-189-2017">Freitas et al. 2017</a>) that was modified to run coupled biosphere-atmosphere simulations (<a href="https://dx.doi.org/10.5194/hess-19-241-2015">Knox et al. 2015</a>; <a href="https://dx.doi.org/10.1016/j.agrformet.2015.07.006">Swann et al. 2015</a>).  <i>Note</i>: This has not been tested in a while, so the code may need updates to work with the most recent version of ED2.
 - <b> Doc: </b> Contains additional ED2 documentation automatically generated with Doxygen.
 - <b> Ramspost: </b> The BRAMS post-processing program, which generates GrADS files. <i>Note</i>: This has not been tested in a while, so the code may need updates to work with the most recent version of ED2.
 - <b> RAPP: </b> This directory contains the NCEP reanalysis pre-processor, that produces meteorological forcing in the ED-friendly format (HDF5) based on the NCEP/NCAR reanalysis (Kalnay et al 1996). The source code (src) and a build directory are included. The run directory contains the namelist and a shell script to help with the downloading process. A brief instruction can be found in the directory too.<i>Note</i>: This has not been tested in a while, and other users have developed scripts to convert more up-to-date reanalyses.
 - <b> R-utils: </b> A collection of R scripts utilities for model pre- and post-processing (mostly called by R scripts located in ED/Template).

## <a name="implementation"> Implementation Notes </a> 
The primary data structure in ED, which can be found in ed_state_vars.F90, is a named, nested array of arrays. Each level of the heirarchy contains many fields of depth one, but the key large scale structure is as follows:
 - <b> grid: </b> The most coarse data in the model. Basically just a simulation book-keeping linking of polygons. 
 - <b> polygon: </b> A collection of sites sharing a meteorology.
 - <b> site: </b> A collection of patches sharing a common soil system and ground hydrology.
 - <b> patch: </b> A collection of cohorts sharing a disturbance history and age.
 - <b> cohort: </b> A collection of plants of identical PFT and height.
 
Note: height and age, being continuous variables, are "binned". "Identical" in the above refers to bin membership.

## <a name="info"> Further Information </a> 

- Most of the existing documentation on how to pre-process, compile, and run the ED2 model is available in our <a href="https://github.com/EDmodel/ED2/wiki">Wiki</a>.
- For the technical description of the various ED2 model features, we suggest looking at this <a href="https://github.com/EDmodel/ED2/wiki/References-for-technical-description-of-ED-2.2">reference list</a> (<b>Tip</b>: do not skip the Supporting Information of these papers, they contain relevant details).
- For a partial list of studies that have used ED or ED2, check <a href="https://github.com/EDmodel/ED2/wiki/Publications">here</a>.

##<a name="stable">Stable version</a>

- The latest stable version of the code (ED-2.2) can be found <a href="https://github.com/EDmodel/ED2/releases">here</a>.  
- For former versions of ED or ED2, check the <a href="http://www.oeb.harvard.edu/faculty/moorcroft/code_and_data/index.html">Moorcroft Lab website</a>.

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
