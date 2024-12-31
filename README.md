
<!-- NOTE: This page uses a combination of HTML and markdown to display properly in both GitHub and Doxygen. -->

## Table of Contents
1. <a href="#overview"> Model Overview </a>
2. <a href="#stable">Current version and stable release versions</a>
3. <a href="#contents"> Repository Contents </a> 
4. <a href="#implementation"> Implementation Notes </a> 
5. <a href="#info"> Further Information </a>
6. <a href="#develop">Code Development, Pull Requests, and Commits</a>
7. <a href="#docker">Using Docker</a>


## <a name="overview"> Model Overview </a>
The Ecosystem Demography Biosphere Model (ED2) is an integrated terrestrial biosphere model incorporating hydrology, land-surface biophysics, vegetation dynamics, and soil carbon and nitrogen biogeochemistry (<a href="https://dx.doi.org/10.5194/gmd-12-4309-2019">Longo et al. 2019</a>;<a href="https://dx.doi.org/10.1029/2008JG000812">Medvigy et al., 2009</a>).  Like its predecessor, ED (<a href="https://dx.doi.org/10.1890/0012-9615(2001)071[0557:AMFSVD]2.0.CO;2">Moorcroft et al., 2001</a>), ED2 uses a set of size- and age-structured partial differential equations that track the changing structure and composition of the plant canopy. With the ED2 model, in contrast to conventional biosphere models in which ecosystems with climatological grid cells are represented in a highly aggregated manner, the state of the aboveground ecosystem is described by the density of trees of different sizes and how this varies across horizontal space for a series of plant functional types. This more detailed description of ecosystem composition and structure enables the ED2 model to make realistic projections of both the fast-timescale exchanges of carbon, water and energy between the land and atmosphere, and long-term vegetation dynamics incorporating effects of ecosystem heterogeneity, including disturbance history and recovery.  

## <a name="stable">Current version and stable release versions</a>

- This code available on GitHub is the current version of the model, updated frequently.
- The latest stable version of the code (ED-2.2) can be found <a href="https://github.com/EDmodel/ED2/releases">here</a>.  
- For former versions of ED or ED2, check the <a href="http://www.oeb.harvard.edu/faculty/moorcroft/code_and_data/index.html">Moorcroft Lab website</a>.
- A summary of the main changes between stable releases is available <a href="https://github.com/EDmodel/ED2/wiki/ED2-release-notes">here</a>.

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

Note: height and age, being continuous variables, are "binned". "Identical" in this context means sufficiently similar to be placed in the same bin.  These bins are dynamically defined, based on the number of classes sought by the user and the similarity along the age and height axes.

## <a name="info"> Further Information </a> 

- Most of the existing documentation on how to pre-process, compile, and run the ED2 model is available in our <a href="https://github.com/EDmodel/ED2/wiki">Wiki</a>.
- For the technical description of the various ED2 model features, we suggest looking at this <a href="https://github.com/EDmodel/ED2/wiki/References-for-technical-description-of-ED-2.2">reference list</a> (<b>Tip</b>: do not skip the Supporting Information of these papers, they contain relevant details).
- For a partial list of studies that have used ED or ED2, check <a href="https://github.com/EDmodel/ED2/wiki/Publications">here</a>.

## <a name="develop"> Code Development, Pull Requests, and Commits</a>

If you plan to develop the code, please refer to the Wiki entries on <a href="https://github.com/EDmodel/ED2/wiki/Code-organization-and-design-philosophy">code organization and design philosophy</a>, to ensure your code developments are consistent with the existing model. Also, make sure that the code is thoroughly tested, and successfully passes the internal GitHub tests.

We strongly encourage that code developments are properly documented.  Please refer to the <a href="https://github.com/EDmodel/ED2/wiki/ED2-Documentation-with-Doxygen">Doxygen</a> instructions, and especially the <a href="https://github.com/EDmodel/ED2/wiki/ED2-Documentation-with-Doxygen#doxygit"> Doxygen and Git commits</a> section, so additional documentation can be automatically generated from the source code comments.

## <a name="docker"> Using Docker</a>

To use ED2 with Docker use either the `Dockerfile.gnu` for the GNU compiler (works with both x86 and arm64) or `Dockerfile.intel`  which compiled ED2 using the intel compiler (x86 only). For example the following command builds the GNU version (run this in the root of the ED2 source code):

```
docker build -t edmodel/ed2:gnu -f Dockerfile.gnu .
```

Once you have build ED2 model you can run it using:

```bash
docker run -ti --rm --ulimit stack=-1 --volume ${PWD}:/data edmodel/ed2:gnu
```

`-ti` : tells docker to use an interactive shell session
`--rm` : will remove the container after it finishes running, make sure outputs are written to /data or current folder
`--ulimit stack=-1` : set the stack to be unlimited in size, this is needed otherwise ED2 will coredump
`--volume ${PWD}:/data` : mounts the current folder to the /data folder. This is where the container starts

If no arguments are given it will run `ed2` in the /data folder. Otherwise you can pass in any arguments, for example you can use `ed2 -f Templates/ED2IN-tonzi.harvest` to run ed2 with the input file Templates/ED2IN-tonzi.harvest.
