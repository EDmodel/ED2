#Makefile include paths.mk

# RAMS root directory.
#

EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/mlongo/repository/EDBRAMS

# MCD: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/mcd/EDBRAMS
# KIM: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/kim/ed-code/EDBRAMS
# RGK: EDBRAMS_ROOT=/Home2ln/rknox/Models/EDBRAMS
# DMM: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/dmm2/ED2/my-edbrams
# MLO: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/mlongo/EDBRAMS

BRAMS_ROOT=$(EDBRAMS_ROOT)/BRAMS
ED_ROOT=$(EDBRAMS_ROOT)/ED

# Versions.

BRAMS_VERSION=4.0.6
ED_VERSION=2.1

# Source directories.

# New paths like rams 6.x
BC=$(BRAMS_ROOT)/src/brams/bc
CATT=$(BRAMS_ROOT)/src/brams/catt
CORE=$(BRAMS_ROOT)/src/brams/core
CUPARM=$(BRAMS_ROOT)/src/brams/cuparm
ED_MIXED=$(BRAMS_ROOT)/src/brams/ed2
FDDA=$(BRAMS_ROOT)/src/brams/fdda
INIT=$(BRAMS_ROOT)/src/brams/init
IO=$(BRAMS_ROOT)/src/brams/io
ISAN=$(BRAMS_ROOT)/src/brams/isan
MASS=$(BRAMS_ROOT)/src/brams/mass
MEMORY=$(BRAMS_ROOT)/src/brams/memory
MICRO=$(BRAMS_ROOT)/src/brams/micro
MKSFC=$(BRAMS_ROOT)/src/brams/mksfc
MPI=$(BRAMS_ROOT)/src/brams/mpi
NESTING=$(BRAMS_ROOT)/src/brams/nesting
OLDGRELL=$(BRAMS_ROOT)/src/brams/oldgrell
RADIATE=$(BRAMS_ROOT)/src/brams/radiate
SIB=$(BRAMS_ROOT)/src/brams/sib
SOIL_MOISTURE=$(BRAMS_ROOT)/src/brams/soil_moisture
SURFACE=$(BRAMS_ROOT)/src/brams/surface
TEB_SPM=$(BRAMS_ROOT)/src/brams/teb_spm
TURB=$(BRAMS_ROOT)/src/brams/turb


# ED directories that will be accessed by BRAMS
ED_DYNAMICS=$(ED_ROOT)/src/dynamics
ED_IO=$(ED_ROOT)/src/io
ED_INIT=$(ED_ROOT)/src/init
ED_MEMORY=$(ED_ROOT)/src/memory
ED_MPI=$(ED_ROOT)/src/mpi
ED_UTILS=$(ED_ROOT)/src/utils


UTILS_LIB=$(BRAMS_ROOT)/src/utils/lib
EFF=$(BRAMS_ROOT)/src/utils/eff
NCARGD=$(BRAMS_ROOT)/src/utils/ncargd
UTILS_INCS=$(BRAMS_ROOT)/src/utils/include
UTILS_MODS=$(BRAMS_ROOT)/src/utils/modules

ISAN=$(BRAMS_ROOT)/src/brams/isan
ISAN_MODS=$(BRAMS_ROOT)/src/brams/isan

# Including 
