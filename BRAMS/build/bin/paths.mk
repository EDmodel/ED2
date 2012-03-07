#Makefile include paths.mk

# RAMS root directory.
#
EDBRAMS_ROOT=/home/rknox/Models/r28X/EDBRAMS

# MCD: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/mcd/EDBRAMS
# KIM: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/kim/ed-code/EDBRAMS
# RGK: EDBRAMS_ROOT=/Home2ln/rknox/Models/EDBRAMS
# DMM: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/dmm2/ED2/my-edbrams
# MLO: EDBRAMS_ROOT=/n/Moorcroft_Lab/Users/mlongo/EDBRAMS

#EDBRAMS_ROOT=/home/rknox/Models/Mainline/EDBRAMS

BRAMS_ROOT=$(EDBRAMS_ROOT)/BRAMS
ED_ROOT=$(EDBRAMS_ROOT)/ED

# Versions.

BRAMS_VERSION=4.0.6
ED_VERSION=2.1

# Source directories.

# New paths like rams 6.x
BC=$(BRAMS_ROOT)/src/bc
CATT=$(BRAMS_ROOT)/src/catt
CORE=$(BRAMS_ROOT)/src/core
CUPARM=$(BRAMS_ROOT)/src/cuparm
ED_MIXED=$(BRAMS_ROOT)/src/ed2
FDDA=$(BRAMS_ROOT)/src/fdda
INIT=$(BRAMS_ROOT)/src/init
IO=$(BRAMS_ROOT)/src/io
ISAN=$(BRAMS_ROOT)/src/isan
MASS=$(BRAMS_ROOT)/src/mass
MEMORY=$(BRAMS_ROOT)/src/memory
MICRO=$(BRAMS_ROOT)/src/micro
MKSFC=$(BRAMS_ROOT)/src/mksfc
MPI=$(BRAMS_ROOT)/src/mpi
NESTING=$(BRAMS_ROOT)/src/nesting
OLDGRELL=$(BRAMS_ROOT)/src/oldgrell
RADIATE=$(BRAMS_ROOT)/src/radiate
SIB=$(BRAMS_ROOT)/src/sib
SOIL_MOISTURE=$(BRAMS_ROOT)/src/soil_moisture
SURFACE=$(BRAMS_ROOT)/src/surface
TEB_SPM=$(BRAMS_ROOT)/src/teb_spm
TURB=$(BRAMS_ROOT)/src/turb
MNTADVEC=$(BRAMS_ROOT)/src/mnt_advec


# ED directories that will be accessed by BRAMS
ED_DYNAMICS=$(ED_ROOT)/src/dynamics
ED_IO=$(ED_ROOT)/src/io
ED_INIT=$(ED_ROOT)/src/init
ED_MEMORY=$(ED_ROOT)/src/memory
ED_MPI=$(ED_ROOT)/src/mpi
ED_UTILS=$(ED_ROOT)/src/utils


UTILS_LIB=$(BRAMS_ROOT)/src/lib
EFF=$(BRAMS_ROOT)/src/eff
NCARGD=$(BRAMS_ROOT)/src/ncargd
UTILS_INCS=$(BRAMS_ROOT)/src/include
UTILS_MODS=$(BRAMS_ROOT)/src/modules

ISAN=$(BRAMS_ROOT)/src/isan
ISAN_MODS=$(BRAMS_ROOT)/src/isan

# Including 
