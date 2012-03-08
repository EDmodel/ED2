#Makefile include paths.mk

# RAMS root directory.

#

ED_ROOT=../..
# MCD: ED_ROOT=/n/Moorcroft_Lab/Users/mcd/EDBRAMS/ED
# KIM: ED_ROOT=/n/Moorcroft_Lab/Users/kim/ed-code/EDBRAMS/ED
# RGK: ED_ROOT=/Home2ln/rknox/Models/EDBRAMS/ED
# DMM: ED_ROOT=/n/Moorcroft_Lab/Users/dmm2/ED2/my-edbrams/ED
# MLO: ED_ROOT=/n/Moorcroft_Lab/Users/mlongo/EDBRAMS/ED

# Versions.
ED_VERSION=2.1

# Source directories.

# New paths like rams 6.x
ED_DRIVER=$(ED_ROOT)/src/driver
ED_DOC=$(ED_ROOT)/src/doc
ED_IO=$(ED_ROOT)/src/io
ED_INIT=$(ED_ROOT)/src/init
ED_MEMORY=$(ED_ROOT)/src/memory
ED_MPI=$(ED_ROOT)/src/mpi
ED_DYNAMICS=$(ED_ROOT)/src/dynamics
ED_POST=$(ED_ROOT)/src/post
ED_PREPROC=$(ED_ROOT)/src/preproc
ED_UTILS=$(ED_ROOT)/src/utils
ED_MEMTEST=$(ED_ROOT)/src/newmemtest

# Includes
ED_INCS=$(ED_ROOT)/src/include


