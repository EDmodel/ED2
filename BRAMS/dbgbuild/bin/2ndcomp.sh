#!/bin/bash

if [ "x${1}" == "x" ]
then
   modroot="$(pwd)/../../.."
else
   modroot="${1}"
fi

modsrc="${modroot}/BRAMS/src/*/*.f90 
        ${modroot}/BRAMS/src/*/*.F90
        ${modroot}/BRAMS/src/*/*.c
        ${modroot}/ED/src/*/*.f90 
        ${modroot}/ED/src/*/*.F90
        ${modroot}/ED/src/*/*.c"

for srcnow in ${modsrc}
do
   srcbase=$(basename ${srcnow})
   srcpref=$(basename ${srcbase} .f90)
   srcpref=$(basename ${srcpref} .F90)
   srcpref=$(basename ${srcpref} .c  )

   /bin/rm -fv ${srcbase} ${srcpref}.o ${srcpref}.mod
done
