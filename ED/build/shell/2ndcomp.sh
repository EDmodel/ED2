#!/bin/bash

if [ "x${1}" == "x" ]
then
   edroot="$(pwd)/../.."
else
   edroot=${1}
fi

edsrc="${edroot}/src/*/*.f90 ${edroot}/src/*/*.F90 ${edroot}/src/*/*.c"

for srcnow in ${edsrc}
do
   srcbase=$(basename ${srcnow})
   srcpref=$(basename ${srcbase} .f90)
   srcpref=$(basename ${srcpref} .F90)
   srcpref=$(basename ${srcpref} .c  )

   /bin/rm -fv ${srcbase} ${srcpref}.o ${srcpref}.mod
done
