#!/bin/sh
if [ 'x'${1} == 'xclean' ]
then
  make OPT=opt clean
else
  make OPT=opt
fi
