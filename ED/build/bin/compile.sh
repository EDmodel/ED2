#!/bin/sh

if [ 'z'${1} == 'zclean' ]
then
  ./install.sh clean
else
  bsub -Is -q interact ./install.sh
fi
