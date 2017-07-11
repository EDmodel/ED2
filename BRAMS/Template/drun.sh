#!/bin/bash
#bsub -q wofsy -J debug.edbrams -n 36 -a openmpi -Is /bin/bash
#bsub -q west1_moorcroft -R span[ptile=12] -J debug.edbrams -n 48 -a openmpi -Is /bin/bash
bsub -q moorcroft_6100 -J debug.edbrams -n 48 -a openmpi -Is /bin/bash
