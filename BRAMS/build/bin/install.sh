#!/bin/bash
# Create ncarg-dummy library **************************
#
# Compile ncar, if you are using the NCAR dummy libraries
#
# by Daniel Merli Lamosa PAD/CPTEC/INPE
# *****************************************************

home_dir=`pwd`

# Set HDF path -----------------------------------------
Find_HDF() {

 echo "Checking HDF......"
 
 TEMP1=`find ${home_dir}/.hdf4_libs -name "*.a" -print 2>/dev/null`
 TEMP=`echo $TEMP1 | grep libsz.a | grep libz.a | grep libjpeg.a |grep libdf.a |grep libmfhdf.a | wc -l | awk '{print $1}'`

 if [ $TEMP -eq 0 ] ; then
   echo
   echo "ERROR! HDF not detected!"
   echo
   echo "Run \"configure\" script to set a corrrect HDF PATH"
   echo
   echo "You need this libs:"
   echo "  - libsz.a"
   echo "  - libz.a"
   echo "  - libjpeg.a"
   echo "  - libdf.a"
   echo "  - libmfhdf.a"
   echo
   echo "For more information see www.cptec.inpe.br/brams"
   echo
   echo "Tank you to use brams"
   echo "Brams team"
 else
   echo "HDF OK!"
 fi

}
# ------------------------------------------------------


Find_include() {

  echo
  echo "Checking configure files....."

  # Check if include.mk.opt is set
  ls include.mk.opt &> /dev/null 

  if [ $? -ne 0 ] ; then
    echo "ERROR! File include.mk.opt in not define!"
    echo
    echo "Run \"configure\" script to set this file"
    echo
    exit 8
  fi 
  echo "files configured!"

}

# Compile Brams source --------------------------
Make() {

make -f Make_model OPT=opt $1

}
# -----------------------------------------------


# Main script -----------------------------------
if [ "$1" == "clean" ] ; then
  ls include.mk.opt &> /dev/null
  if [ $? -eq 0 ] ; then
    Make clean 
  fi
else
  Find_HDF
  Find_include
#  echo "Cleanning old files..."
#  Make clean
  echo "Installing Brams.."
  Make
fi
# -----------------------------------------------
