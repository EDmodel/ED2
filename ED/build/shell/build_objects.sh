#!/bin/bash

# Add your build commands here
SOURCES=$( (cd ../../.. && find ED/src -name "*.[Ff]90" -and -not -path "ED/src/preproc/*"))

# Remove ED/src/driver/edmain.F90 from OBJECTS
SOURCES=$(echo $SOURCES | sed 's/ED\/src\/driver\/edmain.F90//')

echo SOURCES = $SOURCES ED/src/utils/utils_c.c > sources.mk

rm -f dependency.mk

# Build dependency rules for object files
for src in $SOURCES; do
  # Mod depends on the source file building
  SRC_MOD_NAME=$(echo $(basename "${src%.[Ff]90}.mod"))
  echo $SRC_MOD_NAME: ${src%.[Ff]90}.o >> dependency.mk

  # Find modules used in the source file and add as dependency for the object file
  MODS_USED=$(cd ../../.. && grep -i '^\s*use\s\+' $src | awk '{print $2}' | sed 's/,.*//' | sort | uniq)

  # Update MODS_USED to include only modules that has a corresponding
  # source file with the pattern ED/src/*/$mod.[Ff]90
  MODS_USED=$(for mod in $MODS_USED; do
    if [ -f ../../../ED/src/*/$mod.F90 -o -f ../../../ED/src/*/$mod.f90 ]; then
      echo $mod
    fi
  done)

  MODS_USED_MOD=$(for line in $MODS_USED; do echo -n "$line.mod "; done)

  length=$(echo $MODS_USED_MOD | wc -w)
  if [ $length -gt 0 ]; then
    echo ${src%.[Ff]90}.o: $MODS_USED_MOD >> dependency.mk
  fi;
done
