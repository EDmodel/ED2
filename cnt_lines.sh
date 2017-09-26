#!/bin/bash
models="ED BRAMS Ramspost RAPP"

for model in ${models}
do
   echo "========================================================================="
   echo " + Model ${model}: "
   echo "   "
   modellines=0
   direcs=`ls -1 ${model}/src`
   for dir in ${direcs}
   do
      if [ ${dir} != "test_cases" ] && [ ${dir} != "doc" ] && [ ${dir} != "preproc" ] && [ ${dir} != "post" ]
      then
         echo -n "   - Directory ${dir}: "
         files=`ls -1 ${model}/src/${dir}`
         dirlines=0
         for file in ${files}
         do
            nlines=`sed '/^ *$/ d' ${model}/src/${dir}/${file} | wc -l`
            let dirlines=${dirlines}+${nlines}
         done
         echo "${dirlines} lines"
         let modellines=${modellines}+${dirlines}
      fi
   done
   echo "   - Total: ${modellines}"
   echo "========================================================================="
   echo "   "
   echo "   "
   echo "   "
done
