#!/bin/sh
if [ 'x'${1} = 'xall' ]
then
   rm -fvr analy
   rm -fvr histo
   rm -fvr output
   mkdir analy
   mkdir histo
   mkdir output
   rm -fv core.* fort.* *_out.out *_lsf.out *_out.err 
   rm -fv budget_state_* thermo_state_* photo_state_*
else
   rm -fvr analy
   rm -fvr histo
   rm -fvr output
   mkdir analy
   mkdir histo
   mkdir output
   rm -fv core.* fort.* *_out.out *_lsf.out *_out.err
   rm -fv budget_* thermo_state_* photo_state_*
fi
