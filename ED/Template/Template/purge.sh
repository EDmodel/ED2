#!/bin/bash
rm -fvr analy
rm -fvr histo
rm -fvr output
rm -fvr rdata_*
mkdir analy
mkdir histo
mkdir output
rm -fv core.* fort.* *_out.out *_lsf.out *_out.err
rm -fv budget_* thermo_state_* photo_state_*
