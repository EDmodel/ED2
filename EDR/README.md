# ED2 Radiative Transfer Model (RTM) code
**Corresponding author**
Toni Viskari
Brookhaven National Laboratory (BNL)
tviskari@bnl.gov


ED2-RTM a single-step model version of baseline ED2 which both reads in radiative properties and writes out the radiative outputs over a given wavelength spectrum. As ED2 divides shortwave radiation to broad categories of PAR and NIR for its calculatations, both the input and output spectra are divided similarly. Currently the ED2-RTM is only able to read in leaf reflectance and transmittance ovat, and outputs canopy reflectance over a spectrum, but other parameters and variables are easy to add to the code in a similar manner than those two.

The model is compiled in ED2/EDR/build/bin similarly to the baseline ED2 with install.sh. ED2-RTM uses the history files from the baseline ED2 model from the time of interest. At the time of writing, there have been some issues in getting ED2 to write the history files for a given hour, so all the history files are for a date starting at mid-night. Since there is no short-wave radiation then from the model perspective, it won't actually run. Thus the history file name needs to be changed from -000000 to -120000. The EDR ED2IN should be similar to the regular ED2IN file used for the run and have the clock time in the history file readin to reflect the given time.

 The runs requires the follwing files, which should all be vector files without commas:

 -lenghts.dat: This contains the lengths of the PAR and NIR parameter vectors in that order. The values are separated by just a space.

 -trans_par.dat, reflect_par.dat, trans_nir.dat and reflect_nir.dat: These all contain the values of the parameters in the file name. The PAR and NIR range vectors need to have same lengths as given in the lengths.dat file. All values are again separated by a space.

 The output files for albedo are albedo_par.dat and albedo_nir.dat.

 Both the input and output files should be located in the same folder where the model is ran from.

**Acknowledgements**
Development funded by NASA grant NNX14AH65G
