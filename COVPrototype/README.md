# Macros_COVProto

Some macros used to analyze COV proto simulations.

- _PeakFitter.C_ \
Hard coded fit of spectra : plot calibration curve, resolution curve.\
Function to propagate data calibration systematic errors to the simulation.

- _Recalib.C_ \
Macro to Recalibrate already calibrated data with Eddy's code.

- _Shift.C_ \
Macro to shift a spectrum by Drawing random values in the spectrum to shift.

- _Spec_Analysis.cc_ \
Basis code to analyse the simulations (Scale in cps and dru, form the data, etc).\
Calls Shift.C and Recal.C.\
To display the list of functions:\
`>> .L Spec_Analysis.cc`\
`>> help()`

- _FitNorm.cc_ \
Macro to fit the simulation spectra to the data\
Some functions from Spec_Analysis copy-pasted\
To display the list of functionss:\
`>> .L FitNorm.cc` \
`>> help()`

- _MyCOVcodes contains old macro used to analyze orsay data with the bin_to_bin method. This is now depreciated: use COVLab instead. 
