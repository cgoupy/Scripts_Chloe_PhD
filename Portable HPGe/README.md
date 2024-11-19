# Scripts used for the HPGe portable detector analysis

### DSCPEFileConverter.C
Code to convert the .dat files of the DAp spectrometer into root histograms. The main function DSPECFileConverter() takes two arguments: the input file name and the target histogram Bin width in keV

### Normalize_Simu_and_Conv.C
Macro to normalize in rate the simulation obtained with EnerGe:
Apply_resol_Ge() applies the resolution of the HPGe detector, it is called by Conv_Simu()
Conv_Simu() apply the resolution to the simulation obtained with EnerGe
Normalize_Simu() can be called on the file produced by Conv_Simu() to normalize it in DRU. 
!!!! WARNING !!!! you need to specify the correct flux and PLane size of the contribution in the body of the function !!!! WARNING !!!!
 
### Fit.C
The main function is Fit(TString fileData_name, TString file_simu_U238, TString file_simu_K40, TString file_simu_Th232), it fits a data file (root file produced with DSCPEFileConverter.C) with the simulations of U238, K40 and Th232, including the resolution of the detector, produced with Conv_Simu() from Normalize_Simu_and_Conv.C)

### Fit_withMu.C does the same thing as Fit.C but includes as well the muon, neutron and Rn simulation files. In the function Fit(), it is possible to fix one or multiple parameters.

### PlotSpectra.C plots the data spectra obtained with the HPGe detector and converted in root by  DSCPEFileConverter.C

### Sum_generator.C is a small macro to sum the gamma spectra obtained with the fill Gun of the VNS
