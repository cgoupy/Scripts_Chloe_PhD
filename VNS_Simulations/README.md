# Scripts used to normalize VNS simulations

### Histograms_HighStats_COVThr.C Normalizes the Choozerent simulations with different COV thresholds

There are two main functions:
PlotEventsCryoDet(TString SimuFilename = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", Double_t Phi=1.523, Double_t Pane_size=120*120, Bool_t NeutronQuenched=FALSE, TString Algo="Z")

PlotEventsCryoDet_Z_T(TString SimuFilenameZ = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", TString SimuFilenameT = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", Double_t Phi=1.523, Double_t Pane_size=120*120,  Double_t Pane_size_simuT=120*120, Bool_t NeutronQuenched=FALSE)

The first one takes a simulation file as an argument (analyzed with T or Z, you can tune with the last parameter), the flux of the contribution, the plane size, it is possible to quench the simulation by enabling the option.

The second one adds up two simulations files analyzed with Z and T (which can not be automatically done with hadd).

### plots_Contam_Simu.ipynb
Python script that plots the pie charts of the contamination in the VNS
