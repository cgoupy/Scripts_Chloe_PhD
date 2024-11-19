# Scripts for LBR Data-MC comparison 

## Procedure to do Data-MC comparison - C.G. 10/2023
1- You have the simulation files (starting with Analyzer_*) that are in deposited energy and not normalized (it basically gathers all the simulated events).

2- You normalize all the simulations by the equivalent time of the simulation (computed from the tangent plane size, the number of primaries and the flux) this is what is done by an individual script, and it produces the file started by Histos... This is what the macro Histograms_HighStats_LBR.C is doing. 

Gammas are separated in three contributions U, Th, K that are waited by different fluxes (see my slides about the gamma measurement in the CM) and the tangent plane is 120*120 m²

Muons have a flux of 0.019/2.2 in the UGL, tangent plane is also 120*120 m²

Contamination is already normalized in the Histos file, it is the sum of 100 individually simulated contributions.

3- You have to compute the resolution function of the COV by fitting the peak in the data. In extracted somehow a linear resolution law but might be useful to check if it is really valid.

4- You have to fold the simulation spectrum with the resolution function, This is also not easy, because the normalization is not easy to control. But ok this is done by the macro Normalize_Simu_and_Conv., the function Conv_Simu()

*2- and 4- could be done together by applying the reoslution events by events when filling th  histogram, this what is done for example for the BS analysis, this is not what is done here, but the method used here can create some artefacts in the folding or some normalization issues...

5- And then finally you can compare with the data. If you are not in DRU then it is really important to have the exact same binning.

In the end, you need simulation files (Analyzer_".root), normalized files (Histos".root) and folded with convolution files (ConvLin_".root) and macro scripts.

If you just want to compare the data with the simulation you can directly use the files were the simulation is folded with a resolution (ConvLin_*)

