# Scripts and data of the Bonner Sphere measurement analysis

### Most of these are root notebooks, to be launched with:
`root --notebook`

>> BS_Analysis_GoodCalib.ipynb does the full Bonner sphere analysis (Surface and VNS) => construct the model from the simulations and fit the data. 

>> BS_Calib.C (root script) defines set of functions and classes to calibrate from the data-MC comparison the BS spectrometer with Saclay measurements, main function is BS_Calib(). iIt plots the calibration plots of the BS from the fit of the resolution and the calibration of the simulation on data spectra. 
    After updating the paths in the main function (lines 1300 - 1306) can be used by launching:
    `root -l BS_Calib.C`

>> BS_Calibration.ipynb is a notebook doing the calibration of the BS with the positon of the end of the source spectra. 

>> BS_Compare_BS_Simulation.ipynb notebook doing an aboslute data-MC comparison of bckg and source measurements at Saclay. 
