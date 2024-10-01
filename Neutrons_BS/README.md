# Scripts and data of the Bonner Sphere measurement analysis

Most of these are root notebooks, to be launched with:
`root --notabook`

>> BS_Analysis_GoodCalib.ipynb does the full Bonner sphere analysis (Surface and VNS) => construct the model from the simulations and fit the data. 

>> BS_Calib.C (root script) defines set of functions and classes to calibrate from the data-MC comparison the BS spectrometer with Saclay measurements, main funciton is BS_Calib()

>> BS_Calibration.ipynb is a notebook doing the calibration of the BS with the positon of the end of the source spectra. 

>> BS_Compare_BS_Simulation.ipynb notebook doing an aboslute data-MC comparison of bckg and source measurements at Saclay. 

