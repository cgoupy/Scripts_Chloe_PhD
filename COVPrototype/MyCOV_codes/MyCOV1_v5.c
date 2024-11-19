#include<iostream>
#include<fstream>
#include<algorithm>
#include<cmath>
#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"

using namespace std;
unsigned short int ByteSwap(unsigned short int *x);

Double_t Ge_Pulse(Double_t *x, Double_t *par) {

  Double_t Pulse = 0.;

  if (x[0] < par[3]) {Pulse = par[0] + par[1]*x[0];}
  else {Pulse = par[0] + par[1]*x[0] + par[2]*par[4]*(exp(-(x[0] - par[3])/par[4]) - exp(-(x[0] - par[3])/par[5]))/(par[4] - par[5]);}

  return Pulse;
}

float ReverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

void MyCOV1()

{

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Open output file for histograms                                         // 
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////
    
  char myfile_name[ ] = "MyCOV1.root";
  Double_t Cal1 =  .0009533;
  Double_t Cal2 =  .0006513;
  Double_t Off1 = -.001038;
  Double_t Off2 = -.000835;
  Double_t CalNTD =  .0002361;
  Double_t OffNTD =  .0011361;
  Double_t DTacc = 30.;
  Int_t Tbin_1kHz_1 = 1;
  Int_t Tbin_1kHz_2 = 0;


  cout << " Open root file: " << myfile_name << " " << endl;
  TFile fHistos(myfile_name,"recreate");

  TH1F *h1_RunTime, *h1_Thresholds, *h1_Gains;
  TH1F *h1_Raw1, *h1_Raw2, *h1_Raw3, *h1_Raw4, *h1_Raw5;
  TH1F *h1_COV1, *h1_COV2, *h1_DIF1, *h1_DIF2, *h1_ASM1, *h1_ASM2, *h1_DIF_ASM1, *h1_DIF_ASM2;
  TH1F *h1_Peak1, *h1_Peak2, *h1_Peak1_MeV, *h1_Peak2_MeV, *h1_Peak1_4MeV, *h1_Peak2_4MeV, *h1_Peak1_4MeV_AVG, *h1_Peak2_4MeV_AVG;
  TH1F *h1_BSLM1, *h1_BSLM2, *h1_BSLS1, *h1_BSLS2, *h1_BSL_NTD_M, *h1_BSL_NTD_S;
  TH1F *h1_bsl_M1, *h1_bsl_S1, *h1_bsl_M2, *h1_bsl_S2;
  TH1F *h1_DeltaTCOV, *h1_DeltaT_NTD_COV1, *h1_DeltaT_NTD_COV2, *h1_DeltaT_NTD_COV1and2;
  TH1F *h1_Peak1_keV, *h1_Peak2_keV, *h1_Peak1_C_keV, *h1_Peak2_C_keV, *h1_Peak1_C_250keV, *h1_Peak2_C_250keV, *h1_Peak1_C_MeV, *h1_Peak2_C_MeV;
  TH1F *h1_Peak1_100mV, *h1_Peak2_100mV, *h1_Peak1_100keV, *h1_Peak2_100keV;
  TH1F *h1_DIF_NTD, *h1_NTD, *h1_NTD_Amplitude, *h1_NTD_Amplitude_keV;
  TH2F *h2_Width_Amplitude, *h2_NTD_T1_T2, *h2_DeltaT_NTD_COV1_P_NTD, *h2_DeltaT_NTD_COV2_P_NTD, *h2_ASM1_Amplitude, *h2_ASM2_Amplitude;
  TH1F *h1_NTD_1_MeV, *h1_NTD_2_MeV, *h1_NTD_3_MeV, *h1_NTD_4_MeV, *h1_NTD_5_MeV, *h1_NTD_6_MeV;
  TH1F *h1_T_NTD;
  TH1F *h1_HPF, *h1_LPF;
  TH1F *h1_COV2_Test, *h1_Fit;
  TH1F *h1_COV1_AVG, *h1_DIF1_AVG, *h1_Peak1_100mV_AVG, *h1_Peak1_100keV_AVG, *h1_Peak1_100mV_1kHz_AVG, *h1_Peak1_100keV_1kHz_AVG;
  TH1F *h1_COV2_AVG, *h1_DIF2_AVG, *h1_Peak2_100mV_AVG, *h1_Peak2_100keV_AVG, *h1_Peak2_100mV_1kHz_AVG, *h1_Peak2_100keV_1kHz_AVG;
  TH1F *h1_BSLM1_AVG, *h1_BSLS1_AVG, *h1_BSLM2_AVG, *h1_BSLS2_AVG;
  TH1F *h1_P_COV1, *h1_P_COV2, *h1_P_noCOV1, *h1_P_noCOV2;

  Int_t Ndiv = 2;

  h1_RunTime = new TH1F("h1_RunTime"," ",1,0.,1.);
  h1_Thresholds = new TH1F("h1_Thresholds"," ",3,0.,3.);
  h1_Gains = new TH1F("h1_Gains"," ",3,0.,3.);
  h1_Raw1 = new TH1F("h1_Raw1"," ",100000,0.,100000.);
  h1_Raw2 = new TH1F("h1_Raw2"," ",100000,0.,100000.);
  h1_Raw3 = new TH1F("h1_Raw3"," ",100000,0.,100000.);
  h1_Raw4 = new TH1F("h1_Raw4"," ",100000,0.,100000.);
  h1_Raw5 = new TH1F("h1_Raw5"," ",100000,0.,100000.);
  h1_COV1 = new TH1F("h1_COV1"," ",100000,0.,100000.);
  h1_COV2 = new TH1F("h1_COV2"," ",100000,0.,100000.);
  h1_DIF1 = new TH1F("h1_DIF1"," ",100000,0.,100000.);
  h1_DIF2 = new TH1F("h1_DIF2"," ",100000,0.,100000.);
  h1_ASM1 = new TH1F("h1_ASM1"," ",100000,0.,100000.);
  h1_ASM2 = new TH1F("h1_ASM2"," ",100000,0.,100000.);
  h1_DIF_ASM1 = new TH1F("h1_DIF_ASM1"," ",100000,0.,100000.);
  h1_DIF_ASM2 = new TH1F("h1_DIF_ASM2"," ",100000,0.,100000.);
  h1_COV1_AVG = new TH1F("h1_COV1_AVG"," ",100000/Ndiv,0.,100000.);
  h1_COV2_AVG = new TH1F("h1_COV2_AVG"," ",100000/Ndiv,0.,100000.);
  h1_DIF1_AVG = new TH1F("h1_DIF1_AVG"," ",100000/Ndiv,0.,100000.);
  h1_DIF2_AVG = new TH1F("h1_DIF2_AVG"," ",100000/Ndiv,0.,100000.);
  h1_DIF_NTD = new TH1F("h1_DIF_NTD"," ",100000,0.,100000.);
  h1_Peak1 = new TH1F("h1_Peak1"," ",10000,0.,20.);
  h1_Peak2 = new TH1F("h1_Peak2"," ",10000,0.,20.);
  h1_Peak1_keV = new TH1F("h1_Peak1_keV"," ",10000,Off1,20./Cal1);
  h1_Peak2_keV = new TH1F("h1_Peak2_keV"," ",10000,Off2,20./Cal2);
  h1_Peak1_C_keV = new TH1F("h1_Peak1_C_keV"," ",2000,0.,4000.);
  h1_Peak2_C_keV = new TH1F("h1_Peak2_C_keV"," ",2000,0.,4000.);
  h1_Peak1_C_250keV = new TH1F("h1_Peak1_C_100keV"," ",250,0.,250.);
  h1_Peak2_C_250keV = new TH1F("h1_Peak2_C_100keV"," ",250,0.,250.);
  h1_Peak1_MeV = new TH1F("h1_Peak1_MeV"," ",2000,0.,20.);
  h1_Peak2_MeV = new TH1F("h1_Peak2_MeV"," ",2000,0.,20.);
  h1_Peak1_4MeV = new TH1F("h1_Peak1_4MeV"," ",2000,0.,4.);
  h1_Peak2_4MeV = new TH1F("h1_Peak2_4MeV"," ",2000,0.,4.);
  h1_Peak1_C_MeV = new TH1F("h1_Peak1_C_MeV"," ",2000,0.,20.);
  h1_Peak2_C_MeV = new TH1F("h1_Peak2_C_MeV"," ",2000,0.,20.);
  h1_Peak1_100mV = new TH1F("h1_Peak1_100mV"," ", 300,0.,0.1);
  h1_Peak2_100mV = new TH1F("h1_Peak2_100mV"," ", 300,0.,0.1);
  h1_Peak1_100keV = new TH1F("h1_Peak1_100keV"," ",300,Off1,0.1/Cal1);
  h1_Peak2_100keV = new TH1F("h1_Peak2_100keV"," ",300,Off2,0.1/Cal2);
  h1_Peak1_100mV_AVG = new TH1F("h1_Peak1_100mV_AVG"," ", 300,0.,0.1);
  h1_Peak2_100mV_AVG = new TH1F("h1_Peak2_100mV_AVG"," ", 300,0.,0.1);
  h1_Peak1_100keV_AVG = new TH1F("h1_Peak1_100keV_AVG"," ",300,Off1,0.1/Cal1);
  h1_Peak2_100keV_AVG = new TH1F("h1_Peak2_100keV_AVG"," ",300,Off1,0.1/Cal2);
  h1_Peak1_100mV_1kHz_AVG = new TH1F("h1_Peak1_100mV_1kHz_AVG"," ", 300,0.,0.1);
  h1_Peak2_100mV_1kHz_AVG = new TH1F("h1_Peak2_100mV_1kHz_AVG"," ", 300,0.,0.1);
  h1_Peak1_100keV_1kHz_AVG = new TH1F("h1_Peak1_100keV_1kHz_AVG"," ",300,Off1,0.1/Cal1);
  h1_Peak2_100keV_1kHz_AVG = new TH1F("h1_Peak2_100keV_1kHz_AVG"," ",300,Off1,0.1/Cal2);
  h1_Peak1_4MeV_AVG = new TH1F("h1_Peak1_4MeV_AVG"," ",2000,0.,4.);
  h1_Peak2_4MeV_AVG = new TH1F("h1_Peak2_4MeV_AVG"," ",2000,0.,4.);
  h1_BSLM1 = new TH1F("h1_BSLM1"," ",  200,-5.0,5.0);
  h1_BSLM2 = new TH1F("h1_BSLM2"," ",  200,-5.0,5.0);
  h1_BSLS1 = new TH1F("h1_BSLS1"," ",  100, 0.0,5.0);
  h1_BSLS2 = new TH1F("h1_BSLS2"," ",  100, 0.0,5.0);
  h1_BSLM1_AVG = new TH1F("h1_BSLM1_AVG"," ",  200,-5.0,5.0);
  h1_BSLS1_AVG = new TH1F("h1_BSLS1_AVG"," ",  100, 0.0,5.0);
  h1_BSLM2_AVG = new TH1F("h1_BSLM2_AVG"," ",  200,-5.0,5.0);
  h1_BSLS2_AVG = new TH1F("h1_BSLS2_AVG"," ",  100, 0.0,5.0);
  h1_bsl_M1 = new TH1F("h1_bsl_M1"," ",  200,-5.0,5.0);
  h1_bsl_M2 = new TH1F("h1_bsl_M2"," ",  200,-5.0,5.0);
  h1_bsl_S1 = new TH1F("h1_bsl_S1"," ",  100, 0.0,5.0);
  h1_bsl_S2 = new TH1F("h1_bsl_S2"," ",  100, 0.0,5.0);
  h1_BSL_NTD_M = new TH1F("h1_BSL_NTD_M"," ",  200,-2.0,2.0);
  h1_BSL_NTD_S = new TH1F("h1_BSL_NTD_S"," ",  200, 0.,10.0);
  h1_DeltaTCOV = new TH1F("h1_DeltaTCOV"," ",  100, -50.,50.);
  h1_NTD = new TH1F("h1_NTD"," ",100000,0.,100000.);
  h1_NTD_Amplitude      = new TH1F("h1_NTD_Amplitude"    ," ",1000,      0.,        2.);
  h1_NTD_Amplitude_keV  = new TH1F("h1_NTD_Amplitude_keV"," ",1000,  OffNTD, 2./CalNTD);
  h1_T_NTD = new TH1F("h1_T_NTD"," ",2000,0.,100000.);
  h1_DeltaT_NTD_COV1 = new TH1F("h1_DeltaT_NTD_COV1"," ", 250, -25.,25.);
  h1_DeltaT_NTD_COV2 = new TH1F("h1_DeltaT_NTD_COV2"," ", 250, -25.,25.);
  h1_DeltaT_NTD_COV1and2 = new TH1F("h1_DeltaT_NTD_COV1and2"," ", 250, -25.,25.);
  h1_NTD_1_MeV  = new TH1F("h1_NTD_1_MeV"," ",2000,  0., 20.);
  h1_NTD_2_MeV  = new TH1F("h1_NTD_2_MeV"," ",2000,  0., 20.);
  h1_NTD_3_MeV  = new TH1F("h1_NTD_3_MeV"," ",2000,  0., 20.);
  h1_NTD_4_MeV  = new TH1F("h1_NTD_4_MeV"," ",2000,  0., 20.);
  h1_NTD_5_MeV  = new TH1F("h1_NTD_5_MeV"," ",2000,  0., 20.);
  h1_NTD_6_MeV  = new TH1F("h1_NTD_6_MeV"," ",2000,  0., 20.);
  h1_P_COV1 = new TH1F("h1_P_COV1"," ",2000,0.,20.);
  h1_P_COV2 = new TH1F("h1_P_COV2"," ",2000,0.,20.);
  h1_P_noCOV1 = new TH1F("h1_P_noCOV1"," ",2000,0.,20.);
  h1_P_noCOV2 = new TH1F("h1_P_noCOV2"," ",2000,0.,20.);
  h2_Width_Amplitude  = new TH2F("h2_Width_Amplitude"," ",50,  0.,50.,1000.,0.,5.0);
  h2_ASM1_Amplitude  = new TH2F("h2_ASM1_Amplitude"," ",100, 0.,20.,100.,-0.50,0.50);
  h2_ASM2_Amplitude  = new TH2F("h2_ASM2_Amplitude"," ",100, 0.,20.,100.,-0.50,0.50);
  h2_NTD_T1_T2  = new TH2F("h2_NTD_T1_T2"," ",100, -50.,50.,100.,-50.,50.);
  h2_DeltaT_NTD_COV1_P_NTD  = new TH2F("h2_DeltaT_NTD_COV1_P_NTD"," ",100,  -2., 2.,20,-10.,10.);
  h2_DeltaT_NTD_COV2_P_NTD  = new TH2F("h2_DeltaT_NTD_COV2_P_NTD"," ",100,  -2., 2.,20,-10.,10.);
  h1_HPF = new TH1F("h1_HPF"," ",100000,0.,100000.);
  h1_LPF = new TH1F("h1_LPF"," ",100000,0.,100000.);
  h1_COV2_Test = new TH1F("h1_COV2_Test"," ",100000,0.,100000.);
  h1_Fit = new TH1F("h1_Fit"," ",100000,0.,100000.);

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

  Bool_t Flag_Little_Endian = true;

  float Raw1[100000], Raw2[100000], Raw3[100000], Raw4[100000], Raw5[100000];
  float COV1[100000], COV2[100000], Sum_DIF1[100000], Sum_DIF2[100000];
  float ASM1[100000], ASM2[100000];
  // Sum_DIF_ASM1[100000];
  //, Sum_DIF_ASM2[100000];
  float COV1_AVG[50000], COV2_AVG[50000], Sum_DIF1_AVG[50000], Sum_DIF2_AVG[50000];
  float Peak_NTD[100000];
  Double_t T_COV1[10000], T_COV2[10000], T_NTD[10000], T_COV1and2[10000];
  float P_COV1[10000], P_COV2[10000], P_NTD[10000];
  float DIF1_bias[1000], DIF2_bias[1000], DIF1_AVG_bias[1000], DIF2_AVG_bias[1000];
  float Tab1_E[50],Tab2_E[50];
  Int_t Tab1_K[50], Tab2_K[50], Tab1_OK[50], Tab2_OK[50];
  Int_t itab1, itab2;
  Int_t nbsl;
  Double_t bsl_M, bsl_S;
   
  for (int k = 0; k < 10000; k++) {
    T_COV1[k] = -1.;
    T_COV2[k] = -1.;
    T_NTD[k]  = -1.;
    T_COV1and2[k] = -1.;
  }
   
  float DIF1, DIF2;
  Double_t Sn, Sx, Sy, Sxx, Sxy;
  Double_t aal, bbl, aau, bbu;

  signed short ampl;
  int i;
  unsigned iw, iw_ms, iw_sel, iw_NTD, nw1, nw2;
  unsigned iw_NTD_sel = 0;
  Bool_t Flag_NTD = true;
  Int_t NS = 100000;
  Int_t NS_NTD = 10000000;
  float TPeak1[10000], TPeak2[10000]; 
  float Delta_T;

  iw_sel = 30;

  cout << " Start MyCOV1 !!! " << endl;

// Read Ge detector data

  ifstream rf1("/sps/lbno/edoardo/Actuator/Logged Data_2021_05_12_17_27_27.Bin1", ios::in | ios::binary);
  if (!rf1) {cout << "Cannot open input file 1!" << endl;}

  ifstream rf2("/sps/lbno/edoardo/Actuator/Logged Data_2021_05_12_17_27_27.Bin2", ios::in | ios::binary);
  if (!rf2) {cout << "Cannot open input file 2!" << endl;}

  ifstream rf3("/sps/lbno/edoardo/Actuator/Logged Data_2021_05_12_17_27_27.Bin3", ios::in | ios::binary);
  if (!rf3) {cout << "Cannot open input file 3!" << endl;}

  ifstream rf4("/sps/lbno/edoardo/Actuator/Logged Data_2021_05_12_17_27_27.Bin4", ios::in | ios::binary);
  if (!rf4) {cout << "Cannot open input file 4!" << endl;}

  ifstream rf5("/sps/lbno/edoardo/Actuator/Logged Data_2021_05_12_17_27_27.NTD", ios::in | ios::binary);
  if (!rf5) {cout << "Cannot open input file 5!" << endl;}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

  ofstream wf1("Logged Data_2021_05_12_17_27_27.Raw1", ios::out | ios::binary);
  if (!wf1) {cout << "Cannot open output file 1!" << endl;}

  ofstream wf2("Logged Data_2021_05_12_17_27_27.Raw2", ios::out | ios::binary);
  if (!wf2) {cout << "Cannot open output file 2!" << endl;}

  ofstream wf3("Logged Data_2021_05_12_17_27_27.Raw3", ios::out | ios::binary);
  if (!wf3) {cout << "Cannot open output file 3!" << endl;}

  ofstream wf4("Logged Data_2021_05_12_17_27_27.Raw4", ios::out | ios::binary);
  if (!wf4) {cout << "Cannot open output file 4!" << endl;}

  ofstream wf5("Output_2021_05_12_17_27_27.NTD", ios::out | ios::binary);
  if (!wf5) {cout << "Cannot open output file 5!" << endl;}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

  unsigned short int x, hibyte, lobyte;
  unsigned kt;
  float rx;

  Bool_t Flag_Peak1[100000], Flag_Peak2[100000], Flag_DIF1[100000], Flag_DIF2[100000];
  Bool_t Flag_Peak1_AVG[50000], Flag_Peak2_AVG[50000], Flag_DIF1_AVG[50000], Flag_DIF2_AVG[50000];
  Bool_t Flag_Peak1_1kHz[100000], Flag_Peak2_1kHz[100000], Flag_Peak1_1kHz_AVG[50000], Flag_Peak2_1kHz_AVG[50000];

  Double_t RT1_Threshold, RT2_Threshold;
  Double_t Threshold_keV_1 =  4.0;
  Double_t Threshold_keV_2 = 10.0;
  Double_t Peak1_Threshold =  4.0; 
  Double_t Peak2_Threshold = 10.0; 
  Double_t Peak1_Upper_Threshold =  7.0; 
  Double_t Peak2_Upper_Threshold = 15.0; 
  Double_t NTD_Threshold_keV = 25.0;
  //  Double_t NTD_Threshold = 0.0060;
  Double_t NTD_Threshold = NTD_Threshold_keV*CalNTD;
  Double_t Gain_1 = 1000.;
  Double_t Gain_2 = 1000.;
  Double_t NTD_Gain = 8000.;

  h1_Thresholds->SetBinContent(1,Threshold_keV_1);
  h1_Thresholds->SetBinContent(2,Threshold_keV_2);
  h1_Thresholds->SetBinContent(3,NTD_Threshold_keV);

  h1_Gains->SetBinContent(1,Gain_1);
  h1_Gains->SetBinContent(2,Gain_2);
  h1_Gains->SetBinContent(3,NTD_Gain);

  Int_t max_bin1, max_bin2;
  Int_t np1 = 0;
  Int_t np2 = 0;
  Int_t nt1 = 0;
  Int_t nt2 = 0;
  Int_t np_NTD = 0;
  Int_t ip_NTD = 0;
  Int_t ip1_ms = 0;
  Int_t ip2_ms = 0;
  Int_t ip1and2_ms = 0;
  Int_t ip_COV1, ip_COV2, ip_COV1and2;  

  Int_t j1, j2, i_NTD, j_NTD, k_min, k_max;
  Int_t nBSL1, nBSL2, nBSL_NTD;
  Double_t BSL_M1, BSL_S1, BSL_M2, BSL_S2, BSL_NTD_M, BSL_NTD_S;
  Int_t nBSL1_AVG, nBSL2_AVG;
  Double_t BSL_M1_AVG, BSL_S1_AVG, BSL_M2_AVG, BSL_S2_AVG;
  float NTD_Avg, DIF_NTD, NTD[100000], HPF[100000], LPF[100000];
 
  Double_t Sum1, Sum2, Sum3, Sum_ASM1, Sum_ASM2,Sum_ASM3;  

  Bool_t Flag_HPF_COV = false;
  Double_t alpha_COV = 1.;
  Double_t HPFreq_COV = 100.;
  if (Flag_HPF_COV) alpha_COV = (1./HPFreq_COV)/((1./HPFreq_COV) + .00001);

  Bool_t Flag_HPF_NTD = false;
  Double_t alpha_NTD = 1.;
  Double_t HPFreq_NTD = 1.;
  if (Flag_HPF_NTD) alpha_NTD = (1./HPFreq_NTD)/((1./HPFreq_NTD) + .001);

  Double_t beta_NTD = 1./(1000. + 1.);

  RT1_Threshold = Threshold_keV_1*Cal1;
  RT2_Threshold = Threshold_keV_2*Cal2;

  i = 0;
  nw1 = 0;
  nw2 = 0;
  iw = 0;
  iw_ms = 0;
  ip1_ms = 0;
  ip2_ms = 0;

  i_NTD = 0;
  j_NTD = 0;
 
  Int_t nnn = 0;

  while (!rf1.eof() && !rf2.eof() && !rf3.eof() && !rf4.eof() && !rf5.eof()) {
    //while (!rf1.eof() && !rf2.eof() && !rf3.eof() && !rf4.eof() && !rf5.eof() && nw1 <  40*NS) {
    rf1.read((char *) &x, sizeof x);
    hibyte = (x & 0xff00);
    lobyte = (x & 0xff); 
    ampl = (lobyte << 8) | (hibyte >> 8); 
    //wf1.write((char *) &ampl, sizeof ampl);
    if (i < NS) Raw1[i] = ampl;
    rf2.read((char *) &x, sizeof x);
    hibyte = (x & 0xff00);
    lobyte = (x & 0xff); 
    ampl = (lobyte << 8) | (hibyte >> 8); 
    //wf2.write((char *) &ampl, sizeof ampl);
    if (i < NS) {
      Raw2[i] = ampl;
      COV1[i] = 20.*(Raw2[i] - Raw1[i])/65536.;
      ASM1[i] = 20.*(Raw1[i] + Raw2[i])/65536.;
    }
    rf3.read((char *) &x, sizeof x);
    hibyte = (x & 0xff00);
    lobyte = (x & 0xff); 
    ampl = (lobyte << 8) | (hibyte >> 8); 
    //wf3.write((char *) &ampl, sizeof ampl);
    if (i < NS) Raw3[i] = ampl;
    rf4.read((char *) &x, sizeof x);
    hibyte = (x & 0xff00);
    lobyte = (x & 0xff); 
    ampl = (lobyte << 8) | (hibyte >> 8); 
    //wf4.write((char *) &ampl, sizeof ampl);
    if (i < NS) { 
      Raw4[i] = ampl;
      COV2[i] = 20.*(Raw4[i] - Raw3[i])/65536.;
      ASM2[i] = 20.*(Raw3[i] + Raw4[i])/65536.;
    }
    rf5.read((char *) &rx, sizeof rx);
    if (i < NS) { 
      Raw5[i] = ReverseFloat(rx);
    }

//
// NTD : move from 100kHz to 1kHz sampling
// 
    if (i_NTD < NS_NTD) {
      if ((i_NTD + 1)%100 == 0) {
	NTD_Avg = NTD_Avg + Raw5[i]/100.;
        NTD[j_NTD] = NTD_Avg;
        j_NTD++;
        NTD_Avg = 0.;
      } 
      else {NTD_Avg = NTD_Avg + Raw5[i]/100.;}
    } 

    i++;
    iw_ms++;
    i_NTD++;

    if (nw1 < 4000000000) {nw1++;}
    else {nw2++;}
    
    if (i == NS) {
      i = 0;
      if ((nw1%10000000 == 0 && nw2 == 0) || (nw2%10000000 == 0 && nw2 > 0)) cout << " *********** " << nw1 << " " << nw2 << " " << np1 << " " << np2 << " " << endl;
      for (int ii = 0; ii < NS; ii++) {
        h1_Fit->SetBinContent(ii + 1, COV1[ii]);
      }
      
      for (int ii = 0; ii < NS/Ndiv; ii++) {
        COV1_AVG[ii] = 0.;
        for (int idiv = 0; idiv < Ndiv; idiv++) {
	  COV1_AVG[ii] = COV1_AVG[ii] + COV1[2*ii + idiv]/float(Ndiv);
        }
        COV2_AVG[ii] = 0.;
        for (int idiv = 0; idiv < Ndiv; idiv++) {
	  COV2_AVG[ii] = COV2_AVG[ii] + COV2[2*ii + idiv]/float(Ndiv);
        }
      }

//
// Apply HP filter
//
      if (Flag_HPF_COV) {
        for (int ii = 0; ii < NS; ii++) {
          if (ii == 0) HPF[ii] = COV1[ii];                    
          if (ii > 0) HPF[ii] = alpha_COV*HPF[ii - 1] + alpha_COV*(COV1[ii] - COV1[ii - 1]); 
        }         
        for (int ii = 0; ii < NS; ii++) {
          COV1[ii] = HPF[ii];
        }         
        for (int ii = 0; ii < NS; ii++) {
          if (ii == 0) HPF[ii] = COV2[ii];                    
          if (ii > 0) HPF[ii] = alpha_COV*HPF[ii - 1] + alpha_COV*(COV2[ii] - COV2[ii - 1]); 
        }         
        for (int ii = 0; ii < NS; ii++) {
          COV2[ii] = HPF[ii];
        }
      }         
      
      DIF1 = 0.;
      DIF2 = 0.;
       
      for (int k = 0; k < NS; k++) {
        Flag_Peak1[k] = false;
        Flag_Peak2[k] = false;
        Flag_Peak1_1kHz[k] = false;
        Flag_Peak2_1kHz[k] = false;
        Flag_DIF1[k] = false;
        Flag_DIF2[k] = false;
      } 
      for (int k = 0; k < NS/Ndiv; k++) {
        Flag_DIF1_AVG[k] = false;
        Flag_DIF2_AVG[k] = false;
        Flag_Peak1_1kHz_AVG[k] = false;
        Flag_Peak2_1kHz_AVG[k] = false;
      }

      nBSL1 = 0;
      nBSL2 = 0;
      BSL_M1 = 0.;
      BSL_S1 = 0.;
      BSL_M2 = 0.;
      BSL_S2 = 0.;      
      for (int k = 0; k < NS; k++) {
        if (k > 0) {
          DIF1 = COV1[k] - COV1[k - 1];
          DIF2 = COV2[k] - COV2[k - 1];
          if (DIF1 > RT1_Threshold) {Flag_DIF1[k] = true;}
          else {
            BSL_M1 = BSL_M1 + DIF1;
            BSL_S1 = BSL_S1 + DIF1*DIF1;
            nBSL1++;
          }
          if (DIF2 > RT2_Threshold) {Flag_DIF2[k] = true;}
          else {
            BSL_M2 = BSL_M2 + DIF2;
            BSL_S2 = BSL_S2 + DIF2*DIF2;
            nBSL2++;
          }
        }
        if ((k + 1)%100 == 0) {
          BSL_M1 = BSL_M1/max(1,nBSL1); 
          BSL_S1 = sqrt(BSL_S1/max(1,nBSL1) - BSL_M1*BSL_M1);
          if (nBSL1 > 1) {
            h1_BSLM1->Fill(BSL_M1/Cal1);
            h1_BSLS1->Fill(BSL_S1/Cal1);
          }
          DIF1_bias[k/100] = BSL_M1;
//
          BSL_M2 = BSL_M2/max(1,nBSL2); 
          BSL_S2 = sqrt(BSL_S2/max(1,nBSL2) - BSL_M2*BSL_M2);
          if (nBSL2 > 1) {
            h1_BSLM2->Fill(BSL_M2/Cal2);
            h1_BSLS2->Fill(BSL_S2/Cal2);
          }
          DIF2_bias[k/100] = BSL_M2;
//
          nBSL1 = 0;
          BSL_M1 = 0.;
          BSL_S1 = 0.;
//
          nBSL2 = 0;
          BSL_M2 = 0.;
          BSL_S2 = 0.;
        } 
      }

      DIF1 = 0.;
      DIF2 = 0.;
      nBSL1_AVG = 0;
      nBSL2_AVG = 0;
      BSL_M1_AVG = 0.;
      BSL_S1_AVG = 0.;
      BSL_M2_AVG = 0.;
      BSL_S2_AVG = 0.;
      for (int k = 0; k < NS/Ndiv; k++) {
        if (k > 0) {
          DIF1 = COV1_AVG[k] - COV1_AVG[k - 1];
          if (DIF1 > RT1_Threshold) {Flag_DIF1_AVG[k] = true;}
          else {
            BSL_M1_AVG = BSL_M1_AVG + DIF1;
            BSL_S1_AVG = BSL_S1_AVG + DIF1*DIF1;
            nBSL1_AVG++;
          }
        }
        if (k > 0) {
          DIF2 = COV2_AVG[k] - COV2_AVG[k - 1];
          if (DIF2 > RT2_Threshold) {Flag_DIF2_AVG[k] = true;}
          else {
            BSL_M2_AVG = BSL_M2_AVG + DIF2;
            BSL_S2_AVG = BSL_S2_AVG + DIF2*DIF2;
            nBSL2_AVG++;
          }
        }

        if ((k + 1)%(100/Ndiv) == 0) {
          BSL_M1_AVG = BSL_M1_AVG/max(1,nBSL1_AVG); 
          BSL_S1_AVG = sqrt(BSL_S1_AVG/max(1,nBSL1_AVG) - BSL_M1_AVG*BSL_M1_AVG);
          if (nBSL1_AVG > 1) {
            h1_BSLM1_AVG->Fill(BSL_M1_AVG/Cal1);
            h1_BSLS1_AVG->Fill(BSL_S1_AVG/Cal1);
          }
          DIF1_AVG_bias[k/(100/Ndiv)] = BSL_M1_AVG;
//
          nBSL1_AVG = 0;
          BSL_M1_AVG = 0.;
          BSL_S1_AVG = 0.;
//
          BSL_M2_AVG = BSL_M2_AVG/max(1,nBSL2_AVG); 
          BSL_S2_AVG = sqrt(BSL_S2_AVG/max(1,nBSL2_AVG) - BSL_M2_AVG*BSL_M2_AVG);
          if (nBSL2_AVG > 1) {
            h1_BSLM2_AVG->Fill(BSL_M2_AVG/Cal2);
            h1_BSLS2_AVG->Fill(BSL_S2_AVG/Cal2);
          }
          DIF2_AVG_bias[k/(100/Ndiv)] = BSL_M2_AVG;

          nBSL2_AVG = 0;
          BSL_M2_AVG = 0.;
          BSL_S2_AVG = 0.;
        }
      }
//
      itab1 = 0;
      itab2 = 0;
      for (int k = 1; k < NS; k++) {
        Flag_Peak1[k] = false; 
        Flag_Peak1_1kHz[k] = false; 
        Sum_DIF1[k] = 0.;
        if (Flag_DIF1[k]) {
          j1 = k;
          for (int j = min(k + 1, NS - 1); j < NS; j++) {
            if (j > k && Flag_DIF1[j]) {Flag_DIF1[j] = false;}
	    else {
              j2 = j; 
              break;
            }
          } 

          nbsl = 0;
          bsl_M = 0.;
          bsl_S = 0.;

          Sum1 = 0.;
          for (int kk = j1 - 1; kk < j2 + 1; kk++) {Sum1 = Sum1 + COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100];}
          Sum2 = 0.;
          for (int kk = j1 - 21; kk < j1 - 1; kk++) {
            Sum2 = Sum2 + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100])*float(j2 - j1 + 2)/40.;
            bsl_M = bsl_M + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100]);
            bsl_S = bsl_S + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100])*(COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100]);
            nbsl++;
          }          
          Sum3 = 0.;
          for (int kk = j2 + 2; kk < j2 + 22; kk++) {
            Sum3 = Sum3 + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100])*float(j2 - j1 + 2)/40.;    
            bsl_M = bsl_M + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100]);
            bsl_S = bsl_S + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100])*(COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100]);
            nbsl++;
     	  }
          bsl_M = bsl_M/max(1,nbsl); 
          bsl_S = sqrt(bsl_S/max(1,nbsl) - bsl_M*bsl_M);
          if (nbsl > 1) {
            h1_bsl_M1->Fill(bsl_M/Cal1);
            h1_bsl_S1->Fill(bsl_S/Cal1);
          }          
          Sum_DIF1[k] = Sum1 - Sum2 - Sum3;
          if ((Sum_DIF1[k] > Peak1_Threshold*Cal1) && ((k%100 != 2 && Sum_DIF1[k] <= Peak1_Upper_Threshold*Cal1) || Sum_DIF1[k] > Peak1_Upper_Threshold*Cal1)) Flag_Peak1[k] = true;

          Sum_ASM1 = 0.;
          for (int kk = j1 - 1; kk < j2 + 1; kk++)           {Sum_ASM1 = Sum_ASM1 + ASM1[kk] - ASM1[kk - 1];}
          Sum_ASM2 = 0.;
          for (int kk = j1 - 20 - 1; kk < j1 - 1; kk++) {Sum_ASM2 = Sum_ASM2 + (ASM1[kk] - ASM1[kk - 1])*float(j2 - j1 + 2)/40.;}
          Sum_ASM3 = 0.;
	  for (int kk = j2 + 2; kk < j2 + 20 + 2; kk++) {Sum_ASM3 = Sum_ASM3 + (ASM1[kk] - ASM1[kk - 1])*float(j2 - j1 + 2)/40.;}

          if (Flag_Peak1[k]) h2_ASM1_Amplitude->Fill(Sum_DIF1[k],Sum_ASM1 - Sum_ASM2 - Sum_ASM3);  
	 	
	}
      }

      for (int k = 1; k < NS/Ndiv; k++) {
       Flag_Peak1_AVG[k] = false; 
       Flag_Peak1_1kHz_AVG[k] = false; 
       Sum_DIF1_AVG[k] = 0.;
        if (Flag_DIF1_AVG[k]) {
          j1 = k;
          for (int j = min(k + 1, NS/Ndiv - 1); j < NS/Ndiv; j++) {
            if (j > k && Flag_DIF1_AVG[j]) {Flag_DIF1_AVG[j] = false;}
	    else {
              j2 = j; 
              break;
            }
          } 
          Sum1 = 0.;
          for (int kk = j1 - 1; kk < j2 + 1; kk++) {Sum1 = Sum1 + COV1_AVG[kk] - COV1_AVG[kk - 1] - DIF1_AVG_bias[kk/(100/Ndiv)];}
          Sum2 = 0.;
          for (int kk = j1 - 20/Ndiv - 1; kk < j1 - 1; kk++) {Sum2 = Sum2 + (COV1_AVG[kk] - COV1_AVG[kk - 1] - DIF1_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}
          Sum3 = 0.;
	  for (int kk = j2 + 2; kk < j2 + 20/Ndiv + 2; kk++) {Sum3 = Sum3 + (COV1_AVG[kk] - COV1_AVG[kk - 1] - DIF1_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}   
          Sum_DIF1_AVG[k] = Sum1 - Sum2 - Sum3;
          if (Sum_DIF1_AVG[k] > Peak1_Threshold*Cal1) {
            if ((k%50 != Tbin_1kHz_1 && Sum_DIF1_AVG[k] <= Peak1_Upper_Threshold*Cal1) || Sum_DIF1_AVG[k] > Peak1_Upper_Threshold*Cal1) Flag_Peak1_AVG[k] = true;
            if (k%50 == Tbin_1kHz_1 && Sum_DIF1_AVG[k] <= Peak1_Upper_Threshold*Cal1) Flag_Peak1_1kHz_AVG[k] = true;
          }
	}
      }

      for (int k = 1; k < NS; k++) {
        Flag_Peak2[k] = false; 
        Flag_Peak2_1kHz[k] = false; 
        Sum_DIF2[k] = 0.;
        if (Flag_DIF2[k]) {
          j1 = k;
          for (int j = min(k + 1, NS - 1); j < NS; j++) {
            if (j > k && Flag_DIF2[j]) {Flag_DIF2[j] = false;}
	    else {
              j2 = j; 
              break;
            }
          }

          nbsl = 0;
          bsl_M = 0.;
          bsl_S = 0.;

          Sum1 = 0.;
          for (int kk = j1 - 1; kk < j2 + 1; kk++) {Sum1 = Sum1 + COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100];}
          Sum2 = 0.;
          for (int kk = j1 - 21; kk < j1 - 1; kk++) {
            Sum2 = Sum2 + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100])*float(j2 - j1 + 2)/40.;
            bsl_M = bsl_M + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100]);
            bsl_S = bsl_S + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100])*(COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100]);
            nbsl++;
          }          
          Sum3 = 0.;
          for (int kk = j2 + 2; kk < j2 + 22; kk++) {
            Sum3 = Sum3 + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100])*float(j2 - j1 + 2)/40.;    
            bsl_M = bsl_M + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100]);
            bsl_S = bsl_S + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100])*(COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100]);
            nbsl++;
     	  }

          bsl_M = bsl_M/max(1,nbsl); 
          bsl_S = sqrt(bsl_S/max(1,nbsl) - bsl_M*bsl_M);
          if (nbsl > 1) {
            h1_bsl_M2->Fill(bsl_M/Cal2);
            h1_bsl_S2->Fill(bsl_S/Cal2);
          }          
          Sum_DIF2[k] = Sum1 - Sum2 - Sum3;
          if ((Sum_DIF2[k] > Peak2_Threshold*Cal2) && ((k%100 != 1 && Sum_DIF2[k] <= Peak2_Upper_Threshold*Cal2) || Sum_DIF2[k] > Peak2_Upper_Threshold*Cal2)) Flag_Peak2[k] = true;

          Sum_ASM1 = 0.;
          for (int kk = j1 - 1; kk < j2 + 1; kk++)           {Sum_ASM1 = Sum_ASM1 + ASM2[kk] - ASM2[kk - 1];}
          Sum_ASM2 = 0.;
          for (int kk = j1 - 20 - 1; kk < j1 - 1; kk++) {Sum_ASM2 = Sum_ASM2 + (ASM2[kk] - ASM2[kk - 1])*float(j2 - j1 + 2)/40.;}
          Sum_ASM3 = 0.;
	  for (int kk = j2 + 2; kk < j2 + 20 + 2; kk++) {Sum_ASM3 = Sum_ASM3 + (ASM2[kk] - ASM2[kk - 1])*float(j2 - j1 + 2)/40.;}

          if (Flag_Peak2[k]) h2_ASM2_Amplitude->Fill(Sum_DIF2[k],Sum_ASM1 - Sum_ASM2 - Sum_ASM3);
                    
       	}
      }

      for (int k = 1; k < NS/Ndiv; k++) {
        Flag_Peak2_AVG[k] = false; 
        Flag_Peak2_1kHz_AVG[k] = false; 
        Sum_DIF2_AVG[k] = 0.;
        if (Flag_DIF2_AVG[k]) {
          j1 = k;
          for (int j = min(k + 1, NS/Ndiv - 1); j < NS/Ndiv; j++) {
            if (j > k && Flag_DIF2_AVG[j]) {Flag_DIF2_AVG[j] = false;}
	    else {
              j2 = j; 
              break;
            }
          } 
          Sum1 = 0.;
          for (int kk = j1 - 1; kk < j2 + 1; kk++) {Sum1 = Sum1 + COV2_AVG[kk] - COV2_AVG[kk - 1] - DIF2_AVG_bias[kk/(100/Ndiv)];}
          Sum2 = 0.;
          for (int kk = j1 - 20/Ndiv - 1; kk < j1 - 1; kk++) {Sum2 = Sum2 + (COV2_AVG[kk] - COV2_AVG[kk - 1] - DIF2_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}
          Sum3 = 0.;
	  for (int kk = j2 + 2; kk < j2 + 20/Ndiv + 2; kk++) {Sum3 = Sum3 + (COV2_AVG[kk] - COV2_AVG[kk - 1] - DIF2_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}   
          Sum_DIF2_AVG[k] = Sum1 - Sum2 - Sum3;
          if (Sum_DIF2_AVG[k] > Peak2_Threshold*Cal2) {
            if ((k%50 != Tbin_1kHz_2 && Sum_DIF2_AVG[k] <= Peak2_Upper_Threshold*Cal2) || Sum_DIF2_AVG[k] > Peak2_Upper_Threshold*Cal2) Flag_Peak2_AVG[k] = true;
            if (k%50 == Tbin_1kHz_2 && Sum_DIF2_AVG[k] <= Peak2_Upper_Threshold*Cal2) Flag_Peak2_1kHz_AVG[k] = true;
          }
	}
      }
          
      for (int k = 0; k < NS; k++) {
        kt = iw*1000 + k/100;
        if (Flag_Peak1[k]) {
          if (nt1 < 10000) {
            TPeak1[nt1] = float(iw)*NS + k;
            nt1++;
          }
          h1_Peak1->Fill(Sum_DIF1[k]);
          h1_Peak1_100mV->Fill(Sum_DIF1[k]);
          h1_Peak1_MeV->Fill(Sum_DIF1[k]/Cal1/1000.);
          h1_Peak1_4MeV->Fill(Sum_DIF1[k]/Cal1/1000.);
          T_COV1[ip1_ms] = kt;
          //T_COV1[ip1_ms] = float(iw*1000) + float(k)/100.;
          //cout << " kt " << kt << " " << float(iw*1000) + float(k)/100. << " " << endl;
          P_COV1[ip1_ms] = Sum_DIF1[k]/Cal1;
          np1++;
          ip1_ms++;
        }      
        if (Flag_Peak2[k]) {
          if (nt2 < 10000) {
            TPeak2[nt2] = float(iw)*NS + k;
            nt2++;
          }
          h1_Peak2->Fill(Sum_DIF2[k]);
          h1_Peak2_100mV->Fill(Sum_DIF2[k]);
          h1_Peak2_MeV->Fill(Sum_DIF2[k]/Cal2/1000.);
          h1_Peak2_4MeV->Fill(Sum_DIF2[k]/Cal2/1000.);
          T_COV2[ip2_ms] = kt;
          //T_COV2[ip2_ms] = float(iw*1000) + float(k)/100.;
          P_COV2[ip2_ms] = Sum_DIF2[k]/Cal2;
          np2++;
          ip2_ms++;
        }      
        if (Flag_Peak1[k] && Flag_Peak2[k]) {
          h1_Peak1_C_keV->Fill(Sum_DIF1[k]/Cal1);
          h1_Peak2_C_keV->Fill(Sum_DIF2[k]/Cal2);
          if (Sum_DIF1[k]/Cal1 <= 250.) h1_Peak1_C_250keV->Fill(Sum_DIF1[k]/Cal1);
          if (Sum_DIF2[k]/Cal2 <= 250.) h1_Peak2_C_250keV->Fill(Sum_DIF2[k]/Cal2);
          h1_Peak1_C_MeV->Fill(Sum_DIF1[k]/Cal1/1000.);
          h1_Peak2_C_MeV->Fill(Sum_DIF2[k]/Cal2/1000.);
          T_COV1and2[ip1and2_ms] = kt;
          ip1and2_ms++;
        }
      }

      for (int k = 0; k < NS/Ndiv; k++) {
        if (Flag_Peak1_AVG[k]) {h1_Peak1_100mV_AVG->Fill(Sum_DIF1_AVG[k]);}      
        if (Flag_Peak2_AVG[k]) {h1_Peak2_100mV_AVG->Fill(Sum_DIF2_AVG[k]);}      
        if (Flag_Peak1_AVG[k]) {h1_Peak1_4MeV_AVG->Fill(Sum_DIF1_AVG[k]/Cal1/1000.);}      
        if (Flag_Peak2_AVG[k]) {h1_Peak2_4MeV_AVG->Fill(Sum_DIF2_AVG[k]/Cal2/1000.);}      
        if (Flag_Peak1_1kHz_AVG[k]) {h1_Peak1_100mV_1kHz_AVG->Fill(Sum_DIF1_AVG[k]);}      
        if (Flag_Peak2_1kHz_AVG[k]) {h1_Peak2_100mV_1kHz_AVG->Fill(Sum_DIF2_AVG[k]);}      
      }

      iw++;

      if (iw == iw_sel) {
        for (int ii = 0; ii < NS; ii++) {
          h1_COV2_Test->SetBinContent(ii + 1, COV2[ii]);
        }

        TH1 *hm2 = 0;
        TVirtualFFT::SetTransform(0);
        hm2 = h1_COV2_Test->FFT(hm2, "MAG");
//Look at the DC component and the Nyquist harmonic:
        Int_t n = 100000;
        Double_t re, im;
//That's the way to get the current transform object:
        TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
//Use the following method to get the full output:
        Double_t *re_full = new Double_t[n];
        Double_t *im_full = new Double_t[n];
        for (int iii = 0; iii < n; iii++) {
          re_full[iii] = 0.;
          im_full[iii] = 0.;
        }
        fft->GetPointsComplex(re_full,im_full);
    
        for (int iii = 0; iii < 100000; iii++) {
          re_full[iii] = re_full[iii]/100000.;
          im_full[iii] = im_full[iii]/100000.;
        }        
        for (int iii = 0; iii < 21; iii++) {
          re_full[(2*iii + 1)*1000] = 0.;
          im_full[(2*iii + 1)*1000] = 0.;
        }

//Now let's make a backward transform:
        TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
        fft_back->SetPointsComplex(re_full,im_full);
        fft_back->Transform();
        TH1 *hb = 0;
//Let's look at the output
        hb = TH1::TransformHisto(fft_back,hb,"Re");

        for (int iii = 0; iii < 100000; iii++) {
          //cout << " iii "  << COV2[iii] << " " << hb->GetBinContent(iii + 1) << " " << endl;
          //COV2[iii] = hb->GetBinContent(iii + 1);
        }

        for (int ii = 0; ii < NS; ii++) {
          h1_Raw1->SetBinContent(ii + 1, Raw1[ii]);
          h1_Raw2->SetBinContent(ii + 1, Raw2[ii]);
          h1_Raw3->SetBinContent(ii + 1, Raw3[ii]);
          h1_Raw4->SetBinContent(ii + 1, Raw4[ii]);
          h1_COV1->SetBinContent(ii + 1, COV1[ii]);
          h1_COV2->SetBinContent(ii + 1, COV2[ii]);
          h1_ASM1->SetBinContent(ii + 1, ASM1[ii]);
          h1_ASM2->SetBinContent(ii + 1, ASM2[ii]);
          if (ii > 0) {  
            h1_DIF1->SetBinContent(ii + 1, COV1[ii] - COV1[ii - 1]);
            h1_DIF2->SetBinContent(ii + 1, COV2[ii] - COV2[ii - 1]);
            h1_DIF_ASM1->SetBinContent(ii + 1, ASM1[ii] - ASM1[ii - 1]);
            h1_DIF_ASM2->SetBinContent(ii + 1, ASM2[ii] - ASM2[ii - 1]);
          }
          h1_Raw5->SetBinContent(ii + 1, Raw5[ii]);
        }
        for (int ii = 0; ii < NS/Ndiv; ii++) {
          h1_COV1_AVG->SetBinContent(ii + 1, COV1_AVG[ii]);
          if (ii > 0) h1_DIF1_AVG->SetBinContent(ii + 1, COV1_AVG[ii] - COV1_AVG[ii - 1]);         
          h1_COV2_AVG->SetBinContent(ii + 1, COV2_AVG[ii]);
          if (ii > 0) h1_DIF2_AVG->SetBinContent(ii + 1, COV2_AVG[ii] - COV2_AVG[ii - 1]);         
        }


        Int_t nBSL2T = 0;
        Double_t BSL_M2T = 0.;
        Double_t BSL_S2T = 0.;      
        for (int kk = 0; kk < NS; kk++) {
          if (kk > 0) {
            DIF2 = COV2[kk] - COV2[kk - 1];
            if (DIF2 <= RT2_Threshold) {
              BSL_M2T = BSL_M2T + DIF2;
              BSL_S2T = BSL_S2T + DIF2*DIF2;
              nBSL2T++;
            }
          }
        } 
        BSL_M2T = BSL_M2T/max(1,nBSL2T);
        BSL_S2T = sqrt(BSL_S2T/max(1,nBSL2T) - BSL_M2T*BSL_M2T);
        //h1_BSLM2T->Fill(BSL_M2T);
        //h1_BSLS2T->Fill(BSL_S2T);
        //cout << " BSL2 " << BSL_M2T << " " << BSL_S2T << " " << endl;
      }
      
      if (iw_ms == NS_NTD) {
        ip_COV1 = ip1_ms;
        ip_COV2 = ip2_ms;
        ip_COV1and2 = ip1and2_ms;       
        iw_ms = 0; 
        ip1_ms = 0;
        ip2_ms = 0;
        ip1and2_ms = 0;
      }
    }

    if (i_NTD == NS_NTD) {
      i_NTD = 0;
      j_NTD = 0; 
      k_min = 0;
      k_max = 0;
      ip_NTD = 0;
//
// Apply HP filter
//
      for (int ii = 0; ii < NS_NTD/100; ii++) {
        if (ii == 0) HPF[ii] = NTD[ii];                    
        if (ii > 0) HPF[ii] = alpha_NTD*HPF[ii - 1] + alpha_NTD*(NTD[ii] - NTD[ii - 1]); 
      }         
      for (int ii = 0; ii < NS_NTD/100; ii++) {
        NTD[ii] = HPF[ii];
      }         

      nBSL_NTD = 0;
      BSL_NTD_M = 0.;
      BSL_NTD_S = 0.;

      for (int k = 1; k < NS_NTD/100; k++) {
        DIF_NTD = NTD[k - 1] - NTD[k];
        if (abs(DIF_NTD) < 0.003) {
          BSL_NTD_M = BSL_NTD_M + DIF_NTD;
          BSL_NTD_S = BSL_NTD_S + DIF_NTD*DIF_NTD;
          nBSL_NTD++;
        }
      }
      BSL_NTD_M = BSL_NTD_M/max(1,nBSL_NTD); 
      BSL_NTD_S = sqrt(BSL_NTD_S/max(1,nBSL_NTD) - BSL_NTD_M*BSL_NTD_M);
      h1_BSL_NTD_M->Fill(BSL_NTD_M/CalNTD);
      h1_BSL_NTD_S->Fill(BSL_NTD_S/CalNTD);

      for (int k = 1; k < NS_NTD/100; k++) {
        if (k > k_max) { 
          DIF_NTD = NTD[k - 1] - NTD[k];
          k_min = k - 1;
          k_max = k - 1;
          if (DIF_NTD < 0.0025) {
            BSL_NTD_M = BSL_NTD_M + DIF_NTD;
            BSL_NTD_S = BSL_NTD_S + DIF_NTD*DIF_NTD;
            nBSL_NTD++;
          }
 
          if (DIF_NTD > 0.) {
            k_max++;
            for (int l = k + 1; l < NS_NTD/100; l++) {
              if (NTD[l - 1] - NTD[l] > 0.) {k_max++;}
              else {break;}              
            }
	    if (NTD[k_min] - NTD[k_max] > 1.5 && k_max - k_min < 5) {
	      //	      if (iw_NTD == 17) {
	      if (iw_NTD == 23) {
                iw_NTD_sel = iw_NTD;
                Flag_NTD = false;
	        //cout << " " << iw_NTD << " " << NTD[k_min] - NTD[k_max] << " " << k_max - k_min << " " << k_min << " " << k_max << " " << NTD[k_min] << " " << NTD[k_max] << " " << endl;
              }
              if (Flag_NTD) {
                iw_NTD_sel = iw_NTD;
                Flag_NTD = false;
	        //cout << " " << iw_NTD << " " << NTD[k_min] - NTD[k_max] << " " << k_max - k_min << " " << k_min << " " << k_max << " " << NTD[k_min] << " " << NTD[k_max] << " " << endl;
              }
            }
            if (k_max - k_min > 5 && k_max < 99990 && NTD[k_min] - NTD[k_max] > NTD_Threshold) {
              if (k_min > 99970) cout << " WARNING !!! " << k_min << " " << k_max << " " << iw_NTD << " " << NTD[k_min] - NTD[k_max] << " " << endl; 
              h1_T_NTD->Fill(float(k_max));
              h1_NTD_Amplitude->Fill(NTD[k_min] - NTD[k_max]);
              //float out_NTD = NTD[k_min] - NTD[k_max];
              //wf5.write((char *) &out_NTD, sizeof out_NTD);
            }
            if (NTD[k_min] - NTD[k_max] > NTD_Threshold) h2_Width_Amplitude->Fill(float(k_max - k_min),NTD[k_min] - NTD[k_max]);
            if (k_max - k_min > 5 && k_max < 99990 && NTD[k_min] - NTD[k_max] > NTD_Threshold) {
              //T_NTD[ip_NTD] = iw_NTD*NS + k_min;
              P_NTD[ip_NTD] = NTD[k_min] - NTD[k_max];

              Int_t ll0 = k_min;
	      Double_t sl_max = -99999.;
              for (int ll = k_min + 1; ll < k_max + 1; ll++) {
                if (NTD[ll - 1] - NTD[ll] > sl_max) {
                  ll0 = ll - 1;
                  sl_max = NTD[ll - 1] - NTD[ll];
                }                  
              }
              Double_t k0 =  ll0 - (NTD[k_min] - NTD[ll0])/sl_max;
              T_NTD[ip_NTD] = iw_NTD*NS + k0 - 1.6;
              ip_NTD++;
              np_NTD++;
            }
          }
        }   
      }

      Double_t alpha = 1000./(1000. + 1.);
      Double_t beta = 1./(1000. + 1.);

      if (iw_NTD == iw_NTD_sel && Flag_NTD == false) {
        for (int ii = 0; ii < NS_NTD/100; ii++) {
          h1_NTD->SetBinContent(ii + 1, NTD[ii]);
	  if (ii > 0) h1_DIF_NTD->SetBinContent(ii + 1, -(NTD[ii] - NTD[ii - 1]));
        }          
        for (int ii = 0; ii < NS_NTD/100; ii++) {
          if (ii == 0) HPF[ii] = NTD[ii];          
          if (ii == 0) LPF[ii] = NTD[ii];          
	  if (ii > 0) HPF[ii] = alpha*HPF[ii - 1] + alpha*(NTD[ii] - NTD[ii - 1]);
	  if (ii > 0) LPF[ii] = beta*NTD[ii] + (1. - beta)*(LPF[ii - 1]);
          h1_HPF->SetBinContent(ii + 1, HPF[ii]);  
          h1_LPF->SetBinContent(ii + 1, LPF[ii]);  
        }         
      }

      //cout << " ip_NTD = " << ip_NTD << " " << endl;

      for (int ip = 0; ip < ip_NTD; ip++) {
        Bool_t Flag_COV1 = false;
        Bool_t Flag_COV2 = false;
        Bool_t Flag_noCOV1 = false;
        Bool_t Flag_noCOV2 = false;
        Bool_t Flag_noCOV1andCOV2 = false;
        Bool_t Flag_COV1andCOV2 = false;
         

          h1_NTD_1_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          for (int jp = 0; jp < ip_COV1; jp++) {
        if (abs(T_NTD[ip] - T_COV1[jp]) <= 50.) {
              if (abs(T_NTD[ip] - T_COV1[jp]) <= 25.) h1_DeltaT_NTD_COV1->Fill(T_NTD[ip] - T_COV1[jp]);
              h2_DeltaT_NTD_COV1_P_NTD->Fill(log10(P_NTD[ip]/CalNTD/1000.),T_NTD[ip] - T_COV1[jp]);
            }
            if ((T_NTD[ip] - T_COV1[jp]) >= -2. && (T_NTD[ip] - T_COV1[jp]) < 4.) {
              Flag_COV1 = true;
              h1_P_COV1->Fill(P_COV1[jp]/1000.);
            }
            if ((T_NTD[ip] - T_COV1[jp] - DTacc) >= -2. && (T_NTD[ip] - T_COV1[jp] - DTacc) < 4.) {
              Flag_noCOV1 = true;
              Flag_noCOV1_M = true;
              h1_P_noCOV1->Fill(P_COV1[jp]/1000.);
              h1_P_noCOV1_DTacc->Fill(P_COV1[jp]/1000.);
            }
        
            if ((T_NTD[ip] - T_COV1[jp] + DTacc) >= -2. && (T_NTD[ip] - T_COV1[jp] + DTacc) < 4.) {
              Flag_noCOV1 = true;
              Flag_noCOV1_P = true;
              h1_P_noCOV1->Fill(P_COV1[jp]/1000.);
              h1_P_noCOV1_DTacc->Fill(P_COV1[jp]/1000.);
            }
       
          }
          for (int jp = 0; jp < ip_COV2; jp++) {
        if (abs(T_NTD[ip] - T_COV2[jp]) <= 50.) {
              if (abs(T_NTD[ip] - T_COV2[jp]) <= 25.) h1_DeltaT_NTD_COV2->Fill(T_NTD[ip] - T_COV2[jp]);
              h2_DeltaT_NTD_COV2_P_NTD->Fill(log10(P_NTD[ip]/CalNTD/1000.),T_NTD[ip] - T_COV2[jp]);
            }
            if ((T_NTD[ip] - T_COV2[jp]) >= -2. && (T_NTD[ip] - T_COV2[jp]) < 4.) {
              Flag_COV2 = true;
              h1_P_COV2->Fill(P_COV2[jp]/1000.);
            }
        if ((T_NTD[ip] - T_COV2[jp] - DTacc) >= -2. && (T_NTD[ip] - T_COV2[jp] - DTacc) < 4.) {
              Flag_noCOV2 = true;
              Flag_noCOV2_M = true;
              h1_P_noCOV2->Fill(P_COV2[jp]/1000.);
              h1_P_noCOV2_DTacc->Fill(P_COV2[jp]/1000.);
            }
        
        if ((T_NTD[ip] - T_COV2[jp] + DTacc) >= -2. && (T_NTD[ip] - T_COV2[jp] + DTacc) < 4.) {
              Flag_noCOV2 = true;
              Flag_noCOV2_P = true;
              h1_P_noCOV2->Fill(P_COV2[jp]/1000.);
              h1_P_noCOV2_DTacc->Fill(P_COV2[jp]/1000.);
            }
       
          }

          if (Flag_COV1 == false && Flag_COV2 == false) h1_NTD_2_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          if ((Flag_noCOV1_M && !Flag_COV2) || (Flag_noCOV2_M && !Flag_COV1)) h1_NTD_2_MeV_DTacc->Fill(P_NTD[ip]/CalNTD/1000.);
          if ((Flag_noCOV1_P && !Flag_COV2) || (Flag_noCOV2_P && !Flag_COV1)) h1_NTD_2_MeV_DTacc->Fill(P_NTD[ip]/CalNTD/1000.);
          if (Flag_noCOV1 == true) h1_NTD_3_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          if (Flag_noCOV2 == true) h1_NTD_4_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          if (Flag_noCOV1 == true && Flag_COV2 == false) h1_NTD_7_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          if (Flag_noCOV2 == true && Flag_COV1 == false) h1_NTD_8_MeV->Fill(P_NTD[ip]/CalNTD/1000.);

          
          for (int jp1 = 0; jp1 < ip_COV1; jp1++) {
        if ((T_NTD[ip] - T_COV1[jp1]) >= -2. && (T_NTD[ip] - T_COV1[jp1]) < 4.) {
              for (int jp2 = 0; jp2 < ip_COV2; jp2++) {
                if ((T_NTD[ip] - T_COV2[jp2]) >= -2. && (T_NTD[ip] - T_COV2[jp2]) < 4.) {
                  Flag_COV1andCOV2 = true;
                }
              }
            }
          }
          if (Flag_COV1andCOV2 == true) h1_NTD_6_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
  //
        }
        iw_NTD++;
      }
    }

    for (int k1 = 0; k1 < nt1; k1++) {
      for (int k2 = 0; k2 < nt2; k2++) {
        Delta_T = TPeak2[k2] - TPeak1[k1];
        if (abs(Delta_T) < 100.) h1_DeltaTCOV->Fill(Delta_T);
      }
    }

  Double_t cc;
  for (int i = 0; i < 10000; i++) {
    cc = h1_Peak1->GetBinContent(i + 1);
    h1_Peak1_keV->SetBinContent(i + 1, cc);
    cc = h1_Peak2->GetBinContent(i + 1);
    h1_Peak2_keV->SetBinContent(i + 1, cc);
  }
  for (int i = 0; i < 300; i++) {
    cc = h1_Peak1_100mV->GetBinContent(i + 1);
    h1_Peak1_100keV->SetBinContent(i + 1, cc);
    cc = h1_Peak2_100mV->GetBinContent(i + 1);
    h1_Peak2_100keV->SetBinContent(i + 1, cc);
    cc = h1_Peak1_100mV_AVG->GetBinContent(i + 1);
    h1_Peak1_100keV_AVG->SetBinContent(i + 1, cc);
    cc = h1_Peak2_100mV_AVG->GetBinContent(i + 1);
    h1_Peak2_100keV_AVG->SetBinContent(i + 1, cc);
    cc = h1_Peak1_100mV_1kHz_AVG->GetBinContent(i + 1);
    h1_Peak1_100keV_1kHz_AVG->SetBinContent(i + 1, cc);
    cc = h1_Peak2_100mV_1kHz_AVG->GetBinContent(i + 1);
    h1_Peak2_100keV_1kHz_AVG->SetBinContent(i + 1, cc);
  }
  for (int i = 0; i < 1000; i++) {
    cc = h1_NTD_Amplitude->GetBinContent(i + 1);
    h1_NTD_Amplitude_keV->SetBinContent(i + 1, cc);
  }

  cout << " " << nw1 << " " << nw2 << " " << np1 << " " << np2 << " " << np_NTD << " " << iw << " " << iw_ms << " " << iw_NTD << " " << endl;
  cout << " " << nnn << " " << endl;
    
  h1_RunTime->SetBinContent(1, (float(nw1) + float(nw2))/float(NS));

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Write histograms to file                                                // 
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  fHistos.Write(); 
  fHistos.Close();
 
  cout << " " << endl;
  cout << " Output root file closed " << endl;	
  cout << " " << endl;

}
