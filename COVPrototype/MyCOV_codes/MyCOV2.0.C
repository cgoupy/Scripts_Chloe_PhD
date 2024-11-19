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
  //Double_t Cal1 =  .0025739;
  //Double_t Cal2 =  .0025739;
  Double_t Cal1 =  .0009533;
  Double_t Cal2 =  .0006513;
  //Double_t Off1 = .0;
  //Double_t Off2 = .0;
  Double_t Off1 = -.001038;
  Double_t Off2 = -.000835;
  //  Double_t CalNTD =  .0002361;
  //  Double_t OffNTD =  .0011361;
  Double_t CalNTD =  .00024025;
  Double_t OffNTD = -.00156369;
  Double_t DTacc = 30.;
  Double_t h1_Alpha_Min = 3500.;
  Double_t h1_Alpha_Max = 5500.;
  Int_t Tbin_1kHz_1 = 1;
  Int_t Tbin_1kHz_2 = 0;
  Int_t kt1, kt2, dkt1_before, dkt1_after, dkt2_before, dkt2_after;

  Int_t N_OVL = 5;
  Double_t E_OVL[5] = { 2.,  2., 2.,2., 2.};
  Double_t RP1[5] = { 0., -1.06058, -1.55123, -1.58302, -1.58355};
  Double_t RP2[5] = { 0.,  1.01586,  1.51078,  1.53515,  1.53509};
  Double_t RP3[5] = { 0., -0.78211, -1.08075, -1.09541, -1.09723};
  Double_t RP4[5] = { 0.,  0.74168,  1.02268,  1.03590,  1.03548};

  cout << " Open root file: " << myfile_name << " " << endl;
  TFile fHistos(myfile_name,"recreate");

  TH1F *h1_RunTime, *h1_Thresholds, *h1_Gains;
  TH1F *h1_Raw1, *h1_Raw2, *h1_Raw3, *h1_Raw4, *h1_Raw5;
  TH1F *h1_COV1, *h1_COV2, *h1_DIF1, *h1_DIF2;
  TH1F *h1_Peak1, *h1_Peak2, *h1_Peak1_MeV, *h1_Peak2_MeV, *h1_Peak1_4MeV, *h1_Peak2_4MeV, *h1_Peak1_4MeV_AVG, *h1_Peak2_4MeV_AVG;
  TH1F *h1_BSLM1, *h1_BSLM2, *h1_BSLS1, *h1_BSLS2, *h1_BSL_NTD_M, *h1_BSL_NTD_S;
  TH1F *h1_bsl_M1, *h1_bsl_S1, *h1_bsl_M2, *h1_bsl_S2;
  TH1F *h1_DeltaTCOV;
  TH1F *h1_Peak1_keV, *h1_Peak2_keV, *h1_Peak1_C_keV, *h1_Peak2_C_keV, *h1_Peak1_C_250keV, *h1_Peak2_C_250keV, *h1_Peak1_C_MeV, *h1_Peak2_C_MeV;
  TH1F *h1_Peak1_100mV, *h1_Peak2_100mV, *h1_Peak1_100keV, *h1_Peak2_100keV;
  TH1F *h1_DIF_NTD, *h1_NTD, *h1_NTD_Amplitude, *h1_NTD_Amplitude_keV;
  TH2F *h2_Width_Amplitude, *h2_NTD_T1_T2, *h2_DeltaT_NTD_COV1_P_NTD, *h2_DeltaT_NTD_COV2_P_NTD;
  TH1F *h1_NTD_1_MeV, *h1_NTD_2_MeV, *h1_NTD_3_MeV, *h1_NTD_4_MeV, *h1_NTD_5_MeV, *h1_NTD_6_MeV, *h1_NTD_7_MeV, *h1_NTD_8_MeV, *h1_NTD_9_MeV, *h1_NTD_10_MeV, *h1_NTD_2_MeV_DTacc;
  TH1F *h1_Raw1_Overlay, *h1_Raw2_Overlay, *h1_Raw3_Overlay, *h1_Raw4_Overlay; 
  TH1F *h1_COV1_OVL, *h1_COV2_OVL, *h1_E_OVL;
  TH1F *h1_T_NTD;
  TH1F *h1_COV1_AVG, *h1_DIF1_AVG, *h1_Peak1_100mV_AVG, *h1_Peak1_100keV_AVG, *h1_Peak1_100mV_1kHz_AVG, *h1_Peak1_100keV_1kHz_AVG;
  TH1F *h1_COV2_AVG, *h1_DIF2_AVG, *h1_Peak2_100mV_AVG, *h1_Peak2_100keV_AVG, *h1_Peak2_100mV_1kHz_AVG, *h1_Peak2_100keV_1kHz_AVG;
  TH1F *h1_BSLM1_AVG, *h1_BSLS1_AVG, *h1_BSLM2_AVG, *h1_BSLS2_AVG;
  TH1F *h1_P_COV1, *h1_P_COV2, *h1_P_noCOV1, *h1_P_noCOV2, *h1_P_noCOV1_DTacc, *h1_P_noCOV2_DTacc;
  TH1F *h1_Peak1_200keV_AVG, *h1_Peak2_200keV_AVG;

  Int_t Ndiv = 2;
  Int_t inw = 0;

  h1_RunTime = new TH1F("h1_RunTime"," ",1,0.,1.);
  h1_Thresholds = new TH1F("h1_Thresholds"," ",3,0.,3.);
  h1_Gains = new TH1F("h1_Gains"," ",3,0.,3.);
  h1_E_OVL = new TH1F("h1_E_OVL"," ",200,0.,200.);  
  TH1F *h1_DeltaT_Closest = new TH1F("h1_DeltaT_Closest"," ",500,-50., 50.);
  TH2F *h2_E_NTD_DeltaT_Closest = new TH2F("h2_E_NTD_DeltaT_Closest "," ", 50, 0., 2000., 50,-5., 5.);
  TH1F *h1_BSLM_RAW = new TH1F("h1_BSLM_RAW"," ", 100, -.25, .25);
  TH2F *h2_Time_BSLM_RAW = new TH2F("h2_Time_BSLM_RAW"," ", 60, 0.,60.,50, -.25, 0.25);
  TH2F *h2_BSLM_E_Alpha = new TH2F("h2_BSLM_E_Alpha"," ", 50, -.25, 0.25, 75,5000.,5150.); 
  h1_Raw1 = new TH1F("h1_Raw1"," ",100000,0.,100000.);
  h1_Raw2 = new TH1F("h1_Raw2"," ",100000,0.,100000.);
  h1_Raw3 = new TH1F("h1_Raw3"," ",100000,0.,100000.);
  h1_Raw4 = new TH1F("h1_Raw4"," ",100000,0.,100000.);
  h1_Raw5 = new TH1F("h1_Raw5"," ",100000,0.,100000.);
  h1_COV1 = new TH1F("h1_COV1"," ",100000,0.,100000.);
  h1_COV2 = new TH1F("h1_COV2"," ",100000,0.,100000.);
  h1_DIF1 = new TH1F("h1_DIF1"," ",100000,0.,100000.);
  h1_DIF2 = new TH1F("h1_DIF2"," ",100000,0.,100000.);
  h1_COV1_AVG = new TH1F("h1_COV1_AVG"," ",100000/Ndiv,0.,100000.);
  h1_COV2_AVG = new TH1F("h1_COV2_AVG"," ",100000/Ndiv,0.,100000.);
  h1_DIF1_AVG = new TH1F("h1_DIF1_AVG"," ",100000/Ndiv,0.,100000.);
  h1_DIF2_AVG = new TH1F("h1_DIF2_AVG"," ",100000/Ndiv,0.,100000.);
  h1_DIF_NTD = new TH1F("h1_DIF_NTD"," ",100000,0.,100000.);
  h1_Raw1_Overlay = new TH1F("h1_Raw1_Overlay"," ",5000,0.,5000.);
  h1_Raw2_Overlay = new TH1F("h1_Raw2_Overlay"," ",5000,0.,5000.);
  h1_Raw3_Overlay = new TH1F("h1_Raw3_Overlay"," ",5000,0.,5000.);
  h1_Raw4_Overlay = new TH1F("h1_Raw4_Overlay"," ",5000,0.,5000.);
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
  TH1F *h1_Peak1_MeV_AVG = new TH1F("h1_Peak1_MeV_AVG"," ",2000,0.,20.);
  TH1F *h1_Peak2_MeV_AVG = new TH1F("h1_Peak2_MeV_AVG"," ",2000,0.,20.);
  TH1F *h1_Peak1_C_MeV_AVG = new TH1F("h1_Peak1_C_MeV_AVG"," ",2000,0.,20.);
  TH1F *h1_Peak2_C_MeV_AVG = new TH1F("h1_Peak2_C_MeV_AVG"," ",2000,0.,20.);
  h1_Peak1_4MeV_AVG = new TH1F("h1_Peak1_4MeV_AVG"," ",400,0.,4.);
  h1_Peak2_4MeV_AVG = new TH1F("h1_Peak2_4MeV_AVG"," ",400,0.,4.);
  h1_Peak1_200keV_AVG = new TH1F("h1_Peak1_200keV_AVG"," ",400,0.,200.);
  h1_Peak2_200keV_AVG = new TH1F("h1_Peak2_200keV_AVG"," ",400,0.,200.);
  h1_COV1_OVL = new TH1F("h1_COV1_OVL"," ",400,0.,200.);
  h1_COV2_OVL = new TH1F("h1_COV2_OVL"," ",400,0.,200.);
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
  TH2F *h2_BSLM_E_NTD = new TH2F("h2_BSLM_E_NTD"," ",250, -2.50, 2.50,500,3000.,8000.);
  h1_DeltaTCOV = new TH1F("h1_DeltaTCOV"," ",  100, -50.,50.);
  h1_NTD = new TH1F("h1_NTD"," ",100000,0.,100000.);
  h1_NTD_Amplitude      = new TH1F("h1_NTD_Amplitude"    ," ", 500,      0.,        2.);
  h1_NTD_Amplitude_keV  = new TH1F("h1_NTD_Amplitude_keV"," ", 500,  OffNTD, 2./CalNTD);
  h1_T_NTD = new TH1F("h1_T_NTD"," ",2000,0.,100000.);
  TH1F *h1_DeltaT_NTD_COV1 = new TH1F("h1_DeltaT_NTD_COV1"," ", 250, -25.,25.);
  TH1F *h1_DeltaT_NTD_COV2 = new TH1F("h1_DeltaT_NTD_COV2"," ", 250, -25.,25.);
  TH1F *h1_DeltaT_NTD_COV1and2 = new TH1F("h1_DeltaT_NTD_COV1and2"," ", 250, -25.,25.);
  h1_NTD_1_MeV  = new TH1F("h1_NTD_1_MeV"," ",2000,  0., 20.);
  h1_NTD_2_MeV  = new TH1F("h1_NTD_2_MeV"," ",2000,  0., 20.);
  h1_NTD_3_MeV  = new TH1F("h1_NTD_3_MeV"," ",2000,  0., 20.);
  h1_NTD_4_MeV  = new TH1F("h1_NTD_4_MeV"," ",2000,  0., 20.);
  h1_NTD_5_MeV  = new TH1F("h1_NTD_5_MeV"," ",2000,  0., 20.);
  h1_NTD_6_MeV  = new TH1F("h1_NTD_6_MeV"," ",2000,  0., 20.);
  h1_NTD_7_MeV  = new TH1F("h1_NTD_7_MeV"," ",2000,  0., 20.);
  h1_NTD_8_MeV  = new TH1F("h1_NTD_8_MeV"," ",2000,  0., 20.);
  h1_NTD_9_MeV  = new TH1F("h1_NTD_9_MeV"," ",2000,  0., 20.);
  h1_NTD_10_MeV  = new TH1F("h1_NTD_10_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_11_MeV  = new TH1F("h1_NTD_11_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_12_MeV  = new TH1F("h1_NTD_12_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_13_MeV  = new TH1F("h1_NTD_13_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_21_MeV  = new TH1F("h1_NTD_21_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_22_MeV  = new TH1F("h1_NTD_22_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_23_MeV  = new TH1F("h1_NTD_23_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_24_MeV  = new TH1F("h1_NTD_24_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_25_MeV  = new TH1F("h1_NTD_25_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_26_MeV  = new TH1F("h1_NTD_26_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_27_MeV  = new TH1F("h1_NTD_27_MeV"," ",2000,  0., 20.);
  TH1F *h1_NTD_28_MeV  = new TH1F("h1_NTD_28_MeV"," ",2000,  0., 20.);
  h1_NTD_2_MeV_DTacc  = new TH1F("h1_NTD_2_MeV_DTacc"," ",2000,  0., 20.);
  h1_P_COV1 = new TH1F("h1_P_COV1"," ",2000,0.,20.);
  h1_P_COV2 = new TH1F("h1_P_COV2"," ",2000,0.,20.);
  h1_P_noCOV1 = new TH1F("h1_P_noCOV1"," ",2000,0.,20.);
  h1_P_noCOV2 = new TH1F("h1_P_noCOV2"," ",2000,0.,20.);
  h1_P_noCOV1_DTacc = new TH1F("h1_P_noCOV1_DTacc"," ",2000,0.,20.);
  h1_P_noCOV2_DTacc = new TH1F("h1_P_noCOV2_DTacc"," ",2000,0.,20.);
  h2_Width_Amplitude  = new TH2F("h2_Width_Amplitude"," ",50,  0.,50.,1000.,0.,5.0);
  h2_NTD_T1_T2  = new TH2F("h2_NTD_T1_T2"," ",100, -50.,50.,100.,-50.,50.);
  h2_DeltaT_NTD_COV1_P_NTD  = new TH2F("h2_DeltaT_NTD_COV1_P_NTD"," ",100,  -2., 2.,20,-10.,10.);
  h2_DeltaT_NTD_COV2_P_NTD  = new TH2F("h2_DeltaT_NTD_COV2_P_NTD"," ",100,  -2., 2.,20,-10.,10.);
  TH1F *h1_BslM1_AVG[20], *h1_BslS1_AVG[20], *h1_BslM2_AVG[20], *h1_BslS2_AVG[20];
  for (int ik = 0; ik < 20; ik++) {
    h1_BslM1_AVG[ik] = new TH1F(Form("h1_BslM1_AVG_%02d",ik)," ",  200,-5.0,5.0);
    h1_BslS1_AVG[ik] = new TH1F(Form("h1_BslS1_AVG_%02d",ik)," ",  100, 0.0,5.0);
    h1_BslM2_AVG[ik] = new TH1F(Form("h1_BslM2_AVG_%02d",ik)," ",  200,-5.0,5.0);
    h1_BslS2_AVG[ik] = new TH1F(Form("h1_BslS2_AVG_%02d",ik)," ",  100, 0.0,5.0);
  }
  TH1F *h1_Peak1_Am241_AVG[20];
  for (int ik = 0; ik < 20; ik++) {
    h1_Peak1_Am241_AVG[ik] = new TH1F(Form("h1_Peak1_Am241_AVG_%02d",ik)," ",60, 40.,100.);
  }
  TH1F *h1_Alpha[60];
  for (int ik = 0; ik < 60; ik++) {h1_Alpha[ik] = new TH1F(Form("h1_Alpha_%02d",ik)," ",100,h1_Alpha_Min,h1_Alpha_Max);}
  TH1F *h1_DeltaT_Closest_E[20];
  for (int ik = 0; ik < 20; ik++) {h1_DeltaT_Closest_E[ik] = new TH1F(Form("h1_DeltaT_Closest_E_%02d",ik)," ",500,-50.,50.);}
  TH1F *h1_IE[20];
  for (int ik = 0; ik < 20; ik++) {h1_IE[ik] = new TH1F(Form("h1_IE_%02d",ik)," ",300, 0.,15.);}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

  Int_t nov1 = 0;
  Int_t nov2 = 0;
  Bool_t Flag_Little_Endian = true;
 
  float Raw1[100000], Raw2[100000], Raw3[100000], Raw4[100000], Raw5[100000];
  float COV1[100000], COV2[100000], Sum_DIF1[100000], Sum_DIF2[100000];
  float Raw1_Overlay[5000], Raw2_Overlay[5000], Raw3_Overlay[5000], Raw4_Overlay[5000];
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
  for (int k = 0; k < 5000; k++) {
    Raw1_Overlay[k] = 0.;
    Raw2_Overlay[k] = 0.;
    Raw3_Overlay[k] = 0.;
    Raw4_Overlay[k] = 0.;
  }
  float DIF1, DIF2;
  Double_t Sn, Sx, Sy, Sxx, Sxy;
  Double_t aal, bbl, aau, bbu;

  signed short ampl;
  int i;
  unsigned iw, iw_ms, iw_sel, iw_NTD, nw1, nw2;
  unsigned iw_NTD_sel = 0;
  Bool_t Flag_NTD = true;
  Bool_t Flag_Overlay = false;
  Bool_t Flag_Sum3 = true;
  Int_t NS = 100000;
  Int_t FRed = 100;
  Int_t NS_NTD = FRed*NS;
  float TPeak1[10000], TPeak2[10000]; 
  float Delta_T;

  iw_sel = 2;

  cout << " Start MyCOV1 !!! " << endl;

// Read Ge detector data

  ifstream rf1("/sps/lbno/edoardo/Actuator/2021/Logged Data_2021_05_07_00_34_08.Bin1", ios::in | ios::binary);
  if (!rf1) {cout << "Cannot open input file 1!" << endl;}

  ifstream rf2("/sps/lbno/edoardo/Actuator/2021/Logged Data_2021_05_07_00_34_08.Bin2", ios::in | ios::binary);
  if (!rf2) {cout << "Cannot open input file 2!" << endl;}

  ifstream rf3("/sps/lbno/edoardo/Actuator/2021/Logged Data_2021_05_07_00_34_08.Bin3", ios::in | ios::binary);
  if (!rf3) {cout << "Cannot open input file 3!" << endl;}

  ifstream rf4("/sps/lbno/edoardo/Actuator/2021/Logged Data_2021_05_07_00_34_08.Bin4", ios::in | ios::binary);
  if (!rf4) {cout << "Cannot open input file 4!" << endl;}

  ifstream rf5("/sps/lbno/edoardo/Actuator/2021/Logged Data_2021_05_07_00_34_08.NTD", ios::in | ios::binary);
  if (!rf5) {cout << "Cannot open input file 5!" << endl;}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*
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
*/
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
  Int_t nBSL1, nBSL2, nBSL_NTD, nBSL_RAW;
  Double_t BSL_M1, BSL_S1, BSL_M2, BSL_S2, BSL_NTD_M, BSL_NTD_S, BSL_RAW_M, BSL_RAW_S;
  Int_t nBSL1_AVG, nBSL2_AVG;
  Double_t BSL_M1_AVG, BSL_S1_AVG, BSL_M2_AVG, BSL_S2_AVG;
  float NTD_Avg, DIF_NTD, NTD[100000];
  Double_t BSL_NTD_RAW[1000];
 
  Double_t Sum1, Sum2, Sum3;  

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

  Int_t N_Seconds = 100;
  while (!rf1.eof() && !rf2.eof() && !rf3.eof() && !rf4.eof() && !rf5.eof()) {
    //while (!rf1.eof() && !rf2.eof() && !rf3.eof() && !rf4.eof() && !rf5.eof() && nw1 < N_Seconds*NS) {
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
      //COV1[i] = 20.*((1. + .05*float(inw - 10))*Raw2[i] - (1. - .05*float(inw - 10))*Raw1[i])/65536.;
      //COV1[i] = 20.*(-Raw1[i] - Raw1[i])/65536.;     
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
      //COV2[i] = 20.*(-Raw3[i] - Raw3[i])/65536.;
      //COV2[i] = 20.*((1. + .05*float(inw - 10))*Raw4[i] - (1. - .05*float(inw - 10))*Raw3[i])/65536.;
    }
    rf5.read((char *) &rx, sizeof rx);
    if (i < NS) { 
      Raw5[i] = ReverseFloat(rx);
    }

//
// NTD : move from 100kHz to 100kHz/FRed sampling
// 
    if (i_NTD < NS_NTD) {
      if ((i_NTD + 1)%FRed == 0) {
	NTD_Avg = NTD_Avg + Raw5[i]/float(FRed);
        NTD[j_NTD] = NTD_Avg;
        j_NTD++;
        NTD_Avg = 0.;
      } 
      else {NTD_Avg = NTD_Avg + Raw5[i]/float(FRed);}
    } 

    i++;
    iw_ms++;
    i_NTD++;

    if (nw1 < 4000000000) {
      inw = (nw1/NS)%20;
      nw1++;
    }
    else {
      nw2++;
      inw = (nw2/NS)%20;
    }
    
    if (i == NS) {  
      if (Flag_Overlay) {
        Int_t iww = iw%5;
        if (iww < N_OVL) {
          h1_E_OVL->Fill(E_OVL[iww]); 
          Int_t iww1 = 10000 + iw%100; 
          for (int il = 0; il < 20000; il++) {
            if (il < 5) {
              Raw1[il + iww1] = Raw1[il + iww1] + E_OVL[iww]*RP1[il];  
              Raw2[il + iww1] = Raw2[il + iww1] + E_OVL[iww]*RP2[il];
              Raw3[il + iww1] = Raw3[il + iww1] + E_OVL[iww]*RP3[il];  
              Raw4[il + iww1] = Raw4[il + iww1] + E_OVL[iww]*RP4[il];
            }     
            else {     
  	      Raw1[il + iww1] = Raw1[il + iww1] + E_OVL[iww]*RP1[4]*exp(-(il - 5)*.0004808);  
	      Raw2[il + iww1] = Raw2[il + iww1] + E_OVL[iww]*RP2[4]*exp(-(il - 5)*.0004808);
  	      Raw3[il + iww1] = Raw3[il + iww1] + E_OVL[iww]*RP3[4]*exp(-(il - 5)*.0002887);  
	      Raw4[il + iww1] = Raw4[il + iww1] + E_OVL[iww]*RP4[4]*exp(-(il - 5)*.0002887);
            }           
            COV1[il + iww1] = 20.*(Raw2[il + iww1] - Raw1[il + iww1])/65536.;
            COV2[il + iww1] = 20.*(Raw4[il + iww1] - Raw3[il + iww1])/65536.;
          }
        }     
      }

      i = 0;
      kt1 = 0;
      kt2 = 0;
      dkt1_before = 0;
      dkt2_after = 0; 
      dkt2_before = 0;
      dkt2_after = 0; 
      if ((nw1%10000000 == 0 && nw2 == 0) || (nw2%10000000 == 0 && nw2 > 0)) cout << " *********** " << nw1/NS + nw2/NS << " seconds     " << np1 << " " << np2 << " " << np_NTD << " " << endl;
      //      for (int ii = 0; ii < NS; ii++) {
      //h1_Fit->SetBinContent(ii + 1, COV1[ii]);
      //}
      
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
            h1_BslM1_AVG[inw]->Fill(BSL_M1_AVG/Cal1);
            h1_BslS1_AVG[inw]->Fill(BSL_S1_AVG/Cal1);
 
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
            h1_BslM2_AVG[inw]->Fill(BSL_M2_AVG/Cal2);
            h1_BslS2_AVG[inw]->Fill(BSL_S2_AVG/Cal2);
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
          if (Flag_Sum3) {
            for (int kk = j2 + 2; kk < j2 + 22; kk++) {
              Sum3 = Sum3 + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100])*float(j2 - j1 + 2)/40.;    
              bsl_M = bsl_M + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100]);
              bsl_S = bsl_S + (COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100])*(COV1[kk] - COV1[kk - 1] - DIF1_bias[kk/100]);
              nbsl++;
     	    }
          }
          bsl_M = bsl_M/max(1,nbsl); 
          bsl_S = sqrt(bsl_S/max(1,nbsl) - bsl_M*bsl_M);
          if (nbsl > 1) {
            h1_bsl_M1->Fill(bsl_M/Cal1);
            h1_bsl_S1->Fill(bsl_S/Cal1);
          }          
          Sum_DIF1[k] = Sum1 - Sum2 - Sum3;
	  if ((Sum_DIF1[k] > Peak1_Threshold*Cal1) && ((k%100 != 2 && Sum_DIF1[k] <= Peak1_Upper_Threshold*Cal1) || Sum_DIF1[k] > Peak1_Upper_Threshold*Cal1)) Flag_Peak1[k] = true;

          if (Flag_Peak1[k]) {
            Int_t ll0 = j1;
	    Double_t sl_max = -99999.;
            for (int ll = j1 + 1; ll < j2 + 1; ll++) {
              if (COV1[ll - 1] - COV1[ll] > sl_max) {
                ll0 = ll - 1;
                sl_max = COV1[ll - 1] - COV1[ll];
              }                  
            }
            Double_t k0 =  ll0 - (COV1[j1] - COV1[ll0])/sl_max;
            //cout << " T1 " << j1 << " " << COV1[j1] << " " << j2 << " " << COV1[j2] << " " << k0 << " " << endl; 
          }

          //if (Flag_Peak1[k] == true) cout << " dj ok " << j2 - j1 + 1 << " " << Sum_DIF1[k]/Cal1 << " " << endl;

          if (Flag_Overlay && Flag_Peak1[k]) {
            dkt1_after = k - kt1;
            if (dkt1_after > 17500 && dkt1_before > 17500 && Sum_DIF1[kt1]/Cal1 > 150. && Sum_DIF1[kt1]/Cal1 < 2000.) {
              nov1++;
              cout << " Overlay 1 !!! << " << nov1 << " " << iw << " " << kt1 << " " << dkt1_before << " " << dkt1_after << " " << Sum_DIF1[kt1]/Cal1 << " " << Raw1[kt1] << " " << Raw2[kt1] << " " << endl;
              for (int kk = kt1 - 1; kk < kt1 + 5000 - 1; kk++) {
                Int_t kk0 = kt1 - 1;
                Raw1_Overlay[kk - kk0] = Raw1_Overlay[kk - kk0] + Cal1*(Raw1[kk] - Raw1[kt1 - 1])/Sum_DIF1[kt1];
                Raw2_Overlay[kk - kk0] = Raw2_Overlay[kk - kk0] + Cal1*(Raw2[kk] - Raw2[kt1 - 1])/Sum_DIF1[kt1];                 
              }
            }
            kt1 = k;
            dkt1_before = dkt1_after; 
          }   	
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
          if (Flag_Sum3) {
	    for (int kk = j2 + 2; kk < j2 + 20/Ndiv + 2; kk++) {Sum3 = Sum3 + (COV1_AVG[kk] - COV1_AVG[kk - 1] - DIF1_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}
          }    
          Sum_DIF1_AVG[k] = Sum1 - Sum2 - Sum3;
          if (Sum_DIF1_AVG[k] > Peak1_Threshold*Cal1) {
            if ((k%50 != Tbin_1kHz_1 && Sum_DIF1_AVG[k] <= Peak1_Upper_Threshold*Cal1) || Sum_DIF1_AVG[k] > Peak1_Upper_Threshold*Cal1) Flag_Peak1_AVG[k] = true;
            if (k%50 == Tbin_1kHz_1 && Sum_DIF1_AVG[k] <= Peak1_Upper_Threshold*Cal1) Flag_Peak1_1kHz_AVG[k] = true;
          }
          //if (Flag_Peak1_AVG[k] == true)       }

          //if (Flag_Peak1[k] == true) cout << " dj ok " << j2 - j1 + 1 << " " << Sum_DIF1[k]/Cal1 << " " << endl;

          if (Flag_Overlay && Flag_Peak1[k]) {
            dkt1_after = k - kt1;
            if (dkt1_after > 17500 && dkt1_before > 17500 && Sum_DIF1[kt1]/Cal1 > 150. && Sum_DIF1[kt1]/Cal1 < 2000.) {
              nov1++;
              cout << " Overlay 1 !!! << " << nov1 << " " << iw << " " << kt1 << " " << dkt1_before << " " << dkt1_after << " " << Sum_DIF1[kt1]/Cal1 << " " << Raw1[kt1] << " " << Raw2[kt1] << " " << endl;
              for (int kk = kt1 - 1; kk < kt1 + 5000 - 1; kk++) {
                Int_t kk0 = kt1 - 1;
                Raw1_Overlay[kk - kk0] = Raw1_Overlay[kk - kk0] + Cal1*(Raw1[kk] - Raw1[kt1 - 1])/Sum_DIF1[kt1];
                Raw2_Overlay[kk - kk0] = Raw2_Overlay[kk - kk0] + Cal1*(Raw2[kk] - Raw2[kt1 - 1])/Sum_DIF1[kt1];                 
              }
            }
            kt1 = k;
            dkt1_before = dkt1_after; 
          }   	
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
          if (Flag_Sum3) {
	    for (int kk = j2 + 2; kk < j2 + 20/Ndiv + 2; kk++) {Sum3 = Sum3 + (COV1_AVG[kk] - COV1_AVG[kk - 1] - DIF1_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}
          }    
          Sum_DIF1_AVG[k] = Sum1 - Sum2 - Sum3;
          if (Sum_DIF1_AVG[k] > Peak1_Threshold*Cal1) {
            if ((k%50 != Tbin_1kHz_1 && Sum_DIF1_AVG[k] <= Peak1_Upper_Threshold*Cal1) || Sum_DIF1_AVG[k] > Peak1_Upper_Threshold*Cal1) Flag_Peak1_AVG[k] = true;
            if (k%50 == Tbin_1kHz_1 && Sum_DIF1_AVG[k] <= Peak1_Upper_Threshold*Cal1) Flag_Peak1_1kHz_AVG[k] = true;
          }
          //if (Flag_Peak1_AVG[k] == true) cout << " T2 " << j1 << " " << COV1_AVG[j1] << " " << j2 << " " << COV1_AVG[j2] << " " << endl; 
          //if (Flag_Peak1_AVG[k] == true) cout << " dj ok AVG " << j2 - j1 + 1 << " " << Sum_DIF1_AVG[k]/Cal1 << " " << endl;
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
          if (Flag_Sum3) {
            for (int kk = j2 + 2; kk < j2 + 22; kk++) {
              Sum3 = Sum3 + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100])*float(j2 - j1 + 2)/40.;    
              bsl_M = bsl_M + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100]);
              bsl_S = bsl_S + (COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100])*(COV2[kk] - COV2[kk - 1] - DIF2_bias[kk/100]);
              nbsl++;
     	    }
          }
          bsl_M = bsl_M/max(1,nbsl); 
          bsl_S = sqrt(bsl_S/max(1,nbsl) - bsl_M*bsl_M);
          if (nbsl > 1) {
            h1_bsl_M2->Fill(bsl_M/Cal2);
            h1_bsl_S2->Fill(bsl_S/Cal2);
          }          
          Sum_DIF2[k] = Sum1 - Sum2 - Sum3;
	  if ((Sum_DIF2[k] > Peak2_Threshold*Cal2) && ((k%100 != 1 && Sum_DIF2[k] <= Peak2_Upper_Threshold*Cal2) || Sum_DIF2[k] > Peak2_Upper_Threshold*Cal2)) Flag_Peak2[k] = true;
          //if (Flag_Peak2[k] == true) cout << " dj ok " << j2 - j1 + 1 << " " << Sum_DIF2[k]/Cal2 << " " << endl;

          if (Flag_Overlay && Flag_Peak2[k]) {
            dkt2_after = k - kt2;
            if (dkt2_after > 30000 && dkt2_before > 30000 && Sum_DIF2[kt2]/Cal2 > 150. && Sum_DIF2[kt2]/Cal2 < 2000.) {
              nov2++;
              cout << " Overlay 2 !!! << " << nov2 << " " << iw << " " << kt2 << " " << " " << dkt2_before << " " << dkt2_after << " " << Sum_DIF2[kt2]/Cal2 << " " << Raw3[kt2] << " " << Raw4[kt2] << " " << endl;
              for (int kk = kt2 - 1; kk < kt2 + 5000 - 1; kk++) {
                Int_t kk0 = kt2 - 1;
                Raw3_Overlay[kk - kk0] = Raw3_Overlay[kk - kk0] + Cal2*(Raw3[kk] - Raw3[kt2 - 1])/Sum_DIF2[kt2];
                Raw4_Overlay[kk - kk0] = Raw4_Overlay[kk - kk0] + Cal2*(Raw4[kk] - Raw4[kt2 - 1])/Sum_DIF2[kt2];                 
              }
            }
            kt2 = k;
            dkt2_before = dkt2_after; 
          }   	
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
	  if (Flag_Sum3) {
            for (int kk = j2 + 2; kk < j2 + 20/Ndiv + 2; kk++) {Sum3 = Sum3 + (COV2_AVG[kk] - COV2_AVG[kk - 1] - DIF2_AVG_bias[kk/(100/Ndiv)])*float(j2 - j1 + 2)/(40./float(Ndiv));}
          }   
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
          Int_t iww1 = 10000 + iw%100;
	  if (k == iww1 + 1) h1_COV1_OVL->Fill(Sum_DIF1[k]/Cal1);
          if (nt1 < 10000) {
            TPeak1[nt1] = float(iw)*NS + k;
            nt1++;
          }
          h1_Peak1->Fill(Sum_DIF1[k]);
          h1_Peak1_100mV->Fill(Sum_DIF1[k]);
          h1_Peak1_Am241_AVG[inw]->Fill(Sum_DIF1[k]/Cal1);
          h1_Peak1_MeV->Fill(Sum_DIF1[k]/Cal1/1000.);
          h1_Peak1_4MeV->Fill(Sum_DIF1[k]/Cal1/1000.);
	  /*
          T_COV1[ip1_ms] = 1000.*(iw%FRed) + float(k)/100.;
          P_COV1[ip1_ms] = Sum_DIF1[k]/Cal1;
          np1++;
          ip1_ms++;
	  */
        }      
        if (Flag_Peak2[k]) {
          Int_t iww1 = 10000 + iw%100;
	  if (k == iww1 + 1) h1_COV2_OVL->Fill(Sum_DIF2[k]/Cal2);
          if (nt2 < 10000) {
            TPeak2[nt2] = float(iw)*NS + k;
            nt2++;
          }
          h1_Peak2->Fill(Sum_DIF2[k]);
          h1_Peak2_100mV->Fill(Sum_DIF2[k]);
          h1_Peak2_MeV->Fill(Sum_DIF2[k]/Cal2/1000.);
          h1_Peak2_4MeV->Fill(Sum_DIF2[k]/Cal2/1000.);
	  /*
          T_COV2[ip2_ms] = 1000.*(iw%FRed) + float(k)/100.;
          P_COV2[ip2_ms] = Sum_DIF2[k]/Cal2;
          np2++;
          ip2_ms++;
	  */
        }      
        if (Flag_Peak1[k] && Flag_Peak2[k]) {
          h1_Peak1_C_keV->Fill(Sum_DIF1[k]/Cal1);
          h1_Peak2_C_keV->Fill(Sum_DIF2[k]/Cal2);
          if (Sum_DIF1[k]/Cal1 <= 250.) h1_Peak1_C_250keV->Fill(Sum_DIF1[k]/Cal1);
          if (Sum_DIF2[k]/Cal2 <= 250.) h1_Peak2_C_250keV->Fill(Sum_DIF2[k]/Cal2);
          h1_Peak1_C_MeV->Fill(Sum_DIF1[k]/Cal1/1000.);
          h1_Peak2_C_MeV->Fill(Sum_DIF2[k]/Cal2/1000.);
	  /*
          T_COV1and2[ip1and2_ms] = 1000.*(iw%FRed) + float(k)/100.;;
          ip1and2_ms++;
	  */
        }
      }

      for (int k = 0; k < NS/Ndiv; k++) {

        if (Flag_Peak1_AVG[k]) {
          T_COV1[ip1_ms] = 1000.*(iw%FRed) + float(Ndiv)*float(k)/100.;
          P_COV1[ip1_ms] = Sum_DIF1_AVG[k]/Cal1;
          np1++;
          ip1_ms++;
        }
        if (Flag_Peak2_AVG[k]) {
          T_COV2[ip2_ms] = 1000.*(iw%FRed) + float(Ndiv)*float(k)/100.;
          P_COV2[ip2_ms] = Sum_DIF2_AVG[k]/Cal2;
          np2++;
          ip2_ms++;
        }
        if (Flag_Peak1_AVG[k] && Flag_Peak2_AVG[k]) {
          h1_Peak1_C_MeV_AVG->Fill(Sum_DIF1_AVG[k]/Cal1/1000.);
          h1_Peak2_C_MeV_AVG->Fill(Sum_DIF2_AVG[k]/Cal2/1000.);
          T_COV1and2[ip1and2_ms] = 1000.*(iw%FRed) + float(Ndiv)*float(k)/100.;
          ip1and2_ms++;
        }
        if (Flag_Peak1_AVG[k]) {h1_Peak1_100mV_AVG->Fill(Sum_DIF1_AVG[k]);}      
        if (Flag_Peak2_AVG[k]) {h1_Peak2_100mV_AVG->Fill(Sum_DIF2_AVG[k]);}      
        if (Flag_Peak1_AVG[k]) {h1_Peak1_MeV_AVG->Fill(Sum_DIF1_AVG[k]/Cal1/1000.);}      
        if (Flag_Peak2_AVG[k]) {h1_Peak2_MeV_AVG->Fill(Sum_DIF2_AVG[k]/Cal2/1000.);}      
        if (Flag_Peak1_AVG[k]) {h1_Peak1_4MeV_AVG->Fill(Sum_DIF1_AVG[k]/Cal1/1000.);}      
        if (Flag_Peak2_AVG[k]) {h1_Peak2_4MeV_AVG->Fill(Sum_DIF2_AVG[k]/Cal2/1000.);}      
        if (Flag_Peak1_1kHz_AVG[k]) {h1_Peak1_100mV_1kHz_AVG->Fill(Sum_DIF1_AVG[k]);}      
        if (Flag_Peak2_1kHz_AVG[k]) {h1_Peak2_100mV_1kHz_AVG->Fill(Sum_DIF2_AVG[k]);}      
        if (Flag_Peak1_AVG[k]) {h1_Peak1_200keV_AVG->Fill(Sum_DIF1_AVG[k]/Cal1);}      
        if (Flag_Peak2_AVG[k]) {h1_Peak2_200keV_AVG->Fill(Sum_DIF2_AVG[k]/Cal2);}      
      }

      //      iw++;

      if (iw == iw_sel) {
        for (int ii = 0; ii < NS; ii++) {
          h1_Raw1->SetBinContent(ii + 1, Raw1[ii]);
          h1_Raw2->SetBinContent(ii + 1, Raw2[ii]);
          h1_Raw3->SetBinContent(ii + 1, Raw3[ii]);
          h1_Raw4->SetBinContent(ii + 1, Raw4[ii]);
          h1_COV1->SetBinContent(ii + 1, COV1[ii]);
          h1_COV2->SetBinContent(ii + 1, COV2[ii]);
          if (ii > 0) {  
            h1_DIF1->SetBinContent(ii + 1, COV1[ii] - COV1[ii - 1]);
            h1_DIF2->SetBinContent(ii + 1, COV2[ii] - COV2[ii - 1]);
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
      }

      iw++;
      
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

      nBSL_NTD = 0;
      BSL_NTD_M = 0.;
      BSL_NTD_S = 0.;

      nBSL_RAW = 0;
      BSL_RAW_M = 0.;
      BSL_RAW_S = 0.;

      for (int k = 1; k < NS_NTD/FRed; k++) {
        DIF_NTD = NTD[k - 1] - NTD[k];
        if (abs(DIF_NTD) < NTD_Threshold) {
          BSL_NTD_M = BSL_NTD_M + DIF_NTD;
          BSL_NTD_S = BSL_NTD_S + DIF_NTD*DIF_NTD;
          BSL_RAW_M = BSL_RAW_M + NTD[k];
          BSL_RAW_S = BSL_RAW_S + NTD[k]*NTD[k];
          nBSL_NTD++;
          nBSL_RAW++;
        }
        if ((k + 1)%(NS_NTD/FRed) == 0) {
          BSL_RAW_M = BSL_RAW_M/max(1,nBSL_RAW) - 0.5825; 
          BSL_RAW_S = sqrt(BSL_RAW_S/max(1,nBSL_RAW) - BSL_RAW_M*BSL_RAW_M);
          if (nBSL_RAW > 1) {
            //h1_BSLM_AMP[ik]->Fill(BSL_RAW_M - BSLM_AMP_Shift[ik]);
            //h1_BSLS_AMP[ik]->Fill(BSL_RAW_S);
            //h2_Time_BSLM_NTD->Fill(nw1/NS/360.,BSL_RAW_M - BSLM_AMP_Shift[ik]);
            Int_t ie = (nw1/NS + nw2/NS)/720.;
	    h1_BSLM_RAW->Fill(BSL_RAW_M);
	    h2_Time_BSLM_RAW->Fill(float(ie),BSL_RAW_M);
            BSL_NTD_RAW[ie] = BSL_RAW_M;
          }
	}
        if ((k + 1)%100 == 0) {
          BSL_NTD_M = BSL_NTD_M/max(1,nBSL_NTD); 
          BSL_NTD_S = sqrt(BSL_NTD_S/max(1,nBSL_NTD) - BSL_NTD_M*BSL_NTD_M);
          if (nBSL_NTD > 1) {
            h1_BSL_NTD_M->Fill(BSL_NTD_M/CalNTD);
            h1_BSL_NTD_S->Fill(BSL_NTD_S/CalNTD);
          }
//
          nBSL_NTD = 0;
          BSL_NTD_M = 0.;
          BSL_NTD_S = 0.;
        }
      }

      for (int k = 1; k < NS_NTD/FRed; k++) {
        if (k > k_max) { 
          DIF_NTD = NTD[k - 1] - NTD[k];
          k_min = k - 1;
          k_max = k - 1;
          if (DIF_NTD > 0.) {
            k_max++;
            for (int l = k + 1; l < NS_NTD/FRed; l++) {
              if (NTD[l - 1] - NTD[l] > 0.) {k_max++;}
              else {break;}              
            }
            Int_t ie = (nw1/NS + nw2/NS)/720;
            if (k_max - k_min > 5 && k_max < 99990 && NTD[k_min] - NTD[k_max] > NTD_Threshold) {
              if (k_min > 99970) cout << " WARNING !!! " << k_min << " " << k_max << " " << iw_NTD << " " << NTD[k_min] - NTD[k_max] << " " << endl; 
              h1_T_NTD->Fill(float(k_max));
              h1_NTD_Amplitude->Fill((NTD[k_min] - NTD[k_max])*(1. - .01893*BSL_NTD_RAW[ie]));
            }
            if (NTD[k_min] - NTD[k_max] > NTD_Threshold) h2_Width_Amplitude->Fill(float(k_max - k_min),NTD[k_min] - NTD[k_max]);
            if (k_max - k_min > 5 && k_max < 99990 && NTD[k_min] - NTD[k_max] > NTD_Threshold) {
              //T_NTD[ip_NTD] = iw_NTD*NS + k_min;
              Int_t ie = (nw1/NS + nw2/NS)/720;
              P_NTD[ip_NTD] = (NTD[k_min] - NTD[k_max])*(1. - .01893*BSL_NTD_RAW[ie]);

              Int_t ll0 = k_min;
	      Double_t sl_max = -99999.;
              for (int ll = k_min + 1; ll < k_max + 1; ll++) {
                if (NTD[ll - 1] - NTD[ll] > sl_max) {
                  ll0 = ll - 1;
                  sl_max = NTD[ll - 1] - NTD[ll];
                }                  
              }
              Double_t k0 =  ll0 - (NTD[k_min] - NTD[ll0])/sl_max;
              T_NTD[ip_NTD] = k0*float(FRed)/float(100) - 1.0;
              //Int_t ie = (nw1/NS + nw2/NS)/720;
              if (P_NTD[ip_NTD]/CalNTD > h1_Alpha_Min && P_NTD[ip_NTD]/CalNTD < h1_Alpha_Max) {
                h1_Alpha[ie]->Fill(P_NTD[ip_NTD]/CalNTD);
                h2_BSLM_E_Alpha->Fill(BSL_NTD_RAW[ie],P_NTD[ip_NTD]/CalNTD);
              }
              ip_NTD++;
              np_NTD++;
            }
          }
        }   
      }

     if (iw_NTD == iw_NTD_sel && Flag_NTD == false) {
        for (int ii = 0; ii < NS_NTD/FRed; ii++) {
          h1_NTD->SetBinContent(ii + 1, NTD[ii]);
	  if (ii > 0) h1_DIF_NTD->SetBinContent(ii + 1, -(NTD[ii] - NTD[ii - 1]));
        }          
      }

      //cout << " ip_NTD = " << ip_NTD << " " << endl;

      for (int ip = 0; ip < ip_NTD; ip++) {
        Bool_t Flag_COV1 = false;
        Bool_t Flag_COV2 = false;
        Bool_t Flag_noCOV1 = false;
        Bool_t Flag_noCOV2 = false;
        Bool_t Flag_noCOV1andCOV2_11 = false;
        Bool_t Flag_noCOV1andCOV2_12 = false;
        Bool_t Flag_noCOV1andCOV2_21 = false;
        Bool_t Flag_noCOV1andCOV2_22 = false;
        Bool_t Flag_COV1andCOV2 = false;
        Bool_t Flag_noCOV1_M = false;
        Bool_t Flag_noCOV2_M = false;
        Bool_t Flag_noCOV1_P = false;
        Bool_t Flag_noCOV2_P = false;
         
        Double_t T_Closest;
        Double_t T_Closest_1 = 999999.;
        Double_t T_Closest_2 = 999999.;
        Double_t DT_Min = 999999.;
        Double_t DT1, DT2;
        h1_NTD_1_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        for (int jp = 0; jp < ip_COV1; jp++) {
          DT1 = T_NTD[ip] - T_COV1[jp];
          if (abs(DT1) < DT_Min) {
	    DT_Min = abs(DT1);
            T_Closest_1 = DT1;
	  }
	  if (abs(T_NTD[ip] - T_COV1[jp]) <= 50.) {
            if (abs(T_NTD[ip] - T_COV1[jp]) <= 25.) h1_DeltaT_NTD_COV1->Fill(T_NTD[ip] - T_COV1[jp]);
            h2_DeltaT_NTD_COV1_P_NTD->Fill(log10(P_NTD[ip]/CalNTD/1000.),T_NTD[ip] - T_COV1[jp]);
          } 
          if ((T_NTD[ip] - T_COV1[jp]) >= -2. && (T_NTD[ip] - T_COV1[jp]) < 4.) {
            Flag_COV1 = true;
            h1_P_COV1->Fill(P_COV1[jp]/1000.);
            h1_NTD_21_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          }
          if ((T_NTD[ip] - T_COV1[jp] - DTacc) >= -2. && (T_NTD[ip] - T_COV1[jp] - DTacc) < 4.) {
            Flag_noCOV1 = true;
            Flag_noCOV1_M = true;
            h1_P_noCOV1->Fill(P_COV1[jp]/1000.);
            h1_P_noCOV1_DTacc->Fill(P_COV1[jp]/1000.);
            h1_NTD_22_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
           }
	  
          if ((T_NTD[ip] - T_COV1[jp] + DTacc) >= -2. && (T_NTD[ip] - T_COV1[jp] + DTacc) < 4.) {
            Flag_noCOV1 = true;
            Flag_noCOV1_P = true;
            h1_P_noCOV1->Fill(P_COV1[jp]/1000.);
            h1_P_noCOV1_DTacc->Fill(P_COV1[jp]/1000.);
            h1_NTD_23_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          }
	 
        }
        for (int jp = 0; jp < ip_COV2; jp++) {
          DT2 = T_NTD[ip] - T_COV2[jp];
          if (abs(DT2) < DT_Min) {
	    DT_Min = abs(DT2);
            T_Closest_2 = DT2;
	  }
	  if (abs(T_NTD[ip] - T_COV2[jp]) <= 50.) {
            if (abs(T_NTD[ip] - T_COV2[jp]) <= 25.) h1_DeltaT_NTD_COV2->Fill(T_NTD[ip] - T_COV2[jp]); 
            h2_DeltaT_NTD_COV2_P_NTD->Fill(log10(P_NTD[ip]/CalNTD/1000.),T_NTD[ip] - T_COV2[jp]);
          } 
          if ((T_NTD[ip] - T_COV2[jp]) >= -2. && (T_NTD[ip] - T_COV2[jp]) < 4.) {
            Flag_COV2 = true;
            h1_P_COV2->Fill(P_COV2[jp]/1000.);
            h1_NTD_24_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          }
	  if ((T_NTD[ip] - T_COV2[jp] - DTacc) >= -2. && (T_NTD[ip] - T_COV2[jp] - DTacc) < 4.) {
            Flag_noCOV2 = true;
            Flag_noCOV2_M = true;
            h1_P_noCOV2->Fill(P_COV2[jp]/1000.);
            h1_P_noCOV2_DTacc->Fill(P_COV2[jp]/1000.);
            h1_NTD_25_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          }
	  
	  if ((T_NTD[ip] - T_COV2[jp] + DTacc) >= -2. && (T_NTD[ip] - T_COV2[jp] + DTacc) < 4.) {
            Flag_noCOV2 = true;
            Flag_noCOV2_P = true;
            h1_P_noCOV2->Fill(P_COV2[jp]/1000.);
            h1_P_noCOV2_DTacc->Fill(P_COV2[jp]/1000.);
            h1_NTD_26_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
          }	 
        }
        T_Closest = T_Closest_1;
        if (abs(T_Closest_2) < abs(T_Closest_1)) T_Closest = T_Closest_2;

        h1_DeltaT_Closest->Fill(T_Closest);
        h2_E_NTD_DeltaT_Closest->Fill(P_NTD[ip]/CalNTD,T_Closest);
        Int_t ie = int(2.*log(P_NTD[ip]/CalNTD));
        if (ie < 20) {
          h1_DeltaT_Closest_E[ie]->Fill(T_Closest);
          h1_IE[ie]->Fill(log(P_NTD[ip]/CalNTD));
        }
        if (T_Closest >=  -2. && T_Closest <   4.) h1_NTD_11_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (T_Closest <   -2. || T_Closest >=  4.) h1_NTD_12_MeV->Fill(P_NTD[ip]/CalNTD/1000.);      
        if (T_Closest >=  -9. && T_Closest <  -5.) h1_NTD_13_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (T_Closest >=   5. && T_Closest <   9.) h1_NTD_13_MeV->Fill(P_NTD[ip]/CalNTD/1000.);

        if (Flag_COV1 == false && Flag_COV2 == false) h1_NTD_2_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if ((Flag_noCOV1_M && !Flag_COV2) || (Flag_noCOV2_M && !Flag_COV1)) h1_NTD_2_MeV_DTacc->Fill(P_NTD[ip]/CalNTD/1000.);
        if ((Flag_noCOV1_P && !Flag_COV2) || (Flag_noCOV2_P && !Flag_COV1)) h1_NTD_2_MeV_DTacc->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV1 == true) h1_NTD_3_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV2 == true) h1_NTD_4_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV1 == true && Flag_COV2 == false) h1_NTD_7_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV2 == true && Flag_COV1 == false) h1_NTD_8_MeV->Fill(P_NTD[ip]/CalNTD/1000.);

        for (int jp1 = 0; jp1 < ip_COV1; jp1++) {
          DT1 = T_NTD[ip] - T_COV1[jp1];
          if (abs(DT1) <= 50.) {
            for (int jp2 = 0; jp2 < ip_COV2; jp2++) {
	      DT2 = T_NTD[ip] - T_COV2[jp2]; 
	      if (abs(DT2) <= 50.) h2_NTD_T1_T2->Fill(DT1,DT2);
            }  
          }
        }

        for (int jp = 0; jp < ip_COV1and2; jp++) {
          if (abs(T_NTD[ip] - T_COV1and2[jp]) <= 25.) h1_DeltaT_NTD_COV1and2->Fill(T_NTD[ip] - T_COV1and2[jp]);
        }
//
        for (int jp1 = 0; jp1 < ip_COV1; jp1++) {
	  if (T_NTD[ip] - T_COV1[jp1] - DTacc/2. >= -2. && T_NTD[ip] - T_COV1[jp1] - DTacc/2. < 4.) { 
            for (int jp2 = 0; jp2 < ip_COV2; jp2++) {
              if (T_NTD[ip] - T_COV2[jp2] - DTacc/2. >= -2. && T_NTD[ip] - T_COV2[jp2] - DTacc/2. < 4.) {
                Flag_noCOV1andCOV2_11 = true;
                h1_NTD_27_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
              }
              if (T_NTD[ip] - T_COV2[jp2] + DTacc/2. >= -2. && T_NTD[ip] - T_COV2[jp2] + DTacc/2. < 4.) {
                Flag_noCOV1andCOV2_12 = true;
                h1_NTD_28_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
              }
            }
          }
	  if (T_NTD[ip] - T_COV1[jp1] + DTacc/2. >= -2. && T_NTD[ip] - T_COV1[jp1] + DTacc/2. < 4.) { 
            for (int jp2 = 0; jp2 < ip_COV2; jp2++) {
              if (T_NTD[ip] - T_COV2[jp2] - DTacc/2. >= -2. && T_NTD[ip] - T_COV2[jp2] - DTacc/2. < 4.) { 
                Flag_noCOV1andCOV2_21 = true;
                h1_NTD_28_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
	      }
	      if (T_NTD[ip] - T_COV2[jp2] + DTacc/2. >= -2. && T_NTD[ip] - T_COV2[jp2] + DTacc/2. < 4.) {
                Flag_noCOV1andCOV2_22 = true;
                h1_NTD_27_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
              }
            }
          }	
        }
//
        if (Flag_noCOV1andCOV2_11 || Flag_noCOV1andCOV2_12 || Flag_noCOV1andCOV2_21 || Flag_noCOV1andCOV2_22 ) h1_NTD_5_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV1andCOV2_11 == true) h1_NTD_9_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV1andCOV2_22 == true) h1_NTD_9_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV1andCOV2_12 == true) h1_NTD_10_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
        if (Flag_noCOV1andCOV2_21 == true) h1_NTD_10_MeV->Fill(P_NTD[ip]/CalNTD/1000.);
//
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
    if (i < 5000) {
      cc = Raw1_Overlay[i]/float(nov1);
      if (i < 5) cout << " 1 " << i << " " << cc << " " << endl;
      h1_Raw1_Overlay->SetBinContent(i + 1, cc);
      cc = Raw2_Overlay[i]/float(nov1);
      if (i < 5) cout << " 2 " << i <<  " " <<cc << " " << endl;
      h1_Raw2_Overlay->SetBinContent(i + 1, cc);
      cc = Raw3_Overlay[i]/float(nov2);
      if (i < 5) cout << " 3 " << i << " " << cc << " " << endl;
      h1_Raw3_Overlay->SetBinContent(i + 1, cc);
      cc = Raw4_Overlay[i]/float(nov2);
      if (i < 5) cout << " 4 " << i << " " << cc << " " << endl;
      h1_Raw4_Overlay->SetBinContent(i + 1, cc);
    }
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

  cout << " End processing events " << endl;

  fHistos.Write(); 
  cout << " End writing histograms " << endl;

  fHistos.Close();
 
  cout << " " << endl;
  cout << " Output root file closed " << endl;	
  cout << " " << endl;

}
