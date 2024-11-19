/* Veto.C*/

{

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Open histogram files                                                    //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

  //TFile *f_01 = new TFile("/sps/t2k/edoardo/COV/MyCOV1_Test.root");
  TFile *f_01 = new TFile("../../Data_COV/2021_05_12_17_27_27/2021_05_12_17_27_27_LWOthr-50keV.root");

  f_01->cd();

  Double_t Tbin, LL, UL, XL, P0, P1, P0L, P1L, P0R, P1R;
  Int_t IL;

  Double_t RunTime = h1_RunTime->GetBinContent(1);
  cout << " RunTime = " << RunTime << "s" << endl;
  Double_t N_COV1 = h1_Peak1_MeV_AVG->GetEntries();
  cout << " COV1 candidates = " << N_COV1 << "      Rate = " << float(N_COV1)/RunTime << "/s " << endl;
  Double_t N_COV2 = h1_Peak2_MeV_AVG->GetEntries();
  cout << " COV2 candidates = " << N_COV2 << "      Rate = " << float(N_COV2)/RunTime << "/s " << endl;
  Double_t N_COV1and2 = h1_Peak1_C_MeV_AVG->GetEntries();
  cout << " COV1and2 candidates = " << N_COV1and2 << "      Rate = " << float(N_COV1and2)/RunTime << "/s " << endl;
  cout << "      Rate = " << fixed << setprecision(3) << setw(5) << setfill(' ') << N_COV1and2 << " " << endl;

   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Canvas 001 :                                                            //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 001 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  TCanvas *C_001 = new TCanvas("C_001","C_001: ", 1000,600);
  C_001->Divide(1,1);
  C_001->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000111");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  h1_DeltaT_Closest->GetXaxis()->SetTitle("#DeltaT [ms]");
  h1_DeltaT_Closest->GetYaxis()->SetTitle("Counts / .2 ms");
  h1_DeltaT_Closest->GetXaxis()->SetRangeUser(-25.0, 25.0);
  h1_DeltaT_Closest->GetYaxis()->SetRangeUser(  10.,100000.0);
  h1_DeltaT_Closest->SetLineColor(4);
  h1_DeltaT_Closest->SetLineWidth(2);
    h1_DeltaT_Closest->Rebin(2);
  h1_DeltaT_Closest->Draw();

  TF1 *fit01_L = new TF1("fit01_L","expo",-25.0, -2.0);
  TF1 *fit01_R = new TF1("fit01_R","expo",  4.0, 25.0);

  h1_DeltaT_Closest->Fit("fit01_L","NR");
  h1_DeltaT_Closest->Fit("fit01_R","NR");

  P0 = 0.5*(fit01_L->GetParameter(0) + fit01_R->GetParameter(0));
  P1 = 0.5*(fit01_L->GetParameter(1) - fit01_R->GetParameter(1));

  cout << " P0, P1 = " << P0 << " " << P1 << " " << endl; 

  TF1 *fit01_L0 = new TF1("fit01_L0","expo",-25.0,  0.0);
  TF1 *fit01_R0 = new TF1("fit01_R0","expo", -0.0, 25.0);

  fit01_L0->SetLineStyle(0);
  fit01_L0->SetLineColor(2);
  fit01_L0->SetLineWidth(2);
  fit01_L0->SetParameters(P0,P1);
  fit01_L0->Draw("SAME");

  fit01_R0->SetLineStyle(0);
  fit01_R0->SetLineColor(2);
  fit01_R0->SetLineWidth(2);
  fit01_R0->SetParameters(P0,-P1);
  fit01_R0->Draw("SAME");

  Double_t N_NTD = h1_DeltaT_Closest->GetEntries();
  cout << " NTD candidates = " << N_NTD << " " << endl;

  LL = -2.;
  UL =  3.;

  Double_t Sum = 0.;
    Sum = h1_DeltaT_Closest->Integral(h1_DeltaT_Closest->FindBin(LL), h1_DeltaT_Closest->FindBin(UL));
    //cout << " " << j << " " << h1_DeltaT_Closest->GetBinContent(240 + j + 1) << " " << Sum << " " << endl;
  cout << " NTD candidates with VETO = " << Sum << " " << endl;

  Tbin = h1_DeltaT_Closest->GetBinWidth(1);
  Double_t Sum_L = exp(P0)*(1. - exp( P1*LL))/P1/Tbin;
  Double_t Sum_R = exp(P0)*(1. - exp(-P1*UL))/P1/Tbin;

  cout << " Accidental VETO = " << Sum_L << " + " << Sum_R << " = " << Sum_L + Sum_R << " " << endl; 
  Double_t ACC = Sum_L + Sum_R;  
  cout << " NTD candidates with true VETO = " << Sum -  Sum_L - Sum_R << " " << endl;
  cout << " True VETO rejection = " << (Sum -  Sum_L - Sum_R)/N_NTD << " " << endl;
  
  TLine *Line01_1 = new TLine( -2., 10., -2.,100000.);  
  Line01_1->SetLineStyle(2);
  Line01_1->SetLineColor(2);
  Line01_1->SetLineWidth(4);
  Line01_1->Draw();

  TLine *Line01_2 = new TLine(  4., 10.,  4.,100000.);  
  Line01_2->SetLineStyle(2);
  Line01_2->SetLineColor(2);
  Line01_2->SetLineWidth(4);
  Line01_2->Draw();
  /*
  TLine *Line01_3 = new TLine( -9., 10., -9.,500.);  
  Line01_3->SetLineStyle(2);
  Line01_3->SetLineColor(1);
  Line01_3->SetLineWidth(4);
  Line01_3->Draw();

  TLine *Line01_4 = new TLine( -5., 10., -5.,500.);  
  Line01_4->SetLineStyle(2);
  Line01_4->SetLineColor(1);
  Line01_4->SetLineWidth(4);
  Line01_4->Draw();

  TLine *Line01_5 = new TLine(  5., 10.,  5.,500.);  
  Line01_5->SetLineStyle(2);
  Line01_5->SetLineColor(1);
  Line01_5->SetLineWidth(4);
  Line01_5->Draw();

  TLine *Line01_6 = new TLine(  9., 10.,  9.,500.);  
  Line01_6->SetLineStyle(2);
  Line01_6->SetLineColor(1);
  Line01_6->SetLineWidth(4);
  Line01_6->Draw();
  */
  TLegend *Legend_01 = new TLegend(.15,.70,.45,.85);
  Legend_01->SetFillColor(0);
  Legend_01->SetBorderSize(0);
  Legend_01->AddEntry(h1_DeltaT_NTD_COV1,"T_{LWO}-T_{Closest Veto}","L");
  Legend_01->Draw();
     
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Canvas 002 :                                                            //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 002 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  TCanvas *C_002 = new TCanvas("C_002","C_002: ", 1000,600);
  C_002->Divide(1,1);
  C_002->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000111");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  h1_DeltaT_NTD_COV1and2->GetXaxis()->SetTitle("#DeltaT [ms]");
  h1_DeltaT_NTD_COV1and2->GetYaxis()->SetTitle("Counts / .2 ms");
  h1_DeltaT_NTD_COV1and2->GetXaxis()->SetRangeUser(-25.0, 25.0);
  h1_DeltaT_NTD_COV1and2->GetYaxis()->SetRangeUser(  1.,10000.0);
  h1_DeltaT_NTD_COV1and2->SetLineColor(4);
  h1_DeltaT_NTD_COV1and2->SetLineWidth(2);
  h1_DeltaT_NTD_COV1and2->Draw();

  TF1 *fit02_L = new TF1("fit02_L","expo",-25.0, -2.0);
  TF1 *fit02_R = new TF1("fit02_R","expo",  4.0, 25.0);

  h1_DeltaT_NTD_COV1and2->Fit("fit02_L","NR");
  h1_DeltaT_NTD_COV1and2->Fit("fit02_R","NR");

  P0 = 0.5*(fit02_L->GetParameter(0) + fit02_R->GetParameter(0));
  P1 = 0.5*(fit02_L->GetParameter(1) - fit02_R->GetParameter(1));

  cout << " P0, P1 = " << P0 << " " << P1 << " " << endl; 

  TF1 *fit02_L0 = new TF1("fit02_L0","expo",-50.0,  0.0);
  TF1 *fit02_R0 = new TF1("fit02_R0","expo", -0.0, 50.0);

  fit02_L0->SetLineStyle(0);
  fit02_L0->SetLineColor(2);
  fit02_L0->SetLineWidth(2);
  fit02_L0->SetParameters(P0,P1);
  fit02_L0->Draw("SAME");

  fit02_R0->SetLineStyle(0);
  fit02_R0->SetLineColor(2);
  fit02_R0->SetLineWidth(2);
  fit02_R0->SetParameters(P0,-P1);
  fit02_R0->Draw("SAME");

  Double_t N_NTD_COV1and2 = h1_DeltaT_NTD_COV1and2->GetEntries();
  cout << " NTD candidates = " << N_NTD_COV1and2 << " " << endl;

  
  XL = -25.;
  Tbin = 0.2;
  LL = -2.;
  UL =  4.;
  IL = int((LL - XL)/Tbin + .5);
  
  Sum = 0.;
  for (Int_t j = 0; j < int((UL - LL)/Tbin + .5); j++) {
    Sum = Sum + h1_DeltaT_NTD_COV1and2->GetBinContent(IL + j + 1);
    //cout << " " << j << " " << h1_DeltaT_NTD_COV1and2->GetBinContent(IL + j + 1) << " " << Sum << " " << endl;
  }
  cout << " NTD candidates with COV1and2 = " << Sum << " " << endl;

  Sum_L = exp(P0)*(1. - exp(+P1*LL))/P1/Tbin;
  Sum_R = exp(P0)*(1. - exp(-P1*UL))/P1/Tbin;

  cout << " Accidental VETO = " << Sum_L << " + " << Sum_R << " = " << Sum_L + Sum_R << " " << endl;    
  cout << " NTD candidates with true COV1and2 = " << Sum -  Sum_L - Sum_R << " " << endl;
  cout << " Rate of NTD candidates with true COV1and2 = " << (Sum -  Sum_L - Sum_R)/RunTime << "/s " << endl;
  cout << " True VETO rejection = " << (Sum -  Sum_L - Sum_R)/N_NTD << " " << endl;

  TLine *Line02_1 = new TLine( -2.,  1., -2.,10000.);  
  Line02_1->SetLineStyle(2);
  Line02_1->SetLineColor(2);
  Line02_1->SetLineWidth(4);
  Line02_1->Draw();

  TLine *Line02_2 = new TLine(  4.,  1.,  4.,10000.);  
  Line02_2->SetLineStyle(2);
  Line02_2->SetLineColor(2);
  Line02_2->SetLineWidth(4);
  Line02_2->Draw();

  TLegend *Legend_02 = new TLegend(.15,.70,.45,.85);
  Legend_02->SetFillColor(0);
  Legend_02->SetBorderSize(5);
  Legend_02->AddEntry(h1_DeltaT_NTD_COV1,"T_{NTD}-T_{BOT&TOP}","L");
  Legend_02->Draw();
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Canvas 003 :                                                            //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 003 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  TCanvas *C_003 = new TCanvas("C_003","C_003: ", 1000,600);
  C_003->Divide(1,1);
  C_003->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000111");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  h1_DeltaT_NTD_COV1->GetXaxis()->SetTitle("#DeltaT [ms]");
  h1_DeltaT_NTD_COV1->GetYaxis()->SetTitle("Counts / .2 ms");
  h1_DeltaT_NTD_COV1->GetXaxis()->SetRangeUser(-25.0, 25.0);
  h1_DeltaT_NTD_COV1->GetYaxis()->SetRangeUser( 100.,10000.0);
  h1_DeltaT_NTD_COV1->SetLineColor(4);
  h1_DeltaT_NTD_COV1->SetLineWidth(2);
  h1_DeltaT_NTD_COV1->Draw();

  TF1 *fit03_L = new TF1("fit03_L","expo",-25.0, -2.0);
  TF1 *fit03_R = new TF1("fit03_R","expo",  4.0, 25.0);

  h1_DeltaT_NTD_COV1->Fit("fit03_L","NR");
  h1_DeltaT_NTD_COV1->Fit("fit03_R","NR");

  P0 = 0.5*(fit03_L->GetParameter(0) + fit03_R->GetParameter(0));
  P1 = 0.5*(fit03_L->GetParameter(1) - fit03_R->GetParameter(1));

  cout << " P0, P1 = " << P0 << " " << P1 << " " << endl; 

  TF1 *fit03_L0 = new TF1("fit03_L0","expo",-50.0,  0.0);
  TF1 *fit03_R0 = new TF1("fit03_R0","expo", -0.0, 50.0);

  fit03_L0->SetLineStyle(0);
  fit03_L0->SetLineColor(2);
  fit03_L0->SetLineWidth(2);
  fit03_L0->SetParameters(P0,P1);
  fit03_L0->Draw("SAME");

  fit03_R0->SetLineStyle(0);
  fit03_R0->SetLineColor(2);
  fit03_R0->SetLineWidth(2);
  fit03_R0->SetParameters(P0,-P1);
  fit03_R0->Draw("SAME");

  Double_t N_NTD_COV1 = h1_DeltaT_NTD_COV1->GetEntries();
  cout << " NTD candidates = " << N_NTD_COV1 << " " << endl;

  LL = -2.;
  UL =  4.;

  Sum = 0.;
  for (Int_t j = 0; j < 30; j++) {
    Sum = Sum + h1_DeltaT_NTD_COV1->GetBinContent(115 + j + 1);
    //cout << " " << j << " " << h1_DeltaT_NTD_COV1->GetBinContent(115 + j + 1) << " " << Sum << " " << endl;
  }
  cout << " NTD candidates with VETO = " << Sum << " " << endl;

  Tbin = 0.2;
  Sum_L = exp(P0)*(1. - exp(+P1*LL))/P1/Tbin;
  Sum_R = exp(P0)*(1. - exp(-P1*UL))/P1/Tbin;

  cout << " Accidental VETO = " << Sum_L << " + " << Sum_R << " = " << Sum_L + Sum_R << " " << endl;    
  cout << " NTD candidates with true VETO = " << Sum -  Sum_L - Sum_R << " " << endl;
  cout << " Rate of true coincidences NTD*COV1 = " << (Sum -  Sum_L - Sum_R)/RunTime << "/s " << endl;
  cout << " True VETO rejection = " << (Sum -  Sum_L - Sum_R)/N_NTD << " " << endl;
   
  TLine *Line03_1 = new TLine( -2.,  100., -2.,10000.);  
  Line03_1->SetLineStyle(2);
  Line03_1->SetLineColor(2);
  Line03_1->SetLineWidth(4);
  Line03_1->Draw();

  TLine *Line03_2 = new TLine(  4.,  100.,  4.,10000.);  
  Line03_2->SetLineStyle(2);
  Line03_2->SetLineColor(2);
  Line03_2->SetLineWidth(4);
  Line03_2->Draw();

  TLegend *Legend_03 = new TLegend(.15,.70,.40,.85);
  Legend_03->SetFillColor(0);
  Legend_03->SetBorderSize(5);
  Legend_03->AddEntry(h1_DeltaT_NTD_COV1,"T_{NTD}-T_{BOT}","L");
  Legend_03->Draw();
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Canvas 004 :                                                            //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 004 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  TCanvas *C_004 = new TCanvas("C_004","C_004: ", 1000,600);
  C_004->Divide(1,1);
  C_004->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000111");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  h1_DeltaT_NTD_COV2->GetXaxis()->SetTitle("#DeltaT [ms]");
  h1_DeltaT_NTD_COV2->GetYaxis()->SetTitle("Counts / .2 ms");
  h1_DeltaT_NTD_COV2->GetXaxis()->SetRangeUser(-25.0, 25.0);
  h1_DeltaT_NTD_COV2->GetYaxis()->SetRangeUser(  100.,10000.0);
  h1_DeltaT_NTD_COV2->SetLineColor(4);
  h1_DeltaT_NTD_COV2->SetLineWidth(2);
  h1_DeltaT_NTD_COV2->Draw();

  TF1 *fit04_L = new TF1("fit04_L","expo",-25.0, -2.0);
  TF1 *fit04_R = new TF1("fit04_R","expo",  4.0, 25.0);

  h1_DeltaT_NTD_COV2->Fit("fit04_L","NR");
  h1_DeltaT_NTD_COV2->Fit("fit04_R","NR");

  P0 = 0.5*(fit04_L->GetParameter(0) + fit04_R->GetParameter(0));
  P1 = 0.5*(fit04_L->GetParameter(1) - fit04_R->GetParameter(1));

  cout << " P0, P1 = " << P0 << " " << P1 << " " << endl; 

  TF1 *fit04_L0 = new TF1("fit04_L0","expo",-50.0,  0.0);
  TF1 *fit04_R0 = new TF1("fit04_R0","expo", -0.0, 50.0);

  fit04_L0->SetLineStyle(0);
  fit04_L0->SetLineColor(2);
  fit04_L0->SetLineWidth(2);
  fit04_L0->SetParameters(P0,P1);
  fit04_L0->Draw("SAME");

  fit04_R0->SetLineStyle(0);
  fit04_R0->SetLineColor(2);
  fit04_R0->SetLineWidth(2);
  fit04_R0->SetParameters(P0,-P1);
  fit04_R0->Draw("SAME");

  Double_t N_NTD_COV2 = h1_DeltaT_NTD_COV2->GetEntries();
  cout << " NTD candidates = " << N_NTD_COV2 << " " << endl;

  
  XL = -25.;
  Tbin = 0.2;
  LL = -2.;
  UL =  4.;
  IL = int((LL - XL)/Tbin + .5);
  
  Sum = 0.;
  for (Int_t j = 0; j < int((UL - LL)/Tbin + .5); j++) {
    Sum = Sum + h1_DeltaT_NTD_COV2->GetBinContent(IL + j + 1);
    //cout << " " << j << " " << h1_DeltaT_NTD_COV2->GetBinContent(IL + j + 1) << " " << Sum << " " << endl;
  }
  cout << " NTD candidates with COV2 = " << Sum << " " << endl;

  Sum_L = exp(P0)*(1. - exp(+P1*LL))/P1/Tbin;
  Sum_R = exp(P0)*(1. - exp(-P1*UL))/P1/Tbin;

  cout << " Accidental VETO = " << Sum_L << " + " << Sum_R << " = " << Sum_L + Sum_R << " " << endl;    
  cout << " NTD candidates with true COV2 = " << Sum -  Sum_L - Sum_R << " " << endl;
  cout << " Rate of NTD candidates with true COV2 = " << (Sum -  Sum_L - Sum_R)/RunTime << "/s " << endl;
  cout << " True VETO rejection = " << (Sum -  Sum_L - Sum_R)/N_NTD << " " << endl;

  TLine *Line04_1 = new TLine( -2.,  100., -2.,10000.);  
  Line04_1->SetLineStyle(2);
  Line04_1->SetLineColor(2);
  Line04_1->SetLineWidth(4);
  Line04_1->Draw();

  TLine *Line04_2 = new TLine(  4.,  100.,  4.,10000.);  
  Line04_2->SetLineStyle(2);
  Line04_2->SetLineColor(2);
  Line04_2->SetLineWidth(4);
  Line04_2->Draw();

  TLegend *Legend_04 = new TLegend(.15,.70,.40,.85);
  Legend_04->SetFillColor(0);
  Legend_04->SetBorderSize(5);
  Legend_04->AddEntry(h1_DeltaT_NTD_COV2,"T_{NTD}-T_{TOP}","L");
  Legend_04->Draw();
  /*
  TLatex *Latex_004 = new TLatex( 15.,0.30,"BOT");
  Latex_004->SetTextSize(.080);
  Latex_004->SetTextColor(2);
  Latex_004->SetTextAlign(22);
  Latex_004->Draw();
  */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Canvas 021 :                                                            //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 021 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  TCanvas *C_021 = new TCanvas("C_021","C_021: ", 1000,600);
  C_021->Divide(1,1);
  C_021->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  h1_NTD_1_MeV->GetXaxis()->SetTitle("E [MeV]");
  h1_NTD_1_MeV->GetYaxis()->SetTitle("Counts / 10 keV");
  h1_NTD_1_MeV->GetXaxis()->SetRangeUser(0., .6);
  h1_NTD_1_MeV->GetYaxis()->SetRangeUser(1., 100000.);
  h1_NTD_1_MeV->SetLineColor(4);
  h1_NTD_1_MeV->SetLineWidth(2);
  h1_NTD_1_MeV->Draw();

  Double_t  NTD_1_MeV = h1_NTD_1_MeV->GetEntries();
  cout << " ****** " << NTD_1_MeV << " " << endl;

  h1_NTD_12_MeV->SetLineColor(2);
  h1_NTD_12_MeV->SetFillColor(2);
  h1_NTD_12_MeV->SetLineWidth(2);
  h1_NTD_12_MeV->DrawCopy("SAME");

  Double_t NTD_2_MeV = h1_NTD_2_MeV->GetEntries();
  Double_t NTD_3_MeV = h1_NTD_3_MeV->GetEntries();
  Double_t NTD_4_MeV = h1_NTD_4_MeV->GetEntries();
  Double_t NTD_5_MeV = h1_NTD_5_MeV->GetEntries();
  Double_t NTD_6_MeV = h1_NTD_6_MeV->GetEntries();
  Double_t NTD_7_MeV = h1_NTD_7_MeV->GetEntries();
  Double_t NTD_8_MeV = h1_NTD_8_MeV->GetEntries();
  Double_t NTD_9_MeV = h1_NTD_9_MeV->GetEntries();
  Double_t NTD_10_MeV = h1_NTD_10_MeV->GetEntries();
  Double_t NTD_11_MeV = h1_NTD_11_MeV->GetEntries();
  Double_t NTD_12_MeV = h1_NTD_12_MeV->GetEntries();
  Double_t NTD_13_MeV = h1_NTD_13_MeV->GetEntries();

  cout << " " << NTD_1_MeV << " " << NTD_11_MeV << " " << NTD_12_MeV << " " << NTD_13_MeV << " " << endl; 
  Double_t CORR = ACC/NTD_12_MeV;
  cout << " CORR = " << NTD_12_MeV << " " << ACC << " " << CORR << " " << endl;

  TLegend *Legend_021 = new TLegend(.60,.65,.85,.85);
  Legend_021->SetFillColor(0);
  Legend_021->SetBorderSize(5);
  Legend_021->AddEntry(h1_NTD_1_MeV,"NTD","L");
  Legend_021->AddEntry(h1_NTD_12_MeV,"NTD with Ge Veto" ,"L");
  Legend_021->Draw();

  TLatex *Latex_021 = new TLatex( .18,20000.,"Before correcting for accidentals");
  Latex_021->SetTextSize(.045);
  Latex_021->SetTextColor(1);
  Latex_021->SetTextAlign(22);
  Latex_021->Draw();
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Canvas 022 :                                                            //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 022 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  TCanvas *C_022 = new TCanvas("C_022","C_022: ", 1000,600);
  C_022->Divide(1,1);
  C_022->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  h1_NTD_1_MeV->GetXaxis()->SetTitle("E [MeV]");
  h1_NTD_1_MeV->GetYaxis()->SetTitle("Counts / 10 keV");
  h1_NTD_1_MeV->GetXaxis()->SetRangeUser(0.,   .6);
  h1_NTD_1_MeV->GetYaxis()->SetRangeUser(1., 100000.);
  h1_NTD_1_MeV->SetLineColor(4);
  h1_NTD_1_MeV->SetLineWidth(2);
  h1_NTD_1_MeV->Draw();

  TH1F *h1_NTD_Veto_MeV;

  h1_NTD_Veto_MeV = (TH1F*)h1_NTD_1_MeV->Clone();

  h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_11_MeV, 0., 1.);
  h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_12_MeV, 1., -1*CORR);

  cout << " Veto " << h1_NTD_Veto_MeV->Integral(0,2005) << " " << endl;

  h1_NTD_Veto_MeV->SetLineColor(2);
  h1_NTD_Veto_MeV->SetFillColor(2);
  h1_NTD_Veto_MeV->SetLineWidth(2);
  h1_NTD_Veto_MeV->DrawCopy("SAME");

  TLegend *Legend_022 = new TLegend(.60,.65,.85,.85);
  Legend_022->SetFillColor(0);
  Legend_022->SetBorderSize(0);
  Legend_022->AddEntry(h1_NTD_1_MeV,"NTD","L");
  Legend_022->AddEntry(h1_NTD_Veto_MeV,"NTD with Ge Veto" ,"L");
  Legend_022->Draw();

//  TLatex *Latex_022 = new TLatex( .18 ,20000.,"After correcting for accidentals");
//  Latex_022->SetTextSize(.045);
//  Latex_022->SetTextColor(1);
//  Latex_022->SetTextAlign(22);
//  Latex_022->Draw();
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    //                                                                         //
    // Canvas 023 :                                                            //
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////

      cout << " ++++++++++++++++++++++++ " << endl;
      cout << " ++++++ Canvas 023 ++++++ " << endl;
      cout << " ++++++++++++++++++++++++ " << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
      TCanvas *C_023 = new TCanvas("C_023","C_023: ", 1000,600);
      C_023->SetLogy();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      gStyle->SetOptFit(000000);
      gStyle->SetOptStat("000000001");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    h1_Peak1_MeV->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak1_MeV->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak1_MeV->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak1_MeV->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak1_MeV->SetLineColor(2);
    h1_Peak1_MeV->SetLineWidth(2);
    h1_Peak1_MeV->Draw();
    
    h1_Peak2_MeV->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak2_MeV->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak2_MeV->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak2_MeV->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak2_MeV->SetLineColor(4);
    h1_Peak2_MeV->SetLineWidth(2);
    h1_Peak2_MeV->Draw("SAME");

    h1_Peak2_C_MeV->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak2_C_MeV->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak2_C_MeV->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak2_C_MeV->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak2_C_MeV->SetLineColor(6);
    h1_Peak2_C_MeV->SetLineWidth(2);
    h1_Peak2_C_MeV->Draw("SAME");
    
    h1_Peak1_C_MeV->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak1_C_MeV->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak1_C_MeV->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak1_C_MeV->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak1_C_MeV->SetLineColor(3);
    h1_Peak1_C_MeV->SetLineWidth(2);
    h1_Peak1_C_MeV->Draw("SAME");
     
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///    /////////////////////////////////////////////////////////////////////////////
    //                                                                         //
    // Canvas 024 :                                                            //
    //                                                                         //
    /////////////////////////////////////////////////////////////////////////////

      cout << " ++++++++++++++++++++++++ " << endl;
      cout << " ++++++ Canvas 024 ++++++ " << endl;
      cout << " ++++++++++++++++++++++++ " << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
      TCanvas *C_024 = new TCanvas("C_024","C_024: ", 1000,600);
      C_024->SetLogy();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      gStyle->SetOptFit(000000);
      gStyle->SetOptStat("000000001");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    h1_Peak1_4MeV_AVG->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak1_4MeV_AVG->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak1_4MeV_AVG->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak1_4MeV_AVG->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak1_4MeV_AVG->SetLineColor(2);
    h1_Peak1_4MeV_AVG->SetLineWidth(2);
    h1_Peak1_4MeV_AVG->Draw();
    
    h1_Peak2_4MeV_AVG->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak2_4MeV_AVG->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak2_4MeV_AVG->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak2_4MeV_AVG->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak2_4MeV_AVG->SetLineColor(4);
    h1_Peak2_4MeV_AVG->SetLineWidth(2);
    h1_Peak2_4MeV_AVG->Draw("SAME");

    h1_Peak2_C_4MeV_AVG->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak2_C_4MeV_AVG->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak2_C_4MeV_AVG->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak2_C_4MeV_AVG->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak2_C_4MeV_AVG->SetLineColor(6);
    h1_Peak2_C_4MeV_AVG->SetLineWidth(2);
    h1_Peak2_C_4MeV_AVG->Draw("SAME");
    
    h1_Peak1_C_4MeV_AVG->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
    h1_Peak1_C_4MeV_AVG->GetYaxis()->SetTitle("Counts / (1 keV)");
    h1_Peak1_C_4MeV_AVG->GetXaxis()->SetRangeUser(0.,   20.);
    h1_Peak1_C_4MeV_AVG->GetYaxis()->SetRangeUser(1., 100000.);
    h1_Peak1_C_4MeV_AVG->SetLineColor(3);
    h1_Peak1_C_4MeV_AVG->SetLineWidth(2);
    h1_Peak1_C_4MeV_AVG->Draw("SAME");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
