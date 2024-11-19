/* RecEff_COV.C*/

{

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F *Frame_001 = C_001->DrawFrame(0., 0.,  200., 120);

  Frame_001->GetXaxis()->SetTitle("Reconstructed energy [keVee]");
  Frame_001->GetYaxis()->SetTitle("Efficiency [%]");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Int_t NT = 4000;
  Double_t E[18]  = {   2.,   4.,   6.,   8.,  10.,  12.,  15.,  18.,  20.,  30.,  40.,  60.,  80., 100., 120., 140., 160., 180.};
  Int_t NCOV1[18] = {    1,   34,  561, 2251, 3450, 3848, 3954, 3976, 3981, 3981, 3985, 3985, 3987, 3987, 3987, 3988, 3988, 3989};
  Int_t NCOV2[18] = {    3,    7,   21,   44,  152,  658, 2489, 3697, 3889, 3990, 3991, 3992, 3994, 3995, 3995, 3995, 3995, 3995};

  Double_t Eff1_M[18], Eff1_S[18], Eff2_M[18], Eff2_S[18];

  Int_t nn1 = 18;
  Double_t xx1[18], yy1[18], ex1[18], ey1[18];

  for (Int_t j = 0; j < nn1; j++) {
    
    xx1[j] = E[j];
    ex1[j] =   0.;
    yy1[j] = float(NCOV1[j])/float(NT)*100;
    ey1[j] = sqrt(yy1[j]*(100. - yy1[j])/float(NT));
  }

  TGraphErrors *Eff1 = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
  Eff1->SetMarkerColor(2);
  Eff1->SetMarkerStyle(21);
  Eff1->SetMarkerSize(1.2);
  Eff1->Draw("P");

  for (Int_t j = 0; j < nn1; j++) {
    xx1[j] = E[j];
    ex1[j] =   0.;
    yy1[j] = float(NCOV2[j])/float(NT)*100;
    ey1[j] = sqrt(yy1[j]*(100. - yy1[j])/float(NT));
  }

  TGraphErrors *Eff2 = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
  Eff2->SetMarkerColor(4);
  Eff2->SetMarkerStyle(22);
  Eff2->SetMarkerSize(1.2);
  Eff2->Draw("P SAME");

  TLegend *Legend_001 = new TLegend(.50,.25,.85,.45);
  Legend_001->SetFillColor(0);
  Legend_001->SetBorderSize(0);
  Legend_001->SetTextSize(.04);
  Legend_001->AddEntry(Eff1,"BOT Ge - 5 keV threshold","P");
  Legend_001->AddEntry(Eff2,"TOP Ge - 10 keV threshold","P");
  Legend_001->Draw();

  TF1 *fit_eff1 = new TF1("fit_eff1","[0] +[1]*tanh([2]*x + [3])",  0.,200.);
  fit_eff1->SetParameters(27,27,0.5,-4.0);
  fit_eff1->SetLineColor(2);
  fit_eff1->SetLineWidth(2);
  fit_eff1->SetLineStyle(2);
  Eff1->Fit("fit_eff1","R");

  TF1 *fit_eff2 = new TF1("fit_eff2","[0] +[1]*tanh([2]*x + [3])",  0.,200.);
  fit_eff2->SetParameters(50,50,0.4,-5.0);
  fit_eff2->SetLineColor(4);
  fit_eff2->SetLineWidth(2);
  fit_eff2->SetLineStyle(2);
  Eff2->Fit("fit_eff2","R");

  /*
  TLatex *Latex_001_1 = new TLatex( 511.*Cal1/1000.,1500.,"511 keV");
  Latex_001_1->SetTextSize(.030);
  Latex_001_1->SetTextColor(1);
  Latex_001_1->SetTextAlign(22);
  Latex_001_1->Draw();
  */

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TH1F *Frame_002 = C_002->DrawFrame(0., 0., 1200., 120.);

  Frame_002->GetXaxis()->SetTitle("E [keV]");
  Frame_002->GetYaxis()->SetTitle("Efficiency");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  NT = 320;
  Double_t E_NTD[15]  = {   15.,   20.,   25.,   30.,   40.,    50.,   75.,  100.,  150.,  200.,  300.,  400.,  500.,  750., 1000.};
  Int_t N_NTD[15]     = {     5,     9,    19,    43,   138,    194,   231,   254,   276,   285,   294,   297,   300,   305,   308};

  Double_t Eff_NTD_M[13], Eff_NTD_S[13];

  nn1 = 15;

  for (Int_t j = 0; j < nn1; j++) {
    
    xx1[j] = E_NTD[j];
    ex1[j] =   0.;
    yy1[j] = float(N_NTD[j])/float(NT)*100.;
    ey1[j] = sqrt(yy1[j]*(100. - yy1[j])/float(NT));
  }

  TGraphErrors *Eff3 = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
  Eff3->SetMarkerColor(1);
  Eff3->SetMarkerStyle(20);
  Eff3->SetMarkerSize(1.2);
  Eff3->Draw("P");

  TLegend *Legend_002 = new TLegend(.50,.25,.85,.45);
  Legend_002->SetFillColor(0);
  Legend_002->SetBorderSize(5);
  Legend_002->SetTextSize(.04);
  Legend_002->AddEntry(Eff3,"LWO -  25 keV threshold","P");
  Legend_002->Draw();

  TF1 *fit_eff3 = new TF1("fit_eff3","[0] +[1]*tanh([2]*(x + [3]))",  0.,1200.);
  fit_eff3->SetParameters(51,49,0.006,-50.0);
  fit_eff3->SetLineColor(1);
  fit_eff3->SetLineWidth(2);
  fit_eff3->SetLineStyle(2); 
  Eff3->Fit("fit_eff3","R");

  /*
  TLatex *Latex_002_1 = new TLatex( 511.*Cal1/1000.,1500.,"511 keV");
  Latex_002_1->SetTextSize(.030);
  Latex_002_1->SetTextColor(1);
  Latex_002_1->SetTextAlign(22);
  Latex_002_1->Draw();
  */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
