/* Show_COV.C*/

{

/////////////////////////////////////////////////////////////////////////////
//                                                                         // 
// Open histogram files                                                    //
//                                                                         // 
/////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

  char file_name[ ] = "../../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36.root";
  TFile *f_01 = new TFile(file_name);

  //TFile *f_01 = new TFile("/sps/t2k/edoardo/COV/MyCOV1_BOT_4_4_7_TOP_10_10_15_NTD_00060.root");

  f_01->cd(); 

  Double_t RunTime = h1_RunTime->GetBinContent(1);
  cout << " RunTime = " << RunTime << "s" << endl;
  Double_t Gain_COV1 = h1_Gains->GetBinContent(1);
  Double_t Gain_COV2 = h1_Gains->GetBinContent(2);
  Double_t Threshold_COV1 = h1_Thresholds->GetBinContent(1);
  Double_t Threshold_COV2 = h1_Thresholds->GetBinContent(2);

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
  C_001->Divide(2,2);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  C_001->cd(1);
    
  h1_Raw2->GetXaxis()->SetTitle("T [10#mus]");
  h1_Raw2->GetYaxis()->SetTitle("Amplitude [ADC]");
  h1_Raw2->GetYaxis()->SetTitleOffset(1.5);
  h1_Raw2->SetLineColor(4);
  h1_Raw2->SetLineWidth(2);
  h1_Raw2->Draw();
  
  TLegend *Legend_001_1 = new TLegend(.15,.70,.35,.85);
  Legend_001_1->SetFillColor(0);
  Legend_001_1->SetBorderSize(0);
  Legend_001_1->AddEntry(h1_Raw2,"Ge BOT - 1","L");
  Legend_001_1->Draw();

  C_001->cd(2);
  h1_Raw4->GetXaxis()->SetTitle("T [10#mus]");
  h1_Raw4->GetYaxis()->SetTitle("Amplitude [ADC]");
  h1_Raw4->GetYaxis()->SetTitleOffset(1.5);
  h1_Raw4->SetLineColor(4);
  h1_Raw4->SetLineWidth(2);
  h1_Raw4->Draw();

  TLegend *Legend_001_2 = new TLegend(.15,.70,.35,.85);
  Legend_001_2->SetFillColor(0);
  Legend_001_2->SetBorderSize(0);
  Legend_001_2->AddEntry(h1_Raw4,"Ge TOP - 1","L");
  Legend_001_2->Draw();

  C_001->cd(3);
  h1_Raw1->GetXaxis()->SetTitle("T [10#mus]");
  h1_Raw1->GetYaxis()->SetTitle("Amplitude [ADC]");
  h1_Raw1->GetYaxis()->SetTitleOffset(1.5);
  h1_Raw1->SetLineColor(4);
  h1_Raw1->SetLineWidth(2);
  h1_Raw1->Draw();

  TLegend *Legend_001_3 = new TLegend(.15,.15,.35,.30);
  Legend_001_3->SetFillColor(0);
  Legend_001_3->SetBorderSize(0);
  Legend_001_3->AddEntry(h1_Raw1,"Ge BOT - 2","L");
  Legend_001_3->Draw();

  C_001->cd(4);
  h1_Raw3->GetXaxis()->SetTitle("T [10#mus]");
  h1_Raw3->GetYaxis()->SetTitle("Amplitude [ADC]");
  h1_Raw3->GetYaxis()->SetTitleOffset(1.5);
  h1_Raw3->SetLineColor(4);
  h1_Raw3->SetLineWidth(2);
  h1_Raw3->Draw();

  TLegend *Legend_001_4 = new TLegend(.15,.15,.35,.30);
  Legend_001_4->SetFillColor(0);
  Legend_001_4->SetBorderSize(0);
  Legend_001_4->AddEntry(h1_Raw3,"Ge TOP - 2","L");
  Legend_001_4->Draw();

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
//  C_002->Divide(2,2);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  C_002->cd();
  h1_COV1->GetXaxis()->SetTitle("T [10#mus]");
  h1_COV1->GetYaxis()->SetTitle("Amplitude [V]");
  h1_COV1->SetLineColor(4);
  h1_COV1->SetLineWidth(2);
  h1_COV1->Draw();

  TLegend *Legend_002_1 = new TLegend(.15,.70,.35,.85);
  Legend_002_1->SetFillColor(0);
  Legend_002_1->SetBorderSize(0);
  Legend_002_1->AddEntry(h1_COV1,"Ge BOT","L");
  Legend_002_1->Draw();

    TCanvas *C_002_b = new TCanvas("C_002b","C_002b: ", 1000,600);
  C_002_b->cd();
  h1_COV2->GetXaxis()->SetTitle("T [10#mus]");
  h1_COV2->GetYaxis()->SetTitle("Amplitude [V]");
  h1_COV2->SetLineColor(4);
  h1_COV2->SetLineWidth(2);
  h1_COV2->Draw();

  TLegend *Legend_002_2 = new TLegend(.15,.70,.35,.85);
  Legend_002_2->SetFillColor(0);
  Legend_002_2->SetBorderSize(0);
  Legend_002_2->AddEntry(h1_COV2,"Ge TOP","L");
  Legend_002_2->Draw();

    TCanvas *C_002_c = new TCanvas("C_002c","C_002c: ", 1000,600);
  C_002_c->cd();
  h1_DIF1->GetXaxis()->SetTitle("T [10#mus]");
  h1_DIF1->GetYaxis()->SetTitle("#DeltaAmplitude [V]");
  h1_DIF1->SetLineColor(4);
  h1_DIF1->SetLineWidth(2);
  h1_DIF1->Draw();

  TLegend *Legend_002_3 = new TLegend(.15,.70,.35,.85);
  Legend_002_3->SetFillColor(0);
  Legend_002_3->SetBorderSize(0);
  Legend_002_3->AddEntry(h1_DIF1,"Ge BOT","L");
  Legend_002_3->Draw();

    TCanvas *C_002_d = new TCanvas("C_002d","C_002d: ", 1000,600);
  C_002_d->cd();
  h1_DIF2->GetXaxis()->SetTitle("T [10#mus]");
  h1_DIF2->GetYaxis()->SetTitle("#DeltaAmplitude [V]");
  h1_DIF2->SetLineColor(4);
  h1_DIF2->SetLineWidth(2);
  h1_DIF2->Draw();

  TLegend *Legend_002_4 = new TLegend(.15,.70,.35,.85);
  Legend_002_4->SetFillColor(0);
  Legend_002_4->SetBorderSize(0);
  Legend_002_4->AddEntry(h1_DIF2,"Ge TOP","L");
  Legend_002_4->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Canvas 003bis :                                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
///


    TCanvas *C_002_time = new TCanvas("C_002time","C_002time: ", 1000,600);
    C_002_time->cd();
    TH1D* h1_COV2_b = (TH1D*) h1_COV2->Clone();
    TH1D* h1_COV2_sec = new TH1D("TOP_Ge_pulses_sec", "TOP_Ge_pulses_sec", h1_COV2_b->GetNbinsX(), 0, 1000);
    
    for (int i=1; i<h1_COV2_b->GetNbinsX()+1; i++){
        h1_COV2_sec->SetBinContent(i, h1_COV2_b->GetBinContent(i));
    }
    h1_COV2_sec->GetXaxis()->SetTitle("T [ms]");
    h1_COV2_sec->GetYaxis()->SetTitle("Amplitude [V]");
    h1_COV2_sec->SetLineColor(kGray+2);
    h1_COV2_sec->SetLineWidth(2);
    h1_COV2_sec->Draw();
    
    TH1D* h1_DIF2_b = (TH1D*) h1_DIF2->Clone();

    TH1D* h1_DIF2_sec = new TH1D("TOP_Ge_BtB_sec", "TOP_Ge_BtB_sec", h1_DIF2_b->GetNbinsX(), 0, 1000);
    
    for (int i=1; i<h1_DIF2_b->GetNbinsX()+1; i++){
        h1_DIF2_sec->SetBinContent(i, h1_DIF2_b->GetBinContent(i));
    }
    
    double Min_axis = -1.6;
    double Max_axis = 10.6;
    
    h1_DIF2_sec->GetXaxis()->SetTitle("T [ms]");
    h1_DIF2_sec->GetYaxis()->SetRangeUser(Min_axis, Max_axis);
    h1_DIF2_sec->GetYaxis()->SetTitle("#DeltaAmplitude [V]");
    h1_DIF2_sec->SetLineColor(kRed);
    h1_DIF2_sec->SetLineWidth(2);
    h1_DIF2_sec->Draw("SAME");
    
    TGaxis *axis = new TGaxis(1000, Min_axis+0.3, 1000, Max_axis, Min_axis+0.3, Max_axis, 510, "+L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetTitleColor(kRed);
    axis->SetTextFont(72);
    axis->SetTickLength(0.005);
    axis->SetLabelFont(h1_COV2_sec->GetYaxis()->GetLabelFont());
    axis->SetLabelSize(h1_COV2_sec->GetYaxis()->GetLabelSize());
    axis->SetTitleFont(h1_COV2_sec->GetYaxis()->GetTitleFont());
    axis->SetTitleSize(h1_COV2_sec->GetYaxis()->GetTitleSize());
    axis->SetTitle("#DeltaAmplitude [V]");
    axis->Draw("SAME");
    
    TLegend *Legend_002_time = new TLegend(.15,.70,.35,.85);
    Legend_002_time->SetFillColor(0);
    Legend_002_time->SetBorderSize(0);
    Legend_002_time->AddEntry(h1_COV2_sec,"Pulses TOP Ge","L");
    Legend_002_time->AddEntry(h1_DIF2_sec,"Bin to Bin TOP Ge","L");
    Legend_002_time->Draw("SAME");

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
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  h1_Peak1->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak1->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak1->GetXaxis()->SetRangeUser(  0.,  3.0);
  h1_Peak1->GetYaxis()->SetRangeUser(  1.0, 1000000.);
  h1_Peak1->SetLineColor(4);
  h1_Peak1->SetLineWidth(2);
  h1_Peak1->Draw();

  Double_t N_COV1 = h1_Peak1->GetEntries();

  TLegend *Legend_003 = new TLegend(.55,.65,.85,.85);
  Legend_003->SetFillColor(0);
  Legend_003->SetBorderSize(0);
  Legend_003->AddEntry(h1_Peak1,"Ge BOT","L");
  Legend_003->Draw();

  Double_t Cal1 = .950;

  TLatex *Latex_003_1 = new TLatex( 511.*Cal1/1000.,1500.,"511 keV");
  Latex_003_1->SetTextSize(.030);
  Latex_003_1->SetTextColor(1);
  Latex_003_1->SetTextAlign(22);
  Latex_003_1->Draw();

  TLatex *Latex_003_2 = new TLatex(1461.*Cal1/1000., 500.,"1461 keV");
  Latex_003_2->SetTextSize(.030);
  Latex_003_2->SetTextColor(1);
  Latex_003_2->SetTextAlign(22);
  Latex_003_2->Draw();

  TLatex *Latex_003_3 = new TLatex(2615.*Cal1/1000., 100.,"2615 keV");
  Latex_003_3->SetTextSize(.030);
  Latex_003_3->SetTextColor(1);
  Latex_003_3->SetTextAlign(22);
  Latex_003_3->Draw();

  TLatex *Latex_003_4 = new TLatex( 200.*Cal1/1000.,20000.,"60 keV");
  Latex_003_4->SetTextSize(.030);
  Latex_003_4->SetTextColor(1);
  Latex_003_4->SetTextAlign(22);
  Latex_003_4->Draw();

  cout << " Data entries = " << h1_Peak1->GetEntries() << " " << endl;
  cout << " Ge BOT rate (cps) = " << h1_Peak1->GetEntries()/RunTime << " +/- " << sqrt(h1_Peak1->GetEntries())/RunTime << " " << endl;

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

  TCanvas *C_004 = new TCanvas("C_004","C_004: ",1000, 600);
  C_004->Divide(1,1);
  C_004->SetLogy();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  h1_Peak2->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak2->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak2->GetXaxis()->SetRangeUser(  0.,      3.);
  h1_Peak2->GetYaxis()->SetRangeUser(  1., 100000.);
  h1_Peak2->SetLineColor(4);
  h1_Peak2->SetLineWidth(2);
  h1_Peak2->Draw();

  Double_t N_COV2 = h1_Peak2->GetEntries();

  TLegend *Legend_004 = new TLegend(.55,.65,.85,.85);
  Legend_004->SetFillColor(0);
  Legend_004->SetBorderSize(0);
  Legend_004->AddEntry(h1_Peak2,"Ge TOP","L");
  Legend_004->Draw();

  Double_t Cal2 = .650;

  TLatex *Latex_004_1 = new TLatex( 511.*Cal2/1000.,2000.,"511 keV");
  Latex_004_1->SetTextSize(.030);
  Latex_004_1->SetTextColor(1);
  Latex_004_1->SetTextAlign(22);
  Latex_004_1->Draw();

  TLatex *Latex_004_2 = new TLatex(1461.*Cal2/1000., 800.,"1461 keV");
  Latex_004_2->SetTextSize(.030);
  Latex_004_2->SetTextColor(1);
  Latex_004_2->SetTextAlign(22);
  Latex_004_2->Draw();

  TLatex *Latex_004_3 = new TLatex(2615.*Cal2/1000., 150.,"2615 keV");
  Latex_004_3->SetTextSize(.030);
  Latex_004_3->SetTextColor(1);
  Latex_004_3->SetTextAlign(22);
  Latex_004_3->Draw();

  cout << " Data entries = " << h1_Peak2->GetEntries() << " " << endl;
  cout << " Ge TOP rate (cps) = " << h1_Peak2->GetEntries()/RunTime << " +/- " << sqrt(h1_Peak2->GetEntries())/RunTime << " " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Canvas 005 :                                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 005 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *C_005 = new TCanvas("C_005","C_005: ", 1000,600);
  C_005->Divide(3,2);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Double_t Eg[3]         = {  511.00, 1460.83, 2614.53};
  Double_t COV1_Cal_M[3], COV1_Cal_S[3], COV1_Res_M[3], COV1_Res_S[3];
  Double_t COV2_Cal_M[3], COV2_Cal_S[3], COV2_Res_M[3], COV2_Res_S[3];

  TH1F *h1_Peak1_0511, *h1_Peak1_1461, *h1_Peak1_2615;
  TH1F *h1_Peak2_0511, *h1_Peak2_1461, *h1_Peak2_2615;

  C_005->cd(1);

  h1_Peak1_0511 = (TH1F*)h1_Peak1->Clone();
  h1_Peak1_0511->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak1_0511->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak1_0511->GetXaxis()->SetRangeUser(0.43,0.53);
  h1_Peak1_0511->GetYaxis()->SetRangeUser(0.,800.);
  h1_Peak1_0511->GetYaxis()->SetTitleOffset(1.5);
  h1_Peak1_0511->SetLineColor(4);
  h1_Peak1_0511->SetLineWidth(2);

  TF1 *fit_1_0511 = new TF1("fit_1_0511","expo+gaus(2)",0.43,0.53);

  fit_1_0511->SetParameters(7.,-4.,260.,0.486,0.004);
  h1_Peak1_0511->Fit("fit_1_0511","R");
  COV1_Cal_M[0] = fit_1_0511->GetParameter(3);
  COV1_Cal_S[0]=  fit_1_0511->GetParError(3);
  COV1_Res_M[0] = fit_1_0511->GetParameter(4);
  COV1_Res_S[0]=  fit_1_0511->GetParError(4);

  TLegend *Legend_005_1 = new TLegend(.65,.70,.85,.85);
  Legend_005_1->SetFillColor(0);
  Legend_005_1->SetBorderSize(0);
  Legend_005_1->AddEntry(h1_Peak2,"Ge BOT","L");
  Legend_005_1->Draw();

  TLatex *Latex_005_1 = new TLatex( .45, 640.,"511 keV");
  Latex_005_1->SetTextSize(.050);
  Latex_005_1->SetTextColor(1);
  Latex_005_1->SetTextAlign(22);
  Latex_005_1->Draw();

  C_005->cd(2);

  h1_Peak1_1461 = (TH1F*)h1_Peak1->Clone();
  h1_Peak1_1461->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak1_1461->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak1_1461->GetXaxis()->SetRangeUser(1.25,1.55);
  h1_Peak1_1461->GetYaxis()->SetRangeUser(  0.,200.);
  h1_Peak1_1461->GetYaxis()->SetTitleOffset(1.5);
  h1_Peak1_1461->SetLineColor(4);
  h1_Peak1_1461->SetLineWidth(2);

  TF1 *fit_1_1461 = new TF1("fit_1_1461","expo+gaus(2)",1.25,1.55);

  fit_1_1461->SetParameters(10.,-4.,120.,1.385,0.010);
  h1_Peak1_1461->Fit("fit_1_1461","R");
  COV1_Cal_M[1] = fit_1_1461->GetParameter(3);
  COV1_Cal_S[1]=  fit_1_1461->GetParError(3);
  COV1_Res_M[1] = fit_1_1461->GetParameter(4);
  COV1_Res_S[1]=  fit_1_1461->GetParError(4);

  TLegend *Legend_005_2 = new TLegend(.65,.70,.85,.85);
  Legend_005_2->SetFillColor(0);
  Legend_005_2->SetBorderSize(0);
  Legend_005_2->AddEntry(h1_Peak2,"Ge BOT","L");
  Legend_005_2->Draw();

  TLatex *Latex_005_2 = new TLatex(1.31, 160.,"1461 keV");
  Latex_005_2->SetTextSize(.050);
  Latex_005_2->SetTextColor(1);
  Latex_005_2->SetTextAlign(22);
  Latex_005_2->Draw();

  C_005->cd(3);

  h1_Peak1_2615 = (TH1F*)h1_Peak1->Clone();
  h1_Peak1_2615->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak1_2615->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak1_2615->GetXaxis()->SetRangeUser(2.30,2.60);
  h1_Peak1_2615->GetYaxis()->SetRangeUser(  0., 60.);
  h1_Peak1_2615->GetYaxis()->SetTitleOffset(1.5);
  h1_Peak1_2615->SetLineColor(4);
  h1_Peak1_2615->SetLineWidth(2);

  TF1 *fit_1_2615 = new TF1("fit_1_2615","expo+gaus(2)",2.30,2.60);

  fit_1_2615->SetParameters(4.,-2.,35.,2.490,0.018);
  h1_Peak1_2615->Fit("fit_1_2615","R");
  COV1_Cal_M[2] = fit_1_2615->GetParameter(3);
  COV1_Cal_S[2]=  fit_1_2615->GetParError(3);
  COV1_Res_M[2] = fit_1_2615->GetParameter(4);
  COV1_Res_S[2]=  fit_1_2615->GetParError(4);

  TLegend *Legend_005_3 = new TLegend(.65,.70,.85,.85);
  Legend_005_3->SetFillColor(0);
  Legend_005_3->SetBorderSize(0);
  Legend_005_3->AddEntry(h1_Peak2,"Ge BOT","L");
  Legend_005_3->Draw();

  TLatex *Latex_005_3 = new TLatex(2.36,  48.,"2615 keV");
  Latex_005_3->SetTextSize(.050);
  Latex_005_3->SetTextColor(1);
  Latex_005_3->SetTextAlign(22);
  Latex_005_3->Draw();

  C_005->cd(4);

  h1_Peak2_0511 = (TH1F*)h1_Peak2->Clone();
  h1_Peak2_0511->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak2_0511->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak2_0511->GetXaxis()->SetRangeUser(0.28,0.36);
  h1_Peak2_0511->GetYaxis()->SetRangeUser(  0.,1400.);
  h1_Peak2_0511->GetYaxis()->SetTitleOffset(1.5);
  h1_Peak2_0511->SetLineColor(4);
  h1_Peak2_0511->SetLineWidth(2);

  TF1 *fit_2_0511 = new TF1("fit_2_0511","expo+gaus(2)",0.28,0.36);

  fit_2_0511->SetParameters(7.,-4.,260.,0.325,0.004);
  h1_Peak2_0511->Fit("fit_2_0511","R");
  COV2_Cal_M[0] = fit_2_0511->GetParameter(3);
  COV2_Cal_S[0]=  fit_2_0511->GetParError(3);
  COV2_Res_M[0] = fit_2_0511->GetParameter(4);
  COV2_Res_S[0]=  fit_2_0511->GetParError(4);

  TLegend *Legend_005_4 = new TLegend(.65,.70,.85,.85);
  Legend_005_4->SetFillColor(0);
  Legend_005_4->SetBorderSize(0);
  Legend_005_4->AddEntry(h1_Peak2,"Ge TOP","L");
  Legend_005_4->Draw();

  TLatex *Latex_005_4 = new TLatex(0.296,1120.,"511 keV");
  Latex_005_4->SetTextSize(.050);
  Latex_005_4->SetTextColor(1);
  Latex_005_4->SetTextAlign(22);
  Latex_005_4->Draw();

  C_005->cd(5);

  h1_Peak2_1461 = (TH1F*)h1_Peak2->Clone();
  h1_Peak2_1461->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak2_1461->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak2_1461->GetXaxis()->SetRangeUser(0.82,1.02);
  h1_Peak2_1461->GetYaxis()->SetRangeUser(  0.,400.);
  h1_Peak2_1461->GetYaxis()->SetTitleOffset(1.5);
  h1_Peak2_1461->SetLineColor(4);
  h1_Peak2_1461->SetLineWidth(2);

  TF1 *fit_2_1461 = new TF1("fit_2_1461","expo+gaus(2)",0.82,1.02);

  fit_2_1461->SetParameters( 9.,-5.,215.,0.948,0.006);
  h1_Peak2_1461->Fit("fit_2_1461","R");
  COV2_Cal_M[1] = fit_2_1461->GetParameter(3);
  COV2_Cal_S[1]=  fit_2_1461->GetParError(3);
  COV2_Res_M[1] = fit_2_1461->GetParameter(4);
  COV2_Res_S[1]=  fit_2_1461->GetParError(4);

  TLegend *Legend_005_5 = new TLegend(.65,.70,.85,.85);
  Legend_005_5->SetFillColor(0);
  Legend_005_5->SetBorderSize(0);
  Legend_005_5->AddEntry(h1_Peak2,"Ge TOP","L");
  Legend_005_5->Draw();

  TLatex *Latex_005_5 = new TLatex(0.86, 320.,"1461 keV");
  Latex_005_5->SetTextSize(.050);
  Latex_005_5->SetTextColor(1);
  Latex_005_5->SetTextAlign(22);
  Latex_005_5->Draw();

  C_005->cd(6);

  h1_Peak2_2615 = (TH1F*)h1_Peak2->Clone();
  h1_Peak2_2615->GetXaxis()->SetTitle("Amplitude [V]");
  h1_Peak2_2615->GetYaxis()->SetTitle("Counts / 2 mV");
  h1_Peak2_2615->GetXaxis()->SetRangeUser(1.55,1.80);
  h1_Peak2_2615->GetYaxis()->SetRangeUser(  0.,100.);
  h1_Peak2_2615->GetYaxis()->SetTitleOffset(1.5);
  h1_Peak2_2615->SetLineColor(4);
  h1_Peak2_2615->SetLineWidth(2);

  TF1 *fit_2_2615 = new TF1("fit_2_2615","expo+gaus(2)",1.55,1.80);

  fit_2_2615->SetParameters(4.,-2.,50.,1.700,0.012);
  h1_Peak2_2615->Fit("fit_2_2615","R");
  COV2_Cal_M[2] = fit_2_2615->GetParameter(3);
  COV2_Cal_S[2]=  fit_2_2615->GetParError(3);
  COV2_Res_M[2] = fit_2_2615->GetParameter(4);
  COV2_Res_S[2]=  fit_2_2615->GetParError(4);

  TLegend *Legend_005_6 = new TLegend(.65,.70,.85,.85);
  Legend_005_6->SetFillColor(0);
  Legend_005_6->SetBorderSize(0);
  Legend_005_6->AddEntry(h1_Peak2,"Ge TOP","L");
  Legend_005_6->Draw();

  TLatex *Latex_005_6 = new TLatex(1.60,  80.,"2615 keV");
  Latex_005_6->SetTextSize(.050);
  Latex_005_6->SetTextColor(1);
  Latex_005_6->SetTextAlign(22);
  Latex_005_6->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Canvas 006 :                                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

  cout << " ++++++++++++++++++++++++ " << endl;
  cout << " ++++++ Canvas 006 ++++++ " << endl;
  cout << " ++++++++++++++++++++++++ " << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *C_006 = new TCanvas("C_006","C_006: ", 1000,600);
  C_006->Divide(1,1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(000000);
  gStyle->SetOptStat("000000001");

  TH1F *Frame_006 = C_006->DrawFrame(0., 0., 3000., 3.);

  Frame_006->GetXaxis()->SetTitle("E [keV]");
  Frame_006->GetYaxis()->SetTitle("Amplitude [V]");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Int_t nn1 = 3;
  Double_t xx1[10], yy1[10], ex1[10], ey1[10];

  for (Int_t j = 0; j < nn1; j++) {
    xx1[j] = Eg[j];
    ex1[j] =    0.;
    yy1[j] = COV1_Cal_M[j];
    ey1[j] = COV1_Cal_S[j];
  }

  TGraphErrors *COV1_Cal = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
  COV1_Cal->SetMarkerColor(1);
  COV1_Cal->SetMarkerStyle(20);
  COV1_Cal->SetMarkerSize(1.4);
  COV1_Cal->Draw("P");

  TF1 *fit_01 = new TF1("fit_01","pol1",0.,3000. );
  fit_01->SetRange(0.,3000.);
  fit_01->SetParameters(0.,.001);
  fit_01->SetLineStyle(2);

  COV1_Cal->Fit("fit_01","R");

  Double_t COV1_Conv = fit_01->GetParameter(1);

  for (Int_t j = 0; j < nn1; j++) {
    xx1[j] = Eg[j];
    ex1[j] =    0.;
    yy1[j] = COV2_Cal_M[j];
    ey1[j] = COV2_Cal_S[j];
  }

  TGraphErrors *COV2_Cal = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
  COV2_Cal->SetMarkerColor(4);
  COV2_Cal->SetMarkerStyle(20);
  COV2_Cal->SetMarkerSize(1.4);
  COV2_Cal->Draw("P");

  TF1 *fit_02 = new TF1("fit_02","pol1",0.,3000. );
  fit_02->SetRange(0.,3000.);
  fit_02->SetParameters(0.,.001);
  fit_02->SetLineStyle(2);

  COV2_Cal->Fit("fit_02","R");

  Double_t COV2_Conv = fit_02->GetParameter(1);

  TLegend *Legend_006 = new TLegend(.15,.60,.45,.80);
  Legend_006->SetFillColor(0);
  Legend_006->SetBorderSize(0);
  Legend_006->AddEntry(COV1_Cal,"Ge BOT","P");
  Legend_006->AddEntry(COV2_Cal,"Ge TOP","P");
  Legend_006->Draw();

  TLatex *Latex_006_1 = new TLatex(2000.,2.5,"953 #muV/keV");
  Latex_006_1->SetTextSize(.060);
  Latex_006_1->SetTextColor(1);
  Latex_006_1->SetTextAlign(22);
  Latex_006_1->Draw();

  TLatex *Latex_006_2 = new TLatex(2500.,1.1,"651 #muV/keV");
  Latex_006_2->SetTextSize(.060);
  Latex_006_2->SetTextColor(1);
  Latex_006_2->SetTextAlign(22);
  Latex_006_2->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//// Canvas 007 :                                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 007 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_007 = new TCanvas("C_007","C_007: ", 1000,600);
//  C_007->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
//  TH1F *Frame_007 = C_007->DrawFrame(0., 0., 3000., 1.6);
//
//  Frame_007->GetXaxis()->SetTitle("E [keV]");
//  Frame_007->GetYaxis()->SetTitle("#sigma_{E} / E [%]");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  for (Int_t j = 0; j < nn1; j++) {
//    xx1[j] = Eg[j];
//    ex1[j] =    0.;
//    yy1[j] = 2.*sqrt(log(4.))*COV1_Res_M[j]/COV1_Conv;
//    ey1[j] = 2.*sqrt(log(4.))*COV1_Res_S[j]/COV1_Conv;
//    yy1[j] = 100.*abs(COV1_Res_M[j]/COV1_Conv/Eg[j]);
//    ey1[j] = 100.*COV1_Res_S[j]/COV1_Conv/Eg[j];
//  }
//
//  TGraphErrors *COV1_Res = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
//  COV1_Res->SetMarkerColor(1);
//  COV1_Res->SetMarkerStyle(20);
//  COV1_Res->SetMarkerSize(1.4);
//  COV1_Res->Draw("P");
//
//  for (Int_t j = 0; j < nn1; j++) {
//    xx1[j] = Eg[j];
//    ex1[j] =    0.;
//    yy1[j] = 2.*log(4.)*COV2_Res_M[j]/COV2_Conv;
//    ey1[j] = 2.*log(4.)*COV2_Res_S[j]/COV2_Conv;
//    yy1[j] = 100.*COV2_Res_M[j]/COV2_Conv/Eg[j];
//    ey1[j] = 100.*COV2_Res_S[j]/COV2_Conv/Eg[j];
//  }
//
//  TGraphErrors *COV2_Res = new TGraphErrors(nn1,xx1,yy1,ex1,ey1);
//  COV2_Res->SetMarkerColor(4);
//  COV2_Res->SetMarkerStyle(20);
//  COV2_Res->SetMarkerSize(1.4);
//  COV2_Res->Draw("P");
//
//  TLegend *Legend_007 = new TLegend(.55,.60,.85,.80);
//  Legend_007->SetFillColor(0);
//  Legend_007->SetBorderSize(0);
//  Legend_007->AddEntry(COV1_Res,"Ge BOT","P");
//  Legend_007->AddEntry(COV2_Res,"Ge TOP","P");
//  Legend_007->Draw();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 008 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 008 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_008 = new TCanvas("C_008","C_008: ", 1000,600);
//  C_008->Divide(2,2);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TH1F *h1_COV1_Zoom, *h1_COV2_Zoom, *h1_DIF1_Zoom, *h1_DIF2_Zoom;
//  h1_COV1_Zoom = (TH1F*)h1_COV1->Clone();
//  h1_COV2_Zoom = (TH1F*)h1_COV2->Clone();
//  h1_DIF1_Zoom = (TH1F*)h1_DIF1->Clone();
//  h1_DIF2_Zoom = (TH1F*)h1_DIF2->Clone();
//
//
//  C_008->cd(1);
//  h1_COV1_Zoom->GetXaxis()->SetTitle("T [10#mus]");
//  h1_COV1_Zoom->GetYaxis()->SetTitle("Amplitude [V]");
//  h1_COV1_Zoom->GetXaxis()->SetRangeUser(36000.,48000.);
//  h1_COV1_Zoom->SetLineColor(4);
//  h1_COV1_Zoom->SetLineWidth(2);
//  h1_COV1_Zoom->Draw();
//
//  TLegend *Legend_008_1 = new TLegend(.15,.70,.35,.85);
//  Legend_008_1->SetFillColor(0);
//  Legend_008_1->SetBorderSize(0);
//  Legend_008_1->AddEntry(h1_COV2_Zoom,"Ge BOT","L");
//  Legend_008_1->Draw();
//
//  C_008->cd(2);
//  h1_COV2_Zoom->GetXaxis()->SetTitle("T [10#mus]");
//  h1_COV2_Zoom->GetYaxis()->SetTitle("Amplitude [V]");
//  h1_COV2_Zoom->GetXaxis()->SetRangeUser(70000.,85000.);
//  h1_COV2_Zoom->SetLineColor(4);
//  h1_COV2_Zoom->SetLineWidth(2);
//  h1_COV2_Zoom->Draw();
//
//  TLegend *Legend_008_2 = new TLegend(.65,.70,.85,.85);
//  Legend_008_2->SetFillColor(0);
//  Legend_008_2->SetBorderSize(0);
//  Legend_008_2->AddEntry(h1_COV2_Zoom,"Ge TOP","L");
//  Legend_008_2->Draw();
//
//  C_008->cd(3);
//  h1_DIF1_Zoom->GetXaxis()->SetTitle("T [10#mus]");
//  h1_DIF1_Zoom->GetYaxis()->SetTitle("#DeltaAmplitude [V]");
//  h1_DIF1_Zoom->GetXaxis()->SetRangeUser(40050.,40150.);
//  h1_DIF1_Zoom->SetLineColor(4);
//  h1_DIF1_Zoom->SetLineWidth(2);
//  h1_DIF1_Zoom->Draw();
//
//  TLegend *Legend_008_3 = new TLegend(.15,.70,.35,.85);
//  Legend_008_3->SetFillColor(0);
//  Legend_008_3->SetBorderSize(0);
//  Legend_008_3->AddEntry(h1_COV2,"Ge BOT","L");
//  Legend_008_3->Draw();
//
//  C_008->cd(4);
//  h1_DIF2_Zoom->GetXaxis()->SetTitle("T [10#mus]");
//  h1_DIF2_Zoom->GetYaxis()->SetTitle("#DeltaAmplitude [V]");
//  h1_DIF2_Zoom->GetXaxis()->SetRangeUser(70000.,85000.);
//  h1_DIF2_Zoom->SetLineColor(4);
//  h1_DIF2_Zoom->SetLineWidth(2);
//  h1_DIF2_Zoom->Draw();
//
//  TLegend *Legend_008_4 = new TLegend(.65,.70,.85,.85);
//  Legend_008_4->SetFillColor(0);
//  Legend_008_4->SetBorderSize(0);
//  Legend_008_4->AddEntry(h1_COV2,"Ge TOP","L");
//  Legend_008_4->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 011 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 011 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_011 = new TCanvas("C_011","C_011: ", 1000,600);
//  C_011->Divide(1,1);
//  C_011->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak1_keV->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak1_keV->GetYaxis()->SetTitle("Counts / bin");
//  h1_Peak1_keV->GetXaxis()->SetRangeUser(  0.,  3000.0);
//  h1_Peak1_keV->GetYaxis()->SetRangeUser(  1.0, 1000000.);
//  h1_Peak1_keV->SetLineColor(4);
//  h1_Peak1_keV->SetLineWidth(2);
//  h1_Peak1_keV->Draw();
//
//  TLegend *Legend_011 = new TLegend(.55,.65,.85,.85);
//  Legend_011->SetFillColor(0);
//  Legend_011->SetBorderSize(0);
//  Legend_011->AddEntry(h1_Peak1,"Ge BOT","L");
//  Legend_011->Draw();
//
//  TLatex *Latex_011_1 = new TLatex( 511.,1500.,"511 keV");
//  Latex_011_1->SetTextSize(.030);
//  Latex_011_1->SetTextColor(1);
//  Latex_011_1->SetTextAlign(22);
//  Latex_011_1->Draw();
//
//  TLatex *Latex_011_2 = new TLatex(1461., 500.,"1461 keV");
//  Latex_011_2->SetTextSize(.030);
//  Latex_011_2->SetTextColor(1);
//  Latex_011_2->SetTextAlign(22);
//  Latex_011_2->Draw();
//
//  TLatex *Latex_011_3 = new TLatex(2615., 100.,"2615 keV");
//  Latex_011_3->SetTextSize(.030);
//  Latex_011_3->SetTextColor(1);
//  Latex_011_3->SetTextAlign(22);
//  Latex_011_3->Draw();
//
//  TLatex *Latex_011_4 = new TLatex( 200.,20000.,"60 keV");
//  Latex_011_4->SetTextSize(.030);
//  Latex_011_4->SetTextColor(1);
//  Latex_011_4->SetTextAlign(22);
//  Latex_011_4->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 012 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 012 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_012 = new TCanvas("C_012","C_012: ", 1000,600);
//  C_012->Divide(1,1);
//  C_012->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak2_keV->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak2_keV->GetYaxis()->SetTitle("Counts / bin");
//  h1_Peak2_keV->GetXaxis()->SetRangeUser(  0.,   3000.);
//  h1_Peak2_keV->GetYaxis()->SetRangeUser(  1., 100000.);
//  h1_Peak2_keV->SetLineColor(4);
//  h1_Peak2_keV->SetLineWidth(2);
//  h1_Peak2_keV->Draw();
//
//  TLegend *Legend_012 = new TLegend(.55,.65,.85,.85);
//  Legend_012->SetFillColor(0);
//  Legend_012->SetBorderSize(0);
//  Legend_012->AddEntry(h1_Peak2,"Ge TOP","L");
//  Legend_012->Draw();
//
//  TLatex *Latex_012_1 = new TLatex( 511.,2000.,"511 keV");
//  Latex_012_1->SetTextSize(.030);
//  Latex_012_1->SetTextColor(1);
//  Latex_012_1->SetTextAlign(22);
//  Latex_012_1->Draw();
//
//  TLatex *Latex_012_2 = new TLatex(1461., 800.,"1461 keV");
//  Latex_012_2->SetTextSize(.030);
//  Latex_012_2->SetTextColor(1);
//  Latex_012_2->SetTextAlign(22);
//  Latex_012_2->Draw();
//
//  TLatex *Latex_012_3 = new TLatex(2615., 150.,"2615 keV");
//  Latex_012_3->SetTextSize(.030);
//  Latex_012_3->SetTextColor(1);
//  Latex_012_3->SetTextAlign(22);
//  Latex_012_3->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 013 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 013 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_013 = new TCanvas("C_013","C_013: ", 1000,600);
//  C_013->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TH1F *h1_Peak1_keV_Am;
//  h1_Peak1_keV_Am = (TH1F*)h1_Peak1_keV->Clone();
//
//  h1_Peak1_keV_Am->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak1_keV_Am->GetYaxis()->SetTitle("Counts / bin");
//  h1_Peak1_keV_Am->GetXaxis()->SetRangeUser(  0.,  100.0);
//  h1_Peak1_keV_Am->GetYaxis()->SetRangeUser(  0., 140000.);
//  h1_Peak1_keV_Am->SetLineColor(4);
//  h1_Peak1_keV_Am->SetLineWidth(2);
//  h1_Peak1_keV_Am->Draw();
//
//  TF1 *fit_Am = new TF1("fit_Am","pol1+gaus(2)",45.,75.);
//
//  fit_Am->SetParameters( 9.,-5.,69000.,59.5,3.0);
//  h1_Peak1_keV_Am->Fit("fit_Am","R");
//
//  TLegend *Legend_013 = new TLegend(.15,.65,.45,.85);
//  Legend_013->SetFillColor(0);
//  Legend_013->SetBorderSize(0);
//  Legend_013->AddEntry(h1_Peak1_keV_Am,"Ge BOT","L");
//  Legend_013->Draw();
//
//  TLatex *Latex_013 = new TLatex( 82.,35000.,"FWHM = 6.8 keV");
//  Latex_013->SetTextSize(.05);
//  Latex_013->SetTextColor(1);
//  Latex_013->SetTextAlign(22);
//  Latex_013->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 014 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 014 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_014 = new TCanvas("C_014","C_014: ", 1000,600);
//  C_014->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak1_100keV->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak1_100keV->GetYaxis()->SetTitle("Counts / bin");
//  h1_Peak1_100keV->GetXaxis()->SetRangeUser(  0.,  100.0);
//  h1_Peak1_100keV->GetYaxis()->SetRangeUser(  0., 30000.);
//  h1_Peak1_100keV->SetLineColor(4);
//  h1_Peak1_100keV->SetLineWidth(2);
//  h1_Peak1_100keV->Draw();
//
//  TF1 *fit_100keV = new TF1("fit_100keV","pol1+gaus(2)",50.,70.);
//
//  fit_100keV->SetParameters( 9.,-5.,69000.,59.5,3.0);
//  h1_Peak1_100keV->Fit("fit_100keV","R");
//
//  Double_t aa = fit_100keV->GetParameter(0);
//  Double_t bb = fit_100keV->GetParameter(1);
//
//  Double_t Sum = 0.;
//  for (Int_t j = 0; j < 300; j++) {
//    if (abs(h1_Peak1_100keV->GetXaxis()->GetBinCenter(j + 1) - fit_100keV->GetParameter(3)) < 4.*fit_100keV->GetParameter(4)) {
//      Sum = Sum + h1_Peak1_100keV->GetBinContent(j + 1) - aa - bb*h1_Peak1_100keV->GetXaxis()->GetBinCenter(j + 1);
//    }
//  }
//  cout << " !!!!!!!!!!!!! " << Sum << " " << endl;
//
//  TLegend *Legend_014 = new TLegend(.15,.65,.45,.85);
//  Legend_014->SetFillColor(0);
//  Legend_014->SetBorderSize(0);
//  Legend_014->AddEntry(h1_Peak1_100keV,"Ge BOT","L");
//  Legend_014->Draw();
//
//  TLatex *Latex_014 = new TLatex( 82.,15000.,"FWHM = 3.6 keV");
//  Latex_014->SetTextSize(.05);
//  Latex_014->SetTextColor(1);
//  Latex_014->SetTextAlign(22);
//  Latex_014->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 015 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 015 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_015 = new TCanvas("C_015","C_015: ", 1000,600);
//  C_015->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak1_100keV_AVG->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak1_100keV_AVG->GetYaxis()->SetTitle("Counts / bin");
//  h1_Peak1_100keV_AVG->GetXaxis()->SetRangeUser(  0.,  100.0);
//  h1_Peak1_100keV_AVG->GetYaxis()->SetRangeUser(  0., 30000.);
//  h1_Peak1_100keV_AVG->SetLineColor(4);
//  h1_Peak1_100keV_AVG->SetLineWidth(2);
//  h1_Peak1_100keV_AVG->DrawCopy("SAME");
//
//  TF1 *fit_100keV_AVG = new TF1("fit_100keV_AVG","pol1+gaus(2)",50.,70.);
//
//  fit_100keV_AVG->SetParameters( 9.,-5.,69000.,59.5,3.0);
//  h1_Peak1_100keV_AVG->Fit("fit_100keV_AVG","R");
//
//  aa = fit_100keV_AVG->GetParameter(0);
//  bb = fit_100keV_AVG->GetParameter(1);
//
//  Sum = 0.;
//  for (Int_t j = 0; j < 300; j++) {
//    if (abs(h1_Peak1_100keV_AVG->GetXaxis()->GetBinCenter(j + 1) - fit_100keV_AVG->GetParameter(3)) < 4.*fit_100keV_AVG->GetParameter(4)) {
//      Sum = Sum + h1_Peak1_100keV_AVG->GetBinContent(j + 1) - aa - bb*h1_Peak1_100keV_AVG->GetXaxis()->GetBinCenter(j + 1);
//    }
//  }
//  cout << " !!!!!!!!!!!!! " << Sum << " " << endl;
//  Double_t N_COV1_60keV = Sum;
//
//  h1_Peak1_100keV_1kHz_AVG->SetLineColor(3);
//  h1_Peak1_100keV_1kHz_AVG->SetLineWidth(2);
//  h1_Peak1_100keV_1kHz_AVG->DrawCopy("SAME");;
//
//  TLegend *Legend_015 = new TLegend(.15,.65,.45,.85);
//  Legend_015->SetFillColor(0);
//  Legend_015->SetBorderSize(0);
//  Legend_015->AddEntry(h1_Peak1_100keV_AVG,"Ge BOT","L");
//  Legend_015->AddEntry(h1_Peak1_100keV_1kHz_AVG,"1 kHz X-Talk","L");
//  Legend_015->Draw();
//
//  TLatex *Latex_015 = new TLatex( 82.,15000.,"FWHM = 3.5 keV");
//  Latex_015->SetTextSize(.05);
//  Latex_015->SetTextColor(1);
//  Latex_015->SetTextAlign(22);
//  Latex_015->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 016 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 016 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_016 = new TCanvas("C_016","C_016: ", 1000,600);
//  C_016->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak2_100keV_AVG->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak2_100keV_AVG->GetYaxis()->SetTitle("Counts / bin");
//  h1_Peak2_100keV_AVG->GetXaxis()->SetRangeUser(  0.,  100.);
//  h1_Peak2_100keV_AVG->GetYaxis()->SetRangeUser(  0.,10000.);
//  h1_Peak2_100keV_AVG->SetLineColor(4);
//  h1_Peak2_100keV_AVG->SetLineWidth(2);
//  h1_Peak2_100keV_AVG->Draw();
//
//  //TF1 *fit_100keV_AVG = new TF1("fit_100keV_AVG","pol1+gaus(2)",50.,70.);
//  /*
//  fit_100keV_AVG->SetParameters( 9.,-5.,69000.,59.5,3.0);
//  h1_Peak1_100keV_AVG->Fit("fit_100keV_AVG","R");
//
//  aa = fit_100keV_AVG->GetParameter(0);
//  bb = fit_100keV_AVG->GetParameter(1);
//
//  Sum = 0.;
//  for (Int_t j = 0; j < 300; j++) {
//    if (abs(h1_Peak1_100keV_AVG->GetXaxis()->GetBinCenter(j + 1) - fit_100keV_AVG->GetParameter(3)) < 4.*fit_100keV_AVG->GetParameter(4)) {
//      Sum = Sum + h1_Peak1_100keV_AVG->GetBinContent(j + 1) - aa - bb*h1_Peak1_100keV_AVG->GetXaxis()->GetBinCenter(j + 1);
//    }
//  }
//  cout << " !!!!!!!!!!!!! " << Sum << " " << endl;
//  */
//
//  h1_Peak2_100keV_1kHz_AVG->SetLineColor(3);
//  h1_Peak2_100keV_1kHz_AVG->SetLineWidth(2);
//  h1_Peak2_100keV_1kHz_AVG->DrawCopy("SAME");
//
//  TLegend *Legend_016 = new TLegend(.55,.65,.85,.85);
//  Legend_016->SetFillColor(0);
//  Legend_016->SetBorderSize(0);
//  Legend_016->AddEntry(h1_Peak2_100keV_AVG,"Ge TOP","L");
//  Legend_016->AddEntry(h1_Peak2_100keV_1kHz_AVG,"1 kHz X-Talk","L");
//  Legend_016->Draw();
//  /*
//  TLatex *Latex_016 = new TLatex( 82.,15000.,"FWHM = 3.3 keV");
//  Latex_016->SetTextSize(.05);
//  Latex_016->SetTextColor(1);
//  Latex_016->SetTextAlign(22);
//  Latex_016->Draw();
//  */
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 017 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 017 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_017 = new TCanvas("C_017","C_017: ", 1000,600);
//  C_017->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_BSLS1->GetXaxis()->SetTitle("BSL RMS [keV]");
//  h1_BSLS1->GetYaxis()->SetTitle("Counts / bin");
//  h1_BSLS1->GetXaxis()->SetRangeUser(  0.,       5.);
//  h1_BSLS1->GetYaxis()->SetRangeUser(  0., 6000000.);
//  h1_BSLS1->SetLineColor(4);
//  h1_BSLS1->SetLineWidth(2);
//  h1_BSLS1->Draw();
//
//  Double_t COV1_BSL_RMS = h1_BSLS1->GetMean();
//
//  TLegend *Legend_017 = new TLegend(.55,.65,.85,.85);
//  Legend_017->SetFillColor(0);
//  Legend_017->SetBorderSize(0);
//  Legend_017->AddEntry(h1_BSLS1,"Ge BOT","L");
//  Legend_017->Draw();
//
//  TLatex *Latex_017 = new TLatex(2.8,3000000.,"<RMS> = 1.4 keV");
//  Latex_017->SetTextSize(.05);
//  Latex_017->SetTextColor(1);
//  Latex_017->SetTextAlign(22);
//  Latex_017->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 018 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 018 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_018 = new TCanvas("C_018","C_018: ", 1000,600);
//  C_018->Divide(1,1);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_BSLS2->GetXaxis()->SetTitle("BSL RMS [keV]");
//  h1_BSLS2->GetYaxis()->SetTitle("Counts / bin");
//  h1_BSLS2->GetXaxis()->SetRangeUser(  0.,       5.);
//  h1_BSLS2->GetYaxis()->SetRangeUser(  0., 6000000.);
//  h1_BSLS2->SetLineColor(4);
//  h1_BSLS2->SetLineWidth(2);
//  h1_BSLS2->Draw();
//
//  Double_t COV2_BSL_RMS = h1_BSLS2->GetMean();
//
//  TLegend *Legend_018 = new TLegend(.55,.65,.85,.85);
//  Legend_018->SetFillColor(0);
//  Legend_018->SetBorderSize(0);
//  Legend_018->AddEntry(h1_BSLS2,"Ge TOP","L");
//  Legend_018->Draw();
//
//  TLatex *Latex_018 = new TLatex(3.3,3000000.,"<RMS> = 1.96 keV");
//  Latex_018->SetTextSize(.05);
//  Latex_018->SetTextColor(1);
//  Latex_018->SetTextAlign(22);
//  Latex_018->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 020 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 020 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_020 = new TCanvas("C_020","C_020: ", 1000,600);
//  C_020->Divide(1,1);
//  C_020->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak1_4MeV_AVG->GetXaxis()->SetTitle("E [MeV]");
//  h1_Peak1_4MeV_AVG->GetYaxis()->SetTitle("Counts / 2 keV");
//  h1_Peak1_4MeV_AVG->GetXaxis()->SetRangeUser(0.,      4.);
//  h1_Peak1_4MeV_AVG->GetYaxis()->SetRangeUser(1.,1000000.);
//  h1_Peak1_4MeV_AVG->SetLineColor(4);
//  h1_Peak1_4MeV_AVG->SetLineWidth(1);
//  h1_Peak1_4MeV_AVG->Draw();
//
//  Double_t mysum = 0.;
//  for (Int_t j = 0; j < 40; j++) {mysum = mysum + h1_Peak1_4MeV_AVG->GetBinContent(j + 1);}
//  cout << " zzzzzzzzzz " << mysum << " " << endl;
//  mysum = 0.;
//  for (Int_t j = 0; j < 40; j++) {mysum = mysum + h1_Peak2_4MeV_AVG->GetBinContent(j + 1);}
//  cout << " zzzzzzzzzz " << mysum << " " << endl;
//
//  h1_Peak2_4MeV_AVG->SetLineColor(2);
//  //h1_Peak2_4MeV->SetFillColor(2);
//  h1_Peak2_4MeV_AVG->SetLineWidth(2);
//  h1_Peak2_4MeV_AVG->DrawCopy("SAME");
//
//  TLegend *Legend_020 = new TLegend(.65,.65,.85,.85);
//  Legend_020->SetFillColor(0);
//  Legend_020->SetBorderSize(0);
//  Legend_020->AddEntry(h1_Peak1_4MeV_AVG,"Ge BOT","L");
//  Legend_020->AddEntry(h1_Peak2_4MeV_AVG,"Ge TOP","L");
//  Legend_020->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 021 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 021 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_021 = new TCanvas("C_021","C_021: ", 1000,600);
//  C_021->Divide(1,1);
//  C_021->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak1_MeV->GetXaxis()->SetTitle("E [MeV]");
//  h1_Peak1_MeV->GetYaxis()->SetTitle("Counts / 10 keV");
//  h1_Peak1_MeV->GetXaxis()->SetRangeUser(0.,     20.);
//  h1_Peak1_MeV->GetYaxis()->SetRangeUser(1.,1000000.);
//  h1_Peak1_MeV->SetLineColor(4);
//  h1_Peak1_MeV->SetLineWidth(2);
//  h1_Peak1_MeV->Draw();
//
//  Int_t  N_COV1_T = h1_Peak1_MeV->GetEntries();
//  cout << " ****** " << N_COV1_T << " " << endl;
//  Int_t N_COV1_1 = N_COV1_T;
//  for (int j =  0; j < 3; j++) {
//    N_COV1_1 = N_COV1_1 - h1_Peak1_MeV->GetBinContent(j + 1);
//    cout << " " << j << " " << h1_Peak1_MeV->GetBinContent(j + 1) << " " << N_COV1_1 << " " << endl;
//  }
//
//  h1_Peak1_C_MeV->SetLineColor(2);
//  h1_Peak1_C_MeV->SetFillColor(2);
//  h1_Peak1_C_MeV->SetLineWidth(2);
//  h1_Peak1_C_MeV->DrawCopy("SAME");
//
//  TLegend *Legend_021 = new TLegend(.65,.65,.85,.85);
//  Legend_021->SetFillColor(0);
//  Legend_021->SetBorderSize(0);
//  Legend_021->AddEntry(h1_Peak1_MeV  ,"Ge BOT","L");
//  Legend_021->AddEntry(h1_Peak1_C_MeV,"Ge BOT & TOP" ,"L");
//  Legend_021->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 022 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 022 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_022 = new TCanvas("C_022","C_022: ", 1000,600);
//  C_022->Divide(1,1);
//  C_022->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak2_MeV->GetXaxis()->SetTitle("E [MeV]");
//  h1_Peak2_MeV->GetYaxis()->SetTitle("Counts / 10 keV");
//  h1_Peak2_MeV->GetXaxis()->SetRangeUser(0.,     20.);
//  h1_Peak2_MeV->GetYaxis()->SetRangeUser(1.,1000000.);
//  h1_Peak2_MeV->SetLineColor(4);
//  h1_Peak2_MeV->SetLineWidth(2);
//  h1_Peak2_MeV->Draw();
//
//  Int_t  N_COV2_T = h1_Peak2_MeV->GetEntries();
//  cout << " ****** " << N_COV1_T << " " << endl;
//  Int_t N_COV2_1 = N_COV2_T;
//  for (int j =  0; j < 3; j++) {
//    N_COV2_1 = N_COV2_1 - h1_Peak2_MeV->GetBinContent(j + 1);
//    cout << " " << j << " " << h1_Peak2_MeV->GetBinContent(j + 1) << " " << N_COV2_1 << " " << endl;
//  }
//
//  h1_Peak2_C_MeV->SetLineColor(2);
//  h1_Peak2_C_MeV->SetFillColor(2);
//  h1_Peak2_C_MeV->SetLineWidth(2);
//  h1_Peak2_C_MeV->DrawCopy("SAME");
//
//  TLegend *Legend_022 = new TLegend(.65,.65,.85,.85);
//  Legend_022->SetFillColor(0);
//  Legend_022->SetBorderSize(0);
//  Legend_022->AddEntry(h1_Peak2_MeV  ,"Ge TOP","L");
//  Legend_022->AddEntry(h1_Peak2_C_MeV,"Ge TOP & BOT" ,"L");
//  Legend_022->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 023 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 023 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_023 = new TCanvas("C_023","C_023: ", 1000,600);
//  C_023->Divide(1,1);
//  C_023->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_DeltaTCOV->GetXaxis()->SetTitle("#DeltaT [10 #mus]");
//  h1_DeltaTCOV->GetYaxis()->SetTitle("Counts / 10 #mus");
//  h1_DeltaTCOV->GetXaxis()->SetRangeUser(-10.,  10.);
//  h1_DeltaTCOV->GetYaxis()->SetRangeUser(  0.1, 1000.);
//  h1_DeltaTCOV->SetLineColor(4);
//  h1_DeltaTCOV->SetLineWidth(2);
//  h1_DeltaTCOV->Draw();
//
//  TLegend *Legend_023 = new TLegend(.56,.70,.88,.88);
//  Legend_023->SetFillColor(0);
//  Legend_023->SetBorderSize(0);
//  Legend_023->AddEntry(h1_DeltaTCOV,"T_{TOP}-T_{BOT}","L");
//  Legend_023->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 024 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 024 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_024 = new TCanvas("C_024","C_024: ", 1000,600);
//  C_024->Divide(1,1);
//  C_024->SetLogy();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  h1_Peak1_C_keV->GetXaxis()->SetTitle("E [keV]");
//  h1_Peak1_C_keV->GetYaxis()->SetTitle("Counts / 2 keV");
//  h1_Peak1_C_keV->GetXaxis()->SetRangeUser(0.,4000.);
//  h1_Peak1_C_keV->GetYaxis()->SetRangeUser(1.,1000.);
//  h1_Peak1_C_keV->SetLineColor(4);
//  h1_Peak1_C_keV->SetLineWidth(2);
//  h1_Peak1_C_keV->Draw();
//
//  h1_Peak2_C_keV->SetLineColor(2);
//  //h1_Peak2_C_keV->SetFillColor(2);
//  h1_Peak2_C_keV->SetLineWidth(2);
//  h1_Peak2_C_keV->DrawCopy("SAME");
//
//  TLegend *Legend_024 = new TLegend(.65,.65,.85,.85);
//  Legend_024->SetFillColor(0);
//  Legend_024->SetBorderSize(0);
//  Legend_024->AddEntry(h1_Peak1_C_keV,"Ge BOT","L");
//  Legend_024->AddEntry(h1_Peak2_C_keV,"Ge TOP","L");
//  Legend_024->Draw();
//
//  TLatex *Latex_024_1 = new TLatex( 511.,100.,"511 keV");
//  Latex_024_1->SetTextSize(.030);
//  Latex_024_1->SetTextColor(1);
//  Latex_024_1->SetTextAlign(22);
//  Latex_024_1->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
/////////////////////////////////////////////////////////////////////////////////
//////                                                                         //
////// Canvas 025 :                                                            //
//////                                                                         //
/////////////////////////////////////////////////////////////////////////////////
////
////  cout << " ++++++++++++++++++++++++ " << endl;
////  cout << " ++++++ Canvas 025 ++++++ " << endl;
////  cout << " ++++++++++++++++++++++++ " << endl;
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////  TCanvas *C_025 = new TCanvas("C_025","C_025: ", 1000,600);
////  C_025->Divide(1,1);
////  //C_025->SetLogy();
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////  gStyle->SetOptFit(000000);
////  gStyle->SetOptStat("000000001");
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////  h1_Peak1_C_100keV->GetXaxis()->SetTitle("E [keV]");
////  h1_Peak1_C_100keV->GetYaxis()->SetTitle("Counts / keV");
////  h1_Peak1_C_100keV->GetXaxis()->SetRangeUser(0.,250.);
////  h1_Peak1_C_100keV->GetYaxis()->SetRangeUser(0.,120.);
////  h1_Peak1_C_100keV->SetLineColor(4);
////  h1_Peak1_C_100keV->SetLineWidth(2);
////  h1_Peak1_C_100keV->Draw();
////
////  h1_Peak2_C_100keV->SetLineColor(2);
////  //h1_Peak2_C_keV->SetFillColor(2);
////  h1_Peak2_C_100keV->SetLineWidth(2);
////  h1_Peak2_C_100keV->DrawCopy("SAME");
////
////  TLegend *Legend_025 = new TLegend(.65,.65,.85,.85);
////  Legend_025->SetFillColor(0);
////  Legend_025->SetBorderSize(0);
////  Legend_025->AddEntry(h1_Peak1_C_100keV,"Ge BOT","L");
////  Legend_025->AddEntry(h1_Peak2_C_100keV,"Ge TOP","L");
////  Legend_025->Draw();
////  /*
////  TLatex *Latex_025_1 = new TLatex( 511.,100.,"511 keV");
////  Latex_025_1->SetTextSize(.030);
////  Latex_025_1->SetTextColor(1);
////  Latex_025_1->SetTextAlign(22);
////  Latex_025_1->Draw();
////  */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 031 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 031 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_031 = new TCanvas("C_031","C_031: ", 1000,600);
//  C_031->Divide(2,2);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TH1F *h1_Raw1_Mod, *h1_Raw2_Mod, *h1_Raw3_Mod, *h1_Raw4_Mod;
//
//  h1_Raw2_Mod = (TH1F*)h1_Raw2->Clone();
//
//  C_031->cd(1);
//  h1_Raw2_Mod->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw2_Mod->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw2_Mod->GetXaxis()->SetRangeUser(   0.,  500.);
//  h1_Raw2_Mod->GetYaxis()->SetRangeUser(-100.,  100.);
//  h1_Raw2_Mod->SetLineColor(4);
//  h1_Raw2_Mod->SetLineWidth(2);
//  h1_Raw2_Mod->Draw();
//
//  TLegend *Legend_031_1 = new TLegend(.15,.70,.35,.85);
//  Legend_031_1->SetFillColor(0);
//  Legend_031_1->SetBorderSize(0);
//  Legend_031_1->AddEntry(h1_Raw2_Mod,"Ge BOT - 1","L");
//  Legend_031_1->Draw();
//
//  h1_Raw4_Mod = (TH1F*)h1_Raw4->Clone();
//
//  C_031->cd(2);
//  h1_Raw4_Mod->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw4_Mod->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw4_Mod->GetXaxis()->SetRangeUser(   0.,  500.);
//  h1_Raw4_Mod->GetYaxis()->SetRangeUser(-450., -250.);
//  h1_Raw4_Mod->SetLineColor(4);
//  h1_Raw4_Mod->SetLineWidth(2);
//  h1_Raw4_Mod->Draw();
//
//  TLegend *Legend_031_2 = new TLegend(.15,.70,.35,.85);
//  Legend_031_2->SetFillColor(0);
//  Legend_031_2->SetBorderSize(0);
//  Legend_031_2->AddEntry(h1_Raw4_Mod,"Ge TOP - 1","L");
//  Legend_031_2->Draw();
//
//  h1_Raw1_Mod = (TH1F*)h1_Raw1->Clone();
//
//  C_031->cd(3);
//  h1_Raw1_Mod->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw1_Mod->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw1_Mod->GetXaxis()->SetRangeUser(   0.,  500.);
//  h1_Raw1_Mod->GetYaxis()->SetRangeUser(-100.,  100.);
//  h1_Raw1_Mod->SetLineColor(4);
//  h1_Raw1_Mod->SetLineWidth(2);
//  h1_Raw1_Mod->Draw();
//
//  TLegend *Legend_031_3 = new TLegend(.15,.15,.35,.30);
//  Legend_031_3->SetFillColor(0);
//  Legend_031_3->SetBorderSize(0);
//  Legend_031_3->AddEntry(h1_Raw1_Mod,"Ge BOT - 2","L");
//  Legend_031_3->Draw();
//
//  h1_Raw3_Mod = (TH1F*)h1_Raw3->Clone();
//
//  C_031->cd(4);
//  h1_Raw3_Mod->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw3_Mod->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw3_Mod->GetXaxis()->SetRangeUser(   0.,  500.);
//  h1_Raw3_Mod->GetYaxis()->SetRangeUser(-450., -250.);
//  h1_Raw3_Mod->SetLineColor(4);
//  h1_Raw3_Mod->SetLineWidth(2);
//  h1_Raw3_Mod->Draw();
//
//  TLegend *Legend_031_4 = new TLegend(.15,.15,.35,.30);
//  Legend_031_4->SetFillColor(0);
//  Legend_031_4->SetBorderSize(0);
//  Legend_031_4->AddEntry(h1_Raw3_Mod,"Ge TOP - 2","L");
//  Legend_031_4->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 041 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 041 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_041 = new TCanvas("C_041","C_041: ", 1000,600);
//  C_041->Divide(2,2);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TH1F *h1_Raw1_1 = h1_Raw1->Clone("h1_raw1_1");
//  TH1F *h1_Raw2_1 = h1_Raw2->Clone("h1_raw2_1");
//  TH1F *h1_Raw3_1 = h1_Raw3->Clone("h1_raw3_1");
//  TH1F *h1_Raw4_1 = h1_Raw4->Clone("h1_raw4_1");
//
//  C_041->cd(1);
//
//  h1_Raw2_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw2_1->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw2_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  h1_Raw2_1->GetYaxis()->SetTitleOffset(1.5);
//  h1_Raw2_1->SetLineColor(4);
//  h1_Raw2_1->SetLineWidth(2);
//  h1_Raw2_1->Draw();
//
//  TLegend *Legend_041_1 = new TLegend(.15,.70,.35,.85);
//  Legend_041_1->SetFillColor(0);
//  Legend_041_1->SetBorderSize(0);
//  Legend_041_1->AddEntry(h1_Raw2_1,"Ge BOT - 1","L");
//  Legend_041_1->Draw();
//
//  C_041->cd(2);
//  h1_Raw4_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw4_1->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw4_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  h1_Raw4_1->GetYaxis()->SetTitleOffset(1.5);
//  h1_Raw4_1->SetLineColor(4);
//  h1_Raw4_1->SetLineWidth(2);
//  h1_Raw4_1->Draw();
//
//  TLegend *Legend_041_2 = new TLegend(.15,.70,.35,.85);
//  Legend_041_2->SetFillColor(0);
//  Legend_041_2->SetBorderSize(0);
//  Legend_041_2->AddEntry(h1_Raw4_1,"Ge TOP - 1","L");
//  Legend_041_2->Draw();
//
//  C_041->cd(3);
//  h1_Raw1_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw1_1->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw1_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  h1_Raw1_1->GetYaxis()->SetTitleOffset(1.5);
//  h1_Raw1_1->SetLineColor(4);
//  h1_Raw1_1->SetLineWidth(2);
//  h1_Raw1_1->Draw();
//
//  TLegend *Legend_041_3 = new TLegend(.15,.15,.35,.30);
//  Legend_041_3->SetFillColor(0);
//  Legend_041_3->SetBorderSize(0);
//  Legend_041_3->AddEntry(h1_Raw1_1,"Ge BOT - 2","L");
//  Legend_041_3->Draw();
//
//  C_041->cd(4);
//  h1_Raw3_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_Raw3_1->GetYaxis()->SetTitle("Amplitude [ADC]");
//  h1_Raw3_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  h1_Raw3_1->GetYaxis()->SetTitleOffset(1.5);
//  h1_Raw3_1->SetLineColor(4);
//  h1_Raw3_1->SetLineWidth(2);
//  h1_Raw3_1->Draw();
//
//  TLegend *Legend_041_4 = new TLegend(.15,.15,.35,.30);
//  Legend_041_4->SetFillColor(0);
//  Legend_041_4->SetBorderSize(0);
//  Legend_041_4->AddEntry(h1_Raw3_1,"Ge TOP - 2","L");
//  Legend_041_4->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
////                                                                         //
//// Canvas 042 :                                                            //
////                                                                         //
///////////////////////////////////////////////////////////////////////////////
//
//  cout << " ++++++++++++++++++++++++ " << endl;
//  cout << " ++++++ Canvas 042 ++++++ " << endl;
//  cout << " ++++++++++++++++++++++++ " << endl;
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TCanvas *C_042 = new TCanvas("C_042","C_042: ", 1000,600);
//  C_042->Divide(2,2);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  gStyle->SetOptFit(000000);
//  gStyle->SetOptStat("000000001");
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  TH1F *h1_COV1_1 = h1_COV1->Clone("h1_COV1_1");
//  TH1F *h1_COV2_1 = h1_COV2->Clone("h1_COV2_1");
//  TH1F *h1_DIF1_1 = h1_DIF1->Clone("h1_DIF1_1");
//  TH1F *h1_DIF2_1 = h1_DIF2->Clone("h1_DIF2_1");
//
//  C_042->cd(1);
//  h1_COV1_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_COV1_1->GetYaxis()->SetTitle("Amplitude [V]");
//  h1_COV1_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  //h1_COV1_1->GetYaxis()->SetRangeUser(.0, 1000000.);
//  h1_COV1_1->SetLineColor(4);
//  h1_COV1_1->SetLineWidth(2);
//  h1_COV1_1->Draw();
//
//  TLegend *Legend_042_1 = new TLegend(.15,.70,.35,.85);
//  Legend_042_1->SetFillColor(0);
//  Legend_042_1->SetBorderSize(0);
//  Legend_042_1->AddEntry(h1_COV1_1,"Ge BOT","L");
//  Legend_042_1->Draw();
//
//  C_042->cd(2);
//  h1_COV2_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_COV2_1->GetYaxis()->SetTitle("Amplitude [V]");
//  h1_COV2_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  //h1_COV1_1->GetYaxis()->SetRangeUser(.0, 1000000.);
//  h1_COV2_1->SetLineColor(4);
//  h1_COV2_1->SetLineWidth(2);
//  h1_COV2_1->Draw();
//
//  TLegend *Legend_042_2 = new TLegend(.15,.70,.35,.85);
//  Legend_042_2->SetFillColor(0);
//  Legend_042_2->SetBorderSize(0);
//  Legend_042_2->AddEntry(h1_COV2_1,"Ge TOP","L");
//  Legend_042_2->Draw();
//
//
//  C_042->cd(3);
//  for (Int_t j = 0; j < 100000; j++) {h1_DIF1_1->SetBinContent(j + 1,h1_DIF1_1->GetBinContent(j + 1)/(1./1000.));}
//  h1_DIF1_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_DIF1_1->GetYaxis()->SetTitle("#DeltaAmplitude [mV]");
//  h1_DIF1_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  h1_DIF1_1->GetYaxis()->SetRangeUser(-10., 10.);
//  h1_DIF1_1->SetLineColor(4);
//  h1_DIF1_1->SetLineWidth(2);
//  h1_DIF1_1->Draw();
//
//  TLegend *Legend_042_3 = new TLegend(.15,.70,.35,.85);
//  Legend_042_3->SetFillColor(0);
//  Legend_042_3->SetBorderSize(0);
//  Legend_042_3->AddEntry(h1_DIF1_1,"Ge BOT","L");
//  Legend_042_3->Draw();
//
//  TLine *Line1_1 = new TLine(28000.,-5.*Cal1,29000.,-5.*Cal1);
//  Line1_1->SetLineStyle(2);
//  Line1_1->SetLineColor(2);
//  Line1_1->SetLineWidth(4);
//  Line1_1->Draw();
//
//  TLine *Line1_2 = new TLine(28000., 5.*Cal1,29000., 5.*Cal1);
//  Line1_2->SetLineStyle(2);
//  Line1_2->SetLineColor(2);
//  Line1_2->SetLineWidth(4);
//  Line1_2->Draw();
//
//  C_042->cd(4);
//  for (Int_t j = 0; j < 100000; j++) {h1_DIF2_1->SetBinContent(j + 1,h1_DIF2_1->GetBinContent(j + 1)/(1./1000.));}
//  h1_DIF2_1->GetXaxis()->SetTitle("T [10#mus]");
//  h1_DIF2_1->GetYaxis()->SetTitle("#DeltaAmplitude [mV]");
//  h1_DIF2_1->GetXaxis()->SetRangeUser(28000.,29000.);
//  h1_DIF2_1->GetYaxis()->SetRangeUser(-10., 10.);
//  h1_DIF2_1->SetLineColor(4);
//  h1_DIF2_1->SetLineWidth(2);
//  h1_DIF2_1->Draw();
//
//  TLegend *Legend_042_4 = new TLegend(.15,.70,.35,.85);
//  Legend_042_4->SetFillColor(0);
//  Legend_042_4->SetBorderSize(0);
//  Legend_042_4->AddEntry(h1_DIF2_1,"Ge TOP","L");
//  Legend_042_4->Draw();
//
//  TLine *Line2_1 = new TLine(28000.,-5.*Cal2,29000.,-5.*Cal2);
//  Line2_1->SetLineStyle(2);
//  Line2_1->SetLineColor(2);
//  Line2_1->SetLineWidth(4);
//  Line2_1->Draw();
//
//  TLine *Line2_2 = new TLine(28000., 5.*Cal2,29000., 5.*Cal2);
//  Line2_2->SetLineStyle(2);
//  Line2_2->SetLineColor(2);
//  Line2_2->SetLineWidth(4);
//  Line2_2->Draw();
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    ///
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    /////////////////////////////////////////////////////////////////////////////
//    //                                                                         //
//    // Canvas 043 :                                                            //
//    //                                                                         //
//    /////////////////////////////////////////////////////////////////////////////
//
//      cout << " ++++++++++++++++++++++++ " << endl;
//      cout << " ++++++ Canvas 043 ++++++ " << endl;
//      cout << " ++++++++++++++++++++++++ " << endl;
//
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      TCanvas *C_043 = new TCanvas("C_043","C_043: ", 1000,600);
//
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      gStyle->SetOptFit(000000);
//      gStyle->SetOptStat("000000001");
//
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      h1_NTD_1_MeV->GetXaxis()->SetTitle("Reconstructed energy [MeVee]");
//    h1_NTD_1_MeV->GetYaxis()->SetTitle("Counts / (keV)");
//    h1_NTD_1_MeV->GetXaxis()->SetRangeUser(0.,6.);
//      //h1_COV1_1->GetYaxis()->SetRangeUser(.0, 1000000.);
//    h1_NTD_1_MeV->SetLineColor(3);
//    h1_NTD_1_MeV->SetLineWidth(2);
//
//    x=0.;
//    while (x<0.050){
//        h1_NTD_1_MeV->SetBinContent(h1_NTD_1_MeV->FindBin(x), 0);
//        x+= h1_NTD_1_MeV->GetBinWidth(1);
//    }
//
//    h1_NTD_1_MeV->Draw();
//
//      TLegend *Legend_043_1 = new TLegend(.15,.70,.35,.85);
//    Legend_043_1->SetFillColor(0);
//    Legend_043_1->SetBorderSize(0);
//    Legend_043_1->AddEntry(h1_NTD_1_MeV,"LWO","L");
//    Legend_043_1->Draw();

//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  cout << " " << endl;
//  cout << " " << endl;
//  cout << " *************************** " << endl;
//  cout << " ***** Results Summary ***** " << endl;
//  cout << " *************************** " << endl;
//  cout << " " << endl;
//  cout << " ***********************************************************************************" << endl;
//  cout << " " << endl;
//  cout << " File name : " << file_name << " " << endl;
//  cout << " " << endl;
//  cout << " RunTime (s) = " ;
//  cout << fixed << setprecision(1) << setw(7) << setfill(' ') << RunTime << " (" ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << RunTime/3600. << " h)" << endl;
//  cout << " Ge BOT Threshold (keV) = " << int(Threshold_COV1) << " " << endl;
//  cout << " Ge TOP Threshold (keV) = " << int(Threshold_COV2) << " " << endl;
//  cout << " Ge BOT events above threshold = " << N_COV1 << "  ==> Rate (cps) = " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << N_COV1/RunTime << " +/- " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << sqrt(N_COV1)/RunTime << " " << endl;
//  cout << " Ge BOT events above 30 keV = " << N_COV1_1 << "  ==> Rate (cps) = " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << N_COV1_1/RunTime << " +/- " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << sqrt(N_COV1_1)/RunTime << " " << endl;
//  cout << " Ge BOT 60 keV events = " << N_COV1_60keV << "  ==> Rate (cps) = " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << N_COV1_60keV/RunTime << " +/- " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << sqrt(N_COV1_60keV)/RunTime << " " << endl;
//  cout << " Ge TOP events above threshold = " << N_COV2 << "  ==> Rate (cps) = " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << N_COV2/RunTime << " +/- " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << sqrt(N_COV2)/RunTime << " " << endl;
//  cout << " Ge TOP events above 30 keV = " << N_COV2 << "  ==> Rate (cps) = " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << N_COV2_1/RunTime << " +/- " ;
//  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << sqrt(N_COV2_1)/RunTime << " " << endl;
//  cout << " Ge BOT Calibration Slope (mV/keV) = " ;
//  cout << fixed << setprecision(4) << setw(5) << setfill(' ') << 1000.*fit_01->GetParameter(1) << " +/- " ;
//  cout << fixed << setprecision(4) << setw(5) << setfill(' ') << 1000.*fit_01->GetParError(1) << " " << endl;
//  cout << " Ge BOT Calibration Offset (mV) = " ;
//  cout << fixed << setprecision(3) << setw(5) << setfill(' ') << 1000.*fit_01->GetParameter(0) << " +/- " ;
//  cout << fixed << setprecision(3) << setw(5) << setfill(' ') << 1000.*fit_01->GetParError(0) << " " << endl;
//  cout << " Ge TOP Calibration Slope (mV/keV) = " ;
//  cout << fixed << setprecision(4) << setw(5) << setfill(' ') << 1000.*fit_02->GetParameter(1) << " +/- " ;
//  cout << fixed << setprecision(4) << setw(5) << setfill(' ') << 1000.*fit_02->GetParError(1) << " " << endl;
//  cout << " Ge TOP Calibration Offset (mV) = " ;
//  cout << fixed << setprecision(3) << setw(5) << setfill(' ') << 1000.*fit_02->GetParameter(0) << " +/- " ;
//  cout << fixed << setprecision(3) << setw(5) << setfill(' ') << 1000.*fit_02->GetParError(0) << " " << endl;
//  cout << " Ge BOT Gain = " << int(Gain_COV1) << " " << endl;
//  cout << " Ge TOP Gain = " << int(Gain_COV2) << " " << endl;
//  cout << " Ge BOT Sensitivity (nV/keV) = " ;
//  cout << fixed << setprecision(1) << setw(4) << setfill(' ') << 1000.*fit_01->GetParameter(1)/(Gain_COV1/1000000.) << " " << endl;
//  cout << " Ge TOP Sensitivity (nV/keV) = " ;
//  cout << fixed << setprecision(1) << setw(4) << setfill(' ') << 1000.*fit_02->GetParameter(1)/(Gain_COV2/1000000.) << " " << endl;
//  cout << " Ge BOT BSL RMS (keV) = " ;
//  cout << fixed << setprecision(1) << setw(4) << setfill(' ') << COV1_BSL_RMS << " " << endl;
//  cout << " Ge TOP BSL RMS (keV) = " ;
//  cout << fixed << setprecision(1) << setw(4) << setfill(' ') << COV2_BSL_RMS << " " << endl;
//  cout << " " << endl;
//
//  cout << " ***********************************************************************************" << endl;
//
//
// }
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
