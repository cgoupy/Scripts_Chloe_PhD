double CALIB_COV1_A = 3.635*1E-9; //V/keV^2
double CALIB_COV1_B = 0.0009421;//V/keV
double CALIB_COV1_C = 0.003067;//V

double CALIB_COV2 = 0.0006501; //V/keV
double CALIB_COV2_OR = -0.0006176;//V

double CALIB_NTD = 0.000252; //V/keV
double CALIB_NTD_OR = 0.005886;//V

int DRAWNUM = 1000000000;//

TH1D* RecalibHisto(TH1D* histo, TString Crystal, int Nbins, double Emax){
    double Volt_i, E_pol;
    int Ener_index, Hits_i, NumBin;
    
    NumBin=histo->GetNbinsX();
    TString title = histo->GetTitle();
    TH1D *new_histo = new TH1D("recal_"+title, "recal_"+title, Nbins, 0, Emax);
    
    for (int i=0; i<DRAWNUM; i++){
        Volt_i = histo->GetRandom();
        if (Crystal=="COV1") {E_pol = 1E-3*max((-CALIB_COV1_B+sqrt(CALIB_COV1_B*CALIB_COV1_B-4*CALIB_COV1_A*(CALIB_COV1_C-Volt_i)))/(2*CALIB_COV1_A), (-CALIB_COV1_B-sqrt(CALIB_COV1_B*CALIB_COV1_B-4*CALIB_COV1_A*(CALIB_COV1_C-Volt_i)))/(2*CALIB_COV1_A));}
        else if (Crystal=="COV2") {E_pol = 1E-3*(Volt_i-CALIB_COV2_OR)/CALIB_COV2;}
        else if (Crystal=="NTD") {E_pol = 1E-3*(Volt_i-CALIB_NTD_OR)/CALIB_NTD;}
        new_histo->Fill(E_pol);
    }
    
    //new_histo->Rebin(10);
    new_histo->GetXaxis()->SetTitle("Energy [MeV]");
    new_histo->GetYaxis()->SetTitle("Counts [a.u.]");

    double V_Max;
    if (Crystal=="COV1") { V_Max= CALIB_COV1_A*Emax*Emax*1E6+CALIB_COV1_B*Emax*1E3+CALIB_COV1_C;}
    if (Crystal=="COV2") { V_Max= CALIB_COV2_OR+CALIB_COV2*Emax*1E3;}
    if (Crystal=="NTD") { V_Max= CALIB_NTD_OR+CALIB_NTD*Emax*1E3;}
    
    new_histo->Scale(histo->Integral(0, histo->FindBin(V_Max))/new_histo->Integral(0,new_histo->FindBin(Emax)));
    
    for (int i=1; i<new_histo->GetNbinsX(); i++){
        new_histo->SetBinError(i, sqrt(new_histo->GetBinContent(i)));
    }
    
    return new_histo;
}

void Recalibfile(TString FileinName, TString FileoutName){
    TFile *Filein = new TFile(FileinName, "READ");
    TFile *Fileout = new TFile(FileoutName, "RECREATE");
    
    TH1D *histo_RunTime = (TH1D*) Filein->Get("h1_RunTime");
    histo_RunTime->Write("h1_RunTime");
    
    double Volt_i, E_pol_COV1;
    int Ener_index, Hits_i, NumBin;
    
    //COV1
    TH1D *h_COV1 = (TH1D*) Filein->Get("h1_Peak1");
    h_COV1->GetXaxis()->SetTitle("Amplitude [V]");
    h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    h_COV1->Write("h1_Peak1");
    
    NumBin=h_COV1->GetNbinsX();
    
    TH1D *new_h_COV1_4 = RecalibHisto(h_COV1, "COV1", 2000, 4);
    new_h_COV1_4->Write("h1_Peak1_4MeV");
    new_h_COV1_4->Delete();
    
    TH1D* new_h_COV1_20 = RecalibHisto(h_COV1, "COV1", 2000, 20);
    new_h_COV1_20->Write("h1_Peak1_MeV");
    new_h_COV1_20->Delete();
    
    //COV2
    TH1D *h_COV2 = (TH1D*) Filein->Get("h1_Peak2_MeV");
    h_COV2->Write("h1_Peak2_MeV");
    
    h_COV2 = (TH1D*) Filein->Get("h1_Peak2_4MeV");
    h_COV2->Write("h1_Peak2_4MeV");

    //NTD
    TH1D *h_NTDAmp = (TH1D*) Filein->Get("h1_NTD_Amplitude");
    TH1D* new_h_NTD = RecalibHisto(h_NTDAmp, "NTD", 2000, 20);
    new_h_NTD->Write("h1_NTD_1_MeV_myRecal");
    new_h_NTD->Delete();
    TH1D *h_NTD = (TH1D*) Filein->Get("h1_NTD_1_MeV");
    h_NTD->Write("h1_NTD_1_MeV");
    h_NTD = (TH1D*) Filein->Get("h1_NTD_2_MeV");
    h_NTD->Write("h1_NTD_2_MeV");
    h_NTD = (TH1D*) Filein->Get("h1_NTD_7_MeV");
    h_NTD->Write("h1_NTD_7_MeV");
    h_NTD = (TH1D*) Filein->Get("h1_NTD_8_MeV");
    h_NTD->Write("h1_NTD_8_MeV");
    h_NTD = (TH1D*) Filein->Get("h1_NTD_10_MeV");
    h_NTD->Write("h1_NTD_10_MeV");
    

    Fileout->Close();
    Filein->Close();
}
