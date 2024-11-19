double ENERGY_UNIT = 1e3; //energy unit with respect to keV (GeV =Rate_above_5MeV 1e6 Mev = 1e3, keV = 1, eV = 1e-3,...)

double MASS_Ge = 400*1e-3;
double MASS_LWO = 51.7*1e-3;

int NEW_NBINS = 500;

double veto_lim_min = -2.;
double veto_lim_max = 4.;

void help_form_dataMyCOV1(){
    cout<<"Fit Normalisation of simulations to the COV Prototype data:"<<endl;
    cout<< ""<< endl;
    cout<<"================ CONSTANT VALUES ============="<<endl;
    cout<<"Mass of Ge detector = "<< MASS_Ge<< " kg"<<endl;
    cout<<"Mass of LWO detector = "<< MASS_LWO<< " kg"<<endl;
    cout<<""<<endl;
    cout<<"================ END OF CONSTANT VALUES ============="<<endl;

    cout<<""<<endl;
    cout<<"void Hits_to_cps_data(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, int RunTime=1, int New_Nbins=NEW_NBINS)"<<endl;
    cout<<"\t Normalize in cps the data files"<<endl;
}

// -----------------------------------------
void Hits_to_cps_data(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, int RunTime=1, int New_Nbins=NEW_NBINS){
    
    int Nbins = histo->GetNbinsX();
    double RebinFact = Nbins/New_Nbins;
    //cout<<Nbins<<"/"<<NEW_NBINS<<"="<<RebinFact<<endl;
    histo->Rebin(RebinFact);
    
    Nbins = histo->GetNbinsX();
    double deltaE = histo->GetBinWidth(1)*energy_unit;
    
    double Hits_i, Error_i, cps_i, Error_cps_i;
    for (int i=1; i<=Nbins; i++){
        Hits_i = histo->GetBinContent(i);
        Error_i = histo->GetBinError(i);
        
        cps_i = Hits_i/(RunTime*detect_mass*deltaE);
        Error_cps_i = Error_i/(RunTime*detect_mass*deltaE);
    
        histo->SetBinContent(i, cps_i);
        histo->SetBinError(i, Error_cps_i);
    }
    
    histo->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    //gStyle->SetOptStat(0);
    histo->Draw();
}

// -----------------------------------------
void Write_one_histogram(TH1D* h, TString Crystal = "COV1", Double_t Emax= 4., Double_t Mass = MASS_Ge, int RunTime=10000, int NewNbins=200.){
    
    TString Out_name= Form("%s_Norm_Spectrum_0-%1.fMeV", Crystal.Data(), Emax);
    
    int NumBin = h->FindBin(Emax)-h->FindBin(0);
    TH1D *new_h = new TH1D(Out_name, Out_name, NumBin, 0, Emax);
    
    int Hits_i, Error_i;
    for (int i=1; i<=NumBin; i++){
        Hits_i = h->GetBinContent(i);
        Error_i = h->GetBinError(i);
        
        new_h->SetBinContent(i, Hits_i);
        new_h->SetBinError(i, Error_i);
    }
    
    new_h->GetXaxis()->SetTitle("Energy [MeV]");
    new_h->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h, 1e3, Mass, RunTime, NewNbins);
    
    new_h->Write(Out_name);
    new_h->Delete();
}

// -----------------------------------------
void Veto_histogram(TFile* Filein, Double_t RunTime){
    TH1D *h1_NTD_1_MeV = (TH1D*) Filein->Get("h1_NTD_1_MeV");
    
    TH1D* h1_DeltaT_Closest = (TH1D*) Filein->Get("h1_DeltaT_Closest");
    
    TF1 *fit01_Left = new TF1("fit01_Left","expo",-25.0, -2.0);
    TF1 *fit01_Right = new TF1("fit01_Right","expo",  4.0, 25.0);

    h1_DeltaT_Closest->Fit("fit01_Left","QNR");
    h1_DeltaT_Closest->Fit("fit01_Right","QNR");

    Double_t P0 = 0.5*(fit01_Left->GetParameter(0) + fit01_Right->GetParameter(0)); //Mean normalization
    Double_t P1 = 0.5*(fit01_Left->GetParameter(1) - fit01_Right->GetParameter(1));  //0.5*(pente_left - pente_right) ??

    Double_t Binwidth = h1_DeltaT_Closest->GetBinWidth(1);
    Double_t Sum_L = exp(P0)*(1. - exp( P1*veto_lim_min))/P1/Binwidth;
    Double_t Sum_R = exp(P0)*(1. - exp(-P1*veto_lim_max))/P1/Binwidth;
    
    Double_t N_acc_in_Veto = Sum_L+Sum_R;
    
    TH1D *h1_NTD_11_MeV = (TH1D*) Filein->Get("h1_NTD_11_MeV"); //Vetoed events
    TH1D *h1_NTD_12_MeV = (TH1D*) Filein->Get("h1_NTD_12_MeV"); //Accidentals
    
    Double_t N_veto = h1_NTD_11_MeV->GetEntries();
    Double_t N_tot_accidentals = h1_NTD_12_MeV->GetEntries();

    TH1D* h1_NTD_Veto_MeV = (TH1D*) h1_NTD_1_MeV->Clone();
    TH1D* h1_NTD_Veto_withAcc_MeV = (TH1D*) h1_NTD_1_MeV->Clone();

    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_11_MeV, 1., -1.);
    h1_NTD_Veto_withAcc_MeV->Add(h1_NTD_Veto_withAcc_MeV,h1_NTD_11_MeV, 1., -1.);
    
    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_12_MeV, 1., N_acc_in_Veto/N_tot_accidentals);
    
    h1_NTD_Veto_MeV->GetXaxis()->SetTitle("Energy [MeV]");
    h1_NTD_Veto_MeV->GetYaxis()->SetTitle("Counts [a.u.]");
    h1_NTD_Veto_MeV->Write("LWO-Vetoed_Raw_Spectrum_0-20MeV");
    
    Write_one_histogram(h1_NTD_Veto_MeV, "LWO-Vetoed", 20., MASS_LWO, RunTime, 100.);
    Write_one_histogram(h1_NTD_Veto_MeV, "LWO-Vetoed", 10., MASS_LWO, RunTime, 500.);
    
    h1_NTD_Veto_withAcc_MeV->GetXaxis()->SetTitle("Energy [MeV]");
    h1_NTD_Veto_withAcc_MeV->GetYaxis()->SetTitle("Counts [a.u.]");
    h1_NTD_Veto_withAcc_MeV->Write("LWO-Vetoed_withAcc_Raw_Spectrum_0-20MeV");
    
    Write_one_histogram(h1_NTD_Veto_withAcc_MeV, "LWO-Vetoed_withAcc", 20., MASS_LWO, RunTime, 100.);
    Write_one_histogram(h1_NTD_Veto_withAcc_MeV, "LWO-Vetoed_withAcc", 10., MASS_LWO, RunTime, 500.);
}

// -----------------------------------------
void form_data_fromMyCOV1(TString FileinName, TString FileoutName){
    
    TCanvas *c1 = new TCanvas();
    TFile *Filein = new TFile(FileinName, "READ");
    TFile *Fileout = new TFile(FileoutName, "RECREATE");
    
    TH1D *histo_RunTime = (TH1D*) Filein->Get("h1_RunTime");
    Double_t RunTime=histo_RunTime->GetBinContent(1);
    
    int Hits_i, Error_i, NumBin;
    
    //----------------------  COV1 --------------------------------
    Fileout->mkdir("BotGe");
    Fileout->cd("BotGe");
    TH1D *h = (TH1D*) Filein->Get("h1_Peak1_4MeV");
    
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle("Counts [a.u.]");
    h->Write("BotGe_Raw_Spectrum_0-20MeV");
    
    Write_one_histogram(h, "BotGe", 4., MASS_Ge, RunTime, 200);
    Write_one_histogram(h, "BotGe", 1., MASS_Ge, RunTime, 500);
    
    h = (TH1D*) Filein->Get("h1_Peak1_MeV");
    Write_one_histogram(h, "BotGe", 20., MASS_Ge, RunTime, 100);
    Write_one_histogram(h, "BotGe", 10., MASS_Ge, RunTime, 500);
    //----------------------  end of COV1 --------------------------------
    
    //----------------------  COV2 --------------------------------
    Fileout->cd();
    Fileout->mkdir("TopGe");
    Fileout->cd("TopGe");
    h = (TH1D*) Filein->Get("h1_Peak2_4MeV");
    
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle("Counts [a.u.]");
    h->Write("TopGe_Raw_Spectrum_0-20MeV");
    
    Write_one_histogram(h, "TopGe", 4., MASS_Ge, RunTime, 200);
    Write_one_histogram(h, "TopGe", 1., MASS_Ge, RunTime, 500);
    
    h = (TH1D*) Filein->Get("h1_Peak2_MeV");
    Write_one_histogram(h, "TopGe", 20., MASS_Ge, RunTime, 100);
    Write_one_histogram(h, "TopGe", 10., MASS_Ge, RunTime, 500);
    //----------------------  end of COV2 --------------------------------
    
    //----------------------  NTD --------------------------------
    Fileout->cd();
    Fileout->mkdir("LWO");
    Fileout->cd("LWO");
    h = (TH1D*) Filein->Get("h1_NTD_1_MeV");
    
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle("Counts [a.u.]");
    h->Write("LWO_Raw_Spectrum_0-20MeV");
    
    Write_one_histogram(h, "LWO", 20., MASS_LWO, RunTime, 100.);
    Write_one_histogram(h, "LWO", 10., MASS_LWO, RunTime, 500.);
//    Write_one_histogram(h, "LWO", 1, MASS_LWO, RunTime, 100.);
    
    //Veto
    Veto_histogram(Filein, RunTime);
    //----------------------  end of NTD --------------------------------
    
    Fileout->Close();
    Filein->Close();
    c1->Close();
}

// -----------------------------------------
//---- CovMat
void Save_Cov_Mat_tot(TFile* Filein, TString Crystal, Double_t Emax, Double_t Mass = MASS_Ge, int RunTime=10000){
    
    TString Out_name= Form("%s_CovMat_0-%1.fMeV", Crystal.Data(), Emax);
    
    TH2D* SystCovMat = (TH2D*) Filein->Get(Form("0-%1.fMeV/%s_SystCovMat",  Emax, Crystal.Data()));
    TH2D* StatCovMat = (TH2D*) Filein->Get(Form("0-%1.fMeV/%s_StatCovMat",  Emax, Crystal.Data()));

    TH2D* SumCovMat= (TH2D*) SystCovMat->Clone();
    
    SumCovMat->SetTitle(Out_name);
    SumCovMat->SetName(Out_name);
    SumCovMat->Add(SumCovMat, StatCovMat, 1., 1.);
    
    Double_t deltaE = SumCovMat->GetXaxis()->GetBinWidth(1)*1e3;
    int Nbins = SumCovMat->GetNbinsX();
    Double_t Var_ij, var_cps_ij;
    for (int i=1; i<=Nbins; i++){
        for (int j=1; j<=Nbins; j++){
            Var_ij = SumCovMat->GetBinContent(i, j);
        
            var_cps_ij = Var_ij/((RunTime*Mass*deltaE)*(RunTime*Mass*deltaE));
            if((Out_name=="BotGe_CovMat_0-4MeV")&&(i==100)&&(j==100)){
                cout<<sqrt(var_cps_ij)<<endl;
            }
        
            SumCovMat->SetBinContent(i,j, var_cps_ij);
        }
    }
    
//    SumCovMat->GetZaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    SumCovMat->Write(Out_name);
}
//----

// -----------------------------------------
void form_data_Syst_Stat(TString FileinName, TString FileinName_veto, TString FileoutName){
    
    TCanvas *c1 = new TCanvas();
    TFile *Filein = new TFile(FileinName, "READ");
    TFile *Fileinveto = new TFile(FileinName_veto, "READ");
    TFile *Fileout = new TFile(FileoutName, "RECREATE");
    
    TH1D *histo_RunTime = (TH1D*) Fileinveto->Get("h1_RunTime");
    Double_t RunTime=histo_RunTime->GetBinContent(1);
    
    int Hits_i, Error_i, NumBin;
    
    //----------------------  COV1 --------------------------------
    Fileout->mkdir("BotGe");
    Fileout->cd("BotGe");
    TH1D *h = (TH1D*) Filein->Get("0-4MeV/BotGe_Ref_calib_statErrors");
    cout<<h->GetBinError(100)<<endl;
    Write_one_histogram(h, "BotGe", 4., MASS_Ge, RunTime, 200);
    Save_Cov_Mat_tot(Filein, "BotGe", 4., MASS_Ge, RunTime);
    
    h = (TH1D*) Filein->Get("0-1MeV/BotGe_Ref_calib_statErrors");
    Write_one_histogram(h, "BotGe", 1., MASS_Ge, RunTime, 500);
    Save_Cov_Mat_tot(Filein, "BotGe", 1., MASS_Ge, RunTime);
    
    h = (TH1D*) Filein->Get("0-20MeV/BotGe_Ref_calib_statErrors");
    Write_one_histogram(h, "BotGe", 20., MASS_Ge, RunTime, 100);
    Save_Cov_Mat_tot(Filein, "BotGe", 20., MASS_Ge, RunTime);
    
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle("Counts [a.u.]");
    h->Write("BotGe_Raw_Spectrum_0-20MeV");
    
    h = (TH1D*) Filein->Get("0-10MeV/BotGe_Ref_calib_statErrors");
    Write_one_histogram(h, "BotGe", 10., MASS_Ge, RunTime, 500);
    Save_Cov_Mat_tot(Filein, "BotGe", 10., MASS_Ge, RunTime);
    //----------------------  end of COV1 --------------------------------
    
    //----------------------  COV2 --------------------------------
    Fileout->cd();
    Fileout->mkdir("TopGe");
    Fileout->cd("TopGe");
    h = (TH1D*) Filein->Get("0-4MeV/TopGe_Ref_calib_statErrors");
    Write_one_histogram(h, "TopGe", 4., MASS_Ge, RunTime, 200);
    Save_Cov_Mat_tot(Filein, "TopGe", 4., MASS_Ge, RunTime);
    
    h = (TH1D*) Filein->Get("0-1MeV/TopGe_Ref_calib_statErrors");
    Write_one_histogram(h, "TopGe", 1., MASS_Ge, RunTime, 500);
    Save_Cov_Mat_tot(Filein, "TopGe", 1., MASS_Ge, RunTime);
    
    h = (TH1D*) Filein->Get("0-20MeV/TopGe_Ref_calib_statErrors");
    Write_one_histogram(h, "TopGe", 20., MASS_Ge, RunTime, 100);
    Save_Cov_Mat_tot(Filein, "TopGe", 20., MASS_Ge, RunTime);
    
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle("Counts [a.u.]");
    h->Write("TopGe_Raw_Spectrum_0-20MeV");
    
    h = (TH1D*) Filein->Get("0-10MeV/TopGe_Ref_calib_statErrors");
    Write_one_histogram(h, "TopGe", 10., MASS_Ge, RunTime, 500);
    Save_Cov_Mat_tot(Filein, "TopGe", 10., MASS_Ge, RunTime);
    //----------------------  end of COV2 --------------------------------
    
    //----------------------  NTD --------------------------------
    Fileout->cd();
    Fileout->mkdir("LWO");
    Fileout->cd("LWO");
    h = (TH1D*) Filein->Get("0-20MeV/LWO_Ref_calib_statErrors");
    
    h->GetXaxis()->SetTitle("Energy [MeV]");
    h->GetYaxis()->SetTitle("Counts [a.u.]");
    h->Write("LWO_Raw_Spectrum_0-20MeV");
    Write_one_histogram(h, "LWO", 20., MASS_LWO, RunTime, 100.);
    Save_Cov_Mat_tot(Filein, "LWO", 20., MASS_LWO, RunTime);
    
    h = (TH1D*) Filein->Get("0-10MeV/LWO_Ref_calib_statErrors");
    Write_one_histogram(h, "LWO", 10., MASS_LWO, RunTime, 500.);
    Save_Cov_Mat_tot(Filein, "LWO", 10., MASS_LWO, RunTime);
    
//    Write_one_histogram(h, "LWO", 1, MASS_LWO, RunTime, 100.);
    
    //Veto
    Veto_histogram(Fileinveto, RunTime);
    //----------------------  end of NTD --------------------------------
    
    Fileout->Close();
    Fileinveto->Close();
    Filein->Close();
    c1->Close();
}
