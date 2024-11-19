//Code to fit the simulation of COV proto to simulations
#include "Recalib.C"
#include "Shift.C"
#include "form_data_fromMyCOV1.C"

double FLUX = 0.0134; // events/cm2/s neutrons : 0.0134 ev/cm2/s - muons : 0.019 ev/cm2/s - AmbientGammas : 1.25e6 ev/cm2/d // 3 ev/cm^2/s ?
double PLANE_SURF = 140*140; //in cm^2 //1 for Am (-> Rate in ev/s)
double RATE = FLUX*PLANE_SURF; //events/s

//double ENERGY_UNIT = 1e3; //energy unit with respect to keV (GeV =Rate_above_5MeV 1e6 Mev = 1e3, keV = 1, eV = 1e-3,...)
//
//defined in form_data_fromMyCOV1.C
//double MASS_Ge = 400*1e-3;
//double MASS_LWO = 51.7*1e-3;
//
//int NEW_NBINS = 500;

double THRESH_TOP = 10; //keV
double THRESH_BOT = 4; //keV
double THRESH_LWO = 25; //keV

double NUM_EVENTS = 1e8;
double NUM_E_ANALYZED=2; //1 if no Evis (no resolution applied in the analysis)

double HQF = 1.07; //Heat quenching factor =/= 1 only for neutrons

void help(){
    cout<<"Fit Normalisation of simulations to the COV Prototype data:"<<endl;
    cout<< ""<< endl;
    cout<<"================ CONSTANT VALUES ============="<<endl;
    cout<<"Flux = " <<FLUX<<" events/cm2/s"<<endl;
    cout<<"Surface of plane generating the events = " <<PLANE_SURF<< " cm^2"<<endl;
    cout<<"Rate = "<<RATE<< " events/d"<<endl;
    cout<<"Default energy unit = "<< ENERGY_UNIT << " keV" << endl;
    cout<<"Mass of Ge detector = "<< MASS_Ge<< " kg"<<endl;
    cout<<"Mass of LWO detector = "<< MASS_LWO<< " kg"<<endl;
    cout<<""<<endl;
    cout<<"Thresholds"<<endl;
    cout<<"Threshold top Ge = "<<THRESH_TOP<<" keV"<<endl;
    cout<<"Threshold bottom Ge = "<<THRESH_BOT<<" keV"<<endl;
    cout<<"Threshold LWO = "<<THRESH_LWO<<" keV"<<endl;
    cout<<""<<endl;
    if (NUM_E_ANALYZED==2){cout<<"Edep and Evis analyzed"<<endl;} else if (NUM_E_ANALYZED==1){cout<<"Edep analyzed"<<endl;}
    cout<<""<<endl;
    cout<<"================ END OF CONSTANT VALUES ============="<<endl;
    
    cout<<""<<endl;
    cout<<"double SetNumEvents(TString FileName)"<<endl;
    cout<<"\t return the number of events generated in the file, from the primaries histograms"<<endl;
    
    cout<<""<<endl;
    cout<<"void Hits_to_cps(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE)"<<endl;
    cout<<"\t normalize the histogram histo in count per seconds"<<endl;
    
    cout<<""<<endl;
    cout<<"void Hits_to_dru(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE, int New_Nbins=NEW_NBINS)"<<endl;
    cout<<"\t normalize the histogram histo in count per day per kg"<<endl;
    
    cout<<""<<endl;
    cout<<"TH1D* GoodFormat(TH1D* histo, TString title, int New_Nbins, double Bin_MinVal, double Bin_MaxVal, bool ShiftNeutrons = false)"<<endl;
    cout<<"\t Put the histogram in the New_Nbins number of bins between min and max values, if ShiftNeutrons == true, the histogram is shifted by HQF"<<endl;
    
    cout<<""<<endl;
    cout<<"ROOT::Fit::FitResult Fit_Data_Simu(...)"<<endl;
    cout<<"\t Fitting function of the simulation on the data --  used in Fit_OnRange"<<endl;
    
    cout<<""<<endl;
    cout<<"void plot_with_res(vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron)"<<endl;
    cout<<"\t Plot the fit results with the residuals"<<endl;
    
    cout<<""<<endl;
    cout<<"void plot_with_res_pc(vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron)"<<endl;
    cout<<"\t Plot the fit results with the residuals in %"<<endl;

    cout<<""<<endl;
    cout<<"ROOT::Fit::FitResult Fit_OnRange(vector<Double_t> NumEvents, std::vector<TH1D*> histo_data, std::vector<TH1D*> histo_mu, std::vector<TH1D*> histo_gamma, std::vector<TH1D*> histo_n, std::vector<TH1D*> histo_Am, double Emin, double Emax, Double_t ini_flux_mu, Double_t ini_flux_gamma, Double_t ini_flux_n, Double_t ini_flux_Am, bool mu_fixed = false, bool gamma_fixed=false, bool n_fixed=true)"<<endl;
    cout<<"\t Apply the fit function Fit_Data_Simu on the specified Range and parameters histograms -- used in FitNorm"<<endl;
    
    cout<<""<<endl;
    cout<<"void FitNorm(TString NameFileout, TString fileData_name, TString filein_mu_name, TString filein_gamma_name, TString filein_n_name, TString filein_Am241_name)"<<endl;
    cout<<"\t Main function to fit the simulations on the data and save in a file named NameFileout"<<endl;
    
    cout<<""<<endl;
    cout<<"vector<double> Int_Rate(TH1D* histo, double Estart=0, double Eend=0, double MassTarget = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE)"<<endl;
    cout<<"\t print and return the integrated daily rate in histo between Estart and Eend with its error"<<endl;
    
    cout<<""<<endl;
    cout<<"void Norm_in_cps(TString Filemname, TString NameFileout)"<<endl;
    cout<<"\t Create a root file with the histograms normalized in cps and print the rate in different energy ranges"<<endl;
    
    cout<<""<<endl;
    cout<<"void recalib(TString FileinName, TString FileoutName)"<<endl;
    cout<<"\t uses the function from macro Recalib.C"<<endl;
    
    cout<<""<<endl;
    cout<<"void form_data(TString FileinName, TString FileoutName)"<<endl;
    cout<<"\t Create a root file with the data histograms (from MyCOV Eddy's macro) normalized in cps in the same ranges as the one used in the simulation."<<endl;
    
    cout<<""<<endl;
    cout<<"void Norm_in_cps_simu_veto(TString FileinName, TString FileoutName)"<<endl;
    cout<<"\t --- calls save_CPS_histos()"<<endl;
    cout<<"\t Create a root file with the histograms from OrsayOV analyse normalized in cps and print the rate in different energy ranges"<<endl;
    
    cout<<""<<endl;
    cout<<"void save_CPS_histos(TFile *filein, TFile *fileout, TString Crystal, TString hName, double Crystal_Mass, double thresh_crys, double thresh_crys_fixed, double NumEv, double Rate)"<<endl;
    cout<<"\t Save in fileout the CPS normalized histograms root named -hName- (full histogram path without the energy) of the crystal -Crystal- (Top_Ge, Mid_LWO or Bot_Ge) with its correspond thresholds "<<endl;
    
    cout<<""<<endl;
    cout<<"void Save_Fitted_Simu_range(TFile* Fileout, int Emax, int Emax_simu, TString Trigger, vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am)"<<endl;
    cout<<"\t Save in fileout the CPS normalized histograms corresponding to the defined 'Trigger' ('Triggered', 'Coincidences') from the simulation results"<<endl;
    
    cout<<""<<endl;
    cout<<"void Save_Fitted(TString FileoutName, vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am)"<<endl;
    cout<<"\t Save in FileoutName the whole Simulation/Data fitted. -- call Save_Fitted_Simu_range and called in FitNorm"<<endl;
    
}
void plot_with_res(vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am);
void plot_with_res_pc(vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am, bool Syst=false);
void Save_Fitted(TString FileoutName, vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am, bool = false);

// ---------------------------------------
void Hits_to_cps(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE, int New_Nbins=NEW_NBINS){
    
    int Nbins = histo->GetNbinsX();
    double RebinFact = Nbins/New_Nbins;
    
    //cout<<Nbins<<"/"<<NEW_NBINS<<"="<<RebinFact<<endl;
    
    histo->Rebin(RebinFact);

    Nbins = histo->GetNbinsX();
    double deltaE = histo->GetBinWidth(1)*energy_unit;
    
    double Hits_i, Error_i, dru_i, Error_dru_i;
    for (int i=1; i<=Nbins; i++){
        Hits_i = histo->GetBinContent(i);
        Error_i = histo->GetBinError(i);
        
        dru_i = (Hits_i*rate)/(num_events*detect_mass*deltaE);
        Error_dru_i = (Error_i*rate)/(num_events*detect_mass*deltaE);
        
        histo->SetBinContent(i, dru_i);
        histo->SetBinError(i, Error_dru_i);
    }
    
    histo->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    gStyle->SetOptStat(0);
    //histo->Draw("SAME");
}

// ---------------------------------------
void Hits_to_dru(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE, int New_Nbins=NEW_NBINS){
    
    int Nbins = histo->GetNbinsX();
    double RebinFact = Nbins/New_Nbins;
    
    //cout<<Nbins<<"/"<<NEW_NBINS<<"="<<RebinFact<<endl;
    
    histo->Rebin(RebinFact);

    Nbins = histo->GetNbinsX();
    double deltaE = histo->GetBinWidth(1)*energy_unit;
    
    double Hits_i, Error_i, dru_i, Error_dru_i;
    for (int i=1; i<=Nbins; i++){
        Hits_i = histo->GetBinContent(i);
        Error_i = histo->GetBinError(i);
        
        dru_i = (Hits_i*rate*3600*24)/(num_events*detect_mass*deltaE);
        Error_dru_i = (Error_i*rate*3600*24)/(num_events*detect_mass*deltaE);
        
        histo->SetBinContent(i, dru_i);
        histo->SetBinError(i, Error_dru_i);
    }
    
    histo->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    gStyle->SetOptStat(0);
    //histo->Draw("SAME");
}

// ---------------------------------------
double SetNumEvents(TString FileName = "Analyzer_EP_Nu1g_AtmNSphere_Al2O3_10Mevents.root", bool Am = FALSE){
    
    TFile *file = new TFile(FileName);
    TH1D *histoPrim = (TH1D*) file->Get("PrimariesAlgoDir/Energy/hEinj_log");
    double NumEvents=histoPrim->GetEntries();
    
    if (Am){
        //TH1D *histo = (TH1D*) file->Get("EdepAlgoDir/Vol_Bot_Ge/Edep/nEdep_Vol_Bot_Ge_0-100keV");
        TH1D *histo = (TH1D*) file->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-100keV");
        NumEvents=histo->Integral(histo->FindBin(59.3), histo->FindBin(59.8));
    }
    
    return NumEvents;
}

// ---------------------------------------
TH1D* GoodFormat(TH1D* histo, TString title, int New_Nbins, double Bin_MinVal, double Bin_MaxVal, bool ShiftNeutrons = false, double hqf = HQF){
    
    if (ShiftNeutrons){histo = Shift(histo, hqf);}
    
    int MinBin = histo->FindBin(Bin_MinVal);
    int MaxBin = histo->FindBin(Bin_MaxVal);
    
    TH1D* new_histo = new TH1D(title+"_cut", title+"_cut", MaxBin-MinBin, Bin_MinVal, Bin_MaxVal);
    
    Double_t Hits_i, Error_i;
    for (int i=0; i<=MaxBin-MinBin; i++){
        Hits_i = histo->GetBinContent(i+MinBin-1);
        Error_i = histo->GetBinError(i+MinBin-1);
        
        new_histo->SetBinContent(i, Hits_i);
        new_histo->SetBinError(i, Error_i);
    }
    
    
    int RebinFact = new_histo->GetNbinsX()/New_Nbins;
    new_histo->Rebin(RebinFact);
    new_histo->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    new_histo->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    
//    cout<<title<<" : "<<new_histo->GetNbinsX()<<"/"<<New_Nbins<<endl;
    return new_histo;
}

// ---------------------------------------
void Cut(TH1D* histo, TH2D* Mat, TString title, double Bin_MinVal, double Bin_MaxVal){
    
    int MinBin = histo->FindBin(Bin_MinVal);
    int MaxBin = histo->FindBin(Bin_MaxVal);
    
    TH1D* new_histo = new TH1D(title+"_cut", title+"_cut", MaxBin-MinBin, Bin_MinVal, Bin_MaxVal);
    TH2D* new_Mat = new TH2D("Cov_Mat_"+title+"_cut", "Cov_Mat_"+title+"_cut", MaxBin-MinBin, Bin_MinVal, Bin_MaxVal, MaxBin-MinBin, Bin_MinVal, Bin_MaxVal);
    
    cout<<"MaxBin-MinBin = "<<MaxBin-MinBin<<endl;
    
    Double_t Hits_i, Error_i, Var_ij;
    for (int i=0; i<=MaxBin-MinBin; i++){
        Hits_i = histo->GetBinContent(i+MinBin-1);
        Error_i = histo->GetBinError(i+MinBin-1);
        
        for (int j=0; j<=MaxBin-MinBin; j++){
            Var_ij = Mat->GetBinContent(i+MinBin-1, j+MinBin-1);
            new_Mat->SetBinContent(i, j, Var_ij);
        }
        
        new_histo->SetBinContent(i, Hits_i);
        new_histo->SetBinError(i, Error_i);
    }
    
    new_histo->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    new_histo->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    
    new_Mat->GetYaxis()->SetTitle(Mat->GetYaxis()->GetTitle());
    new_Mat->GetXaxis()->SetTitle(Mat->GetXaxis()->GetTitle());
    new_Mat->GetZaxis()->SetTitle(Mat->GetZaxis()->GetTitle());
    
//    cout<<title<<" : "<<new_histo->GetNbinsX()<<"/"<<New_Nbins<<endl;
    
    *histo= *new_histo;
    *Mat = *new_Mat;
}


// --------------- Fit function
ROOT::Fit::FitResult Fit_Data_Simu(vector<Double_t> NumEvents, TVectorT<Double_t> &E_bins, std::vector<TVectorT<Double_t>> Data_vect, std::vector<TVectorT<Double_t>> Simu_vect_mu, std::vector<TMatrixT<Double_t>> V_Simu_mu, std::vector<TVectorT<Double_t>> Simu_vect_gamma,std::vector<TMatrixT<Double_t>> V_Simu_gamma, std::vector<TVectorT<Double_t>> Simu_vect_n, std::vector<TMatrixT<Double_t>> V_Simu_n, std::vector<TVectorT<Double_t>> Simu_vect_Am241,std::vector<TMatrixT<Double_t>> V_Simu_Am, std::vector<TMatrixT<Double_t>> V_Data, TString &Data_name, Double_t ini_flux_mu, Double_t ini_flux_gamma, Double_t ini_flux_n, Double_t ini_flux_Am, Double_t deltaE, bool mu_fixed = false, bool gamma_fixed=false, bool n_fixed=true, bool Am_fixed=true) {

    //E_bins are the bins centers of the data and simulation histogram
    //Data_vect is the vector of the data values
    //Simu_vect is the vector of the simulation values
    //V_Data and V_Simu are the variance matrix of the data and the simulation
    
    int N_E = E_bins.GetNrows();
    //TVectorT<Double_t> Delta_Data_SimuConv(N_NPE); //difference between Data and alpha*Simu
    
    //Projection on a base where diagonal elements V_Data =/= 0
    int N_zerosCOV1 = 0;
    int N_zerosCOV2 = 0;
    int N_zerosNTD = 0;
    
    for (int i =0; i<N_E; i++){
        N_zerosCOV1+=(Data_vect[0](i)==0); //count the number of zero in Data_vect
        N_zerosCOV2+=(Data_vect[1](i)==0); //count the number of zero in Data_vect
        N_zerosNTD+=(Data_vect[2](i)==0); //count the number of zero in Data_vect
    }
    
    int N_ProjCOV1 = N_E-N_zerosCOV1;
    int N_ProjCOV2 = N_E-N_zerosCOV2;
    int N_ProjNTD = N_E-N_zerosNTD;
    
    TMatrixT<Double_t> ProjCOV1(N_ProjCOV1, N_E);
    TMatrixT<Double_t> ProjCOV2(N_ProjCOV2, N_E);
    TMatrixT<Double_t> ProjNTD(N_ProjNTD, N_E);
    
    int i_ProjCOV1 = 0; int i_ProjCOV2 = 0; int i_ProjNTD = 0;
    for (int j=0; j<N_E; j++){
        if (Data_vect[0](j)!=0){
            ProjCOV1(i_ProjCOV1, j)=1;
            i_ProjCOV1+=1;
        }
        if (Data_vect[1](j)!=0){
            ProjCOV2(i_ProjCOV2, j)=1;
            i_ProjCOV2+=1;
        }
        if (Data_vect[2](j)!=0){
            ProjNTD(i_ProjNTD, j)=1;
            i_ProjNTD+=1;
        }
    }
    
    //---- projection data and simulation spectra
    TVectorT<Double_t> Proj_Data_vectCOV1 = ProjCOV1*Data_vect[0];
    TVectorT<Double_t> Proj_Data_vectCOV2 = ProjCOV2*Data_vect[1];
    TVectorT<Double_t> Proj_Data_vectNTD = ProjNTD*Data_vect[2];
    
    TVectorT<Double_t> Proj_Simu_vect_muCOV1 = ProjCOV1*Simu_vect_mu[0];
    TVectorT<Double_t> Proj_Simu_vect_muCOV2 = ProjCOV2*Simu_vect_mu[1];
    TVectorT<Double_t> Proj_Simu_vect_muNTD = ProjNTD*Simu_vect_mu[2];
    
    TVectorT<Double_t> Proj_Simu_vect_gammaCOV1 = ProjCOV1*Simu_vect_gamma[0];
    TVectorT<Double_t> Proj_Simu_vect_gammaCOV2 = ProjCOV2*Simu_vect_gamma[1];
    TVectorT<Double_t> Proj_Simu_vect_gammaNTD = ProjNTD*Simu_vect_gamma[2];
    
    TVectorT<Double_t> Proj_Simu_vect_nCOV1 = ProjCOV1*Simu_vect_n[0];
    TVectorT<Double_t> Proj_Simu_vect_nCOV2 = ProjCOV2*Simu_vect_n[1];
    TVectorT<Double_t> Proj_Simu_vect_nNTD = ProjNTD*Simu_vect_n[2];
    
    TVectorT<Double_t> Proj_Simu_vect_Am241COV1 = ProjCOV1*Simu_vect_Am241[0];
    TVectorT<Double_t> Proj_Simu_vect_Am241COV2 = ProjCOV2*Simu_vect_Am241[1];
    TVectorT<Double_t> Proj_Simu_vect_Am241NTD = ProjNTD*Simu_vect_Am241[2];
    
    
    //---- projection V_data
    TMatrixT<Double_t> ProjCOV1_T(TMatrixT<Double_t>::kTransposed, ProjCOV1);
    TMatrixT<Double_t> ProjCOV2_T(TMatrixT<Double_t>::kTransposed, ProjCOV2);
    TMatrixT<Double_t> ProjNTD_T(TMatrixT<Double_t>::kTransposed, ProjNTD);
    
    TMatrixT<Double_t> Proj_V_DataCOV1 = ProjCOV1*V_Data[0]*ProjCOV1_T;
    TMatrixT<Double_t> Proj_V_DataCOV2 = ProjCOV2*V_Data[1]*ProjCOV2_T;
    TMatrixT<Double_t> Proj_V_DataNTD = ProjNTD*V_Data[2]*ProjNTD_T;
    
    TMatrixT<Double_t> Proj_V_Simu_muCOV1 = ProjCOV1*V_Simu_mu[0]*ProjCOV1_T;
    TMatrixT<Double_t> Proj_V_Simu_muCOV2 = ProjCOV2*V_Simu_mu[1]*ProjCOV2_T;
    TMatrixT<Double_t> Proj_V_Simu_muNTD = ProjNTD*V_Simu_mu[2]*ProjNTD_T;
    
    TMatrixT<Double_t> Proj_V_Simu_gammaCOV1 = ProjCOV1*V_Simu_gamma[0]*ProjCOV1_T;
    TMatrixT<Double_t> Proj_V_Simu_gammaCOV2 = ProjCOV2*V_Simu_gamma[1]*ProjCOV2_T;
    TMatrixT<Double_t> Proj_V_Simu_gammaNTD = ProjNTD*V_Simu_gamma[2]*ProjNTD_T;
    
    TMatrixT<Double_t> Proj_V_Simu_nCOV1 = ProjCOV1*V_Simu_n[0]*ProjCOV1_T;
    TMatrixT<Double_t> Proj_V_Simu_nCOV2 = ProjCOV2*V_Simu_n[1]*ProjCOV2_T;
    TMatrixT<Double_t> Proj_V_Simu_nNTD = ProjNTD*V_Simu_n[2]*ProjNTD_T;
    
    TMatrixT<Double_t> Proj_V_Simu_AmCOV1 = ProjCOV1*V_Simu_Am[0]*ProjCOV1_T;
    TMatrixT<Double_t> Proj_V_Simu_AmCOV2 = ProjCOV2*V_Simu_Am[1]*ProjCOV2_T;
    TMatrixT<Double_t> Proj_V_Simu_AmNTD = ProjNTD*V_Simu_Am[2]*ProjNTD_T;
    
    
    //---- projection E_bins
    TVectorT<Double_t> Proj_E_binsCOV1 = ProjCOV1*E_bins; TVectorT<Double_t> Proj_E_binsCOV2 = ProjCOV2*E_bins; TVectorT<Double_t> Proj_E_binsNTD = ProjNTD*E_bins;
    
    TMatrixT<Double_t> A_COV1(N_ProjCOV1, N_E); TMatrixT<Double_t> A_COV2(N_ProjCOV2, N_E); TMatrixT<Double_t> A_NTD(N_ProjNTD, N_E);
    
    TMatrixT<Double_t> A_COV1_T(N_E, N_ProjCOV1); TMatrixT<Double_t> A_COV2_T(N_E, N_ProjCOV2); TMatrixT<Double_t> A_NTD_T(N_E, N_ProjNTD);
    
    TMatrixT<Double_t> V_Delta_COV1_i(N_ProjCOV1, N_ProjCOV1); TMatrixT<Double_t> V_Delta_COV2_i(N_ProjCOV2, N_ProjCOV2); TMatrixT<Double_t> V_Delta_NTD_i(N_ProjNTD, N_ProjNTD);
    
    TMatrixT<Double_t> V_Delta_COV1(N_ProjCOV1, N_ProjCOV1); TMatrixT<Double_t> V_Delta_COV2(N_ProjCOV2, N_ProjCOV2); TMatrixT<Double_t> V_Delta_NTD(N_ProjNTD, N_ProjNTD);
    
    TVectorT<Double_t> Delta_Data_Simu_COV1(N_ProjCOV1); TVectorT<Double_t> Delta_Data_Simu_COV2(N_ProjCOV2); TVectorT<Double_t> Delta_Data_Simu_NTD(N_ProjNTD);
    
    //cf. https://root.cern.ch/doc/master/fitCircle_8C.html
    auto chi2Function_Fit = [&](const Double_t *par) {
        //minimisation function computing the sum of squares of chi2
        Double_t flux_mu=par[0];
        Double_t flux_gamma = par[1];
        Double_t flux_n = par[2];
        Double_t flux_Am = par[3];
        
        //Norm from flux
        Double_t Norm_mu = (flux_mu*PLANE_SURF)/(NumEvents[0]*MASS_Ge*deltaE);
        Double_t Norm_gamma = (flux_gamma*PLANE_SURF)/(NumEvents[1]*MASS_Ge*deltaE);
        Double_t Norm_n = (flux_n*PLANE_SURF)/(NumEvents[2]*MASS_Ge*deltaE);
        Double_t Norm_Am = (flux_Am)/(NumEvents[3]*MASS_Ge*deltaE);
        
        /// ---------- COV1
        //Diff Data - normalized simulations
        Delta_Data_Simu_COV1 = Proj_Data_vectCOV1 - Norm_mu*Proj_Simu_vect_muCOV1 - Norm_gamma * Proj_Simu_vect_gammaCOV1 - Norm_n * Proj_Simu_vect_nCOV1 - Norm_Am * Proj_Simu_vect_Am241COV1;
        
        //Variance matrix Data - normalized simulations
        V_Delta_COV1 = Proj_V_DataCOV1 + Norm_mu*Norm_mu*Proj_V_Simu_muCOV1 + Norm_gamma*Norm_gamma*Proj_V_Simu_gammaCOV1 + Norm_n*Norm_n*Proj_V_Simu_nCOV1 + Norm_Am*Norm_Am*Proj_V_Simu_AmCOV1;
        TMatrixT<Double_t> V_Delta_COV1_i(TMatrixT<Double_t>::kInverted, V_Delta_COV1);
        
        /// ------------ COV2
        //Diff Data - normalized simulations
        Delta_Data_Simu_COV2 = Proj_Data_vectCOV2 - Norm_mu*Proj_Simu_vect_muCOV2 - Norm_gamma * Proj_Simu_vect_gammaCOV2 - Norm_n * Proj_Simu_vect_nCOV2 - Norm_Am * Proj_Simu_vect_Am241COV2;
        
        //Variance matrix Data - normalized simulations
        V_Delta_COV2 = Proj_V_DataCOV2 + Norm_mu*Norm_mu*Proj_V_Simu_muCOV2 + Norm_gamma*Norm_gamma*Proj_V_Simu_gammaCOV2 + Norm_n*Norm_n*Proj_V_Simu_nCOV2 + Norm_Am*Norm_Am*Proj_V_Simu_AmCOV2;
        TMatrixT<Double_t> V_Delta_COV2_i(TMatrixT<Double_t>::kInverted, V_Delta_COV2);
        
        /// ------------ LWO
        Norm_mu = (flux_mu*PLANE_SURF)/(NumEvents[0]*MASS_LWO*deltaE);
        Norm_gamma = (flux_gamma*PLANE_SURF)/(NumEvents[1]*MASS_LWO*deltaE);
        Norm_n = (flux_n*PLANE_SURF)/(NumEvents[2]*MASS_LWO*deltaE);
        Norm_Am = (flux_Am)/(NumEvents[3]*MASS_LWO*deltaE);
        
        //Diff Data - normalized simulations
        Delta_Data_Simu_NTD = Proj_Data_vectNTD - Norm_mu*Proj_Simu_vect_muNTD - Norm_gamma * Proj_Simu_vect_gammaNTD - Norm_n * Proj_Simu_vect_nNTD - Norm_Am * Proj_Simu_vect_Am241NTD;
        
        //Variance matrix Data - normalized simulations
        V_Delta_NTD = Proj_V_DataNTD + Norm_mu*Norm_mu*Proj_V_Simu_muNTD + Norm_gamma*Norm_gamma*Proj_V_Simu_gammaNTD + Norm_n*Norm_n*Proj_V_Simu_nNTD + Norm_Am*Norm_Am*Proj_V_Simu_AmNTD;
        TMatrixT<Double_t> V_Delta_NTD_i(TMatrixT<Double_t>::kInverted, V_Delta_NTD);
        
        /// ------------ Total Chi2
        double chi2 = Delta_Data_Simu_COV1 * (V_Delta_COV1_i * Delta_Data_Simu_COV1) + Delta_Data_Simu_COV2 * (V_Delta_COV2_i * Delta_Data_Simu_COV2) + Delta_Data_Simu_NTD * (V_Delta_NTD_i * Delta_Data_Simu_NTD);
        if (!Am_fixed){ chi2 = Delta_Data_Simu_COV1 * (V_Delta_COV1_i * Delta_Data_Simu_COV1);} // if Am not fixed = Am241 fit -> fit only COV1
        return chi2;
    };
    
    // wrap chi2 function in a function object for the fit
    // 4 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn_Fit(chi2Function_Fit, 4);
    ROOT::Fit::Fitter fitter_Fit;
    
    // Initial fit parameters
    Double_t pStart_Fit[4] = {ini_flux_mu, ini_flux_gamma, ini_flux_n, ini_flux_Am};
    fitter_Fit.SetFCN(fcn_Fit, pStart_Fit);

    fitter_Fit.Config().ParSettings(0).SetName("Flux_muon");
    fitter_Fit.Config().ParSettings(1).SetName("Flux_gamma");
    fitter_Fit.Config().ParSettings(2).SetName("Flux_neutrons");
    fitter_Fit.Config().ParSettings(3).SetName("Flux_Am241");
    
    fitter_Fit.Config().ParSettings(0).SetLimits(0.5*ini_flux_mu, 3*ini_flux_mu);
    fitter_Fit.Config().ParSettings(1).SetLimits(0.5*ini_flux_gamma, 3*ini_flux_gamma);
    fitter_Fit.Config().ParSettings(2).SetLimits(0.01, 10*ini_flux_n);
    //fitter_Fit.Config().ParSettings(3).SetLimits(0.5*ini_flux_Am, 3*ini_flux_Am);
    
    if (mu_fixed) fitter_Fit.Config().ParSettings(0).Fix();
    if (gamma_fixed) fitter_Fit.Config().ParSettings(1).Fix();
    if (n_fixed) fitter_Fit.Config().ParSettings(2).Fix();
    if (Am_fixed) fitter_Fit.Config().ParSettings(3).Fix();
    
    // do the fit
    bool ok_Fit = fitter_Fit.FitFCN();
    
    // test fit
    if(!ok_Fit) {
        Error("fit", "%s", Form("%s Fit failed", Data_name.Data()));
    }
    
    const ROOT::Fit::FitResult &result_Fit = fitter_Fit.Result();
    //    result_Fit.Print(std::cout);
    
    //Calcul individual chi2
    Double_t flux_mu = result_Fit.Parameter(0);
    Double_t flux_gamma = result_Fit.Parameter(1);
    Double_t flux_n = result_Fit.Parameter(2);
    Double_t flux_Am = result_Fit.Parameter(3);
    
    Double_t Norm_mu = (flux_mu*PLANE_SURF)/(NumEvents[0]*MASS_Ge*deltaE);
    Double_t Norm_gamma = (flux_gamma*PLANE_SURF)/(NumEvents[1]*MASS_Ge*deltaE);
    Double_t Norm_n = (flux_n*PLANE_SURF)/(NumEvents[2]*MASS_Ge*deltaE);
    Double_t Norm_Am = (flux_Am)/(NumEvents[3]*MASS_Ge*deltaE);
    
    int N_par = result_Fit.NPar();
    
    /// ----------------------- COV1
    int N_dof = Data_vect[0].GetNrows() - N_par;
    
    Delta_Data_Simu_COV1 = Proj_Data_vectCOV1 - Norm_mu*Proj_Simu_vect_muCOV1 - Norm_gamma * Proj_Simu_vect_gammaCOV1 - Norm_n * Proj_Simu_vect_nCOV1 - Norm_Am * Proj_Simu_vect_Am241COV1;
    V_Delta_COV1 = Proj_V_DataCOV1 + Norm_mu*Norm_mu*Proj_V_Simu_muCOV1 + Norm_gamma*Norm_gamma*Proj_V_Simu_gammaCOV1 + Norm_n*Norm_n*Proj_V_Simu_nCOV1 + Norm_Am*Norm_Am*Proj_V_Simu_AmCOV1;
    TMatrixT<Double_t> V_Delta_COV1_i2(TMatrixT<Double_t>::kInverted, V_Delta_COV1);
    cout<<"chi2 COV1/ndf = "<<Delta_Data_Simu_COV1 * (V_Delta_COV1_i2 * Delta_Data_Simu_COV1)<<" / "<<N_dof<<" = "<<(Delta_Data_Simu_COV1 * (V_Delta_COV1_i2 * Delta_Data_Simu_COV1))/N_dof<<endl;
    
    /// ----------------------- COV2
    N_dof = Data_vect[1].GetNrows() - N_par;
    Delta_Data_Simu_COV2 = Proj_Data_vectCOV2 - Norm_mu*Proj_Simu_vect_muCOV2 - Norm_gamma * Proj_Simu_vect_gammaCOV2 - Norm_n * Proj_Simu_vect_nCOV2 - Norm_Am * Proj_Simu_vect_Am241COV2;
    V_Delta_COV2 = Proj_V_DataCOV2 + Norm_mu*Norm_mu*Proj_V_Simu_muCOV2 + Norm_gamma*Norm_gamma*Proj_V_Simu_gammaCOV2 + Norm_n*Norm_n*Proj_V_Simu_nCOV2 + Norm_Am*Norm_Am*Proj_V_Simu_AmCOV2;
    TMatrixT<Double_t> V_Delta_COV2_i2(TMatrixT<Double_t>::kInverted, V_Delta_COV2);
    cout<<"chi2 COV2/ndf = "<<Delta_Data_Simu_COV2 * (V_Delta_COV2_i2 * Delta_Data_Simu_COV2)<<" / "<<N_dof<<" = "<<(Delta_Data_Simu_COV2 * (V_Delta_COV2_i2 * Delta_Data_Simu_COV2))/N_dof<<endl;
    
    /// ----------------------- LWO
    Norm_mu = (flux_mu*PLANE_SURF)/(NumEvents[0]*MASS_LWO*deltaE);
    Norm_gamma = (flux_gamma*PLANE_SURF)/(NumEvents[1]*MASS_LWO*deltaE);
    Norm_n = (flux_n*PLANE_SURF)/(NumEvents[2]*MASS_LWO*deltaE);
    Norm_Am = (flux_Am)/(NumEvents[3]*MASS_LWO*deltaE);
    
    N_dof = Data_vect[2].GetNrows() - N_par;
    Delta_Data_Simu_NTD = Proj_Data_vectNTD - Norm_mu*Proj_Simu_vect_muNTD - Norm_gamma * Proj_Simu_vect_gammaNTD - Norm_n * Proj_Simu_vect_nNTD - Norm_Am * Proj_Simu_vect_Am241NTD;
    V_Delta_NTD = Proj_V_DataNTD + Norm_mu*Norm_mu*Proj_V_Simu_muNTD + Norm_gamma*Norm_gamma*Proj_V_Simu_gammaNTD + Norm_n*Norm_n*Proj_V_Simu_nNTD + Norm_Am*Norm_Am*Proj_V_Simu_AmNTD;
    TMatrixT<Double_t> V_Delta_NTD_i2(TMatrixT<Double_t>::kInverted, V_Delta_NTD);
    cout<<"chi2 NTD/ndf = "<<Delta_Data_Simu_NTD * (V_Delta_NTD_i2 * Delta_Data_Simu_NTD)<<" / "<<N_dof<<" = "<<(Delta_Data_Simu_NTD * (V_Delta_NTD_i2 * Delta_Data_Simu_NTD))/N_dof<<endl;
    cout<<""<<endl;
    return result_Fit;
}

int Maximum(int value1, int value2){
    if (value1>value2) {return value1;}
    else {return value2;}
}

// ---------------------------
ROOT::Fit::FitResult Fit_OnRange(vector<Double_t> NumEvents, std::vector<TH1D*> histo_data, std::vector<TH2D*> Var_data, std::vector<TH1D*> histo_mu, std::vector<TH1D*> histo_gamma, std::vector<TH1D*> histo_n, std::vector<TH1D*> histo_Am, double Emin, double Emax, Double_t ini_flux_mu, Double_t ini_flux_gamma, Double_t ini_flux_n, Double_t ini_flux_Am, bool mu_fixed = false, bool gamma_fixed=false, bool n_fixed=true, bool Am_fixed=true){
    
    // Prepare data vectors
    TH1D* histo_dataCOV1 = histo_data[0];
    TH1D* histo_dataCOV2 = histo_data[1];
    TH1D* histo_dataNTD = histo_data[2];
    
    // Prepare data vectors
    TH2D* V_dataCOV1 = Var_data[0];
    TH2D* V_dataCOV2 = Var_data[1];
    TH2D* V_dataNTD = Var_data[2];
    
    vector<Int_t> MaxBins={histo_dataCOV1->FindBin(Emax), histo_dataCOV2->FindBin(Emax), histo_dataNTD->FindBin(Emax)};
    vector<Int_t> MinBins={histo_dataCOV1->FindBin(Emin), histo_dataCOV2->FindBin(Emin), histo_dataNTD->FindBin(Emin)};
    
    int New_Nbins_COV1=MaxBins[0]-MinBins[0];
    cout<<"New_Nbins_COV1 = "<<New_Nbins_COV1<<endl;
    
    int New_Nbins_COV2=MaxBins[1]-MinBins[1];
    int New_Nbins_NTD=MaxBins[2]-MinBins[2];
    int New_Nbins=Maximum(New_Nbins_COV1, New_Nbins_COV2);
    New_Nbins=Maximum(New_Nbins, New_Nbins_NTD);
    
    double deltaE = histo_dataCOV1->GetBinWidth(1)*1e3;
    
    TString Data_name =histo_dataCOV1->GetName();
    
    TVectorT<Double_t> FormatVector(New_Nbins);
    TMatrixT<Double_t> FormatMatrix(New_Nbins, New_Nbins);
    
    TVectorT<Double_t> E_bins(New_Nbins) ; std::vector<TVectorT<Double_t>> Simu_vect_mu={FormatVector, FormatVector, FormatVector};  std::vector<TVectorT<Double_t>> Simu_vect_gamma={FormatVector, FormatVector, FormatVector}; std::vector<TVectorT<Double_t>> Simu_vect_n={FormatVector, FormatVector, FormatVector}; std::vector<TVectorT<Double_t>> Simu_vect_Am={FormatVector, FormatVector, FormatVector};
    
    std::vector<TVectorT<Double_t>> Data_vect={FormatVector, FormatVector, FormatVector}; std::vector<TMatrixT<Double_t>> V_Data = {FormatMatrix, FormatMatrix, FormatMatrix};
    std::vector<TMatrixT<Double_t>> V_Simu_mu = {FormatMatrix, FormatMatrix, FormatMatrix}; std::vector<TMatrixT<Double_t>> V_Simu_gamma = {FormatMatrix, FormatMatrix, FormatMatrix}; std::vector<TMatrixT<Double_t>> V_Simu_n = {FormatMatrix, FormatMatrix, FormatMatrix}; std::vector<TMatrixT<Double_t>> V_Simu_Am = {FormatMatrix, FormatMatrix, FormatMatrix};
    
    Cut(histo_dataCOV1, V_dataCOV1, histo_dataCOV1->GetTitle(), Emin, Emax);
    Cut(histo_dataCOV2, V_dataCOV2, histo_dataCOV2->GetTitle(), Emin, Emax);
    Cut(histo_dataNTD, V_dataNTD, histo_dataNTD->GetTitle(), Emin, Emax);
    
    cout<<histo_dataCOV1->GetBinError(2)*histo_dataCOV1->GetBinError(2)<<endl;
    cout<<V_dataCOV1->GetBinContent(2, 2)<<endl;
    
    for (int i =0; i<New_Nbins; i++){
        if (i<New_Nbins_COV1){
            E_bins(i) = histo_dataCOV1->GetBinCenter(i+1);
            Data_vect[0](i) = histo_dataCOV1->GetBinContent(i+1);
            for (int j=0; j<New_Nbins_COV1; j++){
                if (j<New_Nbins_COV1){
                    V_Data[0](i,j) = V_dataCOV1->GetBinContent(i+1, j+1);
                }
                else V_Data[0](i,j)=0;
            }
        }
        else {
            Data_vect[0](i) = 0;
            V_Data[0](i,i) = 0;
        }
        if (i<New_Nbins_COV2){
            E_bins(i) = histo_dataCOV1->GetBinCenter(i+1);
            Data_vect[1](i) = histo_dataCOV2->GetBinContent(i+1);
            for (int j=0; j<New_Nbins_COV2; j++){
                if (j<New_Nbins_COV2){
                    V_Data[1](i,j) = V_dataCOV2->GetBinContent(i+1, j+1);
                }
                else V_Data[1](i,j)=0;
            }
        }
        else {
            Data_vect[1](i) = 0;
            V_Data[1](i,i) = 0;
        }
        if (i<New_Nbins_NTD){
            E_bins(i) = histo_dataCOV1->GetBinCenter(i+1);
            Data_vect[2](i) = histo_dataNTD->GetBinContent(i+1);
            for (int j=0; j<New_Nbins_NTD; j++){
                if (j<New_Nbins_NTD){
                    V_Data[2](i,j) = V_dataNTD->GetBinContent(i+1, j+1);
                }
                else V_Data[2](i,j)=0;
            }
        }
        else {
            Data_vect[2](i) = 0;
            V_Data[2](i,i) = 0;
        }
    }
    
    // Prepare simu muons vectors
    TH1D* histo_simu_muCOV1 = histo_mu[0];
    TH1D* histo_simu_muCOV2 = histo_mu[1];
    TH1D* histo_simu_muNTD = histo_mu[2];
    cout<<"open muons"<<endl;

    histo_simu_muCOV1=GoodFormat(histo_simu_muCOV1, "hMuonsCOV1", New_Nbins_COV1, Emin, Emax);
    histo_simu_muCOV2=GoodFormat(histo_simu_muCOV2, "hMuonsCOV2", New_Nbins_COV2, Emin, Emax);
    histo_simu_muNTD=GoodFormat(histo_simu_muNTD, "hMuonsNTD", New_Nbins_NTD, Emin, Emax);

    for (int i =0; i<New_Nbins; i++){
        if (i<New_Nbins_COV1){
            Simu_vect_mu[0](i) = histo_simu_muCOV1->GetBinContent(i+1);
            V_Simu_mu[0](i,i) = histo_simu_muCOV1->GetBinError(i+1)*histo_simu_muCOV1->GetBinError(i+1);
        }
        if (i<New_Nbins_COV2){
            Simu_vect_mu[1](i) = histo_simu_muCOV2->GetBinContent(i+1);
            V_Simu_mu[1](i,i) = histo_simu_muCOV2->GetBinError(i+1)*histo_simu_muCOV2->GetBinError(i+1);
        }
        if (i<New_Nbins_NTD){
            Simu_vect_mu[2](i) = histo_simu_muNTD->GetBinContent(i+1);
            V_Simu_mu[2](i,i) = histo_simu_muNTD->GetBinError(i+1)*histo_simu_muNTD->GetBinError(i+1);
        }
    }
    
    // Prepare simu gammas vectors
    TH1D* histo_simu_gammaCOV1 = histo_gamma[0];
    TH1D* histo_simu_gammaCOV2 = histo_gamma[1];
    TH1D* histo_simu_gammaNTD = histo_gamma[2];
    cout<<"open gammas"<<endl;

    histo_simu_gammaCOV1=GoodFormat(histo_simu_gammaCOV1, "hGammaCOV1", New_Nbins_COV1, Emin, Emax);
    histo_simu_gammaCOV2=GoodFormat(histo_simu_gammaCOV2, "hGammaCOV2", New_Nbins_COV2, Emin, Emax);
    histo_simu_gammaNTD=GoodFormat(histo_simu_gammaNTD, "hGammaNTD", New_Nbins_NTD, Emin, Emax);

    for (int i =0; i<New_Nbins; i++){
        if (i<New_Nbins_COV1){
            Simu_vect_gamma[0](i) = histo_simu_gammaCOV1->GetBinContent(i+1);
            V_Simu_gamma[0](i,i) = histo_simu_gammaCOV1->GetBinError(i+1)*histo_simu_gammaCOV1->GetBinError(i+1);
        }
        if (i<New_Nbins_COV2){
            Simu_vect_gamma[1](i) = histo_simu_gammaCOV2->GetBinContent(i+1);
            V_Simu_gamma[1](i,i) = histo_simu_gammaCOV2->GetBinError(i+1)*histo_simu_gammaCOV2->GetBinError(i+1);
        }
        if (i<New_Nbins_NTD){
            Simu_vect_gamma[2](i) = histo_simu_gammaNTD->GetBinContent(i+1);
            V_Simu_gamma[2](i,i) = histo_simu_gammaNTD->GetBinError(i+1)*histo_simu_gammaNTD->GetBinError(i+1);
        }
    }
    
    // Prepare simu neutrons vectors
    TH1D* histo_simu_nCOV1 = histo_n[0];
    TH1D* histo_simu_nCOV2 = histo_n[1];
    TH1D* histo_simu_nNTD = histo_n[2];

    cout<<"open neutrons"<<endl;

    histo_simu_nCOV1=GoodFormat(histo_simu_nCOV1, "hNeutronCOV1", New_Nbins_COV1, Emin, Emax);
    histo_simu_nCOV2=GoodFormat(histo_simu_nCOV2, "hNeutronCOV2", New_Nbins_COV2, Emin, Emax);
    histo_simu_nNTD=GoodFormat(histo_simu_nNTD, "hNeutronNTD", New_Nbins_NTD, Emin, Emax, true);

    for (int i =0; i<New_Nbins; i++){
        if (i<New_Nbins_COV1){
            Simu_vect_n[0](i) = histo_simu_nCOV1->GetBinContent(i+1);
            V_Simu_n[0](i,i) = histo_simu_nCOV1->GetBinError(i+1)*histo_simu_nCOV1->GetBinError(i+1);
        }
        if (i<New_Nbins_COV2){
            Simu_vect_n[1](i) = histo_simu_nCOV2->GetBinContent(i+1);
            V_Simu_n[1](i,i) = histo_simu_nCOV2->GetBinError(i+1)*histo_simu_nCOV2->GetBinError(i+1);
        }
        if (i<New_Nbins_NTD){
            Simu_vect_n[2](i) = histo_simu_nNTD->GetBinContent(i+1);
            V_Simu_n[2](i,i) = histo_simu_nNTD->GetBinError(i+1)*histo_simu_nNTD->GetBinError(i+1);
        }
    }
    
    // Prepare simu Am241 vectors
    TH1D* histo_simu_AmCOV1 = histo_Am[0];
    TH1D* histo_simu_AmCOV2 = histo_Am[1];
    TH1D* histo_simu_AmNTD = histo_Am[2];

    cout<<"open Am241"<<endl;

    histo_simu_AmCOV1=GoodFormat(histo_simu_AmCOV1, "hAm241COV1", New_Nbins_COV1, Emin, Emax);
    histo_simu_AmCOV2=GoodFormat(histo_simu_AmCOV2, "hAm241COV2", New_Nbins_COV2, Emin, Emax);
    histo_simu_AmNTD=GoodFormat(histo_simu_AmNTD, "hAm241NTD", New_Nbins_NTD, Emin, Emax);
    
    for (int i =0; i<New_Nbins; i++){
        if (i<New_Nbins_COV1){
            Simu_vect_Am[0](i) = histo_simu_AmCOV1->GetBinContent(i+1);
            V_Simu_Am[0](i,i) = histo_simu_AmCOV1->GetBinError(i+1)*histo_simu_AmCOV1->GetBinError(i+1);
        }
        if (i<New_Nbins_COV2){
            Simu_vect_Am[1](i) = histo_simu_AmCOV2->GetBinContent(i+1);
            V_Simu_Am[1](i,i) = histo_simu_AmCOV2->GetBinError(i+1)*histo_simu_AmCOV2->GetBinError(i+1);
        }
        if (i<New_Nbins_NTD){
            Simu_vect_Am[2](i) = histo_simu_AmNTD->GetBinContent(i+1);
            V_Simu_Am[2](i,i) = histo_simu_AmNTD->GetBinError(i+1)*histo_simu_AmNTD->GetBinError(i+1);
        }
    }
    
//    TCanvas *c3 = new TCanvas();
//    histo_dataNTD->SetTitle("NTD");
//    histo_dataNTD->Draw();
//    histo_dataNTD->GetXaxis()->SetRangeUser(4, 6);
//    histo_simu_nNTD->SetLineColor(kRed);
//    Hits_to_cps(histo_simu_muNTD, 1e3, MASS_LWO, NUM_EVENTS, ini_flux_mu*PLANE_SURF);
//    Hits_to_cps(histo_simu_nNTD, 1e3, MASS_LWO, NUM_EVENTS, ini_flux_n*PLANE_SURF);
//    histo_simu_nNTD->Add(histo_simu_muNTD);
//    histo_simu_nNTD->Draw("SAME");
    
    // Do the fit
    ROOT::Fit::FitResult result_Fit = Fit_Data_Simu(NumEvents, E_bins, Data_vect,
                                                    Simu_vect_mu, V_Simu_mu,
                                                    Simu_vect_gamma, V_Simu_gamma,
                                                    Simu_vect_n, V_Simu_n,
                                                    Simu_vect_Am, V_Simu_Am,
                                                    V_Data, Data_name,
                                                    ini_flux_mu, ini_flux_gamma, ini_flux_n, ini_flux_Am,
                                                    deltaE,
                                                    mu_fixed, gamma_fixed, n_fixed, Am_fixed);

    double Flux_muons = result_Fit.Parameter(0);
    double Flux_muons_err = result_Fit.ParError(0);
    double Flux_gamma = result_Fit.Parameter(1);
    double Flux_gamma_err = result_Fit.ParError(1);
    double Flux_neutrons = result_Fit.Parameter(2);
    double Flux_neutrons_err = result_Fit.ParError(2);
    double Flux_Am241 = result_Fit.Parameter(3);
    double Flux_Am241_err = result_Fit.ParError(3);


    double Chi2_Fit = result_Fit.MinFcnValue();
    int N_par = result_Fit.NPar();
    int N_dof = Data_vect[0].GetNrows()+Data_vect[1].GetNrows()+Data_vect[2].GetNrows() - N_par;

    cout << "Fit between "<< Emin << " MeV and "<<Emax<< "MeV : "<<endl;
    cout << "Chi2/dof = " << Chi2_Fit << "/" << N_dof << " = " << Chi2_Fit/N_dof<<endl;
    cout <<" flux muons = " << Flux_muons << " +- " << Flux_muons_err << " ev/cm2/s"<<endl;
    cout <<" flux gamma = " << Flux_gamma << " +- " << Flux_gamma_err << " ev/cm2/s"<<endl;
    cout << " flux neutrons = " << Flux_neutrons << " +- " << Flux_neutrons_err << " ev/cm2/s"<<endl;
    cout << " flux Am241 = " << Flux_Am241 << " +- " << Flux_Am241_err << " ev/s"<<endl;
    cout<<endl;
    cout << "**********************************************************************"<<endl;

//    TCanvas *c3 = new TCanvas();
//    histo_dataNTD->SetTitle("NTD");
//    histo_dataNTD->Draw();
//    histo_dataNTD->GetXaxis()->SetRangeUser(4, 6);
//    histo_simu_nNTD->SetLineColor(kRed);
//    Hits_to_cps(histo_simu_muNTD, 1e3, MASS_LWO, NUM_EVENTS, result_Fit.Parameter(2)*PLANE_SURF, histo_dataNTD->GetNbinsX());
//    Hits_to_cps(histo_simu_nNTD, 1e3, MASS_LWO, NUM_EVENTS, result_Fit.Parameter(2)*PLANE_SURF, histo_dataNTD->GetNbinsX());
//    histo_simu_nNTD->Add(histo_simu_muNTD);
//    histo_simu_nNTD->Draw("SAME");
    
    //Delete spectra
    histo_dataCOV1->Delete();
    histo_dataCOV2->Delete();
    histo_dataNTD->Delete();

    histo_simu_gammaCOV1->Delete();
    histo_simu_gammaCOV2->Delete();
    histo_simu_gammaNTD->Delete();

    histo_simu_muCOV1->Delete();
    histo_simu_muCOV2->Delete();
    histo_simu_muNTD->Delete();

    histo_simu_nCOV1->Delete();
    histo_simu_nCOV2->Delete();
    histo_simu_nNTD->Delete();

    histo_simu_AmCOV1->Delete();
    histo_simu_AmCOV2->Delete();
    histo_simu_AmNTD->Delete();
    
    return result_Fit;
}

//--------------------------------------------------------
// PrintChi2
void PrintChi2(TH1D* histo_data, TH1D* Model, int Npar, Double_t Min, Double_t Max){
    
    double Value_i, E_i, Diff_i, Error_i;
    double chi2=0;
    int ndf =0;
    for (int i = 1; i<histo_data->GetNbinsX(); i++){
        Value_i = histo_data->GetBinContent(i);
        E_i = histo_data->GetBinCenter(i);
        if (Value_i>0&&E_i>Min&&E_i<Max) {
            Diff_i = histo_data->GetBinContent(i)-Model->GetBinContent(i);
            Error_i = sqrt(histo_data->GetBinError(i)*histo_data->GetBinError(i)+Model->GetBinError(i)*Model->GetBinError(i));
           
            chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
            ndf+=1;
        }
    }

    cout<<"\nChi2/ndf = "<<chi2<<"/"<<ndf-Npar<<" = "<<chi2/(ndf-Npar)<<endl;
}

// --------------- Calculate covariance matrix from N spectra
TMatrixD Stat_Covariance_Matrix(TH1D* Spectrum){
    int Nbins = Spectrum->GetNbinsX();
    
    TMatrixD Cov_Mat(Nbins, Nbins);
    
    for (int i=0; i<Nbins; i++){
        Cov_Mat(i,i)=Spectrum->GetBinError(i+1)*Spectrum->GetBinError(i+1);
    }
    return Cov_Mat;
}

// ---------------- TH2D from matrix
//Get TH2D from a TMatrix
TH2D* get_TH2D_from_TMatrixD(string graph_title_, TMatrixD XY_values_, string Z_title_, string Y_title_, string X_title_, double X_min, double X_max, double Y_min, double Y_max) {
    
    auto* th2_histogram = new TH2D(graph_title_.c_str(), graph_title_.c_str(),
                                   XY_values_.GetNcols(), X_min, X_max,
                                   XY_values_.GetNrows(), Y_min, Y_max
                                   );

    for(int i_row = 0 ; i_row < XY_values_.GetNcols() ; i_row++){
        for(int j_col = 0 ; j_col < XY_values_.GetNrows() ; j_col++){
            th2_histogram->SetBinContent(i_row + 1, j_col + 1, (XY_values_)[j_col][i_row]);
        }
    }

    th2_histogram->GetXaxis()->SetTitle(X_title_.c_str());
    th2_histogram->GetYaxis()->SetTitle(Y_title_.c_str());
    th2_histogram->GetZaxis()->SetTitle(Z_title_.c_str());
    th2_histogram->SetStats(0);
    
    return th2_histogram;
    
}

// --------------- Calculate covariance matrix from N spectra
TH2D* Stat_Covariance_TH2D(TH1D* Spectrum){
    int NBins = Spectrum->GetNbinsX();
    Double_t Emin = Spectrum->GetBinLowEdge(1);
    Double_t Emax = Spectrum->GetBinLowEdge(NBins)+Spectrum->GetBinWidth(1);
    TH2D* Cov_histo = get_TH2D_from_TMatrixD(Form("CovMat_%s",Spectrum->GetTitle()), Stat_Covariance_Matrix(Spectrum), "Covariance", "Energy", "Energy", Emin, Emax, Emin, Emax);
    return Cov_histo;
}

void SetErrors_from_CovMat(TH1D* histo, TH2D* CovMat, Double_t Norm =1.){
    int NbinsHisto = histo->GetNbinsX();
    int NbinsMat = CovMat->GetNbinsX();
    
    if (NbinsHisto!=NbinsMat) cerr<<"SetErrors_from_CovMat:: histo and CovMat not compatible : histo has "<<NbinsHisto<<" bins and CovMat has : "<<NbinsMat<<" bins !"<<endl;
    
    for (int i=0; i<NbinsHisto; i++){
        histo->SetBinError(i+1, Norm*sqrt(CovMat->GetBinContent(i+1, i+1)));
    }
}

void FitNorm(TString NameFileout, TString fileData_name, TString filein_mu_name, TString filein_gamma_name, TString filein_n_name, TString filein_Am241_name, bool Syst){
    
    TFile *fileData = new TFile(fileData_name);
    
    /// --------------------------- Muons ----------------------------
    
    //prepare the vectors for the fit
    TH1D* histo_dataCOV1 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV1->SetTitle("histo_dataCOV1");
    TH1D* histo_dataCOV2 = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV2->SetTitle("histo_dataCOV2");
    TH1D* histo_dataNTD = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-20MeV");
    histo_dataNTD->SetTitle("histo_dataNTD");
    
    Double_t Num_mu =SetNumEvents(filein_mu_name);
    cout<<"\nNum events muons = "<<Num_mu<<endl;
    Double_t Num_gamma = SetNumEvents(filein_gamma_name);
    cout<<"Num events gammas = "<<Num_gamma<<endl;
    Double_t Num_n =SetNumEvents(filein_n_name);
    cout<<"Num events neutrons = "<<Num_n<<endl;
    Double_t Num_Am = SetNumEvents(filein_Am241_name);
    cout<<"Num events Am = "<<Num_Am<<endl;
    vector<Double_t> NumEvents = {Num_mu, Num_gamma, Num_n, Num_Am};
    
    TH2D* V_dataCOV1 = Stat_Covariance_TH2D(histo_dataCOV1);
    V_dataCOV1->SetTitle("V_dataCOV1");
    TH2D* V_dataCOV2= Stat_Covariance_TH2D(histo_dataCOV2);
    V_dataCOV2->SetTitle("V_dataCOV2");
    TH2D* V_dataNTD = Stat_Covariance_TH2D(histo_dataNTD);
    V_dataNTD->SetTitle("V_dataNTD");
    
    if (Syst){
        V_dataCOV1 = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-20MeV");
        V_dataCOV1->SetTitle("V_dataCOV1");
        V_dataCOV2 = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-20MeV");
        V_dataCOV2->SetTitle("V_dataCOV2");
        V_dataNTD = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-20MeV");
        V_dataNTD->SetTitle("V_dataNTD");
    }
    
    //Muons
    TFile *filein_mu = new TFile(filein_mu_name);
    
    TH1D*  histo_simu_muCOV1 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D*  histo_simu_muCOV2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D*  histo_simu_muNTD = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    //Gammas
    TFile *filein_gamma = new TFile(filein_gamma_name);
    
    TH1D* histo_simu_gammaCOV1 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_gammaCOV2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_gammaNTD = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    //Neutrons
    TFile *filein_n = new TFile(filein_n_name);
    
    TH1D* histo_simu_nCOV1 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_nCOV2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_nNTD = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");
    histo_simu_nNTD= Shift(histo_simu_nNTD, HQF);
    
    //Am241
    TFile *filein_Am = new TFile(filein_Am241_name);
    
    TH1D* histo_simu_AmCOV1 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_AmCOV2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_AmNTD = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");
    
    std::vector<TH1D*> histos_Data = {histo_dataCOV1, histo_dataCOV2, histo_dataNTD};
    std::vector<TH2D*> Variances_data = {V_dataCOV1, V_dataCOV2, V_dataNTD};
    std::vector<TH1D*> histos_mu = {histo_simu_muCOV1, histo_simu_muCOV2, histo_simu_muNTD};
    std::vector<TH1D*> histos_gamma = {histo_simu_gammaCOV1, histo_simu_gammaCOV2, histo_simu_gammaNTD};
    std::vector<TH1D*> histos_n = {histo_simu_nCOV1, histo_simu_nCOV2, histo_simu_nNTD};
    std::vector<TH1D*> histos_Am = {histo_simu_AmCOV1, histo_simu_AmCOV2, histo_simu_AmNTD};
    
    //Fit
    ROOT::Fit::FitResult fit_muon = Fit_OnRange(NumEvents, histos_Data, Variances_data, histos_mu, histos_gamma, histos_n, histos_Am, 9, 18,  0.018, 0, 0, 0);

    TCanvas *can= new TCanvas("Muons","Muons", 1000, 1000);
    can->Divide(1,3);
    
    cout<<"Prepare plot Muons"<<endl;
    
    histo_dataCOV1 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV1->SetTitle("histo_dataCOV1");
    histo_dataCOV2 = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV2->SetTitle("histo_dataCOV2");
    histo_dataNTD = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-20MeV");
    histo_dataNTD->SetTitle("histo_dataNTD");
    
    if (Syst){
        V_dataCOV1 = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-20MeV");
        V_dataCOV1->SetTitle("V_dataCOV1");
        V_dataCOV2 = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-20MeV");
        V_dataCOV2->SetTitle("V_dataCOV2");
        V_dataNTD = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-20MeV");
        V_dataNTD->SetTitle("V_dataNTD");
        
        SetErrors_from_CovMat(histo_dataCOV1, V_dataCOV1);
        SetErrors_from_CovMat(histo_dataCOV2, V_dataCOV2);
        SetErrors_from_CovMat(histo_dataNTD, V_dataNTD);
    }
    
    //COV1
    can->cd(1);
    TH1D* histo_dataCOV1_muon = (TH1D*) histo_dataCOV1->Clone();
    histo_dataCOV1_muon->SetTitle("COV1");
    gPad->SetRightMargin(0.02);
    histo_dataCOV1_muon->SetTitleSize(0.08);
    histo_dataCOV1_muon->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV1_muon->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV1_muon->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV1_muon->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV1_muon->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV1_muon->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV1_muon->Draw();
    histo_dataCOV1_muon->GetXaxis()->SetRangeUser(6, 18);
    histo_simu_muCOV1=GoodFormat(histo_simu_muCOV1, "histo_simu_muCOV1", histo_dataCOV1->GetNbinsX(), 0, 20);
    histo_simu_muCOV1->SetLineColor(kRed);
    Hits_to_cps(histo_simu_muCOV1, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV1->GetNbinsX());
    histo_simu_muCOV1->Draw("SAME");

    //COV2
    can->cd(2);
    TH1D* histo_dataCOV2_muon = (TH1D*) histo_dataCOV2->Clone();
    histo_dataCOV2_muon->SetTitle("COV2");
    gPad->SetRightMargin(0.02);
    histo_dataCOV2_muon->SetTitleSize(0.08);
    histo_dataCOV2_muon->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV2_muon->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV2_muon->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV2_muon->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV2_muon->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV2_muon->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV2_muon->Draw();
    histo_dataCOV2_muon->GetXaxis()->SetRangeUser(6, 18);
    histo_simu_muCOV2=GoodFormat(histo_simu_muCOV2, "histo_simu_muCOV2", histo_dataCOV2->GetNbinsX(), 0, 20);
    histo_simu_muCOV2->SetLineColor(kRed);
    Hits_to_cps(histo_simu_muCOV2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV2->GetNbinsX());
    histo_simu_muCOV2->Draw("SAME");

    //LWO
    can->cd(3);
    TH1D* histo_dataNTD_muon = (TH1D*) histo_dataNTD->Clone();
    histo_dataNTD_muon->SetTitle("NTD");
    gPad->SetRightMargin(0.02);
    histo_dataNTD_muon->SetTitleSize(0.08);
    histo_dataNTD_muon->GetXaxis()->SetTitleSize(0.06);
    histo_dataNTD_muon->GetXaxis()->SetLabelSize(0.055);
    histo_dataNTD_muon->GetXaxis()->SetTitleOffset(0.75);
    histo_dataNTD_muon->GetYaxis()->SetTitleSize(0.06);
    histo_dataNTD_muon->GetYaxis()->SetLabelSize(0.06);
    histo_dataNTD_muon->GetYaxis()->SetTitleOffset(0.6);
    histo_dataNTD_muon->Draw();
    histo_dataNTD_muon->GetXaxis()->SetRangeUser(6, 18);
    histo_simu_muNTD=GoodFormat(histo_simu_muNTD, "histo_simu_muNTD", histo_dataNTD->GetNbinsX(), 0, 20);
    histo_simu_muNTD->SetLineColor(kRed);
    Hits_to_cps(histo_simu_muNTD, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD->GetNbinsX());
    histo_simu_muNTD->Draw("SAME");
    
    cout<<"Prepare fit Neutrons"<<endl;
    
    /// --------------------------- Neutrons ----------------------------
    //prepare the vectors for the fit
    TH1D* histo_dataCOV1_10MeV = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV1_10MeV->SetTitle("histo_dataCOV1_10MeV");
    TH1D* histo_dataCOV2_10MeV = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV2_10MeV->SetTitle("histo_dataCOV2_10MeV");
    TH1D* histo_dataNTD_10MeV = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_10MeV->SetTitle("histo_dataNTD_10MeV");
    
    TH2D* V_dataCOV1_10MeV = Stat_Covariance_TH2D(histo_dataCOV1_10MeV);
    V_dataCOV1_10MeV->SetTitle("V_dataCOV1_10MeV");
    TH2D* V_dataCOV2_10MeV = Stat_Covariance_TH2D(histo_dataCOV2_10MeV);
    V_dataCOV2_10MeV->SetTitle("V_dataCOV2_10MeV");
    TH2D* V_dataNTD_10MeV = Stat_Covariance_TH2D(histo_dataNTD_10MeV);
    V_dataNTD_10MeV->SetTitle("V_dataNTD_10MeV");
    
    if (Syst){
        V_dataCOV1_10MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-10MeV");
        V_dataCOV1_10MeV->SetTitle("V_dataNTD_10MeV");
        V_dataCOV2_10MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-10MeV");
        V_dataCOV2_10MeV->SetTitle("V_dataCOV1_10MeV");
        V_dataNTD_10MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataNTD_10MeV->SetTitle("V_dataNTD_10MeV");
    }
    
//    histo_dataCOV1_10MeV->Rebin(2.);
//    histo_dataCOV2_10MeV->Rebin(2.);
//    histo_dataNTD_10MeV->Rebin(2.);
//
//    histo_dataCOV1_10MeV->Scale(1/2.);
//    histo_dataCOV2_10MeV->Scale(1/2.);
//    histo_dataNTD_10MeV->Scale(1/2.);
    
    TH1D*  histo_simu_muCOV1_10MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D*  histo_simu_muCOV2_10MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D*  histo_simu_muNTD_10MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_gammaCOV1_10MeV = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_gammaCOV2_10MeV = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_gammaNTD_10MeV= (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_nCOV1_10MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_nCOV2_10MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_nNTD_10MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_AmCOV1_10MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_AmCOV2_10MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_AmNTD_10MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    std::vector<TH1D*> histos_Data_10MeV = {histo_dataCOV1_10MeV, histo_dataCOV2_10MeV, histo_dataNTD_10MeV};
    std::vector<TH2D*> Variances_data_10MeV = {V_dataCOV1_10MeV, V_dataCOV2_10MeV, V_dataNTD_10MeV};
    std::vector<TH1D*> histos_mu_10MeV = {histo_simu_muCOV1_10MeV, histo_simu_muCOV2_10MeV, histo_simu_muNTD_10MeV};
    std::vector<TH1D*> histos_gamma_10MeV = {histo_simu_gammaCOV1_10MeV, histo_simu_gammaCOV2_10MeV, histo_simu_gammaNTD_10MeV};
    std::vector<TH1D*> histos_n_10MeV = {histo_simu_nCOV1_10MeV, histo_simu_nCOV2_10MeV, histo_simu_nNTD_10MeV};
    std::vector<TH1D*> histos_Am_10MeV = {histo_simu_AmCOV1_10MeV, histo_simu_AmCOV2_10MeV, histo_simu_AmNTD_10MeV};

    //Fit
    ROOT::Fit::FitResult fit_neutron = Fit_OnRange(NumEvents, histos_Data_10MeV, Variances_data_10MeV, histos_mu_10MeV, histos_gamma_10MeV, histos_n_10MeV, histos_Am_10MeV, 4.5, 5.5, fit_muon.Parameter(0), 0, 0.015, 0, true, false, false);

    //Plot
    TCanvas *can2= new TCanvas("Neutrons","Neutrons", 1000, 1000);
    can2->Divide(1,3);
    
    histo_dataCOV1_10MeV = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV1_10MeV->SetTitle("histo_dataCOV1_10MeV");
    histo_dataCOV2_10MeV = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV2_10MeV->SetTitle("histo_dataCOV2_10MeV");
    histo_dataNTD_10MeV = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_10MeV->SetTitle("histo_dataNTD_10MeV");
    
    if (Syst){
        V_dataCOV1_10MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-10MeV");
        V_dataCOV1_10MeV->SetTitle("V_dataCOV1_10MeV");
        V_dataCOV2_10MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-10MeV");
        V_dataCOV2_10MeV->SetTitle("V_dataCOV2_10MeV");
        V_dataNTD_10MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataNTD_10MeV->SetTitle("V_dataNTD_10MeV");
        
        SetErrors_from_CovMat(histo_dataCOV1_10MeV, V_dataCOV1_10MeV);
        SetErrors_from_CovMat(histo_dataCOV2_10MeV, V_dataCOV2_10MeV);
        SetErrors_from_CovMat(histo_dataNTD_10MeV, V_dataNTD_10MeV);
    }
    
    //COV1
    can2->cd(1);
    TH1D* histo_dataCOV1_neutrons = (TH1D*) histo_dataCOV1_10MeV->Clone();
    histo_dataCOV1_neutrons->SetTitle("COV1");
    gPad->SetRightMargin(0.02);
    histo_dataCOV1_neutrons->SetTitleSize(0.08);
    histo_dataCOV1_neutrons->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV1_neutrons->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV1_neutrons->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV1_neutrons->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV1_neutrons->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV1_neutrons->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV1_neutrons->Draw();
    histo_dataCOV1_neutrons->GetXaxis()->SetRangeUser(4.5, 5.5);
    histo_simu_muCOV1_10MeV=GoodFormat(histo_simu_muCOV1_10MeV, "histo_simu_muCOV1_10MeV", histo_dataCOV1_10MeV->GetNbinsX(), 0, 10);
    histo_simu_nCOV1_10MeV=GoodFormat(histo_simu_nCOV1_10MeV, "histo_simu_nCOV1_10MeV", histo_dataCOV1_10MeV->GetNbinsX(), 0, 10);
    histo_simu_nCOV1_10MeV->SetLineColor(kRed);
    Hits_to_cps(histo_simu_muCOV1_10MeV, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_10MeV->GetNbinsX());
    Hits_to_cps(histo_simu_nCOV1_10MeV, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataNTD_10MeV->GetNbinsX());
    histo_simu_nCOV1_10MeV->Add(histo_simu_muCOV1_10MeV);
    histo_simu_nCOV1_10MeV->Draw("SAME");

    //COV2
    can2->cd(2);
    TH1D* histo_dataCOV2_neutrons = (TH1D*) histo_dataCOV2_10MeV->Clone();
    histo_dataCOV2_neutrons->SetTitle("COV2");
    gPad->SetRightMargin(0.02);
    histo_dataCOV2_neutrons->SetTitleSize(0.08);
    histo_dataCOV2_neutrons->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV2_neutrons->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV2_neutrons->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV2_neutrons->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV2_neutrons->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV2_neutrons->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV2_neutrons->Draw();
    histo_dataCOV2_neutrons->GetXaxis()->SetRangeUser(4.5, 5.5);
    histo_simu_muCOV2_10MeV=GoodFormat(histo_simu_muCOV2_10MeV, "histo_simu_muCOV2_10MeV", histo_dataCOV2_10MeV->GetNbinsX(), 0, 10);
    histo_simu_nCOV2_10MeV=GoodFormat(histo_simu_nCOV2_10MeV, "histo_simu_nCOV2_10MeV", histo_dataCOV2_10MeV->GetNbinsX(), 0, 10);
    histo_simu_nCOV2_10MeV->SetLineColor(kRed);
    Hits_to_cps(histo_simu_muCOV2_10MeV, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_10MeV->GetNbinsX());
    Hits_to_cps(histo_simu_nCOV2_10MeV, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF, histo_dataNTD_10MeV->GetNbinsX());
    histo_simu_nCOV2_10MeV->Add(histo_simu_muCOV2_10MeV);
    histo_simu_nCOV2_10MeV->Draw("SAME");

    //COV3
    can2->cd(3);
    TH1D* histo_dataNTD_neutrons = (TH1D*) histo_dataNTD_10MeV->Clone();
    histo_dataNTD_neutrons->SetTitle("NTD");
    gPad->SetRightMargin(0.02);
    histo_dataNTD_neutrons->SetTitleSize(0.08);
    histo_dataNTD_neutrons->GetXaxis()->SetTitleSize(0.06);
    histo_dataNTD_neutrons->GetXaxis()->SetLabelSize(0.055);
    histo_dataNTD_neutrons->GetXaxis()->SetTitleOffset(0.75);
    histo_dataNTD_neutrons->GetYaxis()->SetTitleSize(0.06);
    histo_dataNTD_neutrons->GetYaxis()->SetLabelSize(0.06);
    histo_dataNTD_neutrons->GetYaxis()->SetTitleOffset(0.6);
    histo_dataNTD_neutrons->Draw();
    histo_dataNTD_neutrons->GetXaxis()->SetRangeUser(4.5, 5.5);
    histo_simu_muNTD_10MeV=GoodFormat(histo_simu_muNTD_10MeV, "histo_simu_muNTD_10MeV", histo_dataNTD_10MeV->GetNbinsX(), 0, 10);
    histo_simu_nNTD_10MeV=GoodFormat(histo_simu_nNTD_10MeV, "histo_simu_nNTD_10MeV", histo_dataNTD_10MeV->GetNbinsX(), 0, 10, true);
    Hits_to_cps(histo_simu_muNTD_10MeV, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_10MeV->GetNbinsX());
    Hits_to_cps(histo_simu_nNTD_10MeV, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF, histo_dataNTD_10MeV->GetNbinsX());
    histo_simu_nNTD_10MeV->SetLineColor(kRed);
    histo_simu_nNTD_10MeV->Add(histo_simu_muNTD_10MeV);
    histo_simu_nNTD_10MeV->Draw("SAME");
    
    cout<<"Fit gammas..."<<endl;
    /// --------------------------- Gammas ----------------------------
    //prepare the vectors for the fit
    TH1D* histo_dataCOV1_4MeV = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV1_4MeV->SetTitle("histo_dataCOV1_4MeV");
    TH1D* histo_dataCOV2_4MeV = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV2_4MeV->SetTitle("histo_dataCOV2_4MeV");
    TH1D* histo_dataNTD_4MeV = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_4MeV->SetTitle("histo_dataNTD_4MeV");

    TH2D* V_dataCOV1_4MeV = Stat_Covariance_TH2D(histo_dataCOV1_4MeV);
    V_dataCOV1_4MeV->SetTitle("V_dataCOV1_4MeV");
    TH2D* V_dataCOV2_4MeV= Stat_Covariance_TH2D(histo_dataCOV2_4MeV);
    V_dataCOV2_4MeV->SetTitle("V_dataCOV2_4MeV");
    TH2D* V_dataNTD_4MeV = Stat_Covariance_TH2D(histo_dataNTD_4MeV);
    V_dataNTD_4MeV->SetTitle("V_dataNTD_4MeV");
    
    if (Syst){
        V_dataCOV1_4MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-4MeV");
        V_dataCOV1_4MeV->SetTitle("V_dataCOV1_4MeV");
        V_dataCOV2_4MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-4MeV");
        histo_dataCOV2_4MeV->SetTitle("V_dataCOV2_4MeV");
        V_dataCOV1_4MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataCOV1_4MeV->SetTitle("V_dataNTD_4MeV");
    }
    
    TH1D*  histo_simu_muCOV1_4MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D*  histo_simu_muCOV2_4MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D*  histo_simu_muNTD_4MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_gammaCOV1_4MeV = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_gammaCOV2_4MeV = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_gammaNTD_4MeV= (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_nCOV1_4MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_nCOV2_4MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_nNTD_4MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_AmCOV1_4MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_AmCOV2_4MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_AmNTD_4MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

//    histo_dataCOV1_4MeV=GoodFormat(histo_dataCOV1_4MeV, "histo_dataCOV1_4MeV", 200, 0, 4);
//    histo_dataCOV2_4MeV=GoodFormat(histo_dataCOV2_4MeV, "histo_dataCOV2_4MeV", 200, 0, 4);
//    histo_dataNTD_4MeV=GoodFormat(histo_dataNTD_4MeV, "histo_dataNTD_4MeV", 200, 0, 4);

    int NBins_COV1_4MeV =histo_dataCOV1_4MeV->FindBin(4.)-histo_dataCOV1_4MeV->FindBin(0);
    int NBins_COV2_4MeV =histo_dataCOV2_4MeV->FindBin(4.)-histo_dataCOV2_4MeV->FindBin(0);
    int NBins_NTD_4MeV =histo_dataNTD_4MeV->FindBin(4.)-histo_dataNTD_4MeV->FindBin(0);
    
    histo_simu_muCOV1_4MeV=GoodFormat(histo_simu_muCOV1_4MeV, "histo_simu_muCOV1_4MeV", NBins_COV1_4MeV, 0, 4);
    histo_simu_muCOV2_4MeV=GoodFormat(histo_simu_muCOV2_4MeV, "histo_simu_muCOV2_4MeV", NBins_COV2_4MeV, 0, 4);
    histo_simu_muNTD_4MeV=GoodFormat(histo_simu_muNTD_4MeV, "histo_simu_muNTD_4MeV", NBins_NTD_4MeV, 0, 4);

    histo_simu_gammaCOV1_4MeV=GoodFormat(histo_simu_gammaCOV1_4MeV, "histo_simu_gammaCOV1_4MeV", NBins_COV1_4MeV, 0, 4);
    histo_simu_gammaCOV2_4MeV=GoodFormat(histo_simu_gammaCOV2_4MeV, "histo_simu_gammaCOV2_4MeV", NBins_COV2_4MeV, 0, 4);
    histo_simu_gammaNTD_4MeV=GoodFormat(histo_simu_gammaNTD_4MeV, "histo_simu_gammaNTD_4MeV", NBins_NTD_4MeV, 0, 4);

    histo_simu_nCOV1_4MeV=GoodFormat(histo_simu_nCOV1_4MeV, "histo_simu_nCOV1_4MeV", NBins_COV1_4MeV, 0, 4);
    histo_simu_nCOV2_4MeV=GoodFormat(histo_simu_nCOV2_4MeV, "histo_simu_nCOV2_4MeV", NBins_COV2_4MeV, 0, 4);
    histo_simu_nNTD_4MeV=GoodFormat(histo_simu_nNTD_4MeV, "histo_simu_nNTD_4MeV", NBins_NTD_4MeV, 0, 4);

    histo_simu_AmCOV1_4MeV=GoodFormat(histo_simu_AmCOV1_4MeV, "histo_simu_AmCOV1_4MeV", NBins_COV1_4MeV, 0, 4);
    histo_simu_AmCOV2_4MeV=GoodFormat(histo_simu_AmCOV2_4MeV, "histo_simu_AmCOV2_4MeV", NBins_COV2_4MeV, 0, 4);
    histo_simu_AmNTD_4MeV=GoodFormat(histo_simu_AmNTD_4MeV, "histo_simu_AmNTD_4MeV", NBins_NTD_4MeV, 0, 4);

    std::vector<TH1D*> histos_Data_4MeV = {histo_dataCOV1_4MeV, histo_dataCOV2_4MeV, histo_dataNTD_4MeV};
    std::vector<TH2D*> Variances_data_4MeV = {V_dataCOV1_4MeV, V_dataCOV2_4MeV, V_dataNTD_4MeV};
    std::vector<TH1D*> histos_mu_4MeV = {histo_simu_muCOV1_4MeV, histo_simu_muCOV2_4MeV, histo_simu_muNTD_4MeV};
    std::vector<TH1D*> histos_gamma_4MeV = {histo_simu_gammaCOV1_4MeV, histo_simu_gammaCOV2_4MeV, histo_simu_gammaNTD_4MeV};
    std::vector<TH1D*> histos_n_4MeV = {histo_simu_nCOV1_4MeV, histo_simu_nCOV2_4MeV, histo_simu_nNTD_4MeV};
    std::vector<TH1D*> histos_Am_4MeV = {histo_simu_AmCOV1_4MeV, histo_simu_AmCOV2_4MeV, histo_simu_AmNTD_4MeV};

    //Fit
    ROOT::Fit::FitResult fit_gamma = Fit_OnRange(NumEvents, histos_Data_4MeV, Variances_data_4MeV, histos_mu_4MeV, histos_gamma_4MeV, histos_n_4MeV, histos_Am_4MeV, 0.1, 4., fit_muon.Parameter(0), 3.3, fit_neutron.Parameter(2), 0, true, false, true);

    //Plot
    TCanvas *can3= new TCanvas("Gammas","Gammas", 1400, 600);
    can3->Divide(1,3);
    
    histo_dataCOV1_4MeV = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-4MeV");
    histo_dataCOV1_4MeV->SetTitle("histo_dataCOV1_4MeV");
    histo_dataCOV1_4MeV->SetName("histo_dataCOV1_4MeV");
    histo_dataCOV2_4MeV = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-4MeV");
    histo_dataCOV2_4MeV->SetTitle("histo_dataCOV2_4MeV");
    histo_dataCOV2_4MeV->SetName("histo_dataCOV2_4MeV");
    histo_dataNTD_4MeV = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_4MeV->SetTitle("histo_dataNTD_4MeV");
    histo_dataNTD_4MeV->SetName("histo_dataNTD_4MeV");
    
    if (Syst){
        V_dataCOV1_4MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-4MeV");
        V_dataCOV1_4MeV->SetTitle("V_dataCOV1_4MeV");
        V_dataCOV2_4MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-4MeV");
        V_dataCOV2_4MeV->SetTitle("V_dataCOV2_4MeV");
        V_dataNTD_4MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataNTD_4MeV->SetTitle("V_dataNTD_4MeV");
        
        SetErrors_from_CovMat(histo_dataCOV1_4MeV, V_dataCOV1_4MeV);
        SetErrors_from_CovMat(histo_dataCOV2_4MeV, V_dataCOV2_4MeV);
        SetErrors_from_CovMat(histo_dataNTD_4MeV, V_dataNTD_4MeV);
    }
    
//    histo_dataCOV1_4MeV=GoodFormat(histo_dataCOV1_4MeV, "histo_dataCOV1_4MeV", 200, 0, 4);
//    histo_dataCOV2_4MeV=GoodFormat(histo_dataCOV2_4MeV, "histo_dataCOV2_4MeV", 200, 0, 4);
//    histo_dataNTD_4MeV=GoodFormat(histo_dataNTD_4MeV, "histo_dataNTD_4MeV", 200, 0, 4);
    
    //COV1
    can3->cd(1);
    gPad->SetLogy();
    histo_dataCOV1_4MeV->SetTitle("COV1");
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.9);
    histo_dataCOV1_4MeV->SetTitleSize(0.08);
    histo_dataCOV1_4MeV->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV1_4MeV->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV1_4MeV->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV1_4MeV->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV1_4MeV->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV1_4MeV->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV1_4MeV->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1_4MeV->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV1_4MeV->Draw();
    histo_dataCOV1_4MeV->GetXaxis()->SetRangeUser(0, 4);
    Hits_to_cps(histo_simu_muCOV1_4MeV, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, NBins_COV1_4MeV);
    Hits_to_cps(histo_simu_nCOV1_4MeV, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,NBins_COV1_4MeV);
    Hits_to_cps(histo_simu_gammaCOV1_4MeV, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,NBins_COV1_4MeV);
    histo_simu_gammaCOV1_4MeV->SetLineColor(kRed);
    histo_simu_gammaCOV1_4MeV->Add(histo_simu_muCOV1_4MeV);
    histo_simu_gammaCOV1_4MeV->Add(histo_simu_nCOV1_4MeV);
    histo_simu_gammaCOV1_4MeV->Draw("SAME");
    
    //COV2
    can3->cd(2);
    gPad->SetLogy();
    histo_dataCOV2_4MeV->SetTitle("COV2");
    gPad->SetRightMargin(0.02);
    histo_dataCOV2_4MeV->SetTitleSize(0.08);
    histo_dataCOV2_4MeV->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV2_4MeV->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV2_4MeV->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV2_4MeV->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV2_4MeV->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV2_4MeV->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV2_4MeV->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2_4MeV->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2_4MeV->Draw();
    histo_dataCOV2_4MeV->GetXaxis()->SetRangeUser(0, 4);
    Hits_to_cps(histo_simu_muCOV2_4MeV, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, NBins_COV2_4MeV);
    Hits_to_cps(histo_simu_nCOV2_4MeV, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,NBins_COV2_4MeV);
    Hits_to_cps(histo_simu_gammaCOV2_4MeV, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,NBins_COV2_4MeV);
    histo_simu_gammaCOV2_4MeV->SetLineColor(kRed);
    histo_simu_gammaCOV2_4MeV->Add(histo_simu_muCOV2_4MeV);
    histo_simu_gammaCOV2_4MeV->Add(histo_simu_nCOV2_4MeV);
    histo_simu_gammaCOV2_4MeV->Draw("SAME");

    //COV3
    can3->cd(3);
    gPad->SetLogy();
    histo_dataNTD_4MeV->SetTitle("NTD");
    gPad->SetRightMargin(0.02);
    histo_dataNTD_4MeV->SetTitleSize(0.08);
    histo_dataNTD_4MeV->GetXaxis()->SetTitleSize(0.06);
    histo_dataNTD_4MeV->GetXaxis()->SetLabelSize(0.055);
    histo_dataNTD_4MeV->GetXaxis()->SetTitleOffset(0.75);
    histo_dataNTD_4MeV->GetYaxis()->SetTitleSize(0.06);
    histo_dataNTD_4MeV->GetYaxis()->SetLabelSize(0.06);
    histo_dataNTD_4MeV->GetYaxis()->SetTitleOffset(0.6);
    histo_dataNTD_4MeV->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataNTD_4MeV->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataNTD_4MeV->Draw();
    histo_dataNTD_4MeV->GetXaxis()->SetRangeUser(0, 4);
    histo_dataNTD_4MeV->GetYaxis()->SetRangeUser(1e-5, 1);
    Hits_to_cps(histo_simu_muNTD_4MeV, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, NBins_NTD_4MeV);
    Hits_to_cps(histo_simu_nNTD_4MeV, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,NBins_NTD_4MeV);
    Hits_to_cps(histo_simu_gammaNTD_4MeV, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,NBins_NTD_4MeV);
    TH1D* histo_simu_sumNTD_4MeV =(TH1D*) histo_simu_gammaNTD_4MeV->Clone();
    histo_simu_sumNTD_4MeV->SetLineColor(kRed);
    histo_simu_sumNTD_4MeV->Add(histo_simu_muNTD_4MeV);
    histo_simu_sumNTD_4MeV->Add(histo_simu_nNTD_4MeV);
    histo_simu_sumNTD_4MeV->Draw("SAME");
    
    /// --------------------------- Am241 ----------------------------
    cout<<"Fit Am241..."<<endl;
    //prepare the vectors for the fit
    TH1D* histo_dataCOV1_1MeV = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-1MeV");
    histo_dataCOV1_1MeV->SetTitle("histo_dataCOV1_1MeV");
    histo_dataCOV1_1MeV->SetName("histo_dataCOV1_1MeV");
    TH1D* histo_dataCOV2_1MeV = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-1MeV");
    histo_dataCOV2_1MeV->SetTitle("histo_dataCOV2_1MeV");
    histo_dataCOV2_1MeV->SetName("histo_dataCOV2_1MeV");
    TH1D* histo_dataNTD_1MeV = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_1MeV->SetTitle("histo_dataNTD_1MeV");
    histo_dataNTD_1MeV->SetName("histo_dataNTD_1MeV");
    
    TH2D* V_dataCOV1_1MeV = Stat_Covariance_TH2D(histo_dataCOV1_1MeV);
    V_dataCOV1_1MeV->SetTitle("V_dataCOV1_1MeV");
    TH2D* V_dataCOV2_1MeV= Stat_Covariance_TH2D(histo_dataCOV2_1MeV);
    V_dataCOV2_1MeV->SetTitle("V_dataCOV2_1MeV");
    TH2D* V_dataNTD_1MeV = Stat_Covariance_TH2D(histo_dataNTD_1MeV);
    V_dataNTD_1MeV->SetTitle("V_dataNTD_1MeV");
    
    if (Syst){
        V_dataCOV1_1MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-1MeV");
        V_dataCOV1_1MeV->SetTitle("V_dataCOV1_1MeV");
        V_dataCOV1_1MeV->SetName("V_dataCOV1_1MeV");
        V_dataCOV2_1MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-1MeV");
        V_dataCOV2_1MeV->SetTitle("V_dataCOV2_1MeV");
        V_dataCOV2_1MeV->SetName("V_dataCOV2_1MeV");
        V_dataNTD_1MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataNTD_1MeV->SetTitle("V_dataNTD_1MeV");
        V_dataNTD_1MeV->SetName("V_dataNTD_1MeV");
    }
    
    TH1D*  histo_simu_muCOV1_1MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-1MeV");
    TH1D*  histo_simu_muCOV2_1MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-1MeV");
    TH1D*  histo_simu_muNTD_1MeV = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-1MeV");

    TH1D* histo_simu_gammaCOV1_1MeV = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-1MeV");
    TH1D* histo_simu_gammaCOV2_1MeV = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-1MeV");
    TH1D* histo_simu_gammaNTD_1MeV= (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-1MeV");

    TH1D* histo_simu_nCOV1_1MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-1MeV");
    TH1D* histo_simu_nCOV2_1MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-1MeV");
    TH1D* histo_simu_nNTD_1MeV = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-1MeV");

    TH1D* histo_simu_AmCOV1_1MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-1MeV");
    TH1D* histo_simu_AmCOV2_1MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-1MeV");
    TH1D* histo_simu_AmNTD_1MeV = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-1MeV");

    int NBins_COV1_1MeV =histo_dataCOV1_1MeV->FindBin(1.)-histo_dataCOV1_1MeV->FindBin(0);
    int NBins_COV2_1MeV =histo_dataCOV2_1MeV->FindBin(1.)-histo_dataCOV2_1MeV->FindBin(0);
    int NBins_NTD_1MeV =histo_dataNTD_1MeV->FindBin(1.)-histo_dataNTD_1MeV->FindBin(0);
    
    
    histo_simu_muCOV1_1MeV=GoodFormat(histo_simu_muCOV1_1MeV, "histo_simu_muCOV1_1MeV", NBins_COV1_1MeV, 0, 1);
    histo_simu_muCOV2_1MeV=GoodFormat(histo_simu_muCOV2_1MeV, "histo_simu_muCOV2_1MeV", NBins_COV2_1MeV, 0, 1);
    histo_simu_muNTD_1MeV=GoodFormat(histo_simu_muNTD_1MeV, "histo_simu_muNTD_1MeV", NBins_NTD_1MeV, 0, 1);

    histo_simu_gammaCOV1_1MeV=GoodFormat(histo_simu_gammaCOV1_1MeV, "histo_simu_gammaCOV1_1MeV", NBins_COV1_1MeV, 0, 1);
    histo_simu_gammaCOV2_1MeV=GoodFormat(histo_simu_gammaCOV2_1MeV, "histo_simu_gammaCOV2_1MeV", NBins_COV2_1MeV, 0, 1);
    histo_simu_gammaNTD_1MeV=GoodFormat(histo_simu_gammaNTD_1MeV, "histo_simu_gammaNTD_1MeV", NBins_NTD_1MeV, 0, 1);

    histo_simu_nCOV1_1MeV=GoodFormat(histo_simu_nCOV1_1MeV, "histo_simu_nCOV1_1MeV", NBins_COV1_1MeV, 0, 1);
    histo_simu_nCOV2_1MeV=GoodFormat(histo_simu_nCOV2_1MeV, "histo_simu_nCOV2_1MeV", NBins_COV2_1MeV, 0, 1);
    histo_simu_nNTD_1MeV=GoodFormat(histo_simu_nNTD_1MeV, "histo_simu_nNTD_1MeV", NBins_NTD_1MeV, 0, 1);

    histo_simu_AmCOV1_1MeV=GoodFormat(histo_simu_AmCOV1_1MeV, "histo_simu_AmCOV1_1MeV", NBins_COV1_1MeV, 0, 1);
    histo_simu_AmCOV2_1MeV=GoodFormat(histo_simu_AmCOV2_1MeV, "histo_simu_AmCOV2_1MeV", NBins_COV2_1MeV, 0, 1);
    histo_simu_AmNTD_1MeV=GoodFormat(histo_simu_AmNTD_1MeV, "histo_simu_AmNTD_1MeV", NBins_NTD_1MeV, 0, 1);

    std::vector<TH1D*> histos_Data_1MeV = {histo_dataCOV1_1MeV, histo_dataCOV2_1MeV, histo_dataNTD_1MeV};
    std::vector<TH2D*> Variances_data_1MeV = {V_dataCOV1_1MeV, V_dataCOV2_1MeV, V_dataNTD_1MeV};
    std::vector<TH1D*> histos_mu_1MeV = {histo_simu_muCOV1_1MeV, histo_simu_muCOV2_1MeV, histo_simu_muNTD_1MeV};
    std::vector<TH1D*> histos_gamma_1MeV = {histo_simu_gammaCOV1_1MeV, histo_simu_gammaCOV2_1MeV, histo_simu_gammaNTD_1MeV};
    std::vector<TH1D*> histos_n_1MeV = {histo_simu_nCOV1_1MeV, histo_simu_nCOV2_1MeV, histo_simu_nNTD_1MeV};
    std::vector<TH1D*> histos_Am_1MeV = {histo_simu_AmCOV1_1MeV, histo_simu_AmCOV2_1MeV, histo_simu_AmNTD_1MeV};

    //Fit
    ROOT::Fit::FitResult fit_Am = Fit_OnRange(NumEvents, histos_Data_1MeV, Variances_data_1MeV, histos_mu_1MeV, histos_gamma_1MeV, histos_n_1MeV, histos_Am_1MeV, 0.05, 0.1, fit_muon.Parameter(0), fit_gamma.Parameter(1), fit_neutron.Parameter(2), 0.02, true, true, true, false);

    //Plot
    TCanvas *can4= new TCanvas("Am","Am", 1400, 600);
    can4->Divide(1,3);

    TH1D* histo_dataCOV1_1MeV_plot = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-1MeV");
    histo_dataCOV1_1MeV_plot->SetName("histo_dataCOV1_1MeV_plot");
    TH1D* histo_dataCOV2_1MeV_plot = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-1MeV");
    histo_dataCOV2_1MeV_plot->SetName("histo_dataCOV2_1MeV_plot");
    TH1D* histo_dataNTD_1MeV_plot = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_1MeV_plot->SetName("histo_dataNTD_1MeV_plot");
    
    if (Syst){
        V_dataCOV1_1MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-1MeV");
        V_dataCOV1_1MeV->SetTitle("V_dataCOV1_1MeV");
        V_dataCOV2_1MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-1MeV");
        V_dataCOV2_1MeV->SetTitle("V_dataCOV2_1MeV");
        V_dataNTD_1MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataNTD_1MeV->SetTitle("V_dataNTD_1MeV");
        
        SetErrors_from_CovMat(histo_dataCOV1_1MeV_plot, V_dataCOV1_1MeV);
        SetErrors_from_CovMat(histo_dataCOV2_1MeV_plot, V_dataCOV2_1MeV);
        SetErrors_from_CovMat(histo_dataNTD_1MeV_plot, V_dataNTD_1MeV);
    }
    
    //COV1
    can4->cd(1);
    gPad->SetLogy();
    histo_dataCOV1_1MeV_plot->SetTitle("COV1");
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.9);
    histo_dataCOV1_1MeV_plot->SetTitleSize(0.08);
    histo_dataCOV1_1MeV_plot->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV1_1MeV_plot->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV1_1MeV_plot->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV1_1MeV_plot->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV1_1MeV_plot->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV1_1MeV_plot->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV1_1MeV_plot->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1_1MeV_plot->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV1_1MeV_plot->Draw();
    histo_dataCOV1_1MeV_plot->GetXaxis()->SetRangeUser(0, 0.1);
    Hits_to_cps(histo_simu_muCOV1_1MeV, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, NBins_COV1_1MeV);
    Hits_to_cps(histo_simu_nCOV1_1MeV, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,NBins_COV1_1MeV);
    Hits_to_cps(histo_simu_gammaCOV1_1MeV, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,NBins_COV1_1MeV);
    Hits_to_cps(histo_simu_AmCOV1_1MeV, 1e3, MASS_Ge, NumEvents[3], fit_Am.Parameter(3),NBins_COV1_1MeV);
    histo_simu_gammaCOV1_1MeV->SetLineColor(kRed);
    histo_simu_gammaCOV1_1MeV->Add(histo_simu_muCOV1_1MeV);
    histo_simu_gammaCOV1_1MeV->Add(histo_simu_nCOV1_1MeV);
    histo_simu_gammaCOV1_1MeV->Add(histo_simu_AmCOV1_1MeV);
    histo_simu_gammaCOV1_1MeV->Draw("SAME");

    //COV2
    can4->cd(2);
    gPad->SetLogy();
    histo_dataCOV2_1MeV_plot->SetTitle("COV2");
    gPad->SetRightMargin(0.02);
    histo_dataCOV2_1MeV_plot->SetTitleSize(0.08);
    histo_dataCOV2_1MeV_plot->GetXaxis()->SetTitleSize(0.06);
    histo_dataCOV2_1MeV_plot->GetXaxis()->SetLabelSize(0.055);
    histo_dataCOV2_1MeV_plot->GetXaxis()->SetTitleOffset(0.75);
    histo_dataCOV2_1MeV_plot->GetYaxis()->SetTitleSize(0.06);
    histo_dataCOV2_1MeV_plot->GetYaxis()->SetLabelSize(0.06);
    histo_dataCOV2_1MeV_plot->GetYaxis()->SetTitleOffset(0.6);
    histo_dataCOV2_1MeV_plot->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2_1MeV_plot->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2_1MeV_plot->Draw();
    histo_dataCOV2_1MeV_plot->GetXaxis()->SetRangeUser(0, 0.1);
    Hits_to_cps(histo_simu_muCOV2_1MeV, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, NBins_COV2_1MeV);
    Hits_to_cps(histo_simu_nCOV2_1MeV, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,NBins_COV2_1MeV);
    Hits_to_cps(histo_simu_gammaCOV2_1MeV, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,NBins_COV2_1MeV);
    histo_simu_gammaCOV2_1MeV->SetLineColor(kRed);
    histo_simu_gammaCOV2_1MeV->Add(histo_simu_muCOV2_1MeV);
    histo_simu_gammaCOV2_1MeV->Add(histo_simu_nCOV2_1MeV);
    histo_simu_gammaCOV2_1MeV->Draw("SAME");

    //COV3
    can4->cd(3);
    gPad->SetLogy();
    histo_dataNTD_1MeV_plot->SetTitle("NTD");
    gPad->SetRightMargin(0.02);
    histo_dataNTD_1MeV_plot->SetTitleSize(0.08);
    histo_dataNTD_1MeV_plot->GetXaxis()->SetTitleSize(0.06);
    histo_dataNTD_1MeV_plot->GetXaxis()->SetLabelSize(0.055);
    histo_dataNTD_1MeV_plot->GetXaxis()->SetTitleOffset(0.75);
    histo_dataNTD_1MeV_plot->GetYaxis()->SetTitleSize(0.06);
    histo_dataNTD_1MeV_plot->GetYaxis()->SetLabelSize(0.06);
    histo_dataNTD_1MeV_plot->GetYaxis()->SetTitleOffset(0.6);
    histo_dataNTD_1MeV_plot->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataNTD_1MeV_plot->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataNTD_1MeV_plot->Draw();
    histo_dataNTD_1MeV_plot->GetXaxis()->SetRangeUser(0, 0.1);
    histo_dataNTD_1MeV_plot->GetYaxis()->SetRangeUser(1e-5, 1);
    Hits_to_cps(histo_simu_muNTD_1MeV, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, NBins_NTD_1MeV);
    Hits_to_cps(histo_simu_nNTD_1MeV, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,NBins_NTD_1MeV);
    Hits_to_cps(histo_simu_gammaNTD_1MeV, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,NBins_NTD_1MeV);
    TH1D* histo_simu_sumNTD_1MeV =(TH1D*) histo_simu_gammaNTD_1MeV->Clone();
    histo_simu_sumNTD_1MeV->SetLineColor(kRed);
    histo_simu_sumNTD_1MeV->Add(histo_simu_muNTD_1MeV);
    histo_simu_sumNTD_1MeV->Add(histo_simu_nNTD_1MeV);
    histo_simu_sumNTD_1MeV->Draw("SAME");
    
    // Plot all fitted contributions
    //plot_with_res(NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron);

    TFile* fileout = new TFile(NameFileout, "RECREATE");
    plot_with_res_pc(NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron, fit_Am, Syst);

    fileout->Close();
    
    Save_Fitted(NameFileout, NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron, fit_Am, Syst);

//    fileData->Close();
//    filein_mu->Close();
//    filein_Am->Close();
//    filein_n->Close();
//    filein_gamma->Close();
}
                 
void save_CPS_histos(TFile *filein, TFile *fileout, TString Crystal, TString hName, double Crystal_Mass, double thresh_crys, double thresh_crys_fixed, double NumEv, double Rate);

void plot_with_res(vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am){
    
    TH1D* histo_dataCOV1_20MeV_2 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV1_20MeV_2->SetTitle("histo_dataCOV1_20MeV_2");
    TH1D* histo_dataCOV2_20MeV_2 = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV2_20MeV_2->SetTitle("histo_dataCOV2_20MeV_2");
    TH1D* histo_dataNTD_20MeV_2 = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-20MeV");
    histo_dataNTD_20MeV_2->SetTitle("histo_dataNTD_20MeV_2");

    TH1D*  histo_simu_muCOV1_20MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D*  histo_simu_muCOV2_20MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D*  histo_simu_muNTD_20MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    TH1D* histo_simu_gammaCOV1_20MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_gammaCOV2_20MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_gammaNTD_20MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    TH1D* histo_simu_nCOV1_20MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_nCOV2_20MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_nNTD_20MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    TH1D* histo_simu_AmCOV1_20MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_AmCOV2_20MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_AmNTD_20MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    histo_simu_muCOV1_20MeV_2=GoodFormat(histo_simu_muCOV1_20MeV_2, "histo_simu_muCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_muCOV2_20MeV_2=GoodFormat(histo_simu_muCOV2_20MeV_2, "histo_simu_muCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_muNTD_20MeV_2=GoodFormat(histo_simu_muNTD_20MeV_2, "histo_simu_muNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20);

    histo_simu_gammaCOV1_20MeV_2=GoodFormat(histo_simu_gammaCOV1_20MeV_2, "histo_simu_gammaCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_gammaCOV2_20MeV_2=GoodFormat(histo_simu_gammaCOV2_20MeV_2, "histo_simu_gammaCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_gammaNTD_20MeV_2=GoodFormat(histo_simu_gammaNTD_20MeV_2, "histo_simu_gammaNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20);

    histo_simu_nCOV1_20MeV_2=GoodFormat(histo_simu_nCOV1_20MeV_2, "histo_simu_nCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_nCOV2_20MeV_2=GoodFormat(histo_simu_nCOV2_20MeV_2, "histo_simu_nCOV2_10MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_nNTD_20MeV_2=GoodFormat(histo_simu_nNTD_20MeV_2, "histo_simu_nNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20, true);

    histo_simu_AmCOV1_20MeV_2=GoodFormat(histo_simu_AmCOV1_20MeV_2, "histo_simu_AmCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_AmCOV2_20MeV_2=GoodFormat(histo_simu_AmCOV2_20MeV_2, "histo_simu_AmCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_AmNTD_20MeV_2=GoodFormat(histo_simu_AmNTD_20MeV_2, "histo_simu_AmNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20);
    
    TCanvas* cCOV1 =  new TCanvas("Bot Ge (COV1)", "Bot Ge (COV1)", 1400, 600);
    cCOV1->Divide(2,1);
    cCOV1->cd(1);
    gPad->SetLogy();
    histo_dataCOV1_20MeV_2->SetTitle("Bot Ge (COV1) - 0-20MeV");
    gPad->SetRightMargin(0.02);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1_20MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV1_20MeV_2->Draw();
    histo_dataCOV1_20MeV_2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetRangeUser(1e-8, 1);
    Hits_to_cps(histo_simu_muCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV1_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_nCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV1_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_gammaCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV1_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_AmCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_Am.Parameter(3),histo_dataCOV1_20MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV1_20MeV_2 =(TH1D*) histo_simu_gammaCOV1_20MeV_2->Clone();
    histo_simu_sumCOV1_20MeV_2->Add(histo_simu_muCOV1_20MeV_2);
    histo_simu_sumCOV1_20MeV_2->Add(histo_simu_nCOV1_20MeV_2);
    histo_simu_sumCOV1_20MeV_2->Add(histo_simu_AmCOV1_20MeV_2);
    
    histo_dataCOV1_20MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV1_20MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV1_20MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    auto histo_residualsCOV1_20MeV= new TRatioPlot(histo_dataCOV1_20MeV_2, histo_simu_sumCOV1_20MeV_2, "diffsig");
    histo_residualsCOV1_20MeV->Draw("confit");
    histo_residualsCOV1_20MeV->GetLowerRefYaxis()->SetRangeUser(-10, 10);
    histo_residualsCOV1_20MeV->GetLowerRefYaxis()->SetLabelSize(0.03);
    histo_residualsCOV1_20MeV->GetLowerRefXaxis()->SetLabelSize(0.03);
    histo_residualsCOV1_20MeV->GetLowerRefXaxis()->SetTitleSize(0.03);
    histo_residualsCOV1_20MeV->GetLowerRefYaxis()->SetTitleOffset(0);
    histo_residualsCOV1_20MeV->GetLowerRefYaxis()->SetTitle("Residuals [#sigma_{Data}]");
    histo_residualsCOV1_20MeV->GetLowerRefYaxis()->SetTitleSize(0.03);
    histo_residualsCOV1_20MeV->GetLowerRefYaxis()->SetTitleOffset(1.7);
    histo_residualsCOV1_20MeV->GetUpperPad()->SetRightMargin(0.02);
    histo_residualsCOV1_20MeV->GetLowerRefGraph()->SetLineColor(kRed);
    histo_residualsCOV1_20MeV->SetConfidenceIntervalColors();
    cCOV1->Update();
    
    gPad = histo_residualsCOV1_20MeV->GetUpperPad();
    histo_simu_muCOV1_20MeV_2->SetLineColor(kGreen);
    histo_simu_muCOV1_20MeV_2->Draw("SAME");
    histo_simu_nCOV1_20MeV_2->SetLineColor(kMagenta);
    histo_simu_nCOV1_20MeV_2->Draw("SAME");
    histo_simu_gammaCOV1_20MeV_2->SetLineColor(kBlack);
    histo_simu_gammaCOV1_20MeV_2->Draw("SAME");
    histo_simu_AmCOV1_20MeV_2->SetLineColor(kBlue);
    histo_simu_AmCOV1_20MeV_2->Draw("SAME");
    histo_simu_sumCOV1_20MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV1_20MeV_2->Draw("SAME");
    
    auto legendCOV1 = new TLegend(0.2,0.75,0.98,0.9);
    legendCOV1->AddEntry(histo_dataCOV1_20MeV_2,"Data bckg - no shielding");
    legendCOV1->AddEntry(histo_simu_muCOV1_20MeV_2,Form("Simu #mu - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_muon.Parameter(0), fit_muon.ParError(0)));
    legendCOV1->AddEntry(histo_simu_gammaCOV1_20MeV_2, Form("Simu #gamma - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_gamma.Parameter(1), fit_gamma.ParError(1)));
    legendCOV1->AddEntry(histo_simu_nCOV1_20MeV_2, Form("Simu n - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_neutron.Parameter(2), fit_neutron.ParError(2)));
    legendCOV1->AddEntry(histo_simu_AmCOV1_20MeV_2, Form("Simu Am - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_Am.Parameter(3), fit_Am.ParError(3)));
    legendCOV1->AddEntry(histo_simu_sumCOV1_20MeV_2, "Simu #mu + #gamma + n");
    legendCOV1->Draw();
    
//    TH1D* histo_residuals = (TH1D*) histo_dataCOV1_20MeV_2->Clone();
//    double Error_i, Diff_i;
//    for (int i = 1; i<histo_residuals->GetNbinsX(); i++){
//        Error_i = histo_dataCOV1_20MeV_2->GetBinError(i);
//        Diff_i = histo_dataCOV1_20MeV_2->GetBinContent(i)-histo_simu_sumCOV1_20MeV_2->GetBinContent(i);
//        cout<<Diff_i<<endl;
//        cout<<Diff_i/Error_i<<endl;
//        histo_residuals->SetBinContent(i, Diff_i/Error_i);
//        histo_residuals->SetBinError(i, 0);
//    }
//
//    cCOV1->cd(3);
//    gPad->SetPad(0.02, 0.05, 0.5, 0.3);
//    gPad->SetGridy();
//    gPad->SetRightMargin(0.02);
//    histo_residuals->SetTitle("Reduced residuals");
//    histo_residuals->GetYaxis()->SetRangeUser(-20, 20);
//    histo_residuals->GetYaxis()->SetTitle("Residuals [#sigma_{Data}]");
//    histo_residuals->Draw();
    
    TCanvas* cCOV2 =  new TCanvas("Top Ge (COV2)", "Top Ge (COV2)", 1400, 600);
    cCOV2->Divide(2,1);
    cCOV2->cd(1);
    gPad->SetLogy();
    histo_dataCOV2_20MeV_2->SetTitle("Top Ge (COV2) - 0-20MeV");
    gPad->SetRightMargin(0.02);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2_20MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2_20MeV_2->Draw();
    histo_dataCOV2_20MeV_2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetRangeUser(1e-8, 1);
    Hits_to_cps(histo_simu_muCOV2_20MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV2_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_nCOV2_20MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV2_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_gammaCOV2_20MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV2_20MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV2_20MeV_2 =(TH1D*) histo_simu_gammaCOV2_20MeV_2->Clone();
    histo_simu_sumCOV2_20MeV_2->Add(histo_simu_muCOV2_20MeV_2);
    histo_simu_sumCOV2_20MeV_2->Add(histo_simu_nCOV2_20MeV_2);
    
    histo_dataCOV2_20MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV2_20MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV2_20MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    auto histo_residualsCOV2_20MeV= new TRatioPlot(histo_dataCOV2_20MeV_2, histo_simu_sumCOV2_20MeV_2, "diffsig");
    histo_residualsCOV2_20MeV->Draw("confit");
    histo_residualsCOV2_20MeV->GetLowerRefYaxis()->SetRangeUser(-10, 10);
    histo_residualsCOV2_20MeV->GetLowerRefYaxis()->SetLabelSize(0.03);
    histo_residualsCOV2_20MeV->GetLowerRefXaxis()->SetLabelSize(0.03);
    histo_residualsCOV2_20MeV->GetLowerRefXaxis()->SetTitleSize(0.03);
    histo_residualsCOV2_20MeV->GetLowerRefYaxis()->SetTitleOffset(0);
    histo_residualsCOV2_20MeV->GetLowerRefYaxis()->SetTitle("Residuals [#sigma_{Data}]");
    histo_residualsCOV2_20MeV->GetLowerRefYaxis()->SetTitleSize(0.03);
    histo_residualsCOV2_20MeV->GetLowerRefYaxis()->SetTitleOffset(1.7);
    histo_residualsCOV2_20MeV->GetUpperPad()->SetRightMargin(0.02);
    histo_residualsCOV2_20MeV->GetLowerRefGraph()->SetLineColor(kRed);
    histo_residualsCOV2_20MeV->SetConfidenceIntervalColors();
    cCOV2->Update();
    
    gPad = histo_residualsCOV2_20MeV->GetUpperPad();
    histo_simu_muCOV2_20MeV_2->SetLineColor(kGreen);
    histo_simu_muCOV2_20MeV_2->Draw("SAME");
    histo_simu_nCOV2_20MeV_2->SetLineColor(kMagenta);
    histo_simu_nCOV2_20MeV_2->Draw("SAME");
    histo_simu_gammaCOV2_20MeV_2->SetLineColor(kBlack);
    histo_simu_gammaCOV2_20MeV_2->Draw("SAME");
    histo_simu_sumCOV2_20MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV2_20MeV_2->Draw("SAME");
    
    auto legendCOV2 = new TLegend(0.2,0.75,0.98,0.9);
    legendCOV2->AddEntry(histo_dataCOV2_20MeV_2,"Data bckg - no shielding");
    legendCOV2->AddEntry(histo_simu_muCOV2_20MeV_2,Form("Simu #mu - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_muon.Parameter(0), fit_muon.ParError(0)));
    legendCOV2->AddEntry(histo_simu_gammaCOV2_20MeV_2, Form("Simu #gamma - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_gamma.Parameter(1), fit_gamma.ParError(1)));
    legendCOV2->AddEntry(histo_simu_nCOV2_20MeV_2, Form("Simu n - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_neutron.Parameter(2), fit_neutron.ParError(2)));
    legendCOV2->AddEntry(histo_simu_sumCOV2_20MeV_2, "Simu #mu + #gamma + n");
    legendCOV2->Draw();
    
    TCanvas* cNTD =  new TCanvas("LWO", "LWO", 1400, 600);
    cNTD->Divide(2,1);
    cNTD->cd(1);
    gPad->SetLogy();
    histo_dataNTD_20MeV_2->SetTitle("Mid LWO (NTD) - 0-20MeV");
    gPad->SetRightMargin(0.02);
    histo_dataNTD_20MeV_2->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataNTD_20MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataNTD_20MeV_2->Draw();
    histo_dataNTD_20MeV_2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataNTD_20MeV_2->GetYaxis()->SetRangeUser(1e-8, 1);
    Hits_to_cps(histo_simu_muNTD_20MeV_2, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_nNTD_20MeV_2, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataNTD_20MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_gammaNTD_20MeV_2, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataNTD_20MeV_2->GetNbinsX());
    TH1D* histo_simu_sumNTD_20MeV_2 =(TH1D*) histo_simu_gammaNTD_20MeV_2->Clone();
    histo_simu_sumNTD_20MeV_2->Add(histo_simu_muNTD_20MeV_2);
    histo_simu_sumNTD_20MeV_2->Add(histo_simu_nNTD_20MeV_2);
    
    histo_dataNTD_20MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataNTD_20MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataNTD_20MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataNTD_20MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataNTD_20MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    auto histo_residualsNTD_20MeV= new TRatioPlot(histo_dataNTD_20MeV_2, histo_simu_sumNTD_20MeV_2, "diffsig");
    histo_residualsNTD_20MeV->Draw("confit");
    histo_residualsNTD_20MeV->GetLowerRefYaxis()->SetRangeUser(-10, 10);
    histo_residualsNTD_20MeV->GetLowerRefYaxis()->SetLabelSize(0.03);
    histo_residualsNTD_20MeV->GetLowerRefXaxis()->SetLabelSize(0.03);
    histo_residualsNTD_20MeV->GetLowerRefXaxis()->SetTitleSize(0.03);
    histo_residualsNTD_20MeV->GetLowerRefYaxis()->SetTitleOffset(0);
    histo_residualsNTD_20MeV->GetLowerRefYaxis()->SetTitle("Residuals [#sigma_{Data}]");
    histo_residualsNTD_20MeV->GetLowerRefYaxis()->SetTitleSize(0.03);
    histo_residualsNTD_20MeV->GetLowerRefYaxis()->SetTitleOffset(1.7);
    histo_residualsNTD_20MeV->GetUpperPad()->SetRightMargin(0.02);
    histo_residualsNTD_20MeV->GetLowerRefGraph()->SetLineColor(kRed);
    histo_residualsNTD_20MeV->SetConfidenceIntervalColors();
    cNTD->Update();
    
    gPad = histo_residualsNTD_20MeV->GetUpperPad();
    histo_simu_muNTD_20MeV_2->SetLineColor(kGreen);
    histo_simu_muNTD_20MeV_2->Draw("SAME");
    histo_simu_nNTD_20MeV_2->SetLineColor(kMagenta);
    histo_simu_nNTD_20MeV_2->Draw("SAME");
    histo_simu_gammaNTD_20MeV_2->SetLineColor(kBlack);
    histo_simu_gammaNTD_20MeV_2->Draw("SAME");
    histo_simu_sumNTD_20MeV_2->SetLineColor(kRed);
    histo_simu_sumNTD_20MeV_2->Draw("SAME");
    
    auto legendNTD = new TLegend(0.2,0.75,0.98,0.9);
    legendNTD->AddEntry(histo_dataNTD_20MeV_2,"Data bckg - no shielding");
    legendNTD->AddEntry(histo_simu_muNTD_20MeV_2,Form("Simu #mu - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_muon.Parameter(0), fit_muon.ParError(0)));
    legendNTD->AddEntry(histo_simu_gammaNTD_20MeV_2, Form("Simu #gamma - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_gamma.Parameter(1), fit_gamma.ParError(1)));
    legendNTD->AddEntry(histo_simu_nNTD_20MeV_2, Form("Simu n - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_neutron.Parameter(2), fit_neutron.ParError(2)));
    legendNTD->AddEntry(histo_simu_sumNTD_20MeV_2, "Simu #mu + #gamma + n");
    legendNTD->Draw();
    
    
    TH1D* histo_dataCOV1_10MeV_2 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV1_10MeV_2->SetTitle("histo_dataCOV1_10MeV_2");
    TH1D* histo_dataCOV2_10MeV_2 = (TH1D*) fileData->Get("TopGe/TopGe_Spectrum_0-10MeV");
    histo_dataCOV2_10MeV_2->SetTitle("histo_dataCOV2_10MeV_2");
    TH1D* histo_dataNTD_10MeV_2 = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_10MeV_2->SetTitle("histo_dataNTD_10MeV_2");

    TH1D*  histo_simu_muCOV1_10MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D*  histo_simu_muCOV2_10MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D*  histo_simu_muNTD_10MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_gammaCOV1_10MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_gammaCOV2_10MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_gammaNTD_10MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_nCOV1_10MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_nCOV2_10MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_nNTD_10MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_AmCOV1_10MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_AmCOV2_10MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_AmNTD_10MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");
    
    histo_simu_muCOV1_10MeV_2=GoodFormat(histo_simu_muCOV1_10MeV_2, "histo_simu_muCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_muCOV2_10MeV_2=GoodFormat(histo_simu_muCOV2_10MeV_2, "histo_simu_muCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_muNTD_10MeV_2=GoodFormat(histo_simu_muNTD_10MeV_2, "histo_simu_muNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10);

    histo_simu_gammaCOV1_10MeV_2=GoodFormat(histo_simu_gammaCOV1_10MeV_2, "histo_simu_gammaCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_gammaCOV2_10MeV_2=GoodFormat(histo_simu_gammaCOV2_10MeV_2, "histo_simu_gammaCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_gammaNTD_10MeV_2=GoodFormat(histo_simu_gammaNTD_10MeV_2, "histo_simu_gammaNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10);

    histo_simu_nCOV1_10MeV_2=GoodFormat(histo_simu_nCOV1_10MeV_2, "histo_simu_nCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_nCOV2_10MeV_2=GoodFormat(histo_simu_nCOV2_10MeV_2, "histo_simu_nCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_nNTD_10MeV_2=GoodFormat(histo_simu_nNTD_10MeV_2, "histo_simu_nNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10, true);

    histo_simu_AmCOV1_10MeV_2=GoodFormat(histo_simu_AmCOV1_10MeV_2, "histo_simu_AmCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_AmCOV2_10MeV_2=GoodFormat(histo_simu_AmCOV2_10MeV_2, "histo_simu_AmCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_AmNTD_10MeV_2=GoodFormat(histo_simu_AmNTD_10MeV_2, "histo_simu_AmNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10);
    
    cCOV1->cd(2);
    gPad->SetLogy();
    histo_dataCOV1_10MeV_2->SetTitle("Bot Ge (COV1) - 0-1MeV");
    gPad->SetRightMargin(0.02);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1_10MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV1_10MeV_2->Draw();
    histo_dataCOV1_10MeV_2->GetXaxis()->SetRangeUser(0, 4);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetRangeUser(1e-6, 1);
    Hits_to_cps(histo_simu_muCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV1_10MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_nCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV1_10MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_gammaCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV1_10MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV1_10MeV_2 =(TH1D*) histo_simu_gammaCOV1_10MeV_2->Clone();
    histo_simu_sumCOV1_10MeV_2->Add(histo_simu_muCOV1_10MeV_2);
    histo_simu_sumCOV1_10MeV_2->Add(histo_simu_nCOV1_10MeV_2);
    
    histo_dataCOV1_10MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV1_10MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV1_10MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    auto histo_residualsCOV1_10MeV= new TRatioPlot(histo_dataCOV1_10MeV_2, histo_simu_sumCOV1_10MeV_2, "diffsig");
    histo_residualsCOV1_10MeV->Draw("confit");
    histo_residualsCOV1_10MeV->GetLowerRefYaxis()->SetRangeUser(-10, 10);
    histo_residualsCOV1_10MeV->GetLowerRefYaxis()->SetLabelSize(0.03);
    histo_residualsCOV1_10MeV->GetLowerRefXaxis()->SetLabelSize(0.03);
    histo_residualsCOV1_10MeV->GetLowerRefXaxis()->SetTitleSize(0.03);
    histo_residualsCOV1_10MeV->GetLowerRefYaxis()->SetTitleOffset(0);
    histo_residualsCOV1_10MeV->GetLowerRefYaxis()->SetTitle("Residuals [#sigma_{Data}]");
    histo_residualsCOV1_10MeV->GetLowerRefYaxis()->SetTitleSize(0.03);
    histo_residualsCOV1_10MeV->GetLowerRefYaxis()->SetTitleOffset(1.7);
    histo_residualsCOV1_10MeV->GetUpperPad()->SetRightMargin(0.02);
    histo_residualsCOV1_10MeV->GetLowerRefGraph()->SetLineColor(kRed);
    histo_residualsCOV1_10MeV->SetConfidenceIntervalColors();
    cCOV1->Update();
    
    gPad = histo_residualsCOV1_10MeV->GetUpperPad();
    histo_simu_muCOV1_10MeV_2->SetLineColor(kGreen);
    histo_simu_muCOV1_10MeV_2->Draw("SAME");
    histo_simu_nCOV1_10MeV_2->SetLineColor(kMagenta);
    histo_simu_nCOV1_10MeV_2->Draw("SAME");
    histo_simu_gammaCOV1_10MeV_2->SetLineColor(kBlack);
    histo_simu_gammaCOV1_10MeV_2->Draw("SAME");
    histo_simu_sumCOV1_10MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV1_10MeV_2->Draw("SAME");
    
    cCOV2->cd(2);
    gPad->SetLogy();
    histo_dataCOV2_10MeV_2->SetTitle("Top Ge (COV2) - 0-1MeV");
    gPad->SetRightMargin(0.02);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetTitle("Differential Rate (.s^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2_10MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2_10MeV_2->Draw();
    histo_dataCOV2_10MeV_2->GetXaxis()->SetRangeUser(0, 4);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetRangeUser(1e-6, 1);
    Hits_to_cps(histo_simu_muCOV2_10MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV2_10MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_nCOV2_10MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV2_10MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_gammaCOV2_10MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV2_10MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV2_10MeV_2 =(TH1D*) histo_simu_gammaCOV2_10MeV_2->Clone();
    histo_simu_sumCOV2_10MeV_2->Add(histo_simu_muCOV2_10MeV_2);
    histo_simu_sumCOV2_10MeV_2->Add(histo_simu_nCOV2_10MeV_2);
    
    histo_dataCOV2_10MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV2_10MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV2_10MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    auto histo_residualsCOV2_10MeV= new TRatioPlot(histo_dataCOV2_10MeV_2, histo_simu_sumCOV2_10MeV_2, "diffsig");
    histo_residualsCOV2_10MeV->Draw("confit");
    histo_residualsCOV2_10MeV->GetLowerRefYaxis()->SetRangeUser(-10, 10);
    histo_residualsCOV2_10MeV->GetLowerRefYaxis()->SetLabelSize(0.03);
    histo_residualsCOV2_10MeV->GetLowerRefXaxis()->SetLabelSize(0.03);
    histo_residualsCOV2_10MeV->GetLowerRefXaxis()->SetTitleSize(0.03);
    histo_residualsCOV2_10MeV->GetLowerRefYaxis()->SetTitleOffset(0);
    histo_residualsCOV2_10MeV->GetLowerRefYaxis()->SetTitle("Residuals [#sigma_{Data}]");
    histo_residualsCOV2_10MeV->GetLowerRefYaxis()->SetTitleSize(0.03);
    histo_residualsCOV2_10MeV->GetLowerRefYaxis()->SetTitleOffset(1.7);
    histo_residualsCOV2_10MeV->GetUpperPad()->SetRightMargin(0.02);
    histo_residualsCOV2_10MeV->GetLowerRefGraph()->SetLineColor(kRed);
    histo_residualsCOV2_10MeV->SetConfidenceIntervalColors();
    cCOV2->Update();
    
    gPad = histo_residualsCOV2_10MeV->GetUpperPad();
    histo_simu_muCOV2_10MeV_2->SetLineColor(kGreen);
    histo_simu_muCOV2_10MeV_2->Draw("SAME");
    histo_simu_nCOV2_10MeV_2->SetLineColor(kMagenta);
    histo_simu_nCOV2_10MeV_2->Draw("SAME");
    histo_simu_gammaCOV2_10MeV_2->SetLineColor(kBlack);
    histo_simu_gammaCOV2_10MeV_2->Draw("SAME");
    histo_simu_sumCOV2_10MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV2_10MeV_2->Draw("SAME");
    
    cNTD->cd(2);
    gPad->SetLogy();
    histo_dataNTD_10MeV_2->SetTitle("Mid LWO (NTD) - 0-6MeV");
    gPad->SetRightMargin(0.02);
    histo_dataNTD_10MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataNTD_10MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataNTD_10MeV_2->Draw();
    histo_dataNTD_10MeV_2->GetXaxis()->SetRangeUser(0, 6);
    histo_dataNTD_10MeV_2->GetYaxis()->SetRangeUser(1e-8, 1);
    Hits_to_cps(histo_simu_muNTD_10MeV_2, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_10MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_nNTD_10MeV_2, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataNTD_10MeV_2->GetNbinsX());
    Hits_to_cps(histo_simu_gammaNTD_10MeV_2, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataNTD_10MeV_2->GetNbinsX());
    TH1D* histo_simu_sumNTD_10MeV_2 =(TH1D*) histo_simu_gammaNTD_10MeV_2->Clone();
    histo_simu_sumNTD_10MeV_2->Add(histo_simu_muNTD_10MeV_2);
    histo_simu_sumNTD_10MeV_2->Add(histo_simu_nNTD_10MeV_2);
    
    histo_dataNTD_10MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataNTD_10MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataNTD_10MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataNTD_10MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataNTD_10MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    auto histo_residualsNTD_10MeV= new TRatioPlot(histo_dataNTD_10MeV_2, histo_simu_sumNTD_10MeV_2, "diffsig");
    histo_residualsNTD_10MeV->Draw("confit");
    histo_residualsNTD_10MeV->GetLowerRefYaxis()->SetRangeUser(-10, 10);
    histo_residualsNTD_10MeV->GetLowerRefYaxis()->SetLabelSize(0.03);
    histo_residualsNTD_10MeV->GetLowerRefXaxis()->SetLabelSize(0.03);
    histo_residualsNTD_10MeV->GetLowerRefXaxis()->SetTitleSize(0.03);
    histo_residualsNTD_10MeV->GetLowerRefYaxis()->SetTitleOffset(0);
    histo_residualsNTD_10MeV->GetLowerRefYaxis()->SetTitle("Residuals [#sigma_{Data}]");
    histo_residualsNTD_10MeV->GetLowerRefYaxis()->SetTitleSize(0.03);
    histo_residualsNTD_10MeV->GetLowerRefYaxis()->SetTitleOffset(1.7);
    histo_residualsNTD_10MeV->GetUpperPad()->SetRightMargin(0.02);
    histo_residualsNTD_10MeV->GetLowerRefGraph()->SetLineColor(kRed);
    histo_residualsNTD_10MeV->SetConfidenceIntervalColors();
    cNTD->Update();
    
    gPad = histo_residualsNTD_10MeV->GetUpperPad();
    histo_simu_muNTD_10MeV_2->SetLineColor(kGreen);
    histo_simu_muNTD_10MeV_2->Draw("SAME");
    histo_simu_nNTD_10MeV_2->SetLineColor(kMagenta);
    histo_simu_nNTD_10MeV_2->Draw("SAME");
    histo_simu_gammaNTD_10MeV_2->SetLineColor(kBlack);
    histo_simu_gammaNTD_10MeV_2->Draw("SAME");
    histo_simu_sumNTD_10MeV_2->SetLineColor(kRed);
    histo_simu_sumNTD_10MeV_2->Draw("SAME");
}
    
void plot_with_res_pc(vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am, bool Syst=false){
    
    TH1D* histo_dataCOV1_20MeV_2 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV1_20MeV_2->SetTitle("histo_dataCOV1_20MeV_2");
    TH1D* histo_dataCOV2_20MeV_2 = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-20MeV");
    histo_dataCOV2_20MeV_2->SetTitle("histo_dataCOV2_20MeV_2");
    TH1D* histo_dataNTD_20MeV_2 = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-20MeV");
    histo_dataNTD_20MeV_2->SetTitle("histo_dataNTD_20MeV_2");

    TH2D* V_dataCOV1_20MeV = Stat_Covariance_TH2D(histo_dataCOV1_20MeV_2);
    V_dataCOV1_20MeV->SetTitle("V_dataCOV1_20MeV_2");
    TH2D* V_dataCOV2_20MeV= Stat_Covariance_TH2D(histo_dataCOV2_20MeV_2);
    V_dataCOV2_20MeV->SetTitle("V_dataCOV2_20MeV_2");
    TH2D* V_dataNTD_20MeV = Stat_Covariance_TH2D(histo_dataNTD_20MeV_2);
    V_dataNTD_20MeV->SetTitle("V_dataNTD_20MeV_2");
    
    if (Syst){
        V_dataCOV1_20MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-20MeV");
        V_dataCOV1_20MeV->SetTitle("V_dataCOV1_20MeV_2");
        V_dataCOV2_20MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-20MeV");
        V_dataCOV2_20MeV->SetTitle("V_dataCOV2_20MeV_2");
        V_dataNTD_20MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-20MeV");
        V_dataNTD_20MeV->SetTitle("V_dataNTD_20MeV_2");
        
        SetErrors_from_CovMat(histo_dataCOV1_20MeV_2, V_dataCOV1_20MeV);
        SetErrors_from_CovMat(histo_dataCOV2_20MeV_2, V_dataCOV2_20MeV);
        SetErrors_from_CovMat(histo_dataNTD_20MeV_2, V_dataNTD_20MeV);
    }
    
    TH1D*  histo_simu_muCOV1_20MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D*  histo_simu_muCOV2_20MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D*  histo_simu_muNTD_20MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    TH1D* histo_simu_gammaCOV1_20MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_gammaCOV2_20MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_gammaNTD_20MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    TH1D* histo_simu_nCOV1_20MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_nCOV2_20MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_nNTD_20MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    TH1D* histo_simu_AmCOV1_20MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-100MeV");
    TH1D* histo_simu_AmCOV2_20MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    TH1D* histo_simu_AmNTD_20MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");

    histo_simu_muCOV1_20MeV_2=GoodFormat(histo_simu_muCOV1_20MeV_2, "histo_simu_muCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_muCOV2_20MeV_2=GoodFormat(histo_simu_muCOV2_20MeV_2, "histo_simu_muCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_muNTD_20MeV_2=GoodFormat(histo_simu_muNTD_20MeV_2, "histo_simu_muNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20);

    histo_simu_gammaCOV1_20MeV_2=GoodFormat(histo_simu_gammaCOV1_20MeV_2, "histo_simu_gammaCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_gammaCOV2_20MeV_2=GoodFormat(histo_simu_gammaCOV2_20MeV_2, "histo_simu_gammaCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_gammaNTD_20MeV_2=GoodFormat(histo_simu_gammaNTD_20MeV_2, "histo_simu_gammaNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20);

    histo_simu_nCOV1_20MeV_2=GoodFormat(histo_simu_nCOV1_20MeV_2, "histo_simu_nCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_nCOV2_20MeV_2=GoodFormat(histo_simu_nCOV2_20MeV_2, "histo_simu_nCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_nNTD_20MeV_2=GoodFormat(histo_simu_nNTD_20MeV_2, "histo_simu_nNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20, true);

    histo_simu_AmCOV1_20MeV_2=GoodFormat(histo_simu_AmCOV1_20MeV_2, "histo_simu_AmCOV1_20MeV_2", histo_dataCOV1_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_AmCOV2_20MeV_2=GoodFormat(histo_simu_AmCOV2_20MeV_2, "histo_simu_AmCOV2_20MeV_2", histo_dataCOV2_20MeV_2->GetNbinsX(), 0, 20);
    histo_simu_AmNTD_20MeV_2=GoodFormat(histo_simu_AmNTD_20MeV_2, "histo_simu_AmNTD_20MeV_2", histo_dataNTD_20MeV_2->GetNbinsX(), 0, 20);
    
    TCanvas* cCOV1 =  new TCanvas("Bot Ge (COV1)", "Bot Ge (COV1)", 1400, 600);
    cCOV1->Divide(2,2);
    
    cCOV1->cd(1);
    gPad->SetPad(0., 0.23, 0.5, 1.);
    gPad->SetLogy();
    gPad->SetRightMargin(0.02);
    
    if (!TString(histo_dataCOV1_20MeV_2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataCOV1_20MeV_2->Scale(24*3600.);
    histo_dataCOV1_20MeV_2->SetTitle("Bot Ge (COV1) - 0-20MeV");
    histo_dataCOV1_20MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1_20MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV1_20MeV_2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    Hits_to_dru(histo_simu_muCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV1_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_nCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV1_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV1_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_AmCOV1_20MeV_2, 1e3, MASS_Ge, NumEvents[3], fit_Am.Parameter(3),histo_dataCOV1_20MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV1_20MeV_2 =(TH1D*) histo_simu_gammaCOV1_20MeV_2->Clone();
    histo_simu_sumCOV1_20MeV_2->Add(histo_simu_muCOV1_20MeV_2);
    histo_simu_sumCOV1_20MeV_2->Add(histo_simu_nCOV1_20MeV_2);
    histo_simu_sumCOV1_20MeV_2->Add(histo_simu_AmCOV1_20MeV_2);
    
    histo_dataCOV1_20MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetTitleOffset(1.4);
    histo_dataCOV1_20MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV1_20MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV1_20MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    histo_dataCOV1_20MeV_2->Draw("E0LF2");
    histo_simu_muCOV1_20MeV_2->SetLineColor(kCyan-5);
    histo_simu_muCOV1_20MeV_2->SetFillColorAlpha(kCyan-5, 0.20);
    histo_simu_muCOV1_20MeV_2->Draw("E0LF2SAME");
    histo_simu_nCOV1_20MeV_2->SetLineColor(kOrange-5);
    histo_simu_nCOV1_20MeV_2->SetFillColorAlpha(kOrange-5, 0.20);
    histo_simu_nCOV1_20MeV_2->Draw("E0LF2SAME");
    histo_simu_gammaCOV1_20MeV_2->SetLineColor(kMagenta-5);
    histo_simu_gammaCOV1_20MeV_2->SetFillColorAlpha(kMagenta-5, 0.20);
    histo_simu_gammaCOV1_20MeV_2->Draw("E0LF2SAME");
    histo_simu_AmCOV1_20MeV_2->SetLineColor(kRed-7);
    histo_simu_AmCOV1_20MeV_2->SetFillColorAlpha(kRed-7, 0.2);
    histo_simu_AmCOV1_20MeV_2->Draw("E0LF2SAME");
    histo_simu_sumCOV1_20MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV1_20MeV_2->Draw("E0LF2SAME");
    
    auto legendCOV1 = new TLegend(0.2,0.75,0.98,0.9);
    legendCOV1->AddEntry(histo_dataCOV1_20MeV_2,"Data bckg - no shielding");
    legendCOV1->AddEntry(histo_simu_muCOV1_20MeV_2,Form("Simu #mu - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_muon.Parameter(0), fit_muon.ParError(0)));
    legendCOV1->AddEntry(histo_simu_gammaCOV1_20MeV_2, Form("Simu #gamma - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_gamma.Parameter(1), fit_gamma.ParError(1)));
    legendCOV1->AddEntry(histo_simu_nCOV1_20MeV_2, Form("Simu n - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_neutron.Parameter(2), fit_neutron.ParError(2)));
    legendCOV1->AddEntry(histo_simu_AmCOV1_20MeV_2, Form("Simu Am - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_Am.Parameter(3), fit_Am.ParError(3)));
    legendCOV1->AddEntry(histo_simu_sumCOV1_20MeV_2, "Simu #mu + #gamma + n + Am");
    legendCOV1->Draw();
    
    TH1D* histo_residualsCOV1_20MeV_pc = (TH1D*) histo_dataCOV1_20MeV_2->Clone();
    double Value_i, Diff_i, Error_i;
    
    double chi2=0;
    for (int i = 1; i<histo_residualsCOV1_20MeV_pc->GetNbinsX(); i++){
        Value_i = histo_dataCOV1_20MeV_2->GetBinContent(i);
        Diff_i = histo_dataCOV1_20MeV_2->GetBinContent(i)-histo_simu_sumCOV1_20MeV_2->GetBinContent(i);
        Error_i = sqrt(histo_dataCOV1_20MeV_2->GetBinError(i)*histo_dataCOV1_20MeV_2->GetBinError(i)+histo_simu_sumCOV1_20MeV_2->GetBinError(i)*histo_simu_sumCOV1_20MeV_2->GetBinError(i));
        histo_residualsCOV1_20MeV_pc->SetBinContent(i, Diff_i/Value_i*100);
        histo_residualsCOV1_20MeV_pc->SetBinError(i, Error_i/Value_i*100);
        chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
    }

    TLatex *Latex_COV1_20MeV = new TLatex(12, 1000, Form("Chi2/dof = %0.2f / %d = %0.2f", chi2, histo_residualsCOV1_20MeV_pc->GetNbinsX()-4, 1.*chi2/(histo_residualsCOV1_20MeV_pc->GetNbinsX()-1)));
    Latex_COV1_20MeV->SetTextSize(.030);
    Latex_COV1_20MeV->SetTextColor(1);
    Latex_COV1_20MeV->SetTextAlign(22);
    Latex_COV1_20MeV->Draw();
    
    cCOV1->cd(3);
    gPad->SetPad(0., 0., 0.5, 0.3);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetGridy();
    
    histo_residualsCOV1_20MeV_pc->GetYaxis()->SetTitle("Residuals [%]");
    histo_residualsCOV1_20MeV_pc->SetTitle("");
    histo_residualsCOV1_20MeV_pc->GetYaxis()->SetTitleSize(0.08);
    histo_residualsCOV1_20MeV_pc->GetYaxis()->SetTitleOffset(0.5);
    histo_residualsCOV1_20MeV_pc->GetXaxis()->SetTitleOffset(1);
    histo_residualsCOV1_20MeV_pc->GetYaxis()->SetLabelSize(0.08);
    histo_residualsCOV1_20MeV_pc->GetXaxis()->SetLabelSize(0.08);
    histo_residualsCOV1_20MeV_pc->GetXaxis()->SetTitleSize(0.08);
    histo_residualsCOV1_20MeV_pc->GetYaxis()->SetRangeUser(-100, 100);
    histo_residualsCOV1_20MeV_pc->SetFillColorAlpha(kBlue+1, 0.2);
    histo_residualsCOV1_20MeV_pc->Draw("E3");
    
//    cCOV1->cd(3);
//    gPad->SetPad(0.02, 0.05, 0.5, 0.3);
//    gPad->SetGridy();
//    gPad->SetRightMargin(0.02);
//    histo_residuals->SetTitle("Reduced residuals");
//    histo_residuals->GetYaxis()->SetRangeUser(-20, 20);
//    histo_residuals->GetYaxis()->SetTitle("Residuals [#sigma_{Data}]");
//    histo_residuals->Draw();
    
    TCanvas* cCOV2 =  new TCanvas("Top Ge (COV2)", "Top Ge (COV2)", 1400, 600);
    cCOV2->Divide(2,2);
    cCOV2->cd(1);
    gPad->SetPad(0., 0.23, 0.5, 1.);
    gPad->SetLogy();
    gPad->SetRightMargin(0.02);
    
    if (!TString(histo_dataCOV2_20MeV_2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataCOV2_20MeV_2->Scale(24*3600.);
    histo_dataCOV2_20MeV_2->SetTitle("Top Ge (COV2) - 0-20MeV");
    histo_dataCOV2_20MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2_20MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2_20MeV_2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    Hits_to_dru(histo_simu_muCOV2_20MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV2_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_nCOV2_20MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV2_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaCOV2_20MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV2_20MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV2_20MeV_2 =(TH1D*) histo_simu_gammaCOV2_20MeV_2->Clone();
    histo_simu_sumCOV2_20MeV_2->Add(histo_simu_muCOV2_20MeV_2);
    histo_simu_sumCOV2_20MeV_2->Add(histo_simu_nCOV2_20MeV_2);
    
    histo_dataCOV2_20MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV2_20MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV2_20MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV2_20MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    histo_dataCOV2_20MeV_2->Draw("E0LF2");
    histo_simu_muCOV2_20MeV_2->SetLineColor(kCyan-5);
    histo_simu_muCOV2_20MeV_2->SetFillColorAlpha(kCyan-5, 0.20);
    histo_simu_muCOV2_20MeV_2->Draw("E0LF2SAME");
    histo_simu_nCOV2_20MeV_2->SetLineColor(kOrange-5);
    histo_simu_nCOV2_20MeV_2->SetFillColorAlpha(kOrange-5, 0.20);
    histo_simu_nCOV2_20MeV_2->Draw("E0LF2SAME");
    histo_simu_gammaCOV2_20MeV_2->SetLineColor(kMagenta-5);
    histo_simu_gammaCOV2_20MeV_2->SetFillColorAlpha(kMagenta-5, 0.20);
    histo_simu_gammaCOV2_20MeV_2->Draw("E0LF2SAME");
    histo_simu_AmCOV2_20MeV_2->SetLineColor(kRed-7);
    histo_simu_AmCOV2_20MeV_2->SetFillColorAlpha(kRed-7, 0.2);
    histo_simu_AmCOV2_20MeV_2->Draw("E0LF2SAME");
    histo_simu_sumCOV2_20MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV2_20MeV_2->Draw("E0LF2SAME");
    
    auto legendCOV2 = new TLegend(0.2,0.75,0.98,0.9);
    legendCOV2->AddEntry(histo_dataCOV2_20MeV_2,"Data bckg - no shielding");
    legendCOV2->AddEntry(histo_simu_muCOV2_20MeV_2,Form("Simu #mu - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_muon.Parameter(0), fit_muon.ParError(0)));
    legendCOV2->AddEntry(histo_simu_gammaCOV2_20MeV_2, Form("Simu #gamma - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_gamma.Parameter(1), fit_gamma.ParError(1)));
    legendCOV2->AddEntry(histo_simu_nCOV2_20MeV_2, Form("Simu n - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_neutron.Parameter(2), fit_neutron.ParError(2)));
    legendCOV2->AddEntry(histo_simu_sumCOV2_20MeV_2, "Simu #mu + #gamma + n");
    legendCOV2->Draw();
    
    TH1D* histo_residualsCOV2_20MeV_pc = (TH1D*) histo_dataCOV2_20MeV_2->Clone();
    chi2=0;
    for (int i = 1; i<histo_residualsCOV2_20MeV_pc->GetNbinsX(); i++){
        Value_i = histo_dataCOV2_20MeV_2->GetBinContent(i);
        Diff_i = histo_dataCOV2_20MeV_2->GetBinContent(i)-histo_simu_sumCOV2_20MeV_2->GetBinContent(i);
        Error_i = sqrt(histo_dataCOV2_20MeV_2->GetBinError(i)*histo_dataCOV2_20MeV_2->GetBinError(i)+histo_simu_sumCOV2_20MeV_2->GetBinError(i)*histo_simu_sumCOV2_20MeV_2->GetBinError(i));
        histo_residualsCOV2_20MeV_pc->SetBinContent(i, Diff_i/Value_i*100);
        histo_residualsCOV2_20MeV_pc->SetBinError(i, Error_i/Value_i*100);
        chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
    }
    
    TLatex *Latex_COV2_20MeV = new TLatex(12, 1000, Form("Chi2/dof = %0.2f / %d = %0.2f", chi2, histo_residualsCOV2_20MeV_pc->GetNbinsX()-4, 1.*chi2/(histo_residualsCOV2_20MeV_pc->GetNbinsX()-1)));
    Latex_COV2_20MeV->SetTextSize(.030);
    Latex_COV2_20MeV->SetTextColor(1);
    Latex_COV2_20MeV->SetTextAlign(22);
    Latex_COV2_20MeV->Draw();
    
    cCOV2->cd(3);
    gPad->SetPad(0., 0., 0.5, 0.3);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetGridy();
    
    histo_residualsCOV2_20MeV_pc->GetYaxis()->SetTitle("Residuals [%]");
    histo_residualsCOV2_20MeV_pc->SetTitle("");
    histo_residualsCOV2_20MeV_pc->GetYaxis()->SetTitleSize(0.08);
    histo_residualsCOV2_20MeV_pc->GetYaxis()->SetTitleOffset(0.5);
    histo_residualsCOV2_20MeV_pc->GetXaxis()->SetTitleOffset(1);
    histo_residualsCOV2_20MeV_pc->GetYaxis()->SetLabelSize(0.08);
    histo_residualsCOV2_20MeV_pc->GetXaxis()->SetLabelSize(0.08);
    histo_residualsCOV2_20MeV_pc->GetXaxis()->SetTitleSize(0.08);
    histo_residualsCOV2_20MeV_pc->GetYaxis()->SetRangeUser(-100, 100);
    histo_residualsCOV2_20MeV_pc->SetFillColorAlpha(kBlue+1, 0.2);
    histo_residualsCOV2_20MeV_pc->Draw("E3");
    
    TCanvas* cNTD =  new TCanvas("LWO", "LWO", 1400, 600);
    cNTD->Divide(2,2);
    cNTD->cd(1);
    gPad->SetPad(0., 0.23, 0.5, 1.);
    gPad->SetLogy();
    gPad->SetRightMargin(0.02);
    
    if (!TString(histo_dataNTD_20MeV_2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataNTD_20MeV_2->Scale(24*3600.);
    histo_dataNTD_20MeV_2->SetTitle("Mid LWO (NTD) - 0-20MeV");
    histo_dataNTD_20MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataNTD_20MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataNTD_20MeV_2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataNTD_20MeV_2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    Hits_to_dru(histo_simu_muNTD_20MeV_2, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_nNTD_20MeV_2, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataNTD_20MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaNTD_20MeV_2, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataNTD_20MeV_2->GetNbinsX());
    TH1D* histo_simu_sumNTD_20MeV_2 =(TH1D*) histo_simu_gammaNTD_20MeV_2->Clone();
    histo_simu_sumNTD_20MeV_2->Add(histo_simu_muNTD_20MeV_2);
    histo_simu_sumNTD_20MeV_2->Add(histo_simu_nNTD_20MeV_2);
    
    histo_dataNTD_20MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataNTD_20MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataNTD_20MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataNTD_20MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataNTD_20MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    histo_dataNTD_20MeV_2->Draw("E0LF2");
    histo_simu_muNTD_20MeV_2->SetLineColor(kCyan-5);
    histo_simu_muNTD_20MeV_2->SetFillColorAlpha(kCyan-5, 0.20);
    histo_simu_muNTD_20MeV_2->Draw("E0LF2SAME");
    histo_simu_nNTD_20MeV_2->SetLineColor(kOrange-5);
    histo_simu_nNTD_20MeV_2->SetFillColorAlpha(kOrange-5, 0.20);
    histo_simu_nNTD_20MeV_2->Draw("E0LF2SAME");
    histo_simu_gammaNTD_20MeV_2->SetLineColor(kMagenta-5);
    histo_simu_gammaNTD_20MeV_2->SetFillColorAlpha(kMagenta-5, 0.20);
    histo_simu_gammaNTD_20MeV_2->Draw("E0LF2SAME");
    histo_simu_AmNTD_20MeV_2->SetLineColor(kRed-7);
    histo_simu_AmNTD_20MeV_2->SetFillColorAlpha(kRed-7, 0.2);
    histo_simu_AmNTD_20MeV_2->Draw("E0LF2SAME");
    histo_simu_sumNTD_20MeV_2->SetLineColor(kRed);
    histo_simu_sumNTD_20MeV_2->Draw("E0LF2SAME");
    
    auto legendNTD = new TLegend(0.2,0.75,0.98,0.9);
    legendNTD->AddEntry(histo_dataNTD_20MeV_2,"Data bckg - no shielding");
    legendNTD->AddEntry(histo_simu_muNTD_20MeV_2,Form("Simu #mu - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_muon.Parameter(0), fit_muon.ParError(0)));
    legendNTD->AddEntry(histo_simu_gammaNTD_20MeV_2, Form("Simu #gamma - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_gamma.Parameter(1), fit_gamma.ParError(1)));
    legendNTD->AddEntry(histo_simu_nNTD_20MeV_2, Form("Simu n - no shielding - #Phi = %0.4f +- %0.4f ev/(s cm^{2})", fit_neutron.Parameter(2), fit_neutron.ParError(2)));
    legendNTD->AddEntry(histo_simu_sumNTD_20MeV_2, "Simu #mu + #gamma + n");
    legendNTD->Draw();
    
    TH1D* histo_residualsNTD_20MeV_pc = (TH1D*) histo_dataNTD_20MeV_2->Clone();
    chi2 =0;
    for (int i = 1; i<histo_residualsNTD_20MeV_pc->GetNbinsX(); i++){
        Value_i = histo_dataNTD_20MeV_2->GetBinContent(i);
        Diff_i = histo_dataNTD_20MeV_2->GetBinContent(i)-histo_simu_sumNTD_20MeV_2->GetBinContent(i);
        Error_i = sqrt(histo_dataNTD_20MeV_2->GetBinError(i)*histo_dataNTD_20MeV_2->GetBinError(i)+histo_simu_sumNTD_20MeV_2->GetBinError(i)*histo_simu_sumNTD_20MeV_2->GetBinError(i));
        histo_residualsNTD_20MeV_pc->SetBinContent(i, Diff_i/Value_i*100);
        histo_residualsNTD_20MeV_pc->SetBinError(i, Error_i/Value_i*100);
        chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
    }

    TLatex *Latex_NTD_20MeV = new TLatex(12, 1000, Form("Chi2/dof = %0.2f / %d = %0.2f", chi2, histo_residualsNTD_20MeV_pc->GetNbinsX()-4, 1.*chi2/(histo_residualsNTD_20MeV_pc->GetNbinsX()-1)));
    Latex_NTD_20MeV->SetTextSize(.030);
    Latex_NTD_20MeV->SetTextColor(1);
    Latex_NTD_20MeV->SetTextAlign(22);
    Latex_NTD_20MeV->Draw();
    
    cNTD->cd(3);
    gPad->SetPad(0., 0., 0.5, 0.3);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetGridy();
    
    histo_residualsNTD_20MeV_pc->GetYaxis()->SetTitle("Residuals [%]");
    histo_residualsNTD_20MeV_pc->SetTitle("");
    histo_residualsNTD_20MeV_pc->GetYaxis()->SetTitleSize(0.08);
    histo_residualsNTD_20MeV_pc->GetYaxis()->SetTitleOffset(0.5);
    histo_residualsNTD_20MeV_pc->GetXaxis()->SetTitleOffset(1);
    histo_residualsNTD_20MeV_pc->GetYaxis()->SetLabelSize(0.08);
    histo_residualsNTD_20MeV_pc->GetXaxis()->SetLabelSize(0.08);
    histo_residualsNTD_20MeV_pc->GetXaxis()->SetTitleSize(0.08);
    histo_residualsNTD_20MeV_pc->GetYaxis()->SetRangeUser(-100, 100);
    histo_residualsNTD_20MeV_pc->SetFillColorAlpha(kBlue+1, 0.2);
    histo_residualsNTD_20MeV_pc->Draw("E3");
    
    
    TH1D* histo_dataCOV1_10MeV_2 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV1_10MeV_2->SetTitle("histo_dataCOV1_10MeV_2");
    TH1D* histo_dataCOV2_10MeV_2 = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-10MeV");
    histo_dataCOV2_10MeV_2->SetTitle("histo_dataCOV2_10MeV_2");
    TH1D* histo_dataNTD_10MeV_2 = (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-10MeV");
    histo_dataNTD_10MeV_2->SetTitle("histo_dataNTD_10MeV_2");

    TH2D* V_dataCOV1_10MeV = Stat_Covariance_TH2D(histo_dataCOV1_10MeV_2);
    V_dataCOV1_10MeV->SetTitle("V_dataCOV1_10MeV_2");
    TH2D* V_dataCOV2_10MeV= Stat_Covariance_TH2D(histo_dataCOV2_10MeV_2);
    V_dataCOV2_10MeV->SetTitle("V_dataCOV2_10MeV_2");
    TH2D* V_dataNTD_10MeV = Stat_Covariance_TH2D(histo_dataNTD_10MeV_2);
    V_dataNTD_10MeV->SetTitle("V_dataNTD_10MeV_2");
    
    if (Syst){
        V_dataCOV1_10MeV = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-10MeV");
        V_dataCOV1_10MeV->SetTitle("V_dataCOV1_10MeV_2");
        V_dataCOV2_10MeV = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-10MeV");
        V_dataCOV2_10MeV->SetTitle("V_dataCOV2_10MeV_2");
        V_dataNTD_10MeV = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-10MeV");
        V_dataNTD_10MeV->SetTitle("V_dataNTD_10MeV_2");
        
        SetErrors_from_CovMat(histo_dataCOV1_10MeV_2, V_dataCOV1_10MeV);
        SetErrors_from_CovMat(histo_dataCOV2_10MeV_2, V_dataCOV2_10MeV);
        SetErrors_from_CovMat(histo_dataNTD_10MeV_2, V_dataNTD_10MeV);
    }
    

    TH1D*  histo_simu_muCOV1_10MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D*  histo_simu_muCOV2_10MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D*  histo_simu_muNTD_10MeV_2 = (TH1D*) filein_mu->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_gammaCOV1_10MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_gammaCOV2_10MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_gammaNTD_10MeV_2 = (TH1D*) filein_gamma->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_nCOV1_10MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_nCOV2_10MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_nNTD_10MeV_2 = (TH1D*) filein_n->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");

    TH1D* histo_simu_AmCOV1_10MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    TH1D* histo_simu_AmCOV2_10MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    TH1D* histo_simu_AmNTD_10MeV_2 = (TH1D*) filein_Am->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");
    
    histo_simu_muCOV1_10MeV_2=GoodFormat(histo_simu_muCOV1_10MeV_2, "histo_simu_muCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_muCOV2_10MeV_2=GoodFormat(histo_simu_muCOV2_10MeV_2, "histo_simu_muCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_muNTD_10MeV_2=GoodFormat(histo_simu_muNTD_10MeV_2, "histo_simu_muNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10);

    histo_simu_gammaCOV1_10MeV_2=GoodFormat(histo_simu_gammaCOV1_10MeV_2, "histo_simu_gammaCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_gammaCOV2_10MeV_2=GoodFormat(histo_simu_gammaCOV2_10MeV_2, "histo_simu_gammaCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_gammaNTD_10MeV_2=GoodFormat(histo_simu_gammaNTD_10MeV_2, "histo_simu_gammaNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10);

    histo_simu_nCOV1_10MeV_2=GoodFormat(histo_simu_nCOV1_10MeV_2, "histo_simu_nCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_nCOV2_10MeV_2=GoodFormat(histo_simu_nCOV2_10MeV_2, "histo_simu_nCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_nNTD_10MeV_2=GoodFormat(histo_simu_nNTD_10MeV_2, "histo_simu_nNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10, true);

    histo_simu_AmCOV1_10MeV_2=GoodFormat(histo_simu_AmCOV1_10MeV_2, "histo_simu_AmCOV1_10MeV_2", histo_dataCOV1_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_AmCOV2_10MeV_2=GoodFormat(histo_simu_AmCOV2_10MeV_2, "histo_simu_AmCOV2_10MeV_2", histo_dataCOV2_10MeV_2->GetNbinsX(), 0, 10);
    histo_simu_AmNTD_10MeV_2=GoodFormat(histo_simu_AmNTD_10MeV_2, "histo_simu_AmNTD_10MeV_2", histo_dataNTD_10MeV_2->GetNbinsX(), 0, 10);
    
    cCOV1->cd(2);
    gPad->SetPad(0.5, 0.23, 1., 1.);
    gPad->SetLogy();
    gPad->SetRightMargin(0.02);
    
    if (!TString(histo_dataCOV1_10MeV_2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataCOV1_10MeV_2->Scale(24*3600.);
    histo_dataCOV1_10MeV_2->SetTitle("Bot Ge (COV1) - 0-4MeV");
    histo_dataCOV1_10MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1_10MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV1_10MeV_2->GetXaxis()->SetRangeUser(0, 4);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    Hits_to_dru(histo_simu_muCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV1_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_nCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV1_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV1_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_AmCOV1_10MeV_2, 1e3, MASS_Ge, NumEvents[3], fit_Am.Parameter(3),histo_dataCOV1_10MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV1_10MeV_2 =(TH1D*) histo_simu_gammaCOV1_10MeV_2->Clone();
    histo_simu_sumCOV1_10MeV_2->Add(histo_simu_muCOV1_10MeV_2);
    histo_simu_sumCOV1_10MeV_2->Add(histo_simu_nCOV1_10MeV_2);
    histo_simu_sumCOV1_10MeV_2->Add(histo_simu_AmCOV1_10MeV_2);
    
    histo_dataCOV1_10MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV1_10MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV1_10MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV1_10MeV_2->GetXaxis()->SetTitleSize(0.03);
    
    histo_dataCOV1_10MeV_2->Draw("E0LF2");
    histo_simu_muCOV1_10MeV_2->SetLineColor(kCyan-5);
    histo_simu_muCOV1_10MeV_2->SetFillColorAlpha(kCyan-5, 0.20);
    histo_simu_muCOV1_10MeV_2->Draw("E0LF2SAME");
    histo_simu_nCOV1_10MeV_2->SetLineColor(kOrange-5);
    histo_simu_nCOV1_10MeV_2->SetFillColorAlpha(kOrange-5, 0.20);
    histo_simu_nCOV1_10MeV_2->Draw("E0LF2SAME");
    histo_simu_gammaCOV1_10MeV_2->SetLineColor(kMagenta-5);
    histo_simu_gammaCOV1_10MeV_2->SetFillColorAlpha(kMagenta-5, 0.20);
    histo_simu_gammaCOV1_10MeV_2->Draw("E0LF2SAME");
    histo_simu_AmCOV1_10MeV_2->SetLineColor(kRed-7);
    histo_simu_AmCOV1_10MeV_2->SetFillColorAlpha(kRed-7, 0.2);
    histo_simu_AmCOV1_10MeV_2->Draw("E0LF2SAME");
    histo_simu_sumCOV1_10MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV1_10MeV_2->Draw("E0LF2SAME");
    
    TH1D* histo_residualsCOV1_10MeV_pc = (TH1D*) histo_dataCOV1_10MeV_2->Clone();
    chi2=0;
    for (int i = 1; i<histo_residualsCOV1_10MeV_pc->GetNbinsX(); i++){
        Value_i = histo_dataCOV1_10MeV_2->GetBinContent(i);
        Diff_i = histo_dataCOV1_10MeV_2->GetBinContent(i)-histo_simu_sumCOV1_10MeV_2->GetBinContent(i);
        Error_i = sqrt(histo_dataCOV1_10MeV_2->GetBinError(i)*histo_dataCOV1_10MeV_2->GetBinError(i)+histo_simu_sumCOV1_10MeV_2->GetBinError(i)*histo_simu_sumCOV1_10MeV_2->GetBinError(i));
        histo_residualsCOV1_10MeV_pc->SetBinContent(i, Diff_i/Value_i*100);
        histo_residualsCOV1_10MeV_pc->SetBinError(i, Error_i/Value_i*100);
        chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
    }

    TLatex *Latex_COV1_10MeV = new TLatex(3, 1000, Form("Chi2/dof = %0.2f / %d = %0.2f", chi2, histo_residualsCOV1_10MeV_pc->GetNbinsX()-4, 1.*chi2/(histo_residualsCOV1_10MeV_pc->GetNbinsX()-1)));
    Latex_COV1_10MeV->SetTextSize(.030);
    Latex_COV1_10MeV->SetTextColor(1);
    Latex_COV1_10MeV->SetTextAlign(22);
    Latex_COV1_10MeV->Draw();
    
    cCOV1->cd(4);
    gPad->SetPad(0.5, 0., 1., 0.3);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetGridy();
    
    histo_residualsCOV1_10MeV_pc->GetYaxis()->SetTitle("Residuals [%]");
    histo_residualsCOV1_10MeV_pc->SetTitle("");
    histo_residualsCOV1_10MeV_pc->GetYaxis()->SetTitleSize(0.08);
    histo_residualsCOV1_10MeV_pc->GetYaxis()->SetTitleOffset(0.5);
    histo_residualsCOV1_10MeV_pc->GetXaxis()->SetTitleOffset(1);
    histo_residualsCOV1_10MeV_pc->GetYaxis()->SetLabelSize(0.08);
    histo_residualsCOV1_10MeV_pc->GetXaxis()->SetLabelSize(0.08);
    histo_residualsCOV1_10MeV_pc->GetXaxis()->SetTitleSize(0.08);
    histo_residualsCOV1_10MeV_pc->GetYaxis()->SetRangeUser(-100, 100);
    histo_residualsCOV1_10MeV_pc->SetFillColorAlpha(kBlue+1, 0.2);
    histo_residualsCOV1_10MeV_pc->Draw("E3");
    
    cCOV1->Write();
    cCOV1->Close();
    
    cCOV2->cd(2);
    gPad->SetPad(0.5, 0.23, 1, 1.);
    gPad->SetLogy();
    gPad->SetRightMargin(0.02);
    
    if (!TString(histo_dataCOV2_10MeV_2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataCOV2_10MeV_2->Scale(24*3600.);
    histo_dataCOV2_10MeV_2->SetTitle("Top Ge (COV2) - 0-4MeV");
    histo_dataCOV2_10MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2_10MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2_10MeV_2->Draw();
    histo_dataCOV2_10MeV_2->GetXaxis()->SetRangeUser(0, 4);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    Hits_to_dru(histo_simu_muCOV2_10MeV_2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV2_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_nCOV2_10MeV_2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV2_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaCOV2_10MeV_2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV2_10MeV_2->GetNbinsX());
    TH1D* histo_simu_sumCOV2_10MeV_2 =(TH1D*) histo_simu_gammaCOV2_10MeV_2->Clone();
    histo_simu_sumCOV2_10MeV_2->Add(histo_simu_muCOV2_10MeV_2);
    histo_simu_sumCOV2_10MeV_2->Add(histo_simu_nCOV2_10MeV_2);
    
    histo_dataCOV2_10MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataCOV2_10MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataCOV2_10MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataCOV2_10MeV_2->GetXaxis()->SetTitleSize(0.03);

    histo_dataCOV2_10MeV_2->Draw("E0LF2");
    histo_simu_muCOV2_10MeV_2->SetLineColor(kCyan-5);
    histo_simu_muCOV2_10MeV_2->SetFillColorAlpha(kCyan-5, 0.20);
    histo_simu_muCOV2_10MeV_2->Draw("E0LF2SAME");
    histo_simu_nCOV2_10MeV_2->SetLineColor(kOrange-5);
    histo_simu_nCOV2_10MeV_2->SetFillColorAlpha(kOrange-5, 0.20);
    histo_simu_nCOV2_10MeV_2->Draw("E0LF2SAME");
    histo_simu_gammaCOV2_10MeV_2->SetLineColor(kMagenta-5);
    histo_simu_gammaCOV2_10MeV_2->SetFillColorAlpha(kMagenta-5, 0.20);
    histo_simu_gammaCOV2_10MeV_2->Draw("E0LF2SAME");
    histo_simu_AmCOV2_10MeV_2->SetLineColor(kRed-7);
    histo_simu_AmCOV2_10MeV_2->SetFillColorAlpha(kRed-7, 0.2);
    histo_simu_AmCOV2_10MeV_2->Draw("E0LF2SAME");
    histo_simu_sumCOV2_10MeV_2->SetLineColor(kRed);
    histo_simu_sumCOV2_10MeV_2->Draw("E0LF2SAME");
    
    TH1D* histo_residualsCOV2_10MeV_pc = (TH1D*) histo_dataCOV2_10MeV_2->Clone();
    chi2=0;
    for (int i = 1; i<histo_residualsCOV2_10MeV_pc->GetNbinsX(); i++){
        Value_i = histo_dataCOV2_10MeV_2->GetBinContent(i);
        Diff_i = histo_dataCOV2_10MeV_2->GetBinContent(i)-histo_simu_sumCOV2_10MeV_2->GetBinContent(i);
        Error_i = sqrt(histo_dataCOV2_10MeV_2->GetBinError(i)*histo_dataCOV2_10MeV_2->GetBinError(i)+histo_simu_sumCOV2_10MeV_2->GetBinError(i)*histo_simu_sumCOV2_10MeV_2->GetBinError(i));
        histo_residualsCOV2_10MeV_pc->SetBinContent(i, Diff_i/Value_i*100);
        histo_residualsCOV2_10MeV_pc->SetBinError(i, Error_i/Value_i*100);
        chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
    }

    TLatex *Latex_COV2_10MeV = new TLatex(3, 1000, Form("Chi2/dof = %0.2f / %d = %0.2f", chi2, histo_residualsCOV2_10MeV_pc->GetNbinsX()-4, 1.*chi2/(histo_residualsCOV2_10MeV_pc->GetNbinsX()-1)));
    Latex_COV2_10MeV->SetTextSize(.030);
    Latex_COV2_10MeV->SetTextColor(1);
    Latex_COV2_10MeV->SetTextAlign(22);
    Latex_COV2_10MeV->Draw();
    
    cCOV2->cd(4);
    gPad->SetPad(0.5, 0., 1, 0.3);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetGridy();
    
    histo_residualsCOV2_10MeV_pc->GetYaxis()->SetTitle("Residuals [%]");
    histo_residualsCOV2_10MeV_pc->SetTitle("");
    histo_residualsCOV2_10MeV_pc->GetYaxis()->SetTitleSize(0.08);
    histo_residualsCOV2_10MeV_pc->GetYaxis()->SetTitleOffset(0.5);
    histo_residualsCOV2_10MeV_pc->GetXaxis()->SetTitleOffset(1);
    histo_residualsCOV2_10MeV_pc->GetYaxis()->SetLabelSize(0.08);
    histo_residualsCOV2_10MeV_pc->GetXaxis()->SetLabelSize(0.08);
    histo_residualsCOV2_10MeV_pc->GetXaxis()->SetTitleSize(0.08);
    histo_residualsCOV2_10MeV_pc->GetYaxis()->SetRangeUser(-100, 100);
    histo_residualsCOV2_10MeV_pc->SetFillColorAlpha(kBlue+1, 0.2);
    histo_residualsCOV2_10MeV_pc->Draw("E3");
    
    cCOV2->Write();
    cCOV2->Close();
    
    cNTD->cd(2);
    gPad->SetPad(0.5, 0.23, 1., 1.);
    gPad->SetLogy();
    gPad->SetRightMargin(0.02);
    
    if (!TString(histo_dataNTD_10MeV_2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataNTD_10MeV_2->Scale(24*3600.);
    histo_dataNTD_10MeV_2->SetTitle("Mid LWO (NTD) - 0-6MeV");
    histo_dataNTD_10MeV_2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataNTD_10MeV_2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataNTD_10MeV_2->Draw();
    histo_dataNTD_10MeV_2->GetXaxis()->SetRangeUser(0, 6);
    histo_dataNTD_10MeV_2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    Hits_to_dru(histo_simu_muNTD_10MeV_2, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_nNTD_10MeV_2, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataNTD_10MeV_2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaNTD_10MeV_2, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataNTD_10MeV_2->GetNbinsX());
    TH1D* histo_simu_sumNTD_10MeV_2 =(TH1D*) histo_simu_gammaNTD_10MeV_2->Clone();
    histo_simu_sumNTD_10MeV_2->Add(histo_simu_muNTD_10MeV_2);
    histo_simu_sumNTD_10MeV_2->Add(histo_simu_nNTD_10MeV_2);
    
    histo_dataNTD_10MeV_2->GetYaxis()->SetTitleSize(0.03);
    histo_dataNTD_10MeV_2->GetYaxis()->SetTitleOffset(1.7);
    histo_dataNTD_10MeV_2->GetYaxis()->SetLabelSize(0.03);
    histo_dataNTD_10MeV_2->GetXaxis()->SetLabelSize(0.03);
    histo_dataNTD_10MeV_2->GetXaxis()->SetTitleSize(0.03);

    histo_dataNTD_10MeV_2->Draw("E0LF2");
    histo_simu_muNTD_10MeV_2->SetLineColor(kCyan-5);
    histo_simu_muNTD_10MeV_2->SetFillColorAlpha(kCyan-5, 0.20);
    histo_simu_muNTD_10MeV_2->Draw("E0LF2SAME");
    histo_simu_nNTD_10MeV_2->SetLineColor(kOrange-5);
    histo_simu_nNTD_10MeV_2->SetFillColorAlpha(kOrange-5, 0.20);
    histo_simu_nNTD_10MeV_2->Draw("E0LF2SAME");
    histo_simu_gammaNTD_10MeV_2->SetLineColor(kMagenta-5);
    histo_simu_gammaNTD_10MeV_2->SetFillColorAlpha(kMagenta-5, 0.20);
    histo_simu_gammaNTD_10MeV_2->Draw("E0LF2SAME");
    histo_simu_AmNTD_10MeV_2->SetLineColor(kRed-7);
    histo_simu_AmNTD_10MeV_2->SetFillColorAlpha(kRed-7, 0.2);
    histo_simu_AmNTD_10MeV_2->Draw("E0LF2SAME");
    histo_simu_sumNTD_10MeV_2->SetLineColor(kRed);
    histo_simu_sumNTD_10MeV_2->Draw("E0LF2SAME");

    
    TH1D* histo_residualsNTD_10MeV_pc = (TH1D*) histo_dataNTD_10MeV_2->Clone();
    chi2=0;
    for (int i = 1; i<histo_residualsNTD_10MeV_pc->GetNbinsX(); i++){
        Value_i = histo_dataNTD_10MeV_2->GetBinContent(i);
        Diff_i = histo_dataNTD_10MeV_2->GetBinContent(i)-histo_simu_sumNTD_10MeV_2->GetBinContent(i);
        Error_i = sqrt(histo_dataNTD_10MeV_2->GetBinError(i)*histo_dataNTD_10MeV_2->GetBinError(i)+histo_simu_sumNTD_10MeV_2->GetBinError(i)*histo_simu_sumNTD_10MeV_2->GetBinError(i));
        if (Value_i!=0){histo_residualsNTD_10MeV_pc->SetBinContent(i, Diff_i/Value_i*100);
            histo_residualsNTD_10MeV_pc->SetBinError(i, Error_i/Value_i*100);
            chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
        }
        else histo_residualsNTD_10MeV_pc->SetBinContent(i, 0);
        
    }

    TLatex *Latex_NTD_10MeV = new TLatex(4, 1000, Form("Chi2/dof = %0.2f / %d = %0.2f", chi2, histo_residualsNTD_10MeV_pc->GetNbinsX()-4, chi2/(histo_residualsNTD_10MeV_pc->GetNbinsX()-1)));
    Latex_NTD_10MeV->SetTextSize(.030);
    Latex_NTD_10MeV->SetTextColor(1);
    Latex_NTD_10MeV->SetTextAlign(22);
    Latex_NTD_10MeV->Draw();
    
    cNTD->cd(4);
    gPad->SetPad(0.5, 0., 1., 0.3);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.2);
    gPad->SetGridy();
    
    histo_residualsNTD_10MeV_pc->GetYaxis()->SetTitle("Residuals [%]");
    histo_residualsNTD_10MeV_pc->SetTitle("");
    histo_residualsNTD_10MeV_pc->GetYaxis()->SetTitleSize(0.08);
    histo_residualsNTD_10MeV_pc->GetYaxis()->SetTitleOffset(0.5);
    histo_residualsNTD_10MeV_pc->GetXaxis()->SetTitleOffset(1);
    histo_residualsNTD_10MeV_pc->GetYaxis()->SetLabelSize(0.08);
    histo_residualsNTD_10MeV_pc->GetXaxis()->SetLabelSize(0.08);
    histo_residualsNTD_10MeV_pc->GetXaxis()->SetTitleSize(0.08);
    histo_residualsNTD_10MeV_pc->GetYaxis()->SetRangeUser(-100, 100);
    histo_residualsNTD_10MeV_pc->SetFillColorAlpha(kBlue+1, 0.2);
    histo_residualsNTD_10MeV_pc->Draw("E3");
    
    cNTD->Write();
    cNTD->Close();
}

void Save_Fitted_Simu_range(TFile* Fileout, int Emax, int Emax_simu, TString Trigger, vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am, bool Syst){
    
    TString Str_Emax = Form("%d", Emax);
    TString Str_Emax_simu = Form("%d", Emax_simu);
    
    TString Str_Trig, histo_name;
    
    if (Trigger=="Triggered"){
        Str_Trig = "OrsayOVAlgoDir/Triggered/";
        histo_name = "nEvis_triggered";
    }
    else if (Trigger=="coinc_Bot_only"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Bot_only/";
        histo_name = "nEvis_coinc_Bot_only";
    }
    else if (Trigger=="coinc_Mid_only"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Mid_only/";
        histo_name = "nEvis_coinc_Mid_only";
    }
    else if (Trigger=="coinc_Mid_&_Bot"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Mid_&_Bot/";
        histo_name = "nEvis_coinc_Mid_&_Bot";
    }
    else if (Trigger=="coinc_Top_only"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Top_only/";
        histo_name = "nEvis_coinc_Top_only";
    }
    else if (Trigger=="coinc_Top_&_Bot"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Top_&_Bot/";
        histo_name = "nEvis_coinc_Top_&_Bot";
    }
    else if (Trigger=="coinc_Top_&_Mid"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Top_&_Mid/";
        histo_name = "nEvis_coinc_Top_&_Mid";
    }
    else if (Trigger=="coinc_Top_&_Mid_&_Bot"){
        Str_Trig = "OrsayOVAlgoDir/Coincidences/coinc_Top_&_Mid_&_Bot/";
        histo_name = "nEvis_coinc_Top_&_Mid_&_Bot";
    }
    else {
        std::cerr << "Wrong trigger name: " +Trigger+ "\n" << std::endl;
    }
    
    TH1D* histo_dataCOV1 = (TH1D*) fileData->Get("BotGe/BotGe_Norm_Spectrum_0-"+Str_Emax+"MeV");
    histo_dataCOV1->SetTitle("BotGe_Norm_Spectrum_0-"+Str_Emax+"MeV");
    TH1D* histo_dataCOV2 = (TH1D*) fileData->Get("TopGe/TopGe_Norm_Spectrum_0-"+Str_Emax+"MeV");
    histo_dataCOV2->SetTitle("TopGe_Norm_Spectrum_0-"+Str_Emax+"MeV");
    
    TH2D* V_dataCOV1 = Stat_Covariance_TH2D(histo_dataCOV1);
    V_dataCOV1->SetTitle("V_dataCOV1");
    TH2D* V_dataCOV2 = Stat_Covariance_TH2D(histo_dataCOV2);
    V_dataCOV2->SetTitle("V_dataCOV2");

    if (Syst){
        V_dataCOV1 = (TH2D*) fileData->Get("BotGe/BotGe_CovMat_0-"+Str_Emax+"MeV");
        V_dataCOV1->SetTitle("V_dataCOV1_"+Str_Emax+"MeV");
        V_dataCOV2 = (TH2D*) fileData->Get("TopGe/TopGe_CovMat_0-"+Str_Emax+"MeV");
        V_dataCOV2->SetTitle("V_dataCOV2_"+Str_Emax+"MeV");
    }

    TH1D*  histo_simu_muCOV1 = (TH1D*) filein_mu->Get(Str_Trig+"Vol_Bot_Ge/Evis/"+histo_name+"_Vol_Bot_Ge_0-"+Str_Emax_simu+"MeV");
    TH1D*  histo_simu_muCOV2 = (TH1D*) filein_mu->Get(Str_Trig+"Vol_Top_Ge/Evis/"+histo_name+"_Vol_Top_Ge_0-"+Str_Emax_simu+"MeV");
    
    TH1D* histo_simu_gammaCOV1 = (TH1D*) filein_gamma->Get(Str_Trig+"Vol_Bot_Ge/Evis/"+histo_name+"_Vol_Bot_Ge_0-"+Str_Emax_simu+"MeV");
    TH1D* histo_simu_gammaCOV2 = (TH1D*) filein_gamma->Get(Str_Trig+"Vol_Top_Ge/Evis/"+histo_name+"_Vol_Top_Ge_0-"+Str_Emax_simu+"MeV");

    TH1D* histo_simu_nCOV1 = (TH1D*) filein_n->Get(Str_Trig+"Vol_Bot_Ge/Evis/"+histo_name+"_Vol_Bot_Ge_0-"+Str_Emax_simu+"MeV");
    TH1D* histo_simu_nCOV2 = (TH1D*) filein_n->Get(Str_Trig+"Vol_Top_Ge/Evis/"+histo_name+"_Vol_Top_Ge_0-"+Str_Emax_simu+"MeV");

    TH1D* histo_simu_AmCOV1 = (TH1D*) filein_Am->Get(Str_Trig+"Vol_Bot_Ge/Evis/"+histo_name+"_Vol_Bot_Ge_0-"+Str_Emax_simu+"MeV");
    TH1D* histo_simu_AmCOV2 = (TH1D*) filein_Am->Get(Str_Trig+"Vol_Top_Ge/Evis/"+histo_name+"_Vol_Top_Ge_0-"+Str_Emax_simu+"MeV");
    
    histo_simu_muCOV1=GoodFormat(histo_simu_muCOV1, "histo_simu_muCOV1", histo_dataCOV1->GetNbinsX(), 0, Emax);
    histo_simu_muCOV2=GoodFormat(histo_simu_muCOV2, "histo_simu_muCOV2", histo_dataCOV2->GetNbinsX(), 0, Emax);

    histo_simu_gammaCOV1=GoodFormat(histo_simu_gammaCOV1, "histo_simu_gammaCOV1", histo_dataCOV1->GetNbinsX(), 0, Emax);
    histo_simu_gammaCOV2=GoodFormat(histo_simu_gammaCOV2, "histo_simu_gammaCOV2", histo_dataCOV2->GetNbinsX(), 0, Emax);

    histo_simu_nCOV1=GoodFormat(histo_simu_nCOV1, "histo_simu_nCOV1", histo_dataCOV1->GetNbinsX(), 0, Emax);
    histo_simu_nCOV2=GoodFormat(histo_simu_nCOV2, "histo_simu_nCOV2", histo_dataCOV2->GetNbinsX(), 0, Emax);

    histo_simu_AmCOV1=GoodFormat(histo_simu_AmCOV1, "histo_simu_AmCOV1", histo_dataCOV1->GetNbinsX(), 0, Emax);
    histo_simu_AmCOV2=GoodFormat(histo_simu_AmCOV2, "histo_simu_AmCOV2", histo_dataCOV2->GetNbinsX(), 0, Emax);
    
    if (!TString(histo_dataCOV1->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataCOV1->Scale(24*3600.);
    histo_dataCOV1->SetTitle("Bot Ge (COV1) - 0-"+Str_Emax+"MeV");
    histo_dataCOV1->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV1->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    
    Hits_to_dru(histo_simu_muCOV1, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV1->GetNbinsX());
    Hits_to_dru(histo_simu_nCOV1, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV1->GetNbinsX());
    Hits_to_dru(histo_simu_gammaCOV1, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV1->GetNbinsX());
    Hits_to_dru(histo_simu_AmCOV1, 1e3, MASS_Ge, NumEvents[3], fit_Am.Parameter(3),histo_dataCOV1->GetNbinsX());
    
    //sum
    TH1D* histo_simu_sumCOV1 =(TH1D*) histo_simu_gammaCOV1->Clone();
    histo_simu_sumCOV1->Add(histo_simu_muCOV1);
    histo_simu_sumCOV1->Add(histo_simu_nCOV1);
    histo_simu_sumCOV1->Add(histo_simu_AmCOV1);
    
    if (Trigger=="Triggered"){
        Fileout->cd("COV1/Data/Triggered/");
        histo_dataCOV1->Write("BotGe_Norm_Spectrum_0-"+Str_Emax+"MeV");
        if (Syst){
            SetErrors_from_CovMat(histo_dataCOV1, V_dataCOV1, 24*3600.);
            histo_dataCOV1->Write("BotGe_Norm_Spectrum_Syst+StatErr_0-"+Str_Emax+"MeV");
            V_dataCOV1->Delete();
        }
        histo_dataCOV1->Delete();
    }
    
    Fileout->cd("COV1/Muons/"+Trigger);
    histo_simu_muCOV1->Write("histo_simu_muCOV1_0-"+Str_Emax+"MeV");
    histo_simu_muCOV1->Delete();
    
    Fileout->cd("COV1/Gammas/"+Trigger);
    histo_simu_gammaCOV1->Write("histo_simu_gammaCOV1_0-"+Str_Emax+"MeV");
    histo_simu_gammaCOV1->Delete();
    
    Fileout->cd("COV1/Neutrons/"+Trigger);
    histo_simu_nCOV1->Write("histo_simu_nCOV1_0-"+Str_Emax+"MeV");
    histo_simu_nCOV1->Delete();
    
    Fileout->cd("COV1/Am/"+Trigger);
    histo_simu_AmCOV1->Write("histo_simu_AmCOV1_0-"+Str_Emax+"MeV");
    histo_simu_AmCOV1->Delete();
    
    Fileout->cd("COV1/Sum/"+Trigger);
    histo_simu_sumCOV1->Write("histo_simu_sumCOV1_0-"+Str_Emax+"MeV");
    histo_simu_sumCOV1->Delete();
    
    if (!TString(histo_dataCOV2->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataCOV2->Scale(24*3600.);
    histo_dataCOV2->SetTitle("Top Ge (COV2) - 0-"+Str_Emax+"MeV");
    histo_dataCOV2->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
    histo_dataCOV2->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
    histo_dataCOV2->GetXaxis()->SetRangeUser(0, 20);
    histo_dataCOV2->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
    
    Hits_to_dru(histo_simu_muCOV2, 1e3, MASS_Ge, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataCOV2->GetNbinsX());
    Hits_to_dru(histo_simu_nCOV2, 1e3, MASS_Ge, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataCOV2->GetNbinsX());
    Hits_to_dru(histo_simu_gammaCOV2, 1e3, MASS_Ge, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataCOV2->GetNbinsX());
    
    TH1D* histo_simu_sumCOV2 =(TH1D*) histo_simu_gammaCOV2->Clone();
    histo_simu_sumCOV2->Add(histo_simu_muCOV2);
    histo_simu_sumCOV2->Add(histo_simu_nCOV2);
    
    if (Trigger=="Triggered"){
        Fileout->cd("COV2/Data/Triggered");
        histo_dataCOV2->Write("TopGe_Norm_Spectrum_0-"+Str_Emax+"MeV");
        if (Syst){
            SetErrors_from_CovMat(histo_dataCOV2, V_dataCOV2, 24*3600.);
            histo_dataCOV2->Write("TopGe_Norm_Spectrum_Syst+StatErr_0-"+Str_Emax+"MeV");
            V_dataCOV2->Delete();
        }
        histo_dataCOV2->Delete();
    }
    Fileout->cd("COV2/Muons/"+Trigger);
    histo_simu_muCOV2->Write("histo_simu_muCOV2_0-"+Str_Emax+"MeV");
    histo_simu_muCOV2->Delete();

    Fileout->cd("COV2/Gammas/"+Trigger);
    histo_simu_gammaCOV2->Write("histo_simu_gammaCOV2_0-"+Str_Emax+"MeV");
    histo_simu_gammaCOV2->Delete();

    Fileout->cd("COV2/Neutrons/"+Trigger);
    histo_simu_nCOV2->Write("histo_simu_nCOV2_0-"+Str_Emax+"MeV");
    histo_simu_nCOV2->Delete();

    Fileout->cd("COV2/Sum/"+Trigger);
    histo_simu_sumCOV2->Write("histo_simu_sumCOV2_0-"+Str_Emax+"MeV");
    histo_simu_sumCOV2->Delete();

    histo_simu_AmCOV2->Delete();
    
    // NTD
    TH1D* histo_dataNTD = new TH1D();
    TH1D* histo_dataNTD_withAcc = new TH1D();
    TH1D* histo_simu_muNTD = new TH1D();
    TH1D* histo_simu_gammaNTD = new TH1D();
    TH1D* histo_simu_nNTD = new TH1D();
    TH1D* histo_simu_AmNTD = new TH1D();
    
    TH2D* V_dataNTD;
    
    if (Emax>=10){
        if (Trigger=="Triggered"){
            histo_dataNTD= (TH1D*) fileData->Get("LWO/LWO_Norm_Spectrum_0-"+Str_Emax+"MeV");
            histo_dataNTD->SetTitle("histo_dataNTD");
            if (Syst){
                V_dataNTD = (TH2D*) fileData->Get("LWO/LWO_CovMat_0-"+Str_Emax+"MeV");
                V_dataNTD->SetTitle("V_dataLWO_"+Str_Emax+"MeV");
            }
        }
        else {histo_dataNTD= (TH1D*) fileData->Get("LWO/LWO-Vetoed_Norm_Spectrum_0-"+Str_Emax+"MeV");
            histo_dataNTD->SetTitle("histo_dataNTD");
            histo_dataNTD_withAcc = (TH1D*) fileData->Get("LWO/LWO-Vetoed_withAcc_Norm_Spectrum_0-"+Str_Emax+"MeV");
            histo_dataNTD_withAcc->SetTitle("histo_dataNTD_withAcc");
            
            if (!TString(histo_dataNTD_withAcc->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataNTD_withAcc->Scale(24*3600.);
            histo_dataNTD_withAcc->SetTitle("Mid LWO (NTD) - 0-"+Str_Emax+"MeV");
            histo_dataNTD_withAcc->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
            histo_dataNTD_withAcc->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
            histo_dataNTD_withAcc->GetXaxis()->SetRangeUser(0, Emax);
            histo_dataNTD_withAcc->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
        }
        histo_simu_muNTD = (TH1D*) filein_mu->Get(Str_Trig+"Vol_Mid_LWO/Evis/"+histo_name+"_Vol_Mid_LWO_0-"+Str_Emax_simu+"MeV");
        histo_simu_gammaNTD = (TH1D*) filein_gamma->Get(Str_Trig+"Vol_Mid_LWO/Evis/"+histo_name+"_Vol_Mid_LWO_0-"+Str_Emax_simu+"MeV");
        histo_simu_nNTD = (TH1D*) filein_n->Get(Str_Trig+"Vol_Mid_LWO/Evis/"+histo_name+"_Vol_Mid_LWO_0-"+Str_Emax_simu+"MeV");
        histo_simu_AmNTD = (TH1D*) filein_Am->Get(Str_Trig+"Vol_Mid_LWO/Evis/"+histo_name+"_Vol_Mid_LWO_0-"+Str_Emax_simu+"MeV");
        
        histo_simu_muNTD=GoodFormat(histo_simu_muNTD, "histo_simu_muNTD", histo_dataNTD->GetNbinsX(), 0, Emax);
        histo_simu_gammaNTD=GoodFormat(histo_simu_gammaNTD, "histo_simu_gammaNTD", histo_dataNTD->GetNbinsX(), 0, Emax);
        histo_simu_nNTD=GoodFormat(histo_simu_nNTD, "histo_simu_nNTD", histo_dataNTD->GetNbinsX(), 0, Emax, true);
        histo_simu_AmNTD=GoodFormat(histo_simu_AmNTD, "histo_simu_AmNTD", histo_dataNTD->GetNbinsX(), 0, Emax);
        
        if (!TString(histo_dataNTD->GetYaxis()->GetTitle()).EqualTo("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})")) histo_dataNTD->Scale(24*3600.);
        histo_dataNTD->SetTitle("Mid LWO (NTD) - 0-"+Str_Emax+"MeV");
        histo_dataNTD->GetYaxis()->SetTitle("Differential Rate (.d^{-1}.kg^{-1}.keV^{-1})");
        histo_dataNTD->GetXaxis()->SetTitle("Reconstructed energy [MeV]");
        histo_dataNTD->GetXaxis()->SetRangeUser(0, Emax);
        histo_dataNTD->GetYaxis()->SetRangeUser(0.5*1e-2, 1e5);
        
        Hits_to_dru(histo_simu_muNTD, 1e3, MASS_LWO, NumEvents[0], fit_muon.Parameter(0)*PLANE_SURF, histo_dataNTD->GetNbinsX());
        Hits_to_dru(histo_simu_nNTD, 1e3, MASS_LWO, NumEvents[2], fit_neutron.Parameter(2)*PLANE_SURF,histo_dataNTD->GetNbinsX());
        Hits_to_dru(histo_simu_gammaNTD, 1e3, MASS_LWO, NumEvents[1], fit_gamma.Parameter(1)*PLANE_SURF,histo_dataNTD->GetNbinsX());
        
        TH1D* histo_simu_sumNTD=(TH1D*) histo_simu_gammaNTD->Clone();
        histo_simu_sumNTD->Add(histo_simu_muNTD);
        histo_simu_sumNTD->Add(histo_simu_nNTD);

        if (Trigger=="Triggered"){
            Fileout->cd("NTD/Data/Triggered");
            histo_dataNTD->Write("LWO_Norm_Spectrum_0-"+Str_Emax+"MeV", TObject::kOverwrite);
            if (Syst){
                SetErrors_from_CovMat(histo_dataNTD, V_dataNTD, 24*3600.);
                histo_dataNTD->Write("LWO_Norm_Spectrum_Syst+StatErr_0-"+Str_Emax+"MeV", TObject::kOverwrite);
                V_dataNTD->Delete();
            }
            histo_dataNTD->Delete();
        }
        else {
            Fileout->cd("NTD/Data/Accepted");
            histo_dataNTD_withAcc->Write("LWO-Vetoed_withAcc_Norm_Spectrum_0-"+Str_Emax+"MeV", TObject::kOverwrite);
            histo_dataNTD_withAcc->Delete();
            histo_dataNTD->Write("LWO-Vetoed_Norm_Spectrum_0-"+Str_Emax+"MeV", TObject::kOverwrite);
            histo_dataNTD->Delete();
        }
        Fileout->cd("NTD/Muons/"+Trigger);
        histo_simu_muNTD->Write("histo_simu_muNTD_0-"+Str_Emax+"MeV");
        histo_simu_muNTD->Delete();

        Fileout->cd("NTD/Gammas/"+Trigger);
        histo_simu_gammaNTD->Write("histo_simu_gammaNTD_0-"+Str_Emax+"MeV");
        histo_simu_gammaNTD->Delete();

        Fileout->cd("NTD/Neutrons/"+Trigger);
        histo_simu_nNTD->Write("histo_simu_nNTD_0-"+Str_Emax+"MeV");
        histo_simu_nNTD->Delete();

        Fileout->cd("NTD/Sum/"+Trigger);
        histo_simu_sumNTD->Write("histo_simu_sumNTD_0-"+Str_Emax+"MeV");
        histo_simu_sumNTD->Delete();

        histo_simu_AmNTD->Delete();
    }
    
    cout<<"Fit result for following coinc: "<<Trigger<<" for E in [0-"+Str_Emax+"] MeV saved!"<<endl;
}

void Save_Fitted(TString FileoutName, vector<Double_t> NumEvents, TFile* fileData, TFile* filein_mu, TFile* filein_gamma, TFile* filein_n, TFile* filein_Am, ROOT::Fit::FitResult fit_muon, ROOT::Fit::FitResult fit_gamma, ROOT::Fit::FitResult fit_neutron, ROOT::Fit::FitResult fit_Am, bool Syst=false){
    
    TFile* Fileout = new TFile(FileoutName, "UPDATE");
    
    std::vector<TString> Crystals = {"COV1", "COV2", "NTD"};
    std::vector<TString> Coincidences = {"Triggered", "coinc_Bot_only", "coinc_Mid_only", "coinc_Mid_&_Bot", "coinc_Top_only", "coinc_Top_&_Bot", "coinc_Top_&_Mid", "coinc_Top_&_Mid_&_Bot"};
    
    TH1D* Norm = new TH1D("NormConstants_fitResult", "NormConstants_fitResult", 4, 0, 4);
    Norm->SetBinContent(1, fit_muon.Parameter(0));
    Norm->SetBinError(1, fit_muon.ParError(0));
    Norm->SetBinContent(2, fit_gamma.Parameter(1));
    Norm->SetBinError(2, fit_gamma.ParError(1));
    Norm->SetBinContent(3, fit_neutron.Parameter(2));
    Norm->SetBinError(3, fit_neutron.ParError(2));
    Norm->SetBinContent(4, fit_Am.Parameter(3));
    Norm->SetBinError(4, fit_Am.ParError(3));
    Norm->GetXaxis()->SetBinLabel(1, "Muons");
    Norm->GetXaxis()->SetBinLabel(2, "Gammas");
    Norm->GetXaxis()->SetBinLabel(3, "Neutrons");
    Norm->GetXaxis()->SetBinLabel(4, "Am");
    Norm->Write();
    
    TString Crystal, Coinc;
    for (int c =0; c<Crystals.size(); c++){
        Crystal=Crystals[c];
        Fileout->mkdir(Crystal);
        Fileout->mkdir(Crystal+"/Data");
        Fileout->mkdir(Crystal+"/Muons");
        Fileout->mkdir(Crystal+"/Gammas");
        Fileout->mkdir(Crystal+"/Neutrons");
        Fileout->mkdir(Crystal+"/Am");
        Fileout->mkdir(Crystal+"/Sum");
        for (int d =0; d<Coincidences.size(); d++){
            Coinc=Coincidences[d];
    
            Fileout->mkdir(Crystal+"/Data/"+Coinc);
            Fileout->mkdir(Crystal+"/Muons/"+Coinc);
            Fileout->mkdir(Crystal+"/Gammas/"+Coinc);
            Fileout->mkdir(Crystal+"/Neutrons/"+Coinc);
            Fileout->mkdir(Crystal+"/Am/"+Coinc);
            Fileout->mkdir(Crystal+"/Sum/"+Coinc);
        }
    }
    Fileout->mkdir("NTD/Data/Accepted");
    
    for (int d=0; d<Coincidences.size(); d++){
        Coinc=Coincidences[d];
        Save_Fitted_Simu_range(Fileout, 20, 100, Coinc, NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron, fit_Am, Syst);
        Save_Fitted_Simu_range(Fileout, 10, 10, Coinc, NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron, fit_Am, Syst);
        Save_Fitted_Simu_range(Fileout, 4, 10, Coinc, NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron, fit_Am, Syst);
        Save_Fitted_Simu_range(Fileout, 1, 1, Coinc, NumEvents, fileData, filein_mu, filein_gamma, filein_n, filein_Am, fit_muon, fit_gamma, fit_neutron, fit_Am, Syst);
    }
    //Fileout->Close();
}

void Int_Rate(TH1D* histo, vector<double>& Int_Rate, double Estart=0, double Eend=0, double MassTarget = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE, bool Verbose_level=1){
    
    double dE = histo->GetBinWidth(1);
    double Nbins = histo->GetNbinsX();
    
    double BinFirst = histo->FindBin(Estart);
    double BinLast= histo->FindBin(Eend);
    TString Unit= histo->GetXaxis()->GetTitle();
    
    if (Estart==0){
        BinFirst= 1;
    }
    
    if (Eend==0) {
        BinLast = Nbins;
    }
    
    double Integral, Error;
    
    Integral = histo->IntegralAndError(BinFirst, BinLast, Error);
    
    double Rate =(Integral*rate)/num_events;
    double Err_Rate =(Error*rate)/num_events;
    
    Int_Rate.push_back(Rate);
    Int_Rate.push_back(Err_Rate);
    
    
    if (Verbose_level){
    cout<<Estart << " - " << Eend << " ("<< Unit<< ") --- "<< Integral << " +- " << Error << " events >>> " << Rate << " +-  "<< Err_Rate << " ev/s ---- " << Rate/MassTarget << " +- " << Err_Rate/MassTarget <<" ev/s/kg" << endl;
    }
    
    return;
}

void Norm_in_cps_simu_veto(TString Filename, TString NameFileout, bool Am = FALSE) {
    
    TString TxtFilename = NameFileout;
    TPMERegexp("\\.root").Substitute(TxtFilename,".txt");
    
    std::streambuf *psbuf, *backup;
    std::ofstream filestr;
    filestr.open (TxtFilename);
    backup = std::cout.rdbuf();     // back up cout's streambuf
    psbuf = filestr.rdbuf();        // get file's streambuf
    std::cout.rdbuf(psbuf);         // assign streambuf to cout
    
    TFile *Fileout = new TFile(NameFileout, "RECREATE");
    Fileout->Close();
    
    TCanvas *c1 = new TCanvas();
    
    TFile *filein = new TFile(Filename, "READ");
    int Hits_i, Error_i, NumBin;
    
    double NumEv = SetNumEvents(Filename, Am);
    if (Am){
        cout << "Num events = " <<NumEv << " in 60keV peak -- expected number of events in 60keV peak: " << FLUX << " ev/s" << endl;
    }
    else {cout << "Num events = " <<NumEv << " -- Flux : " << FLUX << " ev/cm2/s" << endl;
    }
    cout << " " <<endl;
    
    double Tot_Rate_thresh=0;
    double Error_Rate_thresh_sq=0;
    double Tot_Rate_fixed=0;
    double Error_Rate_fixed_sq=0;
    
    Fileout = new TFile(NameFileout, "UPDATE");
    Fileout->mkdir("Accepted");
    Fileout->mkdir("Triggered");
    
    int Numb_AnalyzedE=NUM_E_ANALYZED;
    TString EFolder[]={"Edep", "Evis"};
    
    for (int i=0; i<Numb_AnalyzedE; i++){
        cout << ""<<endl;
        cout << "------------------------------- "+EFolder[i]+" ------------------------------"<<endl;
        cout<<""<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        cout << "                                    Accepted"<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;

        Fileout->cd("Accepted");
        
        //TopGe
        Fileout->mkdir("Accepted/Top_Ge/"+EFolder[i]);
        Fileout->cd("Accepted/Top_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Top_Ge",
                        "OrsayOVAlgoDir/Accepted/Vol_Top_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_accepted_Vol_Top_Ge",
                        MASS_Ge, THRESH_TOP, 30, NumEv, RATE);

        //LWO
        Fileout->mkdir("Accepted/Mid_LWO/"+EFolder[i]);
        Fileout->cd("Accepted/Mid_LWO/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Mid_LWO",
                        "OrsayOVAlgoDir/Accepted/Vol_Mid_LWO/"+EFolder[i]+"/n"+EFolder[i]+"_accepted_Vol_Mid_LWO",
                        MASS_LWO, THRESH_LWO, 50, NumEv, RATE);

        
        //Bot Ge
        Fileout->mkdir("Accepted/Bot_Ge/"+EFolder[i]);
        Fileout->cd("Accepted/Bot_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Bot_Ge",
                        "OrsayOVAlgoDir/Accepted/Vol_Bot_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_accepted_Vol_Bot_Ge",
                        MASS_Ge, THRESH_BOT, 30, NumEv, RATE);
        
        cout<<""<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        cout << "                                    Triggered"<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        
        //TopGe
        Fileout->mkdir("Triggered/Top_Ge/"+EFolder[i]);
        Fileout->cd("Triggered/Top_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Top_Ge",
                        "OrsayOVAlgoDir/Triggered/Vol_Top_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Top_Ge",
                        MASS_Ge, THRESH_TOP, 30, NumEv, RATE);

        //LWO
        Fileout->mkdir("Triggered/Mid_LWO/"+EFolder[i]);
        Fileout->cd("Triggered/Mid_LWO/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Mid_LWO",
                        "OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Mid_LWO",
                        MASS_LWO, THRESH_LWO, 50, NumEv, RATE);
        
        //Bot Ge
        Fileout->mkdir("Triggered/Bot_Ge/"+EFolder[i]);
        Fileout->cd("Triggered/Bot_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Bot_Ge",
                        "OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Bot_Ge",
                        MASS_Ge, THRESH_BOT, 30, NumEv, RATE
                        );
    }

    std::cout.rdbuf(backup);        // restore cout's original streambuf
    filestr.close();
    
    ifstream f(TxtFilename); cout<<f.rdbuf();
    f.close();

    Fileout->Close();
    c1->Close();
    return;
}

void save_CPS_histos(TFile *filein, TFile *fileout, TString Crystal, TString hName, double Crystal_Mass, double thresh_crys, double thresh_crys_fixed, double NumEv, double Rate){
    
    vector<double> Rate_0, Rate_thresh, Rate_fixed, trash;
    
    cout << " " <<endl;
    cout << Crystal + " - Target mass :" << Crystal_Mass << " kg -- Flux : " << FLUX << " ev/cm2/s" << endl;
    cout << "-------------------------------" <<endl;
    
    TH1D *h_cryst = (TH1D*)  filein->Get(hName+"_0-100keV");
    if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 0, 100, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, thresh_crys, 100, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_fixed, thresh_crys_fixed, 100, Crystal_Mass, NumEv, Rate, 1);
    
    Hits_to_cps(h_cryst, 1, Crystal_Mass, NumEv, Rate, 500);
    h_cryst->Write(Crystal+"_0-100keV", TObject::kOverwrite);
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-1MeV");
    if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 0.1, 1, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, 0.1, 1, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, Rate_fixed, 0.1, 1, Crystal_Mass, NumEv, Rate, 0);
    
    Hits_to_cps(h_cryst, 1e3, Crystal_Mass, NumEv, RATE, 500);
    h_cryst->Write(Crystal+"_0-1MeV", TObject::kOverwrite);
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-10MeV");
    if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 1, 10, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, 1, 10, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, Rate_fixed, 1, 10, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, trash, 5, 5.5, Crystal_Mass, NumEv, Rate, 1);
    
    Hits_to_cps(h_cryst, 1e3, Crystal_Mass, NumEv, Rate, 500);
    h_cryst->Write(Crystal+"_0-10MeV", TObject::kOverwrite);
    
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-10MeV");
    if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    double NumBin = h_cryst->FindBin(4)-h_cryst->FindBin(0);
    TH1D* new_h_cryst = new TH1D("new_h_cryst", "new_h_cryst", NumBin, 0, 4);
    
    double Hits_i,Error_i;
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_cryst->GetBinContent(i);
        Error_i = h_cryst->GetBinError(i);
        
        new_h_cryst->SetBinContent(i, Hits_i);
        new_h_cryst->SetBinError(i, Error_i);
    }
    
    new_h_cryst->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_cryst->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_cryst, 1e3, Crystal_Mass, NumEv, Rate, 200);
    new_h_cryst->Write(Crystal+"_0-4MeV");
    
    h_cryst->Delete();
    new_h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-100MeV");
    if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 10, 20, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, 10, 20, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, Rate_fixed, 10, 20, Crystal_Mass, NumEv, Rate, 0);
//    Int_Rate(h_cryst, trash, 4, 20, Crystal_Mass, NumEv, Rate, 1);
    
    Hits_to_cps(h_cryst, 1e3, Crystal_Mass, NumEv, Rate, 500);
    h_cryst->Write(Crystal+"_0-100MeV", TObject::kOverwrite);
    
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-100MeV");
    if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    NumBin = h_cryst->FindBin(20)-h_cryst->FindBin(0);
    new_h_cryst = new TH1D("new_h_cryst", "new_h_cryst", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_cryst->GetBinContent(i);
        Error_i = h_cryst->GetBinError(i);
        
        new_h_cryst->SetBinContent(i, Hits_i);
        new_h_cryst->SetBinError(i, Error_i);
    }
    
    new_h_cryst->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_cryst->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_cryst, 1e3, Crystal_Mass, NumEv, Rate, 100);
    new_h_cryst->Write(""+Crystal+"_0-20MeV");
    
    double Tot_Rate=0;
    double Error_Rate_sq=0;
    double Tot_Rate_thresh=0;
    double Error_Rate_thresh_sq=0;
    double Tot_Rate_fixed=0;
    double Error_Rate_fixed_sq=0;

    for (int i=0; i<Rate_thresh.size()/2; i++){
        Tot_Rate+=Rate_0[2*i];
        Error_Rate_sq+=Rate_0[2*i+1]*Rate_0[2*i+1];
        Tot_Rate_thresh+=Rate_thresh[2*i];
        Error_Rate_thresh_sq+=Rate_thresh[2*i+1]*Rate_thresh[2*i+1];
        Tot_Rate_fixed+=Rate_fixed[2*i];
        Error_Rate_fixed_sq+=Rate_fixed[2*i+1]*Rate_thresh[2*i+1];
    }

    double Error_Rate=sqrt(Error_Rate_sq);
    double Error_Rate_thresh=sqrt(Error_Rate_thresh_sq);
    double Error_Rate_fixed=sqrt(Error_Rate_fixed_sq);
    cout << ""<<endl;
    cout<<"Total rates :"<<endl;
    cout<<"0 keV - 100 MeV --- "<< Tot_Rate << " +-  "<< Error_Rate << " ev/s ---- " << Tot_Rate/Crystal_Mass << " +- " << Error_Rate/Crystal_Mass <<" ev/s/kg" << endl;
    cout<<thresh_crys << " keV - 100 MeV --- "<< Tot_Rate_thresh << " +-  "<< Error_Rate_thresh << " ev/s ---- " << Tot_Rate_thresh/Crystal_Mass << " +- " << Error_Rate_thresh/Crystal_Mass <<" ev/s/kg" << endl;
    cout<<thresh_crys_fixed<<" keV - 100 MeV --- "<< Tot_Rate_fixed << " +-  "<< Error_Rate_fixed << " ev/s ---- " << Tot_Rate_fixed/Crystal_Mass << " +- " << Error_Rate_fixed/Crystal_Mass <<" ev/s/kg" << endl;
    cout<<"************************************************************************"<<endl;
    
    Rate_thresh.clear();
    Rate_fixed.clear();
    
    h_cryst->Delete();
    new_h_cryst->Delete();
}

//----- print the integrated rate in the energy range
void Print_int_rate(TH1D* histo_0_10,TH1D* histo_0_20, Double_t Emin, Double_t Emax, bool Simple = true){
    Double_t Error10, Error20;
    if (Simple){
        if ((Emax>=10)&(Emin<=10)){
            cout<<histo_0_10->IntegralAndError(histo_0_10->FindBin(Emin), histo_0_10->FindBin(8), Error10, "width")*1e3+histo_0_20->IntegralAndError(histo_0_20->FindBin(8), histo_0_20->FindBin(Emax), Error20, "width")*1e3<<"; "<<sqrt(Error10*Error10+Error20*Error20)*1e3<<endl;
        }
        else if (Emin>10){
            cout<<histo_0_20->IntegralAndError(histo_0_20->FindBin(Emin), histo_0_20->FindBin(Emax), Error20, "width")*1e3<<" ; "<<Error20*1e3<<endl;
        }
        else if (Emax<10){
            cout<<histo_0_10->IntegralAndError(histo_0_10->FindBin(Emin), histo_0_10->FindBin(Emax), Error10, "width")*1e3<<" ; "<<Error10*1e3<<endl;
        }
    }
    else{
        if ((Emax>=10)&(Emin<=10)){
            cout<<"Integrated rate between "<<Emin<<" end "<<Emax<<" = "<<histo_0_10->IntegralAndError(histo_0_10->FindBin(Emin), histo_0_10->FindBin(8), Error10, "width")*1e3+histo_0_20->IntegralAndError(histo_0_20->FindBin(8), histo_0_20->FindBin(Emax), Error20, "width")*1e3<<" +- "<<sqrt(Error10*Error10+Error20*Error20)*1e3<< " ev/day/kg"<<endl;
        }
        else if (Emin>10){
            cout<<"Integrated rate between "<<Emin<<" end "<<Emax<<" = "<<histo_0_20->IntegralAndError(histo_0_20->FindBin(Emin), histo_0_20->FindBin(Emax), Error20, "width")*1e3<<" +- "<<Error20*1e3<< " ev/day/kg"<<endl;
        }
        else if (Emax<10){
            cout<<"Integrated rate between "<<Emin<<" end "<<Emax<<" = "<<histo_0_10->IntegralAndError(histo_0_10->FindBin(Emin), histo_0_10->FindBin(Emax), Error10, "width")*1e3<<" +- "<<Error10*1e3<< " ev/day/kg"<<endl;
        }
    }
}

//--------- Print all the interesting energy ranges for LWO
void Print_LWO_rates(TString FileinName){
    TFile* filein = new TFile(FileinName);
    cout<<"=============== Triggered Data =================== "<<endl;
    TH1D* histo_0_10 = (TH1D*) filein->Get("NTD/Data/Triggered/LWO-Vetoed_Spectrum_0-10MeV");
    TH1D* histo_0_20 = (TH1D*) filein->Get("NTD/Data/Triggered/LWO-Vetoed_Spectrum_0-20MeV");
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 20);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 3);
    Print_int_rate(histo_0_10, histo_0_20, 10, 20);
    Print_int_rate(histo_0_10, histo_0_20, 3, 10);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 0.5);
    cout<<"\n=============== Accepted Data =================== "<<endl;
    histo_0_10 = (TH1D*) filein->Get("NTD/Data/Accepted/LWO-Vetoed_Spectrum_0-10MeV");
    histo_0_20 = (TH1D*) filein->Get("NTD/Data/Accepted/LWO-Vetoed_Spectrum_0-20MeV");
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 20);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 3);
    Print_int_rate(histo_0_10, histo_0_20, 10, 20);
    Print_int_rate(histo_0_10, histo_0_20, 3, 10);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 0.5);
    
    cout<<"\n=============== Accepted Data with Acc =================== "<<endl;
    histo_0_10 = (TH1D*) filein->Get("NTD/Data/Accepted/LWO-Vetoed_withAcc_Norm_Spectrum_0-10MeV");
    histo_0_20 = (TH1D*) filein->Get("NTD/Data/Accepted/LWO-Vetoed_withAcc_Norm_Spectrum_0-20MeV");
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 20);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 3);
    Print_int_rate(histo_0_10, histo_0_20, 10, 20);
    Print_int_rate(histo_0_10, histo_0_20, 3, 10);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 0.5);
    
    cout<<"\n=============== Triggered Simu =================== "<<endl;
    histo_0_10 = (TH1D*) filein->Get("NTD/Sum/Triggered/histo_simu_sumNTD_0-20MeV");
    histo_0_20 = (TH1D*) filein->Get("NTD/Sum/Triggered/histo_simu_sumNTD_0-10MeV");
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 20);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 3);
    Print_int_rate(histo_0_10, histo_0_20, 10, 20);
    Print_int_rate(histo_0_10, histo_0_20, 3, 10);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 0.5);
    cout<<"\n=============== Accepted Simu =================== "<<endl;
    histo_0_10 = (TH1D*) filein->Get("NTD/Sum/coinc_Mid_only/histo_simu_sumNTD_0-20MeV");
    histo_0_20 = (TH1D*) filein->Get("NTD/Sum/coinc_Mid_only/histo_simu_sumNTD_0-20MeV");
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 20);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 3);
    Print_int_rate(histo_0_10, histo_0_20, 10, 20);
    Print_int_rate(histo_0_10, histo_0_20, 3, 10);
    Print_int_rate(histo_0_10, histo_0_20, 0.05, 0.5);
}
