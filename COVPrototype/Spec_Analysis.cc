//Code to analyse the spectra from choozerent analyzer >> dru.
#include "Recalib.C"
#include "Shift.C"

double FLUX = 0.0178; // events/cm2/s neutrons : 0.0134 ev/cm2/s - muons : 0.019 ev/cm2/s - AmbientGammas : 1.25e6 ev/cm2/d // 3 ev/cm^2/s ?
double PLANE_SURF = 140*140;//10*10; //in cm^2 //1 for Am (-> Rate in ev/s)
double RATE = FLUX*PLANE_SURF; //events/s

double ENERGY_UNIT = 1e3; //energy unit with respect to keV (GeV =Rate_above_5MeV 1e6 Mev = 1e3, keV = 1, eV = 1e-3,...)

double MASS_Ge = 400*1e-3;
double MASS_LWO = 51.7*1e-3;

int NEW_NBINS = 500;

double THRESH_TOP = 10; //keV
double THRESH_BOT = 4; //keV
double THRESH_LWO = 25; //keV

double NUM_EVENTS;
double NUM_E_ANALYZED=2; //1 if no Evis (no resolution applied in the analysis)

double HQF = 1.095; //Heat quenching factor =/= 1 only for neutrons

void help(){
    cout<<"Analysis of EnerGe simulation spectra:"<<endl;
    cout<< ""<< endl;
    cout<<"================ CONSTANT VALUES ============="<<endl;
    cout<<"Flux = " <<FLUX<<" events/cm2/s"<<endl;
    cout<<"Surface of plane generating the events = " <<PLANE_SURF<< " cm^2"<<endl;
    cout<<"Rate = "<<RATE<< " events/d"<<endl;
    cout<<"Default energy unit = "<< ENERGY_UNIT << " keV" << endl;
    cout<<"Mass of Ge detector = "<< MASS_Ge<< " kg"<<endl;
    cout<<"Mass of LWO detector = "<< MASS_LWO<< " kg"<<endl;
    cout<<"Num energies analyzed = "<<NUM_E_ANALYZED<< " (if 1: only Edep, if 2: Edep and Evis)"<<endl;
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
    cout<<"void Hits_to_cps_data(TH1D* histo, double energy_unit = ENERGY_UNIT, double detect_mass = MASS_Ge, int RunTime=1)"<<endl;
    cout<<"\t normalize the histogram histo (from Eddy MyCOV script data) in count per seconds (wiith the RunTime as input)"<<endl;
    
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
    
}

void save_CPS_histos(TFile *filein, TFile *fileout, TString Crystal, TString hName, double Crystal_Mass, double thresh_crys, double thresh_crys_fixed, double NumEv, double Rate);

double SetNumEvents(TString FileName = "Analyzer_EP_Nu1g_AtmNSphere_Al2O3_10Mevents.root", bool Am = FALSE){
    
    cout<<"SetNumEvents"<<endl;
    TFile *file = new TFile(FileName);
    TH1D *histoPrim = (TH1D*) file->Get("PrimariesAlgoDir/Energy/hEinj_log");
    NUM_EVENTS=histoPrim->GetEntries();
    
    if (Am){
        //TH1D *histo = (TH1D*) file->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-100keV");
        TH1D *histo = (TH1D*) file->Get("OrsayOVAlgoDirDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-100keV");
        NUM_EVENTS=histo->Integral(histo->FindBin(59.3), histo->FindBin(59.8));
    }
    cout<<"end of SetNumEvents"<<endl;
    return NUM_EVENTS;
}

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
    histo->Draw();
}

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

vector<double> Int_Rate(TH1D* histo, vector<double>& Int_Rate, double Estart=0, double Eend=0, double MassTarget = MASS_Ge, double num_events = NUM_EVENTS, double rate = RATE, bool Verbose_level=1){
    
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
    
    return Int_Rate;
}

void Norm_in_cps(TString Filename, TString NameFileout, bool Am = FALSE) {
    
    TFile *Fileout = new TFile(NameFileout, "RECREATE");
    Fileout->Close();
    
    TCanvas *c1 = new TCanvas();
    
    TFile *filein = new TFile(Filename);
    int Hits_i, Error_i, NumBin;
    
    double NumEv = SetNumEvents(Filename, Am);
    if (Am){
        cout << "Num events = " <<NumEv << " in 60keV peak -- expected number of events in 60keV peak: " << FLUX << " ev/s" << endl;
    }
    else {cout << "Num events = " <<NumEv << " -- Flux : " << FLUX << " ev/cm2/s" << endl;
    }
    cout << " " <<endl;
    
    //Top Ge
    Fileout = new TFile(NameFileout, "UPDATE");
    Fileout->mkdir("Top_Ge");
    Fileout->cd("Top_Ge");
    
    vector<double> Rate_thresh, Rate_fixed;
    cout << "----------------------------------------------------" <<endl;
    cout << "Top Ge - Target mass :" << MASS_Ge << " kg" <<endl;
    
    TH1D *h_Ge_Top = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Edep/nEdep_triggered_Vol_Top_Ge_0-100keV");
    h_Ge_Top->SetTitle("Top_Ge");
    Int_Rate(h_Ge_Top, Rate_thresh, THRESH_TOP, 100, MASS_Ge, NumEv, RATE);
    Int_Rate(h_Ge_Top, Rate_fixed, 30, 100, MASS_Ge, NumEv, RATE, 1);
    
    Hits_to_cps(h_Ge_Top, 1, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Top->Write("Top_Ge_0-100keV", TObject::kOverwrite);
    h_Ge_Top->Delete();
    
    h_Ge_Top = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Edep/nEdep_triggered_Vol_Top_Ge_0-1MeV");
    h_Ge_Top->SetTitle("Top_Ge");
    Int_Rate(h_Ge_Top, Rate_thresh, 0.1, 1, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Top, Rate_fixed, 0.1, 1, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_Ge_Top, 1e3, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Top->Write("Top_Ge_0-1MeV", TObject::kOverwrite);
    h_Ge_Top->Delete();
    
    h_Ge_Top = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Edep/nEdep_triggered_Vol_Top_Ge_0-10MeV");
    h_Ge_Top->SetTitle("Top_Ge");
    Int_Rate(h_Ge_Top, Rate_thresh, 1, 10, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Top, Rate_fixed, 1, 10, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_Ge_Top, 1e3, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Top->Write("Top_Ge_0-10MeV", TObject::kOverwrite);
    h_Ge_Top->Delete();
    
    h_Ge_Top = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Edep/nEdep_triggered_Vol_Top_Ge_0-10MeV");
    NumBin = h_Ge_Top->FindBin(4)-h_Ge_Top->FindBin(0);
    TH1D* new_h_Ge_Top = new TH1D("new_h_Ge_Top", "new_h_Ge_Top", NumBin, 0, 4);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_Ge_Top->GetBinContent(i);
        Error_i = h_Ge_Top->GetBinError(i);
        
        new_h_Ge_Top->SetBinContent(i, Hits_i);
        new_h_Ge_Top->SetBinError(i, Error_i);
    }
    
    new_h_Ge_Top->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_Ge_Top->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_Ge_Top, 1e3, MASS_Ge, NumEv, RATE, 200);
    new_h_Ge_Top->Write("Top_Ge_0-4MeV");
    
    h_Ge_Top->Delete();
    new_h_Ge_Top->Delete();
    
    h_Ge_Top = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Edep/nEdep_triggered_Vol_Top_Ge_0-100MeV");
    h_Ge_Top->SetTitle("Top_Ge");
    Int_Rate(h_Ge_Top, Rate_thresh, 10, 20, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Top, Rate_fixed, 10, 20, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_Ge_Top, 1e3, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Top->Write("Top_Ge_0-100MeV", TObject::kOverwrite);
    
    h_Ge_Top->Delete();
    
    h_Ge_Top = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Edep/nEdep_triggered_Vol_Top_Ge_0-100MeV");
    NumBin = h_Ge_Top->FindBin(20)-h_Ge_Top->FindBin(0);
    new_h_Ge_Top = new TH1D("new_h_Ge_Top", "new_h_Ge_Top", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_Ge_Top->GetBinContent(i);
        Error_i = h_Ge_Top->GetBinError(i);
        
        new_h_Ge_Top->SetBinContent(i, Hits_i);
        new_h_Ge_Top->SetBinError(i, Error_i);
    }
    
    new_h_Ge_Top->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_Ge_Top->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_Ge_Top, 1e3, MASS_Ge, NumEv, RATE, 100);
    new_h_Ge_Top->Write("Top_Ge_0-20MeV");
    new_h_Ge_Top->Delete();
    
//    cout<<""<<endl;
//    vector<double> Rate_above_5MeV;
//    Int_Rate(h_Ge_Top, Rate_above_5MeV, 5, 20, MASS_Ge, NumEv, RATE, 0);
//    cout<<" 5 MeV - 20 MeV --- "<< Rate_above_5MeV[0] << " +-  "<< Rate_above_5MeV[1] << "ev/s ---- " << Rate_above_5MeV[0]/MASS_Ge << " +- " << Rate_above_5MeV[1]/MASS_Ge <<" ev/s/kg" << endl;
//    cout<<""<<endl;
    
    h_Ge_Top->Delete();
    Fileout->Close();
    
    double Tot_Rate_thresh=0;
    double Error_Rate_thresh_sq=0;
    double Tot_Rate_fixed=0;
    double Error_Rate_fixed_sq=0;
    
    for (int i=0; i<Rate_thresh.size()/2; i++){
        Tot_Rate_thresh+=Rate_thresh[2*i];
        Error_Rate_thresh_sq+=Rate_thresh[2*i+1]*Rate_thresh[2*i+1];
        Tot_Rate_fixed+=Rate_fixed[2*i];
        Error_Rate_fixed_sq+=Rate_fixed[2*i+1]*Rate_thresh[2*i+1];
    }
    double Error_Rate_thresh=sqrt(Error_Rate_thresh_sq);
    double Error_Rate_fixed=sqrt(Error_Rate_fixed_sq);
    
    
    cout<<"Total rates :"<<endl;
    cout<<THRESH_TOP << " keV - 20 MeV --- "<< Tot_Rate_thresh << " +-  "<< Error_Rate_thresh << "ev/s ---- " << Tot_Rate_thresh/MASS_Ge << " +- " << Error_Rate_thresh/MASS_Ge <<" ev/s/kg" << endl;
    cout<<"30 keV - 20 MeV --- "<< Tot_Rate_fixed << " +-  "<< Error_Rate_fixed << "ev/s ---- " << Tot_Rate_fixed/MASS_Ge << " +- " << Error_Rate_fixed/MASS_Ge <<" ev/s/kg" << endl;
    
    //Bot Ge
    Fileout = new TFile(NameFileout, "UPDATE");
    Rate_thresh.clear();
    Rate_fixed.clear();
    
    Fileout->mkdir("Bot_Ge");
    Fileout->cd("Bot_Ge");
    
    cout << "----------------------------------------------------" <<endl;
    cout << "Bot Ge - Target mass :" << MASS_Ge << " kg -- Flux : " << FLUX << " ev/cm2/s" << endl;
    
    TH1D *h_Ge_Bot = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-100keV");
    h_Ge_Bot->SetTitle("Bot_Ge");
    Int_Rate(h_Ge_Bot, Rate_thresh, THRESH_BOT, 100, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Bot, Rate_fixed, 30, 100, MASS_Ge, NumEv, RATE, 1);
    
    Hits_to_cps(h_Ge_Bot, 1, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Bot->Write("Bot_Ge_0-100keV", TObject::kOverwrite);
    h_Ge_Bot->Delete();
    
    h_Ge_Bot = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-1MeV");
    h_Ge_Bot->SetTitle("Bot_Ge");
    Int_Rate(h_Ge_Bot, Rate_thresh, 0.1, 1, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Bot, Rate_fixed, 0.1, 1, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_Ge_Bot, 1e3, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Bot->Write("Bot_Ge_0-1MeV", TObject::kOverwrite);
    h_Ge_Bot->Delete();
    
    h_Ge_Bot = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-10MeV");
    h_Ge_Bot->SetTitle("Bot_Ge");
    Int_Rate(h_Ge_Bot, Rate_thresh, 1, 10, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Bot, Rate_fixed, 1, 10, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_Ge_Bot, 1e3, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Bot->Write("Bot_Ge_0-10MeV", TObject::kOverwrite);
    
    h_Ge_Bot->Delete();
    
    h_Ge_Bot = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-10MeV");
    NumBin = h_Ge_Bot->FindBin(4)-h_Ge_Bot->FindBin(0);
    TH1D* new_h_Ge_Bot = new TH1D("new_h_Ge_Bot", "new_h_Ge_Bot", NumBin, 0, 4);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_Ge_Bot->GetBinContent(i);
        Error_i = h_Ge_Bot->GetBinError(i);
        
        new_h_Ge_Bot->SetBinContent(i, Hits_i);
        new_h_Ge_Bot->SetBinError(i, Error_i);
    }
    
    new_h_Ge_Bot->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_Ge_Bot->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_Ge_Bot, 1e3, MASS_Ge, NumEv, RATE, 200);
    new_h_Ge_Bot->Write("Bot_Ge_0-4MeV");
    
    h_Ge_Bot->Delete();
    new_h_Ge_Bot->Delete();
    
    h_Ge_Bot = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-100MeV");
    h_Ge_Bot->SetTitle("Bot_Ge");
    Int_Rate(h_Ge_Bot, Rate_thresh, 10, 20, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_Ge_Bot, Rate_fixed, 10, 20, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_Ge_Bot, 1e3, MASS_Ge, NumEv, RATE, 500);
    h_Ge_Bot->Write("Bot_Ge_0-100MeV", TObject::kOverwrite);
    
    h_Ge_Bot->Delete();
    
    h_Ge_Bot = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Edep/nEdep_triggered_Vol_Bot_Ge_0-100MeV");
    NumBin = h_Ge_Bot->FindBin(20)-h_Ge_Bot->FindBin(0);
    new_h_Ge_Bot = new TH1D("new_h_Ge_Bot", "new_h_Ge_Bot", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_Ge_Bot->GetBinContent(i);
        Error_i = h_Ge_Bot->GetBinError(i);
        
        new_h_Ge_Bot->SetBinContent(i, Hits_i);
        new_h_Ge_Bot->SetBinError(i, Error_i);
    }
    
    new_h_Ge_Bot->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_Ge_Bot->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_Ge_Bot, 1e3, MASS_Ge, NumEv, RATE, 100);
    new_h_Ge_Bot->Write("Bot_Ge_0-20MeV");
    
    h_Ge_Bot->Delete();
    new_h_Ge_Bot->Delete();
    Fileout->Close();
    
    Tot_Rate_thresh=0;
    Error_Rate_thresh_sq=0;
    Tot_Rate_fixed=0;
    Error_Rate_fixed_sq=0;
    for (int i=0;i<Rate_thresh.size()/2; i++){
        Tot_Rate_thresh+=Rate_thresh[2*i];
        Error_Rate_thresh_sq+=Rate_thresh[2*i+1]*Rate_thresh[2*i+1];
        Tot_Rate_fixed+=Rate_fixed[2*i];
        Error_Rate_fixed_sq+=Rate_fixed[2*i+1]*Rate_fixed[2*i+1];
    }
    Error_Rate_thresh=sqrt(Error_Rate_thresh_sq);
    Error_Rate_fixed=sqrt(Error_Rate_fixed_sq);
    
    cout<<"Total rates :"<<endl;
    cout<<THRESH_BOT << " keV - 20 MeV --- "<< Tot_Rate_thresh << " +-  "<< Error_Rate_thresh << "ev/s ---- " << Tot_Rate_thresh/MASS_Ge << " +- " << Error_Rate_thresh/MASS_Ge <<" ev/s/kg" << endl;
    cout<<"30 keV - 20 MeV --- "<< Tot_Rate_fixed << " +-  "<< Error_Rate_fixed << "ev/s ---- " << Tot_Rate_fixed/MASS_Ge << " +- " << Error_Rate_fixed/MASS_Ge <<" ev/s/kg" << endl;
    
    //LWO
    Rate_thresh.clear();
    Rate_fixed.clear();
    
    Fileout = new TFile(NameFileout, "UPDATE");
    Fileout->mkdir("LWO");
    Fileout->cd("LWO");
    
    cout << "----------------------------------------------------" <<endl;
    cout << "LWO - Target mass :" << MASS_LWO << " kg -- Flux : " << FLUX << " ev/cm2/s" << endl;
    
    TH1D *h_LWO = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Edep/nEdep_triggered_Vol_Mid_LWO_0-100keV");
    h_LWO->SetTitle("Bot_Ge");
    Int_Rate(h_LWO, Rate_thresh, THRESH_LWO, 100, MASS_LWO, NumEv, RATE, 1);
    Int_Rate(h_LWO, Rate_fixed, 50, 100, MASS_LWO, NumEv, RATE, 1);
    
    Hits_to_cps(h_LWO, 1, MASS_LWO, NumEv, RATE, 500);
    h_LWO->Write("LWO_0-100keV", TObject::kOverwrite);
    h_LWO->Delete();
    
    h_LWO = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Edep/nEdep_triggered_Vol_Mid_LWO_0-1MeV");
    h_LWO->SetTitle("LWO");
    Int_Rate(h_LWO, Rate_thresh, 0.1, 1, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_LWO, Rate_fixed, 0.1, 1, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_LWO, 1e3, MASS_LWO, NumEv, RATE, 500);
    h_LWO->Write("LWO_0-1MeV", TObject::kOverwrite);
    h_LWO->Delete();
    
    h_LWO = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Edep/nEdep_triggered_Vol_Mid_LWO_0-10MeV");
    h_LWO->SetTitle("LWO");
    
    Int_Rate(h_LWO, Rate_thresh, 1, 10, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_LWO, Rate_fixed, 1, 10, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_LWO, 1e3, MASS_LWO, NumEv, RATE, 500);
    h_LWO->Write("LWO_0-10MeV", TObject::kOverwrite);
    h_LWO->Delete();
    
    h_LWO = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Edep/nEdep_triggered_Vol_Mid_LWO_0-10MeV");
    NumBin = h_LWO->FindBin(4)-h_LWO->FindBin(0);
    TH1D* new_h_LWO = new TH1D("new_h_LWO", "new_h_LWO", NumBin, 0, 4);
    
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_LWO->GetBinContent(i);
        Error_i = h_LWO->GetBinError(i);
        
        new_h_LWO->SetBinContent(i, Hits_i);
        new_h_LWO->SetBinError(i, Error_i);
    }
    
    new_h_LWO->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_LWO->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_LWO, 1e3, MASS_LWO, NumEv, RATE, 200);
    new_h_LWO->Write("LWO_0-4MeV");
    
    cout<<""<<endl;
    vector<double> PersonalRangeRate;
    Int_Rate(h_LWO, PersonalRangeRate, 4.5, 5.5, MASS_LWO, NumEv, RATE, 0);
    cout<<"4.5 MeV - 5.5 MeV --- "<< PersonalRangeRate[0] << " +-  "<< PersonalRangeRate[1] << "ev/s ---- " << PersonalRangeRate[0]/MASS_LWO << " +- " << PersonalRangeRate[1]/MASS_LWO <<" ev/s/kg" << endl;
    cout<<""<<endl;
    
    h_LWO->Delete();
    new_h_LWO->Delete();
    
    h_LWO = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Edep/nEdep_triggered_Vol_Mid_LWO_0-100MeV");
    h_LWO->SetTitle("LWO");
    Int_Rate(h_LWO, Rate_thresh, 10, 20, MASS_Ge, NumEv, RATE, 1);
    Int_Rate(h_LWO, Rate_fixed, 10, 20, MASS_Ge, NumEv, RATE, 0);
    
    Hits_to_cps(h_LWO, 1e3, MASS_LWO, NumEv, RATE, 500);
    h_LWO->Write("LWO_0-100MeV", TObject::kOverwrite);
    h_LWO->Delete();
    
    h_LWO = (TH1D*)  filein->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Edep/nEdep_triggered_Vol_Mid_LWO_0-100MeV");
    NumBin = h_LWO->FindBin(20)-h_LWO->FindBin(0);
    new_h_LWO = new TH1D("new_h_LWO", "new_h_LWO", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_LWO->GetBinContent(i);
        Error_i = h_LWO->GetBinError(i);
        
        new_h_LWO->SetBinContent(i, Hits_i);
        new_h_LWO->SetBinError(i, Error_i);
    }
    
    new_h_LWO->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_LWO->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps(new_h_LWO, 1e3, MASS_LWO, NumEv, RATE, 100);
    new_h_LWO->Write("LWO_0-20MeV");
    
    h_LWO->Delete();
    new_h_LWO->Delete();
    
    Tot_Rate_thresh=0;
    Error_Rate_thresh_sq=0;
    Tot_Rate_fixed=0;
    Error_Rate_fixed_sq=0;
    for (int i=0; i<Rate_thresh.size()/2; i++){
        Tot_Rate_thresh+=Rate_thresh[2*i];
        Error_Rate_thresh_sq+=Rate_thresh[2*i+1]*Rate_thresh[2*i+1];
        Tot_Rate_fixed+=Rate_fixed[2*i];
        Error_Rate_fixed_sq+=Rate_fixed[2*i+1]*Rate_thresh[2*i+1];
    }
    Error_Rate_thresh=sqrt(Error_Rate_thresh_sq);
    Error_Rate_fixed=sqrt(Error_Rate_fixed_sq);
    
    cout<<"Total rates :"<<endl;
    cout<<THRESH_LWO << " keV - 20 MeV --- "<< Tot_Rate_thresh << " +-  "<< Error_Rate_thresh << " ev/s ---- " << Tot_Rate_thresh/MASS_Ge << " +- " << Error_Rate_thresh/MASS_Ge <<" ev/s/kg" << endl;
    cout<<"50 keV - 20 MeV --- "<< Tot_Rate_fixed << " +-  "<< Error_Rate_fixed << " ev/s ---- " << Tot_Rate_fixed/MASS_Ge << " +- " << Error_Rate_fixed/MASS_Ge <<" ev/s/kg" << endl;
    
    Fileout->Close();
    c1->Close();
}

void form_data(TString FileinName, TString FileoutName){
    TCanvas *c1 = new TCanvas();
    TFile *Filein = new TFile(FileinName, "READ");
    TFile *Fileout = new TFile(FileoutName, "RECREATE");
    
    TH1D *histo_RunTime = (TH1D*) Filein->Get("h1_RunTime");
    int RunTime=histo_RunTime->GetBinContent(1);
    
    int Hits_i, Error_i, NumBin;
    
    //COV1
    TH1D *h_COV1 = (TH1D*) Filein->Get("h1_Peak1_4MeV");
    h_COV1->GetXaxis()->SetTitle("Energy [MeV]");
    h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    h_COV1->Write("COV1_Raw_Spectrum_0-4MeV");
    
    NumBin = h_COV1->FindBin(4)-h_COV1->FindBin(0);
    TH1D *new_h_COV1 = new TH1D("new_h_COV1", "new_h_COV1", NumBin, 0, 4);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV1->GetBinContent(i);
        Error_i = h_COV1->GetBinError(i);
        
        new_h_COV1->SetBinContent(i, Hits_i);
        new_h_COV1->SetBinError(i, Error_i);
    }
    
    new_h_COV1->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV1, 1e3, MASS_Ge, RunTime, 200);
    
    new_h_COV1->Write("COV1_Cal_Spectrum_0-4MeV");
    new_h_COV1->Delete();
    
    NumBin = h_COV1->FindBin(1)-h_COV1->FindBin(0);
    new_h_COV1 = new TH1D("new_h_COV1", "new_h_COV1", NumBin, 0, 1);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV1->GetBinContent(i);
        Error_i = h_COV1->GetBinError(i);
        
        new_h_COV1->SetBinContent(i, Hits_i);
        new_h_COV1->SetBinError(i, Error_i);
    }
    
    h_COV1->Delete();
    new_h_COV1->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV1, 1e3, MASS_Ge, RunTime, 500);
    
    new_h_COV1->Write("COV1_Cal_Spectrum_0-1MeV");
    new_h_COV1->Delete();
    
    h_COV1 = (TH1D*) Filein->Get("h1_Peak1_MeV");
    h_COV1->GetXaxis()->SetTitle("Energy [MeV]");
    h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    h_COV1->Write("COV1_Raw_Spectrum_0-20MeV");
    
    NumBin = h_COV1->FindBin(20)-h_COV1->FindBin(0);
    new_h_COV1 = new TH1D("new_h_COV1", "new_h_COV1", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV1->GetBinContent(i);
        Error_i = h_COV1->GetBinError(i);
        
        new_h_COV1->SetBinContent(i, Hits_i);
        new_h_COV1->SetBinError(i, Error_i);
    }
    
    new_h_COV1->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV1, 1e3, MASS_Ge, RunTime, 100);
    
    new_h_COV1->Write("COV1_Cal_Spectrum_0-20MeV");
    new_h_COV1->Delete();
    
    NumBin = h_COV1->FindBin(10)-h_COV1->FindBin(0);
    new_h_COV1 = new TH1D("new_h_COV1", "new_h_COV1", NumBin, 0, 10);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV1->GetBinContent(i);
        Error_i = h_COV1->GetBinError(i);
        
        new_h_COV1->SetBinContent(i, Hits_i);
        new_h_COV1->SetBinError(i, Error_i);
    }
    
    h_COV1->Delete();
    new_h_COV1->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV1->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV1, 1e3, MASS_Ge, RunTime, 1000);
    new_h_COV1->Write("COV1_Cal_Spectrum_0-10MeV_1000bins");
    
    new_h_COV1->Rebin(2);
    new_h_COV1->Scale(1/2.);
    new_h_COV1->Write("COV1_Cal_Spectrum_0-10MeV");
    
    new_h_COV1->Delete();
    
    //COV2
    TH1D *h_COV2 = (TH1D*) Filein->Get("h1_Peak2_4MeV");
    h_COV2->GetXaxis()->SetTitle("Energy [MeV]");
    h_COV2->GetYaxis()->SetTitle("Counts [a.u.]");
    h_COV2->Write("COV2_Raw_Spectrum_0-4MeV");
    
    NumBin = h_COV2->FindBin(4)-h_COV2->FindBin(0);
    TH1D *new_h_COV2 = new TH1D("new_h_COV2", "new_h_COV2", NumBin, 0, 4);
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV2->GetBinContent(i);
        Error_i = h_COV2->GetBinError(i);
        
        new_h_COV2->SetBinContent(i, Hits_i);
        new_h_COV2->SetBinError(i, Error_i);
    }
    
    new_h_COV2->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV2->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV2, 1e3, MASS_Ge, RunTime, 200);
    new_h_COV2->Write("COV2_Cal_Spectrum_0-4MeV");
    
    new_h_COV2->Rebin(2);
    new_h_COV2->Scale(1./2);
    
    new_h_COV2->Delete();
    
    NumBin = h_COV2->FindBin(1)-h_COV2->FindBin(0);
    new_h_COV2 = new TH1D("new_h_COV2", "new_h_COV2", NumBin, 0, 1);
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV2->GetBinContent(i);
        Error_i = h_COV2->GetBinError(i);
        
        new_h_COV2->SetBinContent(i, Hits_i);
        new_h_COV2->SetBinError(i, Error_i);
    }
    
    h_COV2->Delete();
    new_h_COV2->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV2->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV2, 1e3, MASS_Ge, RunTime, 500);
    new_h_COV2->Write("COV2_Cal_Spectrum_0-1MeV");
    new_h_COV2->Delete();
    
    h_COV2 = (TH1D*) Filein->Get("h1_Peak2_MeV");
    h_COV2->GetXaxis()->SetTitle("Energy [MeV]");
    h_COV2->GetYaxis()->SetTitle("Counts [a.u.]");
    h_COV2->Write("COV2_Raw_Spectrum_0-20MeV");
    
    NumBin = h_COV2->FindBin(20)-h_COV2->FindBin(0);
    new_h_COV2 = new TH1D("new_h_COV2", "new_h_COV2", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV2->GetBinContent(i);
        Error_i = h_COV2->GetBinError(i);
        
        new_h_COV2->SetBinContent(i, Hits_i);
        new_h_COV2->SetBinError(i, Error_i);
    }
    
    new_h_COV2->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV2->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV2, 1e3, MASS_Ge, RunTime, 100);
    new_h_COV2->Write("COV2_Cal_Spectrum_0-20MeV");
    new_h_COV2->Delete();
    
    NumBin = h_COV2->FindBin(10)-h_COV2->FindBin(0);
    new_h_COV2 = new TH1D("new_h_COV2", "new_h_COV2", NumBin, 0, 10);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_COV2->GetBinContent(i);
        Error_i = h_COV2->GetBinError(i);
        
        new_h_COV2->SetBinContent(i, Hits_i);
        new_h_COV2->SetBinError(i, Error_i);
    }
    
    h_COV2->Delete();
    new_h_COV2->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_COV2->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_COV2, 1e3, MASS_Ge, RunTime, 500);
    
    new_h_COV2->Write("COV2_Cal_Spectrum_0-10MeV");
    new_h_COV2->Delete();
    
    //NTD
    TH1D *h_NTD = (TH1D*) Filein->Get("h1_NTD_1_MeV");
    h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    h_NTD->Write("NTD_Raw_Spectrum_0-20MeV");
    
    NumBin = h_NTD->FindBin(20)-h_NTD->FindBin(0);
    TH1D *new_h_NTD = new TH1D("new_h_NTD", "new_h_NTD", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_NTD->GetBinContent(i);
        Error_i = h_NTD->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_NTD, 1e3, MASS_LWO, RunTime, 100);
    new_h_NTD->Write("NTD_Cal_Spectrum_0-20MeV");
    new_h_NTD->Delete();
    
    NumBin = h_NTD->FindBin(10)-h_NTD->FindBin(0);
    new_h_NTD = new TH1D("new_h_NTD", "new_h_NTD", NumBin, 0, 10);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_NTD->GetBinContent(i);
        Error_i = h_NTD->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_NTD, 1e3, MASS_LWO, RunTime, 1000);
    new_h_NTD->Write("NTD_Cal_Spectrum_0-10MeV_1000bins");
    
    new_h_NTD->Rebin(2);
    new_h_NTD->Scale(1/2.);
    
    new_h_NTD->Write("NTD_Cal_Spectrum_0-10MeV");
    new_h_NTD->Delete();
    
    NumBin = h_NTD->FindBin(1)-h_NTD->FindBin(0);
    new_h_NTD = new TH1D("new_h_NTD", "new_h_NTD", NumBin, 0, 1);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h_NTD->GetBinContent(i);
        Error_i = h_NTD->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    //h_NTD->Delete();
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    //Hits_to_cps_data(new_h_NTD, 1e3, MASS_Ge, RunTime, 200);
    //new_h_NTD->Write("NTD_Cal_Spectrum_0-1MeV");

    //Veto
    
//    TH1D *h1_NTD_Veto_MeV = (TH1D*) h_NTD->Clone();
//    TH1D *h1_NTD_2_MeV = (TH1D*) Filein->Get("h1_NTD_2_MeV");
//    TH1D *h1_NTD_3_MeV = (TH1D*) Filein->Get("h1_NTD_3_MeV");
//    TH1D *h1_NTD_4_MeV = (TH1D*) Filein->Get("h1_NTD_4_MeV");
//    TH1D *h1_NTD_5_MeV = (TH1D*) Filein->Get("h1_NTD_5_MeV");
//
//    h1_NTD_Veto_MeV->SetTitle("h_NTD_Veto");
//
//    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_2_MeV, 0., 1.);
//    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_3_MeV, 1., .5);
//    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_4_MeV, 1., .5);
//    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_5_MeV, 1.,-.5);
    
    TH1D *h1_NTD_Veto_MeV = (TH1D*) h_NTD->Clone();
    TH1D *h1_NTD_2_MeV = (TH1D*) Filein->Get("h1_NTD_2_MeV");
    TH1D *h1_NTD_7_MeV = (TH1D*) Filein->Get("h1_NTD_7_MeV");
    TH1D *h1_NTD_8_MeV = (TH1D*) Filein->Get("h1_NTD_8_MeV");
    TH1D *h1_NTD_10_MeV = (TH1D*) Filein->Get("h1_NTD_10_MeV");
    
    h1_NTD_Veto_MeV->SetTitle("h_NTD_Veto");
    
    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_2_MeV, 0., 1.);// h1_NTD_2_MeV = all events veto
    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_7_MeV, 1., .5);
    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_8_MeV, 1., .5);
    h1_NTD_Veto_MeV->Add(h1_NTD_Veto_MeV,h1_NTD_10_MeV, 1., 1.);
    
    h1_NTD_Veto_MeV->Write("NTD-Vetoed_Raw_Spectrum_0-20MeV");
    
    NumBin = h1_NTD_Veto_MeV->FindBin(20)-h1_NTD_Veto_MeV->FindBin(0);
    new_h_NTD = new TH1D("new_h_NTD_veto", "new_h_NTD_veto", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h1_NTD_Veto_MeV->GetBinContent(i);
        Error_i = h1_NTD_Veto_MeV->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_NTD, 1e3, MASS_LWO, RunTime, 100);
    new_h_NTD->Write("NTD-Vetoed_Cal_Spectrum_0-20MeV");
    new_h_NTD->Delete();
    
    NumBin = h1_NTD_Veto_MeV->FindBin(10)-h1_NTD_Veto_MeV->FindBin(0);
    new_h_NTD = new TH1D("new_h_NTD", "new_h_NTD", NumBin, 0, 10);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h1_NTD_Veto_MeV->GetBinContent(i);
        Error_i = h1_NTD_Veto_MeV->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_NTD, 1e3, MASS_LWO, RunTime, 500);
    new_h_NTD->Write("NTD-Vetoed_Cal_Spectrum_0-10MeV");
    new_h_NTD->Delete();

    NumBin = h1_NTD_2_MeV->FindBin(20)-h1_NTD_2_MeV->FindBin(0);
    new_h_NTD = new TH1D("new_h_NTD_veto_withAcc", "new_h_NTD_veto_withAcc", NumBin, 0, 20);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h1_NTD_2_MeV->GetBinContent(i);
        Error_i = h1_NTD_2_MeV->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_NTD, 1e3, MASS_LWO, RunTime, 100);
    new_h_NTD->Write("NTD-VetoedWithAcc_Cal_Spectrum_0-20MeV");
    new_h_NTD->Delete();
    
    NumBin = h1_NTD_2_MeV->FindBin(10)-h1_NTD_2_MeV->FindBin(0);
    new_h_NTD = new TH1D("new_h_NTD_veto_withAcc", "new_h_NTD_veto_withAcc", NumBin, 0, 10);
    
    for (int i=1; i<=NumBin; i++){
        Hits_i = h1_NTD_2_MeV->GetBinContent(i);
        Error_i = h1_NTD_2_MeV->GetBinError(i);
        
        new_h_NTD->SetBinContent(i, Hits_i);
        new_h_NTD->SetBinError(i, Error_i);
    }
    
    new_h_NTD->GetXaxis()->SetTitle("Energy [MeV]");
    new_h_NTD->GetYaxis()->SetTitle("Counts [a.u.]");
    
    Hits_to_cps_data(new_h_NTD, 1e3, MASS_LWO, RunTime, 100);
    new_h_NTD->Write("NTD-VetoedWithAcc_Cal_Spectrum_0-10MeV");
    new_h_NTD->Delete();
    
    Fileout->Close();
    Filein->Close();
    c1->Close();
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
                        "OrsayOVAlgoDirDir/Accepted/Vol_Top_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_accepted_Vol_Top_Ge",
                        MASS_Ge, THRESH_TOP, 30, NumEv, RATE);

        //LWO
        Fileout->mkdir("Accepted/Mid_LWO/"+EFolder[i]);
        Fileout->cd("Accepted/Mid_LWO/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Mid_LWO",
                        "OrsayOVAlgoDirDir/Accepted/Vol_Mid_LWO/"+EFolder[i]+"/n"+EFolder[i]+"_accepted_Vol_Mid_LWO",
                        MASS_LWO, THRESH_LWO, 50, NumEv, RATE);

        
        //Bot Ge
        Fileout->mkdir("Accepted/Bot_Ge/"+EFolder[i]);
        Fileout->cd("Accepted/Bot_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Bot_Ge",
                        "OrsayOVAlgoDirDir/Accepted/Vol_Bot_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_accepted_Vol_Bot_Ge",
                        MASS_Ge, THRESH_BOT, 30, NumEv, RATE);
        
        cout<<""<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        cout << "                                    Triggered"<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        
        //TopGe
        Fileout->mkdir("Triggered/Top_Ge/"+EFolder[i]);
        Fileout->cd("Triggered/Top_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Top_Ge",
                        "OrsayOVAlgoDirDir/Triggered/Vol_Top_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Top_Ge",
                        MASS_Ge, THRESH_TOP, 30, NumEv, RATE);

        //LWO
        Fileout->mkdir("Triggered/Mid_LWO/"+EFolder[i]);
        Fileout->cd("Triggered/Mid_LWO/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Mid_LWO",
                        "OrsayOVAlgoDirDir/Triggered/Vol_Mid_LWO/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Mid_LWO",
                        MASS_LWO, THRESH_LWO, 50, NumEv, RATE);
        
        //Bot Ge
        Fileout->mkdir("Triggered/Bot_Ge/"+EFolder[i]);
        Fileout->cd("Triggered/Bot_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Bot_Ge",
                        "OrsayOVAlgoDirDir/Triggered/Vol_Bot_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Bot_Ge",
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

void Norm_in_cps_coinc(TString Filename, TString NameFileout, bool Am = FALSE) {
    
    TString TxtFilename = NameFileout;
    
//    std::streambuf *psbuf, *backup;
//    std::ofstream filestr;
//    filestr.open(TxtFilename);
//
//    backup = std::cout.rdbuf();     // back up cout's streambuf
//    psbuf = filestr.rdbuf();        // get file's streambuf
//    std::cout.rdbuf(psbuf);         // assign streambuf to cout
//    
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
        
        //LWO
        Fileout->mkdir("Accepted/Mid_LWO/"+EFolder[i]);
        Fileout->cd("Accepted/Mid_LWO/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Mid_LWO",
                        "OrsayOVAlgoDir/Coincidences/coinc_MidLWO_only/Vol_Mid_LWO/"+EFolder[i]+"/n"+EFolder[i]+"_coinc_MidLWO_only_Vol_Mid_LWO",
                        MASS_LWO, THRESH_LWO, 50, NumEv, RATE);

        cout<<""<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        cout << "                                    Triggered"<<endl;
        cout << "-----------------------------------------------------------------------------------------" <<endl;
        
        //LWO
        Fileout->mkdir("Triggered/Mid_LWO/"+EFolder[i]);
        Fileout->cd("Triggered/Mid_LWO/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Mid_LWO",
                        "OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Mid_LWO", MASS_LWO, THRESH_LWO, 50, NumEv, RATE);
        
        //TopGe
        Fileout->mkdir("Triggered/Top_Ge/"+EFolder[i]);
        Fileout->cd("Triggered/Top_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Top_Ge",
                        "OrsayOVAlgoDir/Triggered/Vol_Top_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Top_Ge",
                        MASS_Ge, THRESH_TOP, 30, NumEv, RATE);
        
        //Bot Ge
        Fileout->mkdir("Triggered/Bot_Ge/"+EFolder[i]);
        Fileout->cd("Triggered/Bot_Ge/"+EFolder[i]);
        
        save_CPS_histos(filein, Fileout,"Bot_Ge",
                        "OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/"+EFolder[i]+"/n"+EFolder[i]+"_triggered_Vol_Bot_Ge",
                        MASS_Ge, THRESH_BOT, 30, NumEv, RATE
                        );
    }

//    std::cout.rdbuf(backup);        // restore cout's original streambuf
//    filestr.close();
//    
//    ifstream f(TxtFilename); cout<<f.rdbuf();
//    f.close();

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
    //if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    
    Int_Rate(h_cryst, Rate_0, 0, 1, Crystal_Mass, NumEv, Rate, 1);
    
    Int_Rate(h_cryst, Rate_0, 0, 100, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, thresh_crys, 100, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_fixed, thresh_crys_fixed, 100, Crystal_Mass, NumEv, Rate, 1);
    
    Hits_to_cps(h_cryst, 1, Crystal_Mass, NumEv, Rate, 500);
    h_cryst->Write(Crystal+"_0-100keV", TObject::kOverwrite);
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-1MeV");
    //if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 0.1, 1, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, 0.1, 1, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, Rate_fixed, 0.1, 1, Crystal_Mass, NumEv, Rate, 0);
    
    Hits_to_cps(h_cryst, 1e3, Crystal_Mass, NumEv, RATE, 500);
    h_cryst->Write(Crystal+"_0-1MeV", TObject::kOverwrite);
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-10MeV");
    //if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 1, 10, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, 1, 10, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, Rate_fixed, 1, 10, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, trash, 5, 5.5, Crystal_Mass, NumEv, Rate, 1);
    
    Hits_to_cps(h_cryst, 1e3, Crystal_Mass, NumEv, Rate, 500);
    h_cryst->Write(Crystal+"_0-10MeV", TObject::kOverwrite);
    
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-10MeV");
    //if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
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
    //if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
    h_cryst->SetTitle(Crystal);
    Int_Rate(h_cryst, Rate_0, 10, 20, Crystal_Mass, NumEv, Rate, 1);
    Int_Rate(h_cryst, Rate_thresh, 10, 20, Crystal_Mass, NumEv, Rate, 0);
    Int_Rate(h_cryst, Rate_fixed, 10, 20, Crystal_Mass, NumEv, Rate, 0);
//    Int_Rate(h_cryst, trash, 4, 20, Crystal_Mass, NumEv, Rate, 1);
    
    Hits_to_cps(h_cryst, 1e3, Crystal_Mass, NumEv, Rate, 500);
    h_cryst->Write(Crystal+"_0-100MeV", TObject::kOverwrite);
    
    h_cryst->Delete();
    
    h_cryst = (TH1D*)  filein->Get(hName+"_0-100MeV");
    //if (Crystal=="Mid_LWO") { h_cryst = Shift(h_cryst, HQF);}
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
