//Double_t Phi=1.523 ;//ev/s/cm-2/sec 0.271 K40, 1.523 U238, 0.618 Th232
Double_t MASS_CaWO4_1 = 0.7575*1e-3;
Double_t MASS_CaWO4_0 = 0.7575*1e-3;
Double_t MASS_Al2O3 = 0.74625*1e-3;
Double_t MASS_Ge = 1045*1e-3;

Double_t* GetNEvents(TH1D* &h, Double_t xmin = 0., Double_t xmax = 0.)
{
    double BinFirst = h->FindBin(xmin);
    double BinLast= h->FindBin(xmax);

    double NEvents = h->Integral(BinFirst, BinLast);
    
    Double_t* int_ = new Double_t[3];
    int_[0] = NEvents;
    //Use TFeldmanCousin CL method if low stat.
    if (NEvents < 30) {
    TFeldmanCousins fc{0.68};
    int_[1] = fc.CalculateLowerLimit(NEvents,0);
    int_[2] = fc.CalculateUpperLimit(NEvents,0);
    }
    else {
    int_[1] = NEvents-std::sqrt(NEvents);
    int_[2] = NEvents+std::sqrt(NEvents);
    }
    return &int_[0];
}

void Int_day(TH1D* histo, double Estart, double Eend, double MassTarget, double num_events, double daily_rate){
    
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
    
    //previous method to display incertainties
    double Integral, Error;
    
    Integral = histo->IntegralAndError(BinFirst, BinLast, Error);
    
    double Int_Rate = (Integral*daily_rate)/num_events;
    double Error_Rate = (Error*daily_rate)/num_events;
    
    //New method with TFeldmanCousin
    Double_t* Nevents = new Double_t[3];
    Double_t* Rates = new Double_t[3];
    
    Nevents = GetNEvents(histo, Estart, Eend);
    Rates[0] = (Nevents[0]*daily_rate)/num_events;
    Rates[1] = (Nevents[1]*daily_rate)/num_events;
    Rates[2] = (Nevents[2]*daily_rate)/num_events;
    
    if (Int_Rate>5) cout<<"\t\t"<< Form("%.2f", Int_Rate) << ";" << Form("%.2f", Error_Rate)<<";";
    //Feldman Cousins
    else if (Int_Rate>0) cout<<"\t\t"<< Int_Rate << ";" << + Rates[1] - Rates[2]<<";";
    else cout<<"\t\t <"<< Rates[2]<<";\t";
}

void Rate_Latex(TH1D* histo, double Estart, double Eend, double MassTarget, double num_events, double daily_rate){
    
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
    
    //previous method to display incertainties
    double Integral, Error;
    
    Integral = histo->IntegralAndError(BinFirst, BinLast, Error);
    
    double Int_Rate = (Integral*daily_rate)/num_events;
    double Error_Rate = (Error*daily_rate)/num_events;
    
    //New method with TFeldmanCousin
    Double_t* Nevents = new Double_t[3];
    Double_t* Rates = new Double_t[3];
    
    Nevents = GetNEvents(histo, Estart, Eend);
    Rates[0] = (Nevents[0]*daily_rate)/num_events;
    Rates[1] = (Nevents[1]*daily_rate)/num_events;
    Rates[2] = (Nevents[2]*daily_rate)/num_events;
    
    if (Int_Rate>3) cout<<"$"<< Form("%.2f", Int_Rate) << "\\pm" << Form("%.2f", Error_Rate)<<"$";
    //Feldman Cousins
    else if (Int_Rate>0) cout<<"$"<< Form("%.2f", Int_Rate) << "^{+" <<Form("%.2f", Rates[2]-Int_Rate)<<"}_{-"<<Form("%.2f", abs(Rates[1]-Int_Rate))<<"}$";
    else cout<<"$"<<"<"<< Form("%.2f", Rates[2])<<"$";
}

void Print_Table_Latex_cryodet(std::vector<TString> Detector_list, std::vector<std::vector<TH1D*>> LowE_TH1Ds, std::vector<std::vector<TH1D*>> HighE_TH1Ds, std::vector<Double_t> mass_targets, unsigned long long numevents, Double_t dailyrate){
    
    std::cout<<"\\begin{table}"<<std::endl;
    std::cout<<"\\begin{tabular}{ccccccc}"<<std::endl;
    std::cout<<"counts/4 weeks";
    for (int i=0; i<Detector_list.size(); i++){
        std::cout<<"&"<<"\\multicolumn{2}{c}{"<<Detector_list[i]<<"}";
    }
    std::cout<<"\\\\"<<std::endl;
    
    for (int i=0; i<Detector_list.size(); i++){
        std::cout<<"&"<<"50-150 eV"<<"&"<<"0.05-5 keV";
    }
    std::cout<<"\\\\"<<std::endl;
    
    std::vector<TString> Cuts_Names={"No cut", "MV cut", "COV cut", "MV + COV cut", "In coinc with COV", "In coinc with MV"};
    for (int i=0; i<Cuts_Names.size(); i++){
        std::cout<<Cuts_Names[i];
        std::cout<<"&";
        Rate_Latex(LowE_TH1Ds[i][0], 50, 150, mass_targets[i], numevents, dailyrate*7*4.);
        std::cout<<"&";
        Rate_Latex(HighE_TH1Ds[i][0], 0.5, 5, mass_targets[i], numevents, dailyrate*7*4.);
        std::cout<<"&";
        Rate_Latex(LowE_TH1Ds[i][1], 50, 150, mass_targets[i], numevents, dailyrate*7*4.);
        std::cout<<"&";
        Rate_Latex(HighE_TH1Ds[i][1], 0.5, 5, mass_targets[i], numevents, dailyrate*7*4.);
        std::cout<<"&";
        Rate_Latex(LowE_TH1Ds[i][2], 50, 150, mass_targets[i], numevents, dailyrate*7*4.);
        std::cout<<"&";
        Rate_Latex(HighE_TH1Ds[i][2], 0.5, 5, mass_targets[i], numevents, dailyrate*7*4.);
        std::cout<<"\\\\"<<std::endl;
    }
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{table}"<<std::endl;
}

void Print_Table_Latex_COV(std::vector<TH1D*> LowE_TH1Ds, std::vector<TH1D*> HighE_TH1Ds, Double_t mass, unsigned long long numevents, Double_t dailyrate){

    std::cout<<"\\begin{table}"<<std::endl;
    std::cout<<"\\begin{tabular}{ccc}"<<std::endl;
    
    std::cout<<"counts/4 weeks &"<<"0.03-4 MeV"<<"&"<<"0.03-80 MeV \\\\";
    std::cout<<endl;
    
    std::cout<<"No cut";
    std::cout<<"&";
    Rate_Latex(LowE_TH1Ds[0], 0.03, 4, mass, numevents, dailyrate*7*4.);
    std::cout<<"&";
    Rate_Latex(HighE_TH1Ds[0], 0.03, 80, mass, numevents, dailyrate*7*4.);
    std::cout<<"\\\\"<<std::endl;
    
    std::cout<<"MV cut";
    std::cout<<"&";
    Rate_Latex(LowE_TH1Ds[1], 0.03, 4, mass, numevents, dailyrate*7*4.);
    std::cout<<"&";
    Rate_Latex(HighE_TH1Ds[1], 0.03, 80, mass, numevents, dailyrate*7*4.);
    std::cout<<"\\\\"<<std::endl;
  
    std::cout<<"In coinc with CaWO4_0";
    std::cout<<"&";
    Rate_Latex(LowE_TH1Ds[2], 0.03, 4, mass, numevents, dailyrate*7*4.);
    std::cout<<"&";
    Rate_Latex(HighE_TH1Ds[2], 0.03, 80, mass, numevents, dailyrate*7*4.);
    std::cout<<"\\\\"<<std::endl;
    
    std::cout<<"In coinc with CaWO4_1";
    std::cout<<"&";
    Rate_Latex(LowE_TH1Ds[3], 0.03, 4, mass, numevents, dailyrate*7*4.);
    std::cout<<"&";
    Rate_Latex(HighE_TH1Ds[3], 0.03, 80, mass, numevents, dailyrate*7*4.);
    std::cout<<"\\\\"<<std::endl;
    
    std::cout<<"In coinc with Al2O3";
    std::cout<<"&";
    Rate_Latex(LowE_TH1Ds[4], 0.03, 4, mass, numevents, dailyrate*7*4.);
    std::cout<<"&";
    Rate_Latex(HighE_TH1Ds[4], 0.03, 80, mass, numevents, dailyrate*7*4.);
    std::cout<<"\\\\"<<std::endl;
    
    std::cout<<"In coinc with MV";
    std::cout<<"&";
    Rate_Latex(LowE_TH1Ds[5], 0, 4, mass, numevents, dailyrate*7*4.);
    std::cout<<"&";
    Rate_Latex(HighE_TH1Ds[5], 0, 80, mass, numevents, dailyrate*7*4.);
    std::cout<<"\\\\"<<std::endl;
    
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{table}"<<std::endl;
}



void Int_day_FC(TH1D* histo, double Estart, double Eend, double MassTarget, double num_events, double daily_rate){
    
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
    
    //previous method to display incertainties
    double Integral, Error;
    
    Integral = histo->IntegralAndError(BinFirst, BinLast, Error);
    
    double Int_Rate = (Integral*daily_rate)/num_events;
    double Error_Rate = (Error*daily_rate)/num_events;
    
    //New method with TFeldmanCousin
    Double_t* Nevents = new Double_t[3];
    Double_t* Rates = new Double_t[3];
    Nevents = GetNEvents(histo, Estart, Eend);
    Rates[0] = (Nevents[0]*daily_rate)/num_events;
    Rates[1] = (Nevents[1]*daily_rate)/num_events;
    Rates[2] = (Nevents[2]*daily_rate)/num_events;
    
    cout<< "\t\t[" <<  Rates[1]/MassTarget << ";" <<  Rates[2]/MassTarget<<"]";
}

void Histograms_HighStats_LBR(TString SimuFilename = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", Double_t Phi=1.523, Double_t Plane_size=120*120, Bool_t Apply_GeResol=FALSE){
    
    //Get simulations
    TFile* Simuin = new TFile(SimuFilename, "READ");
    //TTree *np = (TTree*) Simuin->Get("Step2AlgoDir/Pri_Ntuple");
    TTree *np = (TTree*) Simuin->Get("NTupleAlgoDir/Pri_Ntuple");
    Int_t p_entries= np->GetEntries();
    //TTree *ns = (TTree*) Simuin->Get("Step2AlgoDir/Sec_Ntuple");
    TTree *ns = (TTree*) Simuin->Get("NTupleAlgoDir/Sec_Ntuple");
    Int_t s_entries= ns->GetEntries();
    
    Float_t pEDep, En, fp, sPDG, NsPidM, NsPid, EOVDiskTop, ECaWO4_0, ECaWO4_1, EAl2O3, NEvents, NEventsM, EMaxMV;
    np->SetBranchAddress(  "EMaxOV", &EOVDiskTop);
    np->SetBranchAddress(  "EMaxMV", &EMaxMV);
    np->SetBranchAddress(   "C0", &ECaWO4_0);
    np->SetBranchAddress(   "C1", &ECaWO4_1);
    np->SetBranchAddress(  "A0", &EAl2O3);
    np->SetBranchAddress(  "nEventM", &NEventsM);
    np->SetBranchAddress(  "nEvent", &NEvents);
    
    // initialize a map with all the histograms
    
    TH1D* CaWO4_0_Histogram_0_10keV = new TH1D("CaWO4_0_Histogram_0_10keV", "CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* CaWO4_1_Histogram_0_10keV = new TH1D("CaWO4_1_Histogram_0_10keV", "CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* Al2O3_Histogram_0_10keV = new TH1D("Al2O3_Histogram_0_10keV", "Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* CaWO4_0_Histogram_0_10MeV = new TH1D("CaWO4_0_Histogram_0_10MeV", "CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* CaWO4_1_Histogram_0_10MeV = new TH1D("CaWO4_1_Histogram_0_10MeV", "CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* Al2O3_Histogram_0_10MeV = new TH1D("Al2O3_Histogram_0_10MeV", "Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* CaWO4_0_Histogram_0_200eV = new TH1D("CaWO4_0_Histogram_0_200eV", "CaWO4_0_Histogram_0_200eV", 200, 0, 200);
    TH1D* CaWO4_1_Histogram_0_200eV  = new TH1D("CaWO4_1_Histogram_0_200eV", "CaWO4_1_Histogram_0_200eV", 200, 0, 200);
    TH1D* Al2O3_Histogram_0_200eV  = new TH1D("Al2O3_Histogram_0_200eV", "Al2O3_Histogram_0_200eV", 200, 0, 200);
    
    TH1D* COV_Histogram_0_100MeV = new TH1D("COV_Histogram_0_100MeV", "COV_Histogram_0_100MeV", 1000, 0, 100);
    TH1D* COV_Histogram_0_5MeV  = new TH1D("COV_Histogram_0_5MeV", "COV_Histogram_0_5MeV", 5000, 0, 5);
    
    //--- MV cuts
    TH1D* MVCut_CaWO4_0_Histogram_0_10keV = new TH1D("MVCut_CaWO4_0_Histogram_0_10keV", "MVCut_CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* MVCut_CaWO4_1_Histogram_0_10keV = new TH1D("MVCut_CaWO4_1_Histogram_0_10keV", "MVCut_CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* MVCut_Al2O3_Histogram_0_10keV = new TH1D("MVCut_Al2O3_Histogram_0_10keV", "MVCut_Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* MVCut_CaWO4_0_Histogram_0_200eV = new TH1D("MVCut_CaWO4_0_Histogram_0_200eV", "MVCut_CaWO4_0_Histogram_0_200eV", 200, 0, 200);
    TH1D* MVCut_CaWO4_1_Histogram_0_200eV  = new TH1D("MVCut_CaWO4_1_Histogram_0_200eV", "MVCut_CaWO4_1_Histogram_0_200eV", 200, 0, 200);
    TH1D* MVCut_Al2O3_Histogram_0_200eV  = new TH1D("MVCut_Al2O3_Histogram_0_200eV", "MVCut_Al2O3_Histogram_0_200eV", 200, 0, 200);
    
    TH1D* MVCut_COV_Histogram_0_100MeV = new TH1D("MVCut_COV_Histogram_0_100MeV", "MVCut_COV_Histogram_0_100MeV", 1000, 0, 100);
    TH1D* MVCut_COV_Histogram_0_5MeV  = new TH1D("MVCut_COV_Histogram_0_5MeV", "MVCut_COV_Histogram_0_5MeV", 5000, 0, 5);
    
    //--- COV cut
    TH1D* COVCut_CaWO4_0_Histogram_0_10keV = new TH1D("COVCut_CaWO4_0_Histogram_0_10keV", "COVCut_CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* COVCut_CaWO4_1_Histogram_0_10keV = new TH1D("COVCut_CaWO4_1_Histogram_0_10keV", "COVCut_CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* COVCut_Al2O3_Histogram_0_10keV = new TH1D("COVCut_Al2O3_Histogram_0_10keV", "COVCut_Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* COVCut_CaWO4_0_Histogram_0_200eV = new TH1D("COVCut_CaWO4_0_Histogram_0_200eV", "COVCut_CaWO4_0_Histogram_0_200eV", 200, 0, 200);
    TH1D* COVCut_CaWO4_1_Histogram_0_200eV  = new TH1D("COVCut_CaWO4_1_Histogram_0_200eV", "COVCut_CaWO4_1_Histogram_0_200eV", 200, 0, 200);
    TH1D* COVCut_Al2O3_Histogram_0_200eV  = new TH1D("COVCut_Al2O3_Histogram_0_200eV", "COVCut_Al2O3_Histogram_0_200eV", 200, 0, 200);
    
    //--- MV+COV cuts
    TH1D* MVCOVCut_CaWO4_0_Histogram_0_10keV = new TH1D("MVCOVCut_CaWO4_0_Histogram_0_10keV", "MVCOVCut_CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* MVCOVCut_CaWO4_1_Histogram_0_10keV = new TH1D("MVCOVCut_CaWO4_1_Histogram_0_10keV", "MVCOVCut_CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* MVCOVCut_Al2O3_Histogram_0_10keV = new TH1D("MVCOVCut_Al2O3_Histogram_0_10keV", "MVCOVCut_Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* MVCOVCut_CaWO4_0_Histogram_0_200eV = new TH1D("MVCOVCut_CaWO4_0_Histogram_0_200eV", "MVCOVCut_CaWO4_0_Histogram_0_200eV", 200, 0, 200);
    TH1D* MVCOVCut_CaWO4_1_Histogram_0_200eV  = new TH1D("MVCOVCut_CaWO4_1_Histogram_0_200eV", "MVCOVCut_CaWO4_1_Histogram_0_200eV", 200, 0, 200);
    TH1D* MVCOVCut_Al2O3_Histogram_0_200eV  = new TH1D("MVCOVCut_Al2O3_Histogram_0_200eV", "MVCOVCut_Al2O3_Histogram_0_200eV", 200, 0, 200);
    
    //------- Coincidences
    // Cryodet with COV
    TH1D* COVCoinc_CaWO4_0_Histogram_0_10keV = new TH1D("COVCoinc_CaWO4_0_Histogram_0_10keV", "COVCoinc_CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* COVCoinc_CaWO4_1_Histogram_0_10keV = new TH1D("COVCoinc_CaWO4_1_Histogram_0_10keV", "COVCoinc_CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* COVCoinc_Al2O3_Histogram_0_10keV = new TH1D("COVCoinc_Al2O3_Histogram_0_10keV", "COVCoinc_Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* COVCoinc_CaWO4_0_Histogram_0_200eV = new TH1D("COVCoinc_CaWO4_0_Histogram_0_200eV", "COVCoinc_CaWO4_0_Histogram_0_200eV", 200, 0, 200);
    TH1D* COVCoinc_CaWO4_1_Histogram_0_200eV  = new TH1D("COVCoinc_CaWO4_1_Histogram_0_200eV", "COVCoinc_CaWO4_1_Histogram_0_200eV", 200, 0, 200);
    TH1D* COVCoinc_Al2O3_Histogram_0_200eV  = new TH1D("COVCoinc_Al2O3_Histogram_0_200eV", "COVCoinc_Al2O3_Histogram_0_200eV", 200, 0, 200);
    
    // Cryodet with MV
    TH1D* MVCoinc_CaWO4_0_Histogram_0_10keV = new TH1D("MVCoinc_CaWO4_0_Histogram_0_10keV", "MVCoinc_CaWO4_0_Histogram_0_10keV", 500, 0, 10);
    TH1D* MVCoinc_CaWO4_1_Histogram_0_10keV = new TH1D("MVCoinc_CaWO4_1_Histogram_0_10keV", "MVCoinc_CaWO4_1_Histogram_0_10keV", 500, 0, 10);
    TH1D* MVCoinc_Al2O3_Histogram_0_10keV = new TH1D("MVCoinc_Al2O3_Histogram_0_10keV", "MVCoinc_Al2O3_Histogram_0_10keV", 500, 0, 10);
    
    TH1D* MVCoinc_CaWO4_0_Histogram_0_200eV = new TH1D("MVCoinc_CaWO4_0_Histogram_0_200eV", "MVCoinc_CaWO4_0_Histogram_0_200eV", 200, 0, 200);
    TH1D* MVCoinc_CaWO4_1_Histogram_0_200eV  = new TH1D("MVCoinc_CaWO4_1_Histogram_0_200eV", "MVCoinc_CaWO4_1_Histogram_0_200eV", 200, 0, 200);
    TH1D* MVCoinc_Al2O3_Histogram_0_200eV  = new TH1D("MVCoinc_Al2O3_Histogram_0_200eV", "MVCoinc_Al2O3_Histogram_0_200eV", 200, 0, 200);
    
    // COV with cryodet
    TH1D* CaWO4_0Coinc_COV_Histogram_0_100MeV = new TH1D("CaWO4_0Coinc_COV_Histogram_0_100MeV", "CaWO4_0Coinc_COV_Histogram_0_100MeV", 1000, 0, 100);
    TH1D* CaWO4_0Coinc_COV_Histogram_0_5MeV  = new TH1D("CaWO4_0Coinc_COV_Histogram_0_5MeV", "CaWO4_0Coinc_COV_Histogram_0_5MeV", 5000, 0, 5);
    
    TH1D* CaWO4_1Coinc_COV_Histogram_0_100MeV = new TH1D("CaWO4_1Coinc_COV_Histogram_0_100MeV", "CaWO4_1Coinc_COV_Histogram_0_100MeV", 1000, 0, 100);
    TH1D* CaWO4_1Coinc_COV_Histogram_0_5MeV  = new TH1D("CaWO4_1Coinc_COV_Histogram_0_5MeV", "CaWO4_1Coinc_COV_Histogram_0_5MeV", 5000, 0, 5);
    
    TH1D* Al2O3Coinc_COV_Histogram_0_100MeV = new TH1D("Al2O3Coinc_COV_Histogram_0_100MeV", "Al2O3Coinc_COV_Histogram_0_100MeV", 1000, 0, 100);
    TH1D* Al2O3Coinc_COV_Histogram_0_5MeV  = new TH1D("Al2O3Coinc_COV_Histogram_0_5MeV", "Al2O3Coinc_COV_Histogram_0_5MeV", 5000, 0, 5);
    
    TH1D* MVCoinc_COV_Histogram_0_100MeV = new TH1D("MVCoinc_COV_Histogram_0_100MeV", "MVCoinc_COV_Histogram_0_100MeV", 1000, 0, 100);
    TH1D* MVCoinc_COV_Histogram_0_5MeV  = new TH1D("MVCoinc_COV_Histogram_0_5MeV", "MVCoinc_COV_Histogram_0_5MeV", 5000, 0, 5);
    
    cout<<"Read simulation..."<<endl;
    
    unsigned long long Nev=0;
    unsigned long long NevM_prev=0;
    unsigned long long Nev_prev=0;
    for (int i=0; i<p_entries; i++){
        np->GetEntry(i);
        
        if (ECaWO4_0>0) CaWO4_0_Histogram_0_10keV->Fill(ECaWO4_0*1e3);
        if (ECaWO4_1>0) CaWO4_1_Histogram_0_10keV->Fill(ECaWO4_1*1e3);
        if (EAl2O3>0) Al2O3_Histogram_0_10keV->Fill(EAl2O3*1e3);
        
        if (ECaWO4_0>0) CaWO4_0_Histogram_0_10MeV->Fill(ECaWO4_0);
        if (ECaWO4_1>0) CaWO4_1_Histogram_0_10MeV->Fill(ECaWO4_1);
        if (EAl2O3>0) Al2O3_Histogram_0_10MeV->Fill(EAl2O3);
        
        if (ECaWO4_0>0) CaWO4_0_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
        if (ECaWO4_1>0) CaWO4_1_Histogram_0_200eV->Fill(ECaWO4_1*1e6);
        if (EAl2O3>0)  Al2O3_Histogram_0_200eV->Fill(EAl2O3*1e6);
        
        if (EOVDiskTop>0) COV_Histogram_0_5MeV->Fill(EOVDiskTop);
        if (EOVDiskTop>0) COV_Histogram_0_100MeV->Fill(EOVDiskTop);
        
        
        //Veto
        if(EOVDiskTop<0.03){ //30keV
            if (ECaWO4_0>0) COVCut_CaWO4_0_Histogram_0_10keV->Fill(ECaWO4_0*1e3);
            if (ECaWO4_1>0) COVCut_CaWO4_1_Histogram_0_10keV->Fill(ECaWO4_1*1e3);
            if (EAl2O3>0)  COVCut_Al2O3_Histogram_0_10keV->Fill(EAl2O3*1e3);
            
            if (ECaWO4_0>0) COVCut_CaWO4_0_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            if (ECaWO4_1>0) COVCut_CaWO4_1_Histogram_0_200eV->Fill(ECaWO4_1*1e6);
            if (EAl2O3>0)  COVCut_Al2O3_Histogram_0_200eV->Fill(EAl2O3*1e6);
        }
        
        if(EMaxMV<5){
            if (ECaWO4_0>0) MVCut_CaWO4_0_Histogram_0_10keV->Fill(ECaWO4_0*1e3);
            if (ECaWO4_1>0) MVCut_CaWO4_1_Histogram_0_10keV->Fill(ECaWO4_1*1e3);
            if (EAl2O3>0)  MVCut_Al2O3_Histogram_0_10keV->Fill(EAl2O3*1e3);
            
            if (ECaWO4_0>0) MVCut_CaWO4_0_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            if (ECaWO4_1>0) MVCut_CaWO4_1_Histogram_0_200eV->Fill(ECaWO4_1*1e6);
            if (EAl2O3>0)  MVCut_Al2O3_Histogram_0_200eV->Fill(EAl2O3*1e6);
            
            if (EOVDiskTop>0) MVCut_COV_Histogram_0_100MeV->Fill(EOVDiskTop);
            if (EOVDiskTop>0) MVCut_COV_Histogram_0_5MeV->Fill(EOVDiskTop);
        }
        
        if((EMaxMV<5)&&(EOVDiskTop<0.03)){
            if (ECaWO4_0>0) MVCOVCut_CaWO4_0_Histogram_0_10keV->Fill(ECaWO4_0*1e3);
            if (ECaWO4_1>0) MVCOVCut_CaWO4_1_Histogram_0_10keV->Fill(ECaWO4_1*1e3);
            if (EAl2O3>0)  MVCOVCut_Al2O3_Histogram_0_10keV->Fill(EAl2O3*1e3);
            
            if (ECaWO4_0>0) MVCOVCut_CaWO4_0_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            if (ECaWO4_1>0) MVCOVCut_CaWO4_1_Histogram_0_200eV->Fill(ECaWO4_1*1e6);
            if (EAl2O3>0) MVCOVCut_Al2O3_Histogram_0_200eV->Fill(EAl2O3*1e6);
        }
        
        // Coincidences
        // Cryodet with COV
        if(EOVDiskTop>0.03){ //30keV
            //Spectra inside cryodet
            if (ECaWO4_0>0) {
                COVCoinc_CaWO4_0_Histogram_0_10keV->Fill(ECaWO4_0*1e3);
                COVCoinc_CaWO4_0_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            }
            if (ECaWO4_1>0) {
                COVCoinc_CaWO4_1_Histogram_0_10keV->Fill(ECaWO4_1*1e3);
                COVCoinc_CaWO4_1_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            }
            if (EAl2O3>0)  {
                COVCoinc_Al2O3_Histogram_0_10keV->Fill(EAl2O3*1e3);
                COVCoinc_Al2O3_Histogram_0_200eV->Fill(EAl2O3*1e6);
            }
        }
        
        // Cryodet with MV
        if(EMaxMV>5.){ //5MeV
            
            //Spectra inside cryodet
            if (ECaWO4_0>0) {
                MVCoinc_CaWO4_0_Histogram_0_10keV->Fill(ECaWO4_0*1e3);
                MVCoinc_CaWO4_0_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            }
            if (ECaWO4_1>0) {
                MVCoinc_CaWO4_1_Histogram_0_10keV->Fill(ECaWO4_1*1e3);
                MVCoinc_CaWO4_1_Histogram_0_200eV->Fill(ECaWO4_0*1e6);
            }
            if (EAl2O3>0)  {
                MVCoinc_Al2O3_Histogram_0_10keV->Fill(EAl2O3*1e3);
                MVCoinc_Al2O3_Histogram_0_200eV->Fill(EAl2O3*1e6);
            }
        }
        
        // COV with Cryodet
        if((ECaWO4_0>50*1e-6)&&(ECaWO4_0<5*1e-3)){ //50 eV - 5keV
            if (EOVDiskTop>0) CaWO4_0Coinc_COV_Histogram_0_5MeV->Fill(EOVDiskTop);
            if (EOVDiskTop>0) CaWO4_0Coinc_COV_Histogram_0_100MeV->Fill(EOVDiskTop);
        }
        
        if((ECaWO4_1>50*1e-6)&&(ECaWO4_1<5*1e-3)){ //50 eV - 5keV
            if (EOVDiskTop>0) CaWO4_1Coinc_COV_Histogram_0_5MeV->Fill(EOVDiskTop);
            if (EOVDiskTop>0) CaWO4_1Coinc_COV_Histogram_0_100MeV->Fill(EOVDiskTop);
        }
        
        if((EAl2O3>50*1e-6)&&(EAl2O3<5*1e-3)){ //50 eV - 5keV
            if (EOVDiskTop>0) Al2O3Coinc_COV_Histogram_0_5MeV->Fill(EOVDiskTop);
            if (EOVDiskTop>0) Al2O3Coinc_COV_Histogram_0_100MeV->Fill(EOVDiskTop);
        }
        
        if(EMaxMV>5.){ //5 MeV
            if (EOVDiskTop>0) MVCoinc_COV_Histogram_0_5MeV->Fill(EOVDiskTop);
            if (EOVDiskTop>0) MVCoinc_COV_Histogram_0_100MeV->Fill(EOVDiskTop);
        }
        
        if ((NEventsM<NevM_prev)||((NEvents<Nev_prev)&&(NEventsM==0))){ //there is a reset!
            Nev+=NevM_prev*1e6+Nev_prev;
        }
        
        NevM_prev=NEventsM;
        Nev_prev = NEvents;
    }
    np->GetEntry(p_entries-1);
    Nev+=NEventsM*1e6+NEvents;
    
    
    //Print Latex table in a txt
    TString SimuFilename_txt = SimuFilename;
    TPMERegexp("Analyzer_").Substitute(SimuFilename_txt,"LatexTable_");
    TPMERegexp(".root").Substitute(SimuFilename_txt,".txt");
    
    std::ofstream ofs{SimuFilename_txt};
    auto cout_buff = std::cout.rdbuf();
    std::cout.rdbuf(ofs.rdbuf());
    
    // Prepare vectors of hitograms for printing
    std::vector<TH1D*> COV_LowE = {COV_Histogram_0_5MeV, MVCut_COV_Histogram_0_5MeV, CaWO4_0Coinc_COV_Histogram_0_5MeV, CaWO4_1Coinc_COV_Histogram_0_5MeV, Al2O3Coinc_COV_Histogram_0_5MeV, MVCoinc_COV_Histogram_0_5MeV};
    std::vector<TH1D*> COV_HighE = {COV_Histogram_0_100MeV, MVCut_COV_Histogram_0_100MeV, CaWO4_0Coinc_COV_Histogram_0_100MeV, CaWO4_1Coinc_COV_Histogram_0_100MeV, Al2O3Coinc_COV_Histogram_0_100MeV, MVCoinc_COV_Histogram_0_100MeV};
    
    std::vector<std::vector<TH1D*>> Cryodet_LowE = {{CaWO4_0_Histogram_0_200eV, CaWO4_1_Histogram_0_200eV, Al2O3_Histogram_0_200eV}, {MVCut_CaWO4_0_Histogram_0_200eV, MVCut_CaWO4_1_Histogram_0_200eV, MVCut_Al2O3_Histogram_0_200eV}, {COVCut_CaWO4_0_Histogram_0_200eV, COVCut_CaWO4_1_Histogram_0_200eV, COVCut_Al2O3_Histogram_0_200eV}, {MVCOVCut_CaWO4_0_Histogram_0_200eV, MVCOVCut_CaWO4_1_Histogram_0_200eV, MVCOVCut_Al2O3_Histogram_0_200eV}, {COVCoinc_CaWO4_0_Histogram_0_200eV, COVCoinc_CaWO4_1_Histogram_0_200eV, COVCoinc_Al2O3_Histogram_0_200eV}, {MVCoinc_CaWO4_0_Histogram_0_200eV, MVCoinc_CaWO4_1_Histogram_0_200eV, MVCoinc_Al2O3_Histogram_0_200eV}};
    
    std::vector<std::vector<TH1D*>> Cryodet_HighE = {{CaWO4_0_Histogram_0_10keV, CaWO4_1_Histogram_0_10keV, Al2O3_Histogram_0_10keV}, {MVCut_CaWO4_0_Histogram_0_10keV, MVCut_CaWO4_1_Histogram_0_10keV, MVCut_Al2O3_Histogram_0_10keV}, {COVCut_CaWO4_0_Histogram_0_10keV, COVCut_CaWO4_1_Histogram_0_10keV, COVCut_Al2O3_Histogram_0_10keV}, {MVCOVCut_CaWO4_0_Histogram_0_10keV, MVCOVCut_CaWO4_1_Histogram_0_10keV, MVCOVCut_Al2O3_Histogram_0_10keV}, {COVCoinc_CaWO4_0_Histogram_0_10keV, COVCoinc_CaWO4_1_Histogram_0_10keV, COVCoinc_Al2O3_Histogram_0_10keV}, {MVCoinc_CaWO4_0_Histogram_0_10keV, MVCoinc_CaWO4_1_Histogram_0_10keV, MVCoinc_Al2O3_Histogram_0_10keV}};
    
    std::cout<<"Latex table for cryodetectors"<<std::endl;
    Print_Table_Latex_cryodet({"CaWO4 0", "CaWO4 1", "Al2O3"}, Cryodet_LowE, Cryodet_HighE, {MASS_CaWO4_0, MASS_CaWO4_1, MASS_Al2O3}, Nev, (Phi*Plane_size*24*3600.));
    
    std::cout<<std::endl;
    std::cout<<"==========================================================================="<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<"Latex table for COV"<<std::endl;
    Print_Table_Latex_COV(COV_LowE, COV_HighE, MASS_Ge,  Nev, (Phi*Plane_size*24*3600.));
    
    std::cout<<std::endl;
    std::cout<<"==========================================================================="<<std::endl;
    std::cout<<std::endl;
        
        
    // Save in a file
    TString SimuFilename_out = SimuFilename;
    TPMERegexp("Analyzer_").Substitute(SimuFilename_out,"Histos_");
    
    TFile* fileout = new TFile(SimuFilename_out, "RECREATE");
    
    Double_t EqTime = Nev/(Phi*Plane_size); //s
    double EqTimeDay = Nev/(Phi*Plane_size*24*3600.);
    std::cout<<"Equivalent time in days = "<<EqTimeDay<<std::endl;
    std::cout<<"Equivalent time in years = "<<EqTimeDay/365<<std::endl;
    std::cout<<Nev<<std::endl;
    
    //restore buffer and print
    std::cout.rdbuf(cout_buff);
    std::ifstream ifs{SimuFilename_txt};
    std::cout << ifs.rdbuf();
    
    //Write in the file the flux, num of primaries and norm
    TParameter<int>* Prim = new TParameter<int>("NPrim", Nev);
    Prim->Write();
    TParameter<float>* EqTime_s = new TParameter<float>("EqTime_s", EqTime);
    EqTime_s->Write();
    
    //Legend
    CaWO4_0_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    CaWO4_1_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    Al2O3_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    
    CaWO4_0_Histogram_0_10MeV->GetXaxis()->SetTitle("Energy [MeV]");
    CaWO4_1_Histogram_0_10MeV->GetXaxis()->SetTitle("Energy [MeV]");
    Al2O3_Histogram_0_10MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    CaWO4_0_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    CaWO4_1_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    Al2O3_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    
    COV_Histogram_0_5MeV->GetXaxis()->SetTitle("Energy [MeV]");
    COV_Histogram_0_100MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    //--- MV cuts
    MVCut_CaWO4_0_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    MVCut_CaWO4_1_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    MVCut_Al2O3_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    
    MVCut_CaWO4_0_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    MVCut_CaWO4_1_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    MVCut_Al2O3_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    
    MVCut_COV_Histogram_0_100MeV->GetXaxis()->SetTitle("Energy [MeV]");
    MVCut_COV_Histogram_0_5MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    //--- COV cut
    COVCut_CaWO4_0_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    COVCut_CaWO4_1_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    COVCut_Al2O3_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    
    COVCut_CaWO4_0_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    COVCut_CaWO4_1_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    COVCut_Al2O3_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    
    //--- MV+COV cuts
    MVCOVCut_CaWO4_0_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    MVCOVCut_CaWO4_1_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    MVCOVCut_Al2O3_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    
    MVCOVCut_CaWO4_0_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    MVCOVCut_CaWO4_1_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    MVCOVCut_Al2O3_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    
    //------- Coincidences
    // Cryodet with COV
    COVCoinc_CaWO4_0_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    COVCoinc_CaWO4_1_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    COVCoinc_Al2O3_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    
    COVCoinc_CaWO4_0_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    COVCoinc_CaWO4_1_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    COVCoinc_Al2O3_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    
    // Cryodet with MV
    MVCoinc_CaWO4_0_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    MVCoinc_CaWO4_1_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    MVCoinc_Al2O3_Histogram_0_10keV->GetXaxis()->SetTitle("Energy [keV]");
    
    MVCoinc_CaWO4_0_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    MVCoinc_CaWO4_1_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    MVCoinc_Al2O3_Histogram_0_200eV->GetXaxis()->SetTitle("Energy [eV]");
    
    // COV with cryodet
    CaWO4_0Coinc_COV_Histogram_0_100MeV->GetXaxis()->SetTitle("Energy [MeV]");
    CaWO4_0Coinc_COV_Histogram_0_5MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    CaWO4_1Coinc_COV_Histogram_0_100MeV->GetXaxis()->SetTitle("Energy [MeV]");
    CaWO4_1Coinc_COV_Histogram_0_5MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    Al2O3Coinc_COV_Histogram_0_100MeV->GetXaxis()->SetTitle("Energy [MeV]");
    Al2O3Coinc_COV_Histogram_0_5MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    MVCoinc_COV_Histogram_0_100MeV->GetXaxis()->SetTitle("Energy [MeV]");
    MVCoinc_COV_Histogram_0_5MeV->GetXaxis()->SetTitle("Energy [MeV]");
    
    
    //------------- Scale in counts/day
    CaWO4_0_Histogram_0_10keV->Scale(1./EqTimeDay);
    CaWO4_1_Histogram_0_10keV->Scale(1./EqTimeDay);
    Al2O3_Histogram_0_10keV->Scale(1./EqTimeDay);
    
    CaWO4_0_Histogram_0_10MeV->Scale(1./EqTimeDay);
    CaWO4_1_Histogram_0_10MeV->Scale(1./EqTimeDay);
    Al2O3_Histogram_0_10MeV->Scale(1./EqTimeDay);
    
    CaWO4_0_Histogram_0_200eV->Scale(1./EqTimeDay);
    CaWO4_1_Histogram_0_200eV->Scale(1./EqTimeDay);
    Al2O3_Histogram_0_200eV->Scale(1./EqTimeDay);
    
    COV_Histogram_0_100MeV->Scale(1./EqTimeDay);
    COV_Histogram_0_5MeV->Scale(1./EqTimeDay);
    
    COVCut_CaWO4_0_Histogram_0_10keV->Scale(1./EqTimeDay);
    COVCut_CaWO4_1_Histogram_0_10keV->Scale(1./EqTimeDay);
    COVCut_Al2O3_Histogram_0_10keV->Scale(1./EqTimeDay);
    
    COVCut_CaWO4_0_Histogram_0_200eV->Scale(1./EqTimeDay);
    COVCut_CaWO4_1_Histogram_0_200eV->Scale(1./EqTimeDay);
    COVCut_Al2O3_Histogram_0_200eV->Scale(1./EqTimeDay);
    
    MVCut_CaWO4_0_Histogram_0_10keV->Scale(1./EqTimeDay);
    MVCut_CaWO4_1_Histogram_0_10keV->Scale(1./EqTimeDay);
    MVCut_Al2O3_Histogram_0_10keV->Scale(1./EqTimeDay);
    
    MVCut_CaWO4_0_Histogram_0_200eV->Scale(1./EqTimeDay);
    MVCut_CaWO4_1_Histogram_0_200eV->Scale(1./EqTimeDay);
    MVCut_Al2O3_Histogram_0_200eV->Scale(1./EqTimeDay);
    
    MVCut_COV_Histogram_0_100MeV->Scale(1./EqTimeDay);
    MVCut_COV_Histogram_0_5MeV->Scale(1./EqTimeDay);
    
    MVCOVCut_CaWO4_0_Histogram_0_10keV->Scale(1./EqTimeDay);
    MVCOVCut_CaWO4_1_Histogram_0_10keV->Scale(1./EqTimeDay);
    MVCOVCut_Al2O3_Histogram_0_10keV->Scale(1./EqTimeDay);
    
    MVCOVCut_CaWO4_0_Histogram_0_200eV->Scale(1./EqTimeDay);
    MVCOVCut_CaWO4_1_Histogram_0_200eV->Scale(1./EqTimeDay);
    MVCOVCut_Al2O3_Histogram_0_200eV->Scale(1./EqTimeDay);
    
    
    //------- Coincidences
    // Cryodet with COV
    COVCoinc_CaWO4_0_Histogram_0_10keV->Scale(1./EqTimeDay);
    COVCoinc_CaWO4_1_Histogram_0_10keV->Scale(1./EqTimeDay);
    COVCoinc_Al2O3_Histogram_0_10keV->Scale(1./EqTimeDay);
    
    COVCoinc_CaWO4_0_Histogram_0_200eV->Scale(1./EqTimeDay);
    COVCoinc_CaWO4_1_Histogram_0_200eV->Scale(1./EqTimeDay);
    COVCoinc_Al2O3_Histogram_0_200eV->Scale(1./EqTimeDay);
    
    // Cryodet with MV
    MVCoinc_CaWO4_0_Histogram_0_10keV->Scale(1./EqTimeDay);
    MVCoinc_CaWO4_1_Histogram_0_10keV->Scale(1./EqTimeDay);
    MVCoinc_Al2O3_Histogram_0_10keV->Scale(1./EqTimeDay);
    
    MVCoinc_CaWO4_0_Histogram_0_200eV->Scale(1./EqTimeDay);
    MVCoinc_CaWO4_1_Histogram_0_200eV->Scale(1./EqTimeDay);
    MVCoinc_Al2O3_Histogram_0_200eV->Scale(1./EqTimeDay);
    
    // COV with cryodet
    CaWO4_0Coinc_COV_Histogram_0_100MeV->Scale(1./EqTimeDay);
    CaWO4_0Coinc_COV_Histogram_0_5MeV->Scale(1./EqTimeDay);
    
    CaWO4_1Coinc_COV_Histogram_0_100MeV->Scale(1./EqTimeDay);
    CaWO4_1Coinc_COV_Histogram_0_5MeV->Scale(1./EqTimeDay);
    
    Al2O3Coinc_COV_Histogram_0_100MeV->Scale(1./EqTimeDay);
    Al2O3Coinc_COV_Histogram_0_5MeV->Scale(1./EqTimeDay);
    
    MVCoinc_COV_Histogram_0_100MeV->Scale(1./EqTimeDay);
    MVCoinc_COV_Histogram_0_5MeV->Scale(1./EqTimeDay);
    
    
    //------------ Legend
    CaWO4_0_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    CaWO4_1_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    Al2O3_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    CaWO4_0_Histogram_0_10MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    CaWO4_1_Histogram_0_10MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    Al2O3_Histogram_0_10MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    CaWO4_0_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    CaWO4_1_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    Al2O3_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    COV_Histogram_0_100MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    COV_Histogram_0_5MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    COVCut_CaWO4_0_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCut_CaWO4_1_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCut_Al2O3_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    COVCut_CaWO4_0_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCut_CaWO4_1_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCut_Al2O3_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCut_CaWO4_0_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCut_CaWO4_1_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCut_Al2O3_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCut_CaWO4_0_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCut_CaWO4_1_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCut_Al2O3_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCut_COV_Histogram_0_100MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCut_COV_Histogram_0_5MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCOVCut_CaWO4_0_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCOVCut_CaWO4_1_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCOVCut_Al2O3_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCOVCut_CaWO4_0_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCOVCut_CaWO4_1_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCOVCut_Al2O3_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    //------- Coincidences
    // Cryodet with COV
    COVCoinc_CaWO4_0_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCoinc_CaWO4_1_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCoinc_Al2O3_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    COVCoinc_CaWO4_0_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCoinc_CaWO4_1_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    COVCoinc_Al2O3_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    // Cryodet with MV
    MVCoinc_CaWO4_0_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCoinc_CaWO4_1_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCoinc_Al2O3_Histogram_0_10keV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCoinc_CaWO4_0_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCoinc_CaWO4_1_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCoinc_Al2O3_Histogram_0_200eV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    // COV with cryodet
    CaWO4_0Coinc_COV_Histogram_0_100MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    CaWO4_0Coinc_COV_Histogram_0_5MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    CaWO4_1Coinc_COV_Histogram_0_100MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    CaWO4_1Coinc_COV_Histogram_0_5MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    Al2O3Coinc_COV_Histogram_0_100MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    Al2O3Coinc_COV_Histogram_0_5MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    MVCoinc_COV_Histogram_0_100MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    MVCoinc_COV_Histogram_0_5MeV->GetYaxis()->SetTitle("Rate [ev/day]");
    
    CaWO4_0_Histogram_0_10keV->Write();
    CaWO4_1_Histogram_0_10keV->Write();
    Al2O3_Histogram_0_10keV->Write();
    
    CaWO4_0_Histogram_0_10MeV->Write();
    CaWO4_1_Histogram_0_10MeV->Write();
    Al2O3_Histogram_0_10MeV->Write();
    
    CaWO4_0_Histogram_0_200eV->Write();
    CaWO4_1_Histogram_0_200eV->Write();
    Al2O3_Histogram_0_200eV->Write();
    
    COV_Histogram_0_5MeV->Write();
    COV_Histogram_0_100MeV->Write();
    
    fileout->mkdir("COVCut");
    fileout->cd("COVCut");
    
    COVCut_CaWO4_0_Histogram_0_10keV->Write();
    COVCut_CaWO4_1_Histogram_0_10keV->Write();
    COVCut_Al2O3_Histogram_0_10keV->Write();
    
    COVCut_CaWO4_0_Histogram_0_200eV->Write();
    COVCut_CaWO4_1_Histogram_0_200eV->Write();
    COVCut_Al2O3_Histogram_0_200eV->Write();
    
    fileout->cd();
    
    fileout->mkdir("MVCut");
    fileout->cd("MVCut");
    
    MVCut_CaWO4_0_Histogram_0_10keV->Write();
    MVCut_CaWO4_1_Histogram_0_10keV->Write();
    MVCut_Al2O3_Histogram_0_10keV->Write();
    
    MVCut_CaWO4_0_Histogram_0_200eV->Write();
    MVCut_CaWO4_1_Histogram_0_200eV->Write();
    MVCut_Al2O3_Histogram_0_200eV->Write();
    
    MVCut_COV_Histogram_0_5MeV->Write();
    MVCut_COV_Histogram_0_100MeV->Write();
    
    fileout->cd();
    
    fileout->mkdir("MVCOVCut");
    fileout->cd("MVCOVCut");
    
    MVCOVCut_CaWO4_0_Histogram_0_10keV->Write();
    MVCOVCut_CaWO4_1_Histogram_0_10keV->Write();
    MVCOVCut_Al2O3_Histogram_0_10keV->Write();
    
    MVCOVCut_CaWO4_0_Histogram_0_200eV->Write();
    MVCOVCut_CaWO4_1_Histogram_0_200eV->Write();
    MVCOVCut_Al2O3_Histogram_0_200eV->Write();
    
    //------- Coincidences
    // Cryodet with COV
    fileout->mkdir("COVcoinc");
    fileout->cd("COVcoinc");
    
    COVCoinc_CaWO4_0_Histogram_0_10keV->Write();
    COVCoinc_CaWO4_1_Histogram_0_10keV->Write();
    COVCoinc_Al2O3_Histogram_0_10keV->Write();
    
    COVCoinc_CaWO4_0_Histogram_0_200eV->Write();
    COVCoinc_CaWO4_1_Histogram_0_200eV->Write();
    COVCoinc_Al2O3_Histogram_0_200eV->Write();
    
    // Cryodet with MV
    fileout->mkdir("MVcoinc");
    fileout->cd("MVcoinc");
    MVCoinc_CaWO4_0_Histogram_0_10keV->Write();
    MVCoinc_CaWO4_1_Histogram_0_10keV->Write();
    MVCoinc_Al2O3_Histogram_0_10keV->Write();
    
    MVCoinc_CaWO4_0_Histogram_0_200eV->Write();
    MVCoinc_CaWO4_1_Histogram_0_200eV->Write();
    MVCoinc_Al2O3_Histogram_0_200eV->Write();
    
    MVCoinc_COV_Histogram_0_100MeV->Write();
    MVCoinc_COV_Histogram_0_5MeV->Write();
    
    // COV with cryodet
    fileout->mkdir("Cryodet_coinc");
    fileout->cd("Cryodet_coinc");
    CaWO4_0Coinc_COV_Histogram_0_100MeV->Write();
    CaWO4_0Coinc_COV_Histogram_0_5MeV->Write();
    
    CaWO4_1Coinc_COV_Histogram_0_100MeV->Write();
    CaWO4_1Coinc_COV_Histogram_0_5MeV->Write();
    
    Al2O3Coinc_COV_Histogram_0_100MeV->Write();
    Al2O3Coinc_COV_Histogram_0_5MeV->Write();
    //---
    
    fileout->Close();
}

