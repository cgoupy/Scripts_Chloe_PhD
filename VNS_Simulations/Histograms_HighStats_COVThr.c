Double_t Phi=1.523 ;//ev/s/cm-2/sec 0.271 K40, 1.523 U238, 0.618 Th232
Double_t MASS_CaWO4_array = 9*0.7575*1e-3;
Double_t MASS_Al2O3_array = 9*0.4975*1e-3;

std::map<TString, std::unique_ptr<TH1D>> mh;

Double_t* GetNEvents(TString hName, Double_t xmin = 0., Double_t xmax = 0.)
{
    double BinFirst = mh[hName]->FindBin(xmin);
    double BinLast=  mh[hName]->FindBin(xmax);

    double NEvents =  mh[hName]->Integral(BinFirst, BinLast);
    
    Double_t* int_ = new Double_t[3];
    int_[0] = NEvents;
    //Use TFeldmanCousin CL method if low stat.
    if (NEvents < 30) {
    TFeldmanCousins fc{0.95};
    int_[1] = fc.CalculateLowerLimit(NEvents,0);
    int_[2] = fc.CalculateUpperLimit(NEvents,0);
    }
    else {
    int_[1] = NEvents-std::sqrt(NEvents);
    int_[2] = NEvents+std::sqrt(NEvents);
    }
    return &int_[0];
}

void Int_day(TString hName, double Estart, double Eend, double MassTarget, double eq_time_days){
    
    double dE = mh[hName]->GetBinWidth(1);
    double Nbins = mh[hName]->GetNbinsX();
    double BinFirst = mh[hName]->FindBin(Estart);
    double BinLast= mh[hName]->FindBin(Eend);
    TString Unit= mh[hName]->GetXaxis()->GetTitle();
    
    if (Estart==0){
        BinFirst= 1;
    }
    if (Eend==0) {
        BinLast = Nbins;
    }
    
    //previous method to display incertainties
    double Integral, Error;
    
    Integral = mh[hName]->IntegralAndError(BinFirst, BinLast, Error);
    
    double Int_Rate = Integral*1./eq_time_days;
    double Error_Rate = Error*1./eq_time_days;
    
    //New method with TFeldmanCousin
    Double_t* Nevents = new Double_t[3];
    Double_t* Rates = new Double_t[3];
    Nevents = GetNEvents(hName, Estart, Eend);
    Rates[0] = Nevents[0]*1./eq_time_days;
    Rates[1] = Nevents[1]*1./eq_time_days;
    Rates[2] = Nevents[2]*1./eq_time_days;
    
    if (Int_Rate>0) cout<< Int_Rate/MassTarget << ";" << Error_Rate/MassTarget<<";";
    else cout<<"<"<< Rates[2]/MassTarget<<";;";
    
    
}

void Int_day_FC(TString hName, double Estart, double Eend, double MassTarget, double num_events, double daily_rate){
    
    double dE = mh[hName]->GetBinWidth(1);
    double Nbins = mh[hName]->GetNbinsX();
    
    double BinFirst = mh[hName]->FindBin(Estart);
    double BinLast= mh[hName]->FindBin(Eend);
    TString Unit= mh[hName]->GetXaxis()->GetTitle();
    
    if (Estart==0){
        BinFirst= 1;
    }
    if (Eend==0) {
        BinLast = Nbins;
    }
    
    //previous method to display incertainties
    double Integral, Error;
    
    Integral = mh[hName]->IntegralAndError(BinFirst, BinLast, Error);
    
    double Int_Rate = (Integral*daily_rate)/num_events;
    double Error_Rate = (Error*daily_rate)/num_events;
    
    //New method with TFeldmanCousin
    Double_t* Nevents = new Double_t[3];
    Double_t* Rates = new Double_t[3];
    Nevents = GetNEvents(hName, Estart, Eend);
    Rates[0] = (Nevents[0]*daily_rate)/num_events;
    Rates[1] = (Nevents[1]*daily_rate)/num_events;
    Rates[2] = (Nevents[2]*daily_rate)/num_events;
    
    cout<< "\t\t[" <<  Rates[1]/MassTarget << ";" <<  Rates[2]/MassTarget<<"]";
}

void InitHistograms(Int_t NEdepRange, std::pair<TString,Double_t> *EdepRange, int Ncuts, TString *CutNames, TString Crystal="CaWO4"){
    
    
    // ========= ACCEPTED ===========
    for (UInt_t k=0; k<Ncuts; k++){
        //Energy deposits for accepted events in lin scale and in different Edep ranges
        for (UInt_t p=0; p<NEdepRange; p++){
            TString EdepName = Form("0-%g%s",EdepRange[p].second,EdepRange[p].first.Data());
            TString hName = Form("Edep_%s_array_"+Crystal+"_%s", CutNames[k].Data(), EdepName.Data());
            std::cout<<hName<<std::endl;
            mh[hName] = make_unique<TH1D>("n"+hName,"Edep / crystal array",1000,0,EdepRange[p].second); //units here in MeV !!
            mh[hName]->GetXaxis()->SetTitle(Form("E_{dep} [%s]",EdepRange[p].first.Data()));
        }
    }
    
    //Histograms axis settings
    for (auto it=mh.begin(); it!=mh.end(); ++it) { it->second->GetYaxis()->SetTitle("Rate [counts/day]"); }
    //for (auto it=mh.begin(); it!=mh.end(); ++it) { it->second->GetYaxis()->SetTitle("Hits []"); }
}

void FillHistograms(Double_t Edep, TString CutName, TString Crystal="CaWO4"){
    Double_t convert_coef=1.;
    
    //filling lin scale histograms in different energy deposit ranges
    if (Edep > 0 && Edep/convert_coef/1e3 < 100)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-100GeV",CutName.Data())]->Fill(Edep/convert_coef/1e3); //unit in GeV here
    }
    
    if (Edep > 0 && Edep/convert_coef/1e3 < 10)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-10GeV",CutName.Data())]->Fill(Edep/convert_coef/1e3); //unit in GeV here
    }
    
    if (Edep > 0 && Edep/convert_coef/1e3 < 1)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-1GeV",CutName.Data())]->Fill(Edep/convert_coef/1e3); //unit in GeV here
    }
    
    if (Edep > 0 && Edep/convert_coef < 100)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-100MeV",CutName.Data())]->Fill(Edep/convert_coef); //unit in MeV here
    }
    
    if (Edep > 0 && Edep/convert_coef < 10)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-10MeV",CutName.Data())]->Fill(Edep/convert_coef); //unit in MeV here
    }
    
    if (Edep >0 && Edep/convert_coef < 1)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-1MeV",CutName.Data())]->Fill(Edep/convert_coef); //unit in MeV here
    }
    
    if (Edep >0 && Edep/convert_coef*1e3 < 100)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-100keV",CutName.Data())]->Fill(Edep/convert_coef*1e3); //unit in keV here
    }
    
    if (Edep >0 && Edep/convert_coef*1e3 < 10)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-10keV",CutName.Data())]->Fill(Edep/convert_coef*1e3); //unit in keV here
    }
    
    if (Edep >0 && Edep/convert_coef*1e3 < 1)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-1keV",CutName.Data())]->Fill(Edep/convert_coef*1e3); //unit in keV here
    }
    
    if (Edep >0 && Edep/convert_coef*1e6 < 100)
    {
        mh[Form("Edep_%s_array_"+Crystal+"_0-100eV",CutName.Data())]->Fill(Edep/convert_coef*1e6); //unit in eV here
    }
}

Double_t ConvertCoeff(TString unit_name){
    if (unit_name == "meV"){
        return 1.e9;
    }
    
    if (unit_name == "eV"){
        return 1.e6;
    }
    
    if (unit_name == "keV"){
        return 1.e3;
    }
    
    if (unit_name == "MeV"){
        return 1.;
    }
    
    if (unit_name == "GeV"){
        return 1.e-3;
    }
    return 1.;
}

void SaveHistos(TFile* fileout, Double_t eq_time_d, unsigned long long Nev, Int_t NEdepRange, std::pair<TString,Double_t> *EdepRange, int Ncuts, TString *CutNames, double EqTimeDay, TString Crystal="CaWO4"){
    
    fileout->mkdir(Crystal);
    fileout->cd(Crystal);
    std::cout<<Crystal<<std::endl;
    
    Double_t mass;
    if (Crystal=="CaWO4"){
        mass=MASS_CaWO4_array;
    }
    else{
        mass=MASS_Al2O3_array;
    }
    
    std::cout<<"Cut;\t\t 10-100 eV;   Error; 0.1-1 keV; Error;   1-10 keV; Error;   10-100 keV; Error;   0.1-1 MeV; Error";
    for (UInt_t k=0; k<Ncuts; k++){
        fileout->mkdir(Crystal+"/"+CutNames[k]);
        std::cout<<"\n"<<CutNames[k].Data()<<";";
        //Energy deposits for accepted events in lin scale and in different Edep ranges
        for (UInt_t p=0; p<NEdepRange; p++){
            TString EdepName = Form("0-%g%s",EdepRange[p].second,EdepRange[p].first.Data());
            TString hName = Form("Edep_%s_array_"+Crystal+"_%s", CutNames[k].Data(), EdepName.Data());
            
            if (p==0) Int_day(hName, (10*ConvertCoeff(EdepRange[p].first.Data()))/1e6, EdepRange[p].second, mass, EqTimeDay);
            else Int_day(hName, EdepRange[p-1].second*ConvertCoeff(EdepRange[p].first.Data())/ConvertCoeff(EdepRange[p-1].first.Data()), EdepRange[p].second, mass, EqTimeDay);
            
            //mh[hName]->Write();

            fileout->cd(Crystal+"/"+CutNames[k]);
            mh[hName]->Scale(1./eq_time_d);
            mh[hName]->Write();
            
        }
    }
    std::cout<<"\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
}


void PlotEventsCryoDet(TString SimuFilename = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", Double_t Phi=1.523, Double_t Pane_size=120*120, Bool_t NeutronQuenched=FALSE, TString Algo="Z"){
    
    //Get simulations
    TFile* Simuin = new TFile(SimuFilename, "READ");
    TTree *np, *ns;
    if (Algo=="T") np = (TTree*) Simuin->Get("NTupleAlgoDir/Pri_Ntuple");
    else np = (TTree*) Simuin->Get("Step2AlgoDir/Pri_Ntuple");
    Int_t p_entries= np->GetEntries();
    if (Algo=="T") ns = (TTree*) Simuin->Get("NTupleAlgoDir/Sec_Ntuple");
    else ns = (TTree*) Simuin->Get("Step2AlgoDir/Sec_Ntuple");
    Int_t s_entries= ns->GetEntries();
    
    Float_t pEDep, En, fp, sPDG, NsPidM, NsPid, EMaxOV, EAl2O3, NEvents, NEventsM, EMaxMV, EMaxIV, TgtHits, Delay, EMaxOV_q;
    Float_t C[9], A[9];
    np->SetBranchAddress(  "EMaxOV", &EMaxOV);
    np->SetBranchAddress(  "EMaxMV", &EMaxMV);
    
    np->SetBranchAddress(   "C0", &C[0]);
    np->SetBranchAddress(   "C1", &C[1]);
    np->SetBranchAddress(   "C2", &C[2]);
    np->SetBranchAddress(   "C3", &C[3]);
    np->SetBranchAddress(   "C4", &C[4]);
    np->SetBranchAddress(   "C5", &C[5]);
    np->SetBranchAddress(   "C6", &C[6]);
    np->SetBranchAddress(   "C7", &C[7]);
    np->SetBranchAddress(   "C8", &C[8]);
    
    np->SetBranchAddress(   "A0", &A[0]);
    np->SetBranchAddress(   "A1", &A[1]);
    np->SetBranchAddress(   "A2", &A[2]);
    np->SetBranchAddress(   "A3", &A[3]);
    np->SetBranchAddress(   "A4", &A[4]);
    np->SetBranchAddress(   "A5", &A[5]);
    np->SetBranchAddress(   "A6", &A[6]);
    np->SetBranchAddress(   "A7", &A[7]);
    np->SetBranchAddress(   "A8", &A[8]);
    
    np->SetBranchAddress(  "EAl2O3", &EAl2O3);
    np->SetBranchAddress(  "EMaxIV", &EMaxIV);
    np->SetBranchAddress(  "TgtHits", &TgtHits);
    np->SetBranchAddress(  "nEventM", &NEventsM);
    np->SetBranchAddress(  "nEvent", &NEvents);
    np->SetBranchAddress(  "Delay", &Delay);

    //Initialize histograms
    static const Int_t NEdepRange = 10;
    std::pair<TString,Double_t> EdepRange[NEdepRange] = {
        std::make_pair("eV",100),
        std::make_pair("keV",1),
        std::make_pair("keV",10),
        std::make_pair("keV",100),
        std::make_pair("MeV",1),
        std::make_pair("MeV",10),
        std::make_pair("MeV",100),
        std::make_pair("GeV",1),
        std::make_pair("GeV",10),
        std::make_pair("GeV",100)
    };
    TString CutNames[32]={
        "NoCut",
        "NoCut_Delayed",
        "MVCut",
        "MVCut_Delayed",
        "COVCut_0_5keVee",
        "COVCut_1keVee",
        "COVCut_10keVee",
        "COVCut_5keVee",
        "COVCut_50keVee",
        "COVCut_80keVee",
        "COVCut_100keVee",
        "COVCut_200keVee",
        "COVCut_300keVee",
        "COVCut_400keVee",
        "COVCut_500keVee",
        "COVCut_1keVee_Delayed",
        "COVCut_5keVee_Delayed",
        "COVCut_10keVee_Delayed",
        "COVCut_50keVee_Delayed",
        "MV&COVCuts_1keVee",
        "MV&COVCuts_5keVee",
        "MV&COVCuts_Delayed",
        "AllCuts_0_5keVee",
        "AllCuts_1keVee",
        "AllCuts_2keVee",
        "AllCuts_3keVee",
        "AllCuts_5keVee",
        "AllCuts_10keVee",
        "AllCuts_20keVee",
        "AllCuts_30keVee",
        "AllCuts_50keVee",
        "AllCuts_Delayed"
    };

    
    InitHistograms(NEdepRange, EdepRange, 32, CutNames, "CaWO4");
    InitHistograms(NEdepRange, EdepRange, 32, CutNames, "Al2O3");
    
    cout<<"Read simulation..."<<endl;
    
    unsigned long long Nev=0;
    unsigned long long NevM_prev=0;
    unsigned long long Nev_prev=0;

    for (int i=0; i<p_entries; i++){
        np->GetEntry(i);
        
        if (NeutronQuenched) EMaxOV_q = EMaxOV/5;
        else EMaxOV_q = EMaxOV;
    
        // loop over CaWO4 crystals:
        for (int k=0; k<9; k++){
            
            FillHistograms(C[k], "NoCut");
            FillHistograms(A[k], "NoCut", "Al2O3");
            
            if ((Delay>1)){
                FillHistograms(C[k], "NoCut_Delayed");
                FillHistograms(A[k], "NoCut_Delayed", "Al2O3");
            }
            
            if (EMaxMV<5){
                FillHistograms(C[k], "MVCut");
                FillHistograms(A[k], "MVCut", "Al2O3");
            }
            
            if ((EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "MVCut_Delayed");
                FillHistograms(A[k], "MVCut_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<0.5e-3){
                FillHistograms(C[k], "COVCut_0_5keVee");
                FillHistograms(A[k], "COVCut_0_5keVee", "Al2O3");
            }
            
            if (EMaxOV_q<1e-3){
                FillHistograms(C[k], "COVCut_1keVee");
                FillHistograms(A[k], "COVCut_1keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_1keVee_Delayed");
                FillHistograms(A[k], "COVCut_1keVee_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<10e-3){
                FillHistograms(C[k], "COVCut_10keVee");
                FillHistograms(A[k], "COVCut_10keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<10e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_10keVee_Delayed");
                FillHistograms(A[k], "COVCut_10keVee_Delayed", "Al2O3");
            }
                
            if (EMaxOV_q<5e-3){
                FillHistograms(C[k], "COVCut_5keVee");
                FillHistograms(A[k], "COVCut_5keVee", "Al2O3");
            }
            
            if (EMaxOV_q<80e-3){
                FillHistograms(C[k], "COVCut_80keVee");
                FillHistograms(A[k], "COVCut_80keVee", "Al2O3");
            }
            
            if (EMaxOV_q<100e-3){
                FillHistograms(C[k], "COVCut_100keVee");
                FillHistograms(A[k], "COVCut_100keVee", "Al2O3");
            }
            
            if (EMaxOV_q<200e-3){
                FillHistograms(C[k], "COVCut_200keVee");
                FillHistograms(A[k], "COVCut_200keVee", "Al2O3");
            }
            
            if (EMaxOV_q<300e-3){
                FillHistograms(C[k], "COVCut_300keVee");
                FillHistograms(A[k], "COVCut_300keVee", "Al2O3");
            }
            
            if (EMaxOV_q<400e-3){
                FillHistograms(C[k], "COVCut_400keVee");
                FillHistograms(A[k], "COVCut_400keVee", "Al2O3");
            }
            
            if (EMaxOV_q<500e-3){
                FillHistograms(C[k], "COVCut_500keVee");
                FillHistograms(A[k], "COVCut_500keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<5e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_5keVee_Delayed");
                FillHistograms(A[k], "COVCut_5keVee_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<50e-3){
                FillHistograms(C[k], "COVCut_50keVee");
                FillHistograms(A[k], "COVCut_50keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<50e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_50keVee_Delayed");
                FillHistograms(A[k], "COVCut_50keVee_Delayed", "Al2O3");
            }
            
            if ((EMaxOV_q<5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "MV&COVCuts_5keVee");
                FillHistograms(A[k], "MV&COVCuts_5keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "MV&COVCuts_1keVee");
                FillHistograms(A[k], "MV&COVCuts_1keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "MV&COVCuts_Delayed");
                FillHistograms(A[k], "MV&COVCuts_Delayed", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<0.5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_0_5keVee");
                FillHistograms(A[k], "AllCuts_0_5keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<1e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_1keVee");
                FillHistograms(A[k], "AllCuts_1keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<2e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_2keVee");
                FillHistograms(A[k], "AllCuts_2keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<3e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_3keVee");
                FillHistograms(A[k], "AllCuts_3keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_5keVee");
                FillHistograms(A[k], "AllCuts_5keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<10e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_10keVee");
                FillHistograms(A[k], "AllCuts_10keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<20e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_20keVee");
                FillHistograms(A[k], "AllCuts_20keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<30e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_30keVee");
                FillHistograms(A[k], "AllCuts_30keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<50e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_50keVee");
                FillHistograms(A[k], "AllCuts_50keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<1e-3)&&(EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "AllCuts_Delayed");
                FillHistograms(C[k], "AllCuts_Delayed", "Al2O3");
            }
        }
        
        if ((NEventsM<NevM_prev)||((NEvents<Nev_prev)&&(NEventsM==0))){ //there is a reset!
            Nev+=NevM_prev*1e6+Nev_prev;
        }

        NevM_prev=NEventsM;
        Nev_prev = NEvents;
    }
    
    np->GetEntry(p_entries-1);
    Nev+=NEventsM*1e6+NEvents;
    
//    Nev=4.812e9; // for old T analyzed files
    
    Double_t EqTime = Nev/(Phi*Pane_size); //s
    double EqTimeDay = Nev/(Phi*Pane_size*24*3600.);
    std::cout<<"Equivalent time in days = "<<EqTimeDay<<std::endl;
    std::cout<<"Equivalent time in years = "<<EqTimeDay/365.<<std::endl;
    std::cout<<Nev<<std::endl;
    
    TString SimuFilename_out = SimuFilename;
    TPMERegexp("Analyzer_").Substitute(SimuFilename_out,"HistosCOVThr_");
//    TPMERegexp("FC").Substitute(SimuFilename_out,"Histos_FC");
    
    TFile* fileout = new TFile(SimuFilename_out, "RECREATE");
    
    //Write in the file the flux, num of primaries and norm
    TParameter<int>* Prim = new TParameter<int>("NPrim", Nev);
    Prim->Write();
    TParameter<float>* EqTime_s = new TParameter<float>("EqTime_s", EqTime);
    EqTime_s->Write();
    
    
    SaveHistos(fileout, EqTimeDay, Nev, NEdepRange, EdepRange, 32, CutNames, EqTimeDay, "CaWO4");
    SaveHistos(fileout, EqTimeDay, Nev, NEdepRange, EdepRange, 32, CutNames, EqTimeDay, "Al2O3");
}

void PlotEventsCryoDet_Z_T(TString SimuFilenameZ = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", TString SimuFilenameT = "Analyzer_Z_NUCLEUS_EmptySite_Step2_Ambient_U238_1-2-3.root", Double_t Phi=1.523, Double_t Pane_size=120*120,  Double_t Pane_size_simuT=120*120, Bool_t NeutronQuenched=FALSE){
    
    //Initialize histograms
    static const Int_t NEdepRange = 10;
    std::pair<TString,Double_t> EdepRange[NEdepRange] = {
        std::make_pair("eV",100),
        std::make_pair("keV",1),
        std::make_pair("keV",10),
        std::make_pair("keV",100),
        std::make_pair("MeV",1),
        std::make_pair("MeV",10),
        std::make_pair("MeV",100),
        std::make_pair("GeV",1),
        std::make_pair("GeV",10),
        std::make_pair("GeV",100)
    };
    
    TString CutNames[32]={
        "NoCut",
        "NoCut_Delayed",
        "MVCut",
        "MVCut_Delayed",
        "COVCut_0_5keVee",
        "COVCut_1keVee",
        "COVCut_10keVee",
        "COVCut_5keVee",
        "COVCut_50keVee",
        "COVCut_80keVee",
        "COVCut_100keVee",
        "COVCut_200keVee",
        "COVCut_300keVee",
        "COVCut_400keVee",
        "COVCut_500keVee",
        "COVCut_1keVee_Delayed",
        "COVCut_5keVee_Delayed",
        "COVCut_10keVee_Delayed",
        "COVCut_50keVee_Delayed",
        "MV&COVCuts_1keVee",
        "MV&COVCuts_5keVee",
        "MV&COVCuts_Delayed",
        "AllCuts_0_5keVee",
        "AllCuts_1keVee",
        "AllCuts_2keVee",
        "AllCuts_3keVee",
        "AllCuts_5keVee",
        "AllCuts_10keVee",
        "AllCuts_20keVee",
        "AllCuts_30keVee",
        "AllCuts_50keVee",
        "AllCuts_Delayed"
    };

    
    InitHistograms(NEdepRange, EdepRange, 32, CutNames, "CaWO4");
    InitHistograms(NEdepRange, EdepRange, 32, CutNames, "Al2O3");
    
    //Simulation Z ---------------------------------------------------------------
    TFile* Simuin = new TFile(SimuFilenameZ, "READ");
    TTree *np = (TTree*) Simuin->Get("Step2AlgoDir/Pri_Ntuple");
    Int_t p_entries= np->GetEntries();
    TTree *ns = (TTree*) Simuin->Get("Step2AlgoDir/Sec_Ntuple");
    Int_t s_entries= ns->GetEntries();
    
    Float_t pEDep, En, fp, sPDG, NsPidM, NsPid, EMaxOV, EAl2O3, NEvents, NEventsM, EMaxMV, EMaxIV, TgtHits, Delay, EMaxOV_q;
    Float_t C[9], A[9];
    np->SetBranchAddress(  "EMaxOV", &EMaxOV);
    np->SetBranchAddress(  "EMaxMV", &EMaxMV);
    
    np->SetBranchAddress(   "C0", &C[0]);
    np->SetBranchAddress(   "C1", &C[1]);
    np->SetBranchAddress(   "C2", &C[2]);
    np->SetBranchAddress(   "C3", &C[3]);
    np->SetBranchAddress(   "C4", &C[4]);
    np->SetBranchAddress(   "C5", &C[5]);
    np->SetBranchAddress(   "C6", &C[6]);
    np->SetBranchAddress(   "C7", &C[7]);
    np->SetBranchAddress(   "C8", &C[8]);
    
    np->SetBranchAddress(   "A0", &A[0]);
    np->SetBranchAddress(   "A1", &A[1]);
    np->SetBranchAddress(   "A2", &A[2]);
    np->SetBranchAddress(   "A3", &A[3]);
    np->SetBranchAddress(   "A4", &A[4]);
    np->SetBranchAddress(   "A5", &A[5]);
    np->SetBranchAddress(   "A6", &A[6]);
    np->SetBranchAddress(   "A7", &A[7]);
    np->SetBranchAddress(   "A8", &A[8]);
    
    np->SetBranchAddress(  "EAl2O3", &EAl2O3);
    np->SetBranchAddress(  "EMaxIV", &EMaxIV);
    np->SetBranchAddress(  "TgtHits", &TgtHits);
    np->SetBranchAddress(  "nEventM", &NEventsM);
    np->SetBranchAddress(  "nEvent", &NEvents);
    np->SetBranchAddress(  "Delay", &Delay);

    cout<<"Read simulation..."<<endl;
    
    unsigned long long Nev=0;
    unsigned long long NevM_prev=0;
    unsigned long long Nev_prev=0;

    for (int i=0; i<p_entries; i++){
        np->GetEntry(i);
        
        if (NeutronQuenched) EMaxOV_q = EMaxOV/5;
        else EMaxOV_q = EMaxOV;
        
        // loop over CaWO4 crystals:
        for (int k=0; k<9; k++){
            
            FillHistograms(C[k], "NoCut");
            FillHistograms(A[k], "NoCut", "Al2O3");
            
            if ((Delay>1)){
                FillHistograms(C[k], "NoCut_Delayed");
                FillHistograms(A[k], "NoCut_Delayed", "Al2O3");
            }
            
            if (EMaxMV<5){
                FillHistograms(C[k], "MVCut");
                FillHistograms(A[k], "MVCut", "Al2O3");
            }
            
            if ((EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "MVCut_Delayed");
                FillHistograms(A[k], "MVCut_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<0.5e-3){
                FillHistograms(C[k], "COVCut_0_5keVee");
                FillHistograms(A[k], "COVCut_0_5keVee", "Al2O3");
            }
            
            if (EMaxOV_q<1e-3){
                FillHistograms(C[k], "COVCut_1keVee");
                FillHistograms(A[k], "COVCut_1keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_1keVee_Delayed");
                FillHistograms(A[k], "COVCut_1keVee_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<10e-3){
                FillHistograms(C[k], "COVCut_10keVee");
                FillHistograms(A[k], "COVCut_10keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<10e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_10keVee_Delayed");
                FillHistograms(A[k], "COVCut_10keVee_Delayed", "Al2O3");
            }
                
            if (EMaxOV_q<5e-3){
                FillHistograms(C[k], "COVCut_5keVee");
                FillHistograms(A[k], "COVCut_5keVee", "Al2O3");
            }
            
            if (EMaxOV_q<80e-3){
                FillHistograms(C[k], "COVCut_80keVee");
                FillHistograms(A[k], "COVCut_80keVee", "Al2O3");
            }
            
            if (EMaxOV_q<100e-3){
                FillHistograms(C[k], "COVCut_100keVee");
                FillHistograms(A[k], "COVCut_100keVee", "Al2O3");
            }
            
            if (EMaxOV_q<200e-3){
                FillHistograms(C[k], "COVCut_200keVee");
                FillHistograms(A[k], "COVCut_200keVee", "Al2O3");
            }
            
            if (EMaxOV_q<300e-3){
                FillHistograms(C[k], "COVCut_300keVee");
                FillHistograms(A[k], "COVCut_300keVee", "Al2O3");
            }
            
            if (EMaxOV_q<400e-3){
                FillHistograms(C[k], "COVCut_400keVee");
                FillHistograms(A[k], "COVCut_400keVee", "Al2O3");
            }
            
            if (EMaxOV_q<500e-3){
                FillHistograms(C[k], "COVCut_500keVee");
                FillHistograms(A[k], "COVCut_500keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<5e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_5keVee_Delayed");
                FillHistograms(A[k], "COVCut_5keVee_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<50e-3){
                FillHistograms(C[k], "COVCut_50keVee");
                FillHistograms(A[k], "COVCut_50keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<50e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_50keVee_Delayed");
                FillHistograms(A[k], "COVCut_50keVee_Delayed", "Al2O3");
            }
            
            if ((EMaxOV_q<5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "MV&COVCuts_5keVee");
                FillHistograms(A[k], "MV&COVCuts_5keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "MV&COVCuts_1keVee");
                FillHistograms(A[k], "MV&COVCuts_1keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "MV&COVCuts_Delayed");
                FillHistograms(A[k], "MV&COVCuts_Delayed", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<0.5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_0_5keVee");
                FillHistograms(A[k], "AllCuts_0_5keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<1e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_1keVee");
                FillHistograms(A[k], "AllCuts_1keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<2e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_2keVee");
                FillHistograms(A[k], "AllCuts_2keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<3e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_3keVee");
                FillHistograms(A[k], "AllCuts_3keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_5keVee");
                FillHistograms(A[k], "AllCuts_5keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<10e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_10keVee");
                FillHistograms(A[k], "AllCuts_10keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<20e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_20keVee");
                FillHistograms(A[k], "AllCuts_20keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<30e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_30keVee");
                FillHistograms(A[k], "AllCuts_30keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<50e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_50keVee");
                FillHistograms(A[k], "AllCuts_50keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<1e-3)&&(EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "AllCuts_Delayed");
                FillHistograms(A[k], "AllCuts_Delayed", "Al2O3");
            }
        }
        
        if ((NEventsM<NevM_prev)||((NEvents<Nev_prev)&&(NEventsM==0))){ //there is a reset!
            Nev+=NevM_prev*1e6+Nev_prev;
        }

        NevM_prev=NEventsM;
        Nev_prev = NEvents;
    }
    
    np->GetEntry(p_entries-1);
    Nev+=NEventsM*1e6+NEvents;
    
    Simuin->Close();
    
    //Simulation T ---------------------------------------------------------------
    Simuin = new TFile(SimuFilenameT, "READ");
    np = (TTree*) Simuin->Get("NTupleAlgoDir/Pri_Ntuple");
    p_entries= np->GetEntries();
    ns = (TTree*) Simuin->Get("NTupleAlgoDir/Sec_Ntuple");
    s_entries= ns->GetEntries();
    
    np->SetBranchAddress(  "EMaxOV", &EMaxOV);
    np->SetBranchAddress(  "EMaxMV", &EMaxMV);
    
    np->SetBranchAddress(   "C0", &C[0]);
    np->SetBranchAddress(   "C1", &C[1]);
    np->SetBranchAddress(   "C2", &C[2]);
    np->SetBranchAddress(   "C3", &C[3]);
    np->SetBranchAddress(   "C4", &C[4]);
    np->SetBranchAddress(   "C5", &C[5]);
    np->SetBranchAddress(   "C6", &C[6]);
    np->SetBranchAddress(   "C7", &C[7]);
    np->SetBranchAddress(   "C8", &C[8]);
    
    np->SetBranchAddress(   "A0", &A[0]);
    np->SetBranchAddress(   "A1", &A[1]);
    np->SetBranchAddress(   "A2", &A[2]);
    np->SetBranchAddress(   "A3", &A[3]);
    np->SetBranchAddress(   "A4", &A[4]);
    np->SetBranchAddress(   "A5", &A[5]);
    np->SetBranchAddress(   "A6", &A[6]);
    np->SetBranchAddress(   "A7", &A[7]);
    np->SetBranchAddress(   "A8", &A[8]);
    
    np->SetBranchAddress(  "EAl2O3", &EAl2O3);
    np->SetBranchAddress(  "EMaxIV", &EMaxIV);
    np->SetBranchAddress(  "TgtHits", &TgtHits);
    np->SetBranchAddress(  "nEventM", &NEventsM);
    np->SetBranchAddress(  "nEvent", &NEvents);
    np->SetBranchAddress(  "Delay", &Delay);

    cout<<"Read simulation..."<<endl;

    unsigned long long Nev_T=0;
    NevM_prev=0;
    Nev_prev=0;

    for (int i=0; i<p_entries; i++){
        np->GetEntry(i);
        
        if (NeutronQuenched) EMaxOV_q = EMaxOV*5;
        else EMaxOV_q = EMaxOV;
        
        // loop over CaWO4 crystals:
        for (int k=0; k<9; k++){
            
            FillHistograms(C[k], "NoCut");
            FillHistograms(A[k], "NoCut", "Al2O3");
            
            if ((Delay>1)){
                FillHistograms(C[k], "NoCut_Delayed");
                FillHistograms(A[k], "NoCut_Delayed", "Al2O3");
            }
            
            if (EMaxMV<5){
                FillHistograms(C[k], "MVCut");
                FillHistograms(A[k], "MVCut", "Al2O3");
            }
            
            if ((EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "MVCut_Delayed");
                FillHistograms(A[k], "MVCut_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<0.5e-3){
                FillHistograms(C[k], "COVCut_0_5keVee");
                FillHistograms(A[k], "COVCut_0_5keVee", "Al2O3");
            }
            
            if (EMaxOV_q<1e-3){
                FillHistograms(C[k], "COVCut_1keVee");
                FillHistograms(A[k], "COVCut_1keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_1keVee_Delayed");
                FillHistograms(A[k], "COVCut_1keVee_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<10e-3){
                FillHistograms(C[k], "COVCut_10keVee");
                FillHistograms(A[k], "COVCut_10keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<10e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_10keVee_Delayed");
                FillHistograms(A[k], "COVCut_10keVee_Delayed", "Al2O3");
            }
                
            if (EMaxOV_q<5e-3){
                FillHistograms(C[k], "COVCut_5keVee");
                FillHistograms(A[k], "COVCut_5keVee", "Al2O3");
            }
            
            if (EMaxOV_q<80e-3){
                FillHistograms(C[k], "COVCut_80keVee");
                FillHistograms(A[k], "COVCut_80keVee", "Al2O3");
            }
            
            if (EMaxOV_q<100e-3){
                FillHistograms(C[k], "COVCut_100keVee");
                FillHistograms(A[k], "COVCut_100keVee", "Al2O3");
            }
            
            if (EMaxOV_q<200e-3){
                FillHistograms(C[k], "COVCut_200keVee");
                FillHistograms(A[k], "COVCut_200keVee", "Al2O3");
            }
            
            if (EMaxOV_q<300e-3){
                FillHistograms(C[k], "COVCut_300keVee");
                FillHistograms(A[k], "COVCut_300keVee", "Al2O3");
            }
            
            if (EMaxOV_q<400e-3){
                FillHistograms(C[k], "COVCut_400keVee");
                FillHistograms(A[k], "COVCut_400keVee", "Al2O3");
            }
            
            if (EMaxOV_q<500e-3){
                FillHistograms(C[k], "COVCut_500keVee");
                FillHistograms(A[k], "COVCut_500keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<5e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_5keVee_Delayed");
                FillHistograms(A[k], "COVCut_5keVee_Delayed", "Al2O3");
            }
            
            if (EMaxOV_q<50e-3){
                FillHistograms(C[k], "COVCut_50keVee");
                FillHistograms(A[k], "COVCut_50keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<50e-3)&&(Delay>1)){
                FillHistograms(C[k], "COVCut_50keVee_Delayed");
                FillHistograms(A[k], "COVCut_50keVee_Delayed", "Al2O3");
            }
            
            if ((EMaxOV_q<5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "MV&COVCuts_5keVee");
                FillHistograms(A[k], "MV&COVCuts_5keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "MV&COVCuts_1keVee");
                FillHistograms(A[k], "MV&COVCuts_1keVee", "Al2O3");
            }
            
            if ((EMaxOV_q<1e-3)&&(EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "MV&COVCuts_Delayed");
                FillHistograms(A[k], "MV&COVCuts_Delayed", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<0.5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_0_5keVee");
                FillHistograms(A[k], "AllCuts_0_5keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<1e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_1keVee");
                FillHistograms(A[k], "AllCuts_1keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<2e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_2keVee");
                FillHistograms(A[k], "AllCuts_2keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<3e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_3keVee");
                FillHistograms(A[k], "AllCuts_3keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<5e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_5keVee");
                FillHistograms(A[k], "AllCuts_5keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<10e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_10keVee");
                FillHistograms(A[k], "AllCuts_10keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<20e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_20keVee");
                FillHistograms(A[k], "AllCuts_20keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<30e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_30keVee");
                FillHistograms(A[k], "AllCuts_30keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<50e-3)&&(EMaxMV<5)){
                FillHistograms(C[k], "AllCuts_50keVee");
                FillHistograms(A[k], "AllCuts_50keVee", "Al2O3");
            }
            
            if ((TgtHits==1)&&(EMaxIV<30e-6)&&(EMaxOV_q<1e-3)&&(EMaxMV<5)&&(Delay>1)){
                FillHistograms(C[k], "AllCuts_Delayed");
                FillHistograms(C[k], "AllCuts_Delayed", "Al2O3");
            }
        }
        
        if ((NEventsM<NevM_prev)||((NEvents<Nev_prev)&&(NEventsM==0))){ //there is a reset!
            Nev_T+=NevM_prev*1e6+Nev_prev;
        }

        NevM_prev=NEventsM;
        Nev_prev = NEvents;
    }
    
    np->GetEntry(p_entries-1);
    Nev_T+=NEventsM*1e6+NEvents;
    
    //Nev=4.812e9; // for old T analyzed files
    
    Double_t EqTime = Nev/(Phi*Pane_size)+Nev_T/(Phi*Pane_size_simuT); //s
    double EqTimeDay = Nev/(Phi*Pane_size*24*3600.)+Nev_T/(Phi*Pane_size_simuT*24*3600.);
    std::cout<<"Equivalent time in days = "<<EqTimeDay<<std::endl;
    std::cout<<"Equivalent time in years = "<<EqTimeDay/365.<<std::endl;
    std::cout<<Nev<<std::endl;
    
    TString SimuFilename_out = SimuFilenameT;
    TPMERegexp("Analyzer_T").Substitute(SimuFilename_out,"HistosCOVThr_ZT_");
//    TPMERegexp("FC").Substitute(SimuFilename_out,"Histos_FC");
    
    std::cout<<SimuFilename_out<<std::endl;
    TFile* fileout = new TFile(SimuFilename_out, "RECREATE");
    
    //Write in the file the flux, num of primaries and norm
    TParameter<int>* Prim = new TParameter<int>("NPrim", Nev+Nev_T);
    Prim->Write();
    TParameter<float>* EqTime_s = new TParameter<float>("EqTime_s", EqTime);
    EqTime_s->Write();
    
    SaveHistos(fileout, EqTimeDay, Nev+Nev_T, NEdepRange, EdepRange, 32, CutNames, EqTimeDay, "CaWO4");
    SaveHistos(fileout, EqTimeDay, Nev+Nev_T, NEdepRange, EdepRange, 32, CutNames, EqTimeDay, "Al2O3");
    fileout->Close();
}
