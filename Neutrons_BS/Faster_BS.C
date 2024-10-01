// Analysis macro for Bonner Sphere ata taken with FASTER ADC (Saclay)
// C.Goupy 11/05/2023


void Faster_BS(TString dir = "../NeutronMeasurement/Faster@Saclay/", TString filelist = "datafiles.txt"){
    
    Int_t nrj;
    ULong64_t timestamp;
    UShort_t label;
    
    Double_t timestamp_sec = 1e-6; //usec
    Double_t Prev_muonTime=0; //usec
    Int_t Muon_trig_index=0;
    Double_t Event_time=0; //usec
    Double_t Coinc_window=100; //usec
    Int_t threshold = 20; //(mV)
    
    //Calibration
    Double_t Calib_Neutr = 67532./1.03; //ADC/MeV (from muons in LiI detector)
    Double_t Calib_NeutrPeak = 200000/2.9; //ADC/MeV (from NP in LiI detector)
    
    Double_t Calib_Muon = 163633/7.5; //ADC/MeV (from muons in MV detector)
    
    TH1D* Muon_Histogram = new TH1D("Muon_Histogram", "Muon_Histogram", 2e2, 0, 20);
    TH1D* Neutron_Histogram = new TH1D("Neutron_Histogram", "Neutron_Histogram", 1e2, 0, 10);
    TH1D* NeutronMV_Histogram = new TH1D("NeutronMV_Histogram", "NeutronMV_Histogram", 1e2, 0, 10);
    TH1D* VetoedMuons_MV_Histogram = new TH1D("VetoedMuons_MV_Histogram", "VetoedMuons_MV_Histogram", 2e2, 0, 20);
    TH1D* VetoedMuons_BS_Histogram = new TH1D("VetoedMuons_BS_Histogram", "VetoedMuons_BS_Histogram", 1e2, 0, 10);
    
    ifstream sfiles;
    TString thisfile;
    sfiles.open(dir+filelist);
    if(!sfiles.is_open()) { cout << "Error opening text file" << endl;}
    
    Double_t Tot_time_sec=0;
    
    TFile* filein;
    TTree* DatTree;
    
    while (sfiles>>thisfile){
        filein = new TFile(dir+thisfile, "READ");
        DatTree = (TTree*) filein->Get("DataTree");
        
    //    DatTree->Draw("nrj>>Hist_Neutr(1E3, 0, 300000)", "(label==1)");
        int NEntries = DatTree->GetEntries();
        
        DatTree->SetBranchAddress("nrj",&nrj);
        DatTree->SetBranchAddress("time",&timestamp);
        DatTree->SetBranchAddress("label",&label);
        
        for (int i=0; i<NEntries; i++){
            DatTree->GetEntry(i);
            
            if (label==2){ //muon event
                if (nrj>threshold){
                    Muon_Histogram->Fill(nrj/Calib_Muon);
                    Prev_muonTime = timestamp*timestamp_sec;
                    Muon_trig_index=i;
                }
            }
            
            if (label==1){ //neutron event
                Event_time = timestamp*timestamp_sec;
                Neutron_Histogram->Fill(nrj/Calib_Neutr);
                if ((Event_time<Prev_muonTime-0.5*Coinc_window)||(Event_time>Prev_muonTime+0.5*Coinc_window)){ //Not vetoed
                    NeutronMV_Histogram->Fill(nrj/Calib_Neutr);
                }
                else {
                    VetoedMuons_BS_Histogram->Fill(nrj/Calib_Neutr);
                    DatTree->GetEntry(Muon_trig_index);
                    VetoedMuons_MV_Histogram->Fill(nrj/Calib_Muon);
                }
            }
        }

        
        DatTree->GetEntry(NEntries-1);
        if (Tot_time_sec<timestamp_sec*1e-6*timestamp){
            Tot_time_sec = timestamp_sec*1e-6*timestamp;
            cout<<"this file tot time = "<<Tot_time_sec/(3600.)<<" h"<<endl;
        }
    }
    
    //Scale histograms in ev/sec
    Muon_Histogram->Scale(1./(Tot_time_sec*Muon_Histogram->GetBinWidth(1)));
    Neutron_Histogram->Scale(1./(Tot_time_sec*Muon_Histogram->GetBinWidth(1)));
    NeutronMV_Histogram->Scale(1./(Tot_time_sec*Muon_Histogram->GetBinWidth(1)));
    VetoedMuons_BS_Histogram->Scale(1./(Tot_time_sec*Muon_Histogram->GetBinWidth(1)));
    VetoedMuons_MV_Histogram->Scale(1./(Tot_time_sec*Muon_Histogram->GetBinWidth(1)));
    
    //Legends
    Muon_Histogram->GetXaxis()->SetTitle("Energy [MeV]");
    Neutron_Histogram->GetXaxis()->SetTitle("Energy [MeV]");
    NeutronMV_Histogram->GetXaxis()->SetTitle("Energy [MeV]");
    VetoedMuons_BS_Histogram->GetXaxis()->SetTitle("Energy [MeV]");
    VetoedMuons_MV_Histogram->GetXaxis()->SetTitle("Energy [MeV]");
    
    Muon_Histogram->GetYaxis()->SetTitle("Differential Rate [cps/keV]");
    Neutron_Histogram->GetYaxis()->SetTitle("Differential Rate [cps/keV]");
    NeutronMV_Histogram->GetYaxis()->SetTitle("Differential Rate [cps/keV]");
    VetoedMuons_BS_Histogram->GetYaxis()->SetTitle("Differential Rate [cps/keV]");
    VetoedMuons_MV_Histogram->GetYaxis()->SetTitle("Differential Rate [cps/keV]");
    
    TString Outfilename = dir+"20230502_LiI_trig-20mV_BotMV_trig-20mV_HV-1600V_histograms.root";
    TFile* fileout = new TFile(Outfilename, "RECREATE");
    
    TCanvas* cMuonPan = new TCanvas("MuonPanel", "MuonPanel", 1200, 1200);
    cMuonPan->Divide(1,2);
    cMuonPan->cd(1);
    Muon_Histogram->Draw();
    Muon_Histogram->Write();
    cMuonPan->cd(2);
    VetoedMuons_MV_Histogram->Draw();
    VetoedMuons_MV_Histogram->Write();
    
    TCanvas* cBS = new TCanvas("BS", "BS", 1200, 1200);
    cBS->Divide(2,2);
    cBS->cd(1);
    Neutron_Histogram->Draw();
    Neutron_Histogram->Write();
    cBS->cd(2);
    NeutronMV_Histogram->Draw();
    NeutronMV_Histogram->Write();
    cBS->cd(3);
    VetoedMuons_BS_Histogram->Draw();
    VetoedMuons_BS_Histogram->Write();
    cBS->cd(4);
    Neutron_Histogram->Draw();
    NeutronMV_Histogram->Draw("SAME");
    NeutronMV_Histogram->SetLineColor(kRed);
}
void Histo_Simu_Neutrons(TString NeutronsSimuFilename = "../Simulations/Cf252/Analyzer_PT_BS_NewGeom_8inch_SaclayCf252_2.root"){
    //Neutrons' simulations
    TFile* NeutronsSimuin = new TFile(NeutronsSimuFilename, "READ");
    TTree *np = (TTree*) NeutronsSimuin->Get("NTupleAlgoDir/Pri_Ntuple");
    Int_t p_entries= np->GetEntries();
    TTree *ns = (TTree*) NeutronsSimuin->Get("NTupleAlgoDir/Sec_Ntuple");
    Int_t s_entries= ns->GetEntries();
    
    Float_t pEDep, En, fp, sPDG, NsPidM, NsPid;
    np->SetBranchAddress(  "EDep", &pEDep);
    ns->SetBranchAddress(   "PDG", &sPDG);
    ns->SetBranchAddress(  "ProcId", &fp);
    ns->SetBranchAddress(  "PidM", &NsPidM);
    ns->SetBranchAddress(   "Pid", &NsPid);
    
    TRandom gen(0);
    
    TH1D* NeutronsSimu_Neutron_Histogram = new TH1D("NeutronsSimu_Neutron_Histogram", "NeutronsSimu_Neutron_Histogram", 100000, 0, 10000);
    
    cout<<"Read Neutrons simulation..."<<endl;
    
    long Nev = 1000000000; //TUM: 1000000000 //Saclay: 4833333333 (ccin2p3, old geom source), 100000000 (local), 1000000000(ccin2p3 source new geometry)
    Double_t Phi = 500000;//ev/s //TUM:12300 Saclay:500000
    Double_t S = 1; //cm2
     
    int isQF;
    int sk=0;
    int n_from_final_process=0;
    int n_fromEdep=0;
    
    
    for (int i=0; i<p_entries; i++){
        if (i%(p_entries/100) ==0) cout<<1+int((1.*i)/p_entries*100)<<"% of events processed\n";
        np->GetEntry(i);
        En = pEDep;        // MeV
        
        //Quenching
        if (En>4.7 && En<4.8) {
            //isQF=1;
            n_fromEdep++;
        }
        
        ns->GetEntry(sk);
        // Get the "id" of the primary particles associated to the secondary particle sk
        int id = ((Int_t) (NsPidM+0.5))*1000000 + ((Int_t) (NsPid+0.5));
        while (id==i && sk < s_entries){
            if ((sPDG==2112 || sPDG==-2112) && pEDep>4){
                n_from_final_process++;
                isQF=1;
                if (En<4.7 || En>4.8) cout<<pEDep<<endl;
            }
            sk++;
            ns->GetEntry(sk);
            id = ((Int_t) (NsPidM+0.5))*1000000 + ((Int_t) (NsPid+0.5));
        }
            
        if (isQF==1){
            En = 0.6*En;
        }
    
        //Resolution LiI(Eu)
        En = En*(1+4./sqrt(En*1e3)*gen.Gaus()); //sigma = 4*sqrt(E)
        
        NeutronsSimu_Neutron_Histogram->Fill(En*1e3);
    }
    
    cout<<"n_fromEdep = "<<n_fromEdep<<endl;
    cout<<"n_from_final_process = "<<n_from_final_process<<endl;

    TString NeutronsSimuFilename_out = NeutronsSimuFilename;
    TPMERegexp("Analyzer_").Substitute(NeutronsSimuFilename_out,"Histos2_");
    TFile* fileout_neutr = new TFile(NeutronsSimuFilename_out, "RECREATE");
    cout<<NeutronsSimuFilename_out<<endl;
    
    Double_t EqTime = Nev/(Phi*S); //s
    
    double NeutrEqTime = Nev/(Phi*S);
    
    NeutronsSimu_Neutron_Histogram->Scale(1./(NeutrEqTime*NeutronsSimu_Neutron_Histogram->GetBinWidth(1)));
    NeutronsSimu_Neutron_Histogram->SetLineColor(kMagenta);
    NeutronsSimu_Neutron_Histogram->Write("NeutronsSimu_Neutron_Histogram");
    //---
}

void Faster_BS_Simu(TString NeutronsSimuFilename = "../Simulations/Neutrons/Analyzer_PT_BS_NewGeom_8inch_FullMV_AtmNeutrons.000-013.root", TString MuonsSimuFilename = "../Simulations/Muons/CCIN2P3/Analyzer_PT_BS_NewGeom_8inch_FullMV_AtmMuons.000-039.root", TString GammasSimuFilename = "../Simulations/AmbGammas/Analyzer_BT_BS_8inch_FullMV_AmbGammas_2.000-142.root"){
    //Simulation
    
    //Neutrons' simulations
    TFile* NeutronsSimuin = new TFile(NeutronsSimuFilename);
    TTree *np = (TTree*) NeutronsSimuin->Get("NTupleAlgoDir/Pri_Ntuple");
    Int_t p_entries= np->GetEntries();
    
    Float_t pEDep, pEMV0, pEMV1, pEMV2, pEMV3, En, EMV, resolEMV;
    np->SetBranchAddress(  "EDep", &pEDep);
    np->SetBranchAddress("EMV0", &pEMV0);
    
    TRandom gen(0);
    
    TH1D* NeutronsSimu_Neutron_Histogram = new TH1D("NeutronsSimu_Neutron_Histogram", "NeutronsSimu_Neutron_Histogram", 1e2, 0, 10);
    TH1D* NeutronsSimu_Muon_Histogram = new TH1D("NeutronsSimu_Muon_Histogram", "NeutronsSimu_Muon_Histogram", 2e2, 0, 20);
    TH1D* NeutronsSimu_NeutronMV_Histogram = new TH1D("NeutronsSimu_NeutronMV_Histogram", "NeutronsSimu_NeutronMV_Histogram", 1e2, 0, 10);
    TH1D* NeutronsSimu_VetoedMuons_MV_Histogram = new TH1D("NeutronsSimu_VetoedMuons_MV_Histogram", "NeutronsSimu_VetoedMuons_MV_Histogram", 2e2, 0, 20);
    TH1D* NeutronsSimu_VetoedMuons_BS_Histogram = new TH1D("NeutronsSimu_VetoedMuons_BS_Histogram", "NeutronsSimu_VetoedMuons_BS_Histogram", 2e2, 0, 20);
    
    cout<<"Read Neutrons simulation..."<<endl;
    
    int Nev = 1000000000;
    Double_t Phi = 0.0134;//ev/cm2/s
    Double_t S = 40*40; //cm2
     
    int isQF;
    for (int i=0; i<p_entries; i++){
        
        if (i%10000000 ==0) cout<<1+int((1.*i)/Nev*100)<<"% of events processed\n";
        np->GetEntry(i);
        En = pEDep;        // MeV
        EMV = pEMV0;       // MeV
        if (EMV>0) resolEMV = EMV*(1+0.75/sqrt(EMV)*gen.Gaus());
        else resolEMV = 0.;
        
        //Quenching
        if (En>4.7 && En<4.8) isQF=1;
        if (isQF==1){
            En = 0.6*En;
        }
        
        //Resolution LiI(Eu)
        En = En*(1+4./sqrt(En)*gen.Gaus()); //sigma = 4*sqrt(E)
        
//        cout<<"no resol ="<< EMV<<endl;
//        cout<<" with resol = "<< resolEMV<<endl;
        if (resolEMV>1. && En>0.2){ //Vetoed
            NeutronsSimu_VetoedMuons_MV_Histogram->Fill(resolEMV);
            NeutronsSimu_VetoedMuons_BS_Histogram->Fill(En);
        }
        if (resolEMV<=1. && En>0.2){ //Not vetoed
            NeutronsSimu_NeutronMV_Histogram->Fill(En);
        }
        NeutronsSimu_Neutron_Histogram->Fill(En);
        NeutronsSimu_Muon_Histogram->Fill(resolEMV);
    }

    TString NeutronsSimuFilename_out = NeutronsSimuFilename;
    TPMERegexp("Analyzer_").Substitute(NeutronsSimuFilename_out,"/Histos_");
    TFile* fileout_neutr = new TFile(NeutronsSimuFilename_out, "RECREATE");
    cout<<NeutronsSimuFilename_out<<endl;
    
    Double_t EqTime = Nev/(Phi*S); //s
    
    double NeutrEqTime = Nev/(Phi*S);
    
    NeutronsSimu_Neutron_Histogram->Scale(1./(NeutrEqTime*NeutronsSimu_Neutron_Histogram->GetBinWidth(1)));
    NeutronsSimu_Neutron_Histogram->SetLineColor(kMagenta);
    NeutronsSimu_Neutron_Histogram->Write("NeutronsSimu_Neutron_Histogram");
    
    NeutronsSimu_Muon_Histogram->Scale(1./(NeutrEqTime*NeutronsSimu_Muon_Histogram->GetBinWidth(1)));
    NeutronsSimu_Muon_Histogram->SetLineColor(kMagenta);
    NeutronsSimu_Muon_Histogram->Write("NeutronsSimu_Muon_Histogram");
    
    NeutronsSimu_NeutronMV_Histogram->Scale(1./(NeutrEqTime*NeutronsSimu_NeutronMV_Histogram->GetBinWidth(1)));
    NeutronsSimu_NeutronMV_Histogram->SetLineColor(kMagenta);
    NeutronsSimu_NeutronMV_Histogram->Write("NeutronsSimu_NeutronMV_Histogram");
    
    NeutronsSimu_VetoedMuons_BS_Histogram->Scale(1./(NeutrEqTime*NeutronsSimu_VetoedMuons_BS_Histogram->GetBinWidth(1)));
    NeutronsSimu_VetoedMuons_BS_Histogram->SetLineColor(kMagenta);
    NeutronsSimu_VetoedMuons_BS_Histogram->Write("NeutronsSimu_VetoedMuons_BS_Histogram");
    
    NeutronsSimu_VetoedMuons_MV_Histogram->Scale(1./(NeutrEqTime*NeutronsSimu_VetoedMuons_MV_Histogram->GetBinWidth(1)));
    NeutronsSimu_VetoedMuons_MV_Histogram->SetLineColor(kMagenta);
    NeutronsSimu_VetoedMuons_MV_Histogram->Write("NeutronsSimu_VetoedMuons_MV_Histogram");
    //---
    
    //Muons' simulations
    TFile* MuonsSimuin = new TFile(MuonsSimuFilename);
    TTree *mp = (TTree*) MuonsSimuin->Get("NTupleAlgoDir/Pri_Ntuple");
    p_entries= mp->GetEntries();
    
    mp->SetBranchAddress(  "EDep", &pEDep);
    mp->SetBranchAddress("EMV0", &pEMV0);
    
    TH1D* MuonSimu_Neutron_Histogram = new TH1D("MuonSimu_Neutron_Histogram", "MuonSimu_Neutron_Histogram", 1e2, 0, 10);
    TH1D* MuonSimu_Muon_Histogram = new TH1D("MuonSimu_Muon_Histogram", "MuonSimu_Muon_Histogram", 2e2, 0, 20);
    TH1D* MuonSimu_NeutronMV_Histogram = new TH1D("MuonSimu_NeutronMV_Histogram", "MuonSimu_NeutronMV_Histogram", 1e2, 0, 10);
    TH1D* MuonSimu_VetoedMuons_MV_Histogram = new TH1D("MuonSimu_VetoedMuons_MV_Histogram", "MuonSimu_VetoedMuons_MV_Histogram", 2e2, 0, 20);
    TH1D* MuonSimu_VetoedMuons_BS_Histogram = new TH1D("MuonSimu_VetoedMuons_BS_Histogram", "MuonSimu_VetoedMuons_BS_Histogram", 2e2, 0, 20);
    
    cout<<"Read Muon simulation..."<<endl;
    
    int Nev_mu = 100000000;
    Double_t Phi_mu = 0.019;//ev/cm2/s
    
    for (int i=0; i<p_entries; i++){
        
        if (i%10000000 ==0) cout<<1+int((1.*i)/Nev_mu*100)<<"% of events processed\n";
        mp->GetEntry(i);
        En = pEDep;        // MeV
        EMV = pEMV0;       // MeV
        if (EMV>0) resolEMV = EMV*(1+0.75/sqrt(EMV)*gen.Gaus());
        else resolEMV = 0.;
//        cout<<"no resol ="<< EMV<<endl;
//        cout<<" with resol = "<< resolEMV<<endl;
        
        //Resolution LiI(Eu)
        En = En*(1+4./sqrt(En)*gen.Gaus()); //sigma = 4*sqrt(E)
        
        if (resolEMV>1. && En>0.2){ //Vetoed
            MuonSimu_VetoedMuons_MV_Histogram->Fill(resolEMV);
            MuonSimu_VetoedMuons_BS_Histogram->Fill(En);
        }
        if (resolEMV<=1. && En>0.2){ //Not vetoed
            MuonSimu_NeutronMV_Histogram->Fill(En);
        }
        MuonSimu_Neutron_Histogram->Fill(En);
        MuonSimu_Muon_Histogram->Fill(resolEMV);
    }

    TString MuonsSimuFilename_out = MuonsSimuFilename;
    TPMERegexp("Analyzer_").Substitute(MuonsSimuFilename_out,"/Histos_");
    TFile* fileout = new TFile(MuonsSimuFilename_out, "RECREATE");
    cout<<MuonsSimuFilename_out<<endl;

    double MuEqTime = Nev_mu/(Phi_mu*S);
    
    MuonSimu_Neutron_Histogram->Scale(1./(MuEqTime*MuonSimu_Neutron_Histogram->GetBinWidth(1)));
    MuonSimu_Neutron_Histogram->SetLineColor(kMagenta);
    MuonSimu_Neutron_Histogram->Write("MuonSimu_Neutron_Histogram");
    
    MuonSimu_Muon_Histogram->Scale(1./(MuEqTime*MuonSimu_Muon_Histogram->GetBinWidth(1)));
    MuonSimu_Muon_Histogram->SetLineColor(kMagenta);
    MuonSimu_Muon_Histogram->Write("MuonSimu_Muon_Histogram");
    
    MuonSimu_NeutronMV_Histogram->Scale(1./(MuEqTime*MuonSimu_NeutronMV_Histogram->GetBinWidth(1)));
    MuonSimu_NeutronMV_Histogram->SetLineColor(kMagenta);
    MuonSimu_NeutronMV_Histogram->Write("MuonSimu_NeutronMV_Histogram");
    
    MuonSimu_VetoedMuons_BS_Histogram->Scale(1./(MuEqTime*MuonSimu_VetoedMuons_BS_Histogram->GetBinWidth(1)));
    MuonSimu_VetoedMuons_BS_Histogram->SetLineColor(kMagenta);
    MuonSimu_VetoedMuons_BS_Histogram->Write("MuonSimu_VetoedMuons_BS_Histogram");
    
    MuonSimu_VetoedMuons_MV_Histogram->Scale(1./(MuEqTime*MuonSimu_VetoedMuons_MV_Histogram->GetBinWidth(1)));
    MuonSimu_VetoedMuons_MV_Histogram->SetLineColor(kMagenta);
    MuonSimu_VetoedMuons_MV_Histogram->Write("MuonSimu_VetoedMuons_MV_Histogram");
    //---
    
    //Gammas' simulations
    TFile* GammasSimuin = new TFile(GammasSimuFilename);
    TTree *gamp = (TTree*) GammasSimuin->Get("NTupleAlgoDir/Pri_Ntuple");
    Int_t gamp_entries= gamp->GetEntries();
    
    gamp->SetBranchAddress(  "EDep", &pEDep);
    gamp->SetBranchAddress("EMV0", &pEMV0);
    
    int Nev_gamma =1000000000; //100000040; //1000000000;
    Double_t PhiGam = 10.;//ev/cm2/s
    
    TH1D* GamSimu_Neutron_Histogram = new TH1D("GamSimu_Neutron_Histogram", "GamSimu_Neutron_Histogram", 1e2, 0, 10);
    TH1D* GamSimu_Muon_Histogram = new TH1D("GamSimu_Muon_Histogram", "GamSimu_Muon_Histogram", 2e2, 0, 20);
    TH1D* GamSimu_NeutronMV_Histogram = new TH1D("GamSimu_NeutronMV_Histogram", "GamSimu_NeutronMV_Histogram", 1e2, 0, 10);
    TH1D* GamSimu_VetoedMuons_MV_Histogram = new TH1D("GamSimu_VetoedMuons_MV_Histogram", "GamSimu_VetoedMuons_MV_Histogram", 2e2, 0, 20);
    TH1D* GamSimu_VetoedMuons_BS_Histogram = new TH1D("GamSimu_VetoedMuons_BS_Histogram", "GamSimu_VetoedMuons_BS_Histogram", 2e2, 0, 20);
    
    cout<<"Read gamma simulation..."<<endl;
    for (int j=0; j<gamp_entries; j++){
        
        if (j%10000000 ==0) cout<<1+int((1.*j)/Nev_gamma*100)<<"% of events processed\n";
        
        gamp->GetEntry(j);
        En = pEDep;        // MeV
        EMV = pEMV0;       // MeV
        if (EMV>0) resolEMV = EMV*(1+0.75/sqrt(EMV)*gen.Gaus());
        else resolEMV = 0.;
        
        //Resolution LiI(Eu)
        En = En*(1+4./sqrt(En)*gen.Gaus()); //sigma = 4*sqrt(E)
        
        if (resolEMV>1. && En>0.2){ //Vetoed
            GamSimu_VetoedMuons_MV_Histogram->Fill(resolEMV);
            GamSimu_VetoedMuons_BS_Histogram->Fill(En);
        }
        if (resolEMV<=1. && En>0.2){ //Not vetoed
            GamSimu_NeutronMV_Histogram->Fill(En);
        }
        GamSimu_Neutron_Histogram->Fill(En);
        GamSimu_Muon_Histogram->Fill(resolEMV);
    }

    double GamEqTime = Nev_gamma/(PhiGam*S);
    
    TString GammasSimuFilename_out = GammasSimuFilename;
    TPMERegexp("Analyzer_").Substitute(GammasSimuFilename_out,"/Histos_");
    TFile* GammasSimuout = new TFile(GammasSimuFilename_out, "RECREATE");
    
    GamSimu_Neutron_Histogram->Scale(1./(GamEqTime*GamSimu_Neutron_Histogram->GetBinWidth(1)));
    GamSimu_Neutron_Histogram->SetLineColor(kGreen);
    GamSimu_Neutron_Histogram->Write("GamSimu_Neutron_Histogram");
    
    GamSimu_Muon_Histogram->Scale(1./(GamEqTime*GamSimu_Muon_Histogram->GetBinWidth(1)));
    GamSimu_Muon_Histogram->SetLineColor(kGreen);
    GamSimu_Muon_Histogram->Write("GamSimu_Muon_Histogram");
    
    GamSimu_NeutronMV_Histogram->Scale(1./(GamEqTime*GamSimu_NeutronMV_Histogram->GetBinWidth(1)));
    GamSimu_NeutronMV_Histogram->SetLineColor(kGreen);
    GamSimu_NeutronMV_Histogram->Write("GamSimu_NeutronMV_Histogram");
    
    GamSimu_VetoedMuons_BS_Histogram->Scale(1./(GamEqTime*GamSimu_VetoedMuons_BS_Histogram->GetBinWidth(1)));
    GamSimu_VetoedMuons_BS_Histogram->SetLineColor(kGreen);
    GamSimu_VetoedMuons_BS_Histogram->Write("GamSimu_VetoedMuons_BS_Histogram");
    GamSimu_VetoedMuons_MV_Histogram->Scale(1./(GamEqTime*GamSimu_VetoedMuons_MV_Histogram->GetBinWidth(1)));
    GamSimu_VetoedMuons_MV_Histogram->SetLineColor(kGreen);
    GamSimu_VetoedMuons_MV_Histogram->Write("GamSimu_VetoedMuons_MV_Histogram");
    //---
    
    fileout->Close();
    GammasSimuin->Close();
    GammasSimuout->Close();
}

void Show_BS_Simu_Data(TString filename = "../NeutronMeasurement/Faster@Saclay/20230502_LiI_trig-20mV_BotMV_trig-20mV_HV-1600V_histograms.root",TString NeutronsSimuFilename = "../Simulations/Neutrons/Histos_PT_BS_NewGeom_8inch_FullMV_AtmNeutrons.000-013.root", TString MuonsSimuFilename = "../Simulations/Muons/CCIN2P3/Histos_PT_BS_NewGeom_8inch_FullMV_AtmMuons.000-039.root", TString GammasSimuFilename = "../Simulations/AmbGammas/Histos_PT_BS_NewGeom_8inch_FullMV_AmbGammas.000-142.root"){
    
    //-------- Data
    cout<<"Data"<<endl;
    TFile* filein_data = new TFile(filename, "READ");
    
    TH1D* Muon_Histogram = (TH1D*) filein_data->Get("Muon_Histogram");
    TH1D* Neutron_Histogram = (TH1D*) filein_data->Get("Neutron_Histogram");
    TH1D* NeutronMV_Histogram = (TH1D*) filein_data->Get("NeutronMV_Histogram");
    TH1D* VetoedMuons_MV_Histogram = (TH1D*) filein_data->Get("VetoedMuons_MV_Histogram");
    TH1D* VetoedMuons_BS_Histogram = (TH1D*) filein_data->Get("VetoedMuons_BS_Histogram");
    
    //-------- Simulation Neutrons
    cout<<"Neutrons"<<endl;
    TFile* filein_neutrons = new TFile(NeutronsSimuFilename, "READ");
    
    TH1D* NeutronsSimu_Neutron_Histogram = (TH1D*) filein_neutrons->Get("NeutronsSimu_Neutron_Histogram");
    TH1D* NeutronsSimu_Muon_Histogram = (TH1D*) filein_neutrons->Get("NeutronsSimu_Muon_Histogram");
    TH1D* NeutronsSimu_NeutronMV_Histogram = (TH1D*) filein_neutrons->Get("NeutronsSimu_NeutronMV_Histogram");
    TH1D* NeutronsSimu_VetoedMuons_BS_Histogram = (TH1D*) filein_neutrons->Get("NeutronsSimu_VetoedMuons_BS_Histogram");
    TH1D* NeutronsSimu_VetoedMuons_MV_Histogram = (TH1D*) filein_neutrons->Get("NeutronsSimu_VetoedMuons_MV_Histogram");
    
    //-------- Simulation Gammas
    cout<<"Gammas"<<endl;
    TFile* filein_gammas = new TFile(GammasSimuFilename, "READ");
    
    TH1D* GamSimu_Neutron_Histogram = (TH1D*) filein_gammas->Get("GamSimu_Neutron_Histogram");
    TH1D* GamSimu_Muon_Histogram = (TH1D*) filein_gammas->Get("GamSimu_Muon_Histogram");
    TH1D* GamSimu_NeutronMV_Histogram = (TH1D*) filein_gammas->Get("GamSimu_NeutronMV_Histogram");
    TH1D* GamSimu_VetoedMuons_BS_Histogram = (TH1D*) filein_gammas->Get("GamSimu_VetoedMuons_BS_Histogram");
    TH1D* GamSimu_VetoedMuons_MV_Histogram = (TH1D*) filein_gammas->Get("GamSimu_VetoedMuons_MV_Histogram");
    
    //-------- Simulation Muons
    cout<<"Muons"<<endl;
    TFile* filein_muons = new TFile(MuonsSimuFilename, "READ");
    
    TH1D* MuonSimu_Neutron_Histogram = (TH1D*) filein_muons->Get("MuonSimu_Neutron_Histogram");
    TH1D* MuonSimu_Muon_Histogram = (TH1D*) filein_muons->Get("MuonSimu_Muon_Histogram");
    TH1D* MuonSimu_NeutronMV_Histogram = (TH1D*) filein_muons->Get("MuonSimu_NeutronMV_Histogram");
    TH1D* MuonSimu_VetoedMuons_BS_Histogram = (TH1D*) filein_muons->Get("MuonSimu_VetoedMuons_BS_Histogram");
    TH1D* MuonSimu_VetoedMuons_MV_Histogram = (TH1D*) filein_muons->Get("MuonSimu_VetoedMuons_MV_Histogram");
    
    //-----
    
    TH1D* SumSimu_Neutron_Histogram = (TH1D*) GamSimu_Neutron_Histogram->Clone();
    TH1D* SumSimu_Muon_Histogram = (TH1D*) GamSimu_Muon_Histogram->Clone();
    TH1D* SumSimu_NeutronMV_Histogram = (TH1D*) GamSimu_NeutronMV_Histogram->Clone();
    TH1D* SumSimu_VetoedMuons_MV_Histogram = (TH1D*) GamSimu_VetoedMuons_MV_Histogram->Clone();
    TH1D* SumSimu_VetoedMuons_BS_Histogram = (TH1D*) GamSimu_VetoedMuons_BS_Histogram->Clone();
    
    SumSimu_Neutron_Histogram->Add(SumSimu_Neutron_Histogram, MuonSimu_Neutron_Histogram, 1, 1);
    SumSimu_Muon_Histogram->Add(SumSimu_Muon_Histogram, MuonSimu_Muon_Histogram, 1, 1);
    SumSimu_NeutronMV_Histogram->Add(SumSimu_NeutronMV_Histogram, MuonSimu_NeutronMV_Histogram, 1, 1);
    SumSimu_VetoedMuons_BS_Histogram->Add(SumSimu_VetoedMuons_BS_Histogram, MuonSimu_VetoedMuons_BS_Histogram, 1, 1);
    SumSimu_VetoedMuons_MV_Histogram->Add(SumSimu_VetoedMuons_MV_Histogram, MuonSimu_VetoedMuons_MV_Histogram, 1, 1);
    
    SumSimu_Neutron_Histogram->Add(SumSimu_Neutron_Histogram, NeutronsSimu_Neutron_Histogram, 1, 1);
    SumSimu_Muon_Histogram->Add(SumSimu_Muon_Histogram, NeutronsSimu_Muon_Histogram, 1, 1);
    SumSimu_NeutronMV_Histogram->Add(SumSimu_NeutronMV_Histogram, NeutronsSimu_NeutronMV_Histogram, 1, 1);
    SumSimu_VetoedMuons_BS_Histogram->Add(SumSimu_VetoedMuons_BS_Histogram, NeutronsSimu_VetoedMuons_BS_Histogram, 1, 1);
    SumSimu_VetoedMuons_MV_Histogram->Add(SumSimu_VetoedMuons_MV_Histogram, NeutronsSimu_VetoedMuons_MV_Histogram, 1, 1);
    
    SumSimu_Neutron_Histogram->SetLineColor(kRed);
    SumSimu_Muon_Histogram->SetLineColor(kRed);
    SumSimu_NeutronMV_Histogram->SetLineColor(kRed);
    SumSimu_VetoedMuons_BS_Histogram->SetLineColor(kRed);
    SumSimu_VetoedMuons_MV_Histogram->SetLineColor(kRed);
    
    
    TCanvas* cMuonPan = new TCanvas("MuonPanel", "MuonPanel", 1200, 1200);
    cMuonPan->Divide(1,2);
    cMuonPan->cd(1);
    Muon_Histogram->Draw();
    MuonSimu_Muon_Histogram->Draw("SAME");
    NeutronsSimu_Muon_Histogram->Draw("SAME");
    GamSimu_Muon_Histogram->Draw("SAME");
    SumSimu_Muon_Histogram->Draw("SAME");
    cMuonPan->cd(2);
    VetoedMuons_MV_Histogram->Draw();
    MuonSimu_VetoedMuons_MV_Histogram->Draw("SAME");
    NeutronsSimu_VetoedMuons_MV_Histogram->Draw("SAME");
    GamSimu_VetoedMuons_MV_Histogram->Draw("SAME");
    SumSimu_VetoedMuons_MV_Histogram->Draw("SAME");
    
    TCanvas* cBS = new TCanvas("BS", "BS", 1200, 1200);
    cBS->Divide(2,2);
    cBS->cd(1);
    Neutron_Histogram->Draw();
    MuonSimu_Neutron_Histogram->Draw("SAME");
    NeutronsSimu_Neutron_Histogram->Draw("SAME");
    GamSimu_Neutron_Histogram->Draw("SAME");
    SumSimu_Neutron_Histogram->Draw("SAME");
    cBS->cd(2);
    NeutronMV_Histogram->Draw();
    MuonSimu_NeutronMV_Histogram->Draw("SAME");
    NeutronsSimu_NeutronMV_Histogram->Draw("SAME");
    GamSimu_NeutronMV_Histogram->Draw("SAME");
    SumSimu_NeutronMV_Histogram->Draw("SAME");
    cBS->cd(3);
    VetoedMuons_BS_Histogram->Draw();
    MuonSimu_VetoedMuons_BS_Histogram->Draw("SAME");
    NeutronsSimu_VetoedMuons_BS_Histogram->Draw("SAME");
    GamSimu_VetoedMuons_BS_Histogram->Draw("SAME");
    SumSimu_VetoedMuons_BS_Histogram->Draw("SAME");
    cBS->cd(4);
    Neutron_Histogram->Draw();
    NeutronMV_Histogram->SetLineColor(kCyan);
//    Simu_Neutron_Histogram->Draw("SAME");
//    Simu_NeutronMV_Histogram->Draw("SAME");
//    Simu_NeutronMV_Histogram->SetLineColor(kOrange);
}
