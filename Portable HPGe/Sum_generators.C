// Small macro to sum the gamma spectra obtained with the fill Gun

void Sum_generator(){
    Double_t Room_surf = 1125200; //cm-2
    //Tables to calculate the yield
    Double_t Nucleide_conc[3] = {60, 35, 392}; // Ra, Th, K (in Bq/kg)
    Double_t Mass_sources[4] = {5.97E5, 5.34E5, 2.48E7, 9.8E4}; //Basement, VNS walls, Ground, PartWalls
    const int NFiles=12;
    TString files[NFiles] = {
        "U238/Analyzer_PS_U238_FillGun_Basementwalls.root",
        "U238/Analyzer_PS_U238_FillGun_VNSwalls.root",
        "U238/Analyzer_PS_U238_FillGun_Ground.root",
        "U238/Analyzer_PS_U238_FillGun_PartitionWalls.root",
        "Th232/Analyzer_PS_Th232_FillGun_Basementwalls.root",
        "Th232/Analyzer_PS_Th232_FillGun_VNSwalls.root",
        "Th232/Analyzer_PS_Th232_FillGun_Ground.root",
        "Th232/Analyzer_PS_Th232_FillGun_PartitionWalls.root",
        "K40/Analyzer_PS_K40_FillGun_BasementWalls.root",
        "K40/Analyzer_PS_K40_FillGun_VNSwalls.root",
        "K40/Analyzer_PS_K40_FillGun_Ground.root",
        "K40/Analyzer_PS_K40_FillGun_PartitionWalls.root"};
    
    TFile* fileout = new TFile("../GammaGenerator/Gammas_in_VNS_U238_Th232_K40.root", "RECREATE");
    fileout->mkdir("U238");
    fileout->mkdir("Th232");
    fileout->mkdir("K40");
    
    TH1D* h_Sums[7];
    
    for (int i=0; i<NFiles; i++){
        TFile* filein = new TFile(Form("../GammaGenerator/%s",files[i].Data()), "READ");
        cout<<files[i]<<" open"<<endl;
        
        TH1D* histoPrim = (TH1D*) filein->Get("PrimariesAlgoDir/Energy/hEinj");
        int Nprim = histoPrim->GetEntries();
        
        TH1D* histo_gammas = (TH1D*) filein->Get("SecondariesAlgoDir/Gammas/nEinc_gamma");
        histo_gammas->Rebin(1);
        
        Double_t Scaling_fact = (Nucleide_conc[i/4]*Mass_sources[i-((i/4)*4)])/Nprim;
        histo_gammas->Scale(Scaling_fact/(histo_gammas->GetBinWidth(1)*Room_surf));
        
        histo_gammas->GetYaxis()->SetTitle("Differential flux [MeV^{-1} cm^{-2} s^{-1}]");
        
        TString Name =files[i];
        TPMERegexp("\\.root").Substitute(Name, "");
        TPMERegexp("^\\w*\\/\\w*\\_").Substitute(Name, "");
        
        cout<<"ok"<<endl;
        switch (i/4){
            case 0:
                fileout->cd("U238");
                if (i==0) h_Sums[i/4] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i/4]->Add(h_Sums[i/4], histo_gammas, 1., 1.);
                break;
            
            case 1:
                fileout->cd("Th232");
                if (i==4) h_Sums[i/4] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i/4]->Add(h_Sums[i/4], histo_gammas, 1., 1.);
                break;
            
            case 2:
                fileout->cd("K40");
                if (i==8) h_Sums[i/4] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i/4]->Add(h_Sums[i/4], histo_gammas, 1., 1.);
                break;
        }
        
        switch (i-((i/4)*4)){
            case 0:
                if (i==0) h_Sums[i-((i/4)*4)+3] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i-((i/4)*4)+3]->Add(h_Sums[i-((i/4)*4)+3], histo_gammas, 1., 1.);
                break;
            
            case 1:
                if (i==1) h_Sums[i-((i/4)*4)+3] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i-((i/4)*4)+3]->Add(h_Sums[i-((i/4)*4)+3], histo_gammas, 1., 1.);
                break;
            
            case 2:
                if (i==2) h_Sums[i-((i/4)*4)+3] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i-((i/4)*4)+3]->Add(h_Sums[i-((i/4)*4)+3], histo_gammas, 1., 1.);
                break;
                
            case 3:
                if (i==3) h_Sums[i-((i/4)*4)+3] = (TH1D*) histo_gammas->Clone();
                else h_Sums[i-((i/4)*4)+3]->Add(h_Sums[i-((i/4)*4)+3], histo_gammas, 1., 1.);
                break;
        }
        
        histo_gammas->Write(Name);
        fileout->cd();
        cout<<Name<<" written"<<endl;
    }
    
    TH1D* h_Sums_tot = (TH1D*) h_Sums[0]->Clone();
    h_Sums_tot->Add(h_Sums_tot, h_Sums[1], 1., 1.);
    h_Sums_tot->Add(h_Sums_tot, h_Sums[2], 1., 1.);
    
    fileout->cd("U238");
    h_Sums[0]->Write("Sum_U238");
    
    fileout->cd("Th232");
    h_Sums[1]->Write("Sum_Th232");
    
    fileout->cd("K40");
    h_Sums[2]->Write("Sum_K40");
    
    fileout->cd();
    h_Sums[3]->Write("Sum_BasementWalls");
    h_Sums[4]->Write("Sum_VNSWalls");
    h_Sums[5]->Write("Sum_Ground");
    h_Sums[6]->Write("Sum_PartitionWalls");
    
    h_Sums_tot->Write("Sum_Tot");
    
}


void PrintDat_Gamma(){
    TFile* filein = new TFile("../GammaGenerator/Gammas_in_VNS_U238_Th232_K40.root");
    TH1D* histos[3];
    histos[0] = (TH1D*) filein->Get("U238/Sum_U238");
    histos[1] = (TH1D*) filein->Get("Th232/Sum_Th232");
    histos[2] = (TH1D*) filein->Get("K40/Sum_K40");
    
    TString txtnames[3] = {"AmbientU238Spectrum.dat", "AmbientTh232Spectrum.dat", "AmbientK40Spectrum.dat" };
    ofstream myfile;
    
    for (int k=0; k<3; k++){
        myfile.open(txtnames[k]);
    
        int Nbins = histos[k]->GetNbinsX();
        //histo->Rebin(Nbins/600);
        
        Nbins = histos[k]->GetNbinsX();
        
        cout<<txtnames[k]<<endl;
        cout<<"=============================="<<endl;
        myfile <<  0<<"       "<<0<<endl;
        double E_i;
        for (int i=1; i<=Nbins; i++){
            E_i=histos[k]->GetBinCenter(i);
            cout <<  E_i  <<"       "<<histos[k]->GetBinContent(i)<<endl;
            myfile <<  E_i<<"       "<<histos[k]->GetBinContent(i)<<endl;
        }
        cout<<endl;
        
        myfile.close();
    }
}

