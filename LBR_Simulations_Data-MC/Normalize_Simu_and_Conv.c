// Macro to normalize in rate the simulation obtained with EnerGe

TH1D* X_rescale(TH1D* h_spectrum, TString newName, Double_t ScaleFactor){ //Change the X-scale of the histograms by a scaling factor:

    int NBins = h_spectrum->GetNbinsX();
    Double_t maxE = h_spectrum->GetBinLowEdge(NBins)+h_spectrum->GetBinWidth(NBins);
    Double_t maxE_scaled = maxE*ScaleFactor;
    TH1D* h_spectrum_scaled = new TH1D(newName, newName, NBins, 0, maxE_scaled); // keep the same binning
    
    Double_t Countsbin_i, Energybin_i, Errorsbin_i;
    for (int i=0; i<NBins; i++){
        Countsbin_i = h_spectrum->GetBinContent(i+1);
        Errorsbin_i = h_spectrum->GetBinError(i+1);
        Energybin_i = h_spectrum->GetBinLowEdge(i+1);
        h_spectrum_scaled->SetBinContent(i+1, Countsbin_i);
        h_spectrum_scaled->SetBinError(i+1, Errorsbin_i);
    }
    
    h_spectrum_scaled->Rebin(5);
    h_spectrum_scaled->GetYaxis()->SetTitle(h_spectrum->GetYaxis()->GetTitle());
    
    return h_spectrum_scaled;
}

TH1D* Apply_resol_Ge(TH1D* h_spectrum, double A, double B, double C){

    //--------- Put the TH1D in vectors
    int N_Sim = h_spectrum->GetNbinsX();
    TVectorD Simu_vect(N_Sim); //spectrum in a vector
    TVectorD Sim_bins(N_Sim); //Bin centers (energies)
    TMatrixD V_Simu(N_Sim, N_Sim); //variance matrix
    
    //Sim_bins are the bins centers of the simulation histogram
    //Simu_vect is the vector of the simulation values
    //V_Simu is the variance matrix of the simulation
    
    Double_t E_k, err_k, cont_k;
    for (int k=1; k<=N_Sim; k++){
        cont_k = h_spectrum->GetBinContent(k);
        E_k = h_spectrum->GetBinCenter(k);
        err_k = h_spectrum->GetBinError(k);
        Simu_vect(k-1) = cont_k;
        Sim_bins(k-1) = E_k;
        V_Simu(k-1, k-1) = err_k*err_k;
    }
    
    // --------- Convolution
    TMatrixD Trans(N_Sim, N_Sim); //transfer matrix of vector dimension
    TMatrixD V_Conv(N_Sim, N_Sim); //Variance matrix of the convolution
    TVectorD SimuConv(N_Sim); //Convolution result
    
    Double_t E_i, E_j, sigma_j;
    for (int i =0; i<N_Sim ; i++){ //construct the transfert matrix with gaussians of the resolution
        E_i = Sim_bins(i);
        for (int j =0; j< N_Sim; j++){
            E_j=Sim_bins(j);
            sigma_j = sqrt(A*E_j+B+C*E_j*E_j)/100; //devided by 100 because the fit was done in %
            Trans(i,j)=TMath::Gaus(E_j, E_i, sigma_j, true); //TransferMatrix: the bin (i;j) contains the value at E_i of the gaussian centered in E_j with resolution sigma_j. The gaussian is normed (true) to conserve the normalization of the simulation
        }
    }
        
    SimuConv = Trans*Simu_vect;
    
    TMatrixD Trans_T(TMatrixD::kTransposed, Trans); //transpose matrix of Trans
    V_Conv = Trans*V_Simu*Trans_T;
    
    //------------- Put the result in a TH1D*
    TH1D* h_spectrum_Conv = (TH1D*) h_spectrum->Clone();//clone because the X axis is the same
    h_spectrum_Conv->SetTitle(Form("%s_conv", h_spectrum_Conv->GetTitle()));
    for (int p=1; p<=N_Sim; p++){
        h_spectrum_Conv->SetBinContent(p, SimuConv(p-1));
        h_spectrum_Conv->SetBinError(p, sqrt(V_Conv(p-1, p-1)));
    }
    
    return h_spectrum_Conv;
}


TH1D* Apply_const_resol(TH1D* h_spectrum, double A){

    //--------- Put the TH1D in vectors
    int N_Sim = h_spectrum->GetNbinsX();
    TVectorD Simu_vect(N_Sim); //spectrum in a vector
    TVectorD Sim_bins(N_Sim); //Bin centers (energies)
    TMatrixD V_Simu(N_Sim, N_Sim); //variance matrix
    
    //Sim_bins are the bins centers of the simulation histogram
    //Simu_vect is the vector of the simulation values
    //V_Simu is the variance matrix of the simulation
    
    Double_t E_k, err_k, cont_k;
    for (int k=1; k<=N_Sim; k++){
        cont_k = h_spectrum->GetBinContent(k);
        E_k = h_spectrum->GetBinCenter(k);
        err_k = h_spectrum->GetBinError(k);
        Simu_vect(k-1) = cont_k;
        Sim_bins(k-1) = E_k;
        V_Simu(k-1, k-1) = err_k*err_k;
    }
    
    // --------- Convolution
    TMatrixD Trans(N_Sim, N_Sim); //transfer matrix of vector dimension
    TMatrixD V_Conv(N_Sim, N_Sim); //Variance matrix of the convolution
    TVectorD SimuConv(N_Sim); //Convolution result
    
    Double_t E_i, E_j, sigma_j;
    for (int i =0; i<N_Sim ; i++){ //construct the transfert matrix with gaussians of the resolution
        E_i = Sim_bins(i);
        for (int j =0; j< N_Sim; j++){
            E_j=Sim_bins(j);
            sigma_j = A/100*E_j; //devided by 100 because the fit was done in %
            Trans(i,j)=TMath::Gaus(E_j, E_i, sigma_j, true); //TransferMatrix: the bin (i;j) contains the value at E_i of the gaussian centered in E_j with resolution sigma_j.
        }
    }
        
    SimuConv = Trans*Simu_vect;
    
    TMatrixD Trans_T(TMatrixD::kTransposed, Trans); //transpose matrix of Trans
    V_Conv = Trans*V_Simu*Trans_T;
    
    //------------- Put the result in a TH1D*
    TH1D* h_spectrum_Conv = (TH1D*) h_spectrum->Clone();//clone because the X axis is the same
    h_spectrum_Conv->SetTitle(Form("%s_conv", h_spectrum_Conv->GetTitle()));
    for (int p=1; p<=N_Sim; p++){
        h_spectrum_Conv->SetBinContent(p, SimuConv(p-1));
        h_spectrum_Conv->SetBinError(p, sqrt(V_Conv(p-1, p-1)));
    }
    
    return h_spectrum_Conv;
}

TH1D* Apply_lin_resol(TH1D* h_spectrum, double A, double B){

    //--------- Put the TH1D in vectors
    int N_Sim = h_spectrum->GetNbinsX();
    TVectorD Simu_vect(N_Sim); //spectrum in a vector
    TVectorD Sim_bins(N_Sim); //Bin centers (energies)
    TMatrixD V_Simu(N_Sim, N_Sim); //variance matrix
    
    //Sim_bins are the bins centers of the simulation histogram
    //Simu_vect is the vector of the simulation values
    //V_Simu is the variance matrix of the simulation
    
    Double_t E_k, err_k, cont_k;
    for (int k=1; k<=N_Sim; k++){
        cont_k = h_spectrum->GetBinContent(k);
        E_k = h_spectrum->GetBinCenter(k);
        err_k = h_spectrum->GetBinError(k);
        Simu_vect(k-1) = cont_k;
        Sim_bins(k-1) = E_k;
        V_Simu(k-1, k-1) = err_k*err_k;
    }
    
    // --------- Convolution
    TMatrixD Trans(N_Sim, N_Sim); //transfer matrix of vector dimension
    TMatrixD V_Conv(N_Sim, N_Sim); //Variance matrix of the convolution
    TVectorD SimuConv(N_Sim); //Convolution result
    
    Double_t E_i, E_j, sigma_j;
    for (int i =0; i<N_Sim ; i++){ //construct the transfert matrix with gaussians of the resolution
        E_i = Sim_bins(i);
        for (int j =0; j< N_Sim; j++){
            E_j=Sim_bins(j);
            sigma_j = A+B*E_j;
            Trans(i,j)=TMath::Gaus(E_j, E_i, sigma_j, true); //TransferMatrix: the bin (i;j) contains the value at E_i of the gaussian centered in E_j with resolution sigma_j.
        }
    }
        
    SimuConv = Trans*Simu_vect;
    
    TMatrixD Trans_T(TMatrixD::kTransposed, Trans); //transpose matrix of Trans
    V_Conv = Trans*V_Simu*Trans_T;
    
    //------------- Put the result in a TH1D*
    TH1D* h_spectrum_Conv = (TH1D*) h_spectrum->Clone();//clone because the X axis is the same
    h_spectrum_Conv->SetTitle(Form("%s_conv", h_spectrum_Conv->GetTitle()));
    for (int p=1; p<=N_Sim; p++){
        h_spectrum_Conv->SetBinContent(p, SimuConv(p-1));
        h_spectrum_Conv->SetBinError(p, sqrt(V_Conv(p-1, p-1)));
    }
    
    return h_spectrum_Conv;
}

void Conv_1_file(TString Dir_Spectrum, TString Name_spectrum, TFile* filein, TFile* fileout, Double_t Emin, Double_t Emax){
    TH1D* h_spectrum = (TH1D*) filein->Get(Dir_Spectrum+Name_spectrum); //to be changed if needed
    
    //Change the X-scale of the histograms (from MeV to keV):
    h_spectrum->Rebin(1);
    
    TH1D* h_spectrum_keV = X_rescale(h_spectrum, Form("%s_keV", h_spectrum->GetTitle()), 1.e3);
    
    h_spectrum_keV->GetXaxis()->SetTitle("Energy [keV]");
    
    //Apply the resolution
//    TH1D* h_spectrum_keV_conv = Apply_const_resol(h_spectrum_keV, 2.3);
//    TH1D* h_spectrum_keV_conv = Apply_resol_Ge(h_spectrum_keV, 986.387, -165813, 0.0757142);
//    TH1D* h_spectrum_keV_conv = Apply_resol_Ge(h_spectrum_keV, -15834.7, 1.90052e+07, 3.7);
//    TH1D* h_spectrum_keV_conv = Apply_lin_resol(h_spectrum_keV, 3.53434, 0.00530013);
      TH1D* h_spectrum_keV_conv = Apply_lin_resol(h_spectrum_keV, 1., 0.04);
//    h_spectrum_keV_conv->GetYaxis()->SetRangeUser(h_spectrum_keV->GetYaxis()->GetXmin(), h_spectrum_keV->GetYaxis()->GetXmax());
    
    //Scale
    Double_t Integral_tot = h_spectrum_keV->Integral(h_spectrum_keV->FindBin(Emin*1e3), h_spectrum_keV->FindBin(Emax*1e3));
    Double_t Integral_conv = h_spectrum_keV_conv->Integral(h_spectrum_keV_conv->FindBin(Emin*1e3), h_spectrum_keV_conv->FindBin(Emax*1e3));
    h_spectrum_keV_conv->Scale(Integral_tot*1./Integral_conv);
    
    h_spectrum_keV->Write("Edep_"+Name_spectrum);
    h_spectrum_keV_conv->Write("Evis_"+Name_spectrum);
    fileout->Close();
}

void Conv_Simu(TString Filename_in = "../Simulation/AmbientGammas/U238/Analyzer_PE_PortableHPGe_ColdFinger_LeadShielding_Ambient_U238_2.000-132.root"){
    
    TFile* filein = new TFile(Filename_in, "READ");
    
    TString Filename_out = Filename_in;
    TPMERegexp("Histos").Substitute(Filename_out, "ConvGe_Histos"); //add Normalized to the name of the file, to create a copy with normalized histograms and not overwrite on the simulation file.
    TFile* fileout = new TFile(Filename_out, "RECREATE");
    
    //================================ Energy range 0-5 MeV
    std::cout<<"0-5MeV"<<std::endl;
    //------------------------- Singles
    std::cout<<"Singles..."<<std::endl;
    Conv_1_file("", "COV_Histogram_0_5MeV", filein, fileout, 0.5, 4);
    
    //------------------------- Veto
    std::cout<<"Veto..."<<std::endl;
    
    fileout = new TFile(Filename_out, "UPDATE");
    fileout->mkdir("MVCut");
    fileout->cd("MVCut");
    Conv_1_file("MVCut/", "MVCut_COV_Histogram_0_5MeV", filein, fileout, 0.5, 4);
    
    //------------------------- Coinc
    std::cout<<"Coinc..."<<std::endl;

    fileout = new TFile(Filename_out, "UPDATE");
    fileout->mkdir("MVCoinc");
    fileout->cd("MVCoinc");
    Conv_1_file("MVcoinc/", "MVCoinc_COV_Histogram_0_5MeV", filein, fileout, 0.5, 4);
    
    //================================ Energy range 0-100 MeV
    std::cout<<"0-100MeV"<<std::endl;
    //------------------------- Singles
    std::cout<<"Singles..."<<std::endl;
    fileout = new TFile(Filename_out, "UPDATE");
    Conv_1_file("", "COV_Histogram_0_100MeV", filein, fileout, 1, 90);
    
    //------------------------- Veto
    std::cout<<"Veto..."<<std::endl;
    
    fileout = new TFile(Filename_out, "UPDATE");
    fileout->cd("MVCut");
    Conv_1_file("MVCut/", "MVCut_COV_Histogram_0_100MeV", filein, fileout, 1, 90);
    
    //------------------------- Coinc
    std::cout<<"Coinc..."<<std::endl;

    fileout = new TFile(Filename_out, "UPDATE");
    fileout->cd("MVCoinc");
    Conv_1_file("MVcoinc/", "MVCoinc_COV_Histogram_0_100MeV", filein, fileout, 1, 90);
    
    filein->Close();
    cout<<"Done! Histo written in"<< Filename_out<<endl;
}

void Normalize_Simu(TString Filename_in = "../Simulation/AmbientGammas/AmbientGammas/Conv_Analyzer_PE_PortableHPGe_AmbientGammas_2.000-001.root"){
    
    TFile* filein = new TFile(Filename_in, "READ");
    
    TH1D* h_primaries = (TH1D*) filein->Get("Primary_gamma_spectrum");
    TH1D* h_spectrum_keV = (TH1D*) filein->Get("Edep_gamma_spectrum");
    TH1D* h_spectrum_keV_conv= (TH1D*) filein->Get("Evis_gamma_spectrum");
    
    //Get the number of generated primaries:
    int N_entries = h_primaries->GetEntries();
    TH1D* h_primaries_keV = X_rescale(h_primaries, Form("%s_keV", h_primaries->GetTitle()), 1e3);
    h_primaries_keV->GetXaxis()->SetTitle("Energy [keV]");
    
    cout<<N_entries<<endl;
    //Calculate the "equivalent time"
    Double_t flux = 8.9347797e-05; //ev/cm2/s reference value for Radon (60 bd/m3)
    Double_t Plane_size = 120*120;//160*160; //cm2 size of the tangent plane
    Double_t Eq_time_s = N_entries/(flux*Plane_size); //in seconds
//    cout<<Eq_time_s<<endl;
    Double_t Eq_time_d = Eq_time_s/(24*3600.); //in days
    
    //Mass of the detector
    Double_t Mass = 1; //1 kg
    
    TH1D* h_spectrum_keV_rate = (TH1D*) h_spectrum_keV_conv->Clone();
    
    TString Filename_out = Filename_in;
    TPMERegexp("Conv_Analyzer").Substitute(Filename_out, Form("Normalized_Flux_%0.3f_Conv_Analyzer", flux)); //add Normalized to the name of the file, to create a copy with normalized histograms and not overwrite on the simulation file.
    
    TFile* fileout = new TFile(Filename_out, "RECREATE");
    //Copy before normalizing
    h_primaries->Write("Primary_gamma_spectrum");
    h_spectrum_keV->Write("Edep_gamma_spectrum");
    h_spectrum_keV_conv->Write("Evis_gamma_spectrum");
    
    //Normalize the histograms
    h_spectrum_keV_conv->Scale(1./(Mass*Eq_time_d*h_spectrum_keV_conv->GetBinWidth(1))); //in ev/d/kg/keV
    h_spectrum_keV_rate->Scale(1./Eq_time_s); //in ev/s
    h_primaries_keV->Scale(1./(Eq_time_d*h_primaries_keV->GetBinWidth(1))); //in ev/d/keV
    
    h_spectrum_keV_conv->GetYaxis()->SetTitle("Differential Rate [ev kg^{-1} d^{-1} keV^{-1}]");
    h_spectrum_keV_rate->GetYaxis()->SetTitle("Counting Rate [ev s^{-1}]");
    h_primaries_keV->GetYaxis()->SetTitle("Differential Flux [ev cm^{-2} d^{-1} keV^{-1}]");
    
    h_primaries_keV->Write("Normalized_Primary_gamma_spectrum");
    h_spectrum_keV->Write("Normalized_Edep_gamma_spectrum");
    h_spectrum_keV_conv->Write("Normalized_Evis_gamma_spectrum_dru");
    h_spectrum_keV_rate->Write("Normalized_Evis_gamma_spectrum_rate");
    
    fileout->Close();
    filein->Close();
    cout<<"Done! Histo written in"<< Filename_out<<endl;
}
