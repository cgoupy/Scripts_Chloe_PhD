// Macro to normalize in rate the simulation obtained with EnerGe

TH1D* X_rescale(TH1D* h_spectrum, TString newName, Double_t ScaleFactor){ //Change the X-scale of the histograms by a scaling factor:

    int NBins = h_spectrum->GetNbinsX();
    Double_t maxE = h_spectrum->GetBinLowEdge(NBins-1)+h_spectrum->GetBinWidth(NBins-1);
    Double_t maxE_scaled = maxE*ScaleFactor;
    TH1D* h_spectrum_scaled = new TH1D(newName, newName,NBins, 0, maxE_scaled);
    
    Double_t Countsbin_i, Energybin_i;
    for (int i=0; i<NBins; i++){
        Countsbin_i = h_spectrum->GetBinContent(i+1);
        Energybin_i = h_spectrum->GetBinLowEdge(i+1);
        for (int j=0; j<Countsbin_i; j++){
            h_spectrum_scaled->Fill(Energybin_i*ScaleFactor); //this method instead of SetBinContent to keep the stat of the histogram ;)
        }
    }
    return h_spectrum_scaled;
}

void Normalize_Simu(){
    TString Filename_in = "../GammaGenerator/AmbientGammas/Analyzer_PE_PortableHPGe_AmbientGammas_fromDeconv.000.root";
    
    TFile* filein = new TFile(Filename_in, "READ");
    TH1D* h_spectrum = (TH1D*) filein->Get("EdepAlgoDir/Vol_HPGe/Edep/Edep_Vol_HPGe_0-3MeV"); //to be changed if needed
    
    TH1D* h_primaries = (TH1D*) filein->Get("PrimariesAlgoDir/Energy/hEinj"); //get the primary spectrum
    
    //Change the X-scale of the histograms (from MeV to keV):
    h_spectrum->Rebin(10);
    
    TH1D* h_primaries_keV = X_rescale(h_primaries, "primary_gammas", 1e3);
    TH1D* h_spectrum_keV = X_rescale(h_spectrum, "Edep_in_HPGe", 1e3);
    
    h_spectrum_keV->GetXaxis()->SetTitle("E_{dep} [keV]");
    h_primaries_keV->GetXaxis()->SetTitle("Energy [keV]");
    
    //Get the number of generated primaries:
    int N_entries = h_primaries_keV->GetEntries();

    //Calculate the "equivalent time"
    Double_t flux = 4.3; //ev/cm2/s reference value for gammas
    Double_t Plane_size = 10*10; //cm2 size of the tangent plane
    Double_t Eq_time_s = N_entries/(flux*Plane_size); //in seconds
    Double_t Eq_time_d = Eq_time_s/(24*3600.); //in days
    
    //Mass of the detector
    Double_t Volume =(5.5*TMath::Pi()*pow(5.55/2,2)-4.4*TMath::Pi()*pow(1.15/2,2));
    Double_t Mass = Volume*5.323/1000;
    
    TH1D* h_spectrum_keV_rate = (TH1D*) h_spectrum_keV->Clone();
    
    //Normalize the histograms
    h_spectrum_keV->Scale(1./(Mass*Eq_time_d*h_spectrum_keV->GetBinWidth(1))); //in ev/d/kg/keV
    h_spectrum_keV_rate->Scale(1./Eq_time_s); //in ev/s
    h_primaries_keV->Scale(1./(Plane_size*Eq_time_d*h_primaries_keV->GetBinWidth(1))); //in ev/d/keV
    
    h_spectrum_keV->GetYaxis()->SetTitle("Differential Rate [ev kg^{-1} d^{-1} keV^{-1}]");
    h_spectrum_keV_rate->GetYaxis()->SetTitle("Counting Rate [ev s^{-1}]");
    h_primaries_keV->GetYaxis()->SetTitle("Differential Flux [ev cm^{-2} d^{-1} keV^{-1}]");
    
    TString Filename_out = Filename_in;
    TPMERegexp("Analyzer").Substitute(Filename_out, Form("Normalized_Flux_%0.1f_Analyzer", flux)); //add Normalized to the name of the file, to create a copy with normalized histograms and not overwrite on the simulation file.

    TFile* fileout = new TFile(Filename_out, "RECREATE");
    h_primaries_keV->Write("Primary_gamma_spectrum");
    h_spectrum_keV->Write("Edep_gamma_spectrum_dru");
    h_spectrum_keV_rate->Write("Edep_gamma_spectrum_rate");
    
    fileout->Close();
    filein->Close();
    cout<<"Done! Histo written in"<< Filename_out<<endl;
}
