//Fit a simulation on a data histogram
Double_t PLANE_SURF=160*160;
Double_t Volume =(5.5*TMath::Pi()*pow(5.55/2,2)-4.4*TMath::Pi()*pow(1.15/2,2));
Double_t MASS_Ge = Volume*5.323/1000;
int REBIN=4;

//Fuction to cut an histogram between loweredge and upperedge
TH1D* Cut_Histo(TH1D* h_qdc, double QDC_loweredge, double QDC_upperedge) {
    
    
    TString Name = h_qdc->GetName();
    Name+="_cut";
    TString Title = h_qdc->GetTitle();
    Title+="_cut";
    
    TH1D* h_qdc_cut_old = (TH1D*) gROOT->FindObject(Name);
    delete h_qdc_cut_old;
    
    TH1D* h_qdc_cut = new TH1D(Name,
                               Title,
                               h_qdc->FindBin(QDC_upperedge) - (h_qdc->FindBin(QDC_loweredge)),
                               h_qdc->GetBinLowEdge(h_qdc->FindBin(QDC_loweredge)),
                               h_qdc->GetBinLowEdge(h_qdc->FindBin(QDC_upperedge)+1)
                               );
    
    for(int i_bin = 1; i_bin <=h_qdc_cut->GetNbinsX()+1 ; i_bin++) {
        int i_bin_h = h_qdc->FindBin(QDC_loweredge) + i_bin;
        h_qdc_cut->SetBinContent(i_bin, h_qdc->GetBinContent(i_bin_h));
        h_qdc_cut->SetBinError(i_bin, h_qdc->GetBinError(i_bin_h));
        
    }
    
    h_qdc_cut->GetXaxis()->SetTitle(h_qdc->GetXaxis()->GetTitle());
    h_qdc_cut->GetYaxis()->SetTitle(h_qdc->GetYaxis()->GetTitle());
    
    return h_qdc_cut;
    
}

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

void Draw_Vector_TH1D(TVectorT<Double_t> Vector, TString Opt){
    TH1D* histo = new TH1D(Vector);
    histo->Draw(Opt);
}

// --------------- Fit function
ROOT::Fit::FitResult Fit_Data_Simu(std::vector<Double_t> NumEvents, TVectorT<Double_t> E_bins, TVectorT<Double_t> Data_vect, std::vector<TVectorT<Double_t>> Simu_vect, std::vector<TMatrixT<Double_t>> V_Simu, TMatrixT<Double_t> V_Data, std::vector<TString> Names, std::vector<Double_t> ini_flux) {

    //E_bins are the bins centers of the data and simulation histogram
    //Data_vect is the vector of the data values
    //Simu_vect is the vector of the simulation values
    //V_Data and V_Simu are the variance matrix of the data and the simulation
    
    int N_E = E_bins.GetNrows();
    double deltaE = E_bins[1]-E_bins[0];
    //TVectorT<Double_t> Delta_Data_SimuConv(N_NPE); //difference between Data and alpha*Simu
    
    //Projection on a base where diagonal elements V_Data =/= 0
    int N_zeros = 0;
    
    for (int i =0; i<N_E; i++){
        N_zeros+=Data_vect(i)==0; //count the number of zero in Data_vect
    }
    
    int N_Proj = N_E-N_zeros;
    cout<<N_Proj<<endl;
    TMatrixT<Double_t> Proj(N_Proj, N_E);
    TMatrixT<Double_t> Proj_T(N_E, N_Proj);
    
    int i_Proj = 0;
    for (int j=0; j<N_E; j++){
        if (Data_vect(j)!=0){
            Proj(i_Proj, j)=1;
            Proj_T(j, i_Proj)=1;
            i_Proj+=1;
        }
    }

    
    //---- projection data
    cout<<"Project Data..."<<endl;
    TVectorT<Double_t> Proj_Data_vect = Proj*Data_vect;
    
    //---- projection V_Data
    TMatrixT<Double_t> Proj_V_Data = Proj*V_Data*Proj_T;
    cout<<"size Proj_V_Data " <<Proj_V_Data.GetNrows()<<"x"<<Proj_V_Data.GetNcols()<<endl;
    
    //---- projection simulation spectra and V_simu
    int N_simu = Simu_vect.size();
    std::vector<TVectorT<Double_t>> Proj_Simu_vect;
    std::vector<TMatrixT<Double_t>> Proj_V_Simu;
    TVectorT<Double_t> FormatVector(N_Proj);
    TMatrixT<Double_t> FormatMatrix(N_Proj, N_Proj);
    
    for (int i=0; i<N_simu; i++){
        Proj_Simu_vect.push_back(FormatVector);
        Proj_V_Simu.push_back(FormatMatrix);
    }
    
    for (int i=0; i<N_simu; i++){
        cout<<"Project simu n° " <<i<<"..."<<endl;
        Proj_Simu_vect[i] = Proj*Simu_vect[i];
        Proj_V_Simu[i] = Proj*V_Simu[i]*Proj_T;
    }
    
    //---- projection E_bins
    TVectorT<Double_t> Proj_E_bins = Proj*E_bins;
    
    //--- Transfer matrix
//    TMatrixT<Double_t> V_Delta_i(N_Proj, N_Proj);
    TMatrixT<Double_t> V_Delta(N_Proj, N_Proj);
    TVectorT<Double_t> Delta_Data_Simu(N_Proj);

    //---- ease calculation of inverted matrix: diagonal matrices since only statistical errors from data considered
//    for (int j=0; j<N_Proj; j++){
//        V_Delta_i(j,j) = 1./Proj_V_Data(j,j);
//    }
    
    //cf. https://root.cern.ch/doc/master/fitCircle_8C.html
    auto chi2Function_Fit = [&](const Double_t *par) {
        //minimisation function computing the sum of squares of chi2
        int N_par=par[0]; //number of simulation to fit
        Double_t Norm[N_par];
        
        //Norm from flux
        for (int i=0; i<N_par; i++){
            Norm[i]= (par[i+1]*PLANE_SURF*24.*3600.)/(NumEvents[i]*MASS_Ge*deltaE);
        }
        //Diff Data - normalized simulations
        Delta_Data_Simu = Proj_Data_vect;
        V_Delta = Proj_V_Data;
        
        for (int i=0; i<N_par; i++){
            Delta_Data_Simu -= Norm[i]*Proj_Simu_vect[i];
            V_Delta += Norm[i]*Norm[i]*Proj_V_Simu[i];
        }
        
        //Variance matrix Data - normalized simulations
         TMatrixT<Double_t> V_Delta_i(TMatrixT<Double_t>::kInverted, V_Delta);
        
        /// ------------ Total Chi2
        double chi2 = Delta_Data_Simu * (V_Delta_i * Delta_Data_Simu);
        return chi2;
    };
    
    // wrap chi2 function in a function object for the fit
    // 4 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn_Fit(chi2Function_Fit, N_simu+1);
    ROOT::Fit::Fitter fitter_Fit;
    
    // Initial fit parameters
    Double_t pStart_Fit[N_simu+1];
    pStart_Fit[0]=N_simu;
    for (int i=1; i<N_simu+1; i++){
        pStart_Fit[i]=ini_flux[i-1];
    }
    
    fitter_Fit.SetFCN(fcn_Fit, pStart_Fit);
    fitter_Fit.Config().ParSettings(0).Fix();
    fitter_Fit.Config().ParSettings(1).SetLimits(0.2, 10.);
    fitter_Fit.Config().ParSettings(2).SetLimits(0.2, 10.);
    fitter_Fit.Config().ParSettings(3).SetLimits(0.2, 10.);
    
//    fitter_Fit.Config().ParSettings(1).Fix();
//    fitter_Fit.Config().ParSettings(2).Fix();
//    fitter_Fit.Config().ParSettings(3).Fix();
    
    for (int i=0; i<N_simu; i++){
        fitter_Fit.Config().ParSettings(0).SetName(Names[i].Data());
    }
    
    // do the fit
    cout<<"Fit starts"<<endl;
    bool ok_Fit = fitter_Fit.FitFCN();
    cout<<"Fit End"<<endl;
    
    // test fit
    if(!ok_Fit) {
        Error("fit", "%s", "Fit failed");
    }
    
    const ROOT::Fit::FitResult &result_Fit = fitter_Fit.Result();
    //    result_Fit.Print(std::cout);
    
    return result_Fit;
}

unsigned long long SetNumEvents(TString FileName, bool Bugged_mu=false){
    
    TFile *file = new TFile(FileName);
    TH1D *histoPrim = (TH1D*) file->Get("Primary_gamma_spectrum");
    unsigned long long NumEvents=histoPrim->GetEntries();

    
    return NumEvents;
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

// ---------------------------------------
TH1D* GoodFormat(TH1D* histo, TString title, int New_Nbins, double Bin_MinVal, double Bin_MaxVal, bool ShiftNeutrons = false){
    
    int MinBin = histo->FindBin(Bin_MinVal);
    int MaxBin = histo->FindBin(Bin_MaxVal);
    
    TH1D* new_histo = new TH1D(title+"_cut", title+"_cut", MaxBin-MinBin, Bin_MinVal, Bin_MaxVal);
    
    Double_t Hits_i, Error_i;
    for (int i=0; i<=MaxBin-MinBin; i++){
        Hits_i = histo->GetBinContent(i+MinBin);
        Error_i = histo->GetBinError(i+MinBin);
        
        new_histo->SetBinContent(i+1, Hits_i);
        new_histo->SetBinError(i+1, Error_i);
    }
    
    
    int RebinFact = new_histo->GetNbinsX()/New_Nbins;
    new_histo->Rebin(RebinFact);
    new_histo->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    new_histo->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    
//    cout<<title<<" : "<<new_histo->GetNbinsX()<<"/"<<New_Nbins<<endl;
    return new_histo;
}


// ---------------------------
ROOT::Fit::FitResult Fit_OnRange(vector<Double_t> NumEvents, TH1D* histo_data, TH2D* Var_data, std::vector<TH1D*> histo_simu, double Emin, double Emax, std::vector<Double_t> ini_flux){
    
    Int_t MaxBins=histo_data->FindBin(Emax);
    Int_t MinBins=histo_data->FindBin(Emin);
    
    int New_Nbins=MaxBins-MinBins;
    cout<<"New_Nbins = "<<New_Nbins<<endl;
    
    double deltaE = histo_data->GetBinWidth(1);
    
    TString Data_name =histo_data->GetName();
    
    TVectorT<Double_t> FormatVector(New_Nbins);
    TMatrixT<Double_t> FormatMatrix(New_Nbins, New_Nbins);
    
    TVectorT<Double_t> E_bins(New_Nbins) ; std::vector<TVectorT<Double_t>> Simu_vect={FormatVector, FormatVector, FormatVector};
    TVectorT<Double_t> Data_vect(New_Nbins);
    TMatrixT<Double_t> V_Data(New_Nbins, New_Nbins);
    std::vector<TMatrixT<Double_t>> V_Simu = {FormatMatrix, FormatMatrix, FormatMatrix};

    Cut(histo_data, Var_data, histo_data->GetTitle(), Emin, Emax);

    for (int i =0; i<New_Nbins; i++){
        E_bins(i) = histo_data->GetBinCenter(i+1);
        Data_vect(i) = histo_data->GetBinContent(i+1);
        for (int j=0; j<New_Nbins; j++){
            V_Data(i,j) = Var_data->GetBinContent(i+1, j+1);
        }
    }

    // Prepare simu vectors
    for (int j=0; j<histo_simu.size(); j++){
        histo_simu[j]=GoodFormat(histo_simu[j], Form("hSimu_%d",j), New_Nbins, Emin, Emax);
    }
    
    for (int i =0; i<New_Nbins; i++){
        if (i<New_Nbins){
            for (int j=0; j<histo_simu.size(); j++){
                Simu_vect[j](i) = histo_simu[j]->GetBinContent(i+1);
                V_Simu[j](i,i) = histo_simu[j]->GetBinError(i+1)*histo_simu[j]->GetBinError(i+1);
            }
        }
    }
    
    std::vector<TString> Names = {"U238", "K40", "Th232"};
    
    // Do the fit
    ROOT::Fit::FitResult result_Fit = Fit_Data_Simu(NumEvents, E_bins, Data_vect,
                                                    Simu_vect, V_Simu,
                                                    V_Data, Names,
                                                    ini_flux);
    
    return result_Fit;
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

void Normalize_Simu(TString Filename_in, Double_t flux){
    TFile* filein = new TFile(Filename_in, "READ");
    
    TH1D* h_primaries = (TH1D*) filein->Get("Primary_gamma_spectrum");
    TH1D* h_spectrum_keV = (TH1D*) filein->Get("Edep_gamma_spectrum");
    TH1D* h_spectrum_keV_conv= (TH1D*) filein->Get("Evis_gamma_spectrum");
    
    //Get the number of generated primaries:
    unsigned long long N_entries = h_primaries->GetEntries();
//    TH1D* h_primaries_keV = X_rescale(h_primaries, Form("%s_keV", h_primaries->GetTitle()), 1e3);
//    h_primaries_keV->GetXaxis()->SetTitle("Energy [keV]");
    
    cout<<N_entries<<endl;
    //Calculate the "equivalent time"
    Double_t Eq_time_s = N_entries/(flux*PLANE_SURF); //in seconds
    Double_t Eq_time_d = Eq_time_s/(24*3600.); //in days

    TH1D* h_spectrum_keV_rate = (TH1D*) h_spectrum_keV_conv->Clone();
    
    TString Filename_out = Filename_in;
    TPMERegexp("Conv_Analyzer").Substitute(Filename_out, Form("Normalized_Flux_%0.3f_Conv_Analyzer", flux)); //add Normalized to the name of the file, to create a copy with normalized histograms and not overwrite on the simulation file.
    
    TFile* fileout = new TFile(Filename_out, "RECREATE");
    //Copy before normalizing
    h_primaries->Write("Primary_gamma_spectrum");
    h_spectrum_keV->Write("Edep_gamma_spectrum");
    h_spectrum_keV_conv->Write("Evis_gamma_spectrum");
    
    //Normalize the histograms
    h_spectrum_keV_conv->Scale(1./(MASS_Ge*Eq_time_d*h_spectrum_keV_conv->GetBinWidth(1))); //in ev/d/kg/keV
    h_spectrum_keV_rate->Scale(1./Eq_time_s); //in ev/s
//    h_primaries_keV->Scale(1./(PLANE_SURF*Eq_time_d*h_primaries_keV->GetBinWidth(1))); //in ev/d/keV
    
    h_spectrum_keV_conv->GetYaxis()->SetTitle("Differential Rate [ev kg^{-1} d^{-1} keV^{-1}]");
    h_spectrum_keV_rate->GetYaxis()->SetTitle("Counting Rate [ev s^{-1}]");
//    h_primaries_keV->GetYaxis()->SetTitle("Differential Flux [ev cm^{-2} d^{-1} keV^{-1}]");
    
//    h_primaries_keV->Write("Normalized_Primary_gamma_spectrum");
    h_spectrum_keV->Write("Normalized_Edep_gamma_spectrum");
    h_spectrum_keV_conv->Write("Normalized_Evis_gamma_spectrum_dru");
    h_spectrum_keV_rate->Write("Normalized_Evis_gamma_spectrum_rate");
    
    fileout->Close();
    filein->Close();
    cout<<"Done! Histo written in "<< Filename_out<<endl;
}

// Residual histogram
std::vector<TH1D*> Residuals_histos(TH1D* histo_data, TH1D* Model, Double_t Min, Double_t Max, bool pc = true){
    TH1D* hRes = (TH1D*) histo_data->Clone();
    TString Title =histo_data->GetTitle();
    TH1D* hReducedRes = new TH1D("Reduced_res_"+Title,"Reduced_res_"+Title, 100, -100, 100);
    double Value_i, E_i, Diff_i, Error_i;
    double chi2=0;
    for (int i = 1; i<hRes->GetNbinsX(); i++){
        Value_i = hRes->GetBinContent(i);
        E_i = hRes->GetBinCenter(i);
        if (Value_i>0&&E_i>Min&&E_i<Max) {
            Diff_i = histo_data->GetBinContent(i)-Model->GetBinContent(i);
            Error_i = sqrt(histo_data->GetBinError(i)*histo_data->GetBinError(i)+Model->GetBinError(i)*Model->GetBinError(i));
            if (pc) {
                hRes->SetBinContent(i, Diff_i/Value_i*100);
                hRes->SetBinError(i, Error_i/Value_i*100);
                hReducedRes->Fill(Diff_i/Value_i*100);
            }
            else {
                hRes->SetBinContent(i, Diff_i/Error_i);
                hRes->SetBinError(i, 1);
                hReducedRes->Fill(Diff_i/Value_i);
            }
            chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
        }
        else {
            hRes->SetBinContent(i, 0);
            hRes->SetBinError(i, 0);}
    }

    cout<<"\nChi2/ndf = "<<chi2<<"/"<<(Max-Min)/Model->GetBinWidth(1)<<" = "<<chi2/((Max-Min)/Model->GetBinWidth(1))<<endl;
    hRes->SetLineColor(Model->GetLineColor());
    hRes->SetTitle("");
    if (pc) hRes->GetYaxis()->SetRangeUser(-100, 100);
    else hRes->GetYaxis()->SetRangeUser(-20, 20);
    
    std::vector<TH1D*> histos;
        histos.push_back(hRes);
        histos.push_back(hReducedRes);
    
    return histos;
}

// Plot residuals and histograms and save
void Plot_and_Save(TString fileData_name, TString file_simu_U238, TString file_simu_K40, TString file_simu_Th232, TString Name_Fileout, Double_t Emin, Double_t Emax, ROOT::Fit::FitResult fit){

    //Gat fit results
    double Flux_U238 = fit.Parameter(1);
    double Flux_U238_err = fit.ParError(1);
    double Flux_K40 = fit.Parameter(2);
    double Flux_K40_err = fit.ParError(2);
    double Flux_Th232 = fit.Parameter(3);
    double Flux_Th232_err = fit.ParError(3);
    TMatrixD Cov_Matrix(4,4);
    fit.GetCovarianceMatrix(Cov_Matrix);
   
    TFile *fileData = new TFile(fileData_name);
    
    TH1D* histo_data = (TH1D*) fileData->Get("Cal_spectrum_DRU");
    histo_data->SetTitle("histo_data");
    
    //Simu U234
    TFile *fileU238 = new TFile(file_simu_U238);

    TH1D* histo_U238 = (TH1D*) fileU238->Get("Normalized_Evis_gamma_spectrum_dru");
    histo_U238->SetTitle("histo_U238");
    
    //Simu K40
    TFile *fileK40 = new TFile(file_simu_K40);
    
    TH1D* histo_K40 = (TH1D*) fileK40->Get("Normalized_Evis_gamma_spectrum_dru");
    histo_K40->SetTitle("histo_K40");
    
    //Simu Th232
    TFile *fileTh232 = new TFile(file_simu_Th232);

    TH1D* histo_Th232 = (TH1D*) fileTh232->Get("Normalized_Evis_gamma_spectrum_dru");
    histo_Th232->SetTitle("histo_Th232");
    
    TH1D* Sum= (TH1D*) histo_U238->Clone();
    Sum->Add(Sum, histo_K40);
    Sum->Add(Sum, histo_Th232);
    Sum->SetTitle("Sum_simu");
    
    //Set the correct errors on Sum:
    TVectorD Fluxes(4);
    TVectorD Histo_errors(4);
    TVectorD Cov_mat_Fluxes(4);
    Double_t Error_sum;
    Fluxes(0)=0; //first parameter is number of fitted simulations
    
    Double_t Error_stat;
    for (int i=0; i<Sum->GetNbinsX(); i++){
        // Error from normalization
        Fluxes(1)=1./Flux_U238*histo_U238->GetBinContent(i+1);
        Fluxes(2)=1./Flux_K40*histo_K40->GetBinContent(i+1);
        Fluxes(3)=1./Flux_Th232*histo_Th232->GetBinContent(i+1);
        
        Histo_errors(1)=histo_U238->GetBinError(i+1)*1./Flux_U238;
        Histo_errors(2)=histo_K40->GetBinError(i+1)*1./Flux_K40;
        Histo_errors(3)=histo_Th232->GetBinError(i+1)*1./Flux_Th232;
        
        Cov_mat_Fluxes = Cov_Matrix*Fluxes;
        Error_sum=sqrt(Fluxes*Cov_mat_Fluxes+Histo_errors*Histo_errors);
        
        Sum->SetBinError(i+1, Error_sum);
    }
    
    histo_data = Cut_Histo(histo_data, 0., 3000.);
    histo_data->Rebin(REBIN);
    
    histo_U238->Rebin(REBIN);
    histo_K40->Rebin(REBIN);
    histo_Th232->Rebin(REBIN);
    Sum->Rebin(REBIN);
    
    histo_data->Scale(1./REBIN);
    histo_U238->Scale(1./REBIN);
    histo_K40->Scale(1./REBIN);
    histo_Th232->Scale(1./REBIN);
    Sum->Scale(1./REBIN);
    
    TFile *fileout = new TFile(Name_Fileout, "RECREATE");
    histo_data->Write("Data");
    histo_U238->Write("Simu_U238");
    histo_K40->Write("Simu_K40");
    histo_Th232->Write("Simu_Th232");
    Sum->Write("Simu_Sum");
    
    std::vector<TH1D*> Residuals_fR = Residuals_histos(histo_data, Sum, 0, Emax);
    Residuals_fR[0]->Write("Residuals_fullRange");
    Residuals_fR[1]->Write("Projected_Residuals_fullRange");
    
    int New_Nbins= histo_data->FindBin(Emax)-histo_data->FindBin(Emin);
    
    histo_U238=GoodFormat(histo_U238, "Fitted_histo_U238", New_Nbins, Emin, Emax);
    histo_K40=GoodFormat(histo_K40, "Fitted_histo_K40", New_Nbins, Emin, Emax);
    histo_Th232=GoodFormat(histo_Th232, "Fitted_histo_Th232", New_Nbins, Emin, Emax);
    Sum=GoodFormat(Sum, "Fitted_simu_Sum", New_Nbins, Emin, Emax);
    histo_data=GoodFormat(histo_data, "Fitted_histo_data", New_Nbins, Emin, Emax);
    
    histo_data->Write("Data_fit_Range");
    histo_U238->Write("Simu_U238_fit_Range");
    histo_K40->Write("Simu_K40_fit_Range");
    histo_Th232->Write("Simu_Th232_fit_Range");
    Sum->Write("Simu_Sum_fit_Range");
    
    std::vector<TH1D*> Residuals = Residuals_histos(histo_data, Sum, Emin, Emax);
    Residuals[0]->Write("Residuals");
    Residuals[1]->Write("Projected_Residuals");
    
    fileout->Close();
}


void Fit(TString fileData_name, TString file_simu_U238, TString file_simu_K40, TString file_simu_Th232){
    
    TFile *fileData = new TFile(fileData_name);
    
    //prepare the vectors for the fit
    TH1D* histo_data = (TH1D*) fileData->Get("Cal_spectrum_DRU");
    histo_data->SetTitle("histo_data");
    
    Double_t Num_U238 =SetNumEvents(file_simu_U238);
    cout<<"\nNum events U238 = "<<Num_U238<<endl;
    Double_t Num_K40 = SetNumEvents(file_simu_K40);
    cout<<"Num events K40 = "<<Num_K40<<endl;
    Double_t Num_Th232 =SetNumEvents(file_simu_Th232);
    cout<<"Num events Th232 = "<<Num_Th232<<endl;
    vector<Double_t> NumEvents = {Num_U238, Num_K40, Num_Th232};
    
    histo_data = Cut_Histo(histo_data, 0., 3000.);
    histo_data->Rebin(REBIN);
    histo_data->Scale(1./REBIN);
    
    TH2D* V_data = Stat_Covariance_TH2D(histo_data);
    V_data->SetTitle("V_data");
    
    //U238
    TFile *filein_U238 = new TFile(file_simu_U238);
    TH1D*  histo_simu_U238 = (TH1D*) filein_U238->Get("Evis_gamma_spectrum");
    histo_simu_U238->Rebin(REBIN);

    //Gammas
    TFile *filein_K40 = new TFile(file_simu_K40);
    TH1D*  histo_simu_K40 = (TH1D*) filein_K40->Get("Evis_gamma_spectrum");
    histo_simu_K40->Rebin(REBIN);
    
    //Neutrons
    TFile *filein_Th232 = new TFile(file_simu_Th232);
    TH1D*  histo_simu_Th232 = (TH1D*) filein_Th232->Get("Evis_gamma_spectrum");
    histo_simu_Th232->Rebin(REBIN);
    
    std::vector<TH1D*> histos_simu = {histo_simu_U238, histo_simu_K40, histo_simu_Th232};
    
    //Fit
    Double_t Emin=600.;
    Double_t Emax=2800.;
    std::vector<Double_t> ini_flux = {1.209, 0.302, 0.700};
    ROOT::Fit::FitResult fit = Fit_OnRange(NumEvents, histo_data, V_data, histos_simu, Emin, Emax, ini_flux);

    cout<<"N_simu= "<<fit.Parameter(0)<<endl;
    double Flux_U238 = fit.Parameter(1);
    double Flux_U238_err = fit.ParError(1);
    double Flux_K40 = fit.Parameter(2);
    double Flux_K40_err = fit.ParError(2);
    double Flux_Th232 = fit.Parameter(3);
    double Flux_Th232_err = fit.ParError(3);
    TMatrixD Cov_Matrix(4,4);
    fit.GetCovarianceMatrix(Cov_Matrix);
    TMatrixD Corr_Matrix(4,4);
    fit.GetCorrelationMatrix(Corr_Matrix);
    
    double Chi2_Fit = fit.MinFcnValue();
    int N_par = fit.NPar()-1;
    int N_dof = histo_data->FindBin(Emax)-histo_data->FindBin(Emin) - N_par;
    
    cout << "\n**********************************************************************"<<endl;
    cout << "Fit between "<< Emin << " keV and "<<Emax<< " keV : "<<endl;
    cout << "Chi2/dof = " << Chi2_Fit << "/" << N_dof << " = " << Chi2_Fit/N_dof<<endl;
    cout <<" flux U238 = " << Flux_U238 << " +- " << Flux_U238_err << " ev/cm2/s"<<endl;
    cout <<" flux K40 = " << Flux_K40 << " +- " << Flux_K40_err << " ev/cm2/s"<<endl;
    cout << " flux Th232 = " << Flux_Th232 << " +- " << Flux_Th232_err << " ev/cm2/s"<<endl;
    
    //Calculate total flux and propagate uncertainties:
    const Double_t flux_val[4] = {0, Flux_U238, Flux_K40, Flux_Th232};
    TVectorD Fluxes_values(4, flux_val);
//    Fluxes_values.Print();
    TVectorD Cov_mat_times_Fluxes_values = Cov_Matrix*Fluxes_values;
//    Cov_mat_times_Fluxes_values.Print();
    Double_t Error_sum=sqrt(Fluxes_values*Cov_mat_times_Fluxes_values);
//    cout<<"\n Covariance matrix = "<<endl;
//    Cov_Matrix.Print();
//    cout<<"\n Correlation matrix = "<<endl;
//    Corr_Matrix.Print();
    cout << "\n total flux = " << Flux_U238 + Flux_K40 + Flux_Th232<< " +- " << Error_sum << " ev/cm2/s"<<endl;
    cout << "**********************************************************************"<<endl;
    
    Normalize_Simu(file_simu_U238, Flux_U238);
    Normalize_Simu(file_simu_Th232, Flux_Th232);
    Normalize_Simu(file_simu_K40, Flux_K40);
    
    TPMERegexp("Conv_Analyzer").Substitute(file_simu_U238, Form("Normalized_Flux_%0.3f_Conv_Analyzer", Flux_U238));
    TPMERegexp("Conv_Analyzer").Substitute(file_simu_Th232, Form("Normalized_Flux_%0.3f_Conv_Analyzer", Flux_Th232));
    TPMERegexp("Conv_Analyzer").Substitute(file_simu_K40, Form("Normalized_Flux_%0.3f_Conv_Analyzer", Flux_K40));
    
    TString Out_name_file = fileData_name;
    TPMERegexp("(\\.*\\/*\\w+\\/)*").Substitute(Out_name_file, "");
    TPMERegexp(".root").Substitute(Out_name_file, "");
    Out_name_file=Form("Fit_%s_U238_%0.3f_K40_%0.3f_Th232_%0.3f.root", Out_name_file.Data(), Flux_U238, Flux_K40, Flux_Th232);
    
    Plot_and_Save(fileData_name, file_simu_U238,file_simu_K40, file_simu_Th232, Out_name_file, Emin, Emax, fit);
    cout<<"Results saved in "<<Out_name_file<<endl;
}

