#include <TMath.h>

double sigma1(double x, TMatrixTSym<double> CovMat){
    double sig1 = CovMat(0, 0) + CovMat(0,1)*x + CovMat(1,0)*x + CovMat(1,1)*x*x;
    return sqrt(sig1);
}

class classsigma1 {
public:
    classsigma1(TMatrixTSym<double> CovMat): CovMat(CovMat){
        cout<<"New sigma1 defined with cov matrix = "<<endl;
        CovMat.Print();
    };
    
    double Evaluate(double* x, double* p={0}){
        double x_val = x[0];
        double sig1 = CovMat(0, 0) + CovMat(0,1)*x_val + CovMat(1,0)*x_val + CovMat(1,1)*x_val*x_val;
        return sqrt(sig1);
    };
    
    double EvaluatePlusF(double* x, double* p){
        double x_val = x[0];
        double b = p[0];
        double a = p[1];
        double func_value = b+a*x_val;
        double sig1 = CovMat(0, 0) + CovMat(0,1)*x_val + CovMat(1,0)*x_val + CovMat(1,1)*x_val*x_val;
        return func_value+sqrt(sig1);
    };
    
    double EvaluateMinusF(double* x, double* p){
        double x_val = x[0];
        double b = p[0];
        double a = p[1];
        double func_value = b+a*x_val;
        double sig1 = CovMat(0, 0) + CovMat(0,1)*x_val + CovMat(1,0)*x_val + CovMat(1,1)*x_val*x_val;
        return func_value-sqrt(sig1);
    };
    
    TMatrixTSym<double> CovMat = TMatrixTSym<double>(2);
};

// Functions used to read data and simulations
//--------------------------------------------------------
TH1D* Calibrate_histo(TH1D* histo, Double_t CalibVal, TString NewTitle, Double_t Offset=0, int NBins_E=512, Double_t Emin=0., Double_t Emax=5000.){
    
    int NBins = histo->GetNbinsX();
    
    TH1D* Calib_histo = new TH1D(NewTitle, Form("%s_keVee", histo->GetTitle()), NBins_E, Emin, Emax);
    
    Int_t MCA_i;
    Double_t E_cal, MCA_i_content, MCA_i_Error;
    
    for(int i=0; i<NBins; i++){
        
        MCA_i = histo->GetBinCenter(i);
        MCA_i_content = histo->GetBinContent(i)*histo->GetBinWidth(i);
        MCA_i_Error = histo->GetBinError(i)*histo->GetBinWidth(i);
        
        E_cal = MCA_i*CalibVal+Offset;
        int Bin_Ecal = Calib_histo->FindBin(E_cal);
        
        Double_t prev_val_Ecal = Calib_histo->GetBinContent(Bin_Ecal);
        Double_t prev_err_Ecal = Calib_histo->GetBinError(Bin_Ecal);
        
        Calib_histo->SetBinContent(Bin_Ecal, prev_val_Ecal+MCA_i_content);
        Calib_histo->SetBinError(Bin_Ecal, sqrt(prev_err_Ecal*prev_err_Ecal+MCA_i_Error*MCA_i_Error));
    }
    
    Calib_histo->Scale(1./Calib_histo->GetBinWidth(1));
    Calib_histo->GetXaxis()->SetTitle("Calibrated Energy [keVee]");
    Calib_histo->GetYaxis()->SetTitle("Differential rate (ev.d^{-1}.keVee^{-1})");

    return Calib_histo;
}

TH1D* Calib(TH1D* histo, Double_t CalibVal, TString NewTitle, Double_t Offset=0, int NBins_E=512, Double_t Emin=0., Double_t Emax=5000.){
    
    int Nbins = histo->GetNbinsX();
    
    double histo_binsize = (histo->GetBinWidth(1));
    double histo_min = histo->GetBinLowEdge(1);
    double histo_max = histo->GetBinLowEdge(Nbins)+histo_binsize;
    
    double New_bins[NBins_E];
    double New_counts[NBins_E];
    double New_errors[NBins_E];
    
    double histo_bin_min = histo_min;
    double Bin_width = (Emax-Emin)/NBins_E;
    int index_bin_E=0;
    double histo_bin_max, histo_bin_max_next, E_bin_max, counts_bin, counts_bin_next, histo_bin_cut, error_bin, error_bin_next, new_error;
    
    //------------ Set New_counts to 0 and define newbins----------------//
    for (int i =0; i< NBins_E; i++){
        New_bins[i]=Emin+i*Bin_width;
        New_counts[i] = 0;
        New_errors[i] = 0;
    }

    for (int i=0; i< Nbins; i++){
        histo_bin_max = histo->GetBinLowEdge(i)+histo_binsize;
        E_bin_max = histo_bin_max*CalibVal+Offset;
        counts_bin = histo->GetBinContent(i);
        error_bin = histo->GetBinError(i);

        if ((E_bin_max>=Emin)&&(E_bin_max<=Emin+index_bin_E*Bin_width)&&(Emin+index_bin_E*Bin_width<Emax)){
            New_counts[index_bin_E] += counts_bin*(histo_bin_max-histo_bin_min)/histo_binsize;
            new_error=error_bin*sqrt((histo_bin_max-histo_bin_min)/histo_binsize);
            New_errors[index_bin_E] = sqrt(New_errors[index_bin_E]*New_errors[index_bin_E] + new_error*new_error);
            
            histo_bin_min = histo_bin_max;
            histo_bin_max_next = histo->GetBinLowEdge(i)+2*histo_binsize;
            histo_bin_cut = (Emin+index_bin_E*Bin_width-Offset)/CalibVal;
            
            if (histo_bin_cut<=histo_bin_max_next){
                counts_bin_next = histo->GetBinContent(i+1);
                error_bin_next = histo->GetBinError(i+1);
                New_counts[index_bin_E] += counts_bin_next*(histo_bin_cut-histo_bin_min)/histo_binsize;
                new_error=error_bin_next*sqrt((histo_bin_cut-histo_bin_min)/histo_binsize);
                New_errors[index_bin_E] = sqrt(New_errors[index_bin_E]*New_errors[index_bin_E] + new_error*new_error);
                index_bin_E+=1;
                histo_bin_min = histo_bin_cut;
                
            }
        }

        else if ((E_bin_max>=Emin)&&(Emin+index_bin_E*Bin_width<Emax)){
            while ((Emin+index_bin_E*Bin_width)<= E_bin_max){
                histo_bin_cut = (Emin+index_bin_E*Bin_width-Offset)/CalibVal;
                index_bin_E+=1;
            }
            New_counts[index_bin_E] += counts_bin*(histo_bin_max-histo_bin_min)/histo_binsize;
            new_error=error_bin_next*sqrt((histo_bin_max-histo_bin_min)/histo_binsize);
            New_errors[index_bin_E] = sqrt(New_errors[index_bin_E]*New_errors[index_bin_E] + new_error*new_error);
            histo_bin_min = histo_bin_max;
        }
    }
    
    TH1D *histo_cal = new TH1D(NewTitle, NewTitle, NBins_E-1, New_bins);
    for (int i =0; i<NBins_E-1; i++){
        histo_cal-> SetBinContent(i+1, New_counts[i]);
        histo_cal-> SetBinError(i+1, New_errors[i]);
    }
    
    histo_cal->Scale(1./histo_cal->GetBinWidth(1));
    histo_cal->GetXaxis()->SetTitle("Calibrated Energy [keVee]");
    histo_cal->GetYaxis()->SetTitle("Differential rate (ev.d^{-1}.keVee^{-1})");
    
    return histo_cal;
}


//---------------------------------------------
std::vector<TVectorD> CalibVector(TVectorD ToCalib, TVectorD TC_bins, TMatrixD Errors, Double_t CalibVal, Double_t Offset, TVectorD Ebins){
    
    int Nbins = ToCalib.GetNrows();
    
    double histo_min = TC_bins(0);
    double histo_binsize = TC_bins(1)-TC_bins(0);
    double histo_max = TC_bins(Nbins-1)+histo_binsize;
    
    int NBins_E = Ebins.GetNrows();
    
    TVectorD New_counts(NBins_E);
    TVectorD New_errors(NBins_E);
    
    double histo_bin_min = histo_min;
    double Bin_width = (Ebins(1)-Ebins(0));
    double Emin = Ebins(0);
    double Emax = Ebins(NBins_E-1)+Bin_width;
    int index_bin_E=0;
    double histo_bin_max, histo_bin_max_next, E_bin_max, counts_bin, counts_bin_next, histo_bin_cut, error_bin, error_bin_next, new_error;
    
    //------------ Set New_counts to 0 and define newbins----------------//
    for (int i =0; i<NBins_E; i++){
        New_counts(i) = 0;
        New_errors(i) = 0;
    }

    for (int i=0; i<Nbins; i++){

        histo_bin_max = TC_bins(i)+histo_binsize;
        E_bin_max = histo_bin_max*CalibVal+Offset;
        
        counts_bin = ToCalib(i);
        error_bin = Errors(i,i);

        if ((E_bin_max>=Emin)&&(E_bin_max<=Emin+index_bin_E*Bin_width)&&(Emin+index_bin_E*Bin_width<Emax)){

            New_counts(index_bin_E) += counts_bin*(histo_bin_max-histo_bin_min)/histo_binsize;
            new_error=error_bin*sqrt((histo_bin_max-histo_bin_min)/histo_binsize);
            New_errors(index_bin_E) = sqrt(New_errors(index_bin_E)*New_errors(index_bin_E) + new_error*new_error);

            histo_bin_min = histo_bin_max;
            histo_bin_max_next = TC_bins(i)+2*histo_binsize;
            histo_bin_cut = (Emin+index_bin_E*Bin_width-Offset)/CalibVal;

            if ((histo_bin_cut<=histo_bin_max_next)&&(i+1<Nbins)){
                counts_bin_next = ToCalib(i+1);
                error_bin_next = Errors(i+1, i+1);
                New_counts(index_bin_E) += counts_bin_next*(histo_bin_cut-histo_bin_min)/histo_binsize;
                new_error=error_bin_next*sqrt((histo_bin_cut-histo_bin_min)/histo_binsize);
                New_errors(index_bin_E) = sqrt(New_errors(index_bin_E)*New_errors(index_bin_E) + new_error*new_error);
                index_bin_E+=1;
                histo_bin_min = histo_bin_cut;

            }
        }

        else if ((E_bin_max>=Emin)&&(Emin+index_bin_E*Bin_width<Emax)){
            while ((Emin+index_bin_E*Bin_width)<= E_bin_max){
                histo_bin_cut = (Emin+index_bin_E*Bin_width-Offset)/CalibVal;
                index_bin_E+=1;
            }
            New_counts(index_bin_E) += counts_bin*(histo_bin_max-histo_bin_min)/histo_binsize;
            new_error=error_bin_next*sqrt((histo_bin_max-histo_bin_min)/histo_binsize);
            New_errors(index_bin_E) = sqrt(New_errors(index_bin_E)*New_errors(index_bin_E) + new_error*new_error);
            histo_bin_min = histo_bin_max;
        }
    }
    
    New_counts=(1./Bin_width)*New_counts;
    New_errors=(1./Bin_width)*New_errors;
    
    std::vector<TVectorD> Cal_v_err;
    Cal_v_err.push_back(New_counts);
    Cal_v_err.push_back(New_errors);
    
    return Cal_v_err;
}

//--------------------------------------------------------
double Get_Sum_of_VectorElmts(TVectorD vector){
    int N = vector.GetNrows();
    double sum = 0;
    for (int i =0 ; i<N ; i++){
        sum+=vector(i);
    }
    return sum;
}

//--------------------------------------------------------
//Standard deviation of a vector
Double_t Sigma(vector<Double_t> Vector){
    int num_val = Vector.size();
    Double_t Mean=0;
    for (int i=0; i<num_val; i++){
        Mean+=Vector[i];
    }
    Mean= Mean/num_val;
    
    Double_t sigma=0;
    for (int i=0; i<num_val; i++){
        sigma+=(Vector[i]-Mean)*(Vector[i]-Mean);
    }
    sigma= sigma/(num_val-1);
    sigma = sqrt(sigma);
    return sigma;
}

//Smooth Step function
//--------------------------------------------------------------
double cdf(double *x, double *par){
    double mean=par[0];
    double sigma = par[1];
    
    double xx =x[0];
    return ROOT::Math::normal_cdf(xx, sigma, mean);
}

//--------------------------------------------------------
TH1D* Put_data_in_histo(TString filename, double real_time, TString Title){
    ifstream in;
    in.open(filename, ios::in);
    cout<<"\n*******************************************************"<<endl;
    cout<<filename<<endl;
    Float_t time=real_time/(24*3600.);
    
    Int_t x,y, MCA_chan;
    TH1D* h1 = new TH1D(Title,Title,512,-0.5,511.5);
    h1->GetXaxis()->SetTitle("Pulse Height [MCA a.u.]");
    h1->GetYaxis()->SetTitle("Counts per day per MCA channel");
    
    while (1) {
        in >> x >> y;
        if (!in.good()) break;
        
        h1->SetBinContent(x+1, y);
        h1->SetBinError(x+1, sqrt(y));
    }
    MCA_chan=h1->GetBinWidth(1);
    h1->Scale(1./(MCA_chan*time));
    return h1;
}

//--------------------------------------------------------
TH1D* Put_simu_in_histo_nq(TString filename, TString Title, int nBin, double Emin, double Emax, Double_t NumEv, Double_t Flux, Double_t MVthr=3.){
    TFile* Infile = new TFile(filename);
    
    Float_t ECr, EMV;
    TRandom gen(0);
    
    Double_t Plane_surf = 40.*40.; //cm2
    
    Double_t Time = NumEv/(Plane_surf*Flux); //s
    
    TTree *InTree = (TTree*) Infile->Get("NTupleAlgoDir/Pri_Ntuple");
    Int_t entries= InTree->GetEntries();
    printf(" Entries:      %10d\n",entries);
    printf(" Generated primaries:  %.0f\n",NumEv);
    printf(" Equivalent time:      %.0f s (%.3f d)\n",Time,Time/86400.);

    TH1D *histo = new TH1D(Title,Title,nBin,Emin,Emax);
    histo->GetXaxis()->SetTitle("Energy [keVee]");
    histo->GetYaxis()->SetTitle("Counts.d^{-1}.keVee^{-1}");

    InTree->SetBranchAddress("EDep", &ECr);
    InTree->SetBranchAddress("EMV0", &EMV);
    for (Int_t pk=0; pk<entries; pk++) {
        InTree->GetEntry(pk);
        
        Float_t Ep = EMV;
        
        if (Ep>0){
            Float_t sp = 0.1008*sqrt(1./Ep*(1.+1./Ep))/2.3548;      // (%)
            Ep *= (1. + sp*gen.Gaus())*1000.;
            if (Ep < MVthr*1000) continue;
        }
        
        Float_t E = ECr*1000.;                            // keV
        Float_t Enew = gen.Gaus(E, 0.077*sqrt(3500)*sqrt(E)); //Neutron peak assumed to be at 3500 keVee
        histo->Fill(Enew);
    }
    Double_t Ei = histo->GetBinWidth(1);
    histo->Scale(1./(Time/86400.)/Ei);
    
    return histo;
}

//--------------------------------------------------------
//This function applies a quenching factor
TH1D* Put_simu_in_histo(TString filename, TString Title, int nBin, double Emin, double Emax, Double_t NEvt, Double_t Flux =0.0134, Double_t MVthr=3., bool Coinc=true, Double_t resol = 0.07){
    
    TFile* Nfile = new TFile(filename);

    Float_t ECr, EMV;
    TRandom gen(0);
    
    Double_t NTime = NEvt/40./40./Flux;     // s
    
    
    TTree *Nnp = (TTree*) Nfile->Get("NTupleAlgoDir/Pri_Ntuple");
    TTree *Nns = (TTree*) Nfile->Get("NTupleAlgoDir/Sec_Ntuple");

    Int_t Np_entries= Nnp->GetEntries();
    Int_t Ns_entries= Nns->GetEntries();
    printf("  Primary Entries:      %10d\n",Np_entries);
    printf("  Secondary Entries:    %10d\n",Ns_entries);
    printf("  Generated primaries:  %.0f\n",NEvt);
    printf("  Equivalent time:      %.0f s (%.3f d)\n",NTime,NTime/86400.);
    
    TH1D *neut = new TH1D(Title,Title,nBin,Emin,Emax);
    
    neut->GetXaxis()->SetTitle("Energy [keVee]");
    neut->GetYaxis()->SetTitle("Counts.d^{-1}.keVee^{-1}");

    Nnp->SetBranchAddress("EDep", &ECr);
    Nnp->SetBranchAddress("EMV0", &EMV);

    Float_t NpNSec,NpDelay, NsPDG,NsPidM,NsPid;
    Nnp->SetBranchAddress(  "NSec", &NpNSec);
    Nnp->SetBranchAddress( "Delay", &NpDelay);
    Nns->SetBranchAddress(   "PDG", &NsPDG);
    Nns->SetBranchAddress(  "PidM", &NsPidM);
    Nns->SetBranchAddress(   "Pid", &NsPid);

    // Calculate a constant q.f. from the position, on the gamma calibrated Energy scale, of the thermal neutron peak (expected at 4.78 MeVee)
    Float_t qf = 3500/4780.; //Quenching factor Emin + 326/512(Position of neutron peak/number of bins)*(Emax-Emin)/4780(Qvalue the neutron capture)
    
    // Loop on primary particles pk
    Int_t pk=0, sk=0, npCheck=0, nsCheck=0, isQF=0, nQuenched=0;
    
    while (pk < Np_entries) {
        
        // Loop on secondary particles sk (defined as those particles entering in the LiI(Eu) crystal)
        if (sk < Ns_entries) Nns->GetEntry(sk);
        
        // Get the "id" of the primary particles associated to the secondary particle sk
        Int_t id = ((Int_t) (NsPidM+0.5))*1000000 + ((Int_t) (NsPid+0.5));
        
         // if the primary particles is pk (the one under examination in primary loop), set isQF=1 in the cases when the
        // secondary particle sk is a neutron or antineutron
        // ==> i.e., if the secondary particles entering in the LiI(Eu) crystal is a neutron, the deposited energy is due
        // to a neutron interaction in LiI(Eu) detector and the quenching factor has to be applied (isQF=1, --> true).
        if(id==pk && sk < Ns_entries) {
            if (NsPDG==2112 || NsPDG==-2112) isQF=1;
            sk++;
            nsCheck++;
        }
        
        // otherwise, no further secondary particles are associated to pk, I can fill the spectrum with the deposited
        // energy in LiI(Eu) detector (branch EDep in primary ntuple), applying or not the quenching factor
        else {
            Nnp->GetEntry(pk);
            // it is only a consinstency check (verify that the number of analysed secondary particles is the expected one)
            // for the pk primary particle under examination
            if (NpNSec != nsCheck) printf(" Primary entry %d: NSec error NSec=%d, nsCheck=%d\n", pk, (int) NpNSec, nsCheck);
    
            // It happens sometime that the primary particles is not associated to any secondary particles (NpNSec=0);
            // Almost always, it is due to very slow neutrons that after entering the LiI(Eu) crystal need long time before
            // having an interaction, and therefore are classified as new, delayed events (NpDelay>0) and loose their correlation
            // with the secondary particle ntuple. Also these events have to be classified as neutron events ==> isQF=1
            // (this point is not easy to explain, maybe we can have a short discussion about it)
            if (NpNSec==0 && NpDelay>0) isQF=1;        // Sicuramente neutroni
            if (isQF == 1) nQuenched++;
            
            // the E deposition in the MV is spread accounting for energy resolution resolution
            Float_t Ep = EMV;
        
            if (Ep>0){
                Float_t sp = 0.1008*sqrt(1./Ep*(1.+1./Ep))/2.3548;      // (%)
                Ep *= (1. + sp*gen.Gaus())*1000.;
            }
            
            if (Coinc){
                if (Ep > MVthr*1000) {
                    Float_t E = ECr*1000.;                                  // keV

                    if (isQF==1) E *= qf;
                    Double_t Eresol = gen.Gaus(E, resol*sqrt(E)); //Neutron peak assumed to be at 3500 keVee                           // keV
                    neut->Fill(Eresol);
                } // else printf("Ep = %f keV\n",EMV);
            }
            
            else {
                Float_t E = ECr*1000.;                                  // keV

                if (isQF==1) E *= qf;
                Double_t Eresol = gen.Gaus(E, resol*sqrt(E)); //Neutron peak assumed to be at 3500 keVee                           // keV
                neut->Fill(Eresol);
            }
            npCheck++;
            nsCheck=0;
            isQF=0;
            pk++;
        }
    }

    printf(" Nb of primary analysed: %d\n", npCheck);
    printf(" Nb of events quenched : %d\n", nQuenched);
    
    Double_t Ei = neut->GetBinWidth(1);
    neut->Scale(1./(NTime/86400.)/Ei);
    return neut;
}

//--------------------------------------------------------
// Residual histogram
TH1D* Residuals(TH1D* histo_data, TF1* Model, Double_t Min, Double_t Max, bool pc = true){
    TH1D* hRes = (TH1D*) histo_data->Clone();
    double Value_i, E_i, Diff_i, Error_i;
    double chi2=0;
    for (int i = 1; i<hRes->GetNbinsX(); i++){
        Value_i = hRes->GetBinContent(i);
        E_i = hRes->GetBinCenter(i);
        if (Value_i>0 && E_i>Min && E_i<Max) {
            Diff_i = histo_data->GetBinContent(i)-Model->Eval(E_i);
            Error_i = histo_data->GetBinError(i);
            if (pc) {
                hRes->SetBinContent(i, Diff_i/Value_i*100);
                hRes->SetBinError(i, Error_i/Value_i*100);
            }
            else {
                hRes->SetBinContent(i, Diff_i/Error_i);
                hRes->SetBinError(i, 1);
            }
            chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
        }
        else
        {hRes->SetBinContent(i, 0);
        hRes->SetBinError(i, 0);}
    }

    //cout<<"\nChi2/ndf = "<<chi2<<"/"<<Model->GetNDF()<<" = "<<chi2/Model->GetNDF()<<endl;
    hRes->SetLineColor(Model->GetLineColor());
    hRes->SetTitle("");
    if (pc) hRes->GetYaxis()->SetRangeUser(-100, 100);
    else hRes->GetYaxis()->SetRangeUser(-20, 20);
    return hRes;
}

//--------------------------------------------------------
// Residual histogram
std::vector<TH1D*> Residuals_histos(TH1D* histo_data, TH1D* Model, Double_t Min, Double_t Max, bool pc = true){
    TH1D* hRes = (TH1D*) histo_data->Clone();
    TString Title = histo_data->GetTitle();
    TH1D* hReducedRes = new TH1D("Reduced_res_"+Title,"Reduced_res_"+Title, 50, -10, 10);
    
    double Value_i, E_i, Diff_i, Error_i;
    double chi2=0;
    int ndf =0;
    for (int i = 1; i<hRes->GetNbinsX(); i++){
        Value_i = hRes->GetBinContent(i);
        E_i = hRes->GetBinCenter(i);
        if (Value_i>0&&E_i>Min&&E_i<Max) {
            Diff_i = histo_data->GetBinContent(i)-Model->GetBinContent(i);
            Error_i = sqrt(histo_data->GetBinError(i)*histo_data->GetBinError(i)+Model->GetBinError(i)*Model->GetBinError(i));
            hReducedRes->Fill(Diff_i/Error_i);
            
            if (pc) {
                hRes->SetBinContent(i, Diff_i/Value_i*100);
                hRes->SetBinError(i, Error_i/Value_i*100);
            }
            else {
                hRes->SetBinContent(i, Diff_i/Error_i);
                hRes->SetBinError(i, 1);
            }
            chi2+=(Diff_i*Diff_i)/(Error_i*Error_i);
            ndf+=1;
        }
        else {
            hRes->SetBinContent(i, 0);
            hRes->SetBinError(i, 0);}
    }

    cout<<"\nChi2/ndf = "<<chi2<<"/"<<ndf-4<<" = "<<chi2/(ndf-4)<<endl;
    hRes->SetLineColor(Model->GetLineColor());
    hRes->SetTitle("");
    
    hReducedRes->SetFillColor(Model->GetLineColor());
        hReducedRes->SetLineColor(Model->GetLineColor());
        hReducedRes->SetTitle("");
    
    if (pc) hRes->GetYaxis()->SetRangeUser(-100, 100);
    else hRes->GetYaxis()->SetRangeUser(-20, 20);
    
    //hReducedRes->GetXaxis()->SetRangeUser(-5, 5);
    std::vector<TH1D*> histos;
    histos.push_back(hRes);
    histos.push_back(hReducedRes);
    
    return histos;
}

//--------------------------------------------------------
//Fuction to cut a vector between loweredge and upperedge
TVectorD Cut_Vector(TVectorD Vector_to_cut, TVectorD X_axis, double X_low, double X_up) {
    
    double Binwidth = X_axis(1)-X_axis(0);
    int NBins = X_axis.GetNrows();
    
    // Find the index of the limits
    int index_low, index_up;
    for (int i=0; i<NBins-1; i++){
        if ((X_low>=X_axis(i))&&(X_low< X_axis(i+1))){
            index_low = i+1;
            break;
        }
    }

    for (int i=NBins-1; i>0; i--){
        if ((X_up<X_axis(i))&&(X_up>=X_axis(i-1))){
            index_up = i;
            break;
        }
    }

    int New_Nbins = index_up-index_low;
    TVectorD Cut_Vector(New_Nbins);
    
    int i_vector_to_cut;
    for(int i = 0; i<New_Nbins ; i++) {
        i_vector_to_cut = index_low + i;
        Cut_Vector(i)=Vector_to_cut(i_vector_to_cut);
    }
    
    return Cut_Vector;
}

//--------------------------------------------------------
//Fuction to cut a Matrix between loweredge and upperedge
TMatrixD Cut_Matrix(TMatrixD Matrix_to_cut, TVectorD X_axis, double X_low, double X_up, TVectorD Y_axis, double Y_low, double Y_up) {
    
    double XBinwidth = X_axis(1)-X_axis(0);
    int XNBins = X_axis.GetNrows();
    
    double YBinwidth = Y_axis(1)-Y_axis(0);
    int YNBins = Y_axis.GetNrows();
    
    // Find the index of the limits
    int X_index_low, X_index_up;
    for (int i=0; i<XNBins-1; i++){
        if ((X_low>=X_axis(i))&&(X_low< X_axis(i+1))){
            X_index_low = i+1;
            break;
        }
    }
    
    for (int i=XNBins-1; i>0; i--){
        if ((X_up<X_axis(i))&&(X_up>=X_axis(i-1))){
            X_index_up = i;
            break;
        }
    }

    // Find the index of the limits
    int Y_index_low, Y_index_up;
    for (int i=0; i<YNBins-1; i++){
        if ((Y_low>=Y_axis(i))&&(Y_low< Y_axis(i+1))){
            Y_index_low = i+1;
            break;
        }
    }
    
    for (int i=YNBins-1; i>0; i--){
        if ((Y_up<Y_axis(i))&&(Y_up>=Y_axis(i-1))){
            Y_index_up = i;
            break;
        }
    }
    
    int XNew_Nbins = X_index_up-X_index_low;
    int YNew_Nbins = Y_index_up-Y_index_low;
    TMatrixD Cut_Matrix(YNew_Nbins, XNew_Nbins);
    
    int i_vector_to_cut, j_vector_to_cut;
    for(int i = 0; i<YNew_Nbins ; i++) {
        i_vector_to_cut = Y_index_low + i;
        for(int j = 0; j<XNew_Nbins ; j++) {
            j_vector_to_cut = Y_index_low + j;
            Cut_Matrix(i, j) = Matrix_to_cut(i_vector_to_cut, j_vector_to_cut);
        }
    }
    
    return Cut_Matrix;
}

//--------------------------------------------------------
//Fit function on vectors
ROOT::Fit::FitResult Do_Conv(TVectorD &MCA_bins, TVectorD &E_bins, TVectorD &Simu_vect, TVectorD &Data_vect, TMatrixD &V_Simu, TMatrixD &V_Data, TString &Simu_name, TString &Data_name, Double_t Emin_fit, Double_t Emax_fit, double Gain_ini = 9., double Norm_ini = 0.8, double Offset_ini =0., double Resol_ini = 4.) {
    
    int N_MCA = MCA_bins.GetNrows();
    int N_E = E_bins.GetNrows();
    //TVectorD Delta_Data_SimuConv(N_MCA); //difference between Data and alpha*Simu
    
    /* Initial guess: find the alpha for normalizing Data1 to Data2 */
    double norm_Data = Get_Sum_of_VectorElmts(Data_vect);
    double norm_Simu = Get_Sum_of_VectorElmts(Simu_vect);
    
    double alpha_Simu_to_Data = norm_Data/norm_Simu;
    //TMatrixD V_Delta_i(N_MCA, N_MCA);
    //TMatrixD V_Delta(N_MCA, N_MCA);
    //TMatrixD A(N_MCA, N_E);

    int N_zeros_final = 0;
    
    TMatrixD Cal(N_E, N_MCA);
    TMatrixD Cal_T(N_MCA, N_E);
    TMatrixD A(N_E, N_E);
    TMatrixD A_T(N_E, N_E);
    
//    TCanvas* c1 = new TCanvas();
//    TGraph* gData = new TGraph(MCA_bins, Data_vect);
//    gData->Draw();
//    TCanvas* c2 = new TCanvas();
//    TGraph* gSimu = new TGraph(E_bins, Simu_vect);
//    gSimu->Draw();

    std::vector<TVectorD> Cal_vect;
    
    //cf. https://root.cern.ch/doc/master/fitCircle_8C.html
    auto chi2Function_ConvFit = [&](const Double_t *par) {
        //minimisation function computing the sum of squares of residuals
        // looping at the graph points
        double alpha = par[0];
        double Gain = par[1];
        double Offset = par[2];
        double Resol = par[3];
        
        Cal_vect = CalibVector(Data_vect, MCA_bins, V_Data, Gain, Offset, E_bins);
        TVectorD Cal_Data_vect = Cal_vect[0];
        TVectorD Err_Cal_Data_vect = Cal_vect[0];
        
        //---- construction de V_data
        TMatrixD Cal_V_Data(N_E, N_E);
        for (int i=0; i<N_E; i++){
            Cal_V_Data(i,i)= Err_Cal_Data_vect(i);
        }
        
        //------ Resolution on the Simulation
        double Binwidth = E_bins(1)-E_bins(0);
        double E_j, E_i, sigma;
        for (int i =0; i<N_E ; i++){
            E_i = E_bins(i)+Binwidth/2;
            for (int j =0; j< N_E; j++){
                E_j=E_bins(j)+Binwidth/2;
                sigma=Resol*sqrt(E_j);
                A(i,j)=TMath::Gaus(E_i, E_j, sigma, true);
            }
        }

        //------ Normalize A
        double Integral_gauss=0;
        double Emin = E_bins(0);
        double Emax = E_bins(N_E-1)+Binwidth;
        TVectorD Colj_A(N_E);

        for (int j=0; j<N_E ; j++){
            E_j=E_bins(j)+Binwidth/2;
            sigma=Resol*sqrt(E_j);
            Integral_gauss=ROOT::Math::normal_cdf(Emax, sigma, E_j)-ROOT::Math::normal_cdf(Emin, sigma, E_j);
            
            for (int i=0; i<N_E; i++){
                Colj_A(i)=A(i,j);
            }
            
            for (int i=0; i<N_E; i++){
                A(i,j)*=Integral_gauss/Get_Sum_of_VectorElmts(Colj_A);
            }
        }

        //Resolution du spectre de simu
        TVectorD Resol_Simu_vect = (A*Simu_vect);

        TMatrixD Af_T(TMatrixD::kTransposed, A);
        TMatrixD V_Resol_Simu = A*V_Simu*Af_T;
        
        //Cut vectors in fit range
        TVectorD Fit_range_Resol_Simu_vect = Cut_Vector(Resol_Simu_vect, E_bins, Emin_fit, Emax_fit);
        TMatrixD Fit_range_V_Resol_Simu = Cut_Matrix(V_Resol_Simu, E_bins, Emin_fit, Emax_fit, E_bins, Emin_fit, Emax_fit);
        TVectorD Fit_range_E_bins = Cut_Vector(E_bins, E_bins, Emin_fit, Emax_fit);
        
        TVectorD Fit_range_Cal_Data_vect = Cut_Vector(Cal_Data_vect, E_bins, Emin_fit, Emax_fit);
        TMatrixD Fit_range_Cal_V_Data = Cut_Matrix(Cal_V_Data, E_bins, Emin_fit, Emax_fit, E_bins, Emin_fit, Emax_fit);
        
        int fit_range_N_E = Fit_range_E_bins.GetNrows();
        
        //----------- Projection in a base where V_data=/=0
        int N_zeros = 0;
        for (int i =0; i<fit_range_N_E; i++){
                N_zeros+=(Fit_range_Cal_Data_vect(i)==0); //count the number of zero in Data_vect
        }

        int N_E_Proj = fit_range_N_E-N_zeros;
        TMatrixD Proj(N_E_Proj, fit_range_N_E);
        
        int i_E_Proj = 0;
        for (int j=0; j<fit_range_N_E; j++){
            if (Fit_range_Cal_Data_vect(j)!=0){
                Proj(i_E_Proj, j)=1;
                i_E_Proj+=1;
            }
        }
        
        //---- Projection du spectre data
        TVectorD Proj_Cal_Data_vect = Proj*Fit_range_Cal_Data_vect;
        
        //---- projection de V_data
        TMatrixD Proj_T(TMatrixD::kTransposed, Proj);
        TMatrixD Proj_Cal_V_Data = Proj*Fit_range_Cal_V_Data*Proj_T;
        
        //---- projection simu
        TVectorD Proj_Simu_vect = Proj*Fit_range_Resol_Simu_vect;
        
        //---- projection de V_simu
        TMatrixD Proj_V_Simu = Proj*Fit_range_V_Resol_Simu*Proj_T;
        
        //---- projection de Energy bins
        TVectorD Proj_E_bins = Proj*Fit_range_E_bins;

        TMatrixD V_Delta(N_E_Proj, N_E_Proj);
        TVectorD Delta_Data_SimuConv(N_E_Proj);
        
        Proj_Simu_vect = alpha * Proj_Simu_vect ;
        
        Delta_Data_SimuConv = Proj_Cal_Data_vect - Proj_Simu_vect ;
        
        V_Delta = Proj_Cal_V_Data + alpha*alpha*Proj_V_Simu;

        TMatrixD V_Delta_i(TMatrixD::kInverted, V_Delta);
        
        int N_dof = Proj_Cal_Data_vect.GetNrows() - 4;
        double chi2 = Delta_Data_SimuConv * (V_Delta_i * Delta_Data_SimuConv);
//        cout<<chi2/N_dof<<endl;
        return chi2/N_dof;
    };
    
    // wrap chi2 function in a function object for the fit
    // 4 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn_ConvFit(chi2Function_ConvFit, 4);
    ROOT::Fit::Fitter fitter_ConvFit;
    
    double pStart_ConvFit[4] = {Norm_ini, Gain_ini, Offset_ini, Resol_ini};
    fitter_ConvFit.SetFCN(fcn_ConvFit, pStart_ConvFit);
    fitter_ConvFit.Config().ParSettings(0).SetName("Norm");
    fitter_ConvFit.Config().ParSettings(0).SetLimits(0.7, 5.);
    //fitter_ConvFit.Config().ParSettings(0).Fix();
    fitter_ConvFit.Config().ParSettings(1).SetName("Gain (keV/MCA)");
//    fitter_ConvFit.Config().ParSettings(1).Fix();
    fitter_ConvFit.Config().ParSettings(1).SetLimits(6., 12.);
    fitter_ConvFit.Config().ParSettings(2).SetName("Offset");
    //fitter_ConvFit.Config().ParSettings(2).Fix();
    fitter_ConvFit.Config().ParSettings(2).SetLimits(-20, 20);
    fitter_ConvFit.Config().ParSettings(3).SetName("Resolution");
    fitter_ConvFit.Config().ParSettings(3).SetLimits(3., 10.);
    //fitter_ConvFit.Config().ParSettings(3).Fix();
    
    // do the fit
    cout<<"First fit ..."<<endl;
    bool ok_ConvFit = fitter_ConvFit.FitFCN();
    if(!ok_ConvFit) {
        Error("fit", "%s", Form("%s / %s Conv Fit failed", Simu_name.Data(), Data_name.Data()));
    }
    
    const ROOT::Fit::FitResult &result_ConvFit = fitter_ConvFit.Result();
    //result_ConvFit.Print(std::cout);
    
    double Norm = result_ConvFit.Parameter(0);
    double Norm_err = result_ConvFit.ParError(0);
    double Gain_val = result_ConvFit.Parameter(1);
    double Gain_err = result_ConvFit.ParError(1);
    double Offset_val = result_ConvFit.Parameter(2);
    double Offset_err = result_ConvFit.ParError(2);
    double Resol_val = result_ConvFit.Parameter(3);
    double Resol_err = result_ConvFit.ParError(3);
    
    //do it again---
    double pStart_ConvFit_2nd[4] = {Norm, Gain_val, Offset_val, Resol_val};

//    fitter_ConvFit.Config().ParSettings(0).SetLimits(0.6, 5.);
    //fitter_ConvFit.Config().ParSettings(0).Fix();
    //fitter_ConvFit.Config().ParSettings(1).Fix();
//    fitter_ConvFit.Config().ParSettings(1).SetLimits(0.8*Gain_val, 1.1*Gain_val);
    //fitter_ConvFit.Config().ParSettings(2).Fix();
//    if (Offset_val>0) {fitter_ConvFit.Config().ParSettings(2).SetLimits(0.8*Offset_val, 1.1*Offset_val);}
//    if (Offset_val<=0) {fitter_ConvFit.Config().ParSettings(2).SetLimits(1.1*Offset_val, 0.8*Offset_val);}
//    fitter_ConvFit.Config().ParSettings(3).SetLimits(3., 6.);
    //fitter_ConvFit.Config().ParSettings(3).Fix();
    
    cout<<"Second fit ..."<<endl;
    ROOT::Math::Functor fcn_ConvFit2(chi2Function_ConvFit, 4);
    ROOT::Fit::Fitter fitter_ConvFit2;
    
    fitter_ConvFit2.SetFCN(fcn_ConvFit2, pStart_ConvFit_2nd);
    fitter_ConvFit2.Config().ParSettings(0).SetName("Norm");
    fitter_ConvFit2.Config().ParSettings(0).SetLimits(0.7, 5.);
    fitter_ConvFit2.Config().ParSettings(1).SetName("Gain (keV/MCA)");
    fitter_ConvFit2.Config().ParSettings(1).SetLimits(6., 12.);
    fitter_ConvFit2.Config().ParSettings(2).SetName("Offset");
    fitter_ConvFit2.Config().ParSettings(2).SetLimits(-20, 20);
    fitter_ConvFit2.Config().ParSettings(3).SetName("Resolution");
    fitter_ConvFit2.Config().ParSettings(3).SetLimits(3., 10.);
    
    ok_ConvFit = fitter_ConvFit2.FitFCN();
    if(!ok_ConvFit) {
        Error("fit", "%s", Form("%s / %s Conv Fit failed", Simu_name.Data(), Data_name.Data()));
    }
    
    const ROOT::Fit::FitResult &result_ConvFit_2nd = fitter_ConvFit2.Result();
    //result_ConvFit_2nd.Print(std::cout);
    
    Norm = result_ConvFit_2nd.Parameter(0);
    Norm_err = result_ConvFit_2nd.ParError(0);
    
    Gain_val = result_ConvFit_2nd.Parameter(1);
    Gain_err = result_ConvFit_2nd.ParError(1);
    
    Offset_val = result_ConvFit_2nd.Parameter(2);
    Offset_err = result_ConvFit_2nd.ParError(2);
    
    Resol_val = result_ConvFit_2nd.Parameter(3);
    Resol_err = result_ConvFit_2nd.ParError(3);
    
    // Reproduce values
    Cal_vect = CalibVector(Data_vect, MCA_bins, V_Data, Gain_val, Offset_val, E_bins);
    TVectorD Cal_Data_vect = Cal_vect[0];
    TVectorD Err_Cal_Data_vect = Cal_vect[0];

    //---- construction de V_data
    TMatrixD Cal_V_Data(N_E, N_E);
    for (int i=0; i<N_E; i++){
        Cal_V_Data(i,i)= Err_Cal_Data_vect(i);
    }
    
    //------ Resolution on the Simulation
    double Binwidth = E_bins(1)-E_bins(0);
    double E_j, E_i, sigma;
    for (int i =0; i<N_E ; i++){
        E_i = E_bins(i)+Binwidth/2;
        for (int j =0; j< N_E; j++){
            E_j=E_bins(j)+Binwidth/2;
            sigma=Resol_val*sqrt(E_j);
            A(i,j)=TMath::Gaus(E_i, E_j, sigma, true);
        }
    }
    
    //------ Normalize A
    double Integral_gauss=0;
    double Emin = E_bins(0);
    double Emax = E_bins(N_E-1)+Binwidth;
    TVectorD Colj_A(N_E);
    
    for (int j=0; j<N_E ; j++){
        E_j=E_bins(j)+Binwidth/2;
        sigma=Resol_val*sqrt(E_j);
        Integral_gauss=ROOT::Math::normal_cdf(Emax, sigma, E_j)-ROOT::Math::normal_cdf(Emin, sigma, E_j);
        
        for (int i=0; i<N_E; i++){
            Colj_A(i)=A(i,j);
        }
        
        for (int i=0; i<N_E; i++){
            A(i,j)*=Integral_gauss/Get_Sum_of_VectorElmts(Colj_A);
        }
    }
    
   
//    for (int j =0; j<N_E ; j++){
//        double E_j = E_bins(j);
//        for (int i=0; i<N_E; i++){
//            Colj_A(i)=A(i,j);
//        }
//        for (int i=0; i<N_E; i++){
//            cout<<Get_Sum_of_VectorElmts(Colj_A)<<endl;
//        }
//    }

//    cout<<"Integral vector before resol ="<<Get_Sum_of_VectorElmts(Simu_vect)<<endl;
    
    //Resolution du spectre de simu
    TVectorD Resol_Simu_vect = (A*Simu_vect);
    
    TMatrixD Af_T(TMatrixD::kTransposed, A);
    TMatrixD V_Resol_Simu = A*V_Simu*Af_T;
    
//    cout<<"Integral vector after resol ="<<Get_Sum_of_VectorElmts(Resol_Simu_vect)<<endl;
    
//    A.Draw("colz");
//    E_bins.Print();
//    Resol_Simu_vect.Print();
    
    //Cut vectors in fit range
    TVectorD Fit_range_Resol_Simu_vect = Cut_Vector(Resol_Simu_vect, E_bins, Emin_fit, Emax_fit);
    TMatrixD Fit_range_V_Resol_Simu = Cut_Matrix(V_Resol_Simu, E_bins, Emin_fit, Emax_fit, E_bins, Emin_fit, Emax_fit);
    TVectorD Fit_range_E_bins = Cut_Vector(E_bins, E_bins, Emin_fit, Emax_fit);
    
    TVectorD Fit_range_Cal_Data_vect = Cut_Vector(Cal_Data_vect, E_bins, Emin_fit, Emax_fit);
    TMatrixD Fit_range_Cal_V_Data = Cut_Matrix(Cal_V_Data, E_bins, Emin_fit, Emax_fit, E_bins, Emin_fit, Emax_fit);
    
    int fit_range_N_E = Fit_range_E_bins.GetNrows();
    
    //----------- Projection in a base where V_data=/=0
    int N_zeros = 0;
    for (int i =0; i<fit_range_N_E; i++){
            N_zeros+=(Fit_range_Cal_Data_vect(i)==0); //count the number of zero in Data_vect
    }

    int N_E_Proj = fit_range_N_E-N_zeros;
    TMatrixD Proj(N_E_Proj, fit_range_N_E);
    
    int i_E_Proj = 0;
    for (int j=0; j<fit_range_N_E; j++){
        if (Fit_range_Cal_Data_vect(j)!=0){
            Proj(i_E_Proj, j)=1;
            i_E_Proj+=1;
        }
    }
    
//    E_bins.Print();
//    Cal_Data_vect.Print();
//
    //---- Projection du spectre data
    TVectorD Proj_Cal_Data_vect = Proj*Fit_range_Cal_Data_vect;

    //---- projection de V_data
    TMatrixD Proj_T(TMatrixD::kTransposed, Proj);
    TMatrixD Proj_Cal_V_Data = Proj*Fit_range_Cal_V_Data*Proj_T;
    
    //---- projection simu
    TVectorD Proj_Simu_vect = Proj*Fit_range_Resol_Simu_vect;
    
    //---- projection de V_simu
    TMatrixD Proj_V_Simu = Proj*Fit_range_V_Resol_Simu*Proj_T;
    
    
    //---- projection de Energy bins
    TVectorD Proj_E_bins = Proj*Fit_range_E_bins;

    TMatrixD V_Delta(N_E_Proj, N_E_Proj);
    TVectorD Delta_Data_SimuConv(N_E_Proj);
    
    Proj_Simu_vect = Norm * Proj_Simu_vect ;
    
    Delta_Data_SimuConv = Proj_Cal_Data_vect - Proj_Simu_vect ;
    
    V_Delta = Proj_Cal_V_Data + Norm*Norm*Proj_V_Simu;

    TMatrixD V_Delta_i(TMatrixD::kInverted, V_Delta);
    
    double chi2 = Delta_Data_SimuConv * (V_Delta_i * Delta_Data_SimuConv);
    cout<<"chi2 = "<<chi2<<endl;
    TCanvas* c1 = new TCanvas();
    TGraph* gData = new TGraph(Proj_E_bins, Proj_Cal_Data_vect);
    
    TGraph* gSimu = new TGraph(Proj_E_bins, Proj_Simu_vect);
    TGraph* gDiff = new TGraph(Proj_E_bins, Delta_Data_SimuConv);
    gData->Draw("AP");
    gData->SetMarkerStyle(21);
    gData->SetMarkerColor(kBlue);
    gSimu->SetLineColor(kRed);
    gSimu->SetMarkerStyle(21);
    gSimu->SetMarkerColor(kRed);
    gSimu->Draw("PSAME");
    TCanvas* c2 = new TCanvas();
    gDiff->SetMarkerStyle(21);
    gDiff->SetMarkerColor(21);
    gDiff->Draw("AP");
    
    double Chi2_dof = result_ConvFit_2nd.MinFcnValue();
    int N_par = result_ConvFit_2nd.NPar();
    int N_dof = N_E_Proj - N_par;
    
    cout << "************************************************"<<endl;
    cout << Form("Best fit %s convolution to %s : ", Simu_name.Data(), Data_name.Data()) <<endl;
    cout << "Gain = " << Gain_val << " keV/MCA" << " +- " << Gain_err << " keV/MCA"<<endl;
    cout << "Offset = " << Offset_val << " keV" << " +- " << Offset_err << " keV"<<endl;
    cout << "Resolution = (" << Resol_val << " +- " << Resol_err << ") x sqrt(E)"<<endl;
    cout <<"Norm = " << Norm << " +- " << Norm_err <<endl;
    cout<<endl;
    cout << "Chi2/dof = " << Chi2_dof*N_dof << "/" << N_dof << " = " << Chi2_dof<<endl;
    cout << "************************************************"<<endl;
    
    return result_ConvFit_2nd;
}

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
                               h_qdc->FindBin(QDC_upperedge)+1 - (h_qdc->FindBin(QDC_loweredge)+1),
                               h_qdc->GetBinLowEdge(h_qdc->FindBin(QDC_loweredge)+1),
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

//--------------------------------------------------------
// Fit function histos
ROOT::Fit::FitResult Fit_simu_data(TH1D* hSimu, TH1D* hData, Double_t Emin_fit, Double_t Emax_fit, double Gain_ini = 9., double Norm_ini = 0.8, double Offset_ini =0., double Resol_ini = 4.){
    
    //-----Creating the vectors to do the fit
    int Nbins_MCA = hData->GetNbinsX();
    int Nbins_E = hSimu->GetNbinsX();
    TVectorD E_bins(Nbins_E) ; TVectorD Simu_vect(Nbins_E);
    TVectorD MCA_bins(Nbins_MCA) ; TVectorD Data_vect(Nbins_MCA);
    TMatrixD V_Simu(Nbins_E, Nbins_E); TMatrixD V_Data(Nbins_MCA, Nbins_MCA);
    TString Data_name = hData->GetTitle();
    TString Simu_name = hSimu->GetTitle();
    Simu_name = "convolved simu "+Simu_name;
    
    double Integral =0.;
    for (int i =0; i<Nbins_E; i++){
        E_bins(i)= hSimu->GetBinLowEdge(i+1);
        Simu_vect(i) = hSimu->GetBinContent(i+1);
        V_Simu(i,i) = hSimu->GetBinError(i+1)*hSimu->GetBinError(i+1);
//        cout<<"E = "<< E_bins(i)<<"  "<<Simu_vect(i)<<endl;
//        Integral+=hSimu->GetBinContent(i+1);
    }
//    cout<<"Integral histo in vector = "<<Integral<<endl;
    
    for (int i =0; i<Nbins_MCA; i++){
        MCA_bins(i)= hData->GetBinLowEdge(i+1);
        Data_vect(i) = hData->GetBinContent(i+1);
        V_Data(i,i) = hData->GetBinError(i+1)*hData->GetBinError(i+1);
    }
    
    TCanvas* c = new TCanvas();
    TGraph* gData = new TGraph(MCA_bins, Data_vect);
    gData->Draw();
    TCanvas* c2 = new TCanvas();
    TGraph* gSimu = new TGraph(E_bins, Simu_vect);
    gSimu->Draw();
    
    ROOT::Fit::FitResult resultfit = Do_Conv(MCA_bins, E_bins, Simu_vect, Data_vect, V_Simu, V_Data, Simu_name, Data_name, Emin_fit, Emax_fit, Gain_ini, Norm_ini, Offset_ini, Resol_ini);
    return resultfit;
}

ROOT::Fit::FitResult Fit_and_Plot(TString Name, TString FilenameData, Double_t LiveTimeData, TString DataTitle, TString FilenameSimu, Double_t FluxSimu, TString SimuTitle, Double_t NumEvSimu, double Gain_ini = 9., double Norm_ini = 0.8, double Offset_ini =0., double Resol_ini = 4., Double_t Emin_fit = 0., Double_t Emax_fit=5000., int RebinFact = 4, Bool_t CoincMu =false, Double_t MV_thr = 1.){
    
    TH1D* h_Data = Put_data_in_histo(FilenameData, LiveTimeData, DataTitle);
    h_Data->Scale(h_Data->GetBinWidth(1));// in ev
    
    TH1D* h_Simu = Put_simu_in_histo(FilenameSimu, SimuTitle, 512/RebinFact, 0., 5000., NumEvSimu, FluxSimu, MV_thr, CoincMu, 0.);
    TCanvas* c = new TCanvas();
    h_Simu->Draw();
    
//    cout<<"Integral histo before resol ="<<h_Simu_cut->Integral()<<endl;
//    TCanvas* c = new TCanvas();
//    h_SimuMuons_8inch_1MeV->Draw();
    
    ROOT::Fit::FitResult Fit = Fit_simu_data(h_Simu, h_Data, Emin_fit, Emax_fit, Gain_ini, Norm_ini, Offset_ini, Resol_ini);
    double Gain = Fit.Parameter(1);
    double Offset = Fit.Parameter(2);
    double Resol = Fit.Parameter(3);
    double Norm = Fit.Parameter(0);
    
    TH1D* hCal_Data = Calib(h_Data, Gain, Form("%s_Cal_fit", h_Data->GetTitle()), Offset, 512/RebinFact, 0., 5000.);

//    for (int i=1; i<hCal_Data->GetNbinsX(); i++){
//        if ((hCal_Data->GetBinLowEdge(i)<Emax_fit)&&(hCal_Data->GetBinLowEdge(i)>Emin_fit)){
//            cout<<hCal_Data->GetBinLowEdge(i)<<"   "<<hCal_Data->GetBinContent(i)<<endl;
//        }
//    }
    
    TH1D* h_Simu_resol = Put_simu_in_histo(FilenameSimu, SimuTitle, 512/RebinFact, 0., 5000., NumEvSimu, FluxSimu, MV_thr, CoincMu, Resol);
    h_Simu_resol->Scale(Norm);
//    cout<<"\n histo simu resol not cut"<<endl;
//    for (int i=1; i<h_Simu_resol->GetNbinsX()+1; i++){
//        if ((h_Simu_resol->GetBinLowEdge(i)<Emax_fit)&&(h_Simu_resol->GetBinLowEdge(i)>Emin_fit)){
//            cout<<h_Simu_resol->GetBinLowEdge(i)<<"   "<<h_Simu_resol->GetBinContent(i)<<endl;
//        }
//    }
//    TH1D* h_Simu_resol_cut=Cut_Histo(h_Simu_resol, Emin_fit, Emax_fit);
//    cout<<"\n histo simu resol after cut"<<endl;
//    for (int i=1; i<h_Simu_resol_cut->GetNbinsX()+1; i++){
//        cout<<h_Simu_resol_cut->GetBinLowEdge(i)<<"   "<<h_Simu_resol_cut->GetBinContent(i)<<endl;
//    }
//    cout<<"Integral histo after resol ="<<h_Simu_resol_cut->Integral()<<endl;
//    for (int i=1; i<h_Simu_resol_cut->GetNbinsX()+1; i++){
//        if ((h_Simu_resol_cut->GetBinLowEdge(i)<Emax_fit)&&(h_Simu_resol_cut->GetBinLowEdge(i)>Emin_fit)){
//            cout<<h_Simu_resol_cut->GetBinLowEdge(i)<<"   "<<h_Simu_resol_cut->GetBinContent(i)<<endl;
//        }
//    }
    
    TCanvas* C = new TCanvas("c"+Name, "c"+Name, 1200, 800);
    C->Draw();
    C->SetLogy();
    
    TPad* p1 = new TPad("phisto_"+Name, "phisto_"+Name, 0., 0.4, 0.795, 1.);
    p1->Draw();
    p1->cd();
    p1->SetRightMargin(0.03);
    p1->SetBottomMargin(0.1);
    hCal_Data->GetXaxis()->SetRangeUser(0, 4580);
    hCal_Data->Draw();
    h_Simu_resol->SetLineColor(kRed);
    h_Simu_resol->Draw("SAME");

    hCal_Data->SetTitle("");
    gStyle->SetOptStat(0);

    TLegend* leg0 = new TLegend(0.6, 0.7, 0.89, 0.89);
    leg0->AddEntry(hCal_Data, "Data (Calibrated)");
    leg0->AddEntry(h_Simu_resol, "Simulation");
    leg0->Draw("SAME");
    leg0->SetLineWidth(0);

    C->cd();
    TPad* p2 = new TPad("pRes_Muons_8inch", "pRes_Muons_8inch", 0., 0., 1., 0.4);
    p2->Draw();
    p2->cd();
    TPad *p21 = new TPad("p1","p1",0.,0.02,0.795,1.);
    p21->Draw();
    p21->cd();
    p21->SetRightMargin(0.03);
    p21->SetBottomMargin(0.2);

    std::vector<TH1D*> histosRes = Residuals_histos(hCal_Data, h_Simu_resol, Emin_fit, Emax_fit, true);
    TH1D* hRes=histosRes[0];
    hRes->GetYaxis()->SetTitle("Residuals [%]");
    hRes->GetYaxis()->SetTitleSize(1.9*hCal_Data->GetYaxis()->GetTitleSize());
    hRes->GetXaxis()->SetTitleSize(1.9*hCal_Data->GetXaxis()->GetTitleSize());
    hRes->GetYaxis()->SetLabelSize(1.9*hCal_Data->GetYaxis()->GetLabelSize());
    hRes->GetXaxis()->SetLabelSize(1.9*hCal_Data->GetXaxis()->GetLabelSize());
    hRes->GetYaxis()->SetTitleOffset(0.7);
    hRes->Draw("");

    TH1D* h_10pc = (TH1D*) hRes->Clone();
    for (int i =1; i<=h_10pc->GetNbinsX(); i++){
        h_10pc->SetBinContent(i, 0);
        h_10pc->SetBinError(i, 10);
    }
    h_10pc->SetLineColorAlpha(kOrange+1, 0.2);
    h_10pc->SetMarkerColorAlpha(kOrange+1, 0.2);
    h_10pc->SetFillColorAlpha(kOrange+1, 0.2);
    h_10pc->Draw("E3SAME");

    TLegend* leg = new TLegend(0.15, 0.3, 0.30, 0.4);
    leg->AddEntry(h_10pc, "#pm 10%");
    leg->Draw("SAME");
    leg->SetLineWidth(0);

    p2->cd();
    TPad *p22 = new TPad("p2","p2",0.795,0.02,1.,1.);
    p22->Draw();
    p22->cd();
    p22->SetLeftMargin(0.15);
    p22->SetBottomMargin(0.2);

    TH1D* hRedRes=histosRes[1];
    hRedRes->Draw("");
    hRedRes->GetXaxis()->SetTitle("Residuals [#sigma]");
    hRedRes->GetYaxis()->SetTitleSize(3.*hCal_Data->GetYaxis()->GetTitleSize());
    hRedRes->GetXaxis()->SetTitleSize(2.7*hCal_Data->GetXaxis()->GetTitleSize());
    hRedRes->GetXaxis()->SetTitleOffset(0.7);
    hRedRes->GetYaxis()->SetLabelSize(3.*hCal_Data->GetYaxis()->GetLabelSize());
    hRedRes->GetXaxis()->SetLabelSize(2.7*hCal_Data->GetXaxis()->GetLabelSize());
    hRedRes->GetXaxis()->SetNdivisions(5);
    
    return Fit;
}

TGraphErrors * Plot_Calib_withErrors(ROOT::Fit::FitResult fit_res, Double_t Emax, TString hist_name, int lincolor=2, bool first=false, double fac=1){
    
    TF1* cal_fit = new TF1("cal_fit","([0]+[1]*x)",0,20000);
    cal_fit->SetParNames("b [keV]","a [keV/MCA]");
    cal_fit->FixParameter(0,fit_res.Parameter(2));
    cal_fit->FixParameter(1,fit_res.Parameter(1));
    
    TMatrixTSym<double> Cov_tot = TMatrixTSym<double>(4);
    fit_res.GetCovarianceMatrix(Cov_tot);
    auto Cov = Cov_tot.GetSub(1, 2, 1, 2);
    Cov_tot.Print();
    
    Cov(1,1)=fit_res.ParError(2)*fit_res.ParError(2);
    Cov.Print();
    
    TH1D *histo = new TH1D(hist_name, hist_name, Emax, 0, Emax);
    for (int i=0; i<histo->GetNbinsX(); i++){
        double MCA_Val = histo->GetBinCenter(i+1);
        histo->SetBinContent(i+1, cal_fit->Eval(MCA_Val));
//        cout<<"Ereal="<<Ereal<<endl;
//        cout<<"V="<<Func_fit1->Eval(Ereal)<<endl;
//        cout<<"E1="<<Ereal-func_p_Sigma1->GetX(Func_fit1->Eval(Ereal), Ereal-500, Ereal+500)<<endl;
    }
    
    TGraphErrors *histo_err = new TGraphErrors(histo);
    for (int i=0; i<histo->GetNbinsX(); i++){
        double MCA = histo->GetBinCenter(i+1);
        std::cout<<sigma1(MCA, Cov)<<std::endl;
        histo_err->SetPointError(i+1, fac*sigma1(MCA, Cov));
    }
    
    cal_fit->SetLineWidth(3);
    cal_fit->SetLineColor(lincolor);
    cal_fit->SetLineStyle(9);
    if (first) {
        histo_err->SetMarkerColor(0);
        histo_err->SetLineColor(0);
        histo_err->SetFillColorAlpha(lincolor,0.4);
        histo_err->Draw("APE4 same");
        histo_err->GetXaxis()->SetRangeUser(0, 200);
        histo_err->GetXaxis()->SetTitle("Pulse Height [MCA unit]");
        histo_err->GetYaxis()->SetTitle("Calibrated Energy [keVee]");
        histo_err->SetTitle("#times#sigma");
        
        cal_fit->Draw("L same");
    }
        
    else {
        histo_err->SetMarkerColor(0);
        histo_err->SetLineColor(0);
        histo_err->SetFillColorAlpha(lincolor,0.4);
        histo_err->Draw("PE4 same");
        histo_err->GetXaxis()->SetRangeUser(0, 200);
        histo_err->GetXaxis()->SetTitle("Pulse Height [MCA unit]");
        histo_err->GetYaxis()->SetTitle("Calibrated Energy [keVee]");
        histo_err->SetTitle("#times#sigma");
        
        cal_fit->Draw("L same");
    }
    
    return histo_err;
}

int BS_Calib(){
    // Muons
    ROOT::Fit::FitResult Fit_mu = Fit_and_Plot("Muons","Measures@Saclay/2022-11-14_BS-A-06_RT_489309.69_LV_182.41_MuonCoinc.dat", 489309.69, "8inch at Saclay muons in coincidence", "../BS_BkgModel/Atm_Muons/BS_D08/Analyzer_T_run01-12.root", 0.019, "8 inches atmospheric muons", 9.60e8, 9., 0.8, 0., 4., 600., 4000.,  4, true, 1.);

    //Co60
    ROOT::Fit::FitResult Fit_Co = Fit_and_Plot("Co60","Measures@Saclay/2022-11-08_BS-A-06_RT_2114.05_LV_2112.16_Co60.dat", 2114.05, "Co60 source", "Simulations/Co60/Analyzer_T_run01-20.root", 10240./(40*40), "Simu Co60 8 inch", 3e8, 9., 0.8, 0., 4., 600., 1200., 2, false);

    //Ba133
    ROOT::Fit::FitResult Fit_Ba = Fit_and_Plot("Ba133","Measures@Saclay/2022-11-08_noBS_RT_302.42_LV_301.89_Ba133.dat", 302.42, "Ba133 source", "Simulations/Ba133/Analyzer_PT_LIEu_Ba133_2.root", 0.6*290000./(40*40),"Ba133 simulation", 10e7, 9., 1.6, 0., 4., 250., 550., 1, false);

    TVectorD GainVal(3), GainErr(3);
    TVectorD OffsetVal(3), OffsetErr(3);
    TVectorD EVal(3), EErr(3);
    TVectorD Resol_ov_EVal(3), Resol_ov_EErr(3);

    //---
    GainVal(0) = Fit_mu.Parameter(1);
    GainErr(0) = Fit_mu.ParError(1);

    GainVal(1) = Fit_Co.Parameter(1);
    GainErr(1) = Fit_Co.ParError(1);

    GainVal(2) = Fit_Ba.Parameter(1);
    GainErr(2) = Fit_Ba.ParError(1);
    
    //---
    OffsetVal(0) = Fit_mu.Parameter(2);
    OffsetErr(0) = Fit_mu.ParError(2);

    OffsetVal(1) = Fit_Co.Parameter(2);
    OffsetErr(1) = Fit_Co.ParError(2);

    OffsetVal(2) = Fit_Ba.Parameter(2);
    OffsetErr(2) = Fit_Ba.ParError(2);

    //---
    EVal(0) = 1000.;
    EErr(0) = 200.;

    EVal(1) = 900.;
    EErr(1) = 200.;

    EVal(2) = 300.;
    EErr(2) = 100.;

    //---
    Resol_ov_EVal(0) = Fit_mu.Parameter(3);
    Resol_ov_EErr(0) = Fit_mu.ParError(2);
    
    Resol_ov_EVal(1) = Fit_Co.Parameter(3);
    Resol_ov_EErr(1) = Fit_Co.ParError(3);
    std::cout<<"error  on Co resol = "<<Resol_ov_EErr(1)<<std::endl;
    
    Resol_ov_EVal(2) = Fit_Ba.Parameter(3);
    Resol_ov_EErr(2) = Fit_Ba.ParError(3);
    
    //-----
    TGraphErrors* gGain = new TGraphErrors(EVal, GainVal, EErr, GainErr);
    gGain->SetMarkerStyle(21);
    gGain->SetMarkerColor(kRed);
    gGain->GetXaxis()->SetTitle("Energy range [keVee]");
    gGain->GetYaxis()->SetTitle("Gain [keVee/MCA]");

    TCanvas* cgain = new TCanvas("Gains", "Gains");
    gGain->Draw("AP");

    TGraphErrors* gResol = new TGraphErrors(EVal, Resol_ov_EVal, EErr, Resol_ov_EErr);
    gResol->SetMarkerStyle(21);
    gResol->SetMarkerColor(kRed);
    gResol->GetYaxis()->SetTitle("Gain [keVee/MCA]");
    gResol->GetYaxis()->SetTitle("#sigma/E []");

    TCanvas* cresol = new TCanvas("Resol", "Resol");
    gResol->Draw("AP");
//
//    //Trying to plot obtained calibration curves with the errors of the fits:
    TCanvas* cCalibs = new TCanvas("Calibs with errors", "Calibs with errors");
    
    auto Ba_err = Plot_Calib_withErrors(Fit_Ba, 200, "Ba133", 2, true);
    auto Co_err = Plot_Calib_withErrors(Fit_Co, 200, "Co60", 3);
    auto Mu_err = Plot_Calib_withErrors(Fit_mu, 200, "Muons", 4, false, 0.01);
    
    TLegend* leg0 = new TLegend(0.6, 0.7, 0.89, 0.89);
    leg0->AddEntry(Mu_err, "Muons");
    leg0->AddEntry(Co_err, "^{60}Co");
    leg0->AddEntry(Ba_err, "^{133}Ba");
    leg0->Draw("SAME");
    leg0->SetLineWidth(0);
    
    return(1);
}
