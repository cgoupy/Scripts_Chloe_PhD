#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>
#include <TPRegexp.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>

using namespace TMath;

int REBIN_VAL = 1; //or 3

//Declare functions that will be efined in the code after
TSpectrum* FindPeaks(TH1D* histo, TString plot_option);
std::vector<std::vector<Double_t>> FitPeaks(TH1D* histo, TString plot_option, bool Print=true, Double_t ZoomAt=0);
void FitResol(TGraphErrors* graph, int color);
void FitResol_Const(TGraphErrors* graph, int color);
double PeakModel(double *x, double *par);
double ResolutionModel(double *x, double *par);

void PlotSpectra(Double_t ZoomAt = 0){
    
    Double_t M = 1;
    
    std::cout << "HPGe mass: " << M << " kg" << std::endl;
    
    //Files to plot
    const Int_t NFiles = 1; //number of files to plot
    TString FileName[NFiles] = {"/Users/cg264002/research/HighStatSimulations/UGL_Simulations/NoIS/Data-Simu_COV/MonPlotsCOV_019_014.root"}; //path and name of the files
    Int_t Color[NFiles] = {1}; //colors for the plot
    
    TString hLegend[NFiles] = {"COV run 19"};
    
    TFile *RunFile = nullptr; //Declare file ROOT class
    TH1D* h[NFiles]; //Declare histogram ROOT class
    Double_t Runtime[NFiles]; //Declare histogram ROOT class
    
    for (Int_t i = 0; i< NFiles; i++){ //Loop over files
        RunFile = new TFile(FileName[i].Data(),"READ"); //read the file i
        h[i] = (TH1D*) RunFile->Get("cal_spectrum"); //Get the correct histogram in the root file
        h[i]->SetDirectory(0); //plot esthetic
        h[i]->SetStats(0); //plot esthetic
        h[i]->SetTitle(""); //remove the title (for esthetic)
        h[i]->SetLineColor(Color[i]);
        h[i]->SetLineWidth(2);
        h[i]->Rebin(REBIN_VAL);
        Double_t binW = h[i]->GetBinWidth(0);
        
        //Normalization in kg-1 d-1 keV-1
        h[i]->Scale(86400/(binW*M));
        h[i]->GetYaxis()->SetTitle("Differential count rate [kg^{-1} d^{-1} keV^{-1}]");
        h[i]->GetYaxis()->SetRangeUser(1,1e6);
        
        RunFile->Close();
    }
    
    //Draw spectra
    TCanvas *CanHisto = new TCanvas("CanHisto","Background spectra",0,0,1000,700);
    CanHisto->cd();
    CanHisto->SetRightMargin(0.02);
    CanHisto->SetTopMargin(0.03);
    
    Double_t ICR[NFiles]; //Integrated counting rate
    Double_t ErrICR[NFiles]; //Error on integrated counting rate
    
    TLegend* leg = new TLegend(0.6, 0.8, 0.97, 0.96); //position of the legend
    std::vector<std::vector<Double_t>> FittedPeaks; //result of the FitPeaks function
    
    for (Int_t i = 0; i< NFiles; i++){
        ICR[i] = h[i]->IntegralAndError(h[i]->FindBin(0),h[i]->FindBin(2700),ErrICR[i],"width");
        std::cout <<"\n"<< FileName[i] << ", integral count rate: " << ICR[i]/(24*60) << " +/- " << ErrICR[i]/(24*60) << " count/(min.kg)" << std::endl;
        h[i]->SetName(hLegend[i].Data()); //Set the correct Name to the histograms
        leg->AddEntry(h[i], hLegend[i].Data()); //Set the legend
        
        CanHisto->cd();
        if (i==0) FittedPeaks = FitPeaks(h[i],"", true, ZoomAt);
        else FittedPeaks = FitPeaks(h[i],"nodraw", true, ZoomAt); //if it is the first histogram, show the peak finder markers, if not, don't show them (just to improve the read and not have red triangles everywhere)
        
        CanHisto->cd();
        h[i]->Draw("E0 HIST SAME");
    }
    
    gPad->SetLogy(); //Log scale for y
    gPad->SetGridx(); //plot grid
    gPad->SetGridy(); //plot grid
    leg->Draw("SAME"); //Draw the legend
}

void PlotEnergyResolution(){
    
    //Mass of the HPGe detector
    Double_t M =1;
    
    std::cout << "HPGe mass: " << M << " kg" << std::endl;
    
    //Files to plot
    const Int_t NFiles = 1; //number of files to plot
    TString FileName[NFiles] = {"/Users/cg264002/research/HighStatSimulations/UGL_Simulations/NoIS/Data-Simu_COV/MonPlotsCOV_019_014.root"}; //path and name of the files
    Int_t Color[NFiles] = {1}; //colors for the plot
    
    TString hLegend[NFiles] = {"COV run 19"};
    
    TFile *RunFile = nullptr; //Declare file ROOT class
    TH1D* h[NFiles]; //Declare histogram ROOT class
    TGraphErrors* Resol[NFiles]; //Declare TGraphErrors ROOT class
    Double_t Runtime[NFiles]; //Declare Runtimes
    
    for (Int_t i = 0; i< NFiles; i++){ //Loop over files
        RunFile = new TFile(FileName[i].Data(),"READ"); //read the file i
        h[i] = (TH1D*) RunFile->Get("cal_spectrum"); //Get the correct histogram in the root file
        h[i]->SetDirectory(0); //plot esthetic
        h[i]->SetStats(0); //plot esthetic
        h[i]->SetTitle(""); //remove the title (for esthetic)
        h[i]->SetLineColor(Color[i]);
        h[i]->SetLineWidth(2);
        h[i]->Rebin(REBIN_VAL);
        Double_t binW = h[i]->GetBinWidth(0);
        
        //Normalization in kg-1 d-1 keV-1
        h[i]->Scale(86400/(binW*M));
        h[i]->GetYaxis()->SetTitle("Differential count rate [kg^{-1} d^{-1} keV^{-1}]");
        h[i]->GetYaxis()->SetRangeUser(1,1e6);
        
        RunFile->Close();
    }
    
    //Draw spectra
    TCanvas *CanGraphResol = new TCanvas("CanGraphResol","Detector resolution",0,0,1000,700);
    CanGraphResol->cd();
    
    TLegend* leg = new TLegend(0.6, 0.8, 0.97, 0.94); //position of the legend
    std::vector<std::vector<Double_t>> FittedPeaks; //result of the FitPeaks function
    
    for (Int_t i = 0; i< NFiles; i++){
        h[i]->SetName(hLegend[i].Data()); //Set the correct Name to the histograms
        
        if (i==0) FittedPeaks = FitPeaks(h[i],"nodraw", false);
        else FittedPeaks = FitPeaks(h[i],"nodraw", false); //if it is the first histogram, show the peak finder markers, if not, don't show them (just to improve the read and not have red triangles everywhere)
        
        int NFittedPeaks = FittedPeaks.size();
        TVectorD FittedPeakPos(NFittedPeaks);
        TVectorD ErrFittedPeakPos(NFittedPeaks);
        TVectorD FittedPeakSigma(NFittedPeaks);
        TVectorD ErrFittedPeakSigma(NFittedPeaks);
        TVectorD FittedPeakResol(NFittedPeaks); // the resolution is defined as sigma/E (in %)
        TVectorD ErrFittedPeakResol(NFittedPeaks);
        
        for (int k =0; k<NFittedPeaks; k++){ //loop over peak to take the results
            FittedPeakPos(k) = FittedPeaks[k][0];
            ErrFittedPeakPos(k) = FittedPeaks[k][1];
            FittedPeakSigma(k) = FittedPeaks[k][2];
            ErrFittedPeakSigma(k) = FittedPeaks[k][3];
            
            FittedPeakResol(k) = FittedPeakSigma(k)/FittedPeakPos(k)*100;
            ErrFittedPeakResol(k) = sqrt((ErrFittedPeakSigma(k)*ErrFittedPeakSigma(k))/(FittedPeakSigma(k)*FittedPeakSigma(k))
                                         + (ErrFittedPeakPos(k)*ErrFittedPeakPos(k))/(FittedPeakPos(k)*FittedPeakPos(k)))
                                    *FittedPeakResol(k); // Error propagation formula
        }
        Resol[i] = new TGraphErrors(FittedPeakPos, FittedPeakResol, ErrFittedPeakPos, ErrFittedPeakResol);
        Resol[i]->SetTitle("");
        Resol[i]->SetMarkerColor(Color[i]);
        Resol[i]->SetMarkerStyle(21);
        
        Resol[i]->GetXaxis()->SetTitle("Energy [keV]");
        Resol[i]->GetYaxis()->SetTitle("Resolution (#sigma/E) [%]");
        
        CanGraphResol->cd();
        gPad->SetRightMargin(0.02);
        gPad->SetTopMargin(0.03);
        gPad->SetLeftMargin(0.1);
        
        if (i==0) Resol[i]->Draw("AP");
        else Resol[i]->Draw("PSAME");
        
        FitResol(Resol[i], Color[i]);
        //FitResol_Const(Resol[i], Color[i]);
        leg->AddEntry(h[i], hLegend[i].Data()); //Set the legend
    }
    
    leg->Draw("SAME"); //Draw the legend
}


// Find peak with ROOT pFinder
TSpectrum* FindPeaks(TH1D* histo, TString plot_option){
    
    int MaxPeakNumber = 100; // maximal number of peak to find in the histogram
    TSpectrum* pFinder = new TSpectrum(MaxPeakNumber); // TSpectrum class from ROOT: https://root.cern.ch/doc/master/classTSpectrum.html
    pFinder->Search(histo, 1.05, plot_option, 0.005); //Search peak with a threshold equal to 0.005*maximum, and a sigma of 0.05, to be honest, those parameters are set arbitrary with the found peaks, you can play a bit with them
    
    return pFinder; //return the TSpectrum to do future operation on it: see FitPeaks function
}

//-------- Define peak model
double PeakModel(double *x, double *par){ //Peak model
    Double_t Norm = par[0];
    Double_t Gauss_pos = par[1];
    Double_t Gauss_sigma = par[2];
    
    Double_t Bckg_a = par[3];
    Double_t Bckg_b = par[4];
    
    Double_t x_val = x[0];
    return Norm*TMath::Gaus(x_val, Gauss_pos, Gauss_sigma, kTRUE)+Bckg_b+Bckg_a*x_val; //sum of a normalized gaussian (so the Norm parameter is the integral of the gaussian) and a straight background of equation b+a*x
}

//Fit peak found with PeakModel
std::vector<std::vector<Double_t>> FitPeaks(TH1D* histo, TString plot_option, bool Print=true, Double_t ZoomAt=0){
    
    TSpectrum* pFinder = FindPeaks(histo, plot_option);

    int NPeaks = pFinder->GetNPeaks();
    cout<<NPeaks<<" peaks found in "<<histo->GetName()<<endl;
    
    double *PeakPos_initial = pFinder->GetPositionX();
    double *PeakAmp_initial = pFinder->GetPositionY();
    
    /* Sort out the peaks by increasing positions and get the maximum amplitude*/
    for(int p = 0; p < NPeaks; p++){
        for(int l = 0; l < NPeaks-1; l++){
            if(PeakPos_initial[l] > PeakPos_initial[l+1]){
                double pos_upper = PeakPos_initial[l]; double pos_lower = PeakPos_initial[l+1];
                double amp_upper = PeakAmp_initial[l]; double amp_lower = PeakAmp_initial[l+1];
                PeakPos_initial[l] = pos_lower; PeakPos_initial[l+1] = pos_upper;
                PeakAmp_initial[l] = amp_lower; PeakAmp_initial[l+1] = amp_upper;
            }
        }
    }
    
    //Create the canvas to plot all the peak fits
    TCanvas *CanPeakFit = new TCanvas(Form("PeakFitted_%s",histo->GetName()),Form("PeakFitted in %s",histo->GetName()),0,0,1500,1000);
    CanPeakFit->Divide(sqrt(NPeaks)-1, sqrt(NPeaks)-1); //divide the canvas, reduced by one column arbitrarly
    gStyle->SetOptTitle(1);
    
    //Prepare the vector to return from the function
    std::vector<std::vector<Double_t>> FittedPeaks; //contains: {peak position, peak position error, peak width, peak width error, peak amp, peak amp error} for each peak
    
    //Prepare the table to print in the terminal
    if (Print){
        cout<<"===================================================================================================================================================="<<endl;
        cout<<"|  Fitted energy (keV)\t|  Error (keV)\t|  Fitted resolution (keV)\t|  Error (keV)\t|  Rate in peak (counts/day/kg)\t|  Error (counts/day/kg)\t|  Chi2/Dof\t\t|"<<endl;
        cout<<"============================================================================================================================================================"<<endl;
        cout<<"|                                                        "<<histo->GetName()<<"                                                                                         \t|"<<endl;
        cout<<"==========================================================================================================================================================================="<<endl;
    }

    int FittedPeak = 1; //Number of correctly fitted peak (for the plot)
    //loop over peaks found
    for (int i=0; i<NPeaks; i++){
        
        //plot
        CanPeakFit->cd(FittedPeak);
        
        //Clone the histogram (for zooming in the plot)
        TH1D* histo_clone = (TH1D*) histo->Clone(Form("%s peak_@%0.f_keV", histo->GetName(), PeakPos_initial[i]));
        
        //Define the fit model
        
        TF1* f_Gauss = new TF1("Gauss", "gaus(0)", PeakPos_initial[i]-20,PeakPos_initial[i]+20); // fonction defined between 0 and 3600, following the user defined function "PeakModel", with 5 parameters
        f_Gauss->SetParameters(PeakAmp_initial[i], PeakPos_initial[i], 20); //set the initial values of the parameters
        f_Gauss->SetParLimits(1, 0.9*PeakPos_initial[i], 1.1*PeakPos_initial[i]); //set the initial values of the parameters
        f_Gauss->SetParLimits(2, 0.1/100.*PeakPos_initial[i], 20/100.*PeakPos_initial[i]); //set the initial values of the parameters
        histo_clone->Fit(f_Gauss, "QNR");
        
        TF1* f_lin = new TF1("lin", "pol1(0)", PeakPos_initial[i]-50, PeakPos_initial[i]+50); // fonction defined between 0 and 3600, following the user defined function "PeakModel", with 5 parameters
        f_lin->SetParameters(20, -0.001); //set the initial values of the parameters
        histo_clone->Fit(f_lin, "QNR");
        
        
        TF1* f_PeakModel = new TF1("f_PeakModel", PeakModel, PeakPos_initial[i]-4*f_Gauss->GetParameter(2), PeakPos_initial[i]+4*f_Gauss->GetParameter(2), 5); // fonction defined between 0 and 3600, following the user defined function "PeakModel", with 5 parameters
        f_PeakModel->SetLineColor(2);
        f_PeakModel->SetParNames("Norm", "Gauss_pos", "Gauss_sigma", "Bckg_a", "Bckg_b"); //set the name of the parameters
        f_PeakModel->SetParameters(f_Gauss->GetParameter(0), f_Gauss->GetParameter(1), f_Gauss->GetParameter(2), f_lin->GetParameter(1), f_lin->GetParameter(0)); //set the initial values of the parameters
        
        f_PeakModel->FixParameter(3, f_lin->GetParameter(1));
        f_PeakModel->FixParameter(4, f_lin->GetParameter(0));
        //f_PeakModel->FixParameter(2, f_Gauss->GetParameter(2));
        f_PeakModel->SetParLimits(2, 0.8*f_Gauss->GetParameter(2), 1.1*f_Gauss->GetParameter(2));
        
//        Fit twice the data (we fit a second time because the first fit is setting more correctly the initial parameters and the secon fit may be more reliable)
        histo_clone->Fit(f_PeakModel, "RQN", "", PeakPos_initial[i]-40, PeakPos_initial[i]+40);
                        // ^Model     ^Q = quite mode, N= not draw, ^low limit of the fit               ^high limit of the fit (defined from the initial guess peak position and 10 times an assumed peak width of 0.1% of the peak energy (i.e. peak position)
                        // Root documentation about how to fit an histogram >> https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
        TFitResultPtr result = histo_clone->Fit(f_PeakModel, "NQS", "", PeakPos_initial[i]-0.01*PeakPos_initial[i], PeakPos_initial[i]+0.01*PeakPos_initial[i]); //save the result in a pointer to a TFitResult (parameter S)

        //Read the obtained parameters an their errors from the fit
        Double_t Norm_fitted = result->Parameter(0);
        Double_t Gauss_pos_fitted = result->Parameter(1);
        Double_t Gauss_sigma_fitted = result->Parameter(2);
        Double_t Bckg_a_fitted = result->Parameter(3);
        Double_t Bckg_b_fitted = result->Parameter(4);

        Double_t Err_Norm_fitted = result->ParError(0);
        Double_t Err_Gauss_pos_fitted = result->ParError(1);
        Double_t Err_Gauss_sigma_fitted = result->ParError(2);
        Double_t Err_Bckg_a_fitted = result->ParError(3);
        Double_t Err_Bckg_b_fitted = result->ParError(4);

        //Calculate the background value in the peak position
        Double_t Mean_bckg = Bckg_a_fitted*Gauss_pos_fitted + Bckg_b_fitted;

        //prepare the plot in the fitting range
        histo_clone->GetXaxis()->SetRangeUser(PeakPos_initial[i]-60, PeakPos_initial[i]+60);
        histo_clone->GetYaxis()->SetRangeUser(Mean_bckg-0.1*Norm_fitted, Mean_bckg+.1*Norm_fitted);
        histo_clone->GetXaxis()->SetNdivisions(0);
        histo_clone->GetYaxis()->SetLabelSize(0.1);
        histo_clone->GetYaxis()->SetTitleSize(0.1);

        //Add a title
        TLatex *Latex_name = new TLatex(PeakPos_initial[i]-50, Mean_bckg-0.05*Norm_fitted, Form("#splitline{%s}{peak_@%0.f_keV}", histo->GetName(), Gauss_pos_fitted));
        Latex_name->SetTextSize(0.1);

        Double_t Chi2 = result->Chi2(); //Get the Chi2
        Double_t Dof = result->Ndf(); //Get the number of degrees of freedom

        //Plot esthetic
        gPad->SetTopMargin(0.05);
        gPad->SetRightMargin(0.);
        gPad->SetLeftMargin(0.3);
        
//        Write some results in the terminal
        if ((Chi2/Dof>0)&(PeakPos_initial[i]>500)&(PeakPos_initial[i]-Gauss_pos_fitted<50)){ //check if the fit converged (check if the chi2/ndf is correct: cut at 1.5, can be changed

            FittedPeak++;
            
            if (Print){
                histo_clone->Draw();
                f_PeakModel->Draw("SAME");
                
                f_Gauss->SetLineColor(kBlue);
                f_Gauss->Draw("SAME");
                
                f_lin->SetLineColor(kGreen);
                f_lin->Draw("SAME");
                Latex_name->Draw("SAME");
                cout<<"|  "<< Gauss_pos_fitted<<"\t\t|  "<<Gauss_pos_fitted<<"\t|  "<<Gauss_sigma_fitted<<"\t\t\t| "<<Err_Gauss_sigma_fitted<<"\t|  "<<Norm_fitted<<"\t\t\t|  "<< Err_Norm_fitted<<"\t\t|  "<< Chi2/Dof << "\t\t|"<< endl;
            }
            else CanPeakFit->Close();
            
            if ((ZoomAt < Gauss_pos_fitted+0.01*Gauss_pos_fitted) && (ZoomAt>Gauss_pos_fitted-0.01*Gauss_pos_fitted)){
                cout<< "zoom"<<endl;
                
                TCanvas *ZoomPeakFit = new TCanvas(Form("%s_PeakAt_%0.f",histo->GetName(), Gauss_pos_fitted),Form("%s_PeakAt_%0.f",histo->GetName(), Gauss_pos_fitted));
                
                ZoomPeakFit->cd();
            
                histo_clone->GetXaxis()->SetNdivisions(510);
                histo_clone->GetYaxis()->SetLabelSize(0.05);
                histo_clone->GetYaxis()->SetTitleSize(0.05);
                
                gPad->SetTopMargin(0.05);
                gPad->SetRightMargin(0.05);
                gPad->SetLeftMargin(0.18);
//
                histo_clone->Draw();
                f_PeakModel->Draw("SAME");
                Latex_name->Draw("SAME");
            }
                
            FittedPeaks.push_back({Gauss_pos_fitted, Err_Gauss_sigma_fitted, Gauss_sigma_fitted, Err_Gauss_sigma_fitted, Norm_fitted, Err_Norm_fitted});
        }
    }
    return FittedPeaks;
}

//-------- Resolution Model
double ResolutionModel(double *x, double *par){ //Peak model
    Double_t A = par[0];
    Double_t B = par[1];
    Double_t C = par[2];
    
    Double_t x_val = x[0];
    return sqrt(A/(x_val) + B/(x_val*x_val) + C); //from Knoll HPGe resolution
}

//Fit resolution with Resolution Model
void FitResol(TGraphErrors* graph, int color){
    
    //Define the fit model
    TF1* f_ResolModel = new TF1("f_ResolModel", ResolutionModel, 100, 3600, 3); // fonction defined between 100 and 3600 keV, following the user defined function "ResolutionModel", with 3 parameters
    f_ResolModel->SetLineColor(color);
    f_ResolModel->SetParNames("A", "B", "C"); //set the name of the parameters
    f_ResolModel->SetParameters(2, 1200, 7e-5); //set the initial values of the parameters
//    f_ResolModel->SetParLimits(1, 0, 1);
//    f_ResolModel->SetParLimits(2, 0, 1);
//    f_ResolModel->SetParLimits(3, 0, 1);
    
    TFitResultPtr result = graph->Fit(f_ResolModel, "SQR");
            // ^Model
    
    //Read the obtained parameters an their errors from the fit
    Double_t A = result->Parameter(0);
    Double_t B = result->Parameter(1);
    Double_t C = result->Parameter(2);

    Double_t Err_A = result->ParError(0);
    Double_t Err_B = result->ParError(1);
    Double_t Err_C = result->ParError(2);
    
    cout<<"---- fit function: sqrt(A/E + B/E^2 + C)"<<endl;
    cout<<" A = ("<<A<< " +- " << Err_A << ") keV" << endl;
    cout<<" B = ("<<B<< " +- " << Err_B << ") keV^2" << endl;
    cout<<" C = ("<<C<< " +- " << Err_C << ") \n" << endl;
    f_ResolModel->Draw("SAME");
}

//Fit resolution with Resolution Model
void FitResol_Const(TGraphErrors* graph, int color){
    
    //Define the fit model
    TF1* f_ResolModel = new TF1("f_ResolModel_const", "pol0(0)", 100, 3600); // fonction defined between 100 and 3600 keV, following the user defined function "ResolutionModel", with 3 parameters
    f_ResolModel->SetLineColor(color);
    f_ResolModel->SetParNames("A"); //set the name of the parameters
    f_ResolModel->SetParameters(2); //set the initial values of the parameters
//    f_ResolModel->SetParLimits(1, 0, 1);
//    f_ResolModel->SetParLimits(2, 0, 1);
//    f_ResolModel->SetParLimits(3, 0, 1);
    
    TFitResultPtr result = graph->Fit(f_ResolModel, "SQR");
            // ^Model
    
    //Read the obtained parameters an their errors from the fit
    Double_t A = result->Parameter(0);

    Double_t Err_A = result->ParError(0);
    
    cout<<"---- fit function: s/E = A "<<endl;
    cout<<" A = ("<<A<< " +- " << Err_A << ") %" << endl;
    f_ResolModel->Draw("SAME");
}

