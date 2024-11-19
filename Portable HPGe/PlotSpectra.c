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

void PlotSpectra(){

    //Mass of the HPGe detector
    Double_t Volume =(5.5*TMath::Pi()*pow(5.55/2,2)-4.4*TMath::Pi()*pow(1.15/2,2));
    Double_t M = Volume*5.323/1000;
    
    std::cout << "HPGe volume: " << Volume << " cm3" << std::endl;
    std::cout << "HPGe mass: " << M << " kg" << std::endl;
    
    const Int_t NFiles = 2;
    TString FileName[NFiles] = {"Data/LaboF/Background_221020.root",
        "Data/Chooz/Background-VNS-291020.root"
        //"Data/Chooz/Background-S11-291020.root"
    };
    Int_t Color[NFiles] = {1,2};
    
    //"CEA, labo F"
    TString hLegend[NFiles] = {"Labo F - CEA", "Chooz, S05 (VNS) - 1^{st} position"};//, "Chooz, S11 (comp. room)"};
    
    TFile *RunFile = nullptr;
    
    TH1D* h[NFiles];
    
    for (Int_t i = 0; i< NFiles; i++){
        RunFile = new TFile(FileName[i].Data(),"READ");
        h[i] = (TH1D*)gDirectory->Get("Cal_spectrum");
        h[i]->SetDirectory(0);
        h[i]->SetStats(0);
        h[i]->SetTitle("");
        h[i]->SetLineColor(Color[i]);
        h[i]->SetLineWidth(2);
        Double_t binW = h[i]->GetBinWidth(0);
        
        //cout << binW << endl;
        //Normalization in kg-1 d-1 keV-1
        h[i]->Scale(86400/(binW*M));
        h[i]->GetYaxis()->SetTitle("Differential count rate [kg^{-1} d^{-1} keV^{-1}]");
        h[i]->GetYaxis()->SetRangeUser(1,1e6);
        
        RunFile->Close();
    }
    
    //Draw spectra
    TCanvas *CanHisto = new TCanvas("CanHisto","Background spectra",0,0,1000,700);
    CanHisto->cd();
    
    Double_t ICR[NFiles];
    Double_t ErrICR[NFiles];
    
    for (Int_t i = 0; i< NFiles; i++){
        //ICR[i] = h[i]->Integral(h[i]->FindBin(40),h[i]->FindBin(2700),"width")/3600;
        ICR[i] = h[i]->IntegralAndError(h[i]->FindBin(0),h[i]->FindBin(2700),ErrICR[i],"width");
        std::cout << FileName[i] << ", integral count rate: " << ICR[i]/(24*60) << " +/- " << ErrICR[i]/(24*60) << " count/(min.kg)" << std::endl;
        //[i]->SetTitle(Form(hLegend[i]+": %g +/- %g min^{-1} kg^{-1}",ICR[i]/(24*60),ErrICR[i]/(24*60)));
        h[i]->SetTitle(hLegend[i].Data());
        h[i]->Draw("E0 HIST SAME");
    }
    
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    
    //    TH1D* hh = (TH1D*) h->Clone();
    //    //hh->SetName("Gamma spectrum");
    //    hh->SetTitle("Eu^{152}");
    //    //hh->SetTitle("Co^{60}");
    //    hh->SetStats(0);
    //
    
    //
}
