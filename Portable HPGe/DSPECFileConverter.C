#include <Rtypes.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TObjString.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void DSPECFileConverter(TString InputFileName, Double_t Binwidth_keV=1.)
{
    // test input file
    struct stat buf; // to test input file
    if (stat(InputFileName,&buf) == -1) { // test whether the file exist or not, and get its size
        cerr << "***Error: file missing: " << InputFileName << "\nExiting...\n";
        return;
    }
    
    // create input file
    ifstream InputFile;
    InputFile.open(InputFileName.Data());
    if (!InputFile.is_open() || !InputFile.good()) {
        cerr << "***Error: can not read file: " << InputFileName << "\nExiting...\n";
        return;
    }
    cout << "Your input file is " << InputFileName << endl;
    
    // create output file
    TString OutputFileName = InputFileName;
    TPMERegexp("\\.Spe").Substitute(OutputFileName, Form("_Bin%.1fkeV.root", Binwidth_keV));
    TPMERegexp("-ASCII\\.Spe").Substitute(OutputFileName,"");
    // Create the vectors
    TH1::SetDefaultSumw2();
    vector<Double_t> vValue;
    
    // Get the line from the input file
    Bool_t end_header=false, compteur=false, end_file=false;
    Double_t CalOrigin=-666, CalSlope=-666;
    Int_t RealTime=-1, LiveTime=-1;
    UInt_t line_number=0;
    TString str;
    
    TPMERegexp regSpace("\\s+");
    TPMERegexp regBegin("\\$DATA:");
    TPMERegexp regTime("\\$MEAS_TIM:"); // beginning of histo
    TPMERegexp regCal("\\$ENER_FIT:");
    TPMERegexp regEnd("$ROI:");
    
    while (str.ReadLine(InputFile).good()) { // get a line from the file
        if ( InputFile.fail() ) { break; }
        
        if (end_header == false) {
            if ( regBegin.Match(str) ) { end_header = true; }
            if (compteur==true) {
                if ( regSpace.Split(str)==2) {
                    RealTime = regSpace[1].Atoi();
                    LiveTime = regSpace[0].Atoi();
                }
                compteur=false;
            }
            if ( regTime.Match(str) ) { compteur = true;}
            continue;
        }
        
        if (line_number == 0) {
            if ( regSpace.Split(str) == 2 ) {
                line_number = regSpace[1].Atoi();
                vValue.reserve(line_number+1);
            } else {
                cerr << "***Error: can not parse line after header: "<<str<<endl;
                return;
            }
            continue;
        }
        
        if (end_file == false) {
            if (regEnd.Match(str)) {
                end_file = true;
                continue;
            }
            regSpace.Substitute(str,"");
            vValue.push_back(str.Atof());
        }
        
        if (regCal.Match(str)) {
            str.ReadLine(InputFile);
            if ( regSpace.Split(str)==2) {
                CalOrigin = regSpace[0].Atof();
                CalSlope =  1.00*regSpace[1].Atof();
            }
        }
        
    }
    InputFile.close();
    
    cout <<"RealTime: " <<RealTime << " s, LiveTime: "<<LiveTime<<" s\n";
    if(RealTime<0) { cerr << "***Error: real time not found in file: "<<InputFileName<<endl; return; }
    
    
    TFile* OutputFile = new TFile(OutputFileName,"RECREATE");
    if (not OutputFile->IsOpen()) { cerr << "***Error: can not create ROOT file: " << OutputFileName << endl; return; }
    cout << "Your output ROOT file is " << OutputFile->GetName() << endl;
    
    TObjString StrTime = Form("RealTime: %d s",RealTime);
    StrTime.Write("RealTime");
    StrTime = Form("LiveTime: %d s",LiveTime);
    StrTime.Write("LiveTime");
    
    //   TPMERegexp("\\.root").Substitute(OutputFileName,"_time.txt");
    //   ofstream Time(OutputFileName, ios::out | ios::trunc);
    //   Time << LiveTime << endl << RealTime << endl;
    //   Time.close();
    
    TPMERegexp("_time\\.txt").Substitute(OutputFileName," run");
    TH1D* hRaw = new TH1D( "Raw_spectrum", "Raw gamma spectrum", static_cast<Int_t>(line_number),0,static_cast<Double_t>(line_number));
    int NcalBins = floor((CalOrigin+static_cast<Double_t>(line_number)*CalSlope)/Binwidth_keV);
    TH1D* hCal = new TH1D("Cal_spectrum", "Calibrated gamma spectrum", NcalBins,0,(CalOrigin+static_cast<Double_t>(line_number)*CalSlope));
    
    int Ncounts_i;
    Double_t DSPCE_i;
    
    TRandom* rand = new TRandom();
    for (UInt_t i=0; i<=line_number; i++) {
        //     hRaw->Fill(i,vValue.at(i));
        //     hCal->Fill(i,vValue.at(i));
        Ncounts_i = vValue.at(i);
        for (int j=0; j<Ncounts_i; j++){
            DSPCE_i = hRaw->GetBinLowEdge(i+1)+rand->Uniform(hRaw->GetBinWidth(i+1));
            hRaw->Fill(DSPCE_i);
            hCal->Fill((CalOrigin+DSPCE_i*CalSlope));
        }
    }
    
    /*  TCanvas* Can = dynamic_cast<TCanvas*>(gROOT->FindObject("Can"));
     *  if (Can) { delete Can; }
     *  Can = new TCanvas("Can","No Title",0,0,1000,1000);*/
    
    hRaw->GetXaxis()->SetTitle("DSPEC bin [a.u.]");
    hRaw->GetXaxis()->CenterTitle();
    hRaw->GetYaxis()->SetTitle("Count number []");
    hRaw->Write();
    
    TH1D* hCal_DRU = (TH1D*) hCal->Clone();
    
    //Mass of the detector
    Double_t Volume =(5.5*TMath::Pi()*pow(5.55/2,2)-4.4*TMath::Pi()*pow(1.15/2,2));
    Double_t Mass = Volume*5.323/1000;
    
    Double_t Eq_time_d=LiveTime/(24*3600.);
    
    //Normalize the histogram
    //hCal_DRU->Rebin(10);
    hCal_DRU->Scale(1./(Mass*Eq_time_d*hCal_DRU->GetBinWidth(1))); //in ev/d/kg/keV
    hCal_DRU->GetXaxis()->SetTitle("Energy [keV]");
    hCal_DRU->GetXaxis()->CenterTitle();
    hCal_DRU->GetYaxis()->SetTitle("Integrated Rate [ev kg^{-1} s^{-1} keV^{-1}]");
    hCal_DRU->SetTitle("Cal_spectrum_DRU");
    hCal_DRU->Write("Cal_spectrum_DRU");
    
    hCal->Scale(1./LiveTime);
    hCal->GetXaxis()->SetTitle("Energy [keV]");
    hCal->GetXaxis()->CenterTitle();
    hCal->GetYaxis()->SetTitle("Count rate [s^{-1}]");
    hCal->Write();
    
    delete hRaw;
    delete hCal;
    
    OutputFile->Close();
    //delete OutputFile;
    return;
}

void CalFromRoot(TString InputFileName)
{
    // create output file
    TString OutputFileName = InputFileName;
    TPMERegexp("\\.root").Substitute(OutputFileName,"_Cal.root");
    
    TFile* InputFile = new TFile(InputFileName,"READ");
    int RealTime;
    int LiveTime;
    TObjString* ObjStrTime;
    TString StrTime;
    std::stringstream ss;
    double time;
    int NUM_RUNS=48;
    for (int i = 1; i<= NUM_RUNS; i++){
        ObjStrTime= (TObjString*) InputFile->Get(Form("RealTime;%d", i));
        StrTime=ObjStrTime->GetString();
        TPMERegexp("RealTime: ").Substitute(StrTime,"");
        TPMERegexp(" s").Substitute(StrTime,"");
        ss << StrTime;
        ss >> time;
        RealTime+=time;
        ss.clear();
        ObjStrTime= (TObjString*) InputFile->Get(Form("LiveTime;%d", i));
        StrTime=ObjStrTime->GetString();
        TPMERegexp("LiveTime: ").Substitute(StrTime,"");
        TPMERegexp(" s").Substitute(StrTime,"");
        ss << StrTime;
        ss >> time;
        LiveTime+=time;
        ss.clear();
    }

    cout <<"RealTime: " <<RealTime << " s, LiveTime: "<<LiveTime<<" s\n";

    TFile* OutputFile = new TFile(OutputFileName,"RECREATE");
    
    if (not OutputFile->IsOpen()) { cerr << "***Error: can not create ROOT file: " << OutputFileName << endl; return; }
    cout << "Your output ROOT file is " << OutputFile->GetName() << endl;

    TObjString ObjStrTime_out = Form("RealTime: %d s",RealTime);
    ObjStrTime_out.Write("RealTime");
    ObjStrTime_out = Form("LiveTime: %d s",LiveTime);
    ObjStrTime_out.Write("LiveTime");

    
    TH1D* hRaw = (TH1D*) InputFile-> Get("Raw_spectrum");
    
    int Nbins = hRaw->GetNbinsX();
    int FirstBin = hRaw->GetBinLowEdge(1);
    int LastBin = hRaw->GetBinLowEdge(Nbins)+hRaw->GetBinWidth(Nbins);
    double Cal_Slope = 2614./11640.;
    
    TH1D* hCal = new TH1D( "Cal_spectrum", "Calibrated gamma spectrum", Nbins, FirstBin*Cal_Slope, LastBin*Cal_Slope);
    
    for (UInt_t i=1; i<=Nbins; i++) {
        hCal->SetBinContent(i,hRaw->GetBinContent(i));
        hCal->SetBinError(i,hRaw->GetBinError(i));
    }
    hRaw->GetXaxis()->SetTitle("DSPEC bin [a.u.]");
    hRaw->GetXaxis()->CenterTitle();
    hRaw->GetYaxis()->SetTitle("Count number []");
    hRaw->Write();

    hCal->Scale(1./LiveTime);
    hCal->GetXaxis()->SetTitle("Energy [keV]");
    hCal->GetXaxis()->CenterTitle();
    hCal->GetYaxis()->SetTitle("Count rate [s^{-1}]");
    hCal->Write();

    delete hRaw;
    delete hCal;

    OutputFile->Close();
    //delete OutputFile;
    return;
}
