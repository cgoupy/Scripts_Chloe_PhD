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

TString Spec_name = "h1_NTD_1_MeV";

using namespace TMath;

std::vector<std::pair<Double_t,Double_t>> PeakFinder(TString InputFileName, UInt_t Npeaks=15, Double_t Thres = 0.1)
{
    // test input file
    struct stat buf; // to test input file
    if (stat(InputFileName,&buf) == -1) { // test whether the file exist or not, and get its size
        cerr << "***Error: file missing: " << InputFileName << "\nExiting...\n";
        exit(-1);
    }
    
    TFile* RunFile = new TFile(InputFileName,"READ");
    TPMERegexp("\\.root").Substitute(InputFileName,"_analyse.txt");
    ofstream Resultat(InputFileName, ios::out | ios::trunc);
    
    //Retrieve histogram of calibrated energies
    TH1D *hCal=0;
    hCal = (TH1D*) RunFile->Get(Spec_name);
    //hCal->Rebin(10);
    gStyle->SetOptStat(0);
    TPMERegexp("\\_analyse.txt").Substitute(InputFileName,"");
    TPMERegexp("\\w+\\/").Substitute(InputFileName,"");
    hCal->SetTitle(InputFileName);
    
    //Use TSpectrum to find the peak candidates
    TSpectrum *s = new TSpectrum(2*Npeaks+5);
    const Int_t Nfound = s->Search(hCal,1,"new",Thres);
    cout << "Found " << Nfound << " candidate peaks to fit\n";
    //hCal->Draw();
    //hCal->DrawCopy("E0 HIST SAME");

    
    //Results of TSpectrum
    Double_t *Xpeaks = s->GetPositionX();
    Double_t *Ypeaks = s->GetPositionY();
    sort(Xpeaks,Xpeaks+Nfound);
    
    std::vector<std::pair<Double_t,Double_t>> peakData;
    
    //Gross estimate of count rate in peaks + filling std::vector
    for (Int_t p=0;p<Nfound;p++) {
        Double_t xp = Xpeaks[p];
        Double_t yp = hCal->GetBinContent(hCal->GetXaxis()->FindBin(xp));
        peakData.push_back(std::make_pair(xp,yp));
        //Float_t xp = xpeaks[Selected[p]];
        // if(p>nfound) break;
        Double_t Signal=0;
        Double_t Background=0;
        Double_t Total=0;
        Total = hCal->Integral(hCal->FindBin(xp-0.3*xp),hCal->FindBin(xp+0.3*xp));
        Background = 3 * ( hCal->Integral(hCal->FindBin(xp-0.5*xp),hCal->FindBin(xp-0.3*xp)) + hCal->Integral(hCal->FindBin(xp+0.3*xp),hCal->FindBin(xp+0.5*xp)) );
        Signal = Total - Background;
        
        cout<<"Peak #" << p << " location: "<< xp << " - rough integral: " << Signal << " Bq" << endl;
        Resultat << xp << " " << Signal << endl;
    }
    
    //    //Estimate linear background
    //    TF1 *fLin = new TF1("flin","pol0(0)+expo(1)",0,4000);
    //    fLin->SetParameters(0,11.616);
    //    fLin->SetParameters(1,6.62);
    //    fLin->SetParameters(2,-0.002);
    //    hCal->Fit("flin","R");
    //
    //    //Loop on all found peaks and eliminate peaks at the background level
    //    Npeaks = 0;
    //    Int_t Selected[Nfound];
    //    for (Int_t p=0;p<Nfound;p++) {
    //        Float_t xp = Xpeaks[p];
    //        Float_t yp = hCal->GetBinContent(hCal->GetXaxis()->FindBin(xp));
    //        if (xp>50 && (xp < 150 || yp > 12+1.2*(fLin->Eval(xp))) ) {
    //            Selected[Npeaks]=p;
    //            //Xpeaks_saved[p]=xpeaks[p];
    //            Npeaks++;
    //        }
    //        if (xp>1500) {
    //            Selected[Npeaks]=p;
    //            //Xpeaks_saved[p]=xpeaks[p];
    //            Npeaks++;
    //        }
    //    }
    //
    //    cout << "Found " << Npeaks << " useful peaks to fit\n";
    
    
    //    cout << "Now fitting: Be patient\n";
    //
    //    //FIT
    //    TF1 **fit = new TF1*[Nfound];
    //    TF1* fitd[Nfound];
    //
    //    for (Int_t p=0;p<Nfound;p++) {
    //
    //        Float_t xp = Xpeaks[p];
    //        Float_t yp = hCal->GetBinContent(hCal->GetXaxis()->FindBin(xp));
    //        fit[p] = new TF1(Form("fit_%d",p),"gaus(0)+pol1(3)",xp-10,xp+10);
    //        //     fline[p] = new TF1(Form("line_%d",p),"pol1",xp-20,xp+20);
    //        //Fit parameter initialization
    //        Double_t par[5];
    //        par[0] = yp;
    //        par[1] = xp;
    //        par[2] = 1.09;
    //        par[3] = 300;
    //        par[4] = -0.9;
    //        fit[p]->SetParameters(par);
    //        fit[p]->SetParLimits(0,par[0]-5,par[0]+1);
    //        fit[p]->SetParLimits(1,xp-0.1,xp+0.1);
    //        fit[p]->SetParLimits(2,par[2]-0.9, par[2]+0.3);
    //        fit[p]->SetNpx(1000);
    //
    //        hCal->Fit(Form("fit_%d",p),"R0");
    //        cout << "Peak #" << p << " location: " << xp << endl;
    //    }
    
    //TPMERegexp(" run").Substitute(InputFileName,"_analyse.txt");
    
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    Resultat.close();
    return peakData;
}


void PeakFitter(TString FileName = "../Data_COV/Norm_MyCOV1_v3_36.root", Int_t PeakNumber=0){
    
    
    //Mass of the HPGe detector
    //Double_t M = (5.5*TMath::Pi()*pow(5.55/2,2)-4.4*TMath::Pi()*pow(1.15/2,2))*5.323/1000;
    
    std::vector<std::pair<Double_t,Double_t>> Peaks = PeakFinder(FileName,100, 0.00001);
    TH1D* h = (TH1D*)gDirectory->Get(Spec_name);
    
    //Draw full spectrum
    TH1D* hh = (TH1D*) h->Clone();
    //hh->SetName("Gamma spectrum");
    hh->SetTitle("Peak fit");
    //hh->SetTitle("Co^{60}");
    hh->SetStats(0);
    
    const Int_t NPeaks = Peaks.size();
    
    TCanvas *CanHisto = new TCanvas("CanHisto","Canvas of energy histogram",0,0,1000,700);
    CanHisto->cd();
    hh->Draw("E0 HIST");
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    
    //Fit peak
    //Int_t PeakNumber = 14;
    
    Double_t Xp = Peaks.at(PeakNumber).first;
    Double_t Yp = Peaks.at(PeakNumber).second;
    
    TF1* fit = new TF1("Peak fit","gaus(0)+pol1(3)",Xp-0.15*Xp,Xp+0.15*Xp);
    fit->SetParNames("A","E_{rec}","\\sigma","p_{0}","p_{1}");
    
    //Fit parameter initialization
    Double_t par[5];
    par[0] = Yp;
    par[1] = Xp;
    par[2] = 0.001;
    par[3] = 0.2;
    par[4] = -0.8;
    fit->SetParameters(par);
    //fit->SetParLimits(0,Yp-0.6*Yp,Yp+0.6*Yp);
    //fit->SetParLimits(1,Xp-0.6*Xp,Xp+0.6*Xp);
    //fit->SetParLimits(2,0.1,2);
    fit->SetNpx(2000);
    
    //Draw fit
    TCanvas *CanFit = new TCanvas("CanFit","Canvas of peak fit",0,0,1000,900);
    CanFit->cd();
    h->GetXaxis()->SetRangeUser(Xp-0.5*Xp,Xp+0.5*Xp);
    h->Draw("E0 HIST");
    h->SetTitle(Form("COV1 - %g keV", Xp*1000));
    //h->SetTitle(Form("Co^{60} - %g keV", Xp));
    gStyle->SetOptStat(0);
    TFitResultPtr r = h->Fit("Peak fit","RS");
    gStyle->SetOptFit(111);
    fit->Draw("LSAME");
    fit->SetLineWidth(2);
    h->SetLineWidth(2);
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    
    r->Print("V");
    TMatrixDSym cov = r->GetCovarianceMatrix();
    TMatrixDSym covPeak(0,2,cov.GetMatrixArray());
    
    //Estimate Area under peak
    //    cout << fit->GetParameter(0) << endl;
    //    cout << fit->GetParameter(1) << endl;
    //    cout << fit->GetParameter(2) << endl;
    
    TF1* PeakFunc = new TF1("Peak function","gaus",Xp-0.3*Xp,Xp+0.3*Xp);
    PeakFunc->SetParameter(0,fit->GetParameter(0));
    PeakFunc->SetParameter(1,fit->GetParameter(1));
    PeakFunc->SetParameter(2,fit->GetParameter(2));
    
    Double_t A = PeakFunc->Integral(fit->GetParameter(1)-100*fit->GetParameter(2),fit->GetParameter(1)+100*fit->GetParameter(2));
    Double_t Err = PeakFunc->IntegralError(fit->GetParameter(1)-100*fit->GetParameter(2),fit->GetParameter(1)+100*fit->GetParameter(2),PeakFunc->GetParameters(),covPeak.GetMatrixArray());
    
    std::cout << "Count rate for peak #" << PeakNumber << ": " << A << " +/- " << Err << " Bq" << std::endl;
    
}

double sigma1(double x, TMatrixTSym<double> CovMat){
    double sig1 = CovMat(0, 0) + CovMat(0,1)*x + CovMat(1,0)*x + CovMat(1,1)*x*x;
    return sqrt(sig1);
}

class classsigma1 {
public:
    classsigma1(TMatrixTSym<double> CovMat): CovMat(CovMat){
//        cout<<"New sigma1 defined with cov matrix = "<<endl;
//        CovMat.Print();
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

std::vector<TFitResultPtr> PlotEnergyCalibration(){
    
    int Emax = 20000;
    //Calibration points from Eu-152 source
    //Double_t Etrue[14] = {39.52, 121.782, 244.697, 295.697, 344.2785, 411.12, 443.961, 778.90, 867.38, 964.06, 1085.84, 1112.076, 1299.14, 1408.01};
    //Double_t Erec[14] = {39.86, 121.85, 244.71, 295.92, 344.44, 411.14, 443.94, 779.05, 867.55, 964.35, 1086.09, 1112.37, 1299.46, 1408.4};
    //Double_t ChanNum[14]={177, 543, 1090, 1317, 1533, 1830, 1978, 3468, 3862, 4293, 4835, 4952, 5785, 6270};
    
    std::vector<TFitResultPtr> FitResults;
    
    //Calibration points COV Bot
    Double_t Etrue[5] = {59.54, 511, 1460.82, 2614.511, 15690};
    Double_t EtrueErr[5] = {0, 0, 0, 0, 10};
    Double_t ChanNum[5]={0.05703, 0.4866, 1.39, 2.489, 15.49};
    Double_t ChanNumErr[5]={0.00001, 0.0004, 0.001, 0.001, 0.09};
    
    TGraphErrors* gr = new TGraphErrors(5,Etrue,ChanNum, EtrueErr, ChanNumErr);
    //gr->SetTitle("Energy scale linearity from Co60");
    gr->SetTitle("Calibration Bot Ge");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(2);
    gr->SetMarkerSize(0.6);
    gr->GetXaxis()->SetTitle("True energy [keV]");
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetXaxis()->SetRangeUser(0,20000);
    gr->GetYaxis()->SetTitle("Amplitude [V]");
    //gr->GetYaxis()->SetTitle("DESPEC bin [a.u.]");
    gr->GetYaxis()->SetTitleOffset(1.4);
    //gr->GetYaxis()->SetRangeUser(0,1500);
    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1);
    
    
    TCanvas *c1 = new TCanvas("c1","Energy linearity",200,10,1200,1200);
//    gr->GetXaxis()->SetRangeUser(-5, 5);
//    gr->GetYaxis()->SetRangeUser(-5, 5);
    gPad->SetGridx();
    gPad->SetGridy();
    
    
    //Fit
    TF1* COV1lin1 =  new TF1("COV1lin1","([0]+[1]*x)",0,Emax);
    COV1lin1->SetParNames("b [#V]","a [#V/keV]");
    
    //Fit parameter initialization
    COV1lin1->SetParameter(0, 0);
    COV1lin1->SetParameter(1, 1);
    //fit->SetParLimits(0,-0.0001, 0.0001);
    //fit->SetParLimits(1,0.9,1.1);

    TFitResultPtr rlin1 = gr->Fit("COV1lin1","RS0", "", 0, Emax);
    FitResults.push_back(rlin1);
    TMatrixTSym<double> COV1_rlin1_Cov = rlin1->GetCovarianceMatrix();
    COV1_rlin1_Cov.Print();
    
    /*Create a histogram to hold the confidence intervals*/
    classsigma1* sigma1 = new classsigma1(COV1_rlin1_Cov);
    TH1D *COV1_1sig = new TH1D("COV1_1sig",
       "COV1_1sig", Emax, 0, Emax);
    for (int i=0; i<COV1_1sig->GetNbinsX(); i++){
        double x = COV1_1sig->GetBinCenter(i);
        COV1_1sig->SetBinContent(i, COV1lin1->Eval(x));
        double x_vect[1]={x};
        COV1_1sig->SetBinError(i, 100*sigma1->Evaluate(x_vect));
    }
    
    COV1_1sig->SetMarkerSize(0);
    COV1_1sig->SetMarkerColor(0);
    COV1_1sig->SetLineColor(0);
    COV1_1sig->SetStats(false);
    COV1_1sig->SetFillColorAlpha(2,0.4);
    //COV1_1sig->Draw("e3 same");
    COV1_1sig->GetXaxis()->SetTitle("True energy [keV]");
    COV1_1sig->GetYaxis()->SetTitle("Amplitude [V]");
    COV1_1sig->SetTitle("100#times#sigma");
    
    COV1lin1->SetLineWidth(3);
    COV1lin1->SetLineColor(2);
    COV1lin1->SetLineStyle(9);
    COV1lin1->Draw("Same");
    COV1lin1->SetTitle(Form("#splitline{a = %0.1f +/- %0.1f #muV/keV    }{b = %0.1f +/- %0.1f #muV}", 1e6*COV1lin1->GetParameter(1), 1e6*COV1lin1->GetParError(1), 1e6*COV1lin1->GetParameter(0), 1e6*COV1lin1->GetParError(0)));
    
    //gr->Draw("P same");
    auto legend1 = new TLegend();
    legend1->AddEntry(gr);
    legend1->AddEntry(COV1lin1);
    legend1->AddEntry(COV1_1sig);
    legend1->Draw();

//    TF1* fit = new TF1("COV1","[0]+[1]*x+[2]*x*x",0,20000);
//    fit->SetParNames("c", "b", "a");
//
//    //Fit parameter initialization
//    fit->SetParameter(0, 0);
//    fit->SetParameter(1, 1);
//    fit->SetParameter(2, 1);
//    //fit->SetParLimits(0,-0.0001, 0.0001);
//    //fit->SetParLimits(1,0.9,1.1);

//    gStyle->SetOptStat(0);
//    TFitResultPtr r = gr->Fit("COV1","S0", "", 0, 20000);
//    gStyle->SetOptFit(111);
//    fit->Draw("LSAME");
//    fit->SetLineWidth(2);
//    fit->SetLineColor(1);
//    fit->SetLineStyle(9);

   //Calibration points COV2
    Double_t Etrue2[4] = {511, 1460.82, 2614.511, 15700};
    Double_t EtrueErr2[4] = {0, 0, 0, 0};
    Double_t ChanNum2[4]={0.3318, 0.9491, 1.701, 10.26};
    Double_t ChanNumErr2[4]={0.0003, 0.0001, 0.001, 0.04};

    TGraphErrors* gr2 = new TGraphErrors(4,Etrue2,ChanNum2, EtrueErr2, ChanNumErr2);
    //gr->SetTitle("Energy scale linearity from Co60");
    gr2->SetTitle("Calibration Top Ge");
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(6);
    gr2->SetMarkerSize(0.6);
    gr2->GetXaxis()->SetTitle("True energy [keV]");
    gr2->GetXaxis()->SetTitleOffset(1.2);
    gr2->GetXaxis()->SetRangeUser(0,20000);
    gr2->GetYaxis()->SetTitle("Amplitude [V]");
    //gr->GetYaxis()->SetTitle("DESPEC bin [a.u.]");
    gr2->GetYaxis()->SetTitleOffset(1.4);
    //gr->GetYaxis()->SetRangeUser(0,1500);
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerSize(1);

    //gr2->Draw("PSAME");

    //Fit
    TF1* COV2lin = new TF1("COV2","([0]+[1]*x)",0,20000);
    COV2lin->SetParNames("b [#V]","a [#V/keV]");
    Double_t par[2];
    par[0]=0;
    par[1]=1;

    //Fit parameter initialization
    COV2lin->SetParameters(par);
    //fit2->SetParLimits(0,-10,10);
    //fit2->SetParLimits(1,0.9,1.1);

    gStyle->SetOptStat(0);
    TFitResultPtr r2 = gr2->Fit("COV2","RS0");
    FitResults.push_back(r2);
    TMatrixTSym<double> COV2_rlin_Cov = r2->GetCovarianceMatrix();
    COV2_rlin_Cov.Print();
    
    /*Create a histogram to hold the confidence intervals*/
    TH1D *COV2_1sig = new TH1D("COV2_1sig",
       "COV2_1sig", Emax, 0, Emax);
    sigma1= new classsigma1(COV2_rlin_Cov);
    
    for (int i=0; i<COV2_1sig->GetNbinsX(); i++){
        double x = COV2_1sig->GetBinCenter(i);
        COV2_1sig->SetBinContent(i, COV2lin->Eval(x));
        double x_vect[1]={x};
        COV2_1sig->SetBinError(i, 100*sigma1->Evaluate(x_vect));
    }
    
    COV2_1sig->SetMarkerSize(0);
    COV2_1sig->SetMarkerColor(0);
    COV2_1sig->SetLineColor(0);
    COV2_1sig->SetStats(false);
    COV2_1sig->SetFillColorAlpha(6,0.4);
    //COV2_1sig->Draw("e3 same");
    COV2_1sig->SetTitle("100#times#sigma");
    
    COV2lin->SetLineWidth(3);
    COV2lin->SetLineColor(6);
    COV2lin->SetLineStyle(9);
    //COV2lin->Draw("Same");
    COV2lin->SetTitle(Form("#splitline{a = %0.1f +/- %0.1f #muV/keV    }{b = %0.1f +/- %0.1f #muV}", 1e6*COV2lin->GetParameter(1), 1e6*COV2lin->GetParError(1), 1e6*COV2lin->GetParameter(0), 1e6*COV2lin->GetParError(0)));
    
    //gr2->Draw("P same");
    auto legend2 = new TLegend();
    legend2->AddEntry(gr2);
    legend2->AddEntry(COV2lin);
    legend2->AddEntry(COV2_1sig);
    legend2->Draw("same");

    //Calibration points LWO
//     Double_t Etrue3[4] = {511, 1460.82, 2614.511, 15180};
//     Double_t EtrueErr3[4] = {0, 0, 0, 20};
//    Double_t ChanNum3[4]={0.1497, 0.3741, 0.664, 3.989};
//     Double_t ChanNumErr3[4]={0.01, 0.0007, 0.002, 0.1};

//    //No shielding
//     Double_t Etrue3[3] = {1460.82, 2614.511, 15180};
//     Double_t EtrueErr3[3] = {0, 0, 20};
//    Double_t ChanNum3[3]={0.3741, 0.664, 3.989};
//     Double_t ChanNumErr3[3]={0.0007, 0.002, 0.1};
    
    //No shielding
     Double_t Etrue3[3] = {511, 1460.82, 2614};
     Double_t EtrueErr3[3] = {0, 0, 20};
    Double_t ChanNum3[3]={0.1199, 0.3505, 0.6227};
     Double_t ChanNumErr3[3]={0.0012, 0.0003, 0.0011};
    
     TGraphErrors* gr3 = new TGraphErrors(3,Etrue3,ChanNum3, EtrueErr3, ChanNumErr3);
     //gr->SetTitle("Energy scale linearity from Co60");
    gr3->SetTitle("Calibration Mid LWO");
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(8);
    gr3->SetMarkerSize(0.6);
    gr3->GetXaxis()->SetTitle("True energy [keV]");
    gr3->GetXaxis()->SetTitleOffset(1.2);
    gr3->GetXaxis()->SetRangeUser(0,20000);
    gr3->GetYaxis()->SetTitle("Amplitude [V]");
     //gr->GetYaxis()->SetTitle("DESPEC bin [a.u.]");
    gr3->GetYaxis()->SetTitleOffset(1.4);
     //gr->GetYaxis()->SetRangeUser(0,1500);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1);

     gr3->Draw("AP");

     //Fit
     TF1* fitLWO = new TF1("LWO","([0]+[1]*x)",0,20000);
    fitLWO->SetParNames("b [V]","a [V/keV]");

     //Fit parameter initialization
    fitLWO->SetParameters(par);
     //fit2->SetParLimits(0,-10,10);
    //fit2->SetParLimits(1,-1000,1000);

    gStyle->SetOptStat(0);
    TFitResultPtr r3 = gr3->Fit("LWO","RS0");
    FitResults.push_back(r3);

    TMatrixTSym<double> LWO_rlin_Cov = r3->GetCovarianceMatrix();
    LWO_rlin_Cov.Print();
    
    /*Create a histogram to hold the confidence intervals*/
    TH1D *LWO_1sig = new TH1D("LWO_1sig",
       "LWO_1sig", Emax, 0, Emax);
    sigma1= new classsigma1(LWO_rlin_Cov);
    for (int i=0; i<LWO_1sig->GetNbinsX(); i++){
        double x = LWO_1sig->GetBinCenter(i);
        LWO_1sig->SetBinContent(i, fitLWO->Eval(x));
        double x_vect[1]={x};
        LWO_1sig->SetBinError(i, 100*sigma1->Evaluate(x_vect));
    }
    
    LWO_1sig->SetMarkerSize(0);
    LWO_1sig->SetMarkerColor(0);
    LWO_1sig->SetLineColor(0);
    LWO_1sig->SetStats(false);
    LWO_1sig->SetFillColorAlpha(8,0.4);
    LWO_1sig->Draw("e3 same");
    LWO_1sig->SetTitle("10#times#sigma");
    
    fitLWO->SetLineWidth(3);
    fitLWO->SetLineColor(8);
    fitLWO->SetLineStyle(9);
    fitLWO->Draw("Same");
    fitLWO->SetTitle(Form("#splitline{a = %0.1f +/- %0.1f #muV/keV    }{b = %0.1f +/- %0.1f #muV}", 1e6*fitLWO->GetParameter(1), 1e6*fitLWO->GetParError(1), 1e6*fitLWO->GetParameter(0), 1e6*fitLWO->GetParError(0)));
    
    gr3->Draw("P same");
    auto legend3 = new TLegend();
    legend3->AddEntry(gr3);
    legend3->AddEntry(fitLWO);
    legend3->AddEntry(LWO_1sig);
    legend3->Draw("same");

    return FitResults;
}

void Plot_error_Calibration(std::vector<TFitResultPtr> Fit_results){

    int Emax = 20000;

    /*COV1*/
    TF1* Func_fit1 = new TF1("Func_fit1","([0]+[1]*x)",0,20000);
    Func_fit1->SetParNames("b [#V]","a [#V/keV]");
    Func_fit1->FixParameter(0,Fit_results[0]->Parameter(0));
    Func_fit1->FixParameter(1,Fit_results[0]->Parameter(1));
    
    auto COV1rlin_Cov = Fit_results[0]->GetCovarianceMatrix();
    
    classsigma1* sigma1 = new classsigma1(COV1rlin_Cov);
    
    TF1* func_p_Sigma1 = new TF1("func_p_Sigma1", sigma1, &classsigma1::EvaluatePlusF, 0, Emax, 2, "classsigma1", "EvaluatePlusF");
    func_p_Sigma1->SetParNames("b [#V]","a [#V/keV]");
    func_p_Sigma1->FixParameter(0,Fit_results[0]->Parameter(0));
    func_p_Sigma1->FixParameter(1,Fit_results[0]->Parameter(1));
    
    TF1* func_m_Sigma1 = new TF1("func_m_Sigma1", sigma1, &classsigma1::EvaluateMinusF, 0, Emax, 2, "classsigma1", "EvaluateMinusF");
    func_m_Sigma1->SetParNames("b [#V]","a [#V/keV]");
    func_m_Sigma1->FixParameter(0,Fit_results[0]->Parameter(0));
    func_m_Sigma1->FixParameter(1,Fit_results[0]->Parameter(1));
    
    cout<<"b = "<<func_p_Sigma1->GetParameter(0)<<endl;
    cout<<"a = "<<func_p_Sigma1->GetParameter(1)<<endl;
    
    TF1* Sigma1 = new TF1("Sigma1", sigma1, &classsigma1::Evaluate, 0, Emax, 0, "classsigma1", "Evaluate");
    
    TF1* Ereal = new TF1("Ereal","x",0,20000);
    
    TH1D *COV1 = new TH1D("COV1", "COV1", Emax, 0, Emax);
    for (int i=0; i<COV1->GetNbinsX(); i++){
        double Ereal = COV1->GetBinCenter(i+1);
        COV1->SetBinContent(i+1, Ereal);
//        cout<<"Ereal="<<Ereal<<endl;
//        cout<<"V="<<Func_fit1->Eval(Ereal)<<endl;
//        cout<<"E1="<<Ereal-func_p_Sigma1->GetX(Func_fit1->Eval(Ereal), Ereal-500, Ereal+500)<<endl;
    }
    
    TGraphAsymmErrors *COV1err = new TGraphAsymmErrors(COV1);
    for (int i=0; i<COV1->GetNbinsX(); i++){
        double Ereal = COV1->GetBinCenter(i);
        COV1err->SetPointEYlow(i, 100*(Ereal-func_p_Sigma1->GetX(Func_fit1->Eval(Ereal), Ereal-500, Ereal+500)));
        COV1err->SetPointEYhigh(i, 100*(func_m_Sigma1->GetX(Func_fit1->Eval(Ereal), Ereal-500, Ereal+500)-Ereal));
    }
    
    
    TCanvas* c2 = new TCanvas("COV1", "COV1");
    //COV1err->SetMarkerSize(0);
    COV1err->SetMarkerColor(2);
    COV1err->SetLineColor(0);
    COV1err->SetFillColorAlpha(2,0.4);
    COV1err->Draw("APe3");
    COV1err->GetXaxis()->SetRangeUser(0, 20000);
    Ereal->SetLineWidth(3);
    Ereal->SetLineColor(2);
    Ereal->SetLineStyle(9);
    Ereal->Draw("Same");
    COV1err->GetXaxis()->SetTitle("True energy [keV]");
    COV1err->GetYaxis()->SetTitle("Calibrated Energy [keV]");
    COV1err->SetTitle("100#times#sigma");
    
    /*COV2*/
    TF1* Func_fit2 = new TF1("Func_fit2","([0]+[1]*x)",0,20000);
    Func_fit2->SetParNames("b [#V]","a [#V/keV]");
    Func_fit2->FixParameter(0,Fit_results[1]->Parameter(0));
    Func_fit2->FixParameter(1,Fit_results[1]->Parameter(1));
    
    auto COV2rlin_Cov = Fit_results[1]->GetCovarianceMatrix();
    
    sigma1 = new classsigma1(COV2rlin_Cov);
    
    TF1* func_p_Sigma1_2 = new TF1("func_p_Sigma1_2", sigma1, &classsigma1::EvaluatePlusF, 0, Emax, 2, "classsigma1", "EvaluatePlusF");
    func_p_Sigma1_2->SetParNames("b [#V]","a [#V/keV]");
    func_p_Sigma1_2->FixParameter(0,Fit_results[1]->Parameter(0));
    func_p_Sigma1_2->FixParameter(1,Fit_results[1]->Parameter(1));
    
    TF1* func_m_Sigma1_2 = new TF1("func_m_Sigma1_2", sigma1, &classsigma1::EvaluateMinusF, 0, Emax, 2, "classsigma1", "EvaluateMinusF");
    func_m_Sigma1_2->SetParNames("b [#V]","a [#V/keV]");
    func_m_Sigma1_2->FixParameter(0,Fit_results[1]->Parameter(0));
    func_m_Sigma1_2->FixParameter(1,Fit_results[1]->Parameter(1));
    
    cout<<"b = "<<func_m_Sigma1_2->GetParameter(0)<<endl;
    cout<<"a = "<<func_m_Sigma1_2->GetParameter(1)<<endl;
    
    TF1* Sigma1_2 = new TF1("Sigma1_2", sigma1, &classsigma1::Evaluate, 0, Emax, 0, "classsigma1", "Evaluate");
    
    TH1D *COV2 = new TH1D("COV2", "COV2", Emax, 0, Emax);
    for (int i=0; i<COV2->GetNbinsX(); i++){
        double Ereal_val = COV2->GetBinCenter(i+1);
        COV2->SetBinContent(i+1, Ereal_val);
//        cout<<"Ereal="<<Ereal<<endl;
//        cout<<"V="<<Func_fit1->Eval(Ereal)<<endl;
//        cout<<"E1="<<Ereal-func_p_Sigma1->GetX(Func_fit1->Eval(Ereal), Ereal-500, Ereal+500)<<endl;
    }
    
    TGraphAsymmErrors *COV2err = new TGraphAsymmErrors(COV2);
    for (int i=0; i<COV2->GetNbinsX(); i++){
        double Ereal_val = COV2->GetBinCenter(i);
        COV2err->SetPointEYlow(i, 100*(Ereal_val-func_p_Sigma1_2->GetX(Func_fit2->Eval(Ereal_val), Ereal_val-500, Ereal_val+500)));
        COV2err->SetPointEYhigh(i, 100*(func_m_Sigma1_2->GetX(Func_fit2->Eval(Ereal_val), Ereal_val-500, Ereal_val+500)-Ereal_val));
    }
    
    auto Ereal_2 = (TF1*) Ereal->Clone();
    
    TCanvas* c3 = new TCanvas("COV2", "COV2");
    //COV1err->SetMarkerSize(0);
    COV2err->SetMarkerColor(6);
    COV2err->SetLineColor(0);
    COV2err->SetFillColorAlpha(6,0.4);
    COV2err->Draw("APe3");
    COV2err->GetXaxis()->SetRangeUser(0, 20000);
    Ereal_2->SetLineWidth(3);
    Ereal_2->SetLineStyle(9);
    Ereal_2->SetLineColor(6);
    Ereal_2->DrawCopy("Same");
    COV2err->GetXaxis()->SetTitle("True energy [keV]");
    COV2err->GetYaxis()->SetTitle("Calibrated Energy [keV]");
    COV2err->SetTitle("100#times#sigma");
    
    /*NTD*/
    TF1* Func_fit3 = new TF1("Func_fit3","([0]+[1]*x)",0,20000);
    Func_fit3->SetParNames("b [#V]","a [#V/keV]");
    Func_fit3->FixParameter(0,Fit_results[2]->Parameter(0));
    Func_fit3->FixParameter(1,Fit_results[2]->Parameter(1));
    
    auto COV3rlin_Cov = Fit_results[2]->GetCovarianceMatrix();
    
    sigma1 = new classsigma1(COV3rlin_Cov);
    
    TF1* func_p_Sigma1_3 = new TF1("func_p_Sigma1_3", sigma1, &classsigma1::EvaluatePlusF, 0, Emax, 2, "classsigma1", "EvaluatePlusF");
    func_p_Sigma1_3->SetParNames("b [#V]","a [#V/keV]");
    func_p_Sigma1_3->FixParameter(0,Fit_results[2]->Parameter(0));
    func_p_Sigma1_3->FixParameter(1,Fit_results[2]->Parameter(1));
    
    TF1* func_m_Sigma1_3 = new TF1("func_m_Sigma1_3", sigma1, &classsigma1::EvaluateMinusF, 0, Emax, 2, "classsigma1", "EvaluateMinusF");
    func_m_Sigma1_3->SetParNames("b [#V]","a [#V/keV]");
    func_m_Sigma1_3->FixParameter(0,Fit_results[2]->Parameter(0));
    func_m_Sigma1_3->FixParameter(1,Fit_results[2]->Parameter(1));
    
    cout<<"b = "<<func_m_Sigma1_3->GetParameter(0)<<endl;
    cout<<"a = "<<func_m_Sigma1_3->GetParameter(1)<<endl;
    
    TF1* Sigma1_3 = new TF1("Sigma1_3", sigma1, &classsigma1::Evaluate, 0, Emax, 0, "classsigma1", "Evaluate");
    
    TH1D *COV3 = new TH1D("COV3", "COV3", Emax, 0, Emax);
    for (int i=0; i<COV3->GetNbinsX(); i++){
        double Ereal_val = COV3->GetBinCenter(i+1);
        COV3->SetBinContent(i+1, Ereal_val);
//        cout<<"Ereal="<<Ereal<<endl;
//        cout<<"V="<<Func_fit1->Eval(Ereal)<<endl;
//        cout<<"E1="<<Ereal-func_p_Sigma1->GetX(Func_fit1->Eval(Ereal), Ereal-500, Ereal+500)<<endl;
    }
    
    TGraphAsymmErrors *COV3err = new TGraphAsymmErrors(COV3);
    for (int i=0; i<COV3->GetNbinsX(); i++){
        double Ereal_val = COV3->GetBinCenter(i);
        COV3err->SetPointEYlow(i, 10*(Ereal_val-func_p_Sigma1_3->GetX(Func_fit3->Eval(Ereal_val), Ereal_val-500, Ereal_val+500)));
        COV3err->SetPointEYhigh(i, 10*(func_m_Sigma1_3->GetX(Func_fit3->Eval(Ereal_val), Ereal_val-500, Ereal_val+500)-Ereal_val));
    }
    
    auto Ereal_3 = (TF1*) Ereal->Clone();
    
    TCanvas* c4 = new TCanvas("COV3", "COV3");
    //COV1err->SetMarkerSize(0);
    COV3err->SetMarkerColor(8);
    COV3err->SetLineColor(0);
    COV3err->SetFillColorAlpha(8,0.4);
    COV3err->Draw("APe3");
    COV3err->GetXaxis()->SetRangeUser(0, 20000);
    Ereal_3->SetLineWidth(3);
    Ereal_3->SetLineStyle(9);
    Ereal_3->SetLineColor(8);
    Ereal_3->DrawCopy("Same");
    COV3err->GetXaxis()->SetTitle("True energy [keV]");
    COV3err->GetYaxis()->SetTitle("Calibrated Energy [keV]");
    COV3err->SetTitle("10#times#sigma");
}


//-------------- High
double EHigh(double Ereal, TFitResultPtr Fit_results){

    TF1* Func_fit = new TF1("Func_fit","([0]+[1]*x)",0,20000);
    Func_fit->SetParNames("b [#V]","a [#V/keV]");
    Func_fit->FixParameter(0,Fit_results->Parameter(0));
    Func_fit->FixParameter(1,Fit_results->Parameter(1));

    auto COVrlin_Cov = Fit_results->GetCovarianceMatrix();
    
    classsigma1* sigma1 = new classsigma1(COVrlin_Cov);
    
    TF1* func_m_Sigma1 = new TF1("func_m_Sigma1", sigma1, &classsigma1::EvaluateMinusF, 0, 20000, 2, "classsigma1", "EvaluateMinusF");
    func_m_Sigma1->SetParNames("b [#V]","a [#V/keV]");
    func_m_Sigma1->FixParameter(0,Fit_results->Parameter(0));
    func_m_Sigma1->FixParameter(1,Fit_results->Parameter(1));
    
    TF1* Sigma1 = new TF1("Sigma1", sigma1, &classsigma1::Evaluate, 0, 20000, 0, "classsigma1", "Evaluate");
    
    double Ehigh = Ereal-func_m_Sigma1->GetX(Func_fit->Eval(Ereal), Ereal-500, Ereal+500);
    return Ehigh;
    
}

//------------ Low
double ELow(double Ereal, TFitResultPtr Fit_results){

    TF1* Func_fit = new TF1("Func_fit","([0]+[1]*x)",0,20000);
    Func_fit->SetParNames("b [#V]","a [#V/keV]");
    Func_fit->FixParameter(0,Fit_results->Parameter(0));
    Func_fit->FixParameter(1,Fit_results->Parameter(1));

    auto COVrlin_Cov = Fit_results->GetCovarianceMatrix();
    
    classsigma1* sigma1 = new classsigma1(COVrlin_Cov);
    
    TF1* func_p_Sigma1 = new TF1("func_p_Sigma1", sigma1, &classsigma1::EvaluatePlusF, 0, 20000, 2, "classsigma1", "EvaluatePlusF");
    func_p_Sigma1->SetParNames("b [#V]","a [#V/keV]");
    func_p_Sigma1->FixParameter(0,Fit_results->Parameter(0));
    func_p_Sigma1->FixParameter(1,Fit_results->Parameter(1));
    
    TF1* Sigma1 = new TF1("Sigma1", sigma1, &classsigma1::Evaluate, 0, 20000, 0, "classsigma1", "Evaluate");
    
    double Elow = func_p_Sigma1->GetX(Func_fit->Eval(Ereal), Ereal-500, Ereal+500)-Ereal;
    return Elow;
    
}


TGraphAsymmErrors* Calculate_Error_syst(TH1D* h0, TH1D* h1, TH1D* h2){
    TGraphAsymmErrors* h_syst = new TGraphAsymmErrors(h0);
    
    double E1, E2;
    for (int i =0; i<h0->GetNbinsX(); i++){
        E1 = h0->GetBinContent(i+1)-h1->GetBinContent(i+1);
        E2 = h0->GetBinContent(i+1)-h2->GetBinContent(i+1);
        
        if (E1*E2<0){
            if(E1>0){
                h_syst->SetPointEYlow(i, h0->GetBinError(i+1)+E1);
                h_syst->SetPointEYhigh(i, h0->GetBinError(i+1)-E2);
            }
            else{
                h_syst->SetPointEYhigh(i, h0->GetBinError(i+1)-E1);
                h_syst->SetPointEYlow(i, h0->GetBinError(i+1)+E2);
            }
        }
        else {
            if(E1>0){
                if(E1<E2){
                    h_syst->SetPointEYlow(i, h0->GetBinError(i+1));
                    h_syst->SetPointEYhigh(i, h0->GetBinError(i+1)+E2);
                }
                else{
                    h_syst->SetPointEYlow(i, h0->GetBinError(i+1));
                    h_syst->SetPointEYhigh(i, h0->GetBinError(i+1)+E1);
                }
            }
            else{
                if(E1<E2){
                    h_syst->SetPointEYlow(i, h0->GetBinError(i+1)-E1);
                    h_syst->SetPointEYhigh(i, h0->GetBinError(i+1));
                }
                else{
                    h_syst->SetPointEYlow(i, h0->GetBinError(i+1)-E2);
                    h_syst->SetPointEYhigh(i, h0->GetBinError(i+1));
                }
            }
        }
    }
    return h_syst;
}


//---------------------- Shift an histogram with sigma function -------------------
void Shift(TFile* NewFile, TString Folder, TString Name, TH1D* histo, TFitResultPtr Fit_result){
    
    Name = Name+"_shifted";
    
    int Nbins = histo->GetNbinsX();
    cout<<histo->GetBinCenter(Nbins)+histo->GetBinWidth(Nbins)/2.<<endl;
    TH1D* histo_shifted_high = new TH1D(Name+"high", Name+"high", Nbins, 0, histo->GetBinCenter(Nbins)+histo->GetBinWidth(Nbins)/2.);
    TH1D* histo_shifted_low = new TH1D(Name+"low", Name+"low", Nbins, 0, histo->GetBinCenter(Nbins)+histo->GetBinWidth(Nbins)/2.);
    histo_shifted_high->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    histo_shifted_high->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    histo_shifted_low->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    histo_shifted_low->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    
    double E_i;
    double E_shifted_high, E_shifted_low;

    for (int i=0; i<1e5; i++){
        if (i%10000==0){cout<<i<<endl;}
        E_i = histo->GetRandom();
        E_shifted_high = E_i+EHigh(E_i*1e3, Fit_result)*1e-3;
        E_shifted_low = E_i-ELow(E_i*1e3, Fit_result)*1e-3;
        histo_shifted_high->Fill(E_shifted_high);
        histo_shifted_low->Fill(E_shifted_low);
    }

    //histo_shifted->Draw();
    histo_shifted_high->Scale(histo->Integral()/histo_shifted_high->Integral());
    histo_shifted_low->Scale(histo->Integral()/histo_shifted_low->Integral());

    for (int i=1; i<=histo_shifted_low->GetNbinsX(); i++){
        histo_shifted_low->SetBinError(i, 0);
        histo_shifted_high->SetBinError(i, 0);
    }
    
    NewFile->cd(Folder+"Evis_High");
    histo_shifted_high->Write(Name);
    
    NewFile->cd(Folder+"Evis_Low");
    histo_shifted_low->Write(Name);
    
    NewFile->cd(Folder+"Evis");
    TGraphAsymmErrors* h_syst = Calculate_Error_syst(histo, histo_shifted_high, histo_shifted_low);
    h_syst->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    h_syst->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    h_syst->SetMarkerColor(8);
    h_syst->SetFillColorAlpha(8,0.4);
    //h_syst->Draw("APe3");
    h_syst->Write(Name+"Syst");
    
    h_syst->Delete();
    histo_shifted_low->Delete();
    histo_shifted_high->Delete();
}

void PropagateCalibrationError(TString fileName){
    
    TString NewfileName = fileName;
    TPMERegexp("\\.root").Substitute(NewfileName,"_ErrCalib.root");
    
    TString cmd = "cp "+fileName+" "+NewfileName;
    std::cout << "cmd exit status = " << gSystem->Exec(cmd) << std::endl;
    
    TFile* NewFile = new TFile(NewfileName, "UPDATE");
    
    NewFile->mkdir("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis_High");
    NewFile->mkdir("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis_Low");
    
    NewFile->mkdir("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis_High");
    NewFile->mkdir("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis_Low");
    
    NewFile->mkdir("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis_High");
    NewFile->mkdir("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis_Low");

    std::vector<TFitResultPtr> Fit_results = PlotEnergyCalibration();
    
    //COV1
    TH1D* Simu_histoCOV10_10MeV = (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/", "nEvis_triggered_Vol_Bot_Ge_0-10MeV", Simu_histoCOV10_10MeV, Fit_results[0]);
    
    TH1D* Simu_histoCOV10_100MeV= (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-10MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/", "nEvis_triggered_Vol_Bot_Ge_0-100MeV", Simu_histoCOV10_100MeV, Fit_results[0]);
    
    TH1D* Simu_histoCOV10_1MeV= (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/Evis/nEvis_triggered_Vol_Bot_Ge_0-1MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Bot_Ge/", "nEvis_triggered_Vol_Bot_Ge_0-1MeV", Simu_histoCOV10_1MeV, Fit_results[0]);
    
    //COV2
    TH1D* Simu_histoCOV20_10MeV = (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-10MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Top_Ge/", "nEvis_triggered_Vol_Top_Ge_0-10MeV", Simu_histoCOV20_10MeV, Fit_results[1]);

    TH1D* Simu_histoCOV20_100MeV= (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-100MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Top_Ge/", "nEvis_triggered_Vol_Top_Ge_0-100MeV", Simu_histoCOV20_100MeV, Fit_results[1]);

    TH1D* Simu_histoCOV20_1MeV= (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Top_Ge/Evis/nEvis_triggered_Vol_Top_Ge_0-1MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Top_Ge/", "nEvis_triggered_Vol_Top_Ge_0-1MeV", Simu_histoCOV20_1MeV, Fit_results[1]);

    //LWO
    TH1D* Simu_histo0_10MeV = (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-10MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/", "nEvis_triggered_Vol_Mid_LWO_0-10MeV", Simu_histo0_10MeV, Fit_results[2]);

    TH1D* Simu_histo0_100MeV= (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-100MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/", "nEvis_triggered_Vol_Mid_LWO_0-100MeV", Simu_histo0_100MeV, Fit_results[2]);

    TH1D* Simu_histo0_1MeV= (TH1D*) NewFile->Get("OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/Evis/nEvis_triggered_Vol_Mid_LWO_0-1MeV");
    Shift(NewFile, "OrsayOVAlgoDir/Triggered/Vol_Mid_LWO/", "nEvis_triggered_Vol_Mid_LWO_0-1MeV", Simu_histo0_1MeV, Fit_results[2]);
    
    NewFile->Close();
}

double Chi2(TF1* Model, TGraphErrors* Values){
    double chi2 =0;
    double Val, errVal, Mod, E;
    for (int i=0; i<Values->GetN(); i++){
        Val=Values->GetPointY(i);
        E=Values->GetPointX(i);
        errVal=Values->GetErrorY(i);
        Mod = Model->Eval(E);
        chi2 += (Val-Mod)*(Val-Mod)/(errVal*errVal);
    }
    return chi2;
}

void Concatenate(Double_t *s1,Double_t *s2, Double_t *s, int len1, int len2)
{
    int j;
    ///Define k to store the values on Kth address Kstart from 0 to len1+len2;
    int k=0;

        for(j=0;j<len1;j++)
        {
            s[k]=s1[j];
            k++;

        }
        for(j=0;j<len2;j++)
        {
            s[k]=s2[j];
            k++;
        }
}

void PlotEnergyResolution(){

    //NTD
//    Double_t Energy[3] = {1461, 2615, 5233};
//    Double_t FWHM[3] = {53.73/1461.*100, 54.47/2615*100, 55.49/5233*100};
//    Double_t Err_Energy[3] = {4, 6, 6};
//    Double_t Err_FWHM[3] = {3.16/1461*10, 8.23/2615*100, 6.86/5233*100};
//
//    //2021_05_04_14_47_35
//    Double_t Energy_shielding[3] = {511, 1461, 2615};
//    Double_t FWHM_shielding[3] = {12.63/511*100, 11.97/1461*100, 17.1/2615*100};
//    Double_t Err_Energy_shielding[3] = {0, 0, 0};
//    Double_t Err_FWHM_shielding[3] = {3.16/511*100, 1.01/1461*100, 4.9/2615*100};
//
//    //2021_05_12_17_27_27 Eddy
    const int Npoints_NTD = 7;
    Double_t Energy[Npoints_NTD] = {511, 583.19, 911.2, 968.97, 1460.83, 2614.53, 5110.49};
    Double_t FWHM[Npoints_NTD] = {2.01741, 1.97213, 1.18474, 0.844022, 0.80022, 0.514869, 0.310146};
    Double_t Err_Energy[Npoints_NTD]={0., 0., 0., 0., 0., 0., 0.};
    Double_t Err_FWHM[Npoints_NTD]={0.209072, 0.481576, 0.243623, 0.241726, 0.0552832, 0.0722689, 0.00763136};
////
////    //COV1
//    const int Npoints_COV1 = 4;
//    Double_t EnergyCOV1[4] = {59.95, 609, 1461, 2615};
//    Double_t FWHMCOV1[4] = {1.545/59.95*100, 3.1/609*100, 7.7/1461*100, 14.4/2615*100};
//    Double_t Err_EnergyCOV1[4] = {0.01, 0.1, 0.01, 1};
//    Double_t Err_FWHMCOV1[4] = {0.005/59.95*100, 0.5/609*100, 0.5/1461*100, 1./2615*100};
//
//    //2021_05_04_14_47_35
//    const int Npoints_COV1 = 4;
//    Double_t EnergyCOV1_shielding[4] = {59.95, 511, 1461, 2615};
//    Double_t FWHMCOV1_shielding[4] = {1.566/59.95*100, 3.9/511*100, 8.6/1461*100, 14./2615*100};
//    Double_t Err_EnergyCOV1_shielding[4] = {0.01, 0.1, 0.01, 1};
//    Double_t Err_FWHMCOV1_shielding[4] = {0.003/59.95*100, 0.3/511*100, 0.3/1461*100, 1./2615*100};

//    //2021_05_12_17_27_27 Eddy
//    const int Npoints_COV1 = 8;
//    Double_t EnergyCOV1[Npoints_COV1] = {59.5, 511, 583.19, 609.31, 911.2, 968.97, 1460.83, 2614.53};
//    Double_t FWHMCOV1[Npoints_COV1] = {2.80874, 0.759667, 0.730349, 0.577326, 0.638798, 0.612, 0.633324, 0.573677};
//    Double_t Err_EnergyCOV1[Npoints_COV1]={0., 0., 0., 0., 0., 0., 0., 0.};
//    Double_t Err_FWHMCOV1[Npoints_COV1]={0.0048997, 0.0437959, 0.0827713, 0.131065, 0.0956908, 0.131382, 0.0229679, 0.0242575};

//    //COV2
//    Double_t Energy[5] = {510, 609, 1461, 1765, 2615};
//    Double_t FWHM[5] = {2.7/510*100, 3.13/609*100, 8.3/1461*100, 8.8/1765*100, 11./2615*100};
//    Double_t Err_Energy[5] = {0.1, 0.5, 0.1, 1, 1};
//    Double_t Err_FWHM[5] = {1./510*100, 0.4/609*100, 0.5/1461*100, 1.4/1765*100, 2.5/2615*100};

//    //2021_05_12_17_27_27 Eddy
//    const int Npoints_COV2 = 7;
//    Double_t Energy[Npoints_COV2] = {511, 583.19, 609.31, 911.2, 968.97, 1460.83, 2614.53};
//    Double_t FWHM[Npoints_COV2] = { 0.765312, 0.558863, 0.678529, 0.544988, 0.801134, 0.507778, 0.446116};
//    Double_t Err_Energy[Npoints_COV2]={0., 0., 0., 0., 0., 0., 0.};
//    Double_t Err_FWHM[Npoints_COV2]={0.0411712, 0.0677204, 0.117033, 0.0577587, 0.0883557, 0.0205737, 0.0307217};
//
//    TGraphErrors* grCOV1 = new TGraphErrors(Npoints_COV1,EnergyCOV1,FWHMCOV1, Err_EnergyCOV1,Err_FWHMCOV1);
////    TGraphErrors* grCOV1_shielding = new TGraphErrors(4,EnergyCOV1_shielding,FWHMCOV1_shielding, Err_EnergyCOV1_shielding,Err_FWHMCOV1_shielding);
//    TGraphErrors* grCOV2 = new TGraphErrors(Npoints_COV2,Energy,FWHM, Err_Energy,Err_FWHM);
//
//    Double_t Energy_both[Npoints_COV1+Npoints_COV2]; Double_t FWHM_both[Npoints_COV1+Npoints_COV2]; Double_t Err_Energy_both[Npoints_COV1+Npoints_COV2]; Double_t Err_FWHM_both[Npoints_COV1+Npoints_COV2];
//    Concatenate(EnergyCOV1,Energy, Energy_both, Npoints_COV1, Npoints_COV2);
//    Concatenate(FWHMCOV1,FWHM, FWHM_both, Npoints_COV1, Npoints_COV2);
//    Concatenate(Err_EnergyCOV1,Err_Energy, Err_Energy_both, Npoints_COV1, Npoints_COV2);
//    Concatenate(Err_FWHMCOV1,Err_FWHM, Err_FWHM_both, Npoints_COV1, Npoints_COV2);
////
//    TGraphErrors* gr_both = new TGraphErrors(Npoints_COV1+Npoints_COV2,Energy_both, FWHM_both, Err_Energy_both, Err_FWHM_both);
    
    TGraphErrors* grNTD = new TGraphErrors(Npoints_NTD,Energy,FWHM, Err_Energy,Err_FWHM);
//    TGraphErrors* grNTD_shielding = new TGraphErrors(3,Energy_shielding,FWHM_shielding, Err_Energy_shielding,Err_FWHM_shielding);

//    grCOV1->SetTitle("Energy resolution COV1");
//    grCOV1->SetMarkerStyle(21);
//    grCOV1->SetMarkerColor(4);
//    grCOV1->SetMarkerSize(0.6);
//    grCOV1->GetXaxis()->SetTitle("True energy [keV]");
//    grCOV1->GetXaxis()->SetTitleOffset(1.2);
//    //gr->GetXaxis()->SetRangeUser(0,1500);
//    grCOV1->GetYaxis()->SetTitle("#sigma/E [%]");
//    grCOV1->GetYaxis()->SetTitleOffset(1.4);
//
//    //gr->GetYaxis()->SetRangeUser(0.6,2);
//
//    grCOV1->SetMarkerStyle(21);
//    grCOV1->SetMarkerSize(1);
//
    TCanvas *c1 = new TCanvas("c1","Energy resolution",1000,600);
//
//    grCOV2->SetTitle("Energy resolution COV2");
//    grCOV2->SetMarkerStyle(21);
//    grCOV2->SetMarkerColor(3);
//    grCOV2->SetMarkerSize(0.6);
//    grCOV2->GetXaxis()->SetTitle("True energy [keV]");
//    grCOV2->GetXaxis()->SetTitleOffset(1.2);
//    grCOV2->GetXaxis()->SetLimits(0,3000);
//    grCOV2->GetYaxis()->SetLimits(0,0.025);
//    grCOV2->GetYaxis()->SetTitle("#sigma/E [%]");
//    grCOV2->GetYaxis()->SetTitleOffset(1.4);
//
//    grCOV1->SetMarkerStyle(21);
//    grCOV1->SetMarkerSize(1);
//    grCOV2->Draw("APSAME");
//    grCOV1->Draw("PSAME");
//    grCOV1_shielding->Draw("PSAME");
//    grCOV1_shielding->SetMarkerStyle(22);
//    grCOV1_shielding->SetMarkerSize(1);
//    grCOV1_shielding->SetMarkerColor(kRed);

    gPad->SetGridx();
    gPad->SetGridy();


    grNTD->SetTitle("Energy resolution NTD");
    grNTD->SetMarkerStyle(21);
    grNTD->SetMarkerColor(3);
    grNTD->SetMarkerSize(0.6);
    grNTD->GetXaxis()->SetTitle("True energy [keV]");
    grNTD->GetXaxis()->SetTitleOffset(1.2);
    grNTD->GetXaxis()->SetLimits(0,6000);
    grNTD->GetYaxis()->SetLimits(0,0.025);
    grNTD->GetYaxis()->SetTitle("#sigma/E [%]");
    grNTD->GetYaxis()->SetTitleOffset(1.4);

    grNTD->SetMarkerStyle(20);
    grNTD->SetMarkerSize(1.5);
    grNTD->Draw("APSAME");
//
//    grNTD_shielding->SetMarkerStyle(21);
//    grNTD_shielding->SetMarkerSize(1.5);
//    grNTD_shielding->SetMarkerColor(kRed);
//    grNTD_shielding->Draw("PSAME");
//
    //Fit
//    TF1* fit = new TF1("Peak fit","pol2",100,1500);
//    fit->SetParNames("A","B","C");
//
    TF1* fit = new TF1("Rac fit","sqrt([0] + [1]*x)/x",0,6000);
    fit->SetParNames("#sigma_{0}^{2}","p1");
    
//    TF1* fit_2 = new TF1("Knoll fit","sqrt([0]/x + [1]/(x*x)+[2])",0,3000);
//    fit_2->SetParNames("A","B", "C");

    //Fit parameter initialization
    Double_t par[2];
    par[0] = 680000;
    par[1] = 0.3;
    fit->SetParameters(par);
    fit->SetParLimits(1,0,10000);
//    fit->FixParameter(1, 0);

//    Double_t par2[3];
//    par2[0] = 0.01*100;
//    par2[1] = 66*100*100;
//    par2[2] = 0.000001*100;
//    fit_2->SetParameters(par2);
//    fit_2->SetParLimits(0,0,10000);
//    fit_2->SetParLimits(1,0,1000000);
////
//    gStyle->SetOptStat(0);
//    TFitResultPtr rKnoll = gr_both->Fit("Knoll fit","RS+");
//    gStyle->SetOptFit(111);
//    fit_2->SetLineColor(kBlue);
//    fit_2->DrawClone("LSAME");
//
//    cout<<"chi2/Ndf Both = "<< Chi2(fit_2, gr_both)<<"/"<< Npoints_COV1+Npoints_COV2-3<<"="<<Chi2(fit_2, gr_both)/(Npoints_COV1+Npoints_COV2-3)<<endl;
//    cout<<"chi2/Ndf COV1 = "<< Chi2(fit_2, grCOV1)<<"/"<< Npoints_COV1-3<<"="<<Chi2(fit_2, grCOV1)/(Npoints_COV1-3)<<endl;
////    cout<<"chi2 COV1 shielding with COV1 fit = "<< Chi2(fit_2, grCOV1_shielding)<<endl;
//
//////    rKnoll = grCOV1_shielding->Fit("Knoll fit","RS+");
//////    fit_2->Draw("LSAME");
//////    fit_2->SetLineColor(kRed);
//////
//////    cout<<"chi2 COV1 shielding = "<< Chi2(fit_2, grCOV1_shielding)<<endl;
//////    cout<<"chi2 COV1 with COV1 shielding fit = "<< Chi2(fit_2, grCOV1)<<endl;
////    cout<<"chi2 COV2 = "<< Chi2(fit_2, grCOV2)<<endl;
////
////    gStyle->SetOptStat(0);
////    TFitResultPtr rKnollCOV2 = grCOV2->Fit("Knoll fit","RS+");
////    gStyle->SetOptFit(111);
////    fit_2->SetLineColor(kRed);
////    fit_2->Draw("LSAME");
//    cout<<"chi2/Ndf COV2 = "<< Chi2(fit_2, grCOV2)<<"/"<< Npoints_COV2-3<<"="<<Chi2(fit_2, grCOV2)/(Npoints_COV2-3)<<endl;
////    cout<<"chi2 COV2 = "<< Chi2(fit_2, grCOV2)<<endl;
////
    TFitResultPtr rNTD = grNTD->Fit("Rac fit","RS0");
    fit->SetLineColor(3);
    fit->DrawClone("LSAME");
    cout<<"chi2/Ndf NTD = "<< Chi2(fit, grNTD)<<"/"<< Npoints_NTD-2<<"="<<Chi2(fit, grNTD)/(Npoints_NTD-2)<<endl;

//    rNTD = grNTD_shielding->Fit("Rac fit","RS0");
//    fit->Draw("LSAME");
//    fit->SetLineColor(kRed);

}

void PlotEffiency(){

    //Calibration points from Eu-152 source
    Double_t Energy[12] = {121.78,244.7,295.94,344.28,443.96,778.9,867.38,964.06,1085.84,1212.95,1299.14,1408.01};

    Double_t PeakArea[12] = {74.72,13.25,0.71,34.48,3.23,8.13,2.40,7.69,4.70,0.61,0.67,7.96};
    Double_t ErrPeakArea[12] = {0.23,0.09,0.03,0.20,0.05,0.1,0.06,0.08,0.07,0.02,0.02,0.05};

    Double_t Intensity[12] = {28.53,7.55,0.44,25.59,2.83,12.93,4.23,14.51,10.11,1.41,1.63,20.87};

    Double_t Eff[12];
    Double_t ErrEff[12];

    //Computation of "efficiency"
    for (Int_t i=0; i<12; i++){
        Eff[i] = PeakArea[i]/(Intensity[i]/100);
        ErrEff[i] = ErrPeakArea[i]/(Intensity[i]/100);
    }

    Double_t Norm = Eff[0];
    //Normalisation of "efficiency"
    for (Int_t i=0; i<12; i++){
        Eff[i] = Eff[i]/Norm*100;
        ErrEff[i] = ErrEff[i]/Norm*100;
    }

    TGraph* gr = new TGraphErrors(12,Energy,Eff,0,ErrEff);
    gr->SetTitle("Relative efficiency");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(0.6);
    gr->GetXaxis()->SetTitle("Energy [keV]");
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetXaxis()->SetRangeUser(0,1500);
    gr->GetYaxis()->SetTitle("\\epsilon");
    gr->GetYaxis()->SetTitleOffset(1.4);
    //gr->GetYaxis()->SetRangeUser(0.6,2);

    TCanvas *c1 = new TCanvas("c1","Efficiency",200,10,600,500);
    gr->Draw("AP");

    gPad->SetGridx();
    gPad->SetGridy();

//    TF1* fit = new TF1("Peak fit","expo",300,1450);
//    //fit->SetParNames("A","B");
//
//    //Fit parameter initialization
//    Double_t par[2];
//    par[0] = 0.5;
//    par[1] = 0.1;
//    fit->SetParameters(par);
//    //fit->SetParLimits(0,-2,2);
//    //fit->SetParLimits(1,-1,1);
//
//    gStyle->SetOptStat(0);
//    TFitResultPtr r = gr->Fit("Peak fit","RS0");
//    gStyle->SetOptFit(111);
//    fit->Draw("LSAME");
//    fit->SetLineWidth(2);
//    fit->SetLineColor(1);
//    fit->SetLineStyle(9);
}

void Fit_Bremsstrahlung(TString InputFileName, UInt_t Npeaks=10, Double_t Thres = 0.005){
    TFile* RunFile = new TFile(InputFileName,"READ");
    TH1D* h = (TH1D*)RunFile->Get("Cal_spectrum");
    
    h->Draw();
    
    TF1* Fit_func = new TF1("Fit_func", "[0]*TMath::Landau(x, [1], [2])", 50, 600);
    Fit_func->SetParameter(0, h->GetMaximum());
    Fit_func->SetParameter(1, 100);
    Fit_func->SetParameter(2, 50);
    Fit_func->SetParName(0, "Amplitude");
    Fit_func->SetParName(1, "mpv");
    Fit_func->SetParName(2, "sigma");
    
    //Gamma_func->Draw("SAME");
    h->Fit(Fit_func, "RN+");
    Fit_func->Draw("SAME");
    
    TH1D* new_h = (TH1D*) h->Clone();
    for (int i = 1; i<=h->GetNbinsX(); i++){
        if (h->GetBinCenter(i)<50){
            new_h->SetBinContent(i,0);
            new_h->SetBinError(i,0);
        }
        else {new_h->SetBinContent(i, h->GetBinContent(i)-Fit_func->Eval(h->GetBinCenter(i)));}
    }
    TCanvas* c2=new TCanvas();
    new_h->Draw();
    
    //Use TSpectrum to find the peak candidates
    TSpectrum *s = new TSpectrum(2*Npeaks+5);
    const Int_t Nfound = s->Search(new_h,15,"new",Thres);
    cout << "Found " << Nfound << " candidate peaks to fit\n";
    //hCal->Draw();
    //hCal->DrawCopy("E0 HIST SAME");

    
    //Results of TSpectrum
    Double_t *Xpeaks = s->GetPositionX();
    Double_t *Ypeaks = s->GetPositionY();
    sort(Xpeaks,Xpeaks+Nfound);
    
    std::vector<std::pair<Double_t,Double_t>> peakData;
    
    //Gross estimate of count rate in peaks + filling std::vector
    for (Int_t p=0;p<Nfound;p++) {
        Double_t xp = Xpeaks[p];
        Double_t yp = h->GetBinContent(h->GetXaxis()->FindBin(xp));
        peakData.push_back(std::make_pair(xp,yp));
        //Float_t xp = xpeaks[Selected[p]];
        // if(p>nfound) break;
        Double_t Signal=0;
        Double_t Background=0;
        Double_t Total=0;
        Total = h->Integral(h->FindBin(xp-10),h->FindBin(xp+10));
        Background = 2 * ( h->Integral(h->FindBin(xp-10),h->FindBin(xp-6)) + h->Integral(h->FindBin(xp+6),h->FindBin(xp+10)) );
        Signal = Total - Background;
        
        cout<<"Peak #" << p << " location: "<< xp << " - rough integral: " << Signal << " Bq" << endl;
    }
}
