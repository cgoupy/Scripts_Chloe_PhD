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

std::vector<std::pair<Double_t,Double_t>> PeakFinder(TString InputFileName, UInt_t Npeaks=15, Double_t Thres = 0.005)
{
    // test input file
    struct stat buf; // to test input file
    if (stat(InputFileName,&buf) == -1) { // test whether the file exist or not, and get its size
        cerr << "***Error: file missing: " << InputFileName << "\nExiting...\n";
        exit(-1);
    }
    
    TFile* RunFile = new TFile(InputFileName,"READ");
    RunFile->ls();
    
    TPMERegexp("\\.root").Substitute(InputFileName,"_analyse.txt");
    ofstream Resultat(InputFileName, ios::out | ios::trunc);

    //Retrieve histogram of calibrated energies
    //TH1D *hCal= (TH1D*) RunFile->Get("cal_spectrum");
    TH1D *hCal= (TH1D*) RunFile->Get("COV_r21_c4c8_Singles");

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
        Total = hCal->Integral(hCal->FindBin(xp-10),hCal->FindBin(xp+10));
        Background = 2 * ( hCal->Integral(hCal->FindBin(xp-10),hCal->FindBin(xp-6)) + hCal->Integral(hCal->FindBin(xp+6),hCal->FindBin(xp+10)) );
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


void PeakFitter(TString FileName = "Calibration-Eu152.root", Int_t PeakNumber=14){
    
    
    //Mass of the HPGe detector
    //Double_t M = (5.5*TMath::Pi()*pow(5.55/2,2)-4.4*TMath::Pi()*pow(1.15/2,2))*5.323/1000;
    
    std::vector<std::pair<Double_t,Double_t>> Peaks = PeakFinder(FileName,10, 0.001);
    
    const Int_t NPeaks = Peaks.size();
    TH1D* h = (TH1D*) gDirectory->Get("COV_r21_c4c8_Singles");
    
    //Draw full spectrum
    TH1D* hh = (TH1D*) h->Clone();
    //hh->SetName("Gamma spectrum");
    hh->SetTitle("COV");
    //hh->SetTitle("Co^{60}");
    hh->SetStats(0);
    
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
    
    TF1* fit_pol = new TF1("Pol fit","pol1(0)",Xp-100,Xp-50);
    fit_pol->SetParNames("p_{0}","p_{1}");
    
    TFitResultPtr r_pol = h->Fit("Pol fit","RSQ");
    
    TF1* fit_gauss = new TF1("Gauss fit","gaus(0)",0.98*Xp,1.02*Xp);
    fit_gauss->SetParNames("A","E_{rec}","\\sigma");
    
    TFitResultPtr r_gauss = h->Fit("Gauss fit","RS");
    
    TF1* fit = new TF1("Peak fit","gaus(0)+pol1(3)",0.95*Xp,1.05*Xp);
    fit->SetParNames("A","E_{rec}","\\sigma","p_{0}","p_{1}");
    
    //Fit parameter initialization
    Double_t par[5];
    par[0] = fit_gauss->GetParameter(0);
    par[1] = fit_gauss->GetParameter(1);
    par[2] = fit_gauss->GetParameter(2);
    
    par[3] = fit_pol->GetParameter(0);
    par[4] = fit_pol->GetParameter(1);
    
    fit->SetParameters(par);
    fit->SetParLimits(0,0.,1.2*par[0]);
    fit->SetParLimits(1,0.9*par[1],1.1*par[1]);
    fit->SetParLimits(2,0.,1.2*par[2]);
    fit->SetParLimits(3,0.8*par[3],1.2*par[3]);
    fit->SetParLimits(4,0.8*par[4],1.2*par[4]);
    fit->SetNpx(2000);
    
    //Draw fit
    TCanvas *CanFit = new TCanvas("CanFit","Canvas of peak fit",0,0,1000,900);
    CanFit->cd();
    h->Draw("E0 HIST");
    h->SetTitle(Form("Peak - %g keV", Xp));
    gStyle->SetOptStat(0);
    TFitResultPtr r = h->Fit("Peak fit","RS");
    gStyle->SetOptFit(111);
    
    fit_pol->SetLineColor(kGreen);
    fit_pol->Draw("LSAME");
    
    fit_gauss->SetLineColor(kGreen);
    fit_gauss->Draw("LSAME");
    
    fit->Draw("LSAME");
    fit->SetLineWidth(2);
    h->GetXaxis()->SetRangeUser(Xp-200,Xp+200);
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
    
    TF1* PeakFunc = new TF1("Peak function","gaus",Xp-100,Xp+100);
    PeakFunc->SetParameter(0,fit->GetParameter(0));
    PeakFunc->SetParameter(1,fit->GetParameter(1));
    PeakFunc->SetParameter(2,fit->GetParameter(2));
    
    Double_t A = PeakFunc->Integral(fit->GetParameter(1)-100*fit->GetParameter(2),fit->GetParameter(1)+100*fit->GetParameter(2));
    Double_t Err = PeakFunc->IntegralError(fit->GetParameter(1)-100*fit->GetParameter(2),fit->GetParameter(1)+100*fit->GetParameter(2),PeakFunc->GetParameters(),covPeak.GetMatrixArray());
    
    std::cout << "Count rate for peak #" << PeakNumber << ": " << A << " +/- " << Err << " Bq" << std::endl;
    
}

void PlotEnergyCalibration(){
    
    //Calibration points from Eu-152 source
    //Double_t Etrue[14] = {39.52, 121.782, 244.697, 295.697, 344.2785, 411.12, 443.961, 778.90, 867.38, 964.06, 1085.84, 1112.076, 1299.14, 1408.01};
    //Double_t Erec[14] = {39.86, 121.85, 244.71, 295.92, 344.44, 411.14, 443.94, 779.05, 867.55, 964.35, 1086.09, 1112.37, 1299.46, 1408.4};
    //Double_t ChanNum[14]={177, 543, 1090, 1317, 1533, 1830, 1978, 3468, 3862, 4293, 4835, 4952, 5785, 6270};
    
    //Calibration points
    Double_t Etrue[3] = {1173.01, 1332.03, 1460.28};
    Double_t Erec[3] = {1173.228, 1332.492 ,1460.82};
    Double_t ChanNum[3]={5222.5, 5930.5, 6501.5};
    
    TGraph* gr = new TGraph(3,Etrue,ChanNum);
    //gr->SetTitle("Energy scale linearity from Co60");
    gr->SetTitle("Calibration curve from Co60");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(0.6);
    gr->GetXaxis()->SetTitle("True energy [keV]");
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetXaxis()->SetRangeUser(0,1500);
    gr->GetYaxis()->SetTitle("Reconstructed energy [keV]");
    //gr->GetYaxis()->SetTitle("DESPEC bin [a.u.]");
    gr->GetYaxis()->SetTitleOffset(1.4);
    //gr->GetYaxis()->SetRangeUser(0,1500);
    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1);
    
    
    TCanvas *c1 = new TCanvas("c1","Energy linearity",200,10,600,500);
    gr->Draw("AP");
    
    gPad->SetGridx();
    gPad->SetGridy();
    
    
    //Fit
    TF1* fit = new TF1("Peak fit","pol1(0)",0,1500);
    fit->SetParNames("A","B");
    
    //Fit parameter initialization
    Double_t par[2];
    par[0] = 0;
    par[1] = 1;
    fit->SetParameters(par);
    fit->SetParLimits(0,-10,10);
    fit->SetParLimits(1,0.9,1.1);
    
    gStyle->SetOptStat(0);
    TFitResultPtr r = gr->Fit("Peak fit","RS0");
    gStyle->SetOptFit(111);
    fit->Draw("LSAME");
    fit->SetLineWidth(2);
    fit->SetLineColor(1);
    fit->SetLineStyle(9);
    
    
}

void PlotEnergyResolution(){
    
    //Calibration points COV run 19
//    Double_t Energy[5] = {603.7, 1109.75, 1455.44, 1751.05, 2194.8};
//    Double_t FWHM[5] = {22.01, 9.097, 14.8, 19.94, 37.56};
//    Double_t Err_Energy[5] = {3.2, 1.7, 0.9, 1.7, 1.76};
//    Double_t Err_FWHM[5] = {12.06, 1.496, 1.1, 1.84, 3.11};
    
    //Calibration points COV run 19
    Double_t Energy[5] = {509.3, 1122.8, 1461.7, 1762.4, 2617.5};
    Double_t FWHM[5] = {6.78, 10.7, 11.9, 18.14, 17.0};
    Double_t Err_Energy[5] = {0.51, 2.011, 1.03, 2.5, 0.72};
    Double_t Err_FWHM[5] = {0.7, 1.98, 0.89, 2.74, 0.73};
    
    //Calibration points from Co-60 source
//    Double_t Energy[3] = {1173.01, 1332.03, 1460.28};
//    Double_t FWHM[3] = {0.7056, 0.7476, 0.8101};
//    Double_t Err_Energy[3] = {0.0, 0.0, 0.1};
//    Double_t Err_FWHM[3] = {0.0096, 0.0107, 0.1257};
    
    TGraphErrors* gr = new TGraphErrors(5,Energy,FWHM, Err_Energy,Err_FWHM);
    gr->SetTitle("Energy resolution from Eu152");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(0.6);
    gr->GetXaxis()->SetTitle("True energy [keV]");
    gr->GetXaxis()->SetTitleOffset(1.2);
    //gr->GetXaxis()->SetRangeUser(0,1500);
    gr->GetYaxis()->SetTitle("FWHM [keV]");
    gr->GetYaxis()->SetTitleOffset(1.4);
    //gr->GetYaxis()->SetRangeUser(0.6,2);
    
    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1);
    
    TCanvas *c1 = new TCanvas("c1","Energy resolution",200,10,600,500);
    gr->Draw("AP");
    
    gPad->SetGridx();
    gPad->SetGridy();
    
    
    //Fit
//    TF1* fit = new TF1("Peak fit","pol2",100,1500);
//    fit->SetParNames("A","B","C");

//    //Fit parameter initialization
//    Double_t par[3];
//    par[0] = 0.5;
//    par[1] = 0.1;
//    par[2] = 0.1;
//    fit->SetParameters(par);
//    fit->SetParLimits(0,0.001,1);
//    fit->SetParLimits(1,0.001,1);
//    fit->SetParLimits(2,0.001,1);
    
    TF1* fit = new TF1("Peak fit","pol1",0,1500);
    fit->SetParNames("A","B");

    //Fit parameter initialization
    Double_t par[2];
    par[0] = 0.5;
    par[1] = 0.1;
    fit->SetParameters(par);
    //fit->SetParLimits(0,-2,2);
    //fit->SetParLimits(1,-1,1);

    gStyle->SetOptStat(0);
    TFitResultPtr r = gr->Fit("Peak fit","RS0");
    gStyle->SetOptFit(111);
    fit->Draw("LSAME");
    fit->SetLineWidth(2);
    fit->SetLineColor(1);
    fit->SetLineStyle(9);

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
