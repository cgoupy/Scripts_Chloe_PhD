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
#include <Math/GSLRndmEngines.h>
#include "MacrosRudolph_Covariance/Utilities_Statistics.cpp"

TString Spec_name = "h1_LWO_Amplitude";

using namespace TMath;

//Derivative with respect to parameter number k for a polynomial function
double k_der_pol(double x, TVectorD coeffs, int k){
    int Dim = coeffs.GetNrows();
    double der=0;
    for (int i=0; i<Dim; i++){
        if (i!=k){
            der+=coeffs(i)*pow(x, i);
        }
        else {
            der+=pow(x, i);
        }
    }
    return der;
}

//1 sigma for polynomial fit
double sigma1(double x, TVectorD coeffs, TMatrixTSym<double> CovMat){
    int Dim = CovMat.GetNrows();
    double sig1 = 0;
    for (int i=0; i<Dim; i++){
        for (int j=0; j<Dim; j++){
            sig1 += CovMat(i, j)*k_der_pol(x, coeffs, i)*k_der_pol(x, coeffs, j);
        }
    }
    return sqrt(sig1);
}

class classsigma1 {
public:
    classsigma1(TVectorD coeffs, TMatrixTSym<double> CovMat): Dim(CovMat.GetNrows()), coeffs(coeffs), CovMat(CovMat){
//        cout<<"New sigma1 defined with cov matrix = "<<endl;
//        CovMat.Print();
    };
    
    double Evaluate(double* x, double* p={}){
        double x_val = x[0];
        return sigma1(x_val, coeffs, CovMat);
    };
    
    double EvaluatePlusF(double* x, double* p){
        double x_val = x[0];
        double b = p[0];
        double a = p[1];
        double func_value = b+a*x_val;
        return func_value+sigma1(x_val, coeffs, CovMat);
    };
    
    double EvaluateMinusF(double* x, double* p){
        double x_val = x[0];
        double b = p[0];
        double a = p[1];
        double func_value = b+a*x_val;
        return func_value-sigma1(x_val, coeffs, CovMat);
    };
    
    const int Dim;
    TVectorD coeffs = TVectorD(Dim);
    TMatrixTSym<double> CovMat = TMatrixTSym<double>(Dim);
};

//---- Vector of parmeters from fit
TVectorD Parameters_from_fit(TFitResultPtr fitresult){
    int N_par = fitresult->NPar();
    TVectorD Parameters(N_par);
    for (int i = 0; i<N_par; i++){
        Parameters(i)=fitresult->Parameter(i);
    }
    return Parameters;
}


std::vector<TFitResultPtr> PlotEnergyCalibration(Bool_t plot = false){
    
    int Emax = 30000;
    int shift_plot_Top=5;
    int shift_plot_LWO=10;
    std::vector<TFitResultPtr> FitResults;
    
    //---------------------------------- Calibration points COV Bot
    
    //2021_04_29_09_54_36
//    const int Npoints = 6;
//    Double_t Etrue[Npoints] = {0., 59.54, 511, 1460.82, 2614.511, 15690};
//    Double_t EtrueErr[Npoints] = {0., 0, 0, 0, 0, 10};
//    Double_t ChanNum[Npoints]={0., 0.05703, 0.4866, 1.39, 2.489, 15.49};
//    Double_t ChanNumErr[Npoints]={0.0015, 0.001, 0.001, 0.001, 4*0.001, 0.09};
//    //Double_t ChanNumErr[Npoints]={10*0.0015, 0.001, 4*0.001, 10*0.001, 10*0.001, 2*0.09};
//    
//    //2021_04_30_19_46_09
//    const int Npoints = 6;
//    Double_t Etrue[Npoints] = {0., 59.54, 511, 1460.82, 2614.511, 15690};
//    Double_t EtrueErr[Npoints] = {0., 0, 0, 0, 0, 10};
//    Double_t ChanNum[Npoints]={0., 0.05711, 0.4869, 1.392, 2.494, 16.16};
//    Double_t ChanNumErr[Npoints]={0.00114, 0.001, 0.01, 0.001, 0.01, 0.4};
    
    //2021_05_04_14_47_35
//    const int Npoints = 6;
//    Double_t Etrue[Npoints] = {0., 59.54, 511, 1460.82, 2614.511, 15690};
//    Double_t EtrueErr[Npoints] = {0., 0, 0, 0, 0, 10};
//    Double_t ChanNum[Npoints]={0., 0.05697, 0.4857, 1.387, 2.489, 15.8};
//    Double_t ChanNumErr[Npoints]={0.00114, 0.0001, 0.001, 0.001, 0.01, 0.2};

//    //2021_05_12_17_27_27 Eddy
//    const int Npoints = 10;
//    Double_t Etrue[Npoints] = {0., 59.54, 511, 583.19, 609.31, 911.2, 968.97, 1460.82, 2614.511, 15690};
//    Double_t EtrueErr[Npoints] = {0., 0, 0, 0, 0, 10};
//    Double_t ChanNum[Npoints]={0., 0.0566599, 0.485805, 0.554988, 0.580273, 0.868017, 0.921819, 1.39078, 2.49256, 15.71};
//    Double_t ChanNumErr[Npoints]={0.00114, 3.40939e-06, 0.000205742, 0.000444611, 0.000891229, 0.000803777, 0.00102633, 0.000392414, 0.000845872, 2*0.09};

    //NUCLEUS COV 2022_05_23_18_50_11
    const int Npoints = 6;
    Double_t Etrue[Npoints] = {0., 511, 1460.82, 1764, 2614.511, 20004.};
    Double_t EtrueErr[Npoints] = {0., 0, 0, 0, 0, 50};
    Double_t ChanNum[Npoints]={0., 0.34, 1.069, 1.2818, 1.9227, 15.09};
    Double_t ChanNumErr[Npoints]={0.00124, 0.01, 0.02, 0.02, 0.02, 0.1};
    
    TGraphErrors* gr = new TGraphErrors(Npoints,Etrue,ChanNum, EtrueErr, ChanNumErr);
    
    gr->SetTitle("Calibration Bot Ge");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kRed+2);
    gr->SetMarkerSize(1.8);
    gr->GetXaxis()->SetTitle("True energy [keV]");
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetXaxis()->SetRangeUser(0,20000);
    gr->GetYaxis()->SetTitle("Amplitude [V]");
    gr->GetYaxis()->SetTitleOffset(1.4);
    gr->SetMarkerStyle(21);
    gr->SetMarkerSize(1.8);
    
    TCanvas *c1 = new TCanvas("c1","Energy linearity",1000,600);
    gPad->SetGridx();
    gPad->SetGridy();
    
    
    //Fit

//    Linear fit
    TF1* Bot_Ge_fit =  new TF1("Bot_Ge_fit","([0]+[1]*x)",0,Emax);
    Bot_Ge_fit->SetParNames("b [#V]","a [#V/keV]");

    //Fit parameter initialization
    Bot_Ge_fit->SetParameter(0, 0);
    Bot_Ge_fit->SetParameter(1, 1);


    TFitResultPtr rBotGe = gr->Fit("Bot_Ge_fit","RS0", "", 0, Emax);
    FitResults.push_back(rBotGe);
    cout<<"Chi2/Ndf = " << rBotGe->Chi2()<<"/"<<rBotGe->Ndf()<<" = "<<rBotGe->Chi2()/rBotGe->Ndf()<<endl;

    //-- Covariance matrix
    TMatrixTSym<double> Bot_Ge_Cov = rBotGe->GetCovarianceMatrix();
    Bot_Ge_Cov.Print();

    //-- Correlation matrix
    TMatrixTSym<double> Bot_Ge_Cor = rBotGe->GetCorrelationMatrix();
    TCanvas *cCor = new TCanvas("Bot Ge", "Bot Ge");
    TH2D* hBot_Ge_Cor = new TH2D(Bot_Ge_Cor);
    gStyle->SetOptStat(0);
    hBot_Ge_Cor->GetZaxis()->SetRangeUser(-1, 1);
    hBot_Ge_Cor->SetTitle("Bot Ge calibration correlation matrix");
    hBot_Ge_Cor->Draw("colz");
    c1->cd();

//    //pol2 fit
//    TF1* Bot_Ge_fit =  new TF1("Bot_Ge_fit","([0]+[1]*x+[2]*x*x)",0,Emax);
//    Bot_Ge_fit->SetParNames("b [V]","a [V/keV]", "c [V/keV^2]");
//
//    //Fit parameter initialization
//    Bot_Ge_fit->SetParameter(0, 0);
//    Bot_Ge_fit->SetParameter(1, 1);
//    Bot_Ge_fit->SetParameter(2, 1);
//
//    TFitResultPtr rBotGe = gr->Fit("Bot_Ge_fit","RS0", "", 0, Emax);
//    FitResults.push_back(rBotGe);
//    cout<<"Chi2/Ndf = " << rBotGe->Chi2()<<"/"<<rBotGe->Ndf()<<" = "<<rBotGe->Chi2()/rBotGe->Ndf()<<endl;
//
//    // -- Covariance matrix
//    TMatrixTSym<double> Bot_Ge_Cov = rBotGe->GetCovarianceMatrix();
//    Bot_Ge_Cov.Print();
//
//    //-- Covariance matrix
//    TMatrixTSym<double> Bot_Ge_Cor = rBotGe->GetCorrelationMatrix();
//    TCanvas *cCor = new TCanvas("Bot Ge", "Bot Ge");
//    TH2D* hBot_Ge_Cor = new TH2D(Bot_Ge_Cor);
//    gStyle->SetOptStat(0);
//    hBot_Ge_Cor->GetZaxis()->SetRangeUser(-1, 1);
//    hBot_Ge_Cor->SetTitle("Bot Ge calibration correlation matrix");
//    hBot_Ge_Cor->Draw("colz");
//    c1->cd();

    
    /*Create a histogram to hold the confidence intervals*/
    classsigma1* sigma1_Bot  = new classsigma1(Parameters_from_fit(rBotGe), Bot_Ge_Cov);
    TH1D *Bot_Ge_1sig = new TH1D("Bot_Ge_1sig",
       "Bot_Ge_1sig", Emax, 0, Emax);
    
    for (int i=0; i<Bot_Ge_1sig->GetNbinsX(); i++){
        double x = Bot_Ge_1sig->GetBinCenter(i+1);
        Bot_Ge_1sig->SetBinContent(i+1, Bot_Ge_fit->Eval(x));
        double x_vect[1]={x};
        Bot_Ge_1sig->SetBinError(i+1, 10*sigma1_Bot->Evaluate(x_vect));
    }
    
    Bot_Ge_1sig->SetMarkerSize(0);
    Bot_Ge_1sig->SetMarkerColor(0);
    Bot_Ge_1sig->SetLineColor(0);
    Bot_Ge_1sig->SetStats(false);
    Bot_Ge_1sig->SetFillColorAlpha(2,0.4);
    if (plot) Bot_Ge_1sig->Draw("e3 same");
    Bot_Ge_1sig->GetXaxis()->SetTitle("True energy [keV]");
    Bot_Ge_1sig->GetYaxis()->SetTitle("Amplitude [V]");
    Bot_Ge_1sig->SetTitle("10#times#sigma");
    
    Bot_Ge_fit->SetLineWidth(3);
    Bot_Ge_fit->SetLineColor(2);
    Bot_Ge_fit->SetLineStyle(9);
    if (plot) Bot_Ge_fit->Draw("Same");
    Bot_Ge_fit->SetTitle(Form("\n #splitline{#splitline{a = %0.1f +/- %0.1f #muV/keV    }{b = %0.1f +/- %0.1f #muV}}{c = %0.1fe-3 +/- %0.1fe-3 #muV/keV^{2}}", 1e6*Bot_Ge_fit->GetParameter(1), 1e6*Bot_Ge_fit->GetParError(1), 1e6*Bot_Ge_fit->GetParameter(0), 1e6*Bot_Ge_fit->GetParError(0), 1e9*Bot_Ge_fit->GetParameter(2), 1e9*Bot_Ge_fit->GetParError(2)));
    
    if (plot) gr->Draw("P same");
    
    auto legend1 = new TLegend();
    legend1->AddEntry(gr);
    legend1->AddEntry(Bot_Ge_fit);
    legend1->AddEntry(Bot_Ge_1sig);
    if (plot) legend1->Draw();
    
    
    // Residuals
    TCanvas *cRes = new TCanvas("Residuals", "Residuals", 200,10,1200,1200);
    gStyle->SetOptStat(0);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.02);
    cRes->cd();
    
    Double_t Res_value[Npoints]={0, 0, 0, 0, 0};
    Double_t Res_Errors[Npoints]={0, 0, 0, 0, 0};
    double Value_i, Diff_i, Error_i;
    
    for (int i = 0; i<Npoints; i++){
        Value_i = ChanNum[i];
        Res_value[i] = (ChanNum[i]-Bot_Ge_fit->Eval(Etrue[i]))/Bot_Ge_fit->Eval(Etrue[i])*100;
        Res_Errors[i] = ChanNumErr[i]/Bot_Ge_fit->Eval(Etrue[i])*100;
    }
    TGraphErrors* ResidualsBot_pc = new TGraphErrors(Npoints,Etrue,Res_value, EtrueErr, Res_Errors);
    
    gPad->SetLeftMargin(0.18);
    ResidualsBot_pc->SetTitle("");
    ResidualsBot_pc->SetMarkerStyle(21);
    ResidualsBot_pc->SetMarkerColor(kRed+2);
    ResidualsBot_pc->SetLineColor(kRed+2);
    ResidualsBot_pc->SetMarkerSize(1.8);
    ResidualsBot_pc->GetXaxis()->SetTitle("True energy [keV]");
    ResidualsBot_pc->GetXaxis()->SetTitleOffset(1.2);
    ResidualsBot_pc->GetXaxis()->SetRangeUser(-100,20000);
    
    double Min_axis = -2;
    double Max_axis = shift_plot_LWO+5;
    
    ResidualsBot_pc->GetYaxis()->SetRangeUser(Min_axis,Max_axis);
    ResidualsBot_pc->GetYaxis()->SetAxisColor(2);
    ResidualsBot_pc->GetYaxis()->SetLabelColor(2);
    ResidualsBot_pc->GetYaxis()->SetTitleOffset(2.5);
    ResidualsBot_pc->GetYaxis()->SetTitle("Residuals [%]");
    
    if (plot) ResidualsBot_pc->Draw("AP");
    
    TH1D *Res_Bot_1sig = new TH1D("Res_Bot_1sig",
       "Res_Bot_1sig", Emax, 0, Emax);
    double x;
    double x_vect[1];
    for (int i=0; i<Res_Bot_1sig->GetNbinsX(); i++){
        x = Res_Bot_1sig->GetBinCenter(i+1);
        Res_Bot_1sig->SetBinContent(i+1, 0);
        x_vect[0]={x};
        Res_Bot_1sig->SetBinError(i+1, (sigma1_Bot->Evaluate(x_vect)/Bot_Ge_fit->Eval(x))*100);
    }
    
    Res_Bot_1sig->SetMarkerSize(0);
    Res_Bot_1sig->SetMarkerColor(0);
    Res_Bot_1sig->SetLineColor(0);
    Res_Bot_1sig->SetStats(false);
    Res_Bot_1sig->SetFillColorAlpha(2,0.4);
    if (plot) Res_Bot_1sig->Draw("e3 same");
    Res_Bot_1sig->GetXaxis()->SetTitle("True energy [keV]");
    Res_Bot_1sig->GetYaxis()->SetTitle("Amplitude [V]");
    Res_Bot_1sig->SetTitle("1 #sigma");
    
    c1->cd();
    
    
//   //--------------------------------- Calibration points Top_Ge
////    const int Npoints_Top = 5;
//////    //2021_04_29_09_54_36
////    Double_t Etrue2[Npoints_Top] = {0, 511, 1460.82, 2614.511, 15700};
////    Double_t EtrueErr2[Npoints_Top] = {0, 0, 0, 0, 10};
////    Double_t ChanNum2[Npoints_Top]={0., 0.3318, 0.9491, 1.701, 10.26};
////    Double_t ChanNumErr2[Npoints_Top]={0.0002, 0.0003, 2*0.0001, 2*0.001, 2*0.04};
//
//    //2021_05_04_14_47_35
////    Double_t Etrue2[Npoints_Top] = {0, 511, 1460.82, 2614.511, 15700};
////    Double_t EtrueErr2[Npoints_Top] = {0, 0, 0, 0, 10};
////    Double_t ChanNum2[Npoints_Top]={0., 0.3319, 0.9496, 1.701, 10.30};
////    Double_t ChanNumErr2[Npoints_Top]={0.0002, 0.0002, 2*0.001, 2*0.001, 0.06};
//
//    //2021_05_12_17_27_27 Eddy
//    const int Npoints_Top = 9;
//    Double_t Etrue2[Npoints_Top] = {0., 511, 583.19, 609.31, 911.2, 968.97, 1460.82, 2614.511, 15700};
//    Double_t EtrueErr2[Npoints_Top] = {0., 0, 0, 0, 0, 10};
//    Double_t ChanNum2[Npoints_Top]={0., 0.333111, 0.37948, 0.396891, 0.59303, 0.630386, 0.951391, 1.70384, 10.19};
//    Double_t ChanNumErr2[Npoints_Top]={0.0002, 0.000123857, 0.000278143, 0.000464021, 0.000319696, 0.000424584, 0.000179156, 0.000505504, 0.02};
//
//    TGraphErrors* gr2 = new TGraphErrors(Npoints_Top,Etrue2,ChanNum2, EtrueErr2, ChanNumErr2);
//    //gr->SetTitle("Energy scale linearity from Co60");
//    gr2->SetTitle("Calibration Top Ge");
//    gr2->SetMarkerStyle(22);
//    gr2->SetMarkerColor(kBlue+2);
//    gr2->SetMarkerSize(0.6);
//    gr2->GetXaxis()->SetTitle("True energy [keV]");
//    gr2->GetXaxis()->SetTitleOffset(1.2);
//    gr2->GetXaxis()->SetRangeUser(0,20000);
//    gr2->GetYaxis()->SetTitle("Amplitude [V]");
//    //gr->GetYaxis()->SetTitle("DESPEC bin [a.u.]");
//    gr2->GetYaxis()->SetTitleOffset(1.4);
//    //gr->GetYaxis()->SetRangeUser(0,1500);
//    gr2->SetMarkerStyle(22);
//    gr2->SetMarkerSize(1.8);
//
//    if (plot) gr2->Draw("PSAME");
//
//    //Fit
//    TF1* Top_Gelin = new TF1("Top_Ge","([0]+[1]*x)",0,20000);
//    Top_Gelin->SetParNames("b [#V]","a [#V/keV]");
//    Double_t par[2];
//    par[0]=0;
//    par[1]=1;
//
//    //Fit parameter initialization
//    Top_Gelin->SetParameters(par);
//    //fit2->SetParLimits(0,-10,10);
//    //fit2->SetParLimits(1,0.9,1.1);
//
//    //gStyle->SetOptStat(0);
//    TFitResultPtr r2 = gr2->Fit("Top_Ge","RS0");
//    FitResults.push_back(r2);
//    TMatrixTSym<double> Top_Ge_rlin_Cov = r2->GetCovarianceMatrix();
//    cout<<"Chi2/Ndf = " << r2->Chi2()<<"/"<<r2->Ndf()<<" = "<<r2->Chi2()/r2->Ndf()<<endl;
//    Top_Ge_rlin_Cov.Print();
//    
//    TMatrixTSym<double> Top_Ge_rlin_Cor = r2->GetCorrelationMatrix();
//    TH2D* hTop_Ge_rlin_Cor = new TH2D(Top_Ge_rlin_Cor);
//    TCanvas *cCor2 = new TCanvas("Top Ge", "Top Ge");
//    gStyle->SetOptStat(0);
//    hTop_Ge_rlin_Cor->GetZaxis()->SetRangeUser(-1, 1);
//    hTop_Ge_rlin_Cor->SetTitle("Top Ge calibration correlation matrix");
//    if (plot) hTop_Ge_rlin_Cor->Draw("colz");
//    c1->cd();
//    
//    /*Create a histogram to hold the confidence intervals*/
//    TH1D *Top_Ge_1sig = new TH1D("Top_Ge_1sig",
//       "Top_Ge_1sig", Emax, 0, Emax);
//    classsigma1* sigma1_Top = new classsigma1(Parameters_from_fit(r2), Top_Ge_rlin_Cov);
//    
//    for (int i=0; i<Top_Ge_1sig->GetNbinsX(); i++){
//        double x = Top_Ge_1sig->GetBinCenter(i+1);
//        Top_Ge_1sig->SetBinContent(i+1, Top_Gelin->Eval(x));
//        double x_vect[1]={x};
//        Top_Ge_1sig->SetBinError(i+1, 10*sigma1_Top->Evaluate(x_vect));
//    }
//    
//    Top_Ge_1sig->SetMarkerSize(0);
//    Top_Ge_1sig->SetMarkerColor(0);
//    Top_Ge_1sig->SetLineColor(0);
//    Top_Ge_1sig->SetStats(false);
//    Top_Ge_1sig->SetFillColorAlpha(4,0.4);
//    if (plot) Top_Ge_1sig->Draw("e3 same");
//    Top_Ge_1sig->SetTitle("10#times#sigma");
//    
//    Top_Gelin->SetLineWidth(3);
//    Top_Gelin->SetLineColor(4);
//    Top_Gelin->SetLineStyle(9);
//    if (plot) Top_Gelin->Draw("Same");
//    Top_Gelin->SetTitle(Form("#splitline{a = %0.1f +/- %0.1f #muV/keV    }{b = %0.1f +/- %0.1f #muV}", 1e6*Top_Gelin->GetParameter(1), 1e6*Top_Gelin->GetParError(1), 1e6*Top_Gelin->GetParameter(0), 1e6*Top_Gelin->GetParError(0)));
//    
//    if (plot) gr2->Draw("P same");
//    auto legend2 = new TLegend();
//    legend2->AddEntry(gr2);
//    legend2->AddEntry(Top_Gelin);
//    legend2->AddEntry(Top_Ge_1sig);
//    if (plot) legend2->Draw("same");
//    
//    // Residuals
//    cRes->cd();
//    
//    Double_t Res_value_Top[Npoints_Top]={0, 0, 0, 0, 0};
//    Double_t Res_Errors_Top[Npoints_Top]={0, 0, 0, 0, 0};
//    
//    for (int i = 0; i<Npoints_Top; i++){
//        Value_i = ChanNum2[i];
//        Res_value_Top[i] = ((ChanNum2[i]-Top_Gelin->Eval(Etrue2[i]))/Top_Gelin->Eval(Etrue2[i]))*100+shift_plot_Top;
//        Res_Errors_Top[i] = (ChanNumErr2[i]/Top_Gelin->Eval(Etrue2[i]))*100;
//    }
//    TGraphErrors* ResidualsTop_pc = new TGraphErrors(Npoints_Top,Etrue2,Res_value_Top, EtrueErr2, Res_Errors_Top);
//
//    ResidualsTop_pc->SetTitle("Residuals Bot Ge");
//    ResidualsTop_pc->SetMarkerStyle(22);
//    ResidualsTop_pc->SetMarkerColor(kBlue+2);
//    ResidualsTop_pc->SetLineColor(kBlue+2);
//    ResidualsTop_pc->SetMarkerSize(1.8);
//    ResidualsTop_pc->GetXaxis()->SetTitle("True energy [keV]");
//    ResidualsTop_pc->GetXaxis()->SetRangeUser(-100,20000);
//    ResidualsTop_pc->GetYaxis()->SetTitle("Residuals [%]");
//    
//    int Ndiv = ResidualsBot_pc->GetYaxis()->GetNdivisions();
//    
//    TGaxis *axis = new TGaxis(-1100, Min_axis, -1100, Max_axis, Min_axis-shift_plot_Top, Max_axis-shift_plot_Top, Ndiv, "");
//    axis->SetLineColor(4);
//    axis->SetLabelColor(4);
//    axis->SetTextFont(72);
//    axis->SetTickLength(0.005);
//    axis->SetLabelFont(ResidualsTop_pc->GetXaxis()->GetLabelFont());
//    axis->SetLabelSize(ResidualsTop_pc->GetXaxis()->GetLabelSize());
//    if (plot) axis->Draw("SAME");
//    
//    if (plot) ResidualsTop_pc->Draw("PSAME");
//    
//    TH1D *Res_Top_1sig = new TH1D("Res_Top_1sig",
//       "Res_Top_1sig", Emax, 0, Emax);
//    for (int i=0; i<Res_Top_1sig->GetNbinsX(); i++){
//        double x = Res_Top_1sig->GetBinCenter(i+1);
//        Res_Top_1sig->SetBinContent(i+1, shift_plot_Top);
//        double x_vect[1]={x};
//        Res_Top_1sig->SetBinError(i+1, 100*sigma1_Top->Evaluate(x_vect)/Top_Gelin->Eval(x));
//    }
//    
//    Res_Top_1sig->SetMarkerSize(0);
//    Res_Top_1sig->SetMarkerColor(0);
//    Res_Top_1sig->SetLineColor(0);
//    Res_Top_1sig->SetStats(false);
//    Res_Top_1sig->SetFillColorAlpha(4,0.4);
//    if (plot) Res_Top_1sig->Draw("e3 same");
//    Res_Top_1sig->GetXaxis()->SetTitle("True energy [keV]");
//    Res_Top_1sig->GetYaxis()->SetTitle("Amplitude [V]");
//    Res_Top_1sig->SetTitle("1 #sigma");
//    
//    c1->cd();
//    

    //Calibration points LWO
    
    

////    //No shielding 2021_04_29_09_54_36
//    const int Npoints_LWO = 4;
//    Double_t Etrue3[Npoints_LWO] = {0, 1460.82, 2614.511, 15180};
//    Double_t EtrueErr3[Npoints_LWO] = {0, 0, 0, 20};
//    Double_t ChanNum3[Npoints_LWO]={0, 0.3741, 0.664, 3.989};
//    Double_t ChanNumErr3[Npoints_LWO]={0.007, 2*0.0007, 2*0.002, 0.1};

    //With shielding 2021_05_04_14_47_35
//    const int Npoints_LWO = 5;
//    Double_t Etrue3[Npoints_LWO] = {0, 511, 1460.82, 2614.511, 15407};
//    Double_t EtrueErr3[Npoints_LWO] = {0, 0, 0, 0, 100};
//    Double_t ChanNum3[Npoints_LWO]={0, 0.1207, 0.3502, 0.6228, 3.6};
//    Double_t ChanNumErr3[Npoints_LWO]={0.007, 2*0.0007, 2*0.0003, 2*0.0012, 0.05};
        
////    2021_05_12_17_27_27 Eddy
//    const int Npoints_LWO = 8;
//    Double_t Etrue3[Npoints_LWO] = {0., 511, 1460.82, 2614.511, 911.2, 968.97, 583.19, 15407};
//    Double_t EtrueErr3[Npoints_LWO] = {0., 0, 0, 0, 0, 10};
//    Double_t ChanNum3[Npoints_LWO]={0., 0.123094, 0.347698, 0.619407, 0.21697, 0.230071, 0.140852, 3.64945};
//    Double_t ChanNumErr3[Npoints_LWO]={0.007, 0.000221719, 0.000201991, 0.000469123, 0.000630377, 0.000547354, 0.000508196, 0.1};

//    //2021_05_12_17_27_27 Chloé
//    const int Npoints_LWO = 5;
//    Double_t Etrue3[Npoints_LWO] = {0., 511, 1460.82, 2614.511, 15407};
//    Double_t EtrueErr3[Npoints_LWO] = {0., 0, 0, 0, 10};
//    Double_t ChanNum3[Npoints_LWO]={0., 1.22679e-01, 3.53829e-01, 6.30289e-01, 3.64945};
//    Double_t ChanNumErr3[Npoints_LWO]={0.007, 0.001, 0.001, 0.002, 0.1};
//
//     TGraphErrors* gr3 = new TGraphErrors(Npoints_LWO,Etrue3,ChanNum3, EtrueErr3, ChanNumErr3);
//     //gr->SetTitle("Energy scale linearity from Co60");
//    gr3->SetTitle("Calibration Mid LWO");
//    gr3->SetMarkerStyle(20);
//    gr3->SetMarkerColor(kGreen+3);
//    gr3->SetMarkerSize(1.8);
//    gr3->GetXaxis()->SetTitle("True energy [keV]");
//    gr3->GetXaxis()->SetTitleOffset(1.2);
//    gr3->GetXaxis()->SetRangeUser(0,20000);
//    gr3->GetYaxis()->SetTitle("Amplitude [V]");
//     //gr->GetYaxis()->SetTitle("DESPEC bin [a.u.]");
//    gr3->GetYaxis()->SetTitleOffset(1.4);
//     //gr->GetYaxis()->SetRangeUser(0,1500);
//    gr3->SetMarkerStyle(20);
//    gr3->SetMarkerSize(1.8);
//
//    if (plot) gr3->Draw("P same");
//
//     //Fit
//     TF1* fitLWO = new TF1("LWO","([0]+[1]*x)",0,20000);
//    fitLWO->SetParNames("b [V]","a [V/keV]");
//
//     //Fit parameter initialization
//    fitLWO->SetParameters(par);
//     //fit2->SetParLimits(0,-10,10);
//    //fit2->SetParLimits(1,-1000,1000);
//
//    //gStyle->SetOptStat(0);
//    TFitResultPtr r3 = gr3->Fit("LWO","RS0");
//    FitResults.push_back(r3);
//
//    TMatrixTSym<double> LWO_rlin_Cov = r3->GetCovarianceMatrix();
//    cout<<"Chi2/Ndf = " << r3->Chi2()<<"/"<<r3->Ndf()<<" = "<<r3->Chi2()/r3->Ndf()<<endl;
//    LWO_rlin_Cov.Print();
//    
//    TMatrixTSym<double> LWO_rlin_Cor = r3->GetCorrelationMatrix();
//    TH2D* hLWO_rlin_Cor = new TH2D(LWO_rlin_Cor);
//    TCanvas *cCor3 = new TCanvas("LWO", "LWO");
//    gStyle->SetOptStat(0);
//    hLWO_rlin_Cor->GetZaxis()->SetRangeUser(-1, 1);
//    hLWO_rlin_Cor->SetTitle("LWO calibration correlation matrix");
//    if (plot) hLWO_rlin_Cor->Draw("colz");
//    c1->cd();
//    
//    /*Create a histogram to hold the confidence intervals*/
//    TH1D *LWO_1sig = new TH1D("LWO_1sig",
//       "LWO_1sig", Emax, 0, Emax);
//    classsigma1* sigma1_LWO = new classsigma1(Parameters_from_fit(r3), LWO_rlin_Cov);
//    for (int i=0; i<LWO_1sig->GetNbinsX(); i++){
//        double x = LWO_1sig->GetBinCenter(i+1);
//        LWO_1sig->SetBinContent(i+1, fitLWO->Eval(x));
//        double x_vect[1]={x};
//        LWO_1sig->SetBinError(i+1, 1*sigma1_LWO->Evaluate(x_vect));
//    }
//    
//    LWO_1sig->SetMarkerSize(0);
//    LWO_1sig->SetMarkerColor(0);
//    LWO_1sig->SetLineColor(0);
//    LWO_1sig->SetStats(false);
//    LWO_1sig->SetFillColorAlpha(8,0.4);
//    if (plot) LWO_1sig->Draw("e3 same");
//    LWO_1sig->SetTitle("1#times#sigma");
//    
//    fitLWO->SetLineWidth(4);
//    fitLWO->SetLineColor(8);
//    fitLWO->SetLineStyle(9);
//    if (plot) fitLWO->Draw("Same");
//    fitLWO->SetTitle(Form("#splitline{a = %0.1f +/- %0.1f #muV/keV    }{b = %0.1f +/- %0.1f #muV}", 1e6*fitLWO->GetParameter(1), 1e6*fitLWO->GetParError(1), 1e6*fitLWO->GetParameter(0), 1e6*fitLWO->GetParError(0)));
//    
//    if (plot) gr3->Draw("P same");
//    auto legend3 = new TLegend();
//    legend3->AddEntry(gr3);
//    legend3->AddEntry(fitLWO);
//    legend3->AddEntry(LWO_1sig);
//    if (plot) legend3->Draw("same");
//
//    // Residuals
//    cRes->cd();
//    
//    Double_t Res_value_LWO[Npoints_LWO]={0, 0, 0, 0};
//    Double_t Res_Errors_LWO[Npoints_LWO]={0, 0, 0, 0};
//    
//    for (int i = 0; i<Npoints_LWO; i++){
//        Value_i = ChanNum3[i];
//        Res_value_LWO[i] = (ChanNum3[i]-fitLWO->Eval(Etrue3[i]))/fitLWO->Eval(Etrue3[i])*100+shift_plot_LWO;
//        Res_Errors_LWO[i] = ChanNumErr3[i]/fitLWO->Eval(Etrue3[i])*100;
//    }
//    TGraphErrors* ResidualsLWO_pc = new TGraphErrors(Npoints_LWO,Etrue3,Res_value_LWO, EtrueErr3, Res_Errors_LWO);
//    
//    ResidualsLWO_pc->SetTitle("Residuals LWO");
//    ResidualsLWO_pc->SetMarkerStyle(20);
//    ResidualsLWO_pc->SetMarkerColor(kGreen+3);
//    ResidualsLWO_pc->SetLineColor(kGreen+3);
//    ResidualsLWO_pc->SetMarkerSize(1.8);
//    ResidualsLWO_pc->GetXaxis()->SetTitle("True energy [keV]");
//    ResidualsLWO_pc->GetXaxis()->SetTitleOffset(1.2);
//    ResidualsLWO_pc->GetXaxis()->SetRangeUser(0,20000);
//    ResidualsLWO_pc->GetYaxis()->SetTitle("Residuals [%]");
//    
//    TGaxis *axisLWO = new TGaxis(-2000,Min_axis, -2000, Max_axis,Min_axis-shift_plot_LWO,Max_axis-shift_plot_LWO,Ndiv,"");
//    axisLWO->SetLineColor(8);
//    axisLWO->SetLabelColor(8);
//    axisLWO->SetTextFont(72);
//    //axisLWO->SetTickLength(0.005);
//    axisLWO->SetLabelFont(hLWO_rlin_Cor->GetXaxis()->GetLabelFont());
//    axisLWO->SetLabelSize(hLWO_rlin_Cor->GetXaxis()->GetLabelSize());
//    if (plot) axisLWO->Draw("SAME");
//    
//    if (plot) ResidualsLWO_pc->Draw("PSAME");
//    
//    TH1D *Res_LWO_1sig = new TH1D("Res_LWO_1sig",
//       "Res_LWO_1sig", Emax, 0, Emax);
//    for (int i=0; i<Res_LWO_1sig->GetNbinsX(); i++){
//        double x = Res_LWO_1sig->GetBinCenter(i+1);
//        Res_LWO_1sig->SetBinContent(i+1, shift_plot_LWO);
//        double x_vect[1]={x};
//        Res_LWO_1sig->SetBinError(i+1, sigma1_LWO->Evaluate(x_vect)/fitLWO->Eval(x)*100);
//    }
//    
//    Res_LWO_1sig->SetMarkerSize(0);
//    Res_LWO_1sig->SetMarkerColor(0);
//    Res_LWO_1sig->SetLineColor(0);
//    Res_LWO_1sig->SetStats(false);
//    Res_LWO_1sig->SetFillColorAlpha(8,0.4);
//    if (plot) Res_LWO_1sig->Draw("e3 same");
//    Res_LWO_1sig->GetXaxis()->SetTitle("True energy [keV]");
//    Res_LWO_1sig->GetYaxis()->SetTitle("Amplitude [V]");
//    Res_LWO_1sig->SetTitle("1#sigma");
//    
//    auto legend_res = new TLegend(0.75,0.8,0.95,0.95);
//    legend_res->AddEntry(ResidualsTop_pc, "Top Ge");
//    legend_res->AddEntry(ResidualsLWO_pc, "LWO");
//    legend_res->AddEntry(ResidualsBot_pc, "Bot Ge");
//    legend_res->SetLineWidth(0);
//    if (plot) legend_res->Draw("same");
//    
//    if (!plot){
//        c1->Close();
//        cRes->Close();
//        cCor2->Close();
//        cCor->Close();
//        cCor3->Close();
//    }
//    
    return FitResults;
}

//------------------------- Covariance Matrix
double Cov(int index_i, int index_j, TMatrixD N, double mean_i, double mean_j, int Num_spe){
    double cov=0;
    for (int k=0; k<Num_spe; k++){
        cov += (N(k, index_i)-mean_i)*(N(k, index_j)-mean_j);
    }
    return cov/(Num_spe-1);
}

//------------------------- Correlation Matrix
double Corr(int index_i, int index_j, TMatrixD N){
    double var_i = N(index_i, index_i);
    double var_j = N(index_j, index_j);
    return N(index_i, index_j)/sqrt(var_i*var_j);
}

// --------------- Calculate mean spectrum from spectra
TVectorD Mean_spectrum(TMatrixD Spectra){
    
    int Nbins = Spectra.GetNcols();
    int N = Spectra.GetNrows();
    
    TVectorD Means(Nbins);
    
    double mean_j;
    for (int j=0; j<Nbins; j++){
        mean_j=0;
        for (int i=0; i<N; i++){
            mean_j +=Spectra(i,j);
        }
        Means(j)=mean_j/N;
    }
    return Means;
}

// --------------- Calculate covariance matrix from N spectra
TMatrixD Covariance_Matrix(TMatrixD Spectra){
    int Nbins = Spectra.GetNcols();
    int N = Spectra.GetNrows();
    
    TVectorD Means = Mean_spectrum(Spectra);
    TMatrixD Cov_Mat(Nbins, Nbins);
    
    for (int i=0; i<Nbins; i++){
        for (int j=0; j<Nbins; j++){
            Cov_Mat(i,j)=Cov(i, j, Spectra, Means(i), Means(j), N);
        }
    }
    return Cov_Mat;
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
//------------------------- Generate 3x two calib coeff from fit result
vector<Double_t> Generate_a_calib(TFitResultPtr FitResult, ROOT::Math::GSLRandomEngine* rnd, int ErrorFact=1){
    vector<Double_t> Calib_val;
    
    Double_t Center_b = FitResult->Parameter(0);
    Double_t Center_a = FitResult->Parameter(1);
    
    auto CovMat = FitResult->GetCovarianceMatrix();
    CovMat *= ErrorFact; // multiply the errors by ErrorFact (default 1)
    
    Double_t Sigma_b = sqrt(CovMat(0,0));
    Double_t Sigma_a = sqrt(CovMat(1,1));
    Double_t Corr_ab = CovMat(1,0)/(Sigma_b*Sigma_a);
    
    Double_t x=0, y=0;
    
    rnd->Gaussian2D(Sigma_a, Sigma_b, Corr_ab, x, y);
    
    Calib_val.push_back(y+Center_b);
    Calib_val.push_back(x+Center_a);

    return Calib_val;
}

//------------------- Function generate N correlated variables
TVectorD Generate_n_correlated_variables(TMatrixDSym Cov_Mat, TRandom3* rand_gen){
    int N_var = Cov_Mat.GetNrows();
    
    TVectorD random_var(N_var);
    
    // fill random_var with uncorrelated samples in a normal distribution centered in Means
    for (int i=0; i<N_var; i++){
        random_var(i) = rand_gen->Gaus(0, 1);
    }

    //Find the matrix C which verifies CC_T = Cov_Mat
    TDecompChol Cholesky(N_var);
    Cholesky.SetMatrix(Cov_Mat);
    Bool_t managed = Cholesky.Decompose();
    
    TMatrixD sqrt_Cov_Mat = Cholesky.GetU();
    TMatrixD sqrt_Cov_Mat_T(TMatrixD::kTransposed, sqrt_Cov_Mat);
    TMatrixD product = sqrt_Cov_Mat_T*sqrt_Cov_Mat;
    
//    if (managed) {
//        cout<< "Cholesky decomposition worked :"<<endl;
//        sqrt_Cov_Mat.Print();
//        cout<<"The product :"<<endl;
//        product.Print();
//        cout<<"should be equal to :"<<endl;
//        Cov_Mat.Print();
//    }
    
    //Multiply the variables by C
    TVectorD corr_var(N_var);
    corr_var = sqrt_Cov_Mat_T*random_var;
    
    return corr_var;
}

//------------------------------- Generate a calib with Cholesky method
vector<Double_t> Generate_a_calib_Chol(TFitResultPtr FitResult, int N_par, TRandom3* rnd, int ErrorFact=1){
    vector<Double_t> Calib_val;
    
    TVectorD Means(N_par);
    TVectorD Corr_var(N_par);
    TMatrixDSym Cov_Mat(N_par);
    
    Cov_Mat =FitResult->GetCovarianceMatrix();
    Cov_Mat*=ErrorFact;
    
    for (int i=0; i<N_par; i++){
        Means(i) = FitResult->Parameter(i);
    }
    
    Corr_var = Generate_n_correlated_variables(Cov_Mat, rnd);
    
    for (int i=0; i<N_par; i++){
        Calib_val.push_back(Corr_var(i)+Means(i));
    }

    return Calib_val;
}

//------------------------- Re calib the histograms
TH1D* RecalibHisto(TH1D* histo, int Nbins, double Emax, Double_t b, Double_t a, TString Name){
    double Volt_i, E_cal;
    int Ener_index, Hits_i, Err_i, NumBin;
    
    NumBin=histo->GetNbinsX();
    TString title = histo->GetTitle();
    TH1D *new_histo = new TH1D(Name, Name, Nbins, 0, Emax);
    
    for (int i=0; i<10000000; i++){
        Volt_i = histo->GetRandom();
        E_cal = ((Volt_i-b)/a)*1e-3;
        new_histo->Fill(E_cal);
    }
    
    //new_histo->Rebin(10);
    new_histo->GetXaxis()->SetTitle("Energy [MeV]");
    new_histo->GetYaxis()->SetTitle("Counts [a.u.]");

    double V_Max;
    V_Max= b+a*Emax*1E3;
    
    new_histo->Scale(histo->Integral(0, histo->FindBin(V_Max))/new_histo->Integral(0,new_histo->FindBin(Emax)));
    
    for (int i=1; i<new_histo->GetNbinsX(); i++){
        new_histo->SetBinError(i, sqrt(new_histo->GetBinContent(i)));
    }
    
    //new_histo->Rebin(10);
    new_histo->GetXaxis()->SetTitle("Energy [MeV]");
    new_histo->GetYaxis()->SetTitle("Counts [a.u.]");
    return new_histo;
}

//------------------------- Generate a Tree from the histograms
TTree* Tree_from_histos(TH1D* histoTop, TH1D* histoBot, TH1D* histoLWO, int NDraw=100000000){
    
    TTree* T = new TTree("Events_COV", "Events_COV");
    Double_t evTop;
    Double_t evBot;
    Double_t evLWO;
    T->Branch("Bot_Ge", &evTop);
    T->Branch("Top_Ge", &evBot);
    T->Branch("LWO", &evLWO);
    
    for (int i=0; i<NDraw; i++){
        evTop = histoTop->GetRandom();
        evBot = histoBot->GetRandom();
        evLWO = histoLWO->GetRandom();
        T->Fill();
    }

    return T;
}

//------------------------- Re calib the histograms from the Tree
TH1D* RecalibHisto_fromTree(TH1D* histo, TTree* T, int Nbins, double Emin, double Emax, vector<Double_t> Cal_pars, TString Name, TString BranchName){
    int N_par = Cal_pars.size();
    
    double Volt_i, E_cal;
    int Ener_index, Hits_i, Err_i, NumBin;
    
    TH1D *new_histo = new TH1D(Name, Name, Nbins, Emin, Emax);
    
    Double_t ev;
    T->SetBranchStatus("*", false);
    T->SetBranchStatus(BranchName, true);
    
    T->SetBranchAddress(BranchName, &ev);
    
    int NDraw = T->GetEntries();
    
    Double_t Delta, E_1, E_2, a, b, c;
    if (N_par==3) {a = Cal_pars[2];}
    b = Cal_pars[1];
    c = Cal_pars[0];
    for (int i=0; i<NDraw; i++){
        T->GetEntry(i);
        Volt_i = ev;
        if (N_par==2) {E_cal = ((Volt_i-c)/b)*1e-3;} //MeV
        if (N_par==3) {
            Delta = b*b-4*a*(c-Volt_i);
            E_1 = (-b-sqrt(Delta))/(2*a)*1e-3;
            E_2 = (-b+sqrt(Delta))/(2*a)*1e-3;
            if ((E_1>0)&&(E_2<0)){E_cal=E_1;}
            else if ((E_2>0)&&(E_1<0)){E_cal=E_2;}
            else if ((E_2>Emax)&&(E_1>Emin)){E_cal=E_1;}
            else if ((E_1>Emax)&&(E_2>Emin)){E_cal=E_2;}
            else if ((E_2<0)&&(E_1<0)){E_cal=-99;}
            else {cout<<"WARNING in calib : don't know how to choose between energy roots :"<<E_1<<" MeV and "<<E_2<<" MeV both possible : event has been ignored"<<endl;
                E_cal=-99;
            }
        }
        if (E_cal!= -99) new_histo->Fill(E_cal);
    }
    
    //new_histo->Rebin(10);
    new_histo->GetXaxis()->SetTitle("Energy [MeV]");
    new_histo->GetYaxis()->SetTitle("Counts [a.u.]");

    double V_Max;
    if (N_par==2) {V_Max = c+b*Emax*1e3;}
    if (N_par==3) {V_Max = c+b*Emax*1e3+a*Emax*Emax*1e6;}
    
    new_histo->Scale(histo->Integral(0, histo->FindBin(V_Max))/new_histo->Integral(0,new_histo->FindBin(Emax)));
    
    for (int i=1; i<new_histo->GetNbinsX(); i++){
        new_histo->SetBinError(i, sqrt(new_histo->GetBinContent(i)));
    }
    
    //new_histo->Rebin(10);
    new_histo->GetXaxis()->SetTitle("Energy [MeV]");
    new_histo->GetYaxis()->SetTitle("Counts [a.u.]");
//    new_histo->GetYaxis()->SetRangeUser(-10, 10);
    new_histo->GetYaxis()->SetRangeUser(1e-2, 1e6);
    return new_histo;
}

//---------------- Save Mean Calib in file
void Save_histo_statErr(TH1D* histo, TString Name_save){
    for (int i=1; i<histo->GetNbinsX()+1; i++){
        histo->SetBinError(i, sqrt(histo->GetBinContent(i)));
    }
    histo->Write(Name_save);
}

//---------------- Combine the two errors
TH1D* Combine_two_errors(TH1D* histo_err1, TH1D* histo_err2, TString Title){
    int Nbins = histo_err1->GetNbinsX();
    
    TH1D* histo_err1_2 = (TH1D*) histo_err1->Clone();
    histo_err1_2->SetTitle(Title);
    histo_err1_2->SetName(Title);
    
    double Err1, Err2;
    for (int i=1; i<Nbins+1; i++){
        Err1=histo_err1->GetBinError(i);
        Err2=histo_err2->GetBinError(i);
        histo_err1_2->SetBinError(i, sqrt(Err1*Err1+Err2*Err2));
    }
    
    return histo_err1_2;
}

//------------------------- Plot N calib and the mean spectrum and the correlation matrix
void Plot_N_calibs_Spectra(int N, TString File="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2.root", TDirectory* Dir=NULL, int Nbins=200, Double_t Emin= 0, Double_t Emax =5, int ErrorFact=1){
    
    TFile* filein = new TFile(File);
    
    TH1D* histo_0_Bot_Ge = (TH1D*) filein->Get("h1_Peak1");
    TH1D* histo_0_Top_Ge = (TH1D*) filein->Get("h1_Peak2");
    TH1D* histo_0_LWO = (TH1D*) filein->Get("h1_NTD_Amplitude");
    
    vector<TFitResultPtr> Fit_results = PlotEnergyCalibration();
    
    TCanvas* cBot_Ge = new TCanvas();
    cBot_Ge->Draw();
    
    TCanvas* cTop_Ge = new TCanvas();
    cTop_Ge->Draw();
    
    TCanvas* cLWO = new TCanvas();
    cLWO->Draw();
    
    ROOT::Math::GSLRandomEngine* rnd = new ROOT::Math::GSLRandomEngine();
    rnd->Initialize();
    rnd->SetSeed(33566);
    rnd->Initialize();
    
    TRandom3* rnd3 = new TRandom3();
    
    TMatrixD Spectra_Bot_Ge(N, Nbins);
    TMatrixD Spectra_Top_Ge(N, Nbins);
    TMatrixD Spectra_LWO(N, Nbins);
    
    TTree* T;
    TString Namefiletree = File;
    TPMERegexp("\\.root").Substitute(Namefiletree,"_filetree.root");
    
    TFile *filetree = TFile::Open(Namefiletree);
    
    if ((!filetree) || (filetree->IsZombie())){
        delete filetree;
        cout<<"Creating "<<Namefiletree<<"..."<<endl;
        TFile* filetree = new TFile(Namefiletree, "CREATE");
        T = Tree_from_histos(histo_0_Bot_Ge, histo_0_Top_Ge, histo_0_LWO, 10000000);
        T->Write("Events_COV");
        filetree->Close();
        filetree = TFile::Open(Namefiletree);
        T=(TTree*) filetree->Get("Events_COV");
    }
    else {
        T=(TTree*) filetree->Get("Events_COV");
    }
    
    if (Dir!=NULL){
        Dir->cd();
        Dir->mkdir("calibrated_histos");
        Dir->mkdir("calibrated_histos/TopGe");
        Dir->mkdir("calibrated_histos/BotGe");
        Dir->mkdir("calibrated_histos/LWO");
    }
    
    //--------------- Calibration Bot Ge
    cBot_Ge->cd();
    //rnd->Initialize();
    //vector<Double_t> Calib_Bot_Ge = Generate_a_calib(Fit_results[0], rnd); //with 2 variables: use the Gaussian2D function from root
    vector<Double_t> Calib_Bot_Ge = Generate_a_calib_Chol(Fit_results[0], 3, rnd3, ErrorFact);
    
    Calib_Bot_Ge[0]=Fit_results[0]->Parameter(0);
    Calib_Bot_Ge[1]=Fit_results[0]->Parameter(1);
    Calib_Bot_Ge[2]=Fit_results[0]->Parameter(2);
    
    TH1D* calib_hist_Bot_Ge_ref = RecalibHisto_fromTree(histo_0_Bot_Ge, T, Nbins, Emin, Emax, Calib_Bot_Ge, "Bot Ge", "Bot_Ge");
    //TH1D* calib_hist_Bot_Ge = RecalibHisto(histo_0_Bot_Ge, Nbins, 20, Calib_Bot_Ge[0], Calib_Bot_Ge[1], "histo_Bot_Ge_0");
    calib_hist_Bot_Ge_ref->Draw("HIST");
    
    Save_histo_statErr(calib_hist_Bot_Ge_ref, "BotGe_Ref_calib_statErrors");
    // -----------------------------------
    
    //--------------- Calibration Top Ge
    cTop_Ge->cd();
    //rnd->Initialize();
//    vector<Double_t> Calib_Top_Ge = Generate_a_calib(Fit_results[1], rnd);
    vector<Double_t> Calib_Top_Ge = Generate_a_calib_Chol(Fit_results[1], 2, rnd3, ErrorFact);
    
    Calib_Top_Ge[0]=Fit_results[1]->Parameter(0);
    Calib_Top_Ge[1]=Fit_results[1]->Parameter(1);
    
    TH1D* calib_hist_Top_Ge_ref = RecalibHisto_fromTree(histo_0_Top_Ge, T, Nbins, Emin, Emax, Calib_Top_Ge, "Top Ge", "Top_Ge");
//    TH1D* calib_hist_Top_Ge = RecalibHisto(histo_0_Top_Ge, Nbins, 20, Calib_Top_Ge[0], Calib_Top_Ge[1], "histo_Top_Ge_0");
    calib_hist_Top_Ge_ref->Draw("HIST");
    
    Save_histo_statErr(calib_hist_Top_Ge_ref, "TopGe_Ref_calib_statErrors");
    // --------------------------------------
    
    //--------------- Calibration LWO
    cLWO->cd();
    //rnd->Initialize();
//    vector<Double_t> Calib_LWO = Generate_a_calib(Fit_results[2], rnd);
    vector<Double_t> Calib_LWO = Generate_a_calib_Chol(Fit_results[2], 2, rnd3, ErrorFact);
    
    Calib_LWO[0]=Fit_results[2]->Parameter(0);
    Calib_LWO[1]=Fit_results[2]->Parameter(1);
    
    TH1D* calib_hist_LWO_ref = RecalibHisto_fromTree(histo_0_LWO, T, Nbins, Emin, Emax, Calib_LWO, "LWO", "LWO");
    //TH1D* calib_hist_LWO = RecalibHisto(histo_0_LWO, Nbins, 20, Calib_LWO[0], Calib_LWO[1], "histo_LWO_0");
    calib_hist_LWO_ref->Draw("HIST");
    
    Save_histo_statErr(calib_hist_LWO_ref, "LWO_Ref_calib_statErrors");
    // --------------------------------------
    
    for (int k=0; k<Nbins; k++){
        Spectra_Bot_Ge(0, k)=calib_hist_Bot_Ge_ref->GetBinContent(k+1);
        Spectra_Top_Ge(0, k)=calib_hist_Top_Ge_ref->GetBinContent(k+1);
        Spectra_LWO(0, k)=calib_hist_LWO_ref->GetBinContent(k+1);
    }
    gStyle->SetOptStat(0);
    
    TH1D* calib_hist_Bot_Ge = (TH1D*) calib_hist_Bot_Ge_ref->Clone();
    calib_hist_Bot_Ge->SetLineColor(kGray);

    TH1D* calib_hist_Top_Ge = (TH1D*) calib_hist_Top_Ge_ref->Clone();
    calib_hist_Top_Ge->SetLineColor(kGray);
    
    TH1D* calib_hist_LWO = (TH1D*) calib_hist_LWO_ref->Clone();
    calib_hist_LWO->SetLineColor(kGray);
    
    for (int i=1; i<N; i++){
        cout<<"i = "<<i<<endl;
        
        //--------------- Calibration Bot Ge
        cBot_Ge->cd();
        Calib_Bot_Ge = Generate_a_calib_Chol(Fit_results[0], 3, rnd3, ErrorFact);
        
        calib_hist_Bot_Ge = RecalibHisto_fromTree(histo_0_Bot_Ge, T, Nbins, Emin, Emax, Calib_Bot_Ge, Form("histo_Bot_Ge_%d", i), "Bot_Ge");
        calib_hist_Bot_Ge->SetLineColor(kGray);
        calib_hist_Bot_Ge->Draw("SAME HIST");
        
        if (Dir!=NULL){
            Dir->cd("calibrated_histos/BotGe");
            calib_hist_Bot_Ge->Write(Form("histo_Bot_Ge_%d", i));
            Dir->cd();
        }
        //----------------------------------------
        
        //---------------- Calibration Top Ge
        cTop_Ge->cd();
        Calib_Top_Ge = Generate_a_calib_Chol(Fit_results[1], 2, rnd3, ErrorFact);
        
        calib_hist_Top_Ge = RecalibHisto_fromTree(histo_0_Top_Ge, T, Nbins, Emin, Emax, Calib_Top_Ge, Form("histo_Top_Ge_%d", i), "Top_Ge");
        calib_hist_Top_Ge->SetLineColor(kGray);
        calib_hist_Top_Ge->Draw("SAME HIST");
        
        if (Dir!=NULL){
            Dir->cd("calibrated_histos/TopGe");
            calib_hist_Top_Ge->Write(Form("histo_Top_Ge_%d", i));
            Dir->cd();
        }
        //----------------------------------------

        //---------------- Calibration LWO
        cLWO->cd();
        Calib_LWO = Generate_a_calib_Chol(Fit_results[2], 2, rnd3, ErrorFact);
        
        calib_hist_LWO = RecalibHisto_fromTree(histo_0_LWO, T, Nbins, Emin, Emax, Calib_LWO, Form("histo_LWO_%d", i), "LWO");
        calib_hist_LWO->SetLineColor(kGray);
        calib_hist_LWO->Draw("SAME HIST");
        
        if ((Dir!=NULL)&&(i%100==0)){
            Dir->cd("calibrated_histos/LWO");
            calib_hist_LWO->Write(Form("histo_LWO_%d", i));
            Dir->cd();
        }
        //----------------------------------------
        
        for (int k=0; k<Nbins; k++){
            Spectra_Bot_Ge(i, k)=calib_hist_Bot_Ge->GetBinContent(k+1);
            Spectra_Top_Ge(i, k)=calib_hist_Top_Ge->GetBinContent(k+1);
            Spectra_LWO(i, k)=calib_hist_LWO->GetBinContent(k+1);
        }
    }
    
    // Calculate the mean spectrum from all the calculated spectra
    TVectorD Means_Bot_Ge = Mean_spectrum(Spectra_Bot_Ge);
    TVectorD Means_Top_Ge = Mean_spectrum(Spectra_Top_Ge);
    TVectorD Means_LWO = Mean_spectrum(Spectra_LWO);
    
    // Calculate the covariance matrices
    TMatrixD Cov_Mat_Bot = Covariance_Matrix(Spectra_Bot_Ge);
    TMatrixD Cov_Mat_Top = Covariance_Matrix(Spectra_Top_Ge);
    TMatrixD Cov_Mat_LWO = Covariance_Matrix(Spectra_LWO);
    
    // Calculate the stat covariance matrices
    TMatrixD StatCov_Mat_Bot = Stat_Covariance_Matrix(calib_hist_Bot_Ge_ref);
    TMatrixD StatCov_Mat_Top = Stat_Covariance_Matrix(calib_hist_Top_Ge_ref);
    TMatrixD StatCov_Mat_LWO = Stat_Covariance_Matrix(calib_hist_LWO_ref);
    
    
    TH2D* hCov_Mat_Top = get_TH2D_from_TMatrixD("TopGe_SystCovMat", Cov_Mat_Top, "Covariance", "Energy [MeV]", "Energy [MeV]", Emin, Emax, Emin, Emax);
    TH2D* hCov_Mat_Bot = get_TH2D_from_TMatrixD("BotGe_SystCovMat", Cov_Mat_Bot, "Covariance", "Energy [MeV]", "Energy [MeV]", Emin, Emax, Emin, Emax);
    TH2D* hCov_Mat_LWO = get_TH2D_from_TMatrixD("LWO_SystCovMat", Cov_Mat_LWO, "Covariance", "Energy [MeV]", "Energy [MeV]", Emin, Emax, Emin, Emax);
    
    TH2D* hStatCov_Mat_Top = get_TH2D_from_TMatrixD("TopGe_StatCovMat", StatCov_Mat_Top, "Covariance", "Energy [MeV]", "Energy [MeV]", Emin, Emax, Emin, Emax);
    TH2D* hStatCov_Mat_Bot = get_TH2D_from_TMatrixD("BotGe_StatCovMat", StatCov_Mat_Bot, "Covariance", "Energy [MeV]", "Energy [MeV]", Emin, Emax, Emin, Emax);
    TH2D* hStatCov_Mat_LWO = get_TH2D_from_TMatrixD("LWO_StatCovMat", StatCov_Mat_LWO, "Covariance", "Energy [MeV]", "Energy [MeV]", Emin, Emax, Emin, Emax);
    
    if (Dir!=NULL){
        hCov_Mat_Top->Write();
        hCov_Mat_Bot->Write();
        hCov_Mat_LWO->Write();
        
        hStatCov_Mat_Top->Write();
        hStatCov_Mat_Bot->Write();
        hStatCov_Mat_LWO->Write();
    }
    
    // Mean TH1D with syst errors - copy the format
    TH1D* hMeans_Bot_Ge = (TH1D*) calib_hist_Bot_Ge_ref->Clone();
    TH1D* hMeans_Top_Ge = (TH1D*) calib_hist_Top_Ge_ref->Clone();
    TH1D* hMeans_LWO = (TH1D*) calib_hist_LWO_ref->Clone();

    hMeans_Bot_Ge->SetLineColor(kBlue+2);
    hMeans_Top_Ge->SetLineColor(kBlue+2);
    hMeans_LWO->SetLineColor(kBlue+2);
    
    hMeans_Bot_Ge->SetTitle("Bot_Ge_mean");
    hMeans_Top_Ge->SetTitle("Top_Ge_mean");
    hMeans_LWO->SetTitle("LWO_mean");
    
    hMeans_Bot_Ge->SetName("Bot_Ge_mean");
    hMeans_Top_Ge->SetName("Top_Ge_mean");
    hMeans_LWO->SetName("LWO_mean");
    
    for (int j=0; j<Nbins; j++){
        hMeans_Bot_Ge->SetBinContent(j+1, Means_Bot_Ge(j));
        hMeans_Bot_Ge->SetBinError(j+1, sqrt(Cov_Mat_Bot(j,j)));
        
        hMeans_Top_Ge->SetBinContent(j+1, Means_Top_Ge(j));
        hMeans_Top_Ge->SetBinError(j+1, sqrt(Cov_Mat_Top(j,j)));
        
        hMeans_LWO->SetBinContent(j+1, Means_LWO(j));
        hMeans_LWO->SetBinError(j+1, sqrt(Cov_Mat_LWO(j,j)));
    }
    
    // Ref spectrum TH1D with stat errors - copy the format
    TH1D* hRef_Bot_Ge = (TH1D*) calib_hist_Bot_Ge_ref->Clone();
    TH1D* hRef_Top_Ge = (TH1D*) calib_hist_Top_Ge_ref->Clone();
    TH1D* hRef_LWO = (TH1D*) calib_hist_LWO_ref->Clone();
    
    hRef_Bot_Ge->SetLineColor(kBlue+2);
    hRef_Top_Ge->SetLineColor(kBlue+2);
    hRef_LWO->SetLineColor(kBlue+2);
    
    hRef_Bot_Ge->SetTitle("Bot_Ge_ref_cal_err");
    hRef_Top_Ge->SetTitle("Top_Ge_ref_cal_err");
    hRef_LWO->SetTitle("LWO_ref_cal_err");
    
    hRef_Bot_Ge->SetName("Bot_Ge_ref_cal_err");
    hRef_Top_Ge->SetName("Top_Ge_ref_cal_err");
    hRef_LWO->SetName("LWO_ref_cal_err");
    
    // set syst errors from cal
    for (int j=0; j<Nbins; j++){
        hRef_Bot_Ge->SetBinError(j+1, sqrt(Cov_Mat_Bot(j,j)));
        hRef_Top_Ge->SetBinError(j+1, sqrt(Cov_Mat_Top(j,j)));
        hRef_LWO->SetBinError(j+1, sqrt(Cov_Mat_LWO(j,j)));
    }
    
    cBot_Ge->cd();
    hMeans_Bot_Ge->Draw("SAME");
    cTop_Ge->cd();
    hMeans_Top_Ge->Draw("SAME");
    cLWO->cd();
    hMeans_LWO->Draw("SAME");
    
    if (Dir!=NULL){
        hMeans_Bot_Ge->Write("BotGe_Mean_calib_systErrors");
        hMeans_Top_Ge->Write("TopGe_Mean_calib_systErrors");
        hMeans_LWO->Write("LWO_Mean_calib_systErrors");
        
        hRef_Bot_Ge->Write("BotGe_Ref_calib_systErrors");
        hRef_Top_Ge->Write("TopGe_Ref_calib_systErrors");
        hRef_LWO->Write("LWO_Ref_calib_systErrors");
    }
    
    // Ref spectrum both errors - copy the format
    TH1D* hRef_bothErr_Bot_Ge = Combine_two_errors(calib_hist_Bot_Ge_ref, hRef_Bot_Ge, "BotGe_Ref_calib_bothErrors");
    TH1D* hRef_bothErr_Top_Ge = Combine_two_errors(calib_hist_Top_Ge_ref, hRef_Top_Ge, "TopGe_Ref_calib_bothErrors");
    TH1D* hRef_bothErr_LWO = Combine_two_errors(calib_hist_LWO_ref, hRef_LWO, "LWO_Ref_calib_bothErrors");
    
    hRef_bothErr_Bot_Ge->SetLineColor(kBlue+2);
    hRef_bothErr_Top_Ge->SetLineColor(kBlue+2);
    hRef_bothErr_LWO->SetLineColor(kBlue+2);
    
    if (Dir!=NULL){
        hRef_bothErr_Bot_Ge->Write("BotGe_Ref_calib_bothErrors");
        hRef_bothErr_Top_Ge->Write("TopGe_Ref_calib_bothErrors");
        hRef_bothErr_LWO->Write("LWO_Ref_calib_bothErrors");
    }
    
    //------------------ Bot Ge correlation matrix
    TH2D* Corr_Mat_Bot = new TH2D("CorrelationMatrix_Bot", "CorrelationMatrix_Bot",
                                  Nbins, hMeans_Bot_Ge->GetBinLowEdge(1), hMeans_Bot_Ge->GetBinLowEdge(Nbins)+hMeans_Bot_Ge->GetBinWidth(1),
                                  Nbins, hMeans_Bot_Ge->GetBinLowEdge(1), hMeans_Bot_Ge->GetBinLowEdge(Nbins)+hMeans_Bot_Ge->GetBinWidth(1));
    Corr_Mat_Bot->GetXaxis()->SetTitle(hMeans_Bot_Ge->GetXaxis()->GetTitle());
    Corr_Mat_Bot->GetYaxis()->SetTitle(hMeans_Bot_Ge->GetXaxis()->GetTitle());

    //------------------ Top Ge correlation matrix
    TH2D* Corr_Mat_Top = new TH2D("CorrelationMatrix_Top", "CorrelationMatrix_Top",
                                  Nbins, hMeans_Top_Ge->GetBinLowEdge(1), hMeans_Top_Ge->GetBinLowEdge(Nbins)+hMeans_Top_Ge->GetBinWidth(1),
                                  Nbins, hMeans_Top_Ge->GetBinLowEdge(1), hMeans_Top_Ge->GetBinLowEdge(Nbins)+hMeans_Top_Ge->GetBinWidth(1));
    Corr_Mat_Top->GetXaxis()->SetTitle(hMeans_Top_Ge->GetXaxis()->GetTitle());
    Corr_Mat_Top->GetYaxis()->SetTitle(hMeans_Top_Ge->GetXaxis()->GetTitle());

    //------------------ LWO correlation matrix
    TH2D* Corr_Mat_LWO = new TH2D("CorrelationMatrix_LWO", "CorrelationMatrix_LWO",
                                  Nbins, hMeans_LWO->GetBinLowEdge(1), hMeans_LWO->GetBinLowEdge(Nbins)+hMeans_LWO->GetBinWidth(1),
                                  Nbins, hMeans_LWO->GetBinLowEdge(1), hMeans_LWO->GetBinLowEdge(Nbins)+hMeans_LWO->GetBinWidth(1));
    Corr_Mat_LWO->GetXaxis()->SetTitle(hMeans_LWO->GetXaxis()->GetTitle());
    Corr_Mat_LWO->GetYaxis()->SetTitle(hMeans_LWO->GetXaxis()->GetTitle());

    
    for (int i=0; i<Nbins; i++){
        for (int j=0; j<Nbins; j++){
            if (TMath::IsNaN(Corr(i, j, Cov_Mat_Top))){
//                cout<<"Nan correlation between"<<i << " and "<<j<<" =" <<Cov_Mat_Top(i, j)<<"/"<<"sqrt("<<Cov_Mat_Top(i,i)<<"*"<<Cov_Mat_Top(j,j)<<endl;
                Corr_Mat_Top->SetBinContent(i+1, j+1, -99);
            }
            else Corr_Mat_Top->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_Top));
           
            if (TMath::IsNaN(Corr(i, j, Cov_Mat_Bot))){
//                cout<<"Nan correlation between"<<i << " and "<<j<<" =" <<Cov_Mat_Bot(i, j)<<"/"<<"sqrt("<<Cov_Mat_Bot(i,i)<<"*"<<Cov_Mat_Bot(j,j)<<endl;
                Corr_Mat_Bot->SetBinContent(i+1, j+1, -99);
            }
            else Corr_Mat_Bot->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_Bot));
            
            if (TMath::IsNaN(Corr(i, j, Cov_Mat_LWO))){
//                cout<<"Nan correlation between"<<i << " and "<<j<<" = "<<Cov_Mat_LWO(i, j)<<"/"<<"sqrt("<<Cov_Mat_LWO(i,i)<<"*"<<Cov_Mat_LWO(j,j)<<endl;
                Corr_Mat_LWO->SetBinContent(i+1, j+1, -99);
            }
            else Corr_Mat_LWO->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_LWO));
        }
    }

    // Bot Ge correlation matrix
    TCanvas* cMat_Bot_Ge = new TCanvas();
    Corr_Mat_Top->Draw("colz");
    Corr_Mat_Top->GetZaxis()->SetRangeUser(-1, 1);
    
    if (Dir!=NULL) Corr_Mat_Top->Write("TopGe_CorrMat");

    // Top Ge correlation matrix
    TCanvas* cMat_Top_Ge = new TCanvas();
    Corr_Mat_Bot->Draw("colz");
    Corr_Mat_Bot->GetZaxis()->SetRangeUser(-1, 1);

    if (Dir!=NULL) Corr_Mat_Bot->Write("BotGe_CorrMat");

    // LWO correlation matrix
    TCanvas* cMat_LWO = new TCanvas();
    Corr_Mat_LWO->Draw("colz");
    Corr_Mat_LWO->GetZaxis()->SetRangeUser(-1, 1);
    
    if (Dir!=NULL) Corr_Mat_LWO->Write("LWO_CorrMat");
}

//---------------- Use the previous function to save in a file with different energy ranges ----------
void Save_Syst_errors(int N, TString FileinName="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_nonlin.root", TString FileoutName="../Data_COV/2021_04_29_09_54_36/MC_Syst_Stat_2021_04_29_09_54_36_nonlin_BinforFit.root"){
    
    TFile* Fileout= new TFile(FileoutName, "RECREATE");
   
    cout<<"Range = 0-1MeV"<<endl;
    TDirectory* ThisRangeDir = Fileout->mkdir("0-1MeV");
    Plot_N_calibs_Spectra(N, FileinName, ThisRangeDir, 500, 0, 1, 1);
    
    cout<<"Range = 0-4MeV"<<endl;
    ThisRangeDir = Fileout->mkdir("0-4MeV");
    Plot_N_calibs_Spectra(N, FileinName, ThisRangeDir, 200, 0, 4, 1);

    cout<<"\nRange = 0-10MeV"<<endl;
    ThisRangeDir = Fileout->mkdir("0-10MeV");
    Plot_N_calibs_Spectra(N, FileinName, ThisRangeDir, 500, 0, 10, 1);
    
    cout<<"\nRange = 0-20MeV"<<endl;
    ThisRangeDir = Fileout->mkdir("0-20MeV");
    Plot_N_calibs_Spectra(N, FileinName, ThisRangeDir, 100, 0, 20, 1);
    
    Fileout->Close();
}

// Do the normalized difference between ref and mean and project in an histogram
TH1D* diff_h1_h2_proj(TH1D* histo_ref, TH1D* histo_mean, TString Name){
    
    int Nbins = histo_mean->GetNbinsX();
    TH1D* hdist_to_mean = new TH1D(Name, Name, 100, -100, 100);
    
    double Mean_content, Ref_content;
    for (int i=0; i<Nbins; i++){
        Mean_content= histo_mean->GetBinContent(i+1);
        Ref_content= histo_ref->GetBinContent(i+1);
        hdist_to_mean->Fill(100*(Ref_content-Mean_content)/Mean_content);
    }
    hdist_to_mean->GetXaxis()->SetTitle("(N^{ref}_{i} - #mu_{i})/#mu_{i} [%]");
    hdist_to_mean->GetYaxis()->SetTitle("Hits []");
    
    return hdist_to_mean;
}

// Do the normalized difference between ref and mean
TH1D* diff_h1_h2(TH1D* histo_ref, TH1D* histo_mean, TString Name, Double_t Emin =0, Double_t Emax =20){
    
    int Nbins = histo_mean->FindBin(Emax)-histo_mean->FindBin(Emin);
    TH1D* hdist_to_mean = new TH1D(Name, Name, Nbins, Emin, Emax);
    
    double Mean_content, Ref_content;
    for (int i=histo_mean->FindBin(Emin); i<histo_mean->FindBin(Emin)+Nbins; i++){
        Mean_content= histo_mean->GetBinContent(i+1);
        Ref_content= histo_ref->GetBinContent(i+1);
        hdist_to_mean->SetBinContent(i+1, 100*(Ref_content-Mean_content)/Mean_content);
    }
    hdist_to_mean->GetYaxis()->SetTitle("(N^{ref}_{i} - #mu_{i})/#mu_{i} [%]");
    hdist_to_mean->GetXaxis()->SetTitle("Energy [MeV]");
    
    return hdist_to_mean;
}

//Read file from Save_Syst_errors and calculate ref-mu
void Diff_ref_Mu_syst_file(TString FileinName="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2_SystErrors.root"){
    TFile* Filein = new TFile(FileinName, "READ");
    
    //Top
    TH1D* histo_mean = (TH1D*) Filein->Get("0-20MeV/TopGe_Mean_calib_systErrors");
    TH1D* histo_ref = (TH1D*) Filein->Get("0-20MeV/TopGe_Ref_calib_statErrors");
    
    TH1D* hdist_to_mean_proj_Top = diff_h1_h2_proj(histo_ref,histo_mean, "Top Ge proj");
    TH1D* hdist_to_mean_Top = diff_h1_h2(histo_ref,histo_mean, "Top Ge");
    
    //Bot
    histo_mean = (TH1D*) Filein->Get("0-20MeV/BotGe_Mean_calib_systErrors");
    histo_ref = (TH1D*) Filein->Get("0-20MeV/BotGe_Ref_calib_statErrors");
    
    TH1D* hdist_to_mean_proj_Bot = diff_h1_h2_proj(histo_ref,histo_mean, "Bot Ge proj");
    TH1D* hdist_to_mean_Bot = diff_h1_h2(histo_ref,histo_mean, "Bot Ge");
    hdist_to_mean_Bot->SetLineColor(kGreen);
    hdist_to_mean_proj_Bot->SetLineColor(kGreen);
    
    //LWO
    histo_mean = (TH1D*) Filein->Get("0-20MeV/LWO_Mean_calib_systErros");
    histo_ref = (TH1D*) Filein->Get("0-20MeV/LWO_Ref_calib_statErrors");
    
    TH1D* hdist_to_mean_LWO = diff_h1_h2(histo_ref,histo_mean, "LWO");
    TH1D* hdist_to_mean_proj_LWO = diff_h1_h2_proj(histo_ref,histo_mean, "LWO proj");
    hdist_to_mean_proj_LWO->SetLineColor(kRed);
    hdist_to_mean_LWO->SetLineColor(kRed);
    
    TCanvas* c1 = new TCanvas("Spectrum", "Spectrum");
    hdist_to_mean_Top->Draw();
    hdist_to_mean_Bot->Draw("SAME");
    hdist_to_mean_LWO->Draw("SAME");
    
    TCanvas* c2 = new TCanvas("Projection", "Projection");
    hdist_to_mean_proj_Top->Draw();
    hdist_to_mean_proj_Bot->Draw("SAME");
    hdist_to_mean_proj_LWO->Draw("SAME");
}

// ------------------------- Quantify the impact of the MC stat
TH1D* MC_stat(TString FileinName_ref="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2_SystErrors.root", TString FileinName_MC="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2_SystErrors.root", TString Name="MC=500", TString Name_histo_ref="0-20MeV/TopGe_Mean_calib_systErros", TString Name_histo_MC="0-20MeV/TopGe_Mean_calib_systErrors"){
    
    TFile* Filein_ref = new TFile(FileinName_ref, "READ");
    TFile* Filein_MC = new TFile(FileinName_MC, "READ");
    
    //Top
    TH1D* histo_ref = (TH1D*) Filein_ref->Get(Name_histo_ref);
    TH1D* histo_MC = (TH1D*) Filein_MC->Get(Name_histo_MC);
    
    int Nbins = histo_ref->GetNbinsX();
    int Emax = histo_ref->GetBinLowEdge(Nbins)+histo_ref->GetBinWidth(Nbins);
    int Emin = histo_ref->GetBinLowEdge(1);
    TH1D* h_MC_stat = new TH1D(Name, Name, Nbins, Emin, Emax);
    
    double MC_error, Ref_error;
    for (int i=0; i<Nbins; i++){
        MC_error= histo_MC->GetBinError(i+1);
        Ref_error= histo_ref->GetBinError(i+1);
        h_MC_stat->SetBinContent(i+1, 100*(Ref_error-MC_error)/Ref_error);
    }
    h_MC_stat->GetXaxis()->SetTitle("Energy [MeV]");
    
    return h_MC_stat;
}

// projected
TH1D* MC_stat_proj(TString FileinName_ref="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2_SystErrors.root", TString FileinName_MC="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2_SystErrors.root", TString Name="MC=500", TString Name_histo_ref="0-20MeV/TopGe_Mean_calib_systErros", TString Name_histo_MC="0-20MeV/TopGe_Mean_calib_systErrors"){
    
    TFile* Filein_ref = new TFile(FileinName_ref, "READ");
    TFile* Filein_MC = new TFile(FileinName_MC, "READ");
    
    //Top
    TH1D* histo_ref = (TH1D*) Filein_ref->Get(Name_histo_ref);
    TH1D* histo_MC = (TH1D*) Filein_MC->Get(Name_histo_MC);
    
    int Nbins = histo_ref->GetNbinsX();
    TH1D* h_MC_stat = new TH1D(Name, Name, 100, -100, 100);
    
    double MC_error, Ref_error;
    for (int i=0; i<Nbins; i++){
        MC_error= histo_MC->GetBinError(i+1);
        Ref_error= histo_ref->GetBinError(i+1);
        h_MC_stat->Fill(100*(Ref_error-MC_error)/Ref_error);
    }
    h_MC_stat->GetYaxis()->SetTitle("Hits []");
    
    return h_MC_stat;
}

void MC_stat_tot(){
    
    cout<<"Top"<<endl;
    int i=0;
    //--- Top
    cout<<i++<<endl;
    TH1D* h_MC_stat_Top_10 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC10.root", "MC=10", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_10->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    cout<<i++<<endl;
    
    TH1D* h_MC_stat_Top_50 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC50.root", "MC=50", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_50->SetLineColor(kGreen);
    h_MC_stat_Top_50->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    cout<<i++<<endl;
    
    TH1D* h_MC_stat_Top_100 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC100.root", "MC=100", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_100->SetLineColor(kRed);
    h_MC_stat_Top_100->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    cout<<i++<<endl;
    
    TH1D* h_MC_stat_Top_500 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC500.root", "MC=500", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_500->SetLineColor(kBlue);
    h_MC_stat_Top_500->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    cout<<i++<<endl;
    
    TH1D* h_MC_stat_Top_1000 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC1000.root", "MC=1000", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_1000->SetLineColor(kYellow+2);
    h_MC_stat_Top_1000->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    cout<<i++<<endl;
    
    TCanvas* C1= new TCanvas("Top Projected", "Top Projected");
    h_MC_stat_Top_1000->Draw("");
    h_MC_stat_Top_500->Draw("SAME");
    h_MC_stat_Top_100->Draw("SAME");
    h_MC_stat_Top_50->Draw("SAME");
    h_MC_stat_Top_10->Draw("SAME");
    
    //Spectrum
    TH1D* h_MC_stat_Top_10_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC10.root", "MC=10", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_10_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Top_50_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC50.root", "MC=50", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_50_spec->SetLineColor(kGreen);
    h_MC_stat_Top_50_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Top_100_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC100.root", "MC=100", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_100_spec->SetLineColor(kRed);
    h_MC_stat_Top_100_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Top_500_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC500.root", "MC=500", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_500_spec->SetLineColor(kBlue);
    h_MC_stat_Top_500_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Top_1000_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC1000.root", "MC=1000", "0-20MeV/TopGe_Mean_calib_systErrors", "0-20MeV/TopGe_Mean_calib_systErrors");
    h_MC_stat_Top_1000_spec->SetLineColor(kYellow+2);
    h_MC_stat_Top_1000_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TCanvas* C2= new TCanvas("Top Spectrum", "Top Spectrum");
    h_MC_stat_Top_10_spec->Draw("");
    h_MC_stat_Top_50_spec->Draw("SAME");
    h_MC_stat_Top_100_spec->Draw("SAME");
    h_MC_stat_Top_500_spec->Draw("SAME");
    h_MC_stat_Top_1000_spec->Draw("SAME");

    //--- Bot
    cout<<"\nBot"<<endl;
    TH1D* h_MC_stat_Bot_10 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC10.root", "MC=10", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_10->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_50 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC50.root", "MC=50", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_50->SetLineColor(kGreen);
    h_MC_stat_Bot_50->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_100 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC100.root", "MC=100", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_100->SetLineColor(kRed);
    h_MC_stat_Bot_100->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_500 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC500.root", "MC=500", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_500->SetLineColor(kBlue);
    h_MC_stat_Bot_500->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_1000 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC1000.root", "MC=1000", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_1000->SetLineColor(kYellow+2);
    h_MC_stat_Bot_1000->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");

    TCanvas* cBot = new TCanvas("Bot projected", "Bot projected");
    h_MC_stat_Bot_1000->Draw("");
    h_MC_stat_Bot_500->Draw("SAME");
    h_MC_stat_Bot_100->Draw("SAME");
    h_MC_stat_Bot_50->Draw("SAME");
    h_MC_stat_Bot_10->Draw("SAME");
    
    //spectrum
    TH1D* h_MC_stat_Bot_10_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC10.root", "MC=10", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_10_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_50_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC50.root", "MC=50", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_50_spec->SetLineColor(kGreen);
    h_MC_stat_Bot_50_spec->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_100_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC100.root", "MC=100", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_100_spec->SetLineColor(kRed);
    h_MC_stat_Bot_100_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_500_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC500.root", "MC=500", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_500_spec->SetLineColor(kBlue);
    h_MC_stat_Bot_500_spec->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_Bot_1000_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC1000.root", "MC=1000", "0-20MeV/BotGe_Mean_calib_systErrors", "0-20MeV/BotGe_Mean_calib_systErrors");
    h_MC_stat_Bot_1000_spec->SetLineColor(kYellow+2);
    h_MC_stat_Bot_1000_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");

    TCanvas* cBot_2 = new TCanvas("Bot spectrum", "Bot spectrum");
    h_MC_stat_Bot_10_spec->Draw("");
    h_MC_stat_Bot_50_spec->Draw("SAME");
    h_MC_stat_Bot_100_spec->Draw("SAME");
    h_MC_stat_Bot_500_spec->Draw("SAME");
    h_MC_stat_Bot_1000_spec->Draw("SAME");
    
    //-- LWO
    cout<<"\nLWO"<<endl;
    TH1D* h_MC_stat_LWO_10 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC10.root", "MC=10", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_10->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_50 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC50.root", "MC=50", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_50->SetLineColor(kGreen);
    h_MC_stat_LWO_50->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_100 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC100.root", "MC=100", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_100->SetLineColor(kRed);
    h_MC_stat_LWO_100->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_500 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC500.root", "MC=500", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_500->SetLineColor(kBlue);
    h_MC_stat_LWO_500->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_1000 = MC_stat_proj("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC1000.root", "MC=1000", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_1000->SetLineColor(kYellow+2);
    h_MC_stat_LWO_1000->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");

    TCanvas* cLWO = new TCanvas("LWO", "LWO");
    h_MC_stat_LWO_1000->Draw("");
    h_MC_stat_LWO_500->Draw("SAME");
    h_MC_stat_LWO_100->Draw("SAME");
    h_MC_stat_LWO_50->Draw("SAME");
    h_MC_stat_LWO_10->Draw("SAME");
    
    //Spectrum
    TH1D* h_MC_stat_LWO_10_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC10.root", "MC=10", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_10_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_50_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC50.root", "MC=50", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_50_spec->SetLineColor(kGreen);
    h_MC_stat_LWO_50_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_100_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC100.root", "MC=100", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_100_spec->SetLineColor(kRed);
    h_MC_stat_LWO_100_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_500_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC500.root", "MC=500", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_500_spec->SetLineColor(kBlue);
    h_MC_stat_LWO_500_spec->GetYaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");
    
    TH1D* h_MC_stat_LWO_1000_spec = MC_stat("../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC5000.root", "../Data_COV/2021_04_29_09_54_36/MC_Tests/2021_04_29_09_54_36_2_MC1000.root", "MC=1000", "0-20MeV/LWO_Mean_calib_systErros", "0-20MeV/LWO_Mean_calib_systErrors");
    h_MC_stat_LWO_1000_spec->SetLineColor(kYellow+2);
    h_MC_stat_LWO_1000_spec->GetXaxis()->SetTitle("(#sigma^{N=5000}_{i} - #sigma^{N=MC}_{i})/#sigma^{N=5000}_{i} [%]");

    TCanvas* cLWO_2 = new TCanvas("LWO spectrum", "LWO spectrum");
//    h_MC_stat_LWO_10_spec->Draw("");
//    h_MC_stat_LWO_50_spec->Draw("SAME");
//    h_MC_stat_LWO_100_spec->Draw("SAME");
    h_MC_stat_LWO_500_spec->Draw("SAME");
    h_MC_stat_LWO_1000_spec->Draw("SAME");

}

//------------ Plot errors
TH1D* errors_h(TH1D* histo, TString Name){
    
    int Nbins = histo->GetNbinsX();
    int Emax = histo->GetBinLowEdge(Nbins)+histo->GetBinWidth(Nbins);
    int Emin = histo->GetBinLowEdge(1);
    TH1D* h_errors = new TH1D(Name, Name, Nbins, Emin, Emax);
    
    double error, content;
    for (int i=0; i<Nbins; i++){
        content= histo->GetBinContent(i+1);
        error= histo->GetBinError(i+1);
        h_errors->SetBinContent(i+1, 100*error/content);
    }
    h_errors->GetYaxis()->SetTitle("#sigma/#N_{i} [%]");
    h_errors->GetXaxis()->SetTitle("Energy [MeV]");
    
    return h_errors;
}

//Read file from Save_Syst_errors and plot errors
void Plot_errors_syst_stat(TString FileinName="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2_MC5000.root"){
    TFile* Filein = new TFile(FileinName, "READ");
    
    //Top
    TH1D* histo_mean = (TH1D*) Filein->Get("0-20MeV/TopGe_Mean_calib_systErrors");
    TH1D* histo_ref = (TH1D*) Filein->Get("0-20MeV/TopGe_Ref_calib_statErrors");
    
    TH1D* herr_mean_Top = errors_h(histo_mean, "SystErrors Top");
    TH1D* herr_ref_Top = errors_h(histo_ref, "StatErrors Top");
    
    herr_ref_Top->SetLineColor(kRed);
    
    //Bot
    histo_mean = (TH1D*) Filein->Get("0-20MeV/BotGe_Mean_calib_systErrors");
    histo_ref = (TH1D*) Filein->Get("0-20MeV/BotGe_Ref_calib_statErrors");
    
    TH1D* herr_mean_Bot = errors_h(histo_mean, "SystErrors Bot");
    TH1D* herr_ref_Bot= errors_h(histo_ref, "StatErrors Bot");
    
    herr_ref_Bot->SetLineColor(kRed);
    
    //LWO
    histo_mean = (TH1D*) Filein->Get("0-20MeV/LWO_Mean_calib_systErros");
    histo_ref = (TH1D*) Filein->Get("0-20MeV/LWO_Ref_calib_statErrors");
    
    TH1D* herr_mean_LWO = errors_h(histo_mean, "SystErrors LWO");
    TH1D* herr_ref_LWO= errors_h(histo_ref, "StatErrors LWO");
    
    herr_ref_LWO->SetLineColor(kRed);
    
    TCanvas* c1 = new TCanvas("Top", "Top");
    herr_mean_Top->Draw();
    herr_ref_Top->Draw("SAME");
    
    TCanvas* c2 = new TCanvas("Bot", "Bot");
    herr_mean_Bot->Draw();
    herr_ref_Bot->Draw("SAME");
    
    TCanvas* cLWO = new TCanvas("LWO", "LWO");
    herr_mean_LWO->Draw();
    herr_ref_LWO->Draw("SAME");
}

//------------------------- Create a gif with N calib and the mean spectrum and the correlation matrix
void Gif_N_calibs_Spectra(int N, TString File="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2.root", int Nbins=4, Double_t Emin= 0, Double_t Emax =40){
    
    TFile* filein = new TFile(File);
    TH1D* histo_0_Bot_Ge = (TH1D*) filein->Get("h1_Peak1");
    TH1D* histo_0_Top_Ge = (TH1D*) filein->Get("h1_Peak2");
    TH1D* histo_0_LWO = (TH1D*) filein->Get("h1_NTD_Amplitude");
    
    vector<TFitResultPtr> Fit_results = PlotEnergyCalibration();
    
    TCanvas* cBot_Ge = new TCanvas();
    cBot_Ge->SetLogy();
    cBot_Ge->Draw();
    
    TCanvas* cTop_Ge = new TCanvas();
    cTop_Ge->SetLogy();
    cTop_Ge->Draw();
    
    TCanvas* cLWO = new TCanvas();
    cLWO->SetLogy();
    cLWO->Draw();
    
    ROOT::Math::GSLRandomEngine* rnd = new ROOT::Math::GSLRandomEngine();
    rnd->Initialize();
    rnd->SetSeed(33566);
    rnd->Initialize();
    
    TRandom3* rnd3 = new TRandom3();
    
    TMatrixD Spectra_Bot_Ge(N, Nbins);
    TMatrixD Spectra_Top_Ge(N, Nbins);
    TMatrixD Spectra_LWO(N, Nbins);
    
    TTree* T;
    TFile *filetree = TFile::Open("filetree.root");
          
    if ((!filetree) || (filetree->IsZombie())){
        TFile* filetree = new TFile("filetree.root", "CREATE");
        T = Tree_from_histos(histo_0_Bot_Ge, histo_0_Top_Ge, histo_0_LWO, 10000000);
    }
    else {
        T=(TTree*) filetree->Get("Events_COV");
    }
    
    cBot_Ge->cd();
    //rnd->Initialize();
    //vector<Double_t> Calib_Bot_Ge = Generate_a_calib(Fit_results[0], rnd); //with 2 variables: use the Gaussian2D function from root
    vector<Double_t> Calib_Bot_Ge = {Fit_results[0]->Parameter(0), Fit_results[0]->Parameter(1), Fit_results[0]->Parameter(2)};
    
    
    // fix b and c;
    double a_BotGe, b_BotGe, c_BotGe;
    a_BotGe=Calib_Bot_Ge[0];
    b_BotGe=Calib_Bot_Ge[1];
    c_BotGe=Calib_Bot_Ge[2];
    
    TH1D* calib_hist_Bot_Ge = RecalibHisto_fromTree(histo_0_Bot_Ge, T, Nbins, Emin, Emax, Calib_Bot_Ge, "Bot Ge", "Bot_Ge");
    //TH1D* calib_hist_Bot_Ge = RecalibHisto(histo_0_Bot_Ge, Nbins, 20, Calib_Bot_Ge[0], Calib_Bot_Ge[1], "histo_Bot_Ge_0");
    //calib_hist_Bot_Ge->SetLineColor(kGray);
    calib_hist_Bot_Ge->Draw("HIST");
    
    int gifcnt =0;
    char gifTop[28];
    char gifBot[28];
    char gifLWO[28];
    sprintf(gifBot, "GIF/N_spectraBot_%.3d.gif", gifcnt);
    cBot_Ge->SaveAs(gifBot);
    
    cTop_Ge->cd();
    //rnd->Initialize();
//    vector<Double_t> Calib_Top_Ge = Generate_a_calib(Fit_results[1], rnd);
    vector<Double_t> Calib_Top_Ge = {Fit_results[1]->Parameter(0), Fit_results[1]->Parameter(1)};
    
    // fix b;
    double a_TopGe, b_TopGe;
    a_TopGe=Calib_Top_Ge[0];
    b_TopGe=Calib_Top_Ge[1];
    
    TH1D* calib_hist_Top_Ge = RecalibHisto_fromTree(histo_0_Top_Ge, T, Nbins, Emin, Emax, Calib_Top_Ge, "Top Ge", "Top_Ge");
//    TH1D* calib_hist_Top_Ge = RecalibHisto(histo_0_Top_Ge, Nbins, 20, Calib_Top_Ge[0], Calib_Top_Ge[1], "histo_Top_Ge_0");
    //calib_hist_Top_Ge->SetLineColor(kGray);
    calib_hist_Top_Ge->Draw("HIST");
    
    sprintf(gifTop, "GIF/N_spectraTop_%.3d.gif", gifcnt);
    cTop_Ge->SaveAs(gifTop);
    
    cLWO->cd();
    //rnd->Initialize();
//    vector<Double_t> Calib_LWO = Generate_a_calib(Fit_results[2], rnd);
    vector<Double_t> Calib_LWO = {Fit_results[2]->Parameter(0), Fit_results[2]->Parameter(1)};
    
    // fix b;
    double a_LWO, b_LWO;
    a_LWO=Calib_LWO[0];
    b_LWO=Calib_LWO[1];
    
    TH1D* calib_hist_LWO = RecalibHisto_fromTree(histo_0_LWO, T, Nbins, Emin, Emax, Calib_LWO, "LWO", "LWO");
    //TH1D* calib_hist_LWO = RecalibHisto(histo_0_LWO, Nbins, 20, Calib_LWO[0], Calib_LWO[1], "histo_LWO_0");
    //calib_hist_LWO->SetLineColor(kGray);
    calib_hist_LWO->Draw("HIST");
    
    sprintf(gifLWO, "GIF/N_spectraLWO_%.3d.gif", gifcnt);
    cLWO->SaveAs(gifLWO);
    
    for (int k=0; k<Nbins; k++){
        Spectra_Bot_Ge(0, k)=calib_hist_Bot_Ge->GetBinContent(k+1);
        Spectra_Top_Ge(0, k)=calib_hist_Top_Ge->GetBinContent(k+1);
        Spectra_LWO(0, k)=calib_hist_LWO->GetBinContent(k+1);
    }
    gStyle->SetOptStat(0);
    
    for (int i=1; i<N; i++){
        cout<<"i = "<<i<<endl;
        
        cBot_Ge->cd();
//        Calib_Bot_Ge = Generate_a_calib(Fit_results[0], rnd);
        //Calib_Bot_Ge = Generate_a_calib_Chol(Fit_results[0], 3, rnd3);
        
        //use fixed b and c
        Calib_Bot_Ge[0]=a_BotGe;
        Calib_Bot_Ge[1]=(1.+0.1*i)*b_BotGe;
        Calib_Bot_Ge[2]=c_BotGe;
        
        calib_hist_Bot_Ge = RecalibHisto_fromTree(histo_0_Bot_Ge, T, Nbins, Emin, Emax, Calib_Bot_Ge, Form("histo_Bot_Ge_%d", i), "Bot_Ge");
//        TH1D* calib_hist_Bot_Ge = RecalibHisto(histo_0_Bot_Ge, Nbins, 20, Calib_Bot_Ge[0], Calib_Bot_Ge[1], Form("histo_Bot_Ge_%d", i));
        //calib_hist_Bot_Ge->SetLineColor(kGray);
        calib_hist_Bot_Ge->Draw("HIST");
        
        gifcnt++;
        sprintf(gifBot, "GIF/N_spectraBot_%.3d.gif", gifcnt);
        cBot_Ge->SaveAs(gifBot);
        
        cTop_Ge->cd();
//        Calib_Top_Ge = Generate_a_calib(Fit_results[1], rnd);
        //Calib_Top_Ge = Generate_a_calib_Chol(Fit_results[1], 2, rnd3);
        
        //use fixed b
        Calib_Top_Ge[0]=(1.+0.1*i)*a_TopGe;
        Calib_Top_Ge[1]=b_TopGe;
        
        calib_hist_Top_Ge = RecalibHisto_fromTree(histo_0_Top_Ge, T, Nbins, Emin, Emax, Calib_Top_Ge, Form("histo_Top_Ge_%d", i), "Top_Ge");
//        TH1D* calib_hist_Top_Ge = RecalibHisto(histo_0_Top_Ge, Nbins, 20, Calib_Top_Ge[0], Calib_Top_Ge[1], Form("histo_Top_Ge_%d", i));
        //calib_hist_Top_Ge->SetLineColor(kGray);
        calib_hist_Top_Ge->Draw("HIST");
        
        sprintf(gifTop, "GIF/N_spectraTop_%.3d.gif", gifcnt);
        cTop_Ge->SaveAs(gifTop);
        
        cLWO->cd();
//        Calib_LWO = Generate_a_calib(Fit_results[2], rnd);
        //Calib_LWO = Generate_a_calib_Chol(Fit_results[2], 2, rnd3);
        
        //use fixed b
        Calib_LWO[0]=(1.+0.1*i)*a_LWO;
        Calib_LWO[1]=b_LWO;
        
        calib_hist_LWO = RecalibHisto_fromTree(histo_0_LWO, T, Nbins, Emin, Emax, Calib_LWO, Form("histo_LWO_%d", i), "LWO");
//        calib_hist_LWO = RecalibHisto(histo_0_LWO, Nbins, 20, Calib_LWO[0], Calib_LWO[1], Form("histo_LWO_%d", i));
        //calib_hist_LWO->SetLineColor(kGray);
        calib_hist_LWO->Draw("HIST");
        
        sprintf(gifLWO, "GIF/N_spectraLWO_%.3d.gif", gifcnt);
        cLWO->SaveAs(gifLWO);
        
        for (int k=0; k<Nbins; k++){
            Spectra_Bot_Ge(i, k)=calib_hist_Bot_Ge->GetBinContent(k+1);
            Spectra_Top_Ge(i, k)=calib_hist_Top_Ge->GetBinContent(k+1);
            Spectra_LWO(i, k)=calib_hist_LWO->GetBinContent(k+1);
        }
    }
    
    TVectorD Means_Bot_Ge(Nbins);
    TVectorD Means_Top_Ge(Nbins);
    TVectorD Means_LWO(Nbins);
    
    for (int j=0; j<Nbins; j++){
        double mean_j_Bot =0;
        double mean_j_Top =0;
        double mean_j_LWO =0;
        for (int i=0; i<N; i++){
            mean_j_Bot+=Spectra_Bot_Ge(i,j);
            mean_j_Top+=Spectra_Top_Ge(i,j);
            mean_j_LWO+=Spectra_LWO(i,j);
        }
        Means_Bot_Ge(j)=mean_j_Bot/N;
        Means_Top_Ge(j)=mean_j_Top/N;
        Means_LWO(j)=mean_j_LWO/N;
    }
    
    TMatrixD Cov_Mat_Top(Nbins, Nbins);
    TMatrixD Cov_Mat_Bot(Nbins, Nbins);
    TMatrixD Cov_Mat_LWO(Nbins, Nbins);
    
    for (int i=0; i<Nbins; i++){
        for (int j=0; j<Nbins; j++){
            Cov_Mat_Bot(i,j)=Cov(i, j, Spectra_Bot_Ge, Means_Bot_Ge(i), Means_Bot_Ge(j), N);
            Cov_Mat_Top(i,j)=Cov(i, j, Spectra_Top_Ge, Means_Top_Ge(i), Means_Top_Ge(j), N);
            Cov_Mat_LWO(i,j)=Cov(i, j, Spectra_LWO, Means_LWO(i), Means_LWO(j), N);
        }
    }
    
    TH1D* hMeans_Bot_Ge = (TH1D*) calib_hist_Bot_Ge->Clone();
    TH1D* hMeans_Top_Ge = (TH1D*) calib_hist_Top_Ge->Clone();
    TH1D* hMeans_LWO = (TH1D*) calib_hist_LWO->Clone();
    
    hMeans_Bot_Ge->SetLineColor(kBlue+2);
    hMeans_Top_Ge->SetLineColor(kBlue+2);
    hMeans_LWO->SetLineColor(kBlue+2);
    
    hMeans_Bot_Ge->SetTitle("Bot_Ge mean");
    hMeans_Top_Ge->SetTitle("Top_Ge mean");
    hMeans_LWO->SetTitle("LWO mean");
    
    for (int j=0; j<Nbins; j++){
        hMeans_Bot_Ge->SetBinContent(j+1, Means_Bot_Ge(j));
        hMeans_Bot_Ge->SetBinError(j+1, sqrt(Cov_Mat_Bot(j,j)));
        hMeans_Top_Ge->SetBinContent(j+1, Means_Top_Ge(j));
        hMeans_Top_Ge->SetBinError(j+1, sqrt(Cov_Mat_Top(j,j)));
        hMeans_LWO->SetBinContent(j+1, Means_LWO(j));
        hMeans_LWO->SetBinError(j+1, sqrt(Cov_Mat_LWO(j,j)));
    }
    
    gifcnt++;
    cBot_Ge->cd();
    hMeans_Bot_Ge->Draw("");
    sprintf(gifBot, "GIF/N_spectraBot_%.3d.gif", gifcnt);
    cBot_Ge->SaveAs(gifBot);
    
    cTop_Ge->cd();
    hMeans_Top_Ge->Draw("");
    sprintf(gifTop, "GIF/N_spectraTop_%.3d.gif", gifcnt);
    cTop_Ge->SaveAs(gifTop);
    
    cLWO->cd();
    hMeans_LWO->Draw("");
    sprintf(gifLWO, "GIF/N_spectraLWO_%.3d.gif", gifcnt);
    cLWO->SaveAs(gifLWO);
    
    cTop_Ge->Close();
    cBot_Ge->Close();
    cLWO->Close();
    
    TH2D* Corr_Mat_Bot = new TH2D("CorrelationMatrix_Bot", "CorrelationMatrix_Bot",
                                  Nbins, hMeans_Bot_Ge->GetBinLowEdge(1), hMeans_Bot_Ge->GetBinLowEdge(Nbins)+hMeans_Bot_Ge->GetBinWidth(1),
                                  Nbins, hMeans_Bot_Ge->GetBinLowEdge(1), hMeans_Bot_Ge->GetBinLowEdge(Nbins)+hMeans_Bot_Ge->GetBinWidth(1));
    Corr_Mat_Bot->GetXaxis()->SetTitle(hMeans_Bot_Ge->GetXaxis()->GetTitle());
    Corr_Mat_Bot->GetYaxis()->SetTitle(hMeans_Bot_Ge->GetXaxis()->GetTitle());
    
    TH2D* Corr_Mat_Top = new TH2D("CorrelationMatrix_Top", "CorrelationMatrix_Top",
                                  Nbins, hMeans_Top_Ge->GetBinLowEdge(1), hMeans_Top_Ge->GetBinLowEdge(Nbins)+hMeans_Top_Ge->GetBinWidth(1),
                                  Nbins, hMeans_Top_Ge->GetBinLowEdge(1), hMeans_Top_Ge->GetBinLowEdge(Nbins)+hMeans_Top_Ge->GetBinWidth(1));
    Corr_Mat_Top->GetXaxis()->SetTitle(hMeans_Top_Ge->GetXaxis()->GetTitle());
    Corr_Mat_Top->GetYaxis()->SetTitle(hMeans_Top_Ge->GetXaxis()->GetTitle());
    
    TH2D* Corr_Mat_LWO = new TH2D("CorrelationMatrix_LWO", "CorrelationMatrix_LWO",
                                  Nbins, hMeans_LWO->GetBinLowEdge(1), hMeans_LWO->GetBinLowEdge(Nbins)+hMeans_LWO->GetBinWidth(1),
                                  Nbins, hMeans_LWO->GetBinLowEdge(1), hMeans_LWO->GetBinLowEdge(Nbins)+hMeans_LWO->GetBinWidth(1));
    Corr_Mat_LWO->GetXaxis()->SetTitle(hMeans_LWO->GetXaxis()->GetTitle());
    Corr_Mat_LWO->GetYaxis()->SetTitle(hMeans_LWO->GetXaxis()->GetTitle());
    
    for (int i=0; i<Nbins; i++){
        for (int j=0; j<Nbins; j++){
            Corr_Mat_Top->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_Top));
            Corr_Mat_Bot->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_Bot));
            Corr_Mat_LWO->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_LWO));
        }
    }
    
    TCanvas* cMat_Bot_Ge = new TCanvas();
    Corr_Mat_Top->Draw("colz TEXT");
    Corr_Mat_Top->GetZaxis()->SetRangeUser(-1, 1);
    
    TCanvas* cMat_Top_Ge = new TCanvas();
    Corr_Mat_Bot->Draw("colz TEXT");
    Corr_Mat_Bot->GetZaxis()->SetRangeUser(-1, 1);
    
    TCanvas* cMat_LWO = new TCanvas();
    Corr_Mat_LWO->Draw("colz TEXT");
    Corr_Mat_LWO->GetZaxis()->SetRangeUser(-1, 1);
    
//    //make an animated gif file using gifsicle.
//    //for information about gifsicle, see http://www.lcdf.org/~eddietwo/gifsicle/
    gSystem->Exec("gifsicle --delay=100 --loop=1000 GIF/N_spectraTop*.gif > GIF/anim_N_spectraTop.gif");
    gSystem->Exec("rm -f GIF/N_spectraTop*.gif");

    gSystem->Exec("gifsicle --delay=100 --loop=1000 GIF/N_spectraBot*.gif > GIF/anim_N_spectraBot.gif");
    gSystem->Exec("rm -f GIF/N_spectraBot*.gif");

    gSystem->Exec("gifsicle --delay=100 --loop=1000 GIF/N_spectraLWO*.gif > GIF/anim_N_spectraLWO.gif");
    gSystem->Exec("rm -f GIF/N_spectraLWO*.gif");
//
    //You can view the animated file anim.gif with Netscape/IE
    //or with gifview as shown below (finish by typing "q" in the window)
    //for more info about gifsicle, gifview, see above url or do
    // gifsicle --help     gifview --help
//    gSystem->Exec("gifview -a GIF/anim_N_spectraTop.gif");
//    gSystem->Exec("gifview -a GIF/anim_N_spectraBot.gif");
//    gSystem->Exec("gifview -a GIF/anim_N_spectraLWO.gif");

}

//------------------------- Plot the Diff in the bins of N calibs and the mean spectrum and the correlation matrix
void DiffBins_N_calibs_Spectra(int N, TString File="../Data_COV/2021_04_29_09_54_36/2021_04_29_09_54_36_2.root", int Nbins=4, Double_t Emin= 0, Double_t Emax =40){
    
    TFile* filein = new TFile(File);
    TH1D* histo_0_Bot_Ge = (TH1D*) filein->Get("h1_Peak1");
    TH1D* histo_0_Top_Ge = (TH1D*) filein->Get("h1_Peak2");
    TH1D* histo_0_LWO = (TH1D*) filein->Get("h1_NTD_Amplitude");
    
    vector<TFitResultPtr> Fit_results = PlotEnergyCalibration();
    
    TCanvas* cBot_Ge = new TCanvas();
    cBot_Ge->Draw();
    
    TCanvas* cTop_Ge = new TCanvas();
    cTop_Ge->Draw();
    
    TCanvas* cLWO = new TCanvas();
    cLWO->Draw();
    
    ROOT::Math::GSLRandomEngine* rnd = new ROOT::Math::GSLRandomEngine();
    rnd->Initialize();
    rnd->SetSeed(33566);
    rnd->Initialize();
    
    TRandom3* rnd3 = new TRandom3();
    
    TMatrixD Spectra_Bot_Ge(N, Nbins);
    TMatrixD Spectra_Top_Ge(N, Nbins);
    TMatrixD Spectra_LWO(N, Nbins);
    
    TTree* T;
    TFile *filetree = TFile::Open("filetree.root");
          
    if ((!filetree) || (filetree->IsZombie())){
        TFile* filetree = new TFile("filetree.root", "CREATE");
        T = Tree_from_histos(histo_0_Bot_Ge, histo_0_Top_Ge, histo_0_LWO, 10000000);
    }
    else {
        T=(TTree*) filetree->Get("Events_COV");
    }
    
    //------------ Bot
    vector<Double_t> Calib_Bot_Ge = {Fit_results[0]->Parameter(0), Fit_results[0]->Parameter(1), Fit_results[0]->Parameter(2)};
    
    // fix b and c;
    double a_BotGe, b_BotGe, c_BotGe;
    a_BotGe=Calib_Bot_Ge[0];
    b_BotGe=Calib_Bot_Ge[1];
    c_BotGe=Calib_Bot_Ge[2];
    
    vector<TH1D*> histos_BotGe, histos_TopGe, histos_LWO;
    TVectorD N0_BotGe(N), N0_TopGe(N), N0_LWO(N);
    
    for (int i=0; i<Nbins; i++){
        histos_BotGe.push_back(new TH1D(Form("Bin_%d_BotGe",i), Form("Bin_%d_BotGe",i), N, 0, N));
        histos_TopGe.push_back(new TH1D(Form("Bin_%d_TopGe",i), Form("Bin_%d_TopGe",i), N, 0, N));
        histos_LWO.push_back(new TH1D(Form("Bin_%d_LWO",i), Form("Bin_%d_LWO",i), N, 0, N));
    }
    
    TH1D* calib_hist_Bot_Ge = RecalibHisto_fromTree(histo_0_Bot_Ge, T, Nbins, Emin, Emax, Calib_Bot_Ge, "Bot Ge", "Bot_Ge");
    calib_hist_Bot_Ge->Draw("HIST");
    
    for (int i=0; i<Nbins; i++){
        N0_BotGe(i)=calib_hist_Bot_Ge->GetBinContent(i+1);
    }
    
    
    //------------ Top
    vector<Double_t> Calib_Top_Ge = {Fit_results[1]->Parameter(0), Fit_results[1]->Parameter(1)};
    
    // fix b;
    double a_TopGe, b_TopGe;
    a_TopGe=Calib_Top_Ge[0];
    b_TopGe=Calib_Top_Ge[1];
    
    TH1D* calib_hist_Top_Ge = RecalibHisto_fromTree(histo_0_Top_Ge, T, Nbins, Emin, Emax, Calib_Top_Ge, "Top Ge", "Top_Ge");
    
    for (int i=0; i<Nbins; i++){
        N0_TopGe(i)=calib_hist_Top_Ge->GetBinContent(i+1);
    }
    
    calib_hist_Top_Ge->Draw("HIST");
    
    //------------ LWO
    vector<Double_t> Calib_LWO = {Fit_results[2]->Parameter(0), Fit_results[2]->Parameter(1)};
    
    // fix b;
    double a_LWO, b_LWO;
    a_LWO=Calib_LWO[0];
    b_LWO=Calib_LWO[1];
    
    TH1D* calib_hist_LWO = RecalibHisto_fromTree(histo_0_LWO, T, Nbins, Emin, Emax, Calib_LWO, "LWO", "LWO");
    
    for (int i=0; i<Nbins; i++){
        N0_LWO(i)=calib_hist_LWO->GetBinContent(i+1);
    }
    
    calib_hist_LWO->Draw("HIST");
    
    for (int k=0; k<Nbins; k++){
        Spectra_Bot_Ge(0, k)=calib_hist_Bot_Ge->GetBinContent(k+1);
        Spectra_Top_Ge(0, k)=calib_hist_Top_Ge->GetBinContent(k+1);
        Spectra_LWO(0, k)=calib_hist_LWO->GetBinContent(k+1);
    }
    
    for (int i=0; i<Nbins; i++){
        histos_BotGe[i]->SetBinContent(1, 0);
        histos_TopGe[i]->SetBinContent(1, 0);
        histos_LWO[i]->SetBinContent(1, 0);
    }
    
    gStyle->SetOptStat(0);
    
    for (int i=1; i<N; i++){
        cout<<"i = "<<i<<endl;
        
        //------------ Bot
        //Calib_Bot_Ge = Generate_a_calib_Chol(Fit_results[0], 3, rnd3);
        
        //use fixed a and b
        Calib_Bot_Ge[0]= a_BotGe;
        Calib_Bot_Ge[1]= (1.+0.1*i)*b_BotGe;
        Calib_Bot_Ge[2]= c_BotGe;
        
        calib_hist_Bot_Ge = RecalibHisto_fromTree(histo_0_Bot_Ge, T, Nbins, Emin, Emax, Calib_Bot_Ge, Form("histo_Bot_Ge_%d", i), "Bot_Ge");
        calib_hist_Bot_Ge->Draw("HIST");
        
        for (int j=0; j<Nbins; j++){
            histos_BotGe[j]->SetBinContent(i+1, N0_BotGe(j)-calib_hist_Bot_Ge->GetBinContent(j+1));
        }
        
        //------------ Top
        //Calib_Top_Ge = Generate_a_calib_Chol(Fit_results[1], 2, rnd3);
        
        //use fixed b
        Calib_Top_Ge[0]=(1.+0.1*i)*a_TopGe;
        Calib_Top_Ge[1]= b_TopGe;
        
        calib_hist_Top_Ge = RecalibHisto_fromTree(histo_0_Top_Ge, T, Nbins, Emin, Emax, Calib_Top_Ge, Form("histo_Top_Ge_%d", i), "Top_Ge");
        calib_hist_Top_Ge->Draw("HIST");

        for (int j=0; j<Nbins; j++){
            histos_TopGe[j]->SetBinContent(i+1, N0_TopGe(j)-calib_hist_Top_Ge->GetBinContent(j+1));
        }
        
        //------------ LWO
        //Calib_LWO = Generate_a_calib_Chol(Fit_results[2], 2, rnd3);
        
        //use fixed b
        Calib_LWO[0]=(1.+0.1*i)*a_LWO;
        Calib_LWO[1]=b_LWO;
        
        calib_hist_LWO = RecalibHisto_fromTree(histo_0_LWO, T, Nbins, Emin, Emax, Calib_LWO, Form("histo_LWO_%d", i), "LWO");
        calib_hist_LWO->Draw("HIST");
        
        for (int j=0; j<Nbins; j++){
            histos_LWO[j]->SetBinContent(i+1, N0_LWO(j)-calib_hist_LWO->GetBinContent(j+1));
        }
        
        for (int k=0; k<Nbins; k++){
            Spectra_Bot_Ge(i, k)=calib_hist_Bot_Ge->GetBinContent(k+1);
            Spectra_Top_Ge(i, k)=calib_hist_Top_Ge->GetBinContent(k+1);
            Spectra_LWO(i, k)=calib_hist_LWO->GetBinContent(k+1);
        }
    }
    
    TVectorD Means_Bot_Ge(Nbins);
    TVectorD Means_Top_Ge(Nbins);
    TVectorD Means_LWO(Nbins);
    
    for (int j=0; j<Nbins; j++){
        double mean_j_Bot =0;
        double mean_j_Top =0;
        double mean_j_LWO =0;
        
        for (int i=0; i<N; i++){
            mean_j_Bot+=Spectra_Bot_Ge(i,j);
            mean_j_Top+=Spectra_Top_Ge(i,j);
            mean_j_LWO+=Spectra_LWO(i,j);
        }
        
        Means_Bot_Ge(j)=mean_j_Bot/N;
        Means_Top_Ge(j)=mean_j_Top/N;
        Means_LWO(j)=mean_j_LWO/N;
    }
    
    TMatrixD Cov_Mat_Top(Nbins, Nbins);
    TMatrixD Cov_Mat_Bot(Nbins, Nbins);
    TMatrixD Cov_Mat_LWO(Nbins, Nbins);
    
    for (int i=0; i<Nbins; i++){
        for (int j=0; j<Nbins; j++){
            Cov_Mat_Bot(i,j)=Cov(i, j, Spectra_Bot_Ge, Means_Bot_Ge(i), Means_Bot_Ge(j), N);
            Cov_Mat_Top(i,j)=Cov(i, j, Spectra_Top_Ge, Means_Top_Ge(i), Means_Top_Ge(j), N);
            Cov_Mat_LWO(i,j)=Cov(i, j, Spectra_LWO, Means_LWO(i), Means_LWO(j), N);
        }
    }
    
    TH1D* hMeans_Bot_Ge = (TH1D*) calib_hist_Bot_Ge->Clone();
    TH1D* hMeans_Top_Ge = (TH1D*) calib_hist_Top_Ge->Clone();
    TH1D* hMeans_LWO = (TH1D*) calib_hist_LWO->Clone();
    
    hMeans_Bot_Ge->SetLineColor(kBlue+2);
    hMeans_Top_Ge->SetLineColor(kBlue+2);
    hMeans_LWO->SetLineColor(kBlue+2);
    
    hMeans_Bot_Ge->SetTitle("Bot_Ge mean");
    hMeans_Top_Ge->SetTitle("Top_Ge mean");
    hMeans_LWO->SetTitle("LWO mean");
    
    for (int j=0; j<Nbins; j++){
        hMeans_Bot_Ge->SetBinContent(j+1, Means_Bot_Ge(j));
        hMeans_Bot_Ge->SetBinError(j+1, sqrt(Cov_Mat_Bot(j,j)));
        
        hMeans_Top_Ge->SetBinContent(j+1, Means_Top_Ge(j));
        hMeans_Top_Ge->SetBinError(j+1, sqrt(Cov_Mat_Top(j,j)));
        
        hMeans_LWO->SetBinContent(j+1, Means_LWO(j));
        hMeans_LWO->SetBinError(j+1, sqrt(Cov_Mat_LWO(j,j)));
    }
    
    cBot_Ge->cd();
    histos_BotGe[0]->Draw("HIST");
    histos_BotGe[0]->GetXaxis()->SetTitle("Draw Number i");
    histos_BotGe[0]->GetYaxis()->SetTitle("(N_{0}-N_{i})");
    histos_BotGe[0]->GetYaxis()->SetRangeUser(-2, 2);
    
    cTop_Ge->cd();
    histos_TopGe[0]->Draw("HIST");
    histos_TopGe[0]->GetXaxis()->SetTitle("Draw Number i");
    histos_TopGe[0]->GetYaxis()->SetTitle("(N_{0}-N_{i})");
    histos_TopGe[0]->GetYaxis()->SetRangeUser(-2, 2);
    
    cLWO->cd();
    histos_LWO[0]->Draw("HIST");
    histos_LWO[0]->GetXaxis()->SetTitle("Draw Number i");
    histos_LWO[0]->GetYaxis()->SetTitle("(N_{0}-N_{i})");
    histos_LWO[0]->GetYaxis()->SetRangeUser(-2, 2);
    
    for (int i=1; i<Nbins; i++){
        cBot_Ge->cd();
        histos_BotGe[i]->SetLineColor(i);
        histos_BotGe[i]->Draw("HIST SAME");
        
        cTop_Ge->cd();
        histos_TopGe[i]->SetLineColor(i);
        histos_TopGe[i]->Draw("HIST SAME");
        
        cLWO->cd();
        histos_LWO[i]->SetLineColor(i);
        histos_LWO[i]->Draw("HIST SAME");
    }
    
//    cTop_Ge->Close();
//    cBot_Ge->Close();
//    cLWO->Close();
    
    TH2D* Corr_Mat_Bot = new TH2D("CorrelationMatrix_Bot", "CorrelationMatrix_Bot",
                                  Nbins, hMeans_Bot_Ge->GetBinLowEdge(1), hMeans_Bot_Ge->GetBinLowEdge(Nbins)+hMeans_Bot_Ge->GetBinWidth(1),
                                  Nbins, hMeans_Bot_Ge->GetBinLowEdge(1), hMeans_Bot_Ge->GetBinLowEdge(Nbins)+hMeans_Bot_Ge->GetBinWidth(1));
    Corr_Mat_Bot->GetXaxis()->SetTitle(hMeans_Bot_Ge->GetXaxis()->GetTitle());
    Corr_Mat_Bot->GetYaxis()->SetTitle(hMeans_Bot_Ge->GetXaxis()->GetTitle());
    
    TH2D* Corr_Mat_Top = new TH2D("CorrelationMatrix_Top", "CorrelationMatrix_Top",
                                  Nbins, hMeans_Top_Ge->GetBinLowEdge(1), hMeans_Top_Ge->GetBinLowEdge(Nbins)+hMeans_Top_Ge->GetBinWidth(1),
                                  Nbins, hMeans_Top_Ge->GetBinLowEdge(1), hMeans_Top_Ge->GetBinLowEdge(Nbins)+hMeans_Top_Ge->GetBinWidth(1));
    Corr_Mat_Top->GetXaxis()->SetTitle(hMeans_Top_Ge->GetXaxis()->GetTitle());
    Corr_Mat_Top->GetYaxis()->SetTitle(hMeans_Top_Ge->GetXaxis()->GetTitle());
    
    TH2D* Corr_Mat_LWO = new TH2D("CorrelationMatrix_LWO", "CorrelationMatrix_LWO",
                                  Nbins, hMeans_LWO->GetBinLowEdge(1), hMeans_LWO->GetBinLowEdge(Nbins)+hMeans_LWO->GetBinWidth(1),
                                  Nbins, hMeans_LWO->GetBinLowEdge(1), hMeans_LWO->GetBinLowEdge(Nbins)+hMeans_LWO->GetBinWidth(1));
    Corr_Mat_LWO->GetXaxis()->SetTitle(hMeans_LWO->GetXaxis()->GetTitle());
    Corr_Mat_LWO->GetYaxis()->SetTitle(hMeans_LWO->GetXaxis()->GetTitle());
    
    for (int i=0; i<Nbins; i++){
        for (int j=0; j<Nbins; j++){
            Corr_Mat_Top->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_Top));
            Corr_Mat_Bot->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_Bot));
            Corr_Mat_LWO->SetBinContent(i+1, j+1, Corr(i, j, Cov_Mat_LWO));
        }
    }
    
    TCanvas* cMat_Bot_Ge = new TCanvas();
    Corr_Mat_Top->Draw("colz TEXT");
    Corr_Mat_Top->GetZaxis()->SetRangeUser(-1, 1);
    
    TCanvas* cMat_Top_Ge = new TCanvas();
    Corr_Mat_Bot->Draw("colz TEXT");
    Corr_Mat_Bot->GetZaxis()->SetRangeUser(-1, 1);
    
    TCanvas* cMat_LWO = new TCanvas();
    Corr_Mat_LWO->Draw("colz TEXT");
    Corr_Mat_LWO->GetZaxis()->SetRangeUser(-1, 1);

}

void Plot_random_coeff(int ErrorFact=1) {
    
    vector<TFitResultPtr> Fit_results = PlotEnergyCalibration();
    
    ROOT::Math::GSLRandomEngine rnd;
    rnd.Initialize();
    
    Double_t Center_b[3];
    Double_t Center_a[3];
    
    Double_t Sigma_b[3];
    Double_t Sigma_a[3];
    Double_t Corr_ab[3];
    
    TMatrixDSym CovMat(2);
    TVectorD Means(2);
    TVectorD Sigma(2);
    
    std::vector<TMatrixDSym> Cryst_CovMat;
    std::vector<TVectorD> Cryst_Means;
    std::vector<TVectorD> Cryst_Sigma;
    
    
    for (int i=0; i<3; i++){
        Center_b[i] = Fit_results[i]->Parameter(0);
        Center_a[i] = Fit_results[i]->Parameter(1);
        CovMat = Fit_results[i]->GetCovarianceMatrix();
        CovMat*=ErrorFact;
        
        Cryst_CovMat.push_back(CovMat);
        
        Sigma_b[i] = sqrt(CovMat(0,0));
        Sigma_a[i] = sqrt(CovMat(1,1));
        
        Means(0)=Center_b[i];
        Means(1)=Center_a[i];
        
        Cryst_Means.push_back(Means);
        
        Sigma(0)=Sigma_b[i];
        Sigma(1)=Sigma_a[i];
        
        Cryst_Sigma.push_back(Sigma);
        
        Corr_ab[i] = CovMat(1,0)/(Sigma_b[i]*Sigma_a[i]);
        
    }
    
    TH1F* ha_Bot_Ge = new TH1F("Slope Bot Ge", "Slope Bot Ge", 100, Center_a[0]-5*Sigma_a[0], Center_a[0]+5*Sigma_a[0]);
    TH1F* hb_Bot_Ge = new TH1F("Origin Bot Ge", "Origin Bot Ge", 100, Center_b[0]-5*Sigma_b[0], Center_b[0]+5*Sigma_b[0]);
    TH2F* hab_Bot_Ge = new TH2F("Slope vs Origin Bot Ge", "Slope vs Origin Bot Ge", 100,  Center_a[0]-5*Sigma_a[0], Center_a[0]+5*Sigma_a[0], 100, Center_b[0]-5*Sigma_b[0], Center_b[0]+5*Sigma_b[0]);
    ha_Bot_Ge->GetXaxis()->SetTitle("Slope [V/keV]");
    ha_Bot_Ge->GetXaxis()->SetTitleSize(0.06);
    ha_Bot_Ge->GetXaxis()->SetLabelSize(0.06);
    ha_Bot_Ge->GetYaxis()->SetTitleSize(0.06);
    ha_Bot_Ge->GetYaxis()->SetLabelSize(0.06);
    
    hb_Bot_Ge->GetXaxis()->SetTitle("Origin [V]");
    hb_Bot_Ge->GetXaxis()->SetTitleSize(0.06);
    hb_Bot_Ge->GetXaxis()->SetLabelSize(0.06);
    hb_Bot_Ge->GetYaxis()->SetTitleSize(0.06);
    hb_Bot_Ge->GetYaxis()->SetLabelSize(0.06);
    
    hab_Bot_Ge->GetXaxis()->SetTitle("Slope [V/keV]");
    hab_Bot_Ge->GetYaxis()->SetTitle("Origin [V]");
    hab_Bot_Ge->GetXaxis()->SetTitleSize(0.06);
    hab_Bot_Ge->GetXaxis()->SetLabelSize(0.06);
    hab_Bot_Ge->GetYaxis()->SetTitleSize(0.06);
    hab_Bot_Ge->GetYaxis()->SetLabelSize(0.06);

    TH1F* ha_Top_Ge = new TH1F("Slope Top Ge", "Slope Top Ge", 100, Center_a[1]-5*Sigma_a[1], Center_a[1]+5*Sigma_a[1]);
    TH1F* hb_Top_Ge = new TH1F("Origin Top Ge", "Origin Top Ge", 100, Center_b[1]-5*Sigma_b[1], Center_b[1]-5*Sigma_b[1]);
    TH2F* hab_Top_Ge = new TH2F("Slope vs Origin Top Ge", "Slope vs Origin Top Ge", 100, Center_a[1]-5*Sigma_a[1], Center_a[1]+5*Sigma_a[1], 100, Center_b[1]-5*Sigma_b[1], Center_b[1]-5*Sigma_b[1]);
    ha_Top_Ge->GetXaxis()->SetTitle("Slope [V/keV]");
    ha_Top_Ge->GetXaxis()->SetTitleSize(0.06);
    ha_Top_Ge->GetXaxis()->SetLabelSize(0.06);
    ha_Top_Ge->GetYaxis()->SetTitleSize(0.06);
    ha_Top_Ge->GetYaxis()->SetLabelSize(0.06);
    
    hb_Top_Ge->GetXaxis()->SetTitle("Origin [V]");
    hb_Top_Ge->GetXaxis()->SetTitleSize(0.06);
    hb_Top_Ge->GetXaxis()->SetLabelSize(0.06);
    hb_Top_Ge->GetYaxis()->SetTitleSize(0.06);
    hb_Top_Ge->GetYaxis()->SetLabelSize(0.06);
    
    hab_Top_Ge->GetXaxis()->SetTitle("Slope [V/keV]");
    hab_Top_Ge->GetYaxis()->SetTitle("Origin [V]");
    hab_Top_Ge->GetXaxis()->SetTitleSize(0.06);
    hab_Top_Ge->GetXaxis()->SetLabelSize(0.06);
    hab_Top_Ge->GetYaxis()->SetTitleSize(0.06);
    hab_Top_Ge->GetYaxis()->SetLabelSize(0.06);

    TH1F* ha_LWO = new TH1F("Slope LWO", "Slope LWO", 100, Center_a[2]-5*Sigma_a[2], Center_a[2]+5*Sigma_a[2]);
    TH1F* hb_LWO = new TH1F("Origin LWO", "Origin LWO", 100, Center_b[2]-5*Sigma_b[2], Center_b[2]-5*Sigma_b[2]);
    TH2F* hab_LWO = new TH2F("Slope vs Origin LWO", "Slope vs Origin LWO", 100, Center_a[2]-5*Sigma_a[2], Center_a[2]+5*Sigma_a[2], 100, Center_b[2]-5*Sigma_b[2], Center_b[2]-5*Sigma_b[2]);
    
    ha_LWO->GetXaxis()->SetTitle("Slope [V/keV]");
    ha_LWO->GetXaxis()->SetTitleSize(0.06);
    ha_LWO->GetXaxis()->SetLabelSize(0.06);
    ha_LWO->GetYaxis()->SetTitleSize(0.06);
    ha_LWO->GetYaxis()->SetLabelSize(0.06);
    
    hb_LWO->GetXaxis()->SetTitle("Origin [V]");
    hb_LWO->GetXaxis()->SetTitleSize(0.06);
    hb_LWO->GetXaxis()->SetLabelSize(0.06);
    hb_LWO->GetYaxis()->SetTitleSize(0.06);
    hb_LWO->GetYaxis()->SetLabelSize(0.06);
    
    hab_LWO->GetXaxis()->SetTitle("Slope [V/keV]");
    hab_LWO->GetYaxis()->SetTitle("Origin [V]");
    hab_LWO->GetXaxis()->SetTitleSize(0.06);
    hab_LWO->GetXaxis()->SetLabelSize(0.06);
    hab_LWO->GetYaxis()->SetTitleSize(0.06);
    hab_LWO->GetYaxis()->SetLabelSize(0.06);
    
    const int MAX = 10000;
    Double_t x, y;
    TVectorD corr_var_Bot(2);
    TVectorD corr_var_Top(2);
    TVectorD corr_var_LWO(2);
    TRandom3* rand_gen = new TRandom3();
    
    Cryst_Sigma[0].Print();
    for (int evnts = 0; evnts < MAX; ++evnts) {
        
        //rnd.Gaussian2D(Sigma_a[0], Sigma_b[0], Corr_ab[0], x, y);
        corr_var_Bot = Generate_n_correlated_variables(Cryst_CovMat[0], rand_gen);
        ha_Bot_Ge->Fill(Center_a[0]+corr_var_Bot(1));
        hb_Bot_Ge->Fill(Center_b[0]+corr_var_Bot(0));
        hab_Bot_Ge->Fill(Center_a[0]+corr_var_Bot(1), Center_b[0]+corr_var_Bot(0));
//        ha_Bot_Ge->Fill(x+Center_a[0]);
//        hb_Bot_Ge->Fill(y+Center_b[0]);
//        hab_Bot_Ge->Fill(x+Center_a[0], y+Center_b[0]);
        
//        rnd.Gaussian2D(Sigma_a[1], Sigma_b[1], Corr_ab[1], x, y);
        corr_var_Top = Generate_n_correlated_variables(Cryst_CovMat[1], rand_gen);
        ha_Top_Ge->Fill(Center_a[1]+corr_var_Top(1));
        hb_Top_Ge->Fill(Center_b[1]+corr_var_Top(0));
        hab_Top_Ge->Fill(Center_a[1]+corr_var_Top(1), Center_b[1]+corr_var_Top(0));
//        ha_Top_Ge->Fill(x+Center_a[1]);
//        hb_Top_Ge->Fill(y+Center_b[1]);
//        hab_Top_Ge->Fill(x+Center_a[1], y+Center_b[1]);
        
//        rnd.Gaussian2D(Sigma_a[2], Sigma_b[2], Corr_ab[2], x, y);
        corr_var_LWO = Generate_n_correlated_variables(Cryst_CovMat[2], rand_gen);
        ha_LWO->Fill(Center_a[2]+corr_var_LWO(1));
        hb_LWO->Fill(Center_b[2]+corr_var_LWO(0));
        hab_LWO->Fill(Center_a[2]+corr_var_LWO(1), Center_b[2]+corr_var_LWO(0));
//        ha_LWO->Fill(x+Center_a[2]);
//        hb_LWO->Fill(y+Center_b[2]);
//        hab_LWO->Fill(x+Center_a[2], y+Center_b[2]);

    }

    TCanvas* c = new TCanvas("c", "Multivariate gaussian random numbers");
    c->Divide(3, 3);
    c->cd(1);
    ha_Bot_Ge->Draw();
    c->cd(2);
    hb_Bot_Ge->Draw();
    c->cd(3);
    hab_Bot_Ge->Draw("COL");
    
    c->cd(4);
    ha_Top_Ge->Draw();
    c->cd(5);
    hb_Top_Ge->Draw();
    c->cd(6);
    hab_Top_Ge->Draw("COL");
    
    c->cd(7);
    ha_LWO->Draw();
    c->cd(8);
    hb_LWO->Draw();
    c->cd(9);
    hab_LWO->Draw("COL");
}

//------------------------- Plot N calib curves
void Plot_N_calibs(int N = 100){
    
    double Emax=20000;
    
    vector<TFitResultPtr> Fit_results = PlotEnergyCalibration();
    
    ROOT::Math::GSLRandomEngine* rnd = new ROOT::Math::GSLRandomEngine();
    
    rnd->Initialize();

    TCanvas* c = new TCanvas();
    
    vector<Double_t> Calib_Bot_Ge = Generate_a_calib(Fit_results[0], rnd);
    TF1* Bot_Ge_curve =  new TF1("Bot_Ge_curve_0","([0]+[1]*x)",0,Emax);
    Bot_Ge_curve->SetParNames("b [#V]","a [#V/keV]");
    Bot_Ge_curve->SetParameter(0, Calib_Bot_Ge[0]);
    Bot_Ge_curve->SetParameter(1, Calib_Bot_Ge[1]);
    Bot_Ge_curve->SetLineColor(kRed);
    Bot_Ge_curve->Draw("");
    
    vector<Double_t> Calib_Top_Ge = Generate_a_calib(Fit_results[1], rnd);
    TF1* Top_Ge_curve =  new TF1("Top_Ge_curve_0","([0]+[1]*x)",0,Emax);
    Top_Ge_curve->SetParNames("b [#V]","a [#V/keV]");
    Top_Ge_curve->SetParameter(0, Calib_Top_Ge[0]);
    Top_Ge_curve->SetParameter(1, Calib_Top_Ge[1]);
    Top_Ge_curve->SetLineColor(kBlue);
    Top_Ge_curve->Draw("SAME");

    vector<Double_t> Calib_LWO = Generate_a_calib(Fit_results[2], rnd);
    TF1* LWO_curve =  new TF1("LWO_curve_0","([0]+[1]*x)",0,Emax);
    LWO_curve->SetParNames("b [#V]","a [#V/keV]");
    LWO_curve->SetParameter(0, Calib_LWO[0]);
    LWO_curve->SetParameter(1, Calib_LWO[1]);
    LWO_curve->SetLineColor(kGreen);
    LWO_curve->Draw("SAME");
    
    for (int i=1; i<N; i++){
        cout<<"i = "<<i<<endl;
        vector<Double_t> Calib_Bot_Ge = Generate_a_calib(Fit_results[0], rnd);
        TF1* Bot_Ge_curve =  new TF1(Form("Bot_Ge_curve_%d", i),"([0]+[1]*x)",0,Emax);
        Bot_Ge_curve->SetParNames("b [#V]","a [#V/keV]");
        Bot_Ge_curve->SetParameter(0, Calib_Bot_Ge[0]);
        Bot_Ge_curve->SetParameter(1, Calib_Bot_Ge[1]);
        Bot_Ge_curve->SetLineColor(kRed+i/10);
        Bot_Ge_curve->Draw("SAME");

        vector<Double_t> Calib_Top_Ge = Generate_a_calib(Fit_results[1], rnd);
        TF1* Top_Ge_curve =  new TF1(Form("Top_Ge_curve_%d", i),"([0]+[1]*x)",0,Emax);
        Top_Ge_curve->SetParNames("b [#V]","a [#V/keV]");
        Top_Ge_curve->SetParameter(0, Calib_Top_Ge[0]);
        Top_Ge_curve->SetParameter(1, Calib_Top_Ge[1]);
        Top_Ge_curve->SetLineColor(kBlue+i/10);
        Top_Ge_curve->Draw("SAME");

        vector<Double_t> Calib_LWO = Generate_a_calib(Fit_results[2], rnd);
        TF1* LWO_curve =  new TF1(Form("LWO_curve_%d", i),"([0]+[1]*x)",0,Emax);
        LWO_curve->SetParNames("b [#V]","a [#V/keV]");
        LWO_curve->SetParameter(0, Calib_LWO[0]);
        LWO_curve->SetParameter(1, Calib_LWO[1]);
        LWO_curve->SetLineColor(kGreen+i/10);
        LWO_curve->Draw("SAME");
    }
}
