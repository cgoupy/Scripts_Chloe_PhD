TH1D* Shift(TH1D* histo, double Shift_factor){
    TString Name = histo->GetName();
    TH1D* histo_shifted = new TH1D(Name+"_shifted", Name+"_shifted", histo->GetNbinsX(), 0, histo->GetBinWidth(histo->GetNbinsX())+histo->GetBinLowEdge(histo->GetNbinsX()));
    histo_shifted->SetTitle(Form("%s_shifted_%.3f", histo->GetTitle(), Shift_factor));
    histo_shifted->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
    histo_shifted->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());
    
    for (int i=1; i<=histo_shifted->GetNbinsX(); i++){
        histo_shifted->SetBinContent(i, 0);
    }
    
    double E_i;
    double E_shifted;
    
    for (int i=0; i<DRAWNUM; i++){
        E_i = histo->GetRandom();
        E_shifted = E_i*Shift_factor;
        histo_shifted->Fill(E_shifted);
    }

    //histo_shifted->Draw();
    histo_shifted->Scale(histo->Integral()/histo_shifted->Integral());
    //cout<<Name<<endl;
    
    for (int i=1; i<histo_shifted->GetNbinsX(); i++){
        histo_shifted->SetBinError(i, histo->GetBinError(histo->FindBin(histo_shifted->GetBinCenter(i)/Shift_factor)));
    }
    return histo_shifted;
}
