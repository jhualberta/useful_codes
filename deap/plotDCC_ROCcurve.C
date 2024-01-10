{
    TFile *fileAr40_qpe = new TFile("outputCut_Merged525runs_CheckDccMC_qPE_ar40_pmt_204_OFF.root"); 
    TFile *fileU232_qpe = new TFile("outputCut_MergedU232_CheckDcc_qpe_cut90to200.root");
    
    TFile *fileAr40_qnscb = new TFile("outputCut_Merged525runs_CheckDccMC_qnscb_ar40_pmt_204_OFF.root");   
    TFile *fileU232_qnscb = new TFile("outputCut_MergedU232_CheckDcc_qnscb_cut90to200.root");

    HAr40_pmtFmaxpe_qpe = (TH1F*)fileAr40_qpe->Get("HpmtFmaxpe");
    HAr40_pmtFmaxpe_qpeMBR800 = (TH1F*)fileAr40_qpe->Get("HpmtFmaxpeMBR800");
    HAr40_pmtFmaxpe_qpeMBR720 = (TH1F*)fileAr40_qpe->Get("HpmtFmaxpeMBR720");
                                                    
    HAr40_pmtFmaxpe_qPrompt = (TH1F*)fileAr40_qpe->Get("HpmtFmaxpePrompt");
    HAr40_pmtFmaxpe_qPromptMBR720 = (TH1F*)fileAr40_qpe->Get("HpmtFmaxpePromptMBR720");
    HAr40_pmtFmaxpe_qPromptMBR800 = (TH1F*)fileAr40_qpe->Get("HpmtFmaxpePromptMBR800");
   
    HAr40_pmtFmaxpe_qnscb = (TH1F*)fileAr40_qnscb->Get("HpmtFmaxpe");
    HAr40_pmtFmaxpe_qnscbMBR800 = (TH1F*)fileAr40_qnscb->Get("HpmtFmaxpeMBR800");
    HAr40_pmtFmaxpe_qnscbMBR720 = (TH1F*)fileAr40_qnscb->Get("HpmtFmaxpeMBR720");

    HU232_pmtFmaxpe_qpe = (TH1F*)fileU232_qpe->Get("HpmtFmaxpe");
    HU232_pmtFmaxpe_qpeMBR800 = (TH1F*)fileU232_qpe->Get("HpmtFmaxpeMBR800");
    HU232_pmtFmaxpe_qpeMBR720 = (TH1F*)fileU232_qpe->Get("HpmtFmaxpeMBR720");
 
    HU232_pmtFmaxpe_qPrompt = (TH1F*)fileU232_qpe->Get("HpmtFmaxpePrompt");
    HU232_pmtFmaxpe_qPromptMBR800 = (TH1F*)fileU232_qpe->Get("HpmtFmaxpePromptMBR800");
    HU232_pmtFmaxpe_qPromptMBR720 = (TH1F*)fileU232_qpe->Get("HpmtFmaxpePromptMBR720");

    HU232_pmtFmaxpe_qnscb = (TH1F*)fileU232_qnscb->Get("HpmtFmaxpe");
    HU232_pmtFmaxpe_qnscbMBR800 = (TH1F*)fileU232_qnscb->Get("HpmtFmaxpeMBR800");
    HU232_pmtFmaxpe_qnscbMBR720 = (TH1F*)fileU232_qnscb->Get("HpmtFmaxpeMBR720");

//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qpe->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qpe->Clone();

    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qPrompt->Clone();
    HU232test = (TH1F*)HU232_pmtFmaxpe_qPrompt->Clone();

//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qpeMBR800->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qpeMBR800->Clone();

//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qpeMBR720->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qpeMBR720->Clone();
//
//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qPrompt->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qPrompt->Clone();
//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qnscb->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qnscb->Clone();

//
//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qPromptMBR800->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qPromptMBR800->Clone();

//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qPromptMBR720->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qPromptMBR720->Clone();
//
//    HAr40test = (TH1F*)HAr40_pmtFmaxpe_qPromptMBR800->Clone();
//    HU232test = (TH1F*)HU232_pmtFmaxpe_qPromptMBR800->Clone();

    HAr40test->Scale(1./HAr40test->Integral());
    HU232test->Scale(1./HU232test->Integral());

    HAr40test->SetLineColor(kRed);
    HU232test->SetLineColor(kBlue);

    double cut = 0.01;
    double step = 0.01;

    int binMax_ar40 = HAr40test->GetMaximumBin();
    int binMax_u232 = HU232test->GetMaximumBin();

    vector<double> sigEff; // Es = 1 - integral from cut to the last bin, how many signal we will keep 
    vector<double> bkgEff; // Eb = integral from cut to the last bin, how many background we will cut

    int npoints = 0;

    startCut = 0.04;
    for(int i = 0; i<14; i++)
    {
      cut = startCut+i*step;
      int start_bin_ar40 = HAr40test->GetXaxis()->FindBin(cut);
      int start_bin_u232 = HU232test->GetXaxis()->FindBin(cut);
      //cout<<start_bin_ar40<< " "<<start_bin_u232<<endl;
      double count_ar40 = HAr40test->Integral(start_bin_ar40, binMax_ar40);
      double count_u232 = HU232test->Integral(start_bin_u232, binMax_u232);
      double sig = 1 - count_ar40;
      double bkg = count_u232;
      sigEff.push_back(sig*100);
      bkgEff.push_back(bkg*100);
      cout<<"step "<<i<<" cut>"<<cut<<" sig="<<(1-count_ar40)*100<<", bkg="<<count_u232*100<<endl;
      npoints++;
    }

    //cout<<npoints<<" steps."<<endl;
    const int nsize = const_cast<const int&>(npoints); 
    
    double array_sig[30];
    double array_bkg[30];
    for(int i = 0; i<nsize; i++)
    {
      array_sig[i] = sigEff[i];
      array_bkg[i] = bkgEff[i];
    }

    TCanvas *c1 = new TCanvas("c1","",800,600);
    c1->cd();
    TGraph *roc = new TGraph(npoints, array_sig, array_bkg);
    roc->Draw("APL*");
    roc->GetYaxis()->SetTitle("background efficiency(%)");
    roc->GetXaxis()->SetTitle("signal efficiency(%)");

    for(int i = 0; i<npoints; i++)
    {
       TLatex *lnum = new TLatex(array_sig[i], array_bkg[i],Form("%.2f",startCut+i*0.01));
       lnum->SetTextFont(12);
       lnum->SetTextSize(0.03);
       lnum->SetTextAlign(12);
       lnum->SetTextAngle(45);
       lnum->SetTextColor(kRed);
       lnum->Draw("same");
    }
 
    TCanvas *c2 = new TCanvas("c2","",800,600);
    c2->cd();
    HAr40test->Draw();
    HU232test->Draw("sames");


//    Double_t *xpeaks=NULL, *ypeaks=NULL;
//  TPolyMarker *pm=(TPolyMarker *)
//    h-->GetListOfFunctions()->FindObject("TPolyMarker");
//  pm->SetMarkerStyle(32);
//  pm->SetMarkerColor(kGreen);
//  pm->SetMarkerSize(0.4);
//  xpeaks=s-->GetPositionX();
//  ypeaks=s-->GetPositionY();

}
