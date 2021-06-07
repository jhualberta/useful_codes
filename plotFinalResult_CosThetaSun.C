{
//  gStyle->SetTitleH(0.5);	
//  gStyle->SetTitleFont(50, "xy");
  
  TFile *fdata = new TFile("Merged_E4to15_Merged_MP_p1mask_ExtractTres_WaterMP6176general_nhit20_Analysis_r200004to207718_p004.root");
  TTree *T = (TTree*)fdata->Get("T2");
  
  TH1F *hcos_total_E4to15 = new TH1F("hcos_total_E4to15","",40,-1,1);
  T->Project("hcos_total_E4to15","cosThetaToSun","energy>4 && energy<15");
  TH1F *hcos_default_E4to15 = new TH1F("hcos_default_E4to15","",40,-1,1);
  T->Project("hcos_default_E4to15","cosThetaToSun","posRad<5500 && energy>4 && energy<15 && scaleLogL>10 && zfactor>-11 && zfactor<1 && Gtest<1.9 && Gtest>0 && Utest<0.95");
  hcos_default_E4to15->GetXaxis()->SetTitle("cos#theta_{sun}");
  hcos_default_E4to15->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");

  TH1F *hcos_total_E5to15 = new TH1F("hcos_total_E5to15","",40,-1,1);
  T->Project("hcos_total_E5to15","cosThetaToSun","energy>5 && energy<15");
  TH1F *hcos_default_E5to15 = new TH1F("hcos_default_E5to15","",40,-1,1);
  T->Project("hcos_default_E5to15","cosThetaToSun","posRad<5500 && energy>5 && energy<15 && scaleLogL>10 && zfactor>-11 && zfactor<1 && Gtest<1.9 && Gtest>0 && Utest<0.95");
  hcos_default_E5to15->GetXaxis()->SetTitle("cos#theta_{sun}");
  hcos_default_E5to15->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");

  TH1F *hcos_total_E6to15 = new TH1F("hcos_total_E6to15","",40,-1,1);
  T->Project("hcos_total_E6to15","cosThetaToSun","energy>6 && energy<15");
  TH1F *hcos_default_E6to15 = new TH1F("hcos_default_E6to15","",40,-1,1);
  T->Project("hcos_default_E6to15","cosThetaToSun","posRad<5500 && energy>6 && energy<15 && scaleLogL>10 && zfactor>-11 && zfactor<1 && Gtest<1.9 && Gtest>0 && Utest<0.95");
  hcos_default_E6to15->GetXaxis()->SetTitle("cos#theta_{sun}");
  hcos_default_E6to15->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");

  TFile *ff0 = new TFile("analysis_E4to15.root");
  TH1F *hbdtE4to15_sig = (TH1F*)ff0->Get("hBDTSigCosThetaToSun");
  TH1F *hbdtE4to15_bkg = (TH1F*)ff0->Get("hBDTBkgCosThetaToSun");
  TH1F *hMLPE4to15_sig = (TH1F*)ff0->Get("hMLPSigCosThetaToSun");
  TH1F *hMLPE4to15_bkg = (TH1F*)ff0->Get("hMLPBkgCosThetaToSun");
  hbdtE4to15_sig->SetName("hbdtE4to15_sig");
  hbdtE4to15_bkg->SetName("hbdtE4to15_bkg");
  hMLPE4to15_sig->SetName("hmlpE4to15_sig");
  hMLPE4to15_bkg->SetName("hmlpE4to15_bkg");

  hbdtE4to15_sig->GetXaxis()->SetTitle("cos#theta_{sun}");
  hbdtE4to15_sig->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");
  hbdtE4to15_sig->GetXaxis()->SetTitleOffset(0.8);
  hbdtE4to15_sig->GetXaxis()->SetTitleSize(0.05);
  hbdtE4to15_sig->GetYaxis()->SetTitleOffset(0.8);
  hbdtE4to15_sig->GetYaxis()->SetTitleSize(0.05);

  hbdtE4to15_bkg->SetLineColor(kRed);
  hbdtE4to15_bkg->SetLineStyle(2);

  hbdtE4to15_sig->SetLineWidth(2);
  hbdtE4to15_bkg->SetLineWidth(2);

  hMLPE4to15_sig->GetXaxis()->SetTitle("cos#theta_{sun}");
  hMLPE4to15_sig->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");
  hMLPE4to15_sig->GetXaxis()->SetTitleOffset(0.8);
  hMLPE4to15_sig->GetXaxis()->SetTitleSize(0.05);
  hMLPE4to15_sig->GetYaxis()->SetTitleOffset(0.8);
  hMLPE4to15_sig->GetYaxis()->SetTitleSize(0.05);

  hMLPE4to15_bkg->SetLineColor(kRed);
  hMLPE4to15_bkg->SetLineStyle(2);

  hMLPE4to15_sig->SetLineWidth(2);
  hMLPE4to15_bkg->SetLineWidth(2);

  TFile *ff = new TFile("analysis_E5to15.root");
  TH1F *hbdt1 = (TH1F*)ff->Get("hBDTSigCosThetaToSun");
  TH1F *hbdt2 = (TH1F*)ff->Get("hBDTBkgCosThetaToSun");
  TH1F *hMLP1 = (TH1F*)ff->Get("hMLPSigCosThetaToSun");
  TH1F *hMLP2 = (TH1F*)ff->Get("hMLPBkgCosThetaToSun");
  hbdt1->SetName("hbdtE5to15_sig");
  hbdt2->SetName("hbdtE5to15_bkg");
  hMLP1->SetName("hmlpE5to15_sig");
  hMLP2->SetName("hmlpE5to15_bkg");

  hbdt1->GetXaxis()->SetTitle("cos#theta_{sun}");
  hbdt1->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");
  hbdt1->GetXaxis()->SetTitleOffset(0.8);
  hbdt1->GetXaxis()->SetTitleSize(0.05);
  hbdt1->GetYaxis()->SetTitleOffset(0.8);
  hbdt1->GetYaxis()->SetTitleSize(0.05);

  hbdt2->SetLineColor(kRed);
  hbdt2->SetLineStyle(2);

  hbdt1->SetLineWidth(2);
  hbdt2->SetLineWidth(2);

  hMLP1->GetXaxis()->SetTitle("cos#theta_{sun}");
  hMLP1->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");
  hMLP1->GetXaxis()->SetTitleOffset(0.8);
  hMLP1->GetXaxis()->SetTitleSize(0.05);
  hMLP1->GetYaxis()->SetTitleOffset(0.8);
  hMLP1->GetYaxis()->SetTitleSize(0.05);

  hMLP2->SetLineColor(kRed);
  hMLP2->SetLineStyle(2);

  hMLP1->SetLineWidth(2);
  hMLP2->SetLineWidth(2);

  TFile *foutput = new TFile("final_solarResults.root","recreate");
  foutput->cd();
  hbdt1->Write();hbdt2->Write();
  hMLP1->Write();hMLP2->Write();
  hbdtE4to15_sig->Write(); hbdtE4to15_bkg->Write();
  hMLPE4to15_sig->Write(); hMLPE4to15_bkg->Write();

  hcos_total_E4to15->Write();
  hcos_total_E5to15->Write(); 

  hcos_default_E4to15->Write();
  hcos_default_E5to15->Write();
  hcos_default_E6to15->Write();

  TCanvas c("c","",800,600);
//  hMLP1->Draw();
//  hMLP2->Draw("same");

//  hMLPE4to15_sig->Draw();
//  hMLPE4to15_bkg->Draw("same");

  hcos_default_E5to15->Draw();
  hcos_default_E5to15->SetLineStyle(2);
  hMLP1->SetLineColor(kRed);
  hbdt1->Draw("same");
  hMLP1->Draw("same");

}
