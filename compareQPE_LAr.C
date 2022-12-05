{

  TFile *f1 = new TFile("Merged_Daniel_deap_ntp_018831.root");// Thorium source
  TFile *f2 = new TFile("Merged_Daniel_deap_ntp_022394.root");

  TH2F *H2_qpe_fprompt_0 = (TH2F*)f1->Get("H2_qpe_fprompt_0");
  TH2F *H2_qpe_fprompt_0_vac = (TH2F*)f2->Get("H2_qpe_fprompt_0");
  double liveTime = 74506.8; // daniel: 5269088.0; //6.09848148148148184e+01
  double liveTimeVac = 103383;
  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);
  c1->cd(1);  gPad->SetLogz();
  H2_qpe_fprompt_0->Draw("colz");
  c1->cd(2);  gPad->SetLogz();
  H2_qpe_fprompt_0_vac->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","qpe vs fprompt",1000,1000);
  c2->Divide(2,2);
  /// alpha region
  c2->cd(1);//gPad->SetLogy();
  int biny1 = H2_qpe_fprompt_0->GetYaxis()->FindBin(0.2);
  int biny2 = H2_qpe_fprompt_0->GetYaxis()->FindBin(0.35);
  int biny1vac = H2_qpe_fprompt_0_vac->GetYaxis()->FindBin(0.2);
  int biny2vac = H2_qpe_fprompt_0_vac->GetYaxis()->FindBin(0.35);
  hQPE_0p2_0p35 = H2_qpe_fprompt_0->ProjectionX("hQPE_0p2_0p35",biny1, biny2);
  hQPE_0p2_0p35_vac = H2_qpe_fprompt_0_vac->ProjectionX("hQPE_0p2_0p35_vac", biny1vac, biny2vac);
  hQPE_0p2_0p35->Scale(1./liveTime);
  hQPE_0p2_0p35_vac->Scale(1./liveTimeVac);
  hQPE_0p2_0p35->Draw();
  hQPE_0p2_0p35_vac->SetLineColor(kRed);
  hQPE_0p2_0p35->SetTitle("Integral fprompt from 0.35 to 0.45");
  hQPE_0p2_0p35_vac->Draw("same");
  TLegend *legend1= new TLegend(0.1,0.7,0.2,0.9);
  legend1->AddEntry(hQPE_0p2_0p35,"LAr physics","l");
  legend1->AddEntry(hQPE_0p2_0p35_vac,"ThSource LAr run","l");
  legend1->Draw("same");

  c2->cd(2);//gPad->SetLogy();
  int biny11 = H2_qpe_fprompt_0->GetYaxis()->FindBin(0.5);
  int biny21 = H2_qpe_fprompt_0->GetYaxis()->FindBin(0.7);
  int biny11vac = H2_qpe_fprompt_0_vac->GetYaxis()->FindBin(0.5);
  int biny21vac = H2_qpe_fprompt_0_vac->GetYaxis()->FindBin(0.7);
  hQPE_0p5_0p7 = H2_qpe_fprompt_0->ProjectionX("hQPE_0p5_0p7",biny11, biny21);
  hQPE_0p5_0p7_vac = H2_qpe_fprompt_0_vac->ProjectionX("hQPE_0p5_0p7_vac", biny11vac, biny21vac);
  hQPE_0p5_0p7->Scale(1./liveTime);
  hQPE_0p5_0p7_vac->Scale(1./liveTimeVac);
  hQPE_0p5_0p7->Draw();
  hQPE_0p5_0p7_vac->SetLineColor(kRed);
  hQPE_0p5_0p7->SetTitle("Integral fprompt from 0.5 to 0.7");
  hQPE_0p5_0p7_vac->Draw("same");
  TLegend *legend2= new TLegend(0.1,0.7,0.2,0.9);
  legend2->AddEntry(hQPE_0p5_0p7,"LAr physics","l");
  legend2->AddEntry(hQPE_0p5_0p7_vac,"ThSource LAr run","l");
  legend2->Draw("same");

  c2->cd(3);//gPad->SetLogy();
  int biny13 = H2_qpe_fprompt_0->GetYaxis()->FindBin(0.8);
  int biny23 = H2_qpe_fprompt_0->GetYaxis()->FindBin(1.0);
  int biny13vac = H2_qpe_fprompt_0_vac->GetYaxis()->FindBin(0.8);
  int biny23vac = H2_qpe_fprompt_0_vac->GetYaxis()->FindBin(1.0);
  hQPE_0p8_1p0 = H2_qpe_fprompt_0->ProjectionX("hQPE_0p8_1p0",biny13, biny23);
  hQPE_0p8_1p0_vac = H2_qpe_fprompt_0_vac->ProjectionX("hQPE_0p8_1p0_vac", biny13vac, biny23vac);
  hQPE_0p8_1p0->Scale(1./liveTime);
  hQPE_0p8_1p0_vac->Scale(1./liveTimeVac);
  hQPE_0p8_1p0_vac->Draw();
  hQPE_0p8_1p0_vac->SetLineColor(kRed);
  hQPE_0p8_1p0_vac->SetTitle("Integral fprompt from 0.7 to 0.9");
  hQPE_0p8_1p0->Draw("same");
  TLegend *legend3 = new TLegend(0.1,0.7,0.2,0.9);
  legend3->AddEntry(hQPE_0p8_1p0,"LAr physics","l");
  legend3->AddEntry(hQPE_0p8_1p0_vac,"ThSource LAr run","l");
  legend3->Draw("same");

  c2->cd(4);//gPad->SetLogy();
  hQPE = (TH1D*)H2_qpe_fprompt_0->ProjectionX("hQPE",0,100);
  hQPE_vac = (TH1D*)H2_qpe_fprompt_0_vac->ProjectionX("hQPE_vac",10,100);
  hQPE->Scale(1./liveTime);
  hQPE_vac->Scale(1./liveTimeVac);
  hQPE_vac->SetLineColor(kRed);
  hQPE->SetTitle("whole");
  hQPE->Draw();
  hQPE_vac->Draw("same");
  TLegend *legend4= new TLegend(0.1,0.7,0.2,0.9);
  legend4->AddEntry(hQPE,"LAr physics","l");
  legend4->AddEntry(hQPE_vac,"ThSource LAr run","l");
  legend4->Draw("same");


}
