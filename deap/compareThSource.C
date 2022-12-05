{

  TFile *f1 = new TFile("Combined5.root");// Thorium source
  TFile *f2 = new TFile("Combined7.root");

  TTree *tree1 = (TTree*)f1->Get("data_satCorr");
  TTree *tree2 = (TTree*)f2->Get("data_satCorr");

  TH2F *H2_nSCBayes_rprompt60Bayes_1 = new TH2F( "H2_nSCBayes_rprompt60Bayes_1", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0);
  TH2F *H2_qpe_fprompt_1             = new TH2F( "H2_qpe_fprompt_1"            , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0);

  TH2F *H2_nSCBayes_rprompt60Bayes_2 = new TH2F( "H2_nSCBayes_rprompt60Bayes_2", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0);
  TH2F *H2_qpe_fprompt_2             = new TH2F( "H2_qpe_fprompt_2"            , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0);

  tree1->Project("H2_nSCBayes_rprompt60Bayes_1", "rprompt60Bayes:nSCBayes","");
  tree1->Project("H2_qpe_fprompt_1", "fprompt:qPE","");

  tree2->Project("H2_nSCBayes_rprompt60Bayes_2", "rprompt60Bayes:nSCBayes","");
  tree2->Project("H2_qpe_fprompt_2", "fprompt:qPE","");


  double liveTime = 1522431.54; // daniel: 5269088.0; //6.09848148148148184e+01
  double liveTime2 = 1369187.29;
  TCanvas *c1 = new TCanvas("c1","nSCBayes vs rprompt60",1600,600);
  c1->Divide(2,1);
  c1->cd(1);  gPad->SetLogz();
  H2_nSCBayes_rprompt60Bayes_1->Draw("colz");
  c1->cd(2);  gPad->SetLogz();
  H2_nSCBayes_rprompt60Bayes_2->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","qpe vs fprompt",1600,600);
  c2->Divide(2,1);
  c2->cd(1);  gPad->SetLogz();
  H2_qpe_fprompt_1->Draw("colz");
  c2->cd(2);  gPad->SetLogz();
  H2_qpe_fprompt_2->Draw("colz");

  TCanvas *c3 = new TCanvas("c3", "", 1000,1000);
  c3->Divide(2,2); 
  c3->cd(1);//gPad->SetLogy();
  int biny1 = H2_nSCBayes_rprompt60Bayes_1->GetYaxis()->FindBin(0.35);
  int biny2 = H2_nSCBayes_rprompt60Bayes_1->GetYaxis()->FindBin(0.45);
  int biny1vac = H2_nSCBayes_rprompt60Bayes_2->GetYaxis()->FindBin(0.35);
  int biny2vac = H2_nSCBayes_rprompt60Bayes_2->GetYaxis()->FindBin(0.45);
  hnSCBayes_0p35_0p45 = H2_nSCBayes_rprompt60Bayes_1->ProjectionX("hnSCBayes_0p35_0p45",biny1, biny2);
  hnSCBayes_0p35_0p45_vac = H2_nSCBayes_rprompt60Bayes_2->ProjectionX("hnSCBayes_0p35_0p45_vac", biny1vac, biny2vac);
  hnSCBayes_0p35_0p45->Scale(1./liveTime);
  hnSCBayes_0p35_0p45_vac->Scale(1./liveTime2);
  hnSCBayes_0p35_0p45->Draw();
  hnSCBayes_0p35_0p45_vac->SetLineColor(kRed);
  hnSCBayes_0p35_0p45->SetTitle("Integral fprompt from 0.35 to 0.45");
  hnSCBayes_0p35_0p45_vac->Draw("same");
  TLegend *legend1= new TLegend(0.1,0.7,0.2,0.9);
  legend1->AddEntry(hnSCBayes_0p35_0p45,"combine 5 (17.62 live-day)","l");
  legend1->AddEntry(hnSCBayes_0p35_0p45_vac,"combine 7 run (15.85 live-day)","l");
  legend1->Draw("same");

  c3->cd(2);//gPad->SetLogy();
  int biny11 = H2_nSCBayes_rprompt60Bayes_1->GetYaxis()->FindBin(0.5);
  int biny21 = H2_nSCBayes_rprompt60Bayes_1->GetYaxis()->FindBin(0.7);
  int biny11vac = H2_nSCBayes_rprompt60Bayes_2->GetYaxis()->FindBin(0.5);
  int biny21vac = H2_nSCBayes_rprompt60Bayes_2->GetYaxis()->FindBin(0.7);
  hnSCBayes_0p5_0p7 = H2_nSCBayes_rprompt60Bayes_1->ProjectionX("hnSCBayes_0p5_0p7",biny11, biny21);
  hnSCBayes_0p5_0p7_vac = H2_nSCBayes_rprompt60Bayes_2->ProjectionX("hnSCBayes_0p5_0p7_vac", biny11vac, biny21vac);
  hnSCBayes_0p5_0p7->Scale(1./liveTime);
  hnSCBayes_0p5_0p7_vac->Scale(1./liveTime2);
  hnSCBayes_0p5_0p7->Draw();
  hnSCBayes_0p5_0p7_vac->SetLineColor(kRed);
  hnSCBayes_0p5_0p7->SetTitle("Integral fprompt from 0.5 to 0.7");
  hnSCBayes_0p5_0p7_vac->Draw("same");
  TLegend *legend2= new TLegend(0.1,0.7,0.2,0.9);
  legend2->AddEntry(hnSCBayes_0p5_0p7,"combine 5 (17.62 live-day)","l");
  legend2->AddEntry(hnSCBayes_0p5_0p7_vac,"combine 7 run (15.85 live-day)","l");
  legend2->Draw("same");

  c3->cd(3);//gPad->SetLogy();
  int biny13 = H2_nSCBayes_rprompt60Bayes_1->GetYaxis()->FindBin(0.7);
  int biny23 = H2_nSCBayes_rprompt60Bayes_1->GetYaxis()->FindBin(0.9);
  int biny13vac = H2_nSCBayes_rprompt60Bayes_2->GetYaxis()->FindBin(0.7);
  int biny23vac = H2_nSCBayes_rprompt60Bayes_2->GetYaxis()->FindBin(0.9);
  hnSCBayes_0p7_0p9 = H2_nSCBayes_rprompt60Bayes_1->ProjectionX("hnSCBayes_0p7_0p9",biny13, biny23);
  hnSCBayes_0p7_0p9_vac = H2_nSCBayes_rprompt60Bayes_2->ProjectionX("hnSCBayes_0p7_0p9_vac", biny13vac, biny23vac);
  hnSCBayes_0p7_0p9->Scale(1./liveTime);
  hnSCBayes_0p7_0p9_vac->Scale(1./liveTime2);
  hnSCBayes_0p7_0p9_vac->Draw();
  hnSCBayes_0p7_0p9_vac->SetLineColor(kRed);
  hnSCBayes_0p7_0p9_vac->SetTitle("Integral fprompt from 0.7 to 0.9");
  hnSCBayes_0p7_0p9->Draw("same");
  TLegend *legend3 = new TLegend(0.1,0.7,0.2,0.9);
  legend3->AddEntry(hnSCBayes_0p7_0p9,"combine 5 (17.62 live-day)","l");
  legend3->AddEntry(hnSCBayes_0p7_0p9_vac,"combine 7 run (15.85 live-day)","l");
  legend3->Draw("same");

  c3->cd(4);//gPad->SetLogy();
  hnSCBayes = (TH1D*)H2_nSCBayes_rprompt60Bayes_1->ProjectionX("hnSCBayes",0,100);
  hnSCBayes_vac = (TH1D*)H2_nSCBayes_rprompt60Bayes_2->ProjectionX("hnSCBayes_vac",10,100);
  hnSCBayes->Scale(1./liveTime);
  hnSCBayes_vac->Scale(1./liveTime2);
  hnSCBayes_vac->SetLineColor(kRed);
  hnSCBayes->SetTitle("whole");
  hnSCBayes->Draw();
  hnSCBayes_vac->Draw("same");
  TLegend *legend4= new TLegend(0.1,0.7,0.2,0.9);
  legend4->AddEntry(hnSCBayes,"combine 5 (17.62 live-day)","l");
  legend4->AddEntry(hnSCBayes_vac,"combine 7 run (15.85 live-day)","l");
  legend4->Draw("same");


}
