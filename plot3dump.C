{
  TFile *f1 = new TFile("multipath_scintwater_gtID_0_2scintpdf_pmtID19.root");
  TFile *f2 = new TFile("multipath_scintwater_gtID_0_2waterpdf_pmtID19.root");
  TFile *f3 = new TFile("multipath_scintwater_gtID_0_scintwater_pmtID19.root");

  hx1 = (TH1F*)f1->Get("Hx project to PMT");
  hy1 = (TH1F*)f1->Get("Hy project to PMT");
  hz1 = (TH1F*)f1->Get("Hz project to PMT");

  hx2 = (TH1F*)f2->Get("Hx project to PMT");
  hy2 = (TH1F*)f2->Get("Hy project to PMT");
  hz2 = (TH1F*)f2->Get("Hz project to PMT");

  hx3 = (TH1F*)f3->Get("Hx project to PMT");
  hy3 = (TH1F*)f3->Get("Hy project to PMT");
  hz3 = (TH1F*)f3->Get("Hz project to PMT");


  hdx1 = (TH1F*)f1->Get("DHx project to PMT");
  hdy1 = (TH1F*)f1->Get("DHy project to PMT");
  hdz1 = (TH1F*)f1->Get("DHz project to PMT");

  hdx2 = (TH1F*)f2->Get("DHx project to PMT");
  hdy2 = (TH1F*)f2->Get("DHy project to PMT");
  hdz2 = (TH1F*)f2->Get("DHz project to PMT");

  hdx3 = (TH1F*)f3->Get("DHx project to PMT");
  hdy3 = (TH1F*)f3->Get("DHy project to PMT");
  hdz3 = (TH1F*)f3->Get("DHz project to PMT");

  hxall1 = (TH1F*)f1->Get("Hx_2_0_px");
  hyall1 = (TH1F*)f1->Get("Hy_2_0_px");
  hzall1 = (TH1F*)f1->Get("Hz_2_0_px");

  hxall2 = (TH1F*)f2->Get("Hx_2_0_px");
  hyall2 = (TH1F*)f2->Get("Hy_2_0_px");
  hzall2 = (TH1F*)f2->Get("Hz_2_0_px");

  hxall3 = (TH1F*)f3->Get("Hx_2_0_px");
  hyall3 = (TH1F*)f3->Get("Hy_2_0_px");
  hzall3 = (TH1F*)f3->Get("Hz_2_0_px");


  hx1->SetLineColor(kRed); hy1->SetLineColor(kRed); hz1->SetLineColor(kRed);
  hx2->SetLineColor(kBlue);hy2->SetLineColor(kBlue);hz2->SetLineColor(kBlue);
  hx3->SetLineColor(kBlack);hy3->SetLineColor(kBlack);hz3->SetLineColor(kBlack);

  hdx1->SetLineColor(kRed); hdy1->SetLineColor(kRed); hdz1->SetLineColor(kRed);
  hdx2->SetLineColor(kBlue);hdy2->SetLineColor(kBlue);hdz2->SetLineColor(kBlue);
  hdx3->SetLineColor(kBlack);hdy3->SetLineColor(kBlack);hdz3->SetLineColor(kBlack);

  hxall1->SetLineColor(kRed); hyall1->SetLineColor(kRed); hzall1->SetLineColor(kRed);
  hxall2->SetLineColor(kBlue);hyall2->SetLineColor(kBlue);hzall2->SetLineColor(kBlue);
  hxall3->SetLineColor(kBlack);hyall3->SetLineColor(kBlack);hzall3->SetLineColor(kBlack);

  TCanvas c("c","one PMT, LLH",800,600);
  c.Divide(2,2);
  c.cd(1);hx1->GetYaxis()->SetRangeUser(0,15);hx1->Draw();hx2->Draw("sames");hx3->Draw("sames");
  c.cd(2);hy1->GetYaxis()->SetRangeUser(0,15);hy1->Draw();hy2->Draw("sames");hy3->Draw("sames");
  c.cd(3);hz1->GetYaxis()->SetRangeUser(0,15);hz1->Draw();hz2->Draw("sames");hz3->Draw("sames");

  TCanvas c1("c1","one PMT, Derivative",800,600);
  c1.Divide(2,2);
  c1.cd(1);hdx1->GetYaxis()->SetRangeUser(-0.1,0.1);hdx1->Draw();hdx2->Draw("sames");hdx3->Draw("sames");
  c1.cd(2);hdy1->GetYaxis()->SetRangeUser(-0.1,0.1);hdy1->Draw();hdy2->Draw("sames");hdy3->Draw("sames");
  c1.cd(3);hdz1->GetYaxis()->SetRangeUser(-0.1,0.1);hdz1->Draw();hdz2->Draw("sames");hdz3->Draw("sames");

  TCanvas c2("c2","all PMT, LLH",800,600);
  c2.Divide(2,2);
  c2.cd(1);hxall1->GetYaxis()->SetRangeUser(0,20000);hxall1->Draw();hxall2->Draw("sames");hxall3->Draw("sames");
  c2.cd(2);hyall1->GetYaxis()->SetRangeUser(0,20000);hyall1->Draw();hyall2->Draw("sames");hyall3->Draw("sames");
  c2.cd(3);hzall1->GetYaxis()->SetRangeUser(0,20000);hzall1->Draw();hzall2->Draw("sames");hzall3->Draw("sames");

}
