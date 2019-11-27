{

  gStyle->SetPalette(kBird);
  TFile *f1 = new TFile("ClassifierMP_FitPartialMP_newFit_1p42_leadingEdge2_fillbottom.root");
  TFile *f2 = new TFile("ResolMPnew_FitMP_2p5MeVbeta_filltop_level4400.root");
  TFile *f3 = new TFile("ResolMPnew_FitMP_2p5MeVbeta_fullfill_level4400.root");

  hDeltaX_top = (TH1F*)f1->Get("hDeltaX_trig");
  hDeltaX_bot = (TH1F*)f2->Get("hDeltaX_trig");
  hDeltaX_full = (TH1F*)f3->Get("hDeltaX_trig");

  hRZ_top = (TH2F*)f1->Get("hfitRZ_trig");
  hRZ_bot = (TH2F*)f2->Get("hfitRZ_trig");
  hRZ_full = (TH2F*)f3->Get("hfitRZ_trig");

  hRZ_top->GetXaxis()->SetTitle("#sqrt{x^{2}+y^{2}} [mm]");
  hRZ_bot->GetXaxis()->SetTitle("#sqrt{x^{2}+y^{2}} [mm]");
  hRZ_full->GetXaxis()->SetTitle("#sqrt{x^{2}+y^{2}} [mm]");

  hRZ_top->GetYaxis()->SetTitle("z [mm]");
  hRZ_bot->GetYaxis()->SetTitle("z [mm]");
  hRZ_full->GetYaxis()->SetTitle("z [mm]");


  hDeltaX_top->GetXaxis()->SetTitle("x_{fit}-x_{MC} [mm]");
  hDeltaX_bot->GetXaxis()->SetTitle("x_{fit}-x_{MC} [mm]");
  hDeltaX_full->GetXaxis()->SetTitle("x_{fit}-x_{MC} [mm]");

  hDeltaX_top->GetYaxis()->SetTitle("Number of events");
  hDeltaX_bot->GetYaxis()->SetTitle("Number of events");
  hDeltaX_full->GetYaxis()->SetTitle("Number of events");

  hDeltaX_top->GetXaxis()->SetRangeUser(-3000,3000);
  hDeltaX_bot->GetXaxis()->SetRangeUser(-3000,3000);
  hDeltaX_full->GetXaxis()->SetRangeUser(-3000,3000);

  TCanvas c1("c1","",800,600);
  c1.SetGrid();
//  hDeltaX_top->Draw();
  hRZ_top->Draw("colz");

  TCanvas c2("c2","",800,600);
  c2.SetGrid();
//  hDeltaX_bot->Draw();
  hRZ_bot->Draw("colz");

  TCanvas c3("c3","",800,600);
  c3.SetGrid();
//  hDeltaX_full->Draw();
  hRZ_full->Draw("colz");

}
