{
  TFile *f1 = new TFile("ResolMPnew_PartialScintNeck_FillTl208_Scint_r1_s0_p1.root");
  TH2F *h2 = (TH2F*)f1->Get("hfitRZ_trig");
  gStyle->SetPalette(57);
  h2->GetXaxis()->SetRangeUser(0,6000);
  h2->GetYaxis()->SetRangeUser(3000,9000);
  h2->Draw("colz");

}
