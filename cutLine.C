{
  TCanvas *c1 = new TCanvas("c1");
  TH1F h("h","test",100,-4,4);
  h.FillRandom("gaus",10000);
  h->Draw();
  Float_t ymax = h.GetMaximum();
  TLine *line = new TLine(-3,ymax,3,ymax);
  line->SetLineColor(kRed);
  line->Draw();
}
