// Extract early hit PDF and smooth it

{
  TFile * f = new TFile("PMT_selector_rat.root");
  //hpdf = (TH1F*)f->Get("htRes_early");

  hpdf = (TH1F*)f->Get("htRes_early");

  TH1F *hcheck = new TH1F("hcheck","h",1600,-100,300);

  hpdf->Smooth(2);

  for(int i = 0;i<1600;i++)
  {
   double cont = hpdf->GetBinContent(i+1);
   cout<<cont<<", ";
   hcheck->SetBinContent(i+1, cont);
  }
  cout<<endl;
  hcheck->Draw();
}
