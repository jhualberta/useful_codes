{
  TFile *ff = new TFile("SaveTest_JeffMC_107055.root");
  TH1F *hz = (TH1F*)ff->Get("hfitZ");
  hz->Rebin(4);
  TH1F *hz_cut = new TH1F("hz_cut","",200,-4000,4000);
  for(int i = 0;i<200;i++) 
  {
    hz_cut->SetBinContent(i+1, hz->GetBinContent(hz->GetXaxis()->FindBin(-4000)+i));
  }
  hz_cut->Draw();
}
