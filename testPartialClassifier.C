{
  TFile *ff = new TFile("PartialScintNeck_FillTl208_Scint_r1_s0_p3.ntuple.root");
  TH2F* h2 = new TH2F("h2", "MP fitted R vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);
  TH1F* hDeltaZ = new TH1F("hDeltaZ", "MP fitted Z - mc Z, FECD cut", 2000,-9000,9000);

  TH2F* hDeltaXvsSkyShine = new TH2F("hDeltaXvsSkyShine", "MP fitted X - mc X vs SkyShine", 2000,-9000,9000,500,0,5);
  TH2F* hDeltaYvsSkyShine = new TH2F("hDeltaYvsSkyShine", "MP fitted Y - mc Y vs SkyShine", 2000,-9000,9000,500,0,5);
  TH2F* hDeltaZvsSkyShine = new TH2F("hDeltaZvsSkyShine", "MP fitted Z - mc Z vs SkyShine", 2000,-9000,9000,500,0,5);

  TH2F* hDeltaZvsNearAV = new TH2F("hDeltaZvsNearAV", "MP fitted Z - mc Z vs NearAV", 2000,-9000,9000,2000,-100000,100);

  gStyle->SetPalette(57);
  TTree* ou = (TTree*)ff->Get("output");
  h2->GetXaxis()->SetRangeUser(0,6000);
  h2->GetYaxis()->SetRangeUser(3000,9000);
//  output->Project("h2","posz_MultiPath:sqrt(posx_MultiPath*posx_MultiPath+posy_MultiPath*posy_MultiPath)","nhits>0");
//  h2->Draw("colz");
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  output->Project("hDeltaXvsSkyShine","skyShine:posx-mcPosx","nhits>0");
  hDeltaXvsSkyShine->Draw("colz");
  c1->cd(2);
  output->Project("hDeltaYvsSkyShine","skyShine:posy-mcPosy","nhits>0");
  hDeltaYvsSkyShine->Draw("colz");
  c1->cd(3); 
  output->Project("hDeltaZvsSkyShine","skyShine:posz-mcPosz","nhits>0");
  hDeltaZvsSkyShine->Draw("colz");

  TFile *fnew = new TFile("testClassifier.root","recreate");
  fnew->cd();
  hDeltaXvsSkyShine->Write();
  hDeltaYvsSkyShine->Write();
  hDeltaZvsSkyShine->Write();
  fnew->Close();

//   output->Project("hDeltaZvsNearAV","nearAV:posz-mcPosz","nhits>0");
//   hDeltaZvsNearAV->Draw("colz");

//  output->Project("hDeltaZ","posz-mcPosz","nhits>0 && nearAV!=-99999");
//  hDeltaZ->Draw();
//  hDeltaZ->Fit("gaus");
}
