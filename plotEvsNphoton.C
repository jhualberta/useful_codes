{

  TFile *f1 = new TFile("FitMultiWater_1p40_Water_2MeV_10000evt_ISOFill.root");
  TFile *f2 = new TFile("FitMultiWater_1p40_Water_3MeV_10000evt_ISOFill.root");
  TFile *f3 = new TFile("FitMultiWater_1p40_Water_4MeV_10000evt_ISOFill.root");
  TFile *f4 = new TFile("FitMultiWater_1p40_Water_5MeV_10000evt_ISOFill.root");
  TFile *f5 = new TFile("FitMultiWater_1p40_Water_6MeV_10000evt_ISOFill.root");
  TFile *f6 = new TFile("FitMultiWater_1p40_Water_7MeV_10000evt_ISOFill.root");
  TFile *f7 = new TFile("FitMultiWater_1p40_Water_8MeV_10000evt_ISOFill.root");
  TFile *f8 = new TFile("FitMultiWater_1p40_Water_9MeV_10000evt_ISOFill.root");
  TFile *f9 = new TFile("FitMultiWater_1p40_Water_10MeV_10000evt_ISOFill.root");

  TH1F *hNCherPhoton2MeV = new TH1F("hNCherPhoton2MeV","",240,0,6000); //25
  TH1F *hNCherPhoton3MeV = new TH1F("hNCherPhoton3MeV","",240,0,6000); //25
  TH1F *hNCherPhoton4MeV = new TH1F("hNCherPhoton4MeV","",240,0,6000); //25
  TH1F *hNCherPhoton5MeV = new TH1F("hNCherPhoton5MeV","",240,0,6000); //25
  TH1F *hNCherPhoton6MeV = new TH1F("hNCherPhoton6MeV","",240,0,6000); //25
  TH1F *hNCherPhoton7MeV = new TH1F("hNCherPhoton7MeV","",240,0,6000); //25
  TH1F *hNCherPhoton8MeV = new TH1F("hNCherPhoton8MeV","",240,0,6000); //25
  TH1F *hNCherPhoton9MeV = new TH1F("hNCherPhoton9MeV","",240,0,6000); //25
  TH1F *hNCherPhoton10MeV = new TH1F("hNCherPhoton10MeV","",240,0,6000); //25


  TTree* t1 = (TTree*)f1->Get("T");
  TTree* t2 = (TTree*)f2->Get("T");
  TTree* t3 = (TTree*)f3->Get("T");
  TTree* t4 = (TTree*)f4->Get("T");
  TTree* t5 = (TTree*)f5->Get("T");
  TTree* t6 = (TTree*)f6->Get("T");
  TTree* t7 = (TTree*)f7->Get("T");
  TTree* t8 = (TTree*)f8->Get("T");
  TTree* t9 = (TTree*)f9->Get("T");

  t1->Project("hNCherPhoton2MeV","ds.mc.nCherPhotons");
  t2->Project("hNCherPhoton3MeV","ds.mc.nCherPhotons");
  t3->Project("hNCherPhoton4MeV","ds.mc.nCherPhotons");
  t4->Project("hNCherPhoton5MeV","ds.mc.nCherPhotons");
  t5->Project("hNCherPhoton6MeV","ds.mc.nCherPhotons");
  t6->Project("hNCherPhoton7MeV","ds.mc.nCherPhotons");
  t7->Project("hNCherPhoton8MeV","ds.mc.nCherPhotons");
  t8->Project("hNCherPhoton9MeV","ds.mc.nCherPhotons");
  t9->Project("hNCherPhoton10MeV","ds.mc.nCherPhotons");


  TFile *ff = new TFile("saveNCherenPhoton.root","recreate");
  ff->cd();
  hNCherPhoton2MeV->Write();
  hNCherPhoton3MeV->Write();
  hNCherPhoton4MeV->Write();
  hNCherPhoton5MeV->Write();
  hNCherPhoton6MeV->Write();
  hNCherPhoton7MeV->Write();
  hNCherPhoton8MeV->Write();
  hNCherPhoton9MeV->Write();
  hNCherPhoton10MeV->Write();
  ff->Close();
 

/*
  int binTot = 240;
  int binE = 100;
  TH2F* h2d = new TH2F("h2d","",10, 0, 10, binTot, 0,6000);
  int k = 2;
  for(int j = 0; j<10; j++)
  {
    double ncher = 0;	  
    if (j == k)
    { 
      for(int i = 0; i<binTot;i++)
      {
      	if(k == 2) ncher = hNCherPhoton2MeV->GetBinContent(i);
      	if(k == 3) ncher = hNCherPhoton3MeV->GetBinContent(i);
      	if(k == 4) ncher = hNCherPhoton4MeV->GetBinContent(i);
      	if(k == 5) ncher = hNCherPhoton5MeV->GetBinContent(i);
      	if(k == 6) ncher = hNCherPhoton6MeV->GetBinContent(i);
      	if(k == 7) ncher = hNCherPhoton7MeV->GetBinContent(i);
      	if(k == 8) ncher = hNCherPhoton8MeV->GetBinContent(i);
      	if(k == 9) ncher = hNCherPhoton9MeV->GetBinContent(i);
      	if(k == 10) ncher = hNCherPhoton10MeV->GetBinContent(i);
        h2d->SetBinContent(j+1,i+1,ncher);
      }
    }
    //else { h2d->SetBinContent(j+1,i+1,ncher);}
    k = k+1;
  }
  h2d->Draw("colz");
//  hNCherPhoton2MeV->Fit("gaus");hNCherPhoton2MeV->Draw();
//  hNCherPhoton3MeV->Draw("sames");hNCherPhoton3MeV->Fit("gaus");
//  hNCherPhoton4MeV->Draw("sames");hNCherPhoton4MeV->Fit("gaus");
//  hNCherPhoton5MeV->Draw("sames");hNCherPhoton5MeV->Fit("gaus");
//  hNCherPhoton6MeV->Draw("sames");hNCherPhoton6MeV->Fit("gaus");
//  hNCherPhoton7MeV->Draw("sames");hNCherPhoton7MeV->Fit("gaus");
//  hNCherPhoton8MeV->Draw("sames");hNCherPhoton8MeV->Fit("gaus");
//  hNCherPhoton9MeV->Draw("sames");hNCherPhoton9MeV->Fit("gaus");
//  hNCherPhoton10MeV->Draw("sames");hNCherPhoton10MeV->Fit("gaus");
*/

}
