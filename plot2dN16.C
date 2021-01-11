{

  TFile *fmcMPW = new TFile("MPWn16_Merged_FitMPtest_nhit6_AntiNeutrinoMC_WaterN16sourceRun_r107055_s0_p8.root");
  TFile *fdataMPW = new TFile("Merged_MPWn16_WaterMP6716_Calibration_r0000100934.root");

  TFile *fmcRat = new TFile("Merged_RatWater_AntiNeutrinoMC_WaterN16sourceRun_r107055.root");
  TFile *fdataRat = new TFile("Merged_RatWater_Calibration_r0000107055.root");


//  TH1F *hEmcMPW = new TH1F("hEmcMPW","",100,0,100);
//  TH1F *hEdataMPW = new TH1F("hEdataMPW","",100,0,100);
//  TH1F *hEmcRat = new TH1F("hEmcRat","",100,0,100);
//  TH1F *hEdataRat = new TH1F("hEdataRat","",100,0,100);
  int bin1 = 100;
  TH2F *hmcMPWbeta14vsE = new TH2F("hmcMPWbeta14vsE","",400,-2,4,bin1,0,20);
  TH2F *hmcRatbeta14vsE = new TH2F("hmcRatbeta14vsE","",400,-2,4,bin1,0,20);


  TH2F *hmcMPWnhitvsE = new TH2F("hmcMPWnhitvsE","",bin1,0,20,100,0,100,);
  TH2F *hmcRatnhitvsE = new TH2F("hmcRatnhitvsE","",bin1,0,20,100,0,100,);
  hmcRatN16vsE->GetXaxis()->SetTitle("#beta_{14}");
  hmcRatN16vsE->GetYaxis()->SetTitle("E_{fit} [MeV]");

  hmcMPWnhitvsE->GetXaxis()->SetTitle("E_{fit} [MeV]");
  hmcMPWnhitvsE->GetYaxis()->SetTitle("nhits");

  hmcRatnhitvsE->GetXaxis()->SetTitle("E_{fit} [MeV]");
  hmcRatnhitvsE->GetYaxis()->SetTitle("nhits");

  TTree *TmcMPW = (TTree*)fmcMPW->Get("T"); 
  TTree *TdataMPW = (TTree*)fdataMPW->Get("T");
  TTree *TmcRat = (TTree*)fmcRat->Get("T");
  TTree *TdataRat = (TTree*)fdataRat->Get("T");

  TmcMPW->Project("hmcMPWbeta14vsE","energy:beta14");
  hmcMPWbeta14vsE->GetXaxis()->SetTitle("#beta_{14}");
  hmcMPWbeta14vsE->GetYaxis()->SetTitle("E_{fit} [MeV]");

  //  TdataMPW->Project("hEdataMPW","nhits","");

  TmcRat->Project("hmcRatbeta14vsE","energy:beta14");
  hmcRatbeta14vsE->GetXaxis()->SetTitle("#beta_{14}");
  hmcRatbeta14vsE->GetYaxis()->SetTitle("E_{fit} [MeV]");


  TmcMPW->Project("hmcRatnhitvsE","nhits:energy");
  hmcMPWnhitvsE->GetYaxis()->SetTitle("Nhits");
  hmcMPWnhitvsE->GetXaxis()->SetTitle("E_{fit} [MeV]");

  //  TdataMPW->Project("hEdataMPW","nhits","");

  TmcRat->Project("hmcRatnhitvsE","nhits:energy");
  hmcRatnhitvsE->GetYaxis()->SetTitle("Nhits");
  hmcRatnhitvsE->GetXaxis()->SetTitle("E_{fit} [MeV]");

 //  TdataRat->Project("hEdataRat","nhits","");
//  TmcMPW->Project("hEmcMPW","energy","");
//  TdataMPW->Project("hEdataMPW","energy","");
//  TmcRat->Project("hEmcRat","energy","");
//  TdataRat->Project("hEdataRat","energy","");

//  hEmcMPW->SetLineColor(kRed);
//  hEdataMPW->SetLineColor(kRed);
//
//  hEmcMPW->SetLineStyle(2);
//  hEmcRat->SetLineStyle(2);
//
// // hEdataMPW->Draw();
// // hEmcMPW->Draw("sames");
// // hEmcRat->Draw("sames");
// // hEdataRat->Draw("sames");
//hEmcMPW->GetXaxis()->SetTitle("E_{fit} [MeV]");
// hEmcMPW->Draw();
// hEmcRat->Draw("sames");

 TCanvas *c = new TCanvas("c","",800,600);
//  hmcMPWbeta14vsE->Draw("colz");
//  hmcRatbeta14vsE->Draw("colz");

  hmcMPWnhitvsE->Draw("colz");
//  hmcRatnhitvsE->Draw("colz");

}
