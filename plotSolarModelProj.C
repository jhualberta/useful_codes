#include "TH2F.h"

void plotSolarModelProj()
{

  TFile *ff = new TFile("SunData.root");
  TH2F *flux1 = (TH2F*)ff->Get("bps09gs");
  int biny = flux1->GetYaxis()->GetNbins();
  int binx = flux1->GetXaxis()->GetNbins();
  TString s0 = "N_{e}";
  TString s1 = "pp";
  TString s2 = "pep";
  TString s3 = "hep";
  TString s4 = "7Be";
  TString s5 = "8B";
  TString s6 = "13N";
  TString s7 = "15O";
  TString s8 = "17F";
  TH1F *hflux0 = new TH1F("hflux0",s0,binx,0,0.5);
  TH1F *hflux1 = new TH1F("hflux1",s1,binx,0,0.5);
  TH1F *hflux2 = new TH1F("hflux2",s2,binx,0,0.5);
  TH1F *hflux3 = new TH1F("hflux3",s3,binx,0,0.5);
  TH1F *hflux4 = new TH1F("hflux4",s4,binx,0,0.5);
  TH1F *hflux5 = new TH1F("hflux5",s5,binx,0,0.5);
  TH1F *hflux6 = new TH1F("hflux6",s6,binx,0,0.5);
  TH1F *hflux7 = new TH1F("hflux7",s7,binx,0,0.5);
  TH1F *hflux8 = new TH1F("hflux8",s8,binx,0,0.5);
  hflux0->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux1->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux2->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux3->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux4->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux5->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux6->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux7->GetXaxis()->SetTitle("R/R_{#odot}");
  hflux8->GetXaxis()->SetTitle("R/R_{#odot}");

  cout<<biny<<" "<<binx<<endl;
  for(int i = 0;i<1000;i++)
  {
    hflux0->SetBinContent(i+1, flux1->GetBinContent(i+1,1));
    hflux1->SetBinContent(i+1, flux1->GetBinContent(i+1,2));
    hflux2->SetBinContent(i+1, flux1->GetBinContent(i+1,3));
    hflux3->SetBinContent(i+1, flux1->GetBinContent(i+1,4));
    hflux4->SetBinContent(i+1, flux1->GetBinContent(i+1,5));
    hflux5->SetBinContent(i+1, flux1->GetBinContent(i+1,6));
    hflux6->SetBinContent(i+1, flux1->GetBinContent(i+1,7));
    hflux7->SetBinContent(i+1, flux1->GetBinContent(i+1,8));
    hflux8->SetBinContent(i+1, flux1->GetBinContent(i+1,9));
  }

  TFile *ffnew = new TFile("projectFlux.root","recreate");
  ffnew->cd();
  hflux0->Write();
  hflux1->Write();
  hflux2->Write();
  hflux3->Write();
  hflux4->Write();
  hflux5->Write();
  hflux6->Write();
  hflux7->Write();
  hflux8->Write();

}
