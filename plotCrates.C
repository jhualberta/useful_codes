#include <string> 

#include <vector> 

void plotCrates()
{
 using namespace std;
 TFile *ff = new TFile("ResolPartial_Analysis50_r0000250127_s000_p000.root");

 vector<TH2F*> hCrateCardChannel;
 for(int i = 0;i<19;i++)
 {
   TH2F *htemp = new TH2F("htemp","card, channel", 16, 0, 15, 32,0,31);
   hCrateCardChannel.push_back((TH2F*)htemp->Clone(Form("hcrate%u",i)));
   delete htemp;
 }

 TCanvas *c = new TCanvas("c","",800,600);
 for(int i = 0;i<19;i++)
 {
   TString num; num.Form("%d",i);
   TString hh= "hcrate"+num;
   const char *hname = hh;
   hCrateCardChannel[i] = (TH2F*)ff->Get(Form("hcrate%u",i));
 }

  c->Divide(10,2);
  gStyle->SetOptStat(0);
  for(int i = 0;i<19;i++)
  {
    c->cd(i+1);
    hCrateCardChannel[i]->Draw("colz");
  }

}
