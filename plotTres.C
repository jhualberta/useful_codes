#include "TH1F.h"
#include <vector> 

void plotTres()
{  
  ofstream out,outputbadfiles;
  ostringstream oss;
  ifstream in0, in1;
  in0.open("mpsimpleList.dat");
  in1.open("ratList.dat");
  char filenames[1500];
  char filenames2[1500];
  vector<TH1F*> htres_rat;
  vector<TH1F*> htres_mpw;

  vector<TH1F*> hfitT_rat;
  vector<TH1F*> hfitT_mpw;

  vector<TH1F*> hfitZ_rat;
  vector<TH1F*> hfitZ_mpw;

  while(in1>>filenames){
   //cout<<" filenames "<<filenames<<endl;
   TFile *f = new TFile(filenames); 
   TH1F *htres = (TH1F*)f->Get("htRes_trig");
   TH1F *hfitZ = (TH1F*)f->Get("hfitZ_trig");
   TH1F *hfitT = (TH1F*)f->Get("hFitTime_trig");
   htres_rat.push_back(htres);
   hfitZ_rat.push_back(hfitZ);
   hfitT_rat.push_back(hfitT);
   //delete f,htres,hfitZ,hfitT;
  }

  while(in0>>filenames2){
   //cout<<" filenames "<<filenames2<<endl;
   TFile *f1 = new TFile(filenames2);
   TH1F *htres1 = (TH1F*)f1->Get("htRes_trig");
   TH1F *hfitZ1 = (TH1F*)f1->Get("hfitZ_trig");
   TH1F *hfitT1 = (TH1F*)f1->Get("hFitTime_trig");
   htres_mpw.push_back(htres1);
   hfitZ_mpw.push_back(hfitZ1);
   hfitT_mpw.push_back(hfitT1);
   //delete f,htres,hfitZ,hfitT;
  }

  // cout<<htres_rat.size()<<" "<<htres_mpw.size()<<endl;

  TCanvas *c = new TCanvas("c","",800,600);
  c->Divide(5,5);
  for(size_t i = 1;i<htres_rat.size()+1;i++)
  {
   c->cd(i);
   gPad->SetLogy();
   htres_rat[i-1]->SetLineColor(kBlue);
   htres_mpw[i-1]->SetLineColor(kRed);
   htres_mpw[i-1]->Draw();htres_rat[i-1]->Draw("sames");
  }

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(5,5);
  for(size_t i = 1;i<hfitT_rat.size()+1;i++)
  {
   c1->cd(i);
   gPad->SetLogy();
   hfitT_rat[i-1]->SetLineColor(kBlue);
   hfitT_mpw[i-1]->SetLineColor(kRed);
   hfitT_rat[i-1]->Draw();hfitT_mpw[i-1]->Draw("sames");
  }


  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Divide(5,5);
  for(size_t i = 1;i<hfitZ_rat.size()+1;i++)
  {
   c2->cd(i);
   gPad->SetLogy();
   hfitZ_rat[i-1]->SetLineColor(kBlue);
   hfitZ_mpw[i-1]->SetLineColor(kRed);
   hfitZ_mpw[i-1]->Draw();hfitZ_rat[i-1]->Draw("sames");
 }

  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->cd();htres_rat[0]->Draw("PLC PMC");
  for(size_t i = 2;i<htres_rat.size()+1;i++)
  {
   gPad->SetLogy();
   htres_rat[i-1]->Draw("SAME PLC PMC");
  }


  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->cd();htres_mpw[0]->Draw("PLC PMC");
  for(size_t i = 2;i<htres_mpw.size()+1;i++)
  {
   gPad->SetLogy();
   htres_mpw[i-1]->Draw("SAME PLC PMC");
  }

}
