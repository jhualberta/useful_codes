#include "TH1F.h"
#include <unistd.h>
#include<vector> 
// using namespace std::this_thread; // sleep_for, sleep_until
// using namespace std::chrono; // nanoseconds, system_clock, seconds

void anim_tResWindow() {
//
// This script is a slightly modified version of hsum.C.
// When run in batch mode, it produces an animated gif file.
//Authors: Rene Brun, Valeriy Onuchin
  TFile *f = new TFile("ResolMPnew_FitMPScint_fitDirect_2p5MeVbeta_center_+x.root");
  TH1F *htResMC = (TH1F*)f->Get("htResMC");
  const int plotNum = 20;
  c1 = new TCanvas("c1","The HSUM example",200,10,800,400);
  c1->SetGrid();
  c1->Divide(2,1);
  gBenchmark->Start("hsum");

// Create some histograms.
  vector<TH1F*> hcosTheta_cutwin1;
  for(int i =0;i<plotNum;i++) 
  { 
    TH1F *htemp = new TH1F("hcosTheta_cutwin1","cut window 1",200,-1,1);
    hcosTheta_cutwin1.push_back(htemp);
    delete htemp;
  }

  TSlider *slider = 0;
  gSystem->Unlink("htResWindow.gif"); // delete old file
 
// Fill histograms randomly
  const Int_t kUPDATE = 1;
  Int_t gifcnt = 0;
  for ( Int_t i=0; i<plotNum; i++) {
     hcosTheta_cutwin1[i] = (TH1F*)f->Get(Form("hcosTheta_cutwin1_%u",i));
     if ((i%kUPDATE) == 0) {
        //if (i == kUPDATE) 
        {
           c1->cd(1);
           hcosTheta_cutwin1[i]->Draw();
           c1->cd(2);gPad->SetLogy();
           htResMC->GetXaxis()->SetRangeUser(-10,10);
           htResMC->Draw(); 
           TLine *line = new TLine(-6,0,-6,1e5);//c1->GetUymax());
           TLine *line1 = new TLine(-1+i*0.5, 0,-1+i*0.5,1e5);
           line->SetLineColor(kRed);line1->SetLineColor(kRed);
           line->Draw("same");line1->Draw("same");
           c1->Update();
        }
        // if (slider) slider->SetRange(0,Float_t(i)/20.);
        c1->Modified();
        c1->Update();
        if (gROOT->IsBatch()) {
           c1->Print("htResWindow.gif+100");// default gif refresh time = 10ms, now is 50*10 ms = 0.5 s
           //printf("i = %d\n", i);
        }
     }
  }
  //slider->SetRange(0,1);
  hcosTheta_cutwin1[0]->Draw("sameaxis"); // to redraw axis hidden by the fill area
  c1->Modified();
  // make infinite animation by adding "++" to the file name
  if (gROOT->IsBatch()) c1->Print("htResWindow.gif++");
  
  //You can view the animated file hsumanim.gif with Netscape/IE or mozilla
  
  gBenchmark->Show("hsum");
}
