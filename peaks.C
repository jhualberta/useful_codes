// Illustrates how to find peaks in histograms.
// This script generates a random number of gaussian peaks
// on top of a linear background.
// The position of the peaks is found via TSpectrum and injected
// as initial values of parameters to make a global fit.
// The background is computed and drawn on top of the original histogram.
//
// To execute this example, do
//  root > .x peaks.C  (generate 10 peaks by default)
//  root > .x peaks.C++ (use the compiler)
//  root > .x peaks.C++(30) (generates 30 peaks)
//
// To execute only the first part of the script (without fitting)
// specify a negative value for the number of peaks, eg
//  root > .x peaks.C(-20)
//
//Author: Rene Brun
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

Int_t npeaks = 30;
Double_t fpeaks(Double_t *x, Double_t *par) {
      Double_t result;
      Double_t norm  = par[0];
      Double_t mean  = par[1];
      Double_t sigma = par[2];
      result = norm*TMath::Gaus(x[0],mean,sigma);
   return result;
}

void peaks(){
    TFile *ff = new TFile("Merged_Fit_100934.root","read");
	
    TH2F *htResVsPMTposx = (TH2F*)ff->Get("htResVsPMTposx");

    TH2F *htResVsPMTposx_left = (TH2F*)htResVsPMTposx->Clone("htResVsPMTposx_left");
    TH2F *htResVsPMTposx_right = (TH2F*)htResVsPMTposx->Clone("htResVsPMTposx_right");
	
    htResVsPMTposx_left->GetXaxis()->SetRangeUser(-8000,-7000);
    htResVsPMTposx_right->GetXaxis()->SetRangeUser(7000,8000);
	
    TH1F *htRes_left = (TH1F*)htResVsPMTposx_left->ProjectionY(); 
    TH1F *htRes_right =  (TH1F*)htResVsPMTposx_right->ProjectionY(); 
    TH1F *h2 = (TH1F*)htRes_left->Clone("h2");

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900); 
    c1->Divide(1,2);c1->cd(1);gPad->SetLogy();
    TSpectrum *s = new TSpectrum(3);
    Int_t nfound = s->Search(htRes_left,10,"",100);//Search (const TH1 *hist, Double_t sigma=2, Option_t *option="", Double_t threshold=0.05)
    printf("Found %d candidate peaks to fit\n",nfound);
    //Estimate background using TSpectrum::Background
    TH1 *hb = s->Background(htRes_left,20,"same");

    TH1F* htRes_left_sub = (TH1F*)htRes_left->Clone("htResVsPMTposx_left_sub");  
    htRes_left_sub->GetXaxis()->SetRangeUser(65,85);
    TH1F *h22 = (TH1F*)htRes_left_sub->Clone("h22");
    TSpectrum *s1 = new TSpectrum(20);
    Int_t nfound1 = s1->Search(htRes_left_sub,10,"",0.001);//Search (const TH1 *hist, Double_t sigma=2, Option_t *option="", Double_t threshold=0.05)
    printf("Found %d candidate peaks to fit\n",nfound1);
    ///Estimate background using TSpectrum::Background
    TH1 *hb1 = s1->Background(htRes_left_sub,20,"same");

    //htRes_left->Draw();htRes_left_sub->Draw("same");
    htRes_left->Draw();//htRes_left_sub->Draw("same");
     c1->Update();
 
   ////estimate linear background using a fitting method
   c1->cd(2);
   gPad->SetLogy(); 
   TF1 *fline = new TF1("fline","pol1",0,1000);
   htRes_left->Fit("fline","qn");
   //Loop on all found peaks. Eliminate peaks at the background level
   Double_t par[3000];
   par[0] = fline->GetParameter(0);
   par[1] = fline->GetParameter(1);
   npeaks = 0;
   Float_t *xpeaks = s->GetPositionX();
   for (int p=0;p<nfound;p++) {
      Float_t xp = xpeaks[p];
      Int_t bin = htRes_left->GetXaxis()->FindBin(xp);
      Float_t yp = htRes_left->GetBinContent(bin);
      if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
      par[0] = yp;
      par[1] = xp;
      par[2] = 3;
      npeaks++;
   }
   printf("Found %d useful peaks to fit\n",npeaks);
   printf("Now fitting: Be patient\n");
   TF1 *fit = new TF1("fit",fpeaks,-100,100,3);
   ////we may have more than the default 25 parameters
   TVirtualFitter::Fitter(h2,2);
   fit->SetParameters(par);
   fit->SetNpx(1000);
   h2->Fit("fit");  
	//htRes_left->Draw();htResVsPMTposx_left_sub->Draw("same");//htRes_right->Draw("sames");
}
