#include "TH1.h"
#include "TF1.h"
#include "TKDE.h"
#include "TCanvas.h"
/*#include "TStopwatch.h"*/
#include "TRandom.h"
#include "Math/DistFunc.h"
#include "TLegend.h"
#include "TFile.h"
// test TKDE

double breit(double *x, double *par)
{
  double val = par[0]*ROOT::Math::breitwigner_pdf(x[0],par[1],par[2])+par[3];
  return val;
}

Double_t lorentz(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
   + .25*par[1]*par[1]);
}

double gaus1(double *x, double *par)
{
  double gaus = par[0]*TMath::Gaus(x[0],0.665,par[1],true);
  return gaus;
}

double total(double *x, double *par)
{
  double total = par[0]*ROOT::Math::breitwigner_pdf(x[0],par[1],par[2])+par[3]+par[4]*TMath::Gaus(x[0],0.665,par[5],true);
  return total;
}

void fitBreitWignerMore(int n = 200) {
   // generate some gaussian points
   int nbin = 200;
   double xmin = -1;
   double xmax = 1;
   // generate some points with bi- gaussian distribution

   TString s = "Merged_10MeV_Resol.root";
   // read the input file
   TFile* fileName = new TFile(s);

   TH1F *h1 = (TH1F*)fileName->Get("hcosThetaPMT_cut");
   //h1->Rebin(2);
   std::vector<double> data(n);
   for (int i = 0; i < n; i++) {
     data[i] = h1->GetBinContent(i+1);
   }
   // scale histogram
   h1->Scale(1./h1->Integral(),"width" );
   h1->SetTitle("Breit-Wigner");
   //h1->Draw();
   // drawn true normalized density

   Double_t par[6];
  
   TF1 *f1 = new TF1("f1",breit,xmin,xmax,4);
   TF1 *g2 = new TF1("g2",lorentz,0.6,0.7,2);
   TF1 *total = new TF1("total",total,-1,1,6);

   f1->SetParameter(0,0.8);
   f1->SetParameter(1,0.3);
   f1->SetParameter(2,0.665);
   //f1->FixParameter(3,0.2);
   h1->Fit(f1);

   //g2->SetParameters(h1->GetMaximum(),0.3);
   g2->SetParameter(0,h1->GetMaximum());
   h1->Fit(g2,"R+");
   f1->GetParameters(&par[0]);
   g2->GetParameters(&par[4]);
   total->SetParameters(par);
   total->FixParameter(1,0.665);
//   total->FixParameter(5,0.665);
   h1->Fit(total);

}
