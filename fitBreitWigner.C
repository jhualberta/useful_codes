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

Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
   + .25*par[1]*par[1]);
}

double breit(double *x, double *par)
{
  double val = par[0]*ROOT::Math::breitwigner_pdf(x[0],par[1],par[2])+par[3];
  //double val = ROOT::Math::breitwigner_pdf(x[0],par[0],par[1])+par[2];
  return val;
}


double gaus1(double *x, double *par)
{
  double gaus = par[2]*TMath::Gaus(x[0],par[0],par[1],true);
  return gaus;
}


void fitBreitWigner(int n = 200) {
   // generate some gaussian points
   int nbin = 200;
   double xmin = 0;
   double xmax = 1;
   // generate some points with bi- gaussian distribution

   TString s = "results_simple.root";
   // read the input file
   TFile* fileName = new TFile(s);

   TH1F *h1 = (TH1F*)fileName->Get("hcosTheta_prompt");
   h1->Sumw2();
   h1->Rebin(2);
   h1->Scale(1./h1->Integral());
   std::vector<double> data(n);
   for (int i = 0; i < n; i++) {
     data[i] = h1->GetBinContent(i+1);
   }
   // scale histogram
   h1->Scale(1./h1->Integral(),"width" );

   h1->SetStats(true);

   h1->SetTitle("Breit-Wigner");
   h1->Draw();
   // drawn true normalized density

   Double_t par[6];
  
//   TF1 *f1 = new TF1("f1","[0]*ROOT::Math::breitwigner_pdf(x,[1],[2])+[3]",xmin,xmax);
//   // TF1 *f1 = new TF1("f1","[0]*ROOT::Math::cauchy_pdf(x,[1],[2])+[3]",xmin,xmax);
//   f1->SetParameter(0,1);
//   f1->SetParameter(1,0.6);
//   f1->SetParameter(2,0.2);
//   f1->SetParameter(3,0);
// 
//   f1->SetParName(0,"c0");
//   f1->SetParName(1,"#Gamma");
//   f1->SetParName(2,"m");
//   f1->SetParName(3,"c1");
//   h1->Fit(f1,"R");
//   double chi2 = f1->GetChisquare();
//   double ndf = f1->GetNDF()/f1->GetMaximum();
//   cout<<"chi/ndf "<<chi2<<"/"<<ndf<<"="<<chi2/ndf<<endl;

   TF1 *f2 = new TF1("f2","[0]*ROOT::Math::Gaus(x,[1],[2])+[3]",xmin,xmax);
   f2->SetParameter(0,1);
   f2->SetParameter(1,0.6);
   f2->SetParameter(2,0.2);
   f2->SetParameter(3,0);

   f2->SetParName(0,"c0");
   f2->SetParName(1,"#mu");
   f2->SetParName(2,"#sigma");
   f2->SetParName(3,"c1");

   h1->Fit(f2,"R");
   double chi2 = f2->GetChisquare();
   double ndf = f2->GetNDF()/f2->GetMaximum();
   cout<<"chi/ndf "<<chi2<<"/"<<ndf<<"="<<chi2/ndf<<endl;
}
