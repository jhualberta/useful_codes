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

void toyBreitWigner() {
   // generate some gaussian points
   int nbin = 200;
   double xmin = -1;
   double xmax = 1;
   // generate some points with bi- gaussian distribution
   TH1F *h1 = new TH1F("h1","",200,-1,1); 
   TF1 *f1 = new TF1("f1","[0]*ROOT::Math::breitwigner_pdf(x[0],[1],[2])+[3]+[4]*TMath::Gaus(x[0],0.665,[5])",xmin,xmax);
   f1->SetParameters(7.22365e-01,6.65000e-01,7.12693e-01,2.14739e-01,7.83947e-02,8.01957e-02);
   f1->Draw();
   Double_t u1,u2,theta;//random variables
   TRandom1 rand;
   Double_t f;
   vector<Double_t> x;
   int run = 0;
   const int RUN = 10000000;
   TRandom1 rand1, rand2;
   while(run<RUN+1)
   {
   //TRandom3 rand1, rand2;
   double u1 = rand1.Uniform(-1,1);
   double u2 = rand2.Uniform(0,1.23);
   u2=double(rand.Rndm());
   double u3 = f1->Eval(u1);
   //cout<<u1<<" "<<u2<<" "<<u3<<endl;
   if (u3>u2) h1->Fill(u1);
   run++;
   }
   TFile *ff = new TFile("fDump.root","recreate");
   ff->cd();h1->Write();
   ff->Close();
}
