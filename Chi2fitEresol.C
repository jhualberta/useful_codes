#include <iostream>
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <TMinuit.h>
#include <algorithm>
////global data

const int binTot = 80;
double hiBin = 12;
double lowBin = 0;
double binStep = (hiBin-lowBin)/binTot;
// parameters
double sigma;
double deltaE;
double Nscale;
double scalepeak;
double range1 = 3;
double range2 = 8.8;//10;
TFile *fdata = new TFile("../Merged_RatWater_AntiNeutrinoMC_WaterN16sourceRun_r107055_s0.root");//
//Merged_RatWater_Calibration_r0000107055_p009.root");//

double Ns[binTot];
double NsX[binTot];
double Nerror[binTot];//errors in data
double nPhotonTable[]={0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 25.0, 75.0, 125.0, 125.0, 150.0, 200.0, 325.0, 300.0, 275.0, 425.0, 400.0, 450.0, 375.0, 550.0, 450.0, 650.0, 750.0, 775.0, 750.0, 825.0, 850.0, 650.0, 950.0, 750.0, 1075.0, 700.0, 1200.0, 850.0, 1150.0, 1025.0, 1425.0, 1350.0, 1200.0, 1425.0, 1525.0, 1325.0, 1600.0, 1700.0, 1500.0, 1700.0, 1800.0, 1575.0, 1650.0, 1900.0, 1450.0, 2050.0, 2150.0, 2000.0, 1900.0, 2200.0, 1675.0, 2200.0, 1975.150390625, 2500.0, 2475.322021484375, 1950.0, 1824.919189453125, 1775.0548095703125, 2424.69189453125, 2524.664306640625, 2624.63916015625, 2824.593994140625, 2350.0, 2324.7216796875, 2550.0, 2250.0, 2475.322021484375, 3100.0, 2675.372802734375, 3050.0, 2775.395263671875, 2975.43603515625, 3050.0, 3150.0, 3200.0, 3050.0, 3200.0, 3550.109619140625, 3275.487548828125, 3075.454345703125, 3849.7431640625, 3075.454345703125, 3150.0, 3624.462890625, 3775.5556640625, 2924.573486328125, 3575.53076171875, 3575.53076171875, 3424.489990234375, 3200.0, 3824.438720703125, 3575.53076171875, 3924.427490234375, 4124.40673828125};

double eResol(double *x, double *par)
{  // x[0] = Teff,  x[1] = E 
   // par[0] = b, par[1] = deltaE, par[2] = N
   double fP = 1./(sqrt(2*TMath::Pi()*x[1])*par[0])*exp(-((1+par[1])*x[0]-x[1])*((1+par[1])*x[0]-x[1])/(2*par[0]*par[0]*x[1]));
   // cout<<x[0]<<" "<<x[1]<<" "<<funcPosResol<<endl;
   return fP*par[2];
}


double cal_convl(double *x, double *par) // x is Teff
{
  double xe = 0;
  double Teff = x[0];
  double convl=0;
  double sxi=0, rxi=0, xi=0;
  double Te = 0.1, Ps = 0;
  int len = binTot;
  for(int ii=0;ii<len;ii++)
  {
    sxi = Ns[ii];
    xi = xi + binStep;
    //Ps = hSx->GetBinContent(ind);
    double xxx[2] = {Teff,Te};// !!! Te>0
    double Pe = sxi*eResol(xxx,par);
    convl += Pe;
    Te += binStep;
    //cout<<Teff<<" "<<Te<<" "<<convl<<endl;
  }
  return convl;
}

void Chi2fitEresol()
{
  //TH1F *hE_Compton = (TH1F*)fmc->Get("hE_Compton");
  //TH1F *hE_PairProduct = (TH1F*)fmc->Get("hE_PairProduct");
  //TH1F *hE_PEabsorb = (TH1F*)fmc->Get("hE_PEabsorb");
  //TH1F *hEsum= new TH1F("hEsum","summed gamma energy",160,0,4000);
  //TH1F *hEsum = (TH1F*)fmc->Get("hNCherenkov");//hEtot");//hNCherenkov");
  //hE_Compton->Add(hE_PairProduct);
  //hE_Compton->Add(hE_PEabsorb);
  //hEsum->Add(hE_Compton);
  //hEsum->Draw();
  TFile *fmc = new TFile("Merged_ftrack.root");//NewSaveNCheren.root");//5MeV.root");//Merged_fEnergy_Rat6176_test_6and7MeVgamma.root");//Merged_fSpatial_New6176_test_6and7MeVgamma_1e5evts.root");
  TH1F *hNCherenkov = (TH1F*)fmc->Get("hNCherenkovWaterNhit6");
  TH1F* hSx = new TH1F("hSx","rebin hSx",binTot,lowBin,hiBin);//MC distribution
  double p0 = 231.4;//-225.118;
  double p1 = 398.9;// 408.66;

  for (int i=0;i<4000;i++)
  {
      double ycounts = hNCherenkov->GetBinContent(i);
      double nphoton = hNCherenkov->GetXaxis()->GetBinCenter(i);
      /// linear transformation method
      double xnew = (1./p1*nphoton +p0/p1);
      hSx->Fill(xnew,ycounts);
  }

  TTree *Tdata = (TTree*)fdata->Get("T");
  TH1D* hNfit_temp = new TH1D("hNfit_temp","rebin hfitX",binTot,lowBin,hiBin);
  TH1D* hNfit= new TH1D("hNfit","rebin hfitX",binTot,lowBin,hiBin);//fitted energy
  double energy, posRad, itr;
  UInt_t nhits;
  Tdata->SetBranchAddress("energy", &energy);
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("posRad", &posRad);
  Tdata->SetBranchAddress("itr", &itr);

  for(int i =0;i<Tdata->GetEntries();i++)
  {
    Tdata->GetEntry(i);
    //cout<<energy<<" "<<nhits<<endl;
    if(nhits>5 && posRad<5300 && itr>0.55)
    {
      hNfit->Fill(energy);
    }
  }
  hSx->Scale(1./hSx->Integral());//hNfit->GetMaximum()/hSx->GetMaximum());
  //hNfit->Draw();
  /// Fill hSx by Gaussian
  //
  //TRandom3 rndgen;
  //for (int i = 0; i < hEsum->GetEntries(); ++i) hSx->Fill(rndgen.Gaus(5.364,0.8522));//084,0.7782));//rndgen.Gaus(5.117,1.026));

  /// artificially invert hSx
  //for (int i=0;i<binTot;i++)
  //{
  //   hSx->SetBinContent(i+1,hSx->GetBinContent(binTot-i));
  //}
   
   // hSx->Draw();
  //for(int i = 0;i<binTot;i++)
  //{
  //  //double xx = hNfit->GetXaxis()->GetBinCenter(i+1);
  //  //if(xx<range1 || xx>range2)// [3, 10] MeV
  //  // {
  //  //  hSx->SetBinContent(i+1,0);
  //  //  hNfit->SetBinContent(i+1,0);
  //  // }
  //  // else
  //  // {
  //    hNfit->SetBinContent(i+1,hNfit_temp->GetBinContent(i));
  //  // }
  //}
  //scalepeak = hNfit->GetMaximum()/hSx->GetMaximum();
  //hNfit->Scale(1./hNfit->Integral());
  //hSx->Scale(scalepeak);
  //hNfit->Scale(10000./hNfit->Integral());
  ////input data	
  //for(int iter=0;iter<binTot;iter++)
  //{Nfit[iter]=hNfit->GetBinContent(iter);}	
  //// x axis values
  //for(int iter=0;iter<binTot;iter++)
  //{NfitX[iter]=hNfit->GetXaxis()->GetBinCenter(iter);}
  int len = binTot;
  for(int iter=0;iter<len;iter++)
  {
    if(hNfit->GetBinError(iter) == 0) Nerror[iter+1]=0.1;
    else Nerror[iter]=hNfit->GetBinError(iter+1);
  }
  
  for(int iter=0;iter<len;iter++)
  {Ns[iter]=hSx->GetBinContent(iter+1);}
  
  for(int iter=0;iter<len;iter++)
  {NsX[iter]= hSx->GetXaxis()->GetBinCenter(iter+1);}

  TF1 *funMain = new TF1("funMain",cal_convl,range1, range2,3);

  //int rebin = 2;
  //hNfit->Rebin(rebin);hSx->Rebin(rebin);

  double bb = 0.1;
  double deltaE = -0.1;
  double Nscale = hNfit->GetMaximum();
//  TF1 *f1= new TF1("f1",eResol,range1, range2,3);
//  f1->SetParameters(bb,deltaE,Nscale); 

  //funMain->FixParameter(0,bb);
  //funMain->FixParameter(1,deltaE);
  //funMain->FixParameter(2,Nscale);
  TCanvas *c = new TCanvas("c","",800,600);
  c->cd();
  cout<<"set pars "<<Nscale<<endl;
  double peakfunc = funMain->GetMaximum();
  hNfit->Sumw2();
  hNfit->Scale(1./Nscale);
  //Nscale = 10000;

  //hNfit->Scale(1./hNfit->Integral());
  //hNfit->Scale(Nscale/hNfit->Integral());
  funMain->SetParNames("b","#delta_{E}","N");
  funMain->SetParameters(bb,deltaE,Nscale);//SetParameters(bb,deltaE,Nscale);
  funMain->SetLineColor(kRed);
  hNfit->Fit("funMain","MErb");
  hNfit->SetLineWidth(2);
  hSx->Scale(hNfit->GetMaximum()/hSx->GetMaximum());
  //hNfit->Draw();
  hSx->Draw("same");

}
