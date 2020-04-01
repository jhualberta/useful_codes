#include <iostream>
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <TMinuit.h>
#include <algorithm>
#include <iomanip>      // std::setprecision
//global data
TString filename = "SaveTest_JeffMC_107055.root";
TString fileSx = "fSpatial_1e5evts_Etot_allprocess.root";//fSpatial_1e5evts_107055.root";
//fSpatial_1e5evts_Etot.root
double fitRange = 6000;
double lowBin = -10000, hiBin = 10000;
const int iNum = 200;
int binTot = iNum;
//fit region for X mm
double hiCut = 2000; 
double loCut = -hiCut;
double Nfit[iNum];//data
double NfitX[iNum];//data
double NSx[iNum];//Spatial Resolution data
double Nerror[iNum];//errors in data
double Nx[iNum];//x values
double mu_P0;//initial param. mu_P for R(x) 
double scalepeak;

double posResol(double *x, double *par)
{   //par[0]=aE,par[1]=sigma_P, par[2]=mu_P, par[3]=tau_P, par[4]=scale

//Fit 3: for fixed aE, scale //fix a_e = 0.45 and fit for 3 params only, otherwise hard to converge
   double funcPosResol = (1.0-par[4])/(TMath::Sqrt(2*TMath::Pi())*par[0])*TMath::Exp(-0.5*pow((x[0]-par[1])/(par[0]),2))+0.55/(2*par[2])*TMath::Exp(-TMath::Abs(x[0]-par[1])/(par[2]));
   return funcPosResol*par[3];
}

double cal_convl(double *x, double *par)
{
  //Position Resolution
  double xx = x[0];
  //Fit 4 
  int funcpar = 5;
  TF1 *fun1 = new TF1("fun1",posResol,lowBin,hiBin,funcpar); 
  fun1->SetParLimits(0,100,400);
  fun1->SetParLimits(1,-500,500);
  fun1->SetParLimits(2,100,400);
  fun1->SetParLimits(4,0.3,0.8);

  fun1->SetParameters(par);

  //double beginX = (lowBin + lowBin+ (hiBin-lowBin)/binTot)/2;//begin center
  double convl = 0, xi=0, sxi=0, rxi=0; 
 
  for(int ii=0;ii<binTot;ii++) 
  {
     xi = Nx[ii];//beginX+ii*(hiBin-lowBin)/binTot;
     sxi = NSx[ii];
     rxi = fun1->Eval(xx-xi);
     convl = convl+ sxi*rxi;
     //cout<<convl<<" "<<sxi<<" "<<rxi<<endl;
  }
  return convl;
}

void Chi2easyfit()
{
  TFile *ff = new TFile(fileSx);
  TFile *ff2 = new TFile(filename);
  //original histo, Spatial Resolution
  TH1F* hSx_temp = (TH1F*)ff->Get("hSz");
  //original histo, Data
  TH1F* hNfit_temp = (TH1F*)ff2->Get("hfitZ");//copy directly	
  ///shift, scale and rebin
  int rebinX = ( (hiBin-lowBin)/binTot )/10; // origin 10 mm/bin to x mm/bin
  hSx_temp->Rebin(rebinX);hNfit_temp->Rebin(rebinX);
  TH1F* hSx = new TH1F("hSx","rebin hSx",binTot,lowBin,hiBin);
  TH1F* hNfit_temp2 = new TH1F("hNfit_temp2","rebin hfitZ",binTot,lowBin,hiBin);//rebin
  TH1F* hNfit = new TH1F("hNfit","rebin hfitZ",binTot,lowBin,hiBin);//shifted
  double mu_fit = hNfit_temp->GetMean();
  for(int i = 0;i<binTot;i++)
  {
    //cut hSx within [-2000,2000] mm, set others 0
    double xx = hSx_temp->GetXaxis()->GetBinCenter(i+1);
    if(xx<loCut+mu_fit || xx>hiCut+mu_fit) 
    {//2000 mm
      hSx->SetBinContent(i+1,0);
    } 
    else{
      hSx->SetBinContent(i+1,hSx_temp->GetBinContent(i));
    }
    //cut hNfit to lowBin, hiBin
    double xxfit = hNfit_temp2->GetXaxis()->GetBinCenter(i+1);
    if(xxfit<loCut+mu_fit|| xx>hiCut+mu_fit)
      hNfit_temp2->SetBinContent(i+1,0);
    else
    hNfit_temp2->SetBinContent(i,hNfit_temp->GetBinContent(i));
  }
  /// Find rough centers for hSx and hNfit 
  //int max1Bin = hSx->GetMaximumBin();
  hSx_temp->GetXaxis()->SetRangeUser(-20,40);
  int centerBinSx = hSx_temp->GetMinimumBin();//center "dip"
  //std::cout<<"center X Sx "<<hSx->GetBinCenter(centerBinSx)<<std::endl;
  hSx_temp->GetXaxis()->UnZoom();//resume hSx

  int meanBinFitX = hNfit_temp2->GetXaxis()->FindBin(hNfit_temp2->GetMean());  

  //hSx->Scale(1e4/hSx->Integral());
  //hNfit->Scale(1e4/hNfit->Integral());

  int shift = centerBinSx - meanBinFitX;

  //std::cout<<"shift "<<leftBin<<" "<<meanBin<<" "<<rightBin<<" "<<centerBin<<" "<<shift<<std::endl;
  int ishift;
  if(shift<0) ishift = -1*shift;
  else ishift = 0;
//  for (int i=ishift;i<=binTot;i++) 
  for(int i = 0;i<binTot;i++)
   hNfit->SetBinContent(i,hNfit_temp2->GetBinContent(i+0));
//  for(int i=0;i<ishift;i++)
//   hNfit->SetBinContent(0);
 
//  mu_P0 = hNfit->GetBinCenter(centerBin);//after shift, hNfit mean bin is at hSx center
  for(int iter=0;iter<binTot;iter++)
  {NSx[iter]=hSx->GetBinContent(iter+1);}

  for(int iter=0;iter<binTot;iter++)
  {Nx[iter]= hSx->GetXaxis()->GetBinCenter(iter+1);}

  mu_P0 = hNfit->GetBinCenter(meanBinFitX) - hSx->GetBinCenter(centerBinSx);

  std::cout<<"mu_P0 "<<mu_P0<<std::endl;
  double max1 = hSx->GetMaximum();//GetBinContent(centerBinSx);
  double max2 = hNfit->GetBinContent(meanBinFitX);
  hNfit->Scale(max1/max2);//scale to hSx

  //hSx->GetXaxis()->UnZoom();
  hNfit->GetXaxis()->SetRangeUser(loCut+mu_fit,hiCut+mu_fit);
  hNfit->Draw();
  const int Npar = 5; 
  TF1 *funMain = new TF1("fitConvl",cal_convl,-1*fitRange, fitRange,Npar); 
  funMain->FixParameter(4,0.45);
  funMain->SetParameters(180,mu_P0,320,20,0.45);
  funMain->SetParName(0,"#sigma_{P}");
  funMain->SetParName(1,"#mu_{P}");
  funMain->SetParName(2,"#tau_{P}");
  funMain->SetParName(3,"Scale");
  funMain->SetParName(4,"#alpha_{e}");
  funMain->SetLineColor(kRed);

  hNfit->Fit("fitConvl","rb");
  double aE = funMain->GetParameter(4),  aEerr= funMain->GetParError(4);
  double sigma = funMain->GetParameter(0), sigmaerr= funMain->GetParError(0);
  double tau = funMain->GetParameter(2), tauerr= funMain->GetParError(2);
  double mu = funMain->GetParameter(1), muerr= funMain->GetParError(1);
  double ndf = funMain->GetNDF();
  double chi2 = funMain->GetChisquare();
  std::cout<<setprecision(4);
  std::cout<<" & "<<aE<<"$\\pm$"<<aEerr<<" & "<<sigma<<"$\\pm$"<<sigmaerr<<" & "<<tau<<"$\\pm$"<<tauerr<<" & "<<mu<<"$\\pm$"<<muerr<<" & "<<chi2<<"/"<<ndf<<std::endl;
  TString newfile = "Convolved_"+ filename;
  TFile *fconv = new TFile(newfile,"recreate");
  TH1F *hConvl= new TH1F("hConvl","tempt convl",binTot,lowBin,hiBin);

 // for(int m=0; m<binTot; m++) {
 //   hConvl->SetBinContent(m+1,htempG->GetBinContent(m));
 // }
 // hConvl->Scale(scalepeak/hConvl->GetMaximum());

  fconv->cd();
  hSx->Write();hNfit->Write();
  funMain->Write();
//  hConvlOutput->Write();


}
