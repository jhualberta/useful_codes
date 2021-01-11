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
#include <iomanip>      // std::setprecision
//global data
TString filename = "Merged_MPW_nhit6_Calibration_r106930_s0.root";
bool defaultFit = false;

double rangeDefault = 2000;
double fitRange = rangeDefault;
double lowBin = -1, hiBin = 5;// the bin width must be 10 mm/bin to match hSx !!!
const int iNum = 2000;
int binTot = 200;
//fit region for X mm
double hiCut = 3000;
double loCut = -hiCut;
double Nfit[iNum];//data
double NfitX[iNum];//data
double NSx[iNum];//Spatial Resolution data
double Nerror[iNum];//errors in data
double Nx[iNum];//x values
double mu_P0;//initial param. mu_P for R(x) 
double scalepeak;

void calBeta14(TString filename)
{
  TFile *ff = new TFile(filename);
  //original histo, Spatial Resolution
  TTree *Tdata = (TTree*)ff->Get("T");
  TH1F* hNfit_temp = new TH1F("hNfit_temp","rebin ", binTot, lowBin, hiBin);
  double posx, posy, posz, beta14, energy, posRad, itr, dirx, diry, dirz;
  UInt_t nhits;
  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("dirx", &dirx);
  Tdata->SetBranchAddress("diry", &diry);
  Tdata->SetBranchAddress("dirz", &dirz);
  Tdata->SetBranchAddress("posRad", &posRad);
  Tdata->SetBranchAddress("itr", &itr);
  Tdata->SetBranchAddress("beta14", &beta14);

  for(int i =0;i<Tdata->GetEntries();i++)
  {
    Tdata->GetEntry(i);
    //cout<<energy<<" "<<nhits<<endl;
    //if(nhits>5)
    if(nhits>5 && itr>0.55)// NO posRad cuts!!
    {
      hNfit_temp->Fill(beta14);
    }
  }
  double mu0 = hNfit_temp->GetMean();
  double range11 = mu0 - 1;
  double range22 = mu0 + 1;

  TF1 *fgaus = new TF1("fgaus","gaus",range11, range22);
  hNfit_temp->Fit("fgaus","rq");
  double mean = fgaus->GetParameter(1);
  double meanErr= fgaus->GetParError(1);
  double sigma = fgaus->GetParameter(2);
  double sigmaErr = fgaus->GetParError(2);
  double ndf = fgaus->GetNDF();
  double chi2 = fgaus->GetChisquare();
  std::cout << std::setprecision(4)<<std::fixed;
  cout<<mean<<", "<<meanErr<<", "<<sigma<<", "<<sigmaErr<<", "<<chi2<<"/"<<ndf<<", "<<chi2/ndf<<endl;

}
