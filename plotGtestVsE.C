#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
#include <RAT/DataCleaningUtility.hh>
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
#include "Math/MinimizerOptions.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <TMinuit.h>
#include <algorithm>
//!!!!!!!!!!!!!!! must use 6.18.8 for energy correction !!!!!!!!!!!!!!!!!
void plotGtestVsE(const char* fname)
{
  TVector3 srcPos(-5.283,-0.209,-1.057);
  TFile *fdata = new TFile(fname);
  double Elow = 0;
  double Ehigh = 12;
  int binNum = (Ehigh-Elow)/0.05;
  TH1F *h1Erec = new TH1F("h1Erec","",binNum,Elow, Ehigh);
  RAT::DU::ReconCorrector *eCorr = RAT::DU::ReconCorrector::Get();
  TTree *Tdata = (TTree*)fdata->Get("T");
  double posx, posy, posz, dirx, diry, dirz, gtest, energyMC, energyOrigin, posRad, itr, beta14, scaleLogL;
  UInt_t nhits;
  Tdata->SetBranchAddress("energyOrigin", &energyOrigin);////!!!! original energy!!!
  Tdata->SetBranchAddress("energyMC", &energyMC);////!!!! original energy!!!
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("posRad", &posRad);
  Tdata->SetBranchAddress("itr", &itr);
  Tdata->SetBranchAddress("beta14", &beta14);

  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
  Tdata->SetBranchAddress("scaleLogL", &scaleLogL);

  Tdata->SetBranchAddress("dirx", &dirx);
  Tdata->SetBranchAddress("diry", &diry);
  Tdata->SetBranchAddress("dirz", &dirz);
  Tdata->SetBranchAddress("Gtest", &gtest);

  TH2F *hGtestVsE = new TH2F("hGtestVsE","Gtest vs E_{recon} - E_{MC}",400,-1.5,2.5,500,0,5);
  TH2F *hGtestVsE_letaCuts = new TH2F("hGtestVsE_letaCuts","Gtest vs E_{recon} - E_{MC}",500,0,5,400,-1.5,2.5);

  for(int i =0;i<Tdata->GetEntries();i++)
  {
    Tdata->GetEntry(i);
    TVector3 evtPos(posx,posy,posz-108);
    TVector3 Dir(dirx,diry,dirz);
    TVector3 Xdiff = (evtPos - srcPos);
    double distance = (evtPos - srcPos).Mag();
    //cout<<energy<<" "<<nhits<<endl;
    double energy = eCorr->CorrectEnergyRSP(energyOrigin,2);//!!! rat-6.18.x correction

    hGtestVsE->Fill(energy-energyMC, gtest);
    if(nhits>5 && itr>0.55 && beta14<0.95 && beta14>-0.12 && posRad<5300 && scaleLogL>10.5)
    {
      /// LETA cuts!!
      //cout<<energy<<" "<<energyOrigin<<endl;
      if(distance>=700) 
      {
         h1Erec->Fill(energy);
         hGtestVsE_letaCuts->Fill(energy-energyMC, gtest);
      }
      else
      {
        double costheta = Xdiff.Unit()*Dir;
        if(costheta>=sqrt(2)/2 && costheta<=1) 
	{
           h1Erec->Fill(energy);
           hGtestVsE_letaCuts->Fill(energy-energyMC, gtest);
	}
      }
    }
  }
  TString ss(fname);
  TString title("Erecon_");
  TString newfile = title+ss;
  cout<<newfile<<endl;
  TFile *fnew = new TFile(newfile,"recreate");
  fnew->cd();
  h1Erec->Write();
  hGtestVsE->Write();
  hGtestVsE_letaCuts->Write();
}
