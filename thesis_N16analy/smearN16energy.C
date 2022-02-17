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

const int binTot = 120;
double hiBin = 12;
double lowBin = 0;
double binStep = (hiBin-lowBin)/binTot;
// parameters
double Escale = 1.0/100;
double Eresol = 0.0369;

TFile *fMC = new TFile("../mpwScan/Merged_MPWn16new_FitMPtest_nhit6_AntiNeutrinoMC_WaterN16sourceRun_r107055_s0.root");//
TFile *fdata = new TFile("../mpwScan/Merged_MPWn16new_WaterMP6176_nhit6_Calibration_r0000107055_p009.root");

TH1F *extractEnergy(TFile* ff, TString hname, int smearType)
{

  TFile *file = ff;
  TTree *Tdata = (TTree*)file->Get("T");
  double energy, posRad, itr, beta14, scaleLogL, Gtest, Utest, medianProbHit, medianProb, medianDevHit, medianDev;
  UInt_t nhits;
  Tdata->SetBranchAddress("energy", &energy);
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("posRad", &posRad);
  Tdata->SetBranchAddress("itr", &itr);
  Tdata->SetBranchAddress("scaleLogL", &scaleLogL);
  Tdata->SetBranchAddress("beta14", &beta14);
  Tdata->SetBranchAddress("Gtest", &Gtest);
  Tdata->SetBranchAddress("Utest", &Utest);
  Tdata->SetBranchAddress("medianProbHit", &medianProbHit);
  Tdata->SetBranchAddress("medianProb", &medianProb);
  Tdata->SetBranchAddress("medianDevHit", &medianDevHit);
  Tdata->SetBranchAddress("medianDev", &medianDev);
 
  TH1F* hEnergy = new TH1F(hname,"",binTot, lowBin, hiBin);
  TRandom3 rr;
  switch(smearType) // no smearing
  {
   case 0:	  
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);

            //cout<<energy<<" "<<nhits<<endl;
            if(nhits>20 && posRad<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hEnergy->Fill(energy);
            }
          }
	break;
   case 1: // smearing Escale up
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);

            double energySmear = energy*(1+Escale);
            if(nhits>20 && posRad<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hEnergy->Fill(energySmear);
            }
          }
        break;
   case 2: // smearing Escale down
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double energySmear = energy*(1-Escale);
            //cout<<energy<<" "<<nhits<<endl;
            if(nhits>20 && posRad<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
               hEnergy->Fill(energySmear);
            }
          }
        break;
   case 3: // smearing Eresol
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
	    double sigma = sqrt(energy)*sqrt((1+Eresol)*(1+Eresol)-1);
	    //    wrong calculation => double sigma = sqrt(energy*((1+Eresol)*(1+Eresol))-1);
            double energySmear = energy+rr.Gaus(0,sigma); // E convolved with Gauss(0,sigma_smear)   similar: double energySmear = rr.Gaus(energy, sigma);
	    cout<<sigma<<" "<<energy<<" "<<energySmear<<nhits<<endl;
            if(nhits>20 && posRad<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hEnergy->Fill(energySmear);
            }
          }
        break;
   default:
  }
  return hEnergy; 
}

void smearN16energy()
{
  TH1F *hEdata = extractEnergy(fdata,"hEdata",0);
  TH1F *hEmc = extractEnergy(fMC,"hEmc",0);
  TH1F *hEmcScaleUp = extractEnergy(fMC,"hEmcScaleUp",1); 
  TH1F *hEmcScaleDown = extractEnergy(fMC,"hEmcScaleDown",2);
  TH1F *hEmcResolDown = extractEnergy(fMC,"hEmcResolDown",3);

  double scale = hEdata->Integral();
  hEdata->GetXaxis()->SetTitle("E_{fit} [MeV]");
  hEdata->GetYaxis()->SetTitle("scaled");
  TCanvas *c = new TCanvas("c","",800,600);

  hEmc->Scale(scale/hEmc->Integral());hEmc->SetLineColor(kGray+2);hEmc->SetLineWidth(2);hEmc->SetLineStyle(4);

  hEmcScaleUp->Scale(scale/hEmcScaleUp->Integral());
  hEmcScaleDown->Scale(scale/hEmcScaleDown->Integral());
  hEmcResolDown->Scale(scale/hEmcResolDown->Integral());

  c->cd();
  hEdata->Draw();
  hEmc->Draw("same");
  hEmcScaleUp->SetLineColor(kRed);hEmcScaleUp->SetLineStyle(2);
  hEmcScaleDown->SetLineColor(kBlue);hEmcScaleDown->SetLineStyle(3);
  hEmcResolDown->SetLineColor(kGreen+2);hEmcResolDown->SetLineStyle(3);
  hEmcScaleUp->Draw("same");
  hEmcScaleDown->Draw("same");
//  hEmcResolDown->Draw("same");
  TLegend *legend1 = new TLegend(0.1,0.7,0.35,0.9);
  legend1->AddEntry(hEdata,"data","l");
  legend1->AddEntry(hEmc,"MC","l");
  legend1->AddEntry(hEmcScaleUp,"E_{scale} up","l");
  legend1->AddEntry(hEmcScaleDown,"E_{scale} down","l");
  legend1->Draw("same");

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  // hEdata->GetXaxis()->SetRangeUser(0,100);
  hEdata->Draw();
  hEmc->Draw("same");
  hEmcResolDown->SetLineColor(kRed);hEmcResolDown->SetLineStyle(3);
  hEmcResolDown->Draw("same");
//  hEmcResolDown->Draw("same");
  TLegend *legend2 = new TLegend(0.1,0.7,0.35,0.9);
  legend2->AddEntry(hEdata,"data","l");
  legend2->AddEntry(hEmc,"MC","l");
  legend2->AddEntry(hEmcResolDown,"Smeared E_{resol}","l");
  legend2->Draw("same");

}
