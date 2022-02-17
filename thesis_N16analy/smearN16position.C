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

const int binTot = 1200;
double hiBin = -6000;
double lowBin = 6000;
double binStep = (hiBin-lowBin)/binTot;
// parameters
double xsigma = 113.6; 
double ysigma = 90.99;
double zsigma = 145.56;

// radius smear up
double deltaXpos = +0.07;
double deltaYpos = +0.02;
double deltaZpos = +0.08;

// radius smear down
double deltaXminus = -0.06;
double deltaYminus = -0.07;
double deltaZminus = -0.01;

TFile *fMC = new TFile("../mpwScan/Merged_MPWn16new_FitMPtest_nhit6_AntiNeutrinoMC_WaterN16sourceRun_r107055_s0.root");//
TFile *fdata = new TFile("../mpwScan/Merged_MPWn16new_WaterMP6176_nhit6_Calibration_r0000107055_p009.root");

TH1F *extractPosRad(TFile* ff, TString hname, int smearType)
{

  TFile *file = ff;
  TTree *Tdata = (TTree*)file->Get("T");
  /// posRad already has z-108 correction!!
  double posx, posy, posz, energy, posRad, itr, beta14, scaleLogL, Gtest, Utest, medianProbHit, medianProb, medianDevHit, medianDev;
  UInt_t nhits;
  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
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
 
  TH1F* hposx = new TH1F(hname,"",binTot, lowBin, hiBin);
  TH1F* hposy = new TH1F(hname,"",binTot, lowBin, hiBin);
  TH1F* hposz = new TH1F(hname,"",binTot, lowBin, hiBin);
  TH1F* hposRad = new TH1F(hname,"",2000, 0, 6000);

  TRandom3 rr;
  switch(smearType) {
   case 0:// no smearing
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            if(nhits>20 && posRad<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
	      hposRad->Fill(posRad);
            }
          }
	break;
   case 1: // smearing sigmaX 
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double xSmear = posx + rr.Gaus(0,xsigma);
            double posRadSmear = sqrt(xSmear*xSmear + posy*posy + posz*posz);
            if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(xSmear);
              hposy->Fill(posy);
              hposz->Fill(posz);
	      hposRad->Fill(posRadSmear);
	    }
        }
        break;
   case 2: // smearing sigmaY 
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double ySmear = posy + rr.Gaus(0,ysigma);
            double posRadSmear = sqrt(posx*posx + ySmear*ySmear + posz*posz);
            if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(ySmear);
              hposy->Fill(posy);
              hposz->Fill(posz);
              hposRad->Fill(posRadSmear);
            }
        }
        break;
   case 3: // smearing sigmaZ
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double zSmear = posz + rr.Gaus(0,zsigma);
            double posRadSmear = sqrt(posx*posx + posy*posy + zSmear*zSmear);
            if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(zSmear);
              hposRad->Fill(posRadSmear);
	    }
	}
        break;
   case 4: // smearing radius up 
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
	    double r2 = (posx*posx+posy*posy+(posz-108)*(posz-108));
            double deltaRpos = sqrt((posx*posx*deltaXpos*deltaXpos+posy*posy*deltaYpos*deltaYpos+posz*posz*deltaZpos*deltaZpos)/r2);
            double posRadSmear = (1+deltaRpos/100)*posRad;// smear up
            if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
              hposRad->Fill(posRadSmear);
            }
	}
        break;
   case 5: // smearing radius down
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double r2 = (posx*posx+posy*posy+(posz-108)*(posz-108));
            double deltaRminus = sqrt((posx*posx*deltaXminus*deltaXminus+posy*posy*deltaYminus*deltaYminus+posz*posz*deltaZminus*deltaZminus)/r2);
            double posRadSmear = (1+deltaRminus/100)*posRad;// smear up
	    if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
              hposRad->Fill(posRadSmear);
            }
	}
        break;

   default:
  }
  return hposRad; 
}

void smearN16position()
{
  TH1F *hPosRad = extractPosRad(fdata,"posR",0);
  TH1F *hPosRadMC = extractPosRad(fMC,"posR MC",0);
  TH1F *hPosRadposXsmear = extractPosRad(fMC,"hposXsmear",1); 
  TH1F *hPosRadposYsmear = extractPosRad(fMC,"hposYsmear",2);
  TH1F *hPosRadposZsmear = extractPosRad(fMC,"hposZsmear",3);
  TH1F *hPosRadSmearUp = extractPosRad(fMC,"hposRadSmearUp",4);
  TH1F *hPosRadSmearDown = extractPosRad(fMC,"hposRadSmearDown",5);

  double scale = hPosRad->Integral();
  hPosRad->GetXaxis()->SetTitle("radius [mm]");
  hPosRad->GetYaxis()->SetTitle("scaled");
  TCanvas *c = new TCanvas("c","",800,600);

  hPosRadMC->Scale(scale/hPosRadMC->Integral());hPosRadMC->SetLineColor(kGray+2);hPosRadMC->SetLineWidth(2);hPosRadMC->SetLineStyle(4);

  hPosRadMCScaleUp->Scale(scale/hPosRadMCScaleUp->Integral());
  hPosRadMCScaleDown->Scale(scale/hPosRadMCScaleDown->Integral());
  hPosRadMCResolDown->Scale(scale/hPosRadMCResolDown->Integral());

  c->cd();
  hPosRad->Draw();
  hPosRadMC->Draw("same");
  hPosRadMCScaleUp->SetLineColor(kRed);hPosRadMCScaleUp->SetLineStyle(2);
  hPosRadMCScaleDown->SetLineColor(kBlue);hPosRadMCScaleDown->SetLineStyle(3);
  hPosRadMCResolDown->SetLineColor(kGreen+2);hPosRadMCResolDown->SetLineStyle(3);
  hPosRadMCScaleUp->Draw("same");
  hPosRadMCScaleDown->Draw("same");
//  hPosRadMCResolDown->Draw("same");
  TLegend *legend1 = new TLegend(0.1,0.7,0.35,0.9);
  legend1->AddEntry(hPosRad,"data","l");
  legend1->AddEntry(hPosRadMC,"MC","l");
  legend1->AddEntry(hPosRadMCScaleUp,"E_{scale} up","l");
  legend1->AddEntry(hPosRadMCScaleDown,"E_{scale} down","l");
  legend1->Draw("same");

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  // hPosRad->GetXaxis()->SetRangeUser(0,100);
  hPosRad->Draw();
  hPosRadMC->Draw("same");
  hPosRadMCResolDown->SetLineColor(kRed);hPosRadMCResolDown->SetLineStyle(3);
  hPosRadMCResolDown->Draw("same");
//  hPosRadMCResolDown->Draw("same");
  TLegend *legend2 = new TLegend(0.1,0.7,0.35,0.9);
  legend2->AddEntry(hPosRad,"data","l");
  legend2->AddEntry(hPosRadMC,"MC","l");
  legend2->AddEntry(hPosRadMCResolDown,"Smeared E_{resol}","l");
  legend2->Draw("same");

}
