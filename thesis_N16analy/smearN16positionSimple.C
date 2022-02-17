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

const int binTot = 200;
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

TFile *ffsave = new TFile("saveSmearPosition.root","recreate");

TH1F *extractPosRad(TFile* ff, TString hname, int smearType, vector<TH1F*>& saveHist)
{

  TTree *Tdata = (TTree*)ff->Get("T");
  /// !!!!!! posRad already has z-108 correction!!
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
 
  TH1F* hposx = new TH1F("hposx","",binTot, lowBin, hiBin);
  TH1F* hposy = new TH1F("hposy","",binTot, lowBin, hiBin);
  TH1F* hposz = new TH1F("hposz","",binTot, lowBin, hiBin);
  TH1F* hposRad = new TH1F("hposRad","",200, 0, 6000);

  TRandom3 rr;
  switch(smearType) {
   case -1:// no smearing, data
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)//Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
              hposRad->Fill(posRad);
            }
        }
        hposx->SetName("hposx_origin");
        hposy->SetName("hposy_origin");
        hposz->SetName("hposz_origin");
        hposRad->SetName("hposRad_origin");
        saveHist.push_back(hposx);
        saveHist.push_back(hposy);
        saveHist.push_back(hposz);
        saveHist.push_back(hposRad);
        break;

   case 0:// no smearing, MC
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)//Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
	      hposRad->Fill(posRad);
	    }
        }
	hposx->SetName("hposx_MCorigin");
        hposy->SetName("hposy_MCorigin");
        hposz->SetName("hposz_MCorigin");
        hposRad->SetName("hposRad_MCorigin");
        saveHist.push_back(hposx);
        saveHist.push_back(hposy);
        saveHist.push_back(hposz);
        saveHist.push_back(hposRad);
	break;
   case 1: // smearing sigmaX 
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double xSmear = posx + rr.Gaus(0,xsigma);
            double posRadSmear = sqrt(xSmear*xSmear + posy*posy + posz*posz);
            // if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)
	    {
              hposx->Fill(xSmear);
              hposy->Fill(posy);
              hposz->Fill(posz);
	      hposRad->Fill(posRadSmear);
	    }
        }
        hposx->SetName("hposx_smearX");
        hposRad->SetName("hposRad_smearX");
        saveHist.push_back(hposx);
        saveHist.push_back(hposRad);
        break;
   case 2: // smearing sigmaY 
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double ySmear = posy + rr.Gaus(0,ysigma);
            double posRadSmear = sqrt(posx*posx + ySmear*ySmear + posz*posz);
            // if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)
	    {
              hposx->Fill(posx);
              hposy->Fill(ySmear);
              hposz->Fill(posz);
              hposRad->Fill(posRadSmear);
            }
        }
        hposy->SetName("hposy_smearY");
        hposRad->SetName("hposRad_smearY");
        saveHist.push_back(hposy);
        saveHist.push_back(hposRad);
	break;
   case 3: // smearing sigmaZ
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double zSmear = posz + rr.Gaus(0,zsigma);
            double posRadSmear = sqrt(posx*posx + posy*posy + zSmear*zSmear);
            // if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)
	    {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(zSmear);
              hposRad->Fill(posRadSmear);
	    }
	}
        hposz->SetName("hposz_smearZ");
        hposRad->SetName("hposRad_smearZ");
        saveHist.push_back(hposz);
        saveHist.push_back(hposRad);
	break;
   case 4: // smearing radius up 
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
	    double r2 = (posx*posx+posy*posy+(posz-108)*(posz-108)); // = posRad*posRad
            double deltaRpos = sqrt((posx*posx*deltaXpos*deltaXpos+posy*posy*deltaYpos*deltaYpos+(posz-108)*(posz-108)*deltaZpos*deltaZpos)/r2);
            double posRadSmear = (1+deltaRpos/100)*posRad;// smear up
            // if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)
	    {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
              hposRad->Fill(posRadSmear);
            }
	}
	//hposx->SetName("hposx_radiusUp");
        //hposy->SetName("hposy_radiusUp");
        //hposz->SetName("hposz_radiusUp");
        hposRad->SetName("hposRad_radiusUp");
        saveHist.push_back(hposRad);
	break;
   case 5: // smearing radius down
        for(int i =0;i<Tdata->GetEntries();i++)
        {
            Tdata->GetEntry(i);
            double zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
            double r2 = (posx*posx+posy*posy+(posz-108)*(posz-108)); // = posRad*posRad
            double deltaRminus = sqrt((posx*posx*deltaXminus*deltaXminus+posy*posy*deltaYminus*deltaYminus+(posz-108)*(posz-108)*deltaZminus*deltaZminus)/r2);
            double posRadSmear = (1-deltaRminus/100)*posRad;// smear down 
	    // if(nhits>20 && posRadSmear<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && zfactor>-11 && zfactor<1 && energy>1.5)
            if(nhits>5 && posRad<6000 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && energy>3.5)
	    {
              hposx->Fill(posx);
              hposy->Fill(posy);
              hposz->Fill(posz);
              hposRad->Fill(posRadSmear);
            }
	}
	//hposx->SetName("hposx_radiusDown");
        //hposy->SetName("hposy_radiusDown");
        //hposz->SetName("hposz_radiusDown");
        hposRad->SetName("hposRad_radiusDown");
        saveHist.push_back(hposRad);
	break;
   default:
  }
}

void smearN16positionSimple()
{

  vector<TH1F*> saveHistOrigin;
  vector<TH1F*> saveHistOriginMC;
  vector<TH1F*> saveHistSmearX;
  vector<TH1F*> saveHistSmearY;
  vector<TH1F*> saveHistSmearZ;
  vector<TH1F*> saveHistPosRadUp;
  vector<TH1F*> saveHistPosRadDown;


  extractPosRad(fdata,"posR",-1, saveHistOrigin);
  extractPosRad(fMC,"posR",0, saveHistOriginMC);

  extractPosRad(fMC,"posR",1, saveHistSmearX);
  extractPosRad(fMC,"posR",2, saveHistSmearY);
  extractPosRad(fMC,"posR",3, saveHistSmearZ);

  extractPosRad(fMC,"posR",4, saveHistPosRadUp);
  extractPosRad(fMC,"posR",5, saveHistPosRadDown);

  ffsave->cd(); 
  for(vector<TH1F*>::iterator it=saveHistOrigin.begin();it!=saveHistOrigin.end();it++)
  (*it)->Write();

  for(vector<TH1F*>::iterator it=saveHistOriginMC.begin();it!=saveHistOriginMC.end();it++)
  (*it)->Write();

  for(vector<TH1F*>::iterator it=saveHistSmearX.begin();it!=saveHistSmearX.end();it++)
  (*it)->Write();

  for(vector<TH1F*>::iterator it=saveHistSmearY.begin();it!=saveHistSmearY.end();it++)
  (*it)->Write();

  for(vector<TH1F*>::iterator it=saveHistSmearZ.begin();it!=saveHistSmearZ.end();it++)
  (*it)->Write();

  for(vector<TH1F*>::iterator it=saveHistPosRadUp.begin();it!=saveHistPosRadUp.end();it++)
  (*it)->Write();

  for(vector<TH1F*>::iterator it=saveHistPosRadDown.begin();it!=saveHistPosRadDown.end();it++)
  (*it)->Write();
}
