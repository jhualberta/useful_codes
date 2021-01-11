#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
#include <RAT/DataCleaningUtility.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
//DATA
void extractESpdf()
{
  TFile *fMC = new TFile("Merged_GeneralPhysicsMC_WaterSolar_NueRun_r200004to206391_s0_p0.ntuple.root");
  TH1F *hcosPdf = new TH1F("hcosPdf","nue ES pdf",200,-1,1);

  TH1F *hcosPdf_cuts = new TH1F("hcosPdf_cuts","nue ES pdf, FV<5.5 m, itr>0.55",200,-1,1);

  TTree *Tdata = (TTree*)fMC->Get("output");
  Int_t nhits, eventID;
  Int_t uTDays, uTSecs, uTNSecs;
  RAT::DU::ReconCorrector *eCorr = RAT::DU::ReconCorrector::Get();
  ULong64_t dcFlagged, runID;
  bool fitValid, waterFit;
  //  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
  double cosThetaToSun; 
  double posx, posy, posz, dirx, diry, dirz, energyOrigin, energy, posRad, itr, beta14, scaleLogL;
  Tdata->SetBranchAddress("energy", &energyOrigin);////!!!! original energy!!!
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("itr", &itr);
  Tdata->SetBranchAddress("beta14", &beta14);
  Tdata->SetBranchAddress("fitValid",&fitValid);
  Tdata->SetBranchAddress("waterFit",&waterFit);
  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
/// Tdata->SetBranchAddress("scaleLogL", &scaleLogL);
  Tdata->SetBranchAddress("dirx", &dirx);
  Tdata->SetBranchAddress("diry", &diry);
  Tdata->SetBranchAddress("dirz", &dirz);
/// Tdata->SetBranchAddress("cosThetaToSun", &cosThetaToSun);
  Tdata->SetBranchAddress("uTDays",&uTDays);
  Tdata->SetBranchAddress("uTSecs",&uTSecs);
  Tdata->SetBranchAddress("uTNSecs",&uTNSecs);

  for(int i =0;i<Tdata->GetEntries();i++)
  {
    Tdata->GetEntry(i);
    if(fitValid && waterFit) {
       TVector3 directionSun = RAT::SunDirection(uTDays,uTSecs,uTNSecs);
       TVector3 ufit(dirx,diry,dirz);
       TVector3 r(posx,posy,posz-108);
       double cosThetaToSun = -directionSun*ufit;
       energy = eCorr->CorrectEnergyRSP(energyOrigin,2);
       //cout<<energyOrigin<<" "<<energy<<endl;
       hcosPdf->Fill(cosThetaToSun);
       if(r.Mag()<5500 && itr>0.55)
       {
         hcosPdf_cuts->Fill(cosThetaToSun);
       }
    }
  }

  TFile *ff = new TFile("esPdf.root","recreate");
  ff->cd();
  hcosPdf->Write();
  hcosPdf_cuts->Write();
}
