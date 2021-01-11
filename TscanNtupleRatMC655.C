#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>

void TscanNtupleRatMC655(const char* filename)
{
  TFile *fold= new TFile(filename);
  TTree *T = (TTree*)fold->Get("output");
  //((TTreePlayer*)(T->GetPlayer()))->SetScanRedirect(true);
  int bin = 200;
  TH1F *hsolar = new TH1F("hsolar","solar #nu candidates, 5.0 < T < 15 MeV", bin, -1,1);
  TH1F *hsolar_5p5 = new TH1F("hsolar_5p5","solar #nu candidates, 5.5 < T < 15 MeV", bin, -1,1);
  TH1F *hsolar_6 = new TH1F("hsolar_6","solar #nu candidates, 6 < T < 15 MeV",bin, -1,1);
  TH2F *hsolarVsE = new TH2F("hsolarVsE","solar #nu candidates vs energy, 5 < T < 15 MeV", bin, -1,1, 10, 5,15);

  double posx, posy, posz, dirx, diry, dirz, itr,beta14, energy;
  Int_t nhits, eventID, triggerWord, uTDays, uTSecs, uTNSecs;
  ULong64_t dcFlagged, runID;
  bool fitValid, waterFit; 
  RAT::DU::ReconCorrector *eCorr = RAT::DU::ReconCorrector::Get();
  T->SetBranchAddress("posx",&posx);
  T->SetBranchAddress("posy",&posy);
  T->SetBranchAddress("posz",&posz);
  T->SetBranchAddress("dirx",&dirx);
  T->SetBranchAddress("diry",&diry);
  T->SetBranchAddress("dirz",&dirz);
  T->SetBranchAddress("itr",&itr);
  T->SetBranchAddress("beta14",&beta14);
  T->SetBranchAddress("nhits",&nhits);
  T->SetBranchAddress("energy",&energy);
  T->SetBranchAddress("fitValid",&fitValid);
  T->SetBranchAddress("waterFit",&waterFit);
  T->SetBranchAddress("dcFlagged",&dcFlagged);
  T->SetBranchAddress("eventID",&eventID);
  T->SetBranchAddress("triggerWord",&triggerWord);
  T->SetBranchAddress("uTDays",&uTDays);
  T->SetBranchAddress("uTSecs",&uTSecs);
  T->SetBranchAddress("uTNSecs",&uTNSecs);

  //read all entries and fill the histograms
  Long64_t nentries = T->GetEntries();

  TString fname(filename);
  for (Long64_t i=0;i<nentries;i++) {
    T->GetEntry(i);
    if(fitValid && waterFit) {
    //if( (dcFlagged & 0xFB0000017FFE) != 0xFB0000017FFE ) continue;
    //if( (dcFlagged & 0x10000017FFE) != 0x10000017FFE ) continue; // for rat6.5.3
    //else
    {	    
      double rCor = sqrt(posx*posx+posy*posy+(posz-108)*(posz-108));
      
      if( rCor<=5300 && nhits>30 && energy>=5 && energy<=15 && itr>=0.55 && beta14<=0.95 && beta14>=-0.12 && !( !(triggerWord & 0x3F) || (triggerWord & 0xBFF9400)) )
      {
       TVector3 directionSun = RAT::SunDirection(uTDays,uTSecs,uTNSecs);
       double sunDirX = directionSun.X();
       double sunDirY = directionSun.Y();
       double sunDirZ = directionSun.Z();
       double cosThetaToSun = -sunDirX*dirx+-sunDirY*diry+-sunDirZ*dirz;
       hsolar->Fill(cosThetaToSun);
       hsolarVsE->Fill(cosThetaToSun,energy);
       if(energy>5.5) hsolar_5p5->Fill(cosThetaToSun);
       if(energy>6.0) hsolar_6->Fill(cosThetaToSun);
       //cout<<fname(14,6)<<" "<<fname(21,4)<<" "<<fname(26,4)<<
      }
    //T->Scan("runID:eventID:triggerWord:energy:nhits:itr:beta14:posx:posy:posz-108:dirx:diry:dirz:sqrt(posx*posx+posy*posy+(posz-108)*(posz-108)):uTDays:uTSecs:uTNSecs","nhits>30 && itr>=0.55 && beta14>=-0.12 && beta14<=0.95 && sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))<=5300 && !( !(triggerWord & 0x3F) || (triggerWord & 0xBFF9400) ) && ( ((dcApplied & 0xFB0000017FFE) & dcFlagged) == (dcApplied & 0xFB0000017FFE) ) && energy>5 && energy<15");
    }
   }
  }
//  TString ss(filename);
//  ss = ss(0,ss.Index('.'));
//  ss = "tot"+ss+".txt";
//  ((TTreePlayer*)(T->GetPlayer()))->SetScanFileName(ss);
//
  TString fnewName = "saveMCsolar"+fname;
  TFile *ff = new TFile(fnewName, "recreate");  
  ff->cd();
  hsolar->Write();
  hsolarVsE->Write();
  hsolar_5p5->Write();
  hsolar_6->Write();
  ff->Close();
}
