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
#include <TVector3.h>
#include <TNtuple.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
void reprocessNtupleSolarRat()
{
    const char* filename = "Analysis_r0000204288_s000_p003.ntuple.root";
//    double mcPosx, mcPosy, mcPosz, mcEdep;
//    double mcmom1x, mcmom1y, mcmom1z;
    double posx, posy, posz, posRad, time, energy, beta14, thetaij, itr, dirx, diry, dirz;
    double posFOM, Gtest, Utest, medProb, medProbHit, medDev, medDevHit;
    UInt_t posFOM2;
    bool waterFit, fitValid;
    Int_t nhits;
    Int_t triggerWord, eventID, runID;
    ULong64_t dcFlagged, dcApplied;
    Int_t uTDays, uTSecs, uTNSecs;

    TFile *fff = new TFile(filename);
    // Load the RAT file
    TTree *tree = (TTree*)fff->Get("output");
    //tree->SetBranchAddress("mcPosx",&mcPosx0);tree->SetBranchAddress("mcPosy",&mcPosy0);tree->SetBranchAddress("mcPosz",&mcPosz0);
    //tree->SetBranchAddress("mcmom1x",&mcmom1x0);tree->SetBranchAddress("mcmom1y",&mcmom1y0);tree->SetBranchAddress("mcmom1z",&mcmom1z0);
    tree->SetBranchAddress("posx", &posx);tree->SetBranchAddress("posy", &posy);tree->SetBranchAddress("posz", &posz);
    tree->SetBranchAddress("dirx", &dirx);tree->SetBranchAddress("diry", &diry);tree->SetBranchAddress("dirz", &dirz);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("beta14", &beta14);
    tree->SetBranchAddress("thetaij", &thetaij);
    tree->SetBranchAddress("itr", &itr);
    tree->SetBranchAddress("nhitsCleaned", &nhits);

    tree->SetBranchAddress("posFOM", &posFOM);
    tree->SetBranchAddress("posFOM2", &posFOM2);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("energyFOMGtest", &Gtest);
    tree->SetBranchAddress("energyFOMUtest", &Utest);
    tree->SetBranchAddress("energyFOMmedProb", &medProb);
    tree->SetBranchAddress("energyFOMmedProbHit", &medProbHit);
    tree->SetBranchAddress("energyFOMmedDev", &medDev);
    tree->SetBranchAddress("energyFOMmedDevHit", &medDevHit);
    tree->SetBranchAddress("triggerWord", &triggerWord);
    tree->SetBranchAddress("dcFlagged", &dcFlagged);
    tree->SetBranchAddress("dcApplied", &dcApplied);

    tree->SetBranchAddress("uTDays", &uTDays);
    tree->SetBranchAddress("uTSecs", &uTSecs);
    tree->SetBranchAddress("uTNSecs",&uTNSecs);

    tree->SetBranchAddress("fitValid", &fitValid);
    tree->SetBranchAddress("waterFit", &waterFit);
    //MultiPath
    //double posXmp, posYmp, posZmp;

    // Processed new Tree

    double energyCor, sunDirx, sunDiry, sunDirz, cosThetaSun, zfactor, logL, scaledLogL;
    TTree *TT = new TTree("T","solar processed");
    TT->Branch("runID", &runID,"runID/i");
    TT->Branch("eventID", &eventID,"eventID/i");
    TT->Branch("beta14", &beta14, "beta14/D");
    TT->Branch("thetaij", &thetaij, "thetaij/D");
    TT->Branch("energy", &energyCor, "energyCor/D");
    TT->Branch("itr", &itr, "itr/D");
    TT->Branch("nhits", &nhits,"nhits/i");//!!! nhitsCleaned
//    TT->Branch("mcPosx", &mcPosx,"mcPosx/D");TT->Branch("mcPosy", &mcPosy,"mcPosy/D");TT->Branch("mcPosz", &mcPosz,"mcPosz/D");
//    TT->Branch("mcEdep", &mcEdep, "mcEdep");
//    TT->Branch("mcmom1x", &mcmom1x,"mcmom1x/D");TT->Branch("mcmom1y", &mcmom1y,"mcmom1y/D");TT->Branch("mcmom1z", &mcmom1z,"mcmom1z/D");
    TT->Branch("posx", &posx,"posx/D");TT->Branch("posy", &posy,"posy/D");TT->Branch("posz", &posz,"posz/D");
    TT->Branch("dirx", &dirx,"dirx/D");TT->Branch("diry", &diry,"diry/D");TT->Branch("dirz", &dirz,"dirz/D");
    TT->Branch("posLogL", &posFOM,"posLogL/D");
    TT->Branch("scaledLogL", &scaledLogL,"scaledLogL/D");

    TT->Branch("energyFOMGtest", &Gtest,"Gtest/D");
    TT->Branch("energyFOMUtest", &Utest,"Utest/D");
    TT->Branch("zfactor", &zfactor,"zfactor/D");

    TT->Branch("triggerWord", &triggerWord,"triggerWord/i");
    TT->Branch("dcFlagged", &dcFlagged,"dcFlagged/l");
    TT->Branch("dcApplied", &dcApplied,"dcApplied/l");

    TT->Branch("cosThetaSun", &cosThetaSun,"cosThetaSun/D");

    TString ss(filename);
    TString newname = "Processed_"+ss;

    TFile *file=new TFile(newname, "RECREATE");
    RAT::DU::ReconCorrector *eCorr = RAT::DU::ReconCorrector::Get();
    //const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();

    for(int i = 0;i<tree->GetEntries();i++)
    {
       tree->GetEntry(i);
       if( nhits>20 && fitValid && waterFit )
       {
         zfactor = (medProbHit != medProb) ? (1 - 3*(medDevHit + medDev)/(medProbHit - medProb)) : -99999;
         logL = posFOM;
	 scaledLogL = posFOM/posFOM2;
	 logL = posFOM;

	 TVector3 sunDir = RAT::SunDirection(uTDays, uTSecs, uTNSecs);
	 TVector3 evtDir(dirx, diry, dirz);
	 cosThetaSun = evtDir*-sunDir;

	 energyCor = eCorr->CorrectEnergyRSP(energy,2);
	 TT->Fill();
       }
    }
   // Write the histograms to fileName
   file->cd();
   TT->Write();
   file->Close();

}
