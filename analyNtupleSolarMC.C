#include <RAT/DU/DSReader.hh>
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

void analyNtupleSolarMC()
{

  ifstream in0, in1;
  ofstream out,outputbadfiles;
  ostringstream oss;

  in0.open("fList.dat");
  char filename[100];
  while(in0>>filename){
    double mcPosx, mcPosy, mcPosz;
    double mcmom1x, mcmom1y, mcmom1z;
    double posx, posy, posz, posRad, time, energy, beta14, thetaij, dirx, diry, dirz, nhits, sunDirx, sunDiry, sunDirz, cosThetaToSun, itr;
    double posFOM, posFOM2;
    //MultiPath
    //double posXmp, posYmp, posZmp;
    int evtID;
    int runID;
    TTree *tree = new TTree("TsolarMC","solar MC");
    tree->Branch("runID", &runID,"runID/I");
    tree->Branch("eventID", &evtID,"eventID/I");
    tree->Branch("beta14", &beta14, "beta14/D");
    tree->Branch("thetaij", &thetaij, "thetaij/D");
    //tree->Branch("energy", energy, "energy/D");
    tree->Branch("itr", &itr, "itr/D");
    tree->Branch("nhits", &nhits,"nhits/D");

    tree->Branch("mcPosx", &mcPosx,"mcPosx/D");tree->Branch("mcPosy", &mcPosy,"mcPosy/D");tree->Branch("mcPosz", &mcPosz,"mcPosz/D");

    tree->Branch("mcmom1x", &mcmom1x,"mcmom1x/D");tree->Branch("mcmom1y", &mcmom1y,"mcmom1y/D");tree->Branch("mcmom1z", &mcmom1z,"mcmom1z/D");
    tree->Branch("posx", &posx,"posx/D");tree->Branch("posy", &posy,"posy/D");tree->Branch("posz", &posz,"posz/D");
    tree->Branch("dirx", &dirx,"dirx/D");tree->Branch("diry", &diry,"diry/D");tree->Branch("dirz", &dirz,"dirz/D");

    tree->Branch("posFOM", &posFOM,"posFOM/D");
    tree->Branch("posFOM2", &posFOM2,"posFOM2/D");

    Int_t i=0;
    Int_t TriggerType;
    Int_t nentries=0;
    ULong64_t bitword;

    cout<<" filename "<<filename<<endl ;
    double mcPosx0, mcPosy0, mcPosz0;
    double mcmom1x0, mcmom1y0, mcmom1z0;
    double posx0, posy0, posz0, posRad0, time0, energy0, beta140, thetaij0, dirx0, diry0, dirz0, itr0;
    double posFOM0;
    int nhits0;
    UInt_t posFOM20;
    int evtID0;
    int runID0;
    ULong64_t dcFlagged;
    bool waterFit, fitValid;

    TFile *fff = new TFile(filename);
    // Load the RAT file
    TTree *treeTempt = (TTree*)fff->Get("output");
    treeTempt->SetBranchAddress("mcPosx",&mcPosx0);treeTempt->SetBranchAddress("mcPosy",&mcPosy0);treeTempt->SetBranchAddress("mcPosz",&mcPosz0);
    treeTempt->SetBranchAddress("mcmom1x",&mcmom1x0);treeTempt->SetBranchAddress("mcmom1y",&mcmom1y0);treeTempt->SetBranchAddress("mcmom1z",&mcmom1z0);

    treeTempt->SetBranchAddress("posx", &posx0);treeTempt->SetBranchAddress("posy", &posy0);treeTempt->SetBranchAddress("posz", &posz0);
    treeTempt->SetBranchAddress("dirx", &dirx0);treeTempt->SetBranchAddress("diry", &diry0);treeTempt->SetBranchAddress("dirz", &dirz0);
    treeTempt->SetBranchAddress("runID", &runID0);
    treeTempt->SetBranchAddress("eventID", &evtID0);
    treeTempt->SetBranchAddress("beta14", &beta140);
    treeTempt->SetBranchAddress("thetaij", &thetaij0);
    treeTempt->SetBranchAddress("itr", &itr0);
    treeTempt->SetBranchAddress("nhits", &nhits0);
    treeTempt->SetBranchAddress("posFOM", &posFOM0);
    treeTempt->SetBranchAddress("posFOM2", &posFOM20);

    treeTempt->SetBranchAddress("dcFlagged", &dcFlagged);
    treeTempt->SetBranchAddress("fitValid", &fitValid);
    treeTempt->SetBranchAddress("waterFit", &waterFit);

    TString s(filename);
    TString runName = s(s.Index("r20")+1,  s.Index("r20")+1+6);
    TString newname = "MCTl208external_500files_"+runName;
    TFile *file=new TFile(newname, "RECREATE");


    for(int i = 0;i<treeTempt->GetEntries();i++)
    {
       treeTempt->GetEntry(i);  
       if( int(dcFlagged/pow(2,32))==196992 && fitValid && waterFit )//&& sqrt(posx0*posx0+posy0*posy0+(posz0-108)*(posz0-108))<8390)
       {
         runID = runID0;
	 evtID = evtID0;
         beta14 = beta140;
	 thetaij = thetaij0;
	 nhits = nhits0;
	 itr = itr0; 
         posx = posx0; posy = posy0; posz = posz0; dirx = dirx0; diry = diry0; dirz = dirz0;
         mcmom1x = mcmom1x0; mcmom1y = mcmom1y0; mcmom1z = mcmom1z0;
	 mcPosx = mcPosx0; mcPosy = mcPosy0; mcPosz = mcPosz0;
	 posFOM = posFOM0; posFOM2 = posFOM20;
         tree->Fill();  
       }
    }
   delete treeTempt;
    // Write the histograms to fileName
   fff->Close();
 
   file->cd();
   tree->Write();
   file->Close();
   delete tree;
  }

}
