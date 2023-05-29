//create a chain
#ifndef __CINT__ 
// include used ROOT classes below
#include "TROOT.h"
#include "TObjArray.h"
#include "TString.h"
#include "TRegexp.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TLine.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TMath.h"
// include used RAT classes below
#include <RAT/DS/EV.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/FitResult.hh>
// include STL classes below (C++ headers)
#include <iostream>
#include <cmath>
#include <stdint.h>
#include <string>
#include <RAT/GeoUtils.hh>
//using namespace std;

#endif 

#include <stdint.h>

void Analysis_coinc(TString filename, int lines, double FV, TString outputname) 
//filename = file with runs, one per line
//lines = total lines in filename
//FV = analysis volume, in mm
//output = name of the output file + directory

{
  //define common pattern in name of files to be merged, creat a virtual TChain
  char *indir =""; //set the directory of data here
  ifstream input;
  char name[200];
  input.open(filename);	
  int counter = 0;
  TObjString *fn = NULL;
  TObjArray *flist = new TObjArray();
  const char *file;
  TString fname;
  int size = lines;
  int j = 0;

  while(j<size)
    {
      input>>name;
      j++;
      string namen = name;
      //cout<<namen<<endl;
      TString s(namen);
      TString prefix1 = s(0, s.Index("_r00") );
      char *prefix = prefix1.Data(); //Analysis module, back in time it was Analysis30
      char *treename = "output";
      string pattern = prefix+namen;
      
      pattern = pattern + "_*";
      
      if (!indir) 
	{
	  cerr << "Failed to get directory name." << endl;
	  return;
	}

      gSystem->ChangeDirectory(indir);
      cout << "Working in directory : " << gSystem->WorkingDirectory() << endl;
      
      void *dir = gSystem->OpenDirectory(".");
      if (!dir) 
	{
	  cerr << "Couldn't open directory : " << indir << endl;
	  return;
	}
    
      TRegexp regexp(pattern.c_str());
      
      while ((file = gSystem->GetDirEntry(dir))) 
	{
	  fname = file;
	  //ignore subfolders and existing merged files
	  
	  if (!fname.CompareTo(".") || !fname.CompareTo("..") ||  (fname.Contains("_cc"))) continue;
	  if (fname.Index(regexp) != kNPOS) 
	    {
	      flist->Add(new TObjString(fname));
	    }
	}// end of while for add
    } // end of while size

  gSystem->FreeDirectory(dir);
  flist->Print("*");
  
  //loop over array and add file to the TChain
  TChain *tc = new TChain(treename,"");
  TIter next(flist);
  while ((fn = (TObjString*)next()))
    {
      tc->Add((fn->GetString()).Data());
      counter++;
    }
  
  input.close();
  Long_t iEntries = tc->GetEntries();
  
  cout << "Created chain of " << counter << " files for a total of " <<  iEntries << " entries" << endl;
  
  // TString outputname;
  TFile *coinc = new TFile(outputname+=".root", "recreate");
  
  RAT::DU::Utility::Get()->LoadDBAndBeginRun();
  RAT::DB* db = RAT::DB::Get();
  
  
  //  ------------ Histograms ------------------ 

  TH1D *hT = new TH1D("hT","Time difference", 1000,0.,10000);
  hT->SetXTitle("Time [#mus]");
  hT->SetYTitle("events for 10 us");

  TH1D *hDr = new TH1D("hDr","Time difference", 1000,0.,2000);
  hDr->SetXTitle("Position difference [mm]");
  hDr->SetYTitle("events every 2 mm");

  TH1D *hZ0 = new TH1D("hZ0","position tagged Z Po215", 1400,-7,7);
  hZ0->SetXTitle("Position [m]");
  hZ0->SetYTitle("events for 100 mm");

  TH1D *hZ = new TH1D("hZ","position tagged Z Rn219", 1400,-7,7);
  hZ->SetXTitle("Position [m]");
  hZ->SetYTitle("events for 100 mm");

  TH1D *hZS = new TH1D("hZS","position tagged Z for RnPo", 1400,-7,7);
  hZS->SetXTitle("Position [m]");
  hZS->SetYTitle("events for 100 mm");

  //--------

  TH1D *hnhits_1 = new TH1D("hnhits_1","nhits corrected for Po event", 6000,0.,6000);
  hnhits_1->SetXTitle("Nhits");
  hnhits_1->SetYTitle("events for 1Nhit");

  TH1D *hnhits_2 = new TH1D("hnhits_2","nhits corrected for Rn event", 6000,0.,6000);
  hnhits_2->SetXTitle("Nhits");
  hnhits_2->SetYTitle("events for 1Nhit");

  TH1D *Sum = new TH1D("Sum","nhits corrected for RnPo event", 6000,0.,6000);
  Sum->SetXTitle("Nhits");
  Sum->SetYTitle("events for 1Nhit");

  TH1D *hene_1 = new TH1D("hene_1","energy for Po event", 2000,0.,20);
  hene_1->SetXTitle("energy");
  hene_1->SetYTitle("events for 10 keV");

  TH1D *hene_2 = new TH1D("hene_2","energy clean for Rn event", 2000,0.,20);
  hene_2->SetXTitle("energy");
  hene_2->SetYTitle("events for 10 keV");

  TH1D *SumE = new TH1D("SumE","energy for RnPo", 2000,0.,20);
  SumE->SetXTitle("energy");
  SumE->SetYTitle("events for 10 keV");

  //------------

  TH1D *itr_rn = new TH1D("itr_rn","itr for Rn event", 2000,-1.,1);
  itr_rn->SetXTitle("itr");
  itr_rn->SetYTitle("events");

  TH1D *itr_po = new TH1D("itr_po","itr for Po event", 2000,-1.,1);
  itr_po->SetXTitle("itr");
  itr_po->SetYTitle("events");

  //------------

  TH2D *hXY0 = new TH2D("hXY0","X vs Y for Po",1800,-9,9, 1800, -9, 9);
  hXY0->SetXTitle("X [m]");
  hXY0->SetYTitle("Y [m]");

  TH2D *hXY = new TH2D("hXY","X vs Y for Rn",1800,-9,9, 1800, -9, 9);
  hXY->SetXTitle("X [m]");
  hXY->SetYTitle("Y [m]");

  TH2D *hZxyS0 = new TH2D("hZxyS0","Position Z vs #rho^2 for Po",8100,0,81, 1800, -9, 9);
  hZxyS0->SetXTitle("#rho^2 [m^2]");
  hZxyS0->SetYTitle("Pos z [m]");

  TH2D *hZxyS = new TH2D("hZxyS","Position Z vs rho^2 for Rn",8100,0,81, 1800, -9, 9);
  hZxyS->SetXTitle("#rho^2 [m^2]");
  hZxyS->SetYTitle("Pos z [m]");

  TH2D *hZxySSum = new TH2D("hZxySSum","Position Z vs #rho^2 for RnPo",8100,0,81, 1800, -9, 9);
  hZxySSum->SetXTitle("#rho^2 [m^2]");
  hZxySSum->SetYTitle("Pos z [m]");

  TH2D *hZxy0 = new TH2D("hZxy0","Position Z vs #rho for Po",900,0,9, 1800, -9, 9);
  hZxy0->SetXTitle("#rho [m]");
  hZxy0->SetYTitle("Pos z [m]");

  TH2D *hZxy = new TH2D("hZxy","Position Z vs #rho for Rn",900,0,9, 1800, -9, 9);
  hZxy->SetXTitle("rho [m]");
  hZxy->SetYTitle("Pos z [m]");

  TH2D *hZxySum = new TH2D("hZxySum","Position Z vs #rho for RnPo ",900,0,9, 1800, -9, 9);
  hZxySum->SetXTitle("rho [m]");
  hZxySum->SetYTitle("Pos z [m]");

  //-----------

  TH2D *hRhit0 = new TH2D("hRhit0","Nhits corrected as a function of radius for Po", 2000,0.,2, 2000, 0, 2000);
  hRhit0->SetXTitle("R^3/R^_AV");
  hRhit0->SetYTitle("Nhits");

  TH2D *hRhit = new TH2D("hRhit","Nhits corrected as a function of radius for Rn", 2000,0.,2, 2000, 0, 2000);
  hRhit->SetXTitle("R^3/R^_AV");
  hRhit->SetYTitle("Nhits");

  TH2D *hZN = new TH2D("hZN","Nhits corrected vs Z for Po",3000,0,6000, 1800, -9, 9);
  hZN->SetXTitle("N [nhits]");
  hZN->SetYTitle("Z [m]");

  ///////////////////------------ Parameters ----------------//////////////////////

  int events, nhits, nhitsCleaned, tag, k_e, eventID, ID0, ID;

  Bool_t fitValid, scintFit;

  Double_t energy, correctedNhits, nhits1, nhits0, itr0, itr1, itr, beta14, posx, posy, posz, posr, rho, rho0, posx0, posy0, posz0, Dr, ene0, ene1, r3, r30, ratio;

  ULong64_t dcFlagged, dcApplied, clockCount50, time0, time1, Dt;

  tag =0;
  k_e=0;

  //------------
  
  tc->SetBranchAddress("fitValid",&fitValid); 
  tc->SetBranchAddress("scintFit",&scintFit);

  tc->SetBranchAddress("nhits",&nhits);
  tc->SetBranchAddress("nhitsCleaned",&nhitsCleaned);
  tc->SetBranchAddress("correctedNhits",&correctedNhits);
  tc->SetBranchAddress("energy",&energy);

  tc->SetBranchAddress("beta14",&beta14);
  tc->SetBranchAddress("itr",&itr);

  tc->SetBranchAddress("posr",&posr);
  tc->SetBranchAddress("posx",&posx);
  tc->SetBranchAddress("posy",&posy);
  tc->SetBranchAddress("posz",&posz);

  tc->SetBranchAddress("dcFlagged",&dcFlagged);
  tc->SetBranchAddress("dcApplied",&dcApplied);
  tc->SetBranchAddress("clockCount50",&clockCount50);
  tc->SetBranchAddress("eventID",&eventID);


  // Loop through all of them

  for(Long_t i=0; i<iEntries;++i)
    {
      tc->GetEntry(i);

      //Correct for Av shift
      posz = posz - 186.4; 
      posr = sqrt(posx*posx+posy*posy+posz*posz);

      rho = sqrt(posx*posx + posy*posy)/1000;
      r3 = posr*posr*posr/(6005*6005*6005);
      ID = eventID;

      if (nhits>0) ratio = (float(nhits) - float(nhitsCleaned))/float(nhits);

      //common cuts

      if (((0x2100000042C2) & dcFlagged ) != (0x2100000042C2)) continue; //DC mask
      if ( scintFit != 1 ) continue;
      if ( ratio > 0.2 ) continue;
      if ( posr > FV ) continue;

      //Search first for Po215

      if ( energy < 0.6 || energy > 1.5 ) continue;

      //Po215 quantities have the index0
      time0 = clockCount50;
      nhits0 = correctedNhits;
      ene0 = energy;
      posx0 = posx;
      posy0 = posy;
      posz0 = posz; 
      rho0 = rho;
      itr0 = itr;
      ID0 = ID;
      r30=r3;

      bool pair = false; //no idetified coincidence to start

      //look back for the Rn219 event

      for (int k =1; k <= (i-k_e); k++)
	{
          if (pair) break; //if already a coicidence, exit
	  
	  tc->GetEntry(i-k);
	  
	  time1 = clockCount50;
	  Dt = ( ( time0 - time1 ) & 0x7FFFFFFFFFF )*20;
	  
	  if ( Dt > 10e6 ) break; //10 ms
	  if ( Dt < 4000 ) continue; //4 us
	  
	  //Correct for Av shift
	  posz = posz - 186.4;
	  posr = sqrt(posx*posx+posy*posy+posz*posz);

	  rho = sqrt(posx*posx + posy*posy)/1000;
	  r3 = posr*posr*posr/(6005*6005*6005);
	  ID = eventID;

	  if (((0x2100000042C2) & dcFlagged ) != (0x2100000042C2)) continue; //DC mask

	  if ( scintFit != 1 ) continue;
	  if ( posr > FV ) continue;
	  if ( energy < 0.6 || energy > 1.1 ) continue; 

	  nhits1 = correctedNhits;
	  ene1 = energy;
	  itr1 = itr;
	  Dr = sqrt( (posx-posx0)**2 + (posy-posy0)**2 + (posz-posz0)**2 );

	  if ( Dr > 800 ) continue;

	  hnhits_1->Fill( nhits0 );
	  hnhits_2->Fill( nhits1 );
	  hene_1->Fill( ene0 );
	  hene_2->Fill( ene1 );

	  hT->Fill(Dt/1e3); //in us
	  hDr->Fill(Dr);

	  hXY->Fill(posx/1000,posy/1000);
	  hXY0->Fill(posx0/1000,posy0/1000);

	  hZxy->Fill( rho, posz/1000);
	  hZxy0->Fill( rho0, posz0/1000);
	  hZxySum->Fill( rho, posz/1000);
	  hZxySum->Fill( rho0, posz0/1000);

	  hZxyS->Fill( (rho*rho), posz/1000);
	  hZxyS0->Fill( (rho0*rho0), posz0/1000);
	  hZxySSum->Fill( (rho*rho), posz/1000);
	  hZxySSum->Fill( (rho0*rho0), posz0/1000);

	  hRhit->Fill(r3, nhits1);
	  hRhit0->Fill(r30, nhits0);
	  
	  itr_po->Fill(itr0);
	  itr_rn->Fill(itr1);
	  hZN->Fill(nhits0,posz/1000);
	  
	  pair = true; //now that a pair is found set it true and exit
	  
	  k_e = k-1;
	  tag++;
	  
	}//Rn loop

  //reset
  time0 = time1 = 0;
  nhits0 = nhits1 = 0;
  ene0=ene1=0;
  
    }//entry loop

  Sum->Add(hnhits_2,hnhits_1,1,1);
  SumE->Add(hene_2,hene_1,1,1);

  cout<<"Total numbers of entries: "<<i<<" Total tagged events: "<<tag<<endl;
  
  coinc->Write();
  coinc->Close();
  
}

