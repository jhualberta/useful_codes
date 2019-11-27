#include <cstdio>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TVector3.h"

#include "sno_style.H"
using namespace::std;

TCanvas* gCanvas = NULL;
// TVector3 srcPos(-1120.8,1041.4,6172.5);// For N16 source
TVector3 srcPos(-1120.8,1041.4,6108.0-108.0); 

double fvCut = 5900;
UInt_t fecdN16 = 9188;//9188;
double tooEarly = -20; // tRes early time cut
double tooLate = 80; //tRes late time cut

void CalcuN16() {
  SetupSNOStyle();
  SetGradientPalette();
  const char* filename = "Merged_Tree_N16_252308to252314.root";//Tree_FitPartial_SNOP_0000252308.root";
  TFile *fff = new TFile(filename,"read");
  // Make plots square
  Double_t w = 800;
  Double_t h = 600;
  // gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
   
  // Check that a file is attached
  if (!gFile) {
    std::cout << "File must be attached" << std::endl;
    return;
  }
  
  // Variables in the gOutput
  double gNhitClean = 0;
  double cosTheta = 0;
  double itr = 0;
  double distSrcToFitPos = 0;

  TVector3* gFitPosition = new TVector3();
 
  std::vector<double>* gPMTRecoTime = new std::vector<double>;
  std::vector<double>* gPMTX = new std::vector<double>;
  std::vector<double>* gPMTY = new std::vector<double>;
  std::vector<double>* gPMTZ = new std::vector<double>;

  std::vector<double>* gQhs = new std::vector<double>;
  std::vector<double>* gQhl = new std::vector<double>;
  std::vector<UInt_t>* gFECD  = new std::vector<UInt_t>;

  // Access the gOutput
  TTree* gOutput = (TTree*) gFile->Get("gOutput");
  if (!gOutput) {
    std::cout << "The selected event gOutput must be present" << std::endl;
    return;
  }
 
  // Access Branches in the Tree
  gOutput->SetBranchAddress("gFitPosition",&gFitPosition);
  /// double values
  gOutput->SetBranchAddress("cosTheta", &cosTheta);
  //gOutput->SetBranchAddress("gFitTime",&gFitTime);
  //gOutput->SetBranchAddress("gFecdTime",&gFecdTime);
  gOutput->SetBranchAddress("gNhitClean",&gNhitClean);
  gOutput->SetBranchAddress("itr",&itr);
  gOutput->SetBranchAddress("distSrcToFitPos", &distSrcToFitPos);
  /// vector containers 
  // gOutput->SetBranchAddress("gPMTTrueTime",&gPMTTrueTime);
  gOutput->SetBranchAddress("gPMTRecoTime",&gPMTRecoTime);
  gOutput->SetBranchAddress("gQhs",&gQhs);
  gOutput->SetBranchAddress("gQhl",&gQhl);

  gOutput->SetBranchAddress("gFECD",&gFECD);
  gOutput->SetBranchAddress("gPMTX",&gPMTX);
  gOutput->SetBranchAddress("gPMTY",&gPMTY);
  gOutput->SetBranchAddress("gPMTZ",&gPMTZ);
  gOutput->SetBranchAddress("distSrcToFitPos", &distSrcToFitPos);

  TString newfile = "extract_"+ TString(filename);
  TFile *outputFile = new TFile(newfile, "RECREATE");

  /// here save the essential variables we are interested !!
  double cosPMT = 0, cosPMTwQhs = 0, cosPMTwQhl = 0;
  UInt_t nhitsC = 0;
  double distCut = 0;
  double tRes = 0;
  double qhs = 0;
  double qhl = 0;
  double posx = 0, posy = 0, posz = 0;
  double pmtx = 0, pmty = 0, pmtz = 0;
  TTree *Tdump = new TTree("Tdump","dump PMT tree data");
  Tdump->Branch("cosPMT", &cosPMT);
//Tdump->Branch("cosPMTwQhs", &cosPMTwQhs);
//Tdump->Branch("cosPMTwQhl", &cosPMTwQhl);

  Tdump->Branch("nhitsC", &nhitsC);
  Tdump->Branch("distCut", &distCut);
  Tdump->Branch("tRes", &tRes);
  Tdump->Branch("qhs", &qhs);
  Tdump->Branch("qhl", &qhl);
  Tdump->Branch("posx", &posx);
  Tdump->Branch("posy", &posy);
  Tdump->Branch("posz", &posz); // Z is corrected

  Tdump->Branch("pmtx", &pmtx);
  Tdump->Branch("pmty", &pmty);
  Tdump->Branch("pmtz", &pmtz); // Z is corrected

  //TFile* outputFile = NULL;
  TString fileName = gFile->GetName();

  double nhitsMax = 2000;
  
  double startTime = -5.;
  double divisionTime = 1.;
  double endTime = 20.;
  
  // Histograms of interest to fill
  TH1F* hCosAngle = new TH1F("hCosAngle", "(#vec{X}_{evt}-#vec{X}_{source})#cdot(#vec{X}_{PMT}-#vec{X}_{evt});Cos(#theta_{PMT});Entries", 200, -1.0, 1.0);

  TH2F* hcosThetaVsTres = new TH2F("hcosThetaVsTres", "Angle between Directions;Cos(#theta_{PMT}); time residual, fecd == 9188, nhitCut", 200, -1.0, 1.0, 400, -100, 300);

  TH2F* hcosThetaVsTres_wQhs = new TH2F("hcosThetaVsTres_wQhs", "Angle between Directions;Cos(#theta_{PMT}); time residual, fecd == 9188, nhitCut", 200, -1.0, 1.0, 400, -100, 300);

  TH2F* hcosThetaVsTres_wQhl = new TH2F("hcosThetaVsTres_wQhl", "Angle between Directions;Cos(#theta_{PMT}); time residual, fecd == 9188, nhitCut", 200, -1.0, 1.0, 400, -100, 300);

  TH2F* hfitRhoZ = new TH2F("hfitRhoZ", "SNO+  fitted Rho vs Z, FECD cut", 1000, 0,9000,2000,-9000,9000);
  TH2F* hfitXZ = new TH2F("hfitXZ", "SNO+  fitted X vs Z, FECD cut", 2000, -9000,9000,2000,-9000,9000);

  gOutput->Project("hCosAngleROI","cosTheta","gFECD == 9188 && gNhitClean>1000 && distSrcToFitPos>1500");
  gOutput->Project("hfitRhoZ","gFitPosition->Z():gFitPosition->Perp()","gFECD == 9188 && gNhitClean>1200 && distSrcToFitPos>1000");
//  hCosAngleROI->Draw();

  std::cout << "Number of Entries: " << gOutput->GetEntries() << std::endl;
//  
//  // Loop through the gOutput
  for (int iev = 0; iev < gOutput->GetEntries(); iev++) {
    if(iev%5000 == 0) cout<<"processed "<<iev<<" events"<<endl;
    gOutput->GetEntry(iev);
    // FECD cut
    vector<UInt_t>::iterator it;
    it = find( gFECD->begin(), gFECD->end(), fecdN16 );

    vector<UInt_t>::iterator it2;
    it2 = find( gFECD->begin(), gFECD->end(), 9207 );

    if ( it == gFECD->end() && it != gFECD->end() ) { //cout<<"not an N16 tagged event! "<<endl; 
      continue;}
//    // Get angle between true direction and the solar direction
    // Skip events that occur outside the Fiducial volume
    TVector3 pos_cor(gFitPosition->X(),gFitPosition->Y(),gFitPosition->Z()-108); //correct z to AV coordinator
    if( pos_cor.Mag() > fvCut ) continue;
    if ( gPMTRecoTime->size() != gPMTX->size() ) { cout<<"data dimension error!"; continue;}
    // Loop PMTs

    for (unsigned int i = 0; i < gPMTX->size(); i++) {
     // if ( gPMTRecoTime->at(i) < -20 || gPMTRecoTime->at(i) > 100) continue;
     TVector3 PMTPos;
     PMTPos.SetXYZ(gPMTX->at(i), gPMTY->at(i), gPMTZ->at(i));
     TVector3 PosToPMT;
     PosToPMT.SetXYZ(PMTPos.X()-pos_cor.X(), PMTPos.Y()-pos_cor.Y(), PMTPos.Z()-pos_cor.Z() );
     TVector3 SrcToPos;
     SrcToPos.SetXYZ(pos_cor.X() -srcPos.X(), pos_cor.Y()-srcPos.Y(), pos_cor.Z()-srcPos.Z());
     double cosTheta = PosToPMT.Unit()*SrcToPos.Unit(); 
     hcosThetaVsTres->Fill(cosTheta, gPMTRecoTime->at(i));
     if(qhs>0) hcosThetaVsTres_wQhs->Fill(cosTheta, gPMTRecoTime->at(i),qhs);
     if(qhl>0) hcosThetaVsTres_wQhl->Fill(cosTheta, gPMTRecoTime->at(i),qhl);

     /// save into new tree
     cosPMT = cosTheta;
     distCut = distSrcToFitPos;
     nhitsC = gNhitClean;
     tRes = gPMTRecoTime->at(i);
     qhs = gQhs->at(i);
     qhl = gQhl->at(i);
     posx = pos_cor.X(); posy = pos_cor.Y(); posz = pos_cor.Z();
     pmtx = PMTPos.X(); pmty = PMTPos.Y(); pmtz = PMTPos.Z();
     Tdump->Fill();
    }
   }
  // hcosThetaVsTres->Draw("colz");  
  //hcosThetaVsTres->ProjectionX("hpx",-5+100,1+100)->Draw(); 
  outputFile->cd();
  hcosThetaVsTres->Write();
  hcosThetaVsTres_wQhs->Write();
  hcosThetaVsTres_wQhl->Write();
  hfitRhoZ->Write();
//  Tdump->Write();
  outputFile->Write();
  outputFile->Close();

}
    
//    // Check the size of all vectors are the same
//    if (gPMTTrueTime->size() != gPMTRecoTime->size() || gPMTTrueTime->size() != gPMTX->size()) {
//      std::cout << "Issue with vector sizes! Entry: " << e << std::endl;
//      std::cout << "Aborting..." << std::endl;
//      return;
//    }
//    
//    double truehighcos = 0.;
//    double trueearlytime = 0.;
//    double truealltime = 0.;
//    
//    double recohighcos = 0.;
//    double recoearlytime = 0.;
//    double recoalltime = 0.;
//    
//    for (unsigned int i = 0; i < gPMTTrueTime->size(); i++) {
//      // Skip the events outside my histogram
//      if (gPMTTrueTime->at(i) < -20 || gPMTTrueTime->at(i) > 80 ||
//          gPMTRecoTime->at(i) < -20 || gPMTRecoTime->at(i) > 80) continue;
//      TVector3 PMTDir;
//      // Get the position of the hit PMT
//      TVector3 PMTPos;
//      // Calculate the cos between Solar direction and PMT position.
//      double cosPMT = TMath::Cos(PosToPMT.Angle(-1*(*gSolarDirection)));
//      if (gTruenhits > ROIMin && gTruenhits < ROIMax) {
//        hPMTTrueTimeROI->Fill(gPMTTrueTime->at(i), cosPMT);
//        hPMTRecoTimeROI->Fill(gPMTRecoTime->at(i), cosPMT);
//      }
//      hPMTTrueTimevsnhits->Fill(gPMTTrueTime->at(i), cosPMT, gTruenhits);
//      hPMTRecoTimevsnhits->Fill(gPMTRecoTime->at(i), cosPMT, gTruenhits);
//
//      // Calculate Classifiers
//      // Early Time
//      if (gPMTTrueTime->at(i) < divisionTime && gPMTTrueTime->at(i) > startTime) {
//        trueearlytime++;
//        if (cosPMT > cosPMTCut) truehighcos++;
//      }
//      // Total time
//      if (gPMTTrueTime->at(i) > startTime && gPMTTrueTime->at(i) < endTime) truealltime++;
//      
//      // Calculate Classifiers
//      // Early Time
//      if (gPMTRecoTime->at(i) < divisionTime && gPMTRecoTime->at(i) > startTime) {
//        recoearlytime++;
//        if (cosPMT > cosPMTCut) recohighcos++;
//      }
//      // Total time
//      if (gPMTRecoTime->at(i) > startTime && gPMTRecoTime->at(i) < endTime) recoalltime++;
//    } // End of loop through hit PMTs
//    
//    hTrueEarlyLightvsnhits->Fill(trueearlytime/truealltime, gTruenhits);
//    hTrueHighCosvsnhits->Fill(truehighcos/trueearlytime, gTruenhits);
//    hTrueCorrelationvsnhits->Fill(truehighcos/trueearlytime, trueearlytime/truealltime, gTruenhits);
//    
//    hRecoEarlyLightvsnhits->Fill(recoearlytime/recoalltime, gTruenhits);
//    hRecoHighCosvsnhits->Fill(recohighcos/recoearlytime, gTruenhits);
//    hRecoCorrelationvsnhits->Fill(recohighcos/recoearlytime, recoearlytime/recoalltime, gTruenhits);
//    
//    if (gTruenhits > ROIMin && gTruenhits < ROIMax) {
//      hTrueEarlyLightROI->Fill(trueearlytime/truealltime);
//      hTrueHighCosROI->Fill(truehighcos/trueearlytime);
//      hTrueCorrelationROI->Fill(truehighcos/trueearlytime, trueearlytime/truealltime);
//
//      hRecoEarlyLightROI->Fill(recoearlytime/recoalltime);
//      hRecoHighCosROI->Fill(recohighcos/recoearlytime);
//      hRecoCorrelationROI->Fill(recohighcos/recoearlytime, recoearlytime/recoalltime);
//    }
//    
//  } // End of loop through gOutput
//  
//  outputFile->cd();
//  
//  outputFile->Write();
//  outputFile->Close();
//}
