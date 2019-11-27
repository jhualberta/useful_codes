/*
Rat6.5.2, track optical photon
*/
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/TrackNode.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNav.hh>
#include <iostream>
#include <vector>
#include "TH2D.h"
#include "TH1D.h"
#include "TH3D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
//NOTE: in Boulay thesis, he calculates total energy of e-, E = sqrt(0.511^2+Ek^2) MeV
using namespace::std;
void NavOpticalTrack()
{
    using namespace RAT;
    //using namespace std;
    RAT::DU::DSReader dsReader("FitRat_TrackOn_RopeOff_MC_partialscint_5MeVgamma_x0y0z-4000_level0_10evts.root");
    TH1F *hSx = new TH1F("hSx","Sx",2000,-10000,10000);
    TH1F *hSy = new TH1F("hSy","Sy",2000,-10000,10000);
    TH1F *hSz = new TH1F("hSz","Sz",2000,-10000,10000);
    TH1F *hTimeElectron = new TH1F("hTimeElectron","child e- global time",100,0,50);
    TH1F *hE_initialgamma = new TH1F("hE_initialgamma","E gamma at track start",1000,0,10);
    TH1F *hE_Compton = new TH1F("hE_Compton","compton energy",100,0,10); 
    TH1F *hE_PairProduct = new TH1F("hE_PairProduct","pair production",100,0,10);
    TH1F *hE_PEabsorb = new TH1F("hE_PEabsorb","PEabsorb",100,0,10);
    double Eelectron = 0.511; 
    TVector3 pos0; pos0.SetXYZ(0,0,-4000);
    
    int countBetaOnly = 0,countNone = 0, count2gamma = 0, count1gamma = 0, countGammaInside = 0;
    int countCompt = 0, countProd = 0, countPhot = 0;
    unsigned int evtNum = 8;//dsReader.GetEntryCount();

    for (size_t iEntry = 7; iEntry<evtNum; iEntry++)//dsReader.GetEntryCount();iEntry++)
    {
       const RAT::DS::Entry& ds = dsReader.GetEntry(iEntry);
       const RAT::DS::MC& rmc= ds.GetMC();
       int particleNum = rmc.GetMCParticleCount();
       std::cout<<"event "<<iEntry<<"  "<<particleNum<<std::endl;
       RAT::TrackNav nav1(&ds);
       RAT::TrackCursor c = nav1.Cursor(false);
       c.GoChild(0);
       RAT::TrackNode *n = c.Here();
       int i = 0;
       vector<size_t> trackID;
       while(1)
       {
         trackID.clear();
         trackID.push_back(n->GetTrackID());
         c.FindNextParticle("opticalphoton");
         n = c.Here(); TVector3 startPos = n->GetPosition();
         c.GoTrackEnd();
         n = c.Here(); double t0 = n->GetGlobalTime();
         trackID.push_back(n->GetTrackID());
         if(n->GetEndVolume() == "innerPMT_pmt") {
          cout<<endl;
          cout<<"optical photon hit innerPMT_PMT"<<endl;
          cout<<"t0 "<<t0<<endl;
          cout<<"sub track start "<<startPos.X()<<", "<<startPos.Y()<<", "<<startPos.Z()<<endl;
          TVector3 endPos(n->GetPosition().X(), n->GetPosition().Y(), n->GetPosition().Z());
          cout<<"sub track end "<<endPos.X()<<", "<<endPos.Y()<<", "<<endPos.Z()<<endl;
          /// calculate track length
          c.GoTrackStart();
          double trackLength = 0;
          while(!c.IsTrackEnd())
          {
           n = c.Here();
           trackLength += n->GetLength();
           c.GoNext();
          }
          double tof = trackLength/(299./1.38486); 
          cout<<"tof "<<tof<<endl;
          cout<<"tof straight "<<(endPos-startPos).Mag()/(299./1.38486)<<endl;
          }
         else continue;
         i++;
         /// when cursor cannot find more track 
         if(trackID[i] == trackID[i-1]) break;
      }
    }//ensure it is gamma11
}
