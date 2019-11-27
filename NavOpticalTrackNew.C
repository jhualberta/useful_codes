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
void NavOpticalTrackNew()
{
    double n_eff = 1.38486;
    
    //double n_eff = 1.59;
    
    using namespace RAT;
    //using namespace std;
    double wavelength[] = {200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0}; 

    double ri[] = {1.41615, 1.39727, 1.38395, 1.37414, 1.36667, 1.36082, 1.35615, 1.35233, 1.34916, 1.3465, 1.34423, 1.34227, 1.34058, 1.33909, 1.33778, 1.33661, 1.33557, 1.33463, 1.33378, 1.33301, 1.33231, 1.33167, 1.33108, 1.33054, 1.33004, 1.32957, 1.32914, 1.32874, 1.32836, 1.32801, 1.32768};
  
    double wavelength_gr[] = {200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0, 710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0, 750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0, 785.0, 790.0, 795.0, 800.0};

    double grVelocity[] = {182.709, 185.104, 187.309, 189.341, 191.217, 192.951, 194.558, 196.047, 197.432, 198.72, 199.92, 201.04, 202.088, 203.068, 203.986, 204.849, 205.659, 206.422, 207.14, 207.818, 208.458, 209.063, 209.636, 210.179, 210.694, 211.183, 211.648, 212.09, 212.51, 212.911, 213.294, 213.659, 214.008, 214.342, 214.661, 214.967, 215.26, 215.541, 215.81, 216.069, 216.318, 216.557, 216.787, 217.008, 217.222, 217.427, 217.626, 217.817, 218.001, 218.18, 218.352, 218.518, 218.679, 218.835, 218.986, 219.132, 219.274, 219.411, 219.544, 219.673, 219.798, 219.92, 220.038, 220.153, 220.264, 220.372, 220.478, 220.58, 220.68, 220.777, 220.872, 220.964, 221.054, 221.141, 221.226, 221.31, 221.391, 221.47, 221.547, 221.622, 221.696, 221.768, 221.838, 221.907, 221.974, 222.039, 222.104, 222.166, 222.228, 222.288, 222.346, 222.404, 222.46, 222.515, 222.569, 222.622, 222.674, 222.725, 222.775, 222.824, 222.871, 222.918, 222.964, 223.01, 223.054, 223.097, 223.14, 223.182, 223.223, 223.264, 223.303, 223.342, 223.381, 223.418, 223.455, 223.491, 223.527, 223.562, 223.597, 223.631, 223.664};


    RAT::DU::DSReader dsReader("FitRat655_mc_TrackOn_water_5MeVgamma_x0y0z-4000_10evts.root");
    TH1F *hSx = new TH1F("hSx","Sx",2000,-10000,10000);
    TH1F *hSy = new TH1F("hSy","Sy",2000,-10000,10000);
    TH1F *hSz = new TH1F("hSz","Sz",2000,-10000,10000);
    TH1F *hTimeElectron = new TH1F("hTimeElectron","child e- global time",100,0,50);
    TH1F *hE_initialgamma = new TH1F("hE_initialgamma","E gamma at track start",1000,0,10);
    TH1F *hE_Compton = new TH1F("hE_Compton","compton energy",100,0,10); 
    TH1F *hE_PairProduct = new TH1F("hE_PairProduct","pair production",100,0,10);
    TH1F *hE_PEabsorb = new TH1F("hE_PEabsorb","PEabsorb",100,0,10);
    double Eelectron = 0.511; 
    TVector3 pos0; pos0.SetXYZ(0,0,4000);
    
    int countBetaOnly = 0,countNone = 0, count2gamma = 0, count1gamma = 0, countGammaInside = 0;
    int countCompt = 0, countProd = 0, countPhot = 0;
    unsigned int evtNum = 1;//dsReader.GetEntryCount();

    for (size_t iEntry = 0; iEntry<evtNum; iEntry++)//dsReader.GetEntryCount();iEntry++)
    {
       const RAT::DS::Entry& ds = dsReader.GetEntry(iEntry);
       const RAT::DS::MC& rmc= ds.GetMC();
       int particleNum = rmc.GetMCParticleCount();
       std::cout<<"event "<<iEntry<<"  "<<particleNum<<std::endl;
       Int_t nevC = ds.GetEVCount();//!! if retriggered events, nevC == 2
       //if(nevC>1) nevC = 1;//!!! remove retrigger events
       //for(Int_t iev=0;iev<nevC; iev++) 
       {
         RAT::TrackNav nav1(&ds);
         RAT::TrackCursor c = nav1.Cursor(false);
         c.GoChild(0);
         RAT::TrackNode *n = c.Here();
         int i = 0;
         vector<size_t> trackID;
         trackID.clear();
         while(1) // searching tracks
         {
           double searchTrack = n->GetTrackID();
           trackID.push_back(searchTrack);
           c.FindNextParticle("opticalphoton");
           n = c.Here(); TVector3 startPos = n->GetPosition();
           searchTrack = n->GetTrackID();
           double t0 = n->GetGlobalTime();
           double E0 = n->GetKineticEnergy();
           while( !c.IsTrackEnd() )
           { 
             c.GoNext();
             n = c.Here();
             TString checkVolume = n->GetStartVolume();//n->GetEndVolume();
             TString checkVolume1 = n->GetEndVolume();
             double trackLength = 0;
             trackLength = trackLength + n->GetLength();
             if( checkVolume(0,15) == "innerPMT_pmtenv") 
             {
                double tf = n->GetGlobalTime();
                double Ef = n->GetKineticEnergy();
                double waveL = 6.626e-34*3e8/(Ef*1.6e-19*1e6)*1e9;
                TVector3 hitPMTpos = n->GetPosition();
                double cn = trackLength/(tf-t0);
                double cn_straight = (hitPMTpos-startPos).Mag()/(tf-t0);
                cout<<"sub track start "<<startPos.X()<<", "<<startPos.Y()<<", "<<startPos.Z()<<endl;
                cout<<"sub track end "<<hitPMTpos.X()<<", "<<hitPMTpos.Y()<<", "<<hitPMTpos.Z()<<endl;
                cout<<"geometry at the PMT: "<<n->GetStartVolume()<<endl;
                cout<<"t0 "<<t0<<" tf "<<tf<<" tf-t0 "<<tf-t0<<" wavelength0 "<<waveL<<endl;
                cout<<" cn = trackLength/(tf-t0)= "<<cn<<" mm/ns,  n_eff= "<<299./(cn)<<endl;
                cout<<" cn = posDiff/(tf-t0)= "<<cn_straight<<" mm/ns,  n_eff= "<<299./(cn_straight)<<endl;
                //cout<<"tof(fixed to 1.384 or 1.59), use trackLength "<<tof<<endl;
                cout<<"tof(fixed to 1.384 or 1.59), use posDiff "<<(hitPMTpos-startPos).Mag()/(299./n_eff)<<endl;
                int index = (waveL-200)/5;
                int index1 = (waveL-200)/20;
                cout<<"look up grVelocity table: "<<wavelength_gr[index]<<"nm, grVelocity = "<<grVelocity[index]<<" mm/ns"<<endl;
                cout<<"look up refraction index table: "<<wavelength[index1]<<"nm, n = "<<ri[index1]<<endl;
                cout<<endl;
                break;
             }
           }
           i++;
           /// when cursor cannot find more track 
           if(searchTrack == trackID[i-1]) break;
         }
      } //not retriggered event
    }//pick one event
}
