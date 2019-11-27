#include <vector>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TF1.h"
#include "TEllipse.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TPluginManager.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "Math/GoFTest.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DataCleaningUtility.hh>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TNtuple.h>

#include <bitset>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;
using namespace TMVA;
//Get pmtInfo
//const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();

void NeutronSearch (const RAT::DS::EV& fev, const RAT::DU::PMTInfo& pmtInfo, Float_t tstart, Float_t tend, TVector3 pos_mc); 
Int_t GetNhits(Float_t *v, Int_t start_index, Float_t width, Int_t nhits);
Int_t GetCluster (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
		                            Int_t N12, Int_t &ncut, Int_t ncth, Float_t thr); 
Float_t GetWeight (const TVector3 xyz, const Float_t v[3]) ;
Float_t EffCos ( Float_t costh ) ;
Int_t GetLowHits (Int_t *ci, Int_t N12, Float_t *wt, Float_t acceptance) ;
Float_t GetWeightThreshold (const Float_t *w, const Float_t frac) ;
Bool_t CheckBackHits (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
		                            Int_t N12, Int_t &ncut, Float_t angle);
Float_t GetTrms (Float_t *ti, Int_t *flag, Int_t N12);
Float_t GetPhiRms (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N12);
Float_t GetThetaMean (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N12);
Int_t GetNXX(Int_t nhits, Float_t *t, Float_t twin, Float_t tcenter);
Double_t GetTriggeredEventTime(const RAT::DS::EV& tev);

const Float_t C_WATER = 21.6628 ; //speed of light in water
const Int_t   MAXNP   = 100;     // Maximum N12 peaks in 420 ns data
const Int_t   MAXN12  = 200; //Maximum N12
const Int_t   MAXCOMB = 10000;
double timeCutValue = 1000.; //1000 us 
//string fitName = "multipath" ;
string fitName = "waterFitter" ;
Int_t N12TH=4 ;
TTree *theOTree; 

typedef struct
{
	Int_t    np; 
        Int_t    N80M; 
        Float_t  T80M; 
        Int_t    N12[MAXNP];
	Int_t    Nc[MAXNP];
	Int_t    Nlow[MAXNP];
	Int_t    Nback[MAXNP];
	Int_t    N100[MAXNP];
	//store event pair
	Int_t    PrunID[MAXNP]; //for primary event
	Int_t    PsubrunID[MAXNP];
	Int_t    PevtID[MAXNP];
	Int_t    PNhits[MAXNP]; 
	Int_t    DNhits[MAXNP]; 
	Int_t    DevtID[MAXNP]; //for delayed event. assume share the same run and subrun, negligible missing rate
	Float_t    trms[MAXNP];
	Float_t    phirms[MAXNP];
	Float_t    thetam[MAXNP];
	Float_t    distan[MAXNP];
	Float_t    deltaT[MAXNP]; //in us
        Float_t    posX[MAXNP];
        Float_t    posY[MAXNP];
        Float_t    posZ[MAXNP];
        Float_t    energy[MAXNP];
        Float_t    logL[MAXNP];
        Float_t    selectedNhits[MAXNP];
        Float_t    logLscale[MAXNP];
   
	//experimental
	Float_t  dalpha[MAXCOMB];
	Int_t    ncomb; 
        void Clear() 
        {
            np   = 0;
            N80M = 0;
            T80M = -9999.;
            for (Int_t i=0; i<MAXNP; i++) {
                N12[i]   = 0;
                Nc[i]   = 0;
                Nlow[i]   = 0;
                Nback[i]   = 0;
                N100[i]   = 0;
		PrunID[i]   = 0;
		PNhits[i]   = 0;
		DNhits[i]   = 0;
		PsubrunID[i]   = 0;
		PevtID[i]   = 0;
		DevtID[i]   = 0;
                trms[i]   = 0;
                phirms[i]   = 0;
                thetam[i]   = 0;
                distan[i]   = 0;
                deltaT[i]   = 0;
                posX[i]     = 0;
                posY[i]     = 0;
                posZ[i]     = 0;
                logL[i]     = 0;
                selectedNhits[i] = 0;
                logLscale[i] = 0;
	    }
            // Experimental
	    ncomb = 0;
            for (Int_t i=0; i<MAXCOMB; i++) dalpha[i]=0;
         }
}results_t;
results_t res; 


//int main(int argc, char** argv)
void Ambe()
{
//  if ( argc != 3 ) {
//	          cout << "Usage: " << argv[0] << " OUTPUT " << 
//			             " INPUT " << endl;
//	         return -1;
//  }
  // Read in file

//P Prompt 
//D Delay
  string Infname = "109133_s000_p000.root";
  TString OUTfname;
//  Infname = argv[2];
//  OUTfname = argv[1];
  RAT::DU::DSReader dsreader(Infname);
  // To get the run info
  const RAT::DS::Run& run = dsreader.GetRun();
   //Get pmtInfo
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TString savefile = "SaveAmBe_"+TString(Infname);
   TFile *fout = new TFile(savefile, "RECREATE");
   theOTree = new TTree("sno2p2", "SNOP 2.2 MeV");
   theOTree->Branch("np",    &res.np,     "np/I");
   theOTree->Branch("N12",    res.N12,    "N12[np]/I");
   theOTree->Branch("ncomb", &res.ncomb,  "ncomb/I");
   theOTree->Branch("dalpha", res.dalpha, "dalpha[ncomb]/F");
   theOTree->Branch("Nc",    res.Nc,    "Nc[np]/I");
   theOTree->Branch("Nlow",    res.Nlow,    "Nlow[np]/I");
   theOTree->Branch("Nback",    res.Nback,    "Nback[np]/I");
   theOTree->Branch("N100",   res.N100,   "N100[np]/I");
   theOTree->Branch("PrunID",   res.PrunID,   "PrunID[np]/I");
   theOTree->Branch("PsubrunID",   res.PsubrunID,   "PsubrunID[np]/I");
   theOTree->Branch("PevtID",   res.PevtID,   "PevtID[np]/I");
   theOTree->Branch("PNhits",   res.PNhits,   "PNhits[np]/I");
   theOTree->Branch("DNhits",   res.DNhits,   "DNhits[np]/I");
   theOTree->Branch("DevtID",   res.DevtID,   "DevtID[np]/I");
   theOTree->Branch("trms",    res.trms,    "trms[np]/F");
   theOTree->Branch("phirms",    res.phirms,    "phirms[np]/F");
   theOTree->Branch("thetam",    res.thetam,    "thetam[np]/F");
   theOTree->Branch("distan",    res.distan,    "distan[np]/F");
   theOTree->Branch("deltaT",    res.deltaT,    "deltaT[np]/F");

   theOTree->Branch("posX",    res.posX,    "posX[np]/F");
   theOTree->Branch("posY",    res.posY,    "posY[np]/F");
   theOTree->Branch("posZ",    res.posZ,    "posZ[np]/F");
   theOTree->Branch("energy",  res.energy,"energy[np]/F");
   theOTree->Branch("logL",    res.logL,    "logL[np]/F");
   theOTree->Branch("selectedNhits", res.selectedNhits, "PositionSelectedNHit[np]/F");
   theOTree->Branch("logLscale",    res.logLscale,    "logLscale[np]/F");
   //cout << "SKY is great************" << endl;
  
  // Loop through events
  for(size_t i=0; i<dsreader.GetEntryCount();i++){
    if(i%10000==0) cout << "Dealing with " << i << "th event" << endl;
    //if(i>50000) break; //use small sample for test
    //if(i>500) cout << "Dealing with " << i << "th event" << endl;
    const RAT::DS::Entry& rds = dsreader.GetEntry(i);
    Int_t runID = run.GetRunID();
    Int_t subrunID = run.GetSubRunID();
    // get the position from fixed postion
    /*double vx, vy, vz; 
    if(runID==109133||runID==109134||runID==109135){vx=0.; vy=0.;vz=0.;} //mm 
    else if(runID==109137) {vx=0.; vy=440.*10;vz=0.;} 
    else if(runID==109140) {vx=0.; vy=300.*10;vz=0.;} 
    else if(runID==109144) {vx=0.; vy=150.*10;vz=0.;} 
    else if(runID==109147) {vx=0.; vy=-150.*10;vz=0.;} 
    else if(runID==109150) {vx=0.; vy=-300.*10;vz=0.;} 
    else if(runID==109153) {vx=0.; vy=-400.*10;vz=0.;} 
    else if(runID==109156) {vx=0.; vy=-260.*10;vz=-260.*10;} 
    else if(runID==109159) {vx=0.; vy=260.*10;vz=-260.*10;} 
    else if(runID==109162) {vx=0.; vy=0.;vz=-550.*10;} 
    else if(runID==109165) {vx=0.; vy=0.;vz=-500.*10;} 
    else if(runID==109168) {vx=0.; vy=0.;vz=-300.*10;} 
    else if(runID==109171) {vx=0.; vy=0.;vz=-450.*10;} 
    else if(runID==109174) {vx=0.; vy=0.;vz=300.*10;} 
    else if(runID==109178) {vx=0.; vy=0.;vz=450.*10;} 
    else if(runID==109181) {vx=0.; vy=0.;vz=550.*10;} 
    else if(runID==109208) {vx=0.; vy=0.;vz=150.*10;} 
    else if(runID==109211) {vx=0.; vy=0.;vz=-150.*10;} 
    else if(runID==109214) {vx=0.; vy=150.*10.;vz=-150.*10;} 
    else if(runID==109217) {vx=0.; vy=-150.*10.;vz=-150.*10;} 
    else if(runID==109220) {vx=0.; vy=-150.*10.;vz=150.*10;} 
    else if(runID==109223) {vx=0.; vy=150.*10.;vz=150.*10;}  
    else if(runID==109226) {vx=0.; vy=0.;vz=500.*10;} 
    else {vx=0.; vy=0.;vz=0.;}  //set a dummy value 
    
    TVector3 pos_mc(vx, vy, vz); 
    */
    TVector3 pos_mc;
    // triggered event
    int nevC = rds.GetEVCount();
    for(int iev=0;iev<nevC; iev++){
      // Get the Event information
      const RAT::DS::EV& rev = rds.GetEV(iev);
      int PevtID = rev.GetGTID(); //trigger ID
      //if (iev>0) continue; //data dont know if re-trigger or not
      int nHits = rev.GetNhits();
      int PNhits = nHits; //for store 
      //First, search for good prompt event
      double TimePri = GetTriggeredEventTime(rev); //ns 
      bool foundPrimary = false;
      //if(nHits>14 && nHits<50.) foundPrimary = true; //define a primary event
      //lower down the threshold to see the distribution
      if(fitName == "waterFitter"){
           if(nHits>14 && nHits<50.) foundPrimary = true; //define a primary event
	   else continue;
           if( !rev.FitResultExists(fitName)) continue;
           if( !rev.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) continue;
           if( !rev.GetFitResult(fitName).GetVertex(0).ValidPosition() ) continue;
           //if( !rev.GetFitResult(fitName).GetVertex(0).ContainsTime() ) continue;
           //if( !rev.GetFitResult(fitName).GetVertex(0).ValidTime() ) continue;
      }
      else{
           if(nHits>4 && nHits<50.) foundPrimary = true; //define a primary event
           else continue; 
	   if( !rev.FitResultExists(fitName)) continue;
	   if( !rev.GetFitResult(fitName).GetValid()) continue;
      }
      RAT::DS::FitResult fResult = rev.GetFitResult(fitName);
      RAT::DS::FitVertex fVertex = fResult.GetVertex(0); 
      pos_mc = fVertex.GetPosition(); //vertex of primary event, reconstructed case
      double fom = rev.GetFitResult(fitName).GetFOM("PositionLogL");
      double selectedNhits = rev.GetFitResult(fitName).GetFOM("PositionSelectedNHit");
      double nhits = rev.GetNhits();
      double fitenergy = fVertex.GetEnergy();
      //cout << "Good primary event is found KO..." << endl;
      ///////////end for good prompt event
      if(foundPrimary) { //second, search for neutron events
         // cout << "Good primary event is found ..." << endl;
	  bool Beyondtime = false; //see whether still in 1000 us time window
	  for(int kd=i; kd<dsreader.GetEntryCount();kd++){ //search from this event
            if(Beyondtime) break;
	    const RAT::DS::Entry& rdelay = dsreader.GetEntry(kd);
	    for(int ievd=0;ievd<rdelay.GetEVCount(); ievd++){
		  const RAT::DS::EV& revdelay = rdelay.GetEV(ievd);
                  int DevtID = revdelay.GetGTID(); //trigger ID
		  double TimeDel = GetTriggeredEventTime(revdelay);
		  double deltaT = (TimeDel-TimePri)/1000.; //ns to us
		  //cout << "deltaT is " << deltaT << endl;
		  bool foundDelay = false;
		  if(deltaT<=1.e-5) continue; 
		  else if (deltaT>1.e-5 && deltaT<timeCutValue) foundDelay = true; 
		  else { Beyondtime = true; break; } //beyond time
		  if (foundDelay) { //do neutron search if within 1000 us
                       //cout << "Good delayed event is found ..." << endl;
		       nHits = revdelay.GetNhits();
		       int DNhits = nHits;
                       //if(nHits<5 || nHits>20.) continue; //if 25, strange events at high end for non-cleaned events
                       if(nHits<4 || nHits>25.) continue; //if 25, strange events at high end
                       if(fitName == "waterFitter"){
                           if( !revdelay.FitResultExists(fitName)) continue;
                           if( !revdelay.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) continue;
                           if( !revdelay.GetFitResult(fitName).GetVertex(0).ValidPosition() ) continue;
                           //if( !revdelay.GetFitResult(fitName).GetVertex(0).ContainsTime() ) continue;
                           //if( !revdelay.GetFitResult(fitName).GetVertex(0).ValidTime() ) continue;
                           if( !revdelay.GetFitResult(fitName).GetVertex(0).ValidEnergy() ) continue;
                        }
                       else{
                           // make sure the position exist
                           if( !revdelay.FitResultExists(fitName)) continue;
                           if( !revdelay.GetFitResult(fitName).GetValid()) continue; 
                       }
                       //cout << "Good delayed event is found KO..." << endl;
                       // Clear results of previous event
                       res.Clear();
                       //search
                       NeutronSearch(revdelay, pmtInfo, 0., 420., pos_mc);
		       //test for -100, 520. This will incease the bkg level but no increase for signal.
                       //NeutronSearch(revdelay, pmtInfo, -100., 520., pos_mc);
                       ////////// Fill the output tree 
                       //head info
                       for(Int_t npk=0; npk<res.np;npk++){ //only store with np evt
	                  res.PrunID[npk]  = runID ;
	                  res.PNhits[npk]  = PNhits ;
	                  res.DNhits[npk]  = DNhits ;
	                  res.PsubrunID[npk]  = subrunID ;
	                  res.PevtID[npk]  = PevtID ;
	                  res.DevtID[npk]  = DevtID ;
	                  res.deltaT[npk]  = deltaT ;
                          res.posX[npk]    = pos_mc.X();
                          res.posY[npk]    = pos_mc.Y();
                          res.posZ[npk]    = pos_mc.Z();
                          res.logL[npk]    = fom;
                          res.selectedNhits[npk] = selectedNhits;
                          res.logLscale[npk] = fom/selectedNhits;
                          res.energy[npk]  = fitenergy;
                       }
                       theOTree->Fill();
                       ///////////////////fill finish
		  }
	    }
        }
      }//find primary event
    } //triggered event loop
  } //entry loop
  //save result
  fout->Write();
  delete theOTree;
  fout->Close();
  delete fout;
  theOTree = 0;
  fout = 0; 
  cout << "Finished successfully" << endl;
}

void NeutronSearch (const RAT::DS::EV& fev, const RAT::DU::PMTInfo& pmtInfo, Float_t tstart, Float_t tend, TVector3 pos_mc) {
        //cout << "SKYY is great************" << endl;
        RAT::DS::FitResult fResult = fev.GetFitResult(fitName);
        RAT::DS::FitVertex fVertex = fResult.GetVertex(0); 
        TVector3 pos_fit = fVertex.GetPosition();
        TVector3 pos_del = pos_fit - pos_mc ;

	const RAT::DS::CalPMTs& calpmts = fev.GetCalPMTs();
	const Int_t MAXPMT = 9438;
	Float_t wt[MAXPMT]; 

	const Int_t MAXHITS = 10000;
	//cout << "id number is " << calpmts.GetPMT(0).GetID() << endl;
	TVector3 postest = pmtInfo.GetPosition(calpmts.GetPMT(0).GetID());
	//cout << "X position is " << postest.X() << endl;
        Int_t   nhits;
	Int_t   cabiz[MAXHITS], cabiz2[MAXHITS];
	Float_t tiskz[MAXHITS], tiskz2[MAXHITS];
	Float_t qiskz[MAXHITS], qiskz2[MAXHITS];
	TVector3 pos[MAXHITS];
        TVector3 pos2[MAXHITS];
	Int_t   index[MAXHITS];
	Float_t hitv_x[MAXHITS]; 
	Float_t hitv_y[MAXHITS];
	Float_t hitv_z[MAXHITS];
	nhits = 0;
	Int_t i;
	for (i=0; i<calpmts.GetCount(); i++) { 
	    cabiz2[nhits]= calpmts.GetPMT(i).GetID();
	    tiskz2[nhits] = calpmts.GetPMT(i).GetTime() ;
	    pos2[nhits] = pmtInfo.GetPosition(calpmts.GetPMT(i).GetID());
	    //if(i==1) cout << "PMT position is " << pos2[i].X() << endl; 
	    qiskz2[nhits] = calpmts.GetPMT(i).GetQHS();
	    nhits++;
        }
	 
	// TOF
	for (i=0; i<nhits; i++) {
		Float_t tof;
		tof = TMath::Sqrt((pos2[i].X()-pos_mc.X())*(pos2[i].X()-pos_mc.X())
				+(pos2[i].Y()-pos_mc.Y())*(pos2[i].Y()-pos_mc.Y())
				+(pos2[i].Z()-pos_mc.Z())*(pos2[i].Z()-pos_mc.Z()))/C_WATER/10.; //mm to cm
		tiskz2[i] -= tof;
        }
        // Sort hits by TOF-corrected time
	TMath::Sort(nhits, tiskz2, index, kFALSE); // In increasing order
	for (i=0; i<nhits; i++){
		cabiz[i] = cabiz2[ index[i] ];
		tiskz[i] = tiskz2[ index[i] ];
		qiskz[i] = qiskz2[ index[i] ];
		pos[i] = pos2[ index[i] ];
		//cout << "tisk is " << tiskz[i] << endl; 
	}
	// Calculate hit vectors
	for (i=0; i<nhits; i++) {
		Float_t pmt_r;
		pmt_r = TMath::Sqrt((pos[i].X()-pos_mc.X())*(pos[i].X()-pos_mc.X())
			     +(pos[i].Y()-pos_mc.Y())*(pos[i].Y()-pos_mc.Y())
	                     +(pos[i].Z()-pos_mc.Z())*(pos[i].Z()-pos_mc.Z())) ;
		hitv_x[i] = (pos[i].X() - pos_mc.X())/pmt_r;
		hitv_y[i] = (pos[i].Y() - pos_mc.Y())/pmt_r;
		hitv_z[i] = (pos[i].Z() - pos_mc.Z())/pmt_r;
	}
        //
	for (Int_t i=0; i<MAXPMT; i++) {
		Float_t v[3];
		v[0] = pos_mc.X();
		v[1] = pos_mc.Y();
		v[2] = pos_mc.Z();
		wt[i] = GetWeight(pos2[i], v);
	}
	Float_t maxwt = TMath::MaxElement(MAXPMT, wt);
	Float_t tot_wt = 0.;
	for (Int_t i=0; i<MAXPMT; i++){
		 wt[i] = wt[i] / maxwt; //calculate relative wt
		 tot_wt += wt[i];
	}
  
        // Use a 10 ns window to search 2.2MeV candidate
	Float_t uvx[MAXN12], uvy[MAXN12], uvz[MAXN12];
	Float_t ti[MAXN12], qi[MAXN12];
	Int_t   ci[MAXN12]; //cable
	TVector3   pi[MAXN12]; //position
	Int_t   uvf[MAXN12]; //flag: 0=not cut, 1=first cut, 2=second cut, etc.
	Int_t   ncut;
	Bool_t   pre_t0_set = kFALSE;
	Float_t  t0 = -1;
	Int_t    N12i, N12 = 0;
	Int_t    Nc, Nback, Nlow; // Neff, Nc1, NhighQ, NlowQ ;
	Float_t  trms, phirms, thetam, distan;// ratio;
	distan = pos_del.Mag()/10.; //mm to cm, one event one position basis
	for ( i=0; i<nhits; i++) {
		if ( tiskz[i] < tstart ) continue;
		if ( tiskz[i] > tend-12. ) continue;
                // Calculate hits in 12 ns window  
		N12i = GetNhits(tiskz, i, 12., nhits);
		//cout << "N12i is " << N12i << endl; 
		// Only consider candidates with N12 >= N12TH && N12 <= 50
		if ((N12i < N12TH) || (N12i > 50)) continue;
		if ( pre_t0_set && (tiskz[i] - t0 > 20.) ) {
		    res.N12[res.np]    = N12;  //store the previous peak
		    res.Nc[res.np]     = Nc;
		    res.Nlow[res.np]   = Nlow ;
		    res.Nback[res.np]   = Nback ;
		    res.N100[res.np]   = GetNXX(nhits, tiskz, 100., t0+6.);
		    //res.N100[res.np]   = GetNXX(nhits, tiskz, 50., t0+6.);
		    res.trms[res.np]   = trms ;
		    res.phirms[res.np]   = phirms ;
		    res.thetam[res.np]   = thetam ;
		    res.distan[res.np]   = distan ;
		    res.np ++;
		    if ( res.np >= MAXNP - 1 ) break; 
		    N12 = 0; //after a peak is stored, re-initialize N12
                } 
		if ( N12i <= N12 ) continue;
		// Now this is a "better" peak.
		// Calculate discriminating variables for the new peak.

                N12 = N12i;
                t0 = tiskz[i];
                pre_t0_set = kTRUE;
                ncut = 0;
                for (Int_t j=0; j<N12; j++) {
                    uvx[j] = hitv_x[i+j]; // hit vector in 10 ns
                    uvy[j] = hitv_y[i+j];
                    uvz[j] = hitv_z[i+j];
                    uvf[j] = 0; // NOT CUT YET
                    ti[j] =  tiskz[i+j];
                    qi[j]  = qiskz[i+j];
                    ci[j]  = cabiz[i+j];
                    pi[j]  = pos[i+j];
                }
		// Check low hits
		//for (Int_t j=0; j<9; j++) {
			//Nlow[j] = GetLowHits (ci, N12, wt, 0.5+j*0.05); //last number is accep
			Nlow = GetLowHits (ci, N12, wt, 0.7); //
		//}
		// cluster search
		Nc = ncut;
		GetCluster(uvx, uvy, uvz, uvf, N12, ncut, 3, 0.97); // 14.1 degree
		Nc = ncut - Nc;
		//check backward hits
		for (Int_t j=0; j<N12; j++) uvf[j] = 0;
		ncut = 0;
		Nback = ncut;
		CheckBackHits (uvx, uvy, uvz, uvf, N12, ncut, 90.);
		Nback = ncut - Nback;
		// trms
		trms = GetTrms (ti, uvf, N12); 
		// phirms
		phirms = GetPhiRms (uvx, uvy, uvz, uvf, N12);
		// theta mean
		thetam = GetThetaMean (uvx, uvy, uvz, uvf, N12);
		
	}
	if ( N12 > 0 ) { //then N12 is stored by N12i, and N12i has passed the N12TH
		res.N12[res.np]    = N12; //store the last peak. IMPORTANT!!!
		res.Nc[res.np]     = Nc;
		res.Nlow[res.np]   = Nlow ;
		res.Nback[res.np]   = Nback ;
		res.N100[res.np]   = GetNXX(nhits, tiskz, 100., t0+6.) ; //t0+5 is tcenter
		//res.N100[res.np]   = GetNXX(nhits, tiskz, 50., t0+6.) ; //t0+5 is tcenter
		res.trms[res.np]   = trms ;
		res.phirms[res.np]   = phirms ;
		res.thetam[res.np]   = thetam ;
		res.distan[res.np]   = distan ;
		res.np ++;
//////////////////////////////////////////add for test
//                for (Int_t j=0; j<N12; j++) {
//                  for (Int_t k=j+1; k<N12; k++) {
//                       Float_t da =  uvx[j]*uvx[k]                         
//			       + uvy[j]*uvy[k] + uvz[j]*uvz[k];
//                     // all hits
//                     res.dalpha[res.ncomb] = da;
//                     res.ncomb ++;
//                     if ( res.ncomb > 10000 ) {
//                         cout << "ncomb: " << res.ncomb << endl;
//                         break;
//                     }
//                 }
//               }
/////////////////////////////////////////////////////
		N12 = 0; //after a peak is stored, re-initialize N12
	}
}

Int_t GetNhits(Float_t *v, Int_t start_index, Float_t width, Int_t nhits)
{
     Int_t i = start_index;
      while (1) {
	       i++;
	       if((i > nhits-1 ) || (TMath::Abs((v[i]-v[start_index])) > width)) break;
      }
      return TMath::Abs(i - start_index);
}

Int_t GetCluster (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
		                                            Int_t N12, Int_t &ncut, Int_t ncth, Float_t thr)
{
    Int_t index[N12];
    Int_t nc, nc_m;
    for (Int_t i=0; i<N12; i++) {
        if (  flag[i] != 0 ) continue; // skip hits that are already in clusters 
        for (Int_t j=0; j<N12; j++) index[j] = 0;
        index [i] = 1; // first hit in a cluster
        nc_m = 0;
        while (1) {
            for (Int_t j=0; j<N12; j++) {
                if (  flag[j] != 0 ) continue; // skip hits that are already in clusters 
                // scan index
                for (Int_t k=0; k<N12; k++) {
                    if ( k == j ) continue;       // do not compute angle to itself
                    if ( index[k] == 0 ) continue;// skip non-candidate hit
                    Float_t dalpha = ux[j]*ux[k]+uy[j]*uy[k]+uz[j]*uz[k];
                    if ( dalpha > thr ) index[j] = 1;
                }
            }
            // count hits in cluster
            nc = 0;
            for (Int_t j=0; j<N12; j++) {
                if ( index[j] == 1 ) nc ++;
            }
            if ( nc > nc_m ) nc_m = nc;
            else break; // no more hits belong to this cluster, stop
        }
        if ( nc_m >= ncth ) { // this is a cluster
            for (Int_t j=0; j<N12; j++) {
                if ( index[j] == 1 ) {
                    ncut ++; //the initial value of ncut = 0
                    flag[j] = ncut;
                }
            }
        }
    }
    return 1;
}

Float_t GetWeight (const TVector3  xyz, const Float_t v[3])
{ //v is event position; xyz is PMT position
     const Float_t ATT_LEN = 9000.; // temporary
     Float_t costh, r, w;
     r = sqrt( (xyz.X() - v[0])*(xyz.X() - v[0])
		      + (xyz.Y() - v[1])*(xyz.Y() - v[1])
		      + (xyz.Z() - v[2])*(xyz.Z() - v[2]) );
     //costh = 1. ; //temporary set
     costh = (xyz.X()*(xyz.X() - v[0]) + xyz.Y()*(xyz.Y() - v[1])+xyz.Z()*(xyz.Z() - v[2]))/sqrt(xyz.X()*xyz.X()+xyz.Y()*xyz.Y()+xyz.Z()*xyz.Z())/r;
     costh = fabs(costh); //temporary set fabs for debug purpose
     w = EffCos( costh ) * exp(-r/ATT_LEN)/r/r;
     //w =  exp(-r/ATT_LEN)/r/r;
     return w;
}
     
Float_t EffCos ( Float_t costh ) 
{ 
	if ( costh < -0.001 || costh > 1.001 ) {
		cerr << " Invalid cos theta: " << costh << endl;
		exit (0);
	}
	return 0.205349 + 0.523981*costh + 0.389951*costh*costh
		        - 0.131959 *costh*costh*costh;
}

Int_t GetLowHits (Int_t *ci, Int_t N12, Float_t *wt, Float_t acceptance) 
{
	Float_t wlow  = GetWeightThreshold(wt, acceptance);
	Int_t nlow = 0;
	for (Int_t j=0; j<N12; j++) {
		if ( wt[ ci[j] - 1] < wlow  ) nlow ++;
	}
	return nlow;
}

Float_t GetWeightThreshold (const Float_t *w, const Float_t frac) 
{
    const Int_t MAXPM = 9438;
    Float_t w1[MAXPM], w2[MAXPM];
    Int_t   index[MAXPM];
    Float_t tot = 0.;
    for (Int_t i=0; i<MAXPM; i++) {
        w2[i] = w[i];
        tot += w[i];
     }
    TMath::Sort(MAXPM, w2, index, kTRUE); // In decreasing order
    for (Int_t i=0; i<MAXPM; i++) {
        w1[i] = w2[ index[i] ];
    }
    Float_t t = 0;
    for (Int_t i=0; i<MAXPM; i++) {
        t += w1[i];
        if ( t/tot > frac ) {
            //if ( verbosity > 1 ) {
            //    //cout << " X/Y/Z= " << VX << " " << VY << " " << VZ << endl;
            //    cout << " # of high weight PMTs (frac=" << frac << "): " << i+1 << endl;
            //}
            return w1[i];
        }
    }

    return -1.;
}

Bool_t CheckBackHits (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag,
		                                            Int_t N12, Int_t &ncut, Float_t angle)
{
    const Float_t THR = TMath::Cos(angle*3.14159265/180.);

    // estimate dir
    Float_t dir[3];
    dir[0] = dir[1] = dir[2] = 0.;
    for (Int_t i=0; i<N12; i++) {
        if ( flag[i] != 0 ) continue;
        dir[0] += ux[i];
        dir[1] += uy[i];
        dir[2] += uz[i];
    }
    Float_t r = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0] = dir[0] / r;
    dir[1] = dir[1] / r;
    dir[2] = dir[2] / r;
    // cut hits going backward
    for ( Int_t i=0; i< N12; i++ ) {
        if ( flag[i] != 0 ) continue;
        Float_t cross = dir[0] * ux[i] + dir[1] * uy[i] + dir[2] * uz[i];
        if ( cross < THR ) {
            ncut ++;
            flag[i] = ncut;
        }
    }
    return kTRUE;
}

Float_t GetTrms (Float_t *ti, Int_t *flag, Int_t N12)
{
    Int_t   neff = 0;
    Float_t tmean= 0.;
    for (Int_t j=0; j<N12; j++) {
        if ( flag[j] != 0 ) continue;
        tmean += ti[j];
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    tmean = tmean/neff;
    Float_t trms = 0.;
    for (Int_t j=0; j<N12; j++) {
        if ( flag[j] != 0 ) continue;
        trms += (ti[j] - tmean) * (ti[j] - tmean);
    }
    trms = TMath::Sqrt(trms/neff); 
    return trms;
}
Float_t GetPhiRms (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N12)
{
    // Calculate standard deviation of Phi of hit vectors
    Float_t vect_sum[3], angle_x[N12], angle_y[N12], angle_z[N12];
    Int_t neff = 0;
    for (Int_t j=0; j<N12; j++){
        if ( flag[j] != 0 ) continue;
        angle_x[neff] = ux[j];
        angle_y[neff] = uy[j];
        angle_z[neff] = uz[j];    
        neff ++;
    }
    if ( neff < 2 ) return -1.;
    // calculate sum of hit vectors (estimate dir)
    vect_sum[0] = vect_sum[1] = vect_sum[2] = 0.;
    for (Int_t i=0; i<neff; i++) {
        vect_sum[0] += ux[i];
        vect_sum[1] += uy[i];
        vect_sum[2] += uz[i];
    }
    Float_t vr = sqrt(vect_sum[0]*vect_sum[0] + vect_sum[1]*vect_sum[1] 
                      + vect_sum[2]*vect_sum[2]);
    vect_sum[0] = vect_sum[0] / vr;
    vect_sum[1] = vect_sum[1] / vr;
    vect_sum[2] = vect_sum[2] / vr;
    // calculate phi_div
    TVector3 hitv[neff], pv[neff], dirv;
    Float_t  phiv[neff];
    dirv.SetXYZ(vect_sum[0], vect_sum[1], vect_sum[2]); //sum vector
    for (Int_t j=0; j<neff; j++) {
        Float_t cross;
        Float_t scale;
        cross = vect_sum[0] * angle_x[j] + vect_sum[1] * angle_y[j]
            + vect_sum[2] * angle_z[j];
        scale = 1 / cross;
        hitv[j].SetXYZ (angle_x[j] * scale, angle_y[j] * scale, angle_z[j] * scale);
        //perpendicular vector
        pv[j] = hitv[j] - dirv;
        // Rotate vdir to z axis
        pv[j].Rotate(-dirv.Phi(), TVector3(0., 0., 1.));
        pv[j].Rotate(-dirv.Theta(), TVector3(0., 1., 0.));
        phiv[j] = pv[j].Phi();
        if ( phiv[j] < 0. ) phiv[j] += 2.*TMath::Pi();
    }
    Float_t phiv2[neff];
    Int_t index[neff]; 
    for (Int_t j=0; j<neff; j++) phiv2[j] = phiv[j];
    TMath::Sort(neff, phiv2, index, kFALSE); // Sort phi in increasing order   
    for (Int_t j=0; j<neff; j++){
        phiv[j] = phiv2[ index[j] ];
    }
    Float_t avephi = 2.*TMath::Pi() / neff;
    Float_t delta_phi;
    Float_t phi_div = 0.;
    for (Int_t j=1; j<neff; j++) {
        delta_phi = phiv[j] - phiv[j-1];
        phi_div += (delta_phi - avephi) * (delta_phi - avephi);
    }
    delta_phi = phiv[0] + 2.*TMath::Pi() - phiv[neff-1];
    phi_div += (delta_phi - avephi) * (delta_phi - avephi);
    phi_div = TMath::Sqrt(phi_div/neff);
    return phi_div*180/TMath::Pi();
}

Float_t GetThetaMean (Float_t *ux, Float_t *uy, Float_t *uz, Int_t *flag, Int_t N12)
{
    Int_t neff = 0;
    Float_t dir[3];
    dir[0] = dir[1] = dir[2] = 0.;
    for (Int_t i=0; i<N12; i++) {
        if ( flag[i] != 0 ) continue;
        dir[0] += ux[i];
        dir[1] += uy[i];
        dir[2] += uz[i];
        neff ++;
    }
    if ( neff < 2 ) return 0.;
    Float_t vr = sqrt(dir[0]*dir[0] + dir[1]*dir[1] 
                      + dir[2]*dir[2]);
    dir[0] = dir[0] / vr;
    dir[1] = dir[1] / vr;
    dir[2] = dir[2] / vr;
 
    Float_t theta;
    Float_t mean = 0.;
    for (Int_t i=0; i<N12; i++) {
        if ( flag[i] != 0 ) continue;
        theta = dir[0]*ux[i] + dir[1]*uy[i] + dir[2]*uz[i];
        theta = TMath::ACos(theta)*180./TMath::Pi();
        mean += theta;
    }
    mean = mean / neff;

    return mean;
}

Int_t GetNXX(Int_t nhits, Float_t *t, Float_t twin, Float_t tcenter)
{
    Int_t nxx = 0;
    for (Int_t i=0; i<nhits; i++){
        if ( t[i] < tcenter - twin/2. ) continue;
        if ( t[i] > tcenter + twin/2. ) break;
        nxx ++;
    }
    return nxx;
}

// Function that returns triggered time of an event is nanoseconds
Double_t GetTriggeredEventTime(const RAT::DS::EV& tev) {
  Double_t day = tev.GetUniversalTime().GetDays() ;
  Double_t second = tev.GetUniversalTime().GetSeconds() ;
  Double_t nanosecond = tev.GetUniversalTime().GetNanoSeconds() ;
  return day*24*3600*1e9 + second*1e9 + nanosecond;  
}
