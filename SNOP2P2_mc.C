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

const Float_t C_WATER = 21.6628 ; //speed of light in water
const Int_t   MAXNP   = 100;     // Maximum N12 peaks in 420 ns data
const Int_t   MAXN12  = 200; //Maximum N12
const Int_t   MAXCOMB = 10000;
string fitName = "multipath" ;
//Int_t N12TH=5 ;
Int_t N12TH=2 ;
TTree *theOTree; 

typedef struct
{
	Int_t    np; 
        Int_t    N80M; 
        Float_t  T80M; 
        Int_t    N12[MAXNP];
        //Int_t    Nhits[MAXNP];
        Int_t    Nhits;
	Int_t    Nc[MAXNP];
	Int_t    Nlow[MAXNP];
	Int_t    Nback[MAXNP];
	Int_t    N100[MAXNP];
	Float_t    trms[MAXNP];
	Float_t    phirms[MAXNP];
	Float_t    thetam[MAXNP];
	Float_t    distan[MAXNP];
	Float_t   VVX;
	Float_t   VVY;
	Float_t   VVZ;

	//experimental
	Float_t  dalpha[MAXCOMB];
	Int_t    ncomb; 

        void Clear() 
        {
            np   = 0;
            Nhits = 0;
            VVX   = 0;
            VVY   = 0;
            VVZ   = 0;
            N80M = 0;
            T80M = -9999.;
            for (Int_t i=0; i<MAXNP; i++) {
                N12[i]   = 0;
                //Nhits[i]   = 0;
                Nc[i]   = 0;
                Nlow[i]   = 0;
                Nback[i]   = 0;
                N100[i]   = 0;
                trms[i]   = 0;
                phirms[i]   = 0;
                thetam[i]   = 0;
                distan[i]   = 0;
	    }
            // Experimental
	    ncomb = 0;
            for (Int_t i=0; i<MAXCOMB; i++) dalpha[i]=0;
         }
}results_t;
results_t res; 

void SNOP2P2_mc() 
//int main(int argc, char** argv)
{
//  if ( argc != 3 ) {
//	          cout << "Usage: " << argv[0] << " OUTPUT " << 
//			             " INPUT " << endl;
//	         return -1;
//  }
  // Read in file
  string Infname = "2p2_water_001.root";
  TString OUTfname = "Save_"+TString(Infname);
  //Infname = argv[2];
  //OUTfname = argv[1];
  RAT::DU::DSReader dsreader(Infname);
   //Get pmtInfo
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();

   TFile *fout = new TFile(OUTfname, "RECREATE");
   theOTree = new TTree("sno2p2", "SNOP 2.2 MeV");
   theOTree->Branch("np",    &res.np,     "np/I");
   theOTree->Branch("N12",    res.N12,    "N12[np]/I");
   //theOTree->Branch("Nhits",    res.Nhits,    "Nhits[np]/I");
   theOTree->Branch("Nhits",    &res.Nhits,    "Nhits/I");
   theOTree->Branch("ncomb", &res.ncomb,  "ncomb/I");
   theOTree->Branch("dalpha", res.dalpha, "dalpha[ncomb]/F");
   theOTree->Branch("Nc",    res.Nc,    "Nc[np]/I");
   theOTree->Branch("Nlow",    res.Nlow,    "Nlow[np]/I");
   theOTree->Branch("Nback",    res.Nback,    "Nback[np]/I");
   theOTree->Branch("N100",   res.N100,   "N100[np]/I");
   theOTree->Branch("trms",    res.trms,    "trms[np]/F");
   theOTree->Branch("phirms",    res.phirms,    "phirms[np]/F");
   theOTree->Branch("thetam",    res.thetam,    "thetam[np]/F");
   theOTree->Branch("distan",    res.distan,    "distan[np]/F");
   theOTree->Branch("VVX",    &res.VVX,    "VVX/F");
   theOTree->Branch("VVY",    &res.VVY,    "VVY/F");
   theOTree->Branch("VVZ",    &res.VVZ,    "VVZ/F");
   //cout << "SKY is great************" << endl;
  //Get pmtInfo
  //const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  
  // Loop through events
  for(int i=0; i<dsreader.GetEntryCount();i++){
    if(i%2000==0) cout << "Dealing with " << i << "th event" << endl;
    //if(i>500) cout << "Dealing with " << i << "th event" << endl;
    const RAT::DS::Entry& rds = dsreader.GetEntry(i);
    const RAT::DS::MC& rmc = rds.GetMC();
    const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0); //-11 is e^+
    // get the position
    TVector3 pos_mc0 = rmcparticle.GetPosition();
    double vvx= pos_mc0.X() ;
    double vvy=pos_mc0.Y();
    double vvz=pos_mc0.Z(); 
    //do {
	    vvx = gRandom->Gaus(vvx, 300); //300 mm
            vvy = gRandom->Gaus(vvy, 300);
            vvz = gRandom->Gaus(vvz, 300);
    //}while ( sqrt(vvx*vvx + vvy*vvy+vvz*vvz) > 6000.); //generate vertex within AV
    TVector3 pos_mc(vvx, vvy, vvz);
    //cout << "SKY is great************" << endl;
    //double tmc = rmcparticle.GetTime();
    // The quick way to get the radius
    // (or you could do pow((mcpos.X()*mcpos.X()+mcpos.Y()*mcpos.Y()+mcpos.Z()*mcpos.Z()),0.5)
    //double mc_r = mcpos.Mag();

    // now access the events - we'll loop through them in case there is more than one triggered event per
    // simulated event
    int nevC = rds.GetEVCount();
    for(int iev=0;iev<nevC; iev++){
      // Get the Event information
      const RAT::DS::EV& rev = rds.GetEV(iev);
      if (iev>0) continue; //skip re-trigger events
      int nHits = rev.GetNhits();
      //if(i>580) cout << "nhits is " << nHits << endl;
      //if(nHits<5 || nHits>25) continue;
      if(nHits<2 || nHits>25) continue; //for showing purpose
      //const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
      // make sure the position exist
      if( !rev.FitResultExists(fitName)) continue;
      if( !rev.GetFitResult(fitName).GetValid()) continue; 
      //consider the trigger efficiency here
      /*if(nHits==5) { //trigger effi 1.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      if(vx>0.1) continue ;
      }
      else if(nHits==6) { //trigger effi 5.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>5.) continue ;
	      if(vx>2.) continue ;
      }
      else if(nHits==7) { //trigger effi 30.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>30.) continue ;
	      if(vx>10.) continue ;
      }
      else if(nHits==8) { //trigger effi 82.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>82.) continue ;
	      if(vx>62.) continue ;
      }
      else if(nHits==9) { //trigger effi 95.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>95.) continue ;
	      if(vx>75.) continue ;
      }
      else if(nHits==10) { //trigger effi 95.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>99.5) continue ;
	      if(vx>80.) continue ;
      }
      else if(nHits==11) { //trigger effi 95.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>99.5) continue ;
	      if(vx>90.) continue ;
      }
      else if(nHits==12) { //trigger effi 95.%.
	      double vx = gRandom->Uniform(0., 100.); 
	      //if(vx>99.5) continue ;
	      if(vx>92.) continue ;
      }
      else { } */ //do nothing to Nhits>=10
      //if(nHits<=25) { 
      //	      double vx = gRandom->Uniform(0., 100.);
      //	      if(vx>60.) continue ; //this number is data driven
      //}
      // Clear results of previous event
      res.Clear();
      //search
      //NeutronSearch(rev, pmtInfo, 0., 420., pos_mc);
      NeutronSearch(rev, pmtInfo, -50., 420., pos_mc);
      // Fill the output tree 
      res.Nhits = nHits ; //one events one nhits
      res.VVX   = vvx ;
      res.VVY   = vvy ;
      res.VVZ   = vvz ;
      //for(Int_t npk=0; npk<res.np;npk++){
      //	   res.Nhits[npk]  = nHits ;
      //}
      theOTree->Fill();
    } //triggered event loop
  } //entry loop
  //save result
  fout->Write();
  delete theOTree;
  fout->Close();
  delete fout;
  theOTree = 0;
  fout = 0; 
  //need a check flag
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
	Float_t  trms, phirms, thetam, distan, VVX, VVY, VVZ;// ratio;
	distan = pos_del.Mag()/10.; //mm to cm, one event one position basis
	VVX = pos_mc.X(); 
	VVY = pos_mc.Y(); 
	VVZ = pos_mc.Z(); 
	for ( i=0; i<nhits; i++) {
		if ( tiskz[i] < tstart ) continue;
		if ( tiskz[i] > tend-12. ) continue;
                // Calculate hits in 12 ns window  
		N12i = GetNhits(tiskz, i, 12., nhits);
		//cout << "N12i is " << N12i << endl; 
		// Only consider candidates with N12 >= N12TH && N12 <= 50
		if ((N12i < N12TH) || (N12i > 50)) continue;
		//if ( pre_t0_set && (tiskz[i] - t0 > 20.) ) {
		if ( pre_t0_set && (tiskz[i] - t0 > 800.) ) {
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
		    //res.VVX[res.np]   = VVX ;
		    //res.VVY[res.np]   = VVY ;
		    //res.VVZ[res.np]   = VVZ ;
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
		//res.VVX[res.np]   = VVX ;
		//res.VVY[res.np]   = VVY ;
		//res.VVZ[res.np]   = VVZ ;
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
