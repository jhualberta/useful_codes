//For Rat-v6.5.2 SNOP data
//2017-6-25
//#include <RAT/DataCleaningUtility.hh>
#include <RAT/DS/Meta.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <vector>
#include "TH2.h"
#include "TH1.h"
#include <TVector3.h>
#include "TF1.h" 
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = 
fit9 = FTP
*/
using namespace std ;
const double ITRval = 0.55;
class MP
{
public:
  void multiPathLight(double *pmt, double *par, vector<double> &results );
private:
  double ScintExWaterPathCalculation(double *incidentVtx, double *pmtvtx );
  double Factor( double transver, double dh );
  double FactorN( double transver, bool pmtBot );
  double fHeight;
  double fDepth;
  TMatrixT<double> fFactor;
};


double fSpeedOfLight = 299.;
double fWaterLevel = 0;
double fScintRIeff = 1.643;
double fWaterRIeff = 1.38486;
double fScintRI = 1.50;//50;
double fWaterRI = 1.33;

void lightPathMP()
{
   const char* filename = "FitRatMaster_mc_TrackOn_water_5MeVgamma_x0y0z-4000_10evts_WL12000.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos;
   bool scint = false;
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   // RAT::DU::Utility::Get()->GetLightPathCalculator() must be called *after* the RAT::DU::DSReader constructor.
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
   const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
   TH1F *hPMTXvsTof_ratMC = new TH1F("hPMTXvsTof_ratMC","pmt X vs tof, rat MC", 400, -9000, 9000);
   TH1F *hPMTYvsTof_ratMC = new TH1F("hPMTYvsTof_ratMC","pmt Y vs tof, rat MC", 400, -9000, 9000);
   TH1F *hPMTZvsTof_ratMC = new TH1F("hPMTZvsTof_ratMC","pmt Z vs tof, rat MC", 400, -9000, 9000);

   TH1F *hPMTXvsTof_mpMC = new TH1F("hPMTXvsTof_mpMC","pmt X vs tof, mp MC", 400, -9000, 9000);
   TH1F *hPMTYvsTof_mpMC = new TH1F("hPMTYvsTof_mpMC","pmt Y vs tof, mp MC", 400, -9000, 9000);
   TH1F *hPMTZvsTof_mpMC = new TH1F("hPMTZvsTof_mpMC","pmt Z vs tof, mp MC", 400, -9000, 9000);

   TH2F *hPMT_pos = new TH2F("hPMT_pos", "hit pmt pos, cosTheta vs Phi ",200,-1,1,628,-TMath::Pi(), TMath::Pi());
   TH2F *hPMT_biasTof = new TH2F("hPMT_biasTof", "pmt pos with delta tof>0.5 ns for mp and rat",200,-1,1,628,-TMath::Pi(), TMath::Pi());

   TVector3 u_fit, pos_fit;
   Double_t theta_e;
   // Double_t grVelocity = 2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
   // if(scint) grVelocity = 183.503698029;//labppo 
   Double_t rPSUP = 8390;
   string process1 = "Scintillation" ;
   string process2 = "Cerenkov" ;
   string process3 = "OpAbsorption" ;
   Double_t energy, wavelength ;
 
   Double_t countFitValid = 0;//total fitted events
   Double_t countFECDtotal =0;//total FECD==9188 events
   Double_t countSuccess = 0;
   Double_t countTimeWindow = 0;//fitValid && FECD && tFit cuts
   Double_t countClean = 0;
   Double_t countTrig = 0;
   Double_t countTrigFECD = 0;
   Double_t countITRbeta14 = 0;
   Double_t countNoCleanTotal=0; Double_t countNoCleanTrig=0;Double_t countNoCleanFECD=0;
   int trigWord = 0x11;//trig number,1<<6, http://www.snoplus.ca/docs/rat/user_manual/html/node47.html#t:trigword
   bool trigCut = false;//true;
   int nhitCut = 0;//15;
   unsigned int fecdID = 9190;
   TString fitName = "partialFitter";
   //TString specFit = "fit10";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   TVector3 mcPos;
   size_t evtnum = 2;//dsReader.GetEntryCount(); 
   for( size_t iEntry = 1; iEntry < evtnum; iEntry++)
   {
    std:cout << " event ID "<< iEntry <<std::endl ;
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
    Int_t nevC = rDS.GetEVCount();//!! if retriggered events, nevC == 2
    const RAT::DS::MC& rmc= rDS.GetMC();  
    const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
    mcPos = rmcparticle.GetPosition();
    double tmc = rmcparticle.GetTime();
    if(nevC>1) nevC = 1;//!!! remove retrigger events
    for(Int_t iev=0;iev<nevC; iev++) {
     /// Get the Event information    
     const RAT::DS::EV& rev = rDS.GetEV(iev);
     std::vector<std::string> fitname = rev.GetFitNames();
     std::vector<std::string>::iterator it;
     double nhits = rev.GetNhits();
     int trig = rev.GetTrigType();//get triggerType
     const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
     if(rev.FitResultExists("partialFitter")) {//find partialFitter exists
        RAT::DS::FitVertex fitVertex = rev.GetFitResult("partialFitter").GetVertex(0);
        //if(rev.GetFitResult("partialFitter").GetValid())//Global Validity of the fitter !!! NOTE: this is different
        {
          pos_fit=fitVertex.GetPosition(); 
 	  double tfit = fitVertex.GetTime();
 	  double par[4] = {mcPos.X(), mcPos.Y(), mcPos.Z(), tmc};
	  if(fitVertex.ValidPosition())//NOTE: check fitPosition Valid !!
 	  { 
            std::cout<<"fit Pos "<<pos_fit.X()<<", "<<pos_fit.Y()<<", "<<pos_fit.Z()<<std::endl;
            std::cout<<"mc Pos "<<mcPos.X()<<", "<<mcPos.Y()<<", "<<mcPos.Z()<<std::endl;
            for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt+= 1)
            {
              TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
              double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
              cout<<"("<<pmtpos.X()<<", "<<pmtpos.Y()<<", "<<pmtpos.Z()<<", "<<hitTime<<") id "<<calpmts.GetPMT(ipmt).GetID()<<" ipmt "<<ipmt<<endl;
              lightPath.CalcByPositionPartial( mcPos, pmtpos );
              double distInInnerAV = lightPath.GetDistInInnerAV();
              double distInAV = lightPath.GetDistInAV();
              double distInWater = lightPath.GetDistInWater();
              const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater );
              cout<<"mc tof "<<transitTime<<endl;
              lightPath.CalcByPositionPartial( pos_fit, pmtpos );
              double distInInnerAV_fit = lightPath.GetDistInInnerAV();
              double distInUpperTarget = lightPath.GetDistInUpperTarget();
              double distInLowerTarget = lightPath.GetDistInLowerTarget();
              RAT::DU::EffectiveVelocity effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();
              effectiveVelocity.BeginOfRun();
              const double transitTimeFit = effectiveVelocity.CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );
              cout<<"fit tof "<<transitTimeFit<<endl;
              double pmt[4] = {pmtpos.X(),pmtpos.Y(),pmtpos.Z(),hitTime};
              MP fitMP;
              vector<double> mpResults;
              fitMP.multiPathLight(pmt, par, mpResults);
              double tof = mpResults[0];
              hPMTXvsTof_ratMC->Fill(pmtpos.X(), transitTime);
              hPMTYvsTof_ratMC->Fill(pmtpos.Y(), transitTime);
              hPMTZvsTof_ratMC->Fill(pmtpos.Z(), transitTime);
              hPMTXvsTof_mpMC->Fill(pmtpos.X(), tof);
              hPMTYvsTof_mpMC->Fill(pmtpos.Y(), tof);
              hPMTZvsTof_mpMC->Fill(pmtpos.Z(), tof);
              hPMT_pos->Fill(pmtpos.CosTheta(), pmtpos.Phi());
              if( abs(tof-transitTime)>0.5 ) hPMT_biasTof->Fill(pmtpos.CosTheta(), pmtpos.Phi());

              cout<<endl;
             }
           } // valid position 
        } // global valid
      } // fit exists
   }//if triggeretd events
  }//for loop
   TFile *f = new TFile("Results_tof.root","recreate");
   f->cd();
   hPMTXvsTof_ratMC->Write();
   hPMTYvsTof_ratMC->Write();
   hPMTZvsTof_ratMC->Write();
   hPMTXvsTof_mpMC->Write();
   hPMTYvsTof_mpMC->Write();
   hPMTZvsTof_mpMC->Write();
   hPMT_pos->Write(); hPMT_biasTof->Write();
   f->Close();
}

double MP::ScintExWaterPathCalculation(double *incidentVtx, double *pmtvtx )
{
  TVector3 startpos, pmtpos, ndiff;
  startpos.SetXYZ(incidentVtx[0],incidentVtx[1],incidentVtx[2]);
  pmtpos.SetXYZ(pmtvtx[0],pmtvtx[1],pmtvtx[2]);
  ndiff = (pmtvtx - startpos).Unit();
  double pathInScint = 0;
  double sqrVal = pow((startpos * ndiff),2) - (startpos).Mag2() + 6050*6050; 
  if(sqrVal<0) pathInScint = 0;
  else {
    pathInScint = -(startpos * ndiff) + sqrt(sqrVal);
    if(pathInScint>6050)
    pathInScint = -(startpos * ndiff) - sqrt(sqrVal);
 
 }
  return pathInScint;
}

double MP::Factor( double transver, double dh ) {

  /// numerically calculate transverse distance
  /// if PMT above water, vertex below, input (transver,dh) = (transver/h, fDepth/fHeight), then:
  /// TransverseDist = fDepth*tan( theta ) + fHeight*tan(asin( fWaterRI/fScintRI*sin( theta ) ))
  /// if PMT below water, vertex above, input (transver,dh) = (transver/-fDepth, fHeight/fDepth), then:
  /// TransverseDist = (-fHeight)*tan( theta ) + (-fDepth)*tan(asin( fWaterRI/fScintRI*sin( theta ) ))
  /// param: fDepth -- [0]; fWaterRI/fScintRI -- [1];
  /// Find theta (always the angle in water)
  TF1 fTransverseDist( "fTransverseDist", "[0]*tan( x )+ 1*tan( asin( [1]*sin( x ) ) )", 0, 0.850908);
 
  fTransverseDist.SetParameter( 0, dh);
  fTransverseDist.SetParameter( 1, fWaterRI/fScintRI );
  
  /// parameter is the vertical distance between PMT and vertex
  double theta = fTransverseDist.GetX(transver);
  /// starting point is below intersection
  if( fabs( transver ) < 0.001 ) return 0;
  /// fraction of straight line distance to point directly below intersection.
  else return dh*tan( theta )/transver;
}

double MP::FactorN( double transver, bool pmtBot ) {
  /// calculates fraction of transverse distance from fDepth to the surface by interpolating, based on fDepth and fHeight
  double alpha = 0, dh = 0, lx = 0;
  if( fabs(fHeight) < 0.1 || fabs(fDepth) < 0.1 ) {
    if (fabs(fHeight) < 0.1)
    { 
      if( pmtBot == false) alpha = 1;
      else alpha = 0;
      //std::cout<<"case 1 "<<std::endl;
    }
    else if (fabs(fDepth) < 0.1) {
      if (pmtBot == false) alpha = 0;
      else alpha = 1;
      //std::cout<<"case 2 "<<std::endl;
    }
  } 
  else {
    if ( pmtBot == false ) { // if the PMT is above the water
      dh = fDepth/fHeight;
      // transver is projection of vector (vertex to PMT) onto detector z-axis
      // lx is the log of the ratio of transver and fHeight to define the extreme conditions
      lx = log(transver/fHeight);
    }
    else if ( pmtBot == true ) { // if the PMT is below the water
      dh = fHeight/fDepth;
      lx = log(transver/(-fDepth));
    }
    //std::cout<<"case 3 "<<std::endl;
    double ldh = log(dh);
    if( lx < -5 ) alpha = ( fWaterRI*dh + fScintRI )/( 1 + dh ); // sVertex close to water level
    else {
      if( lx>5 ) lx = 5; // vetex too low
      if( ldh<-5 ) ldh = -5; // water level too low
      else if( ldh > 5 ) ldh = 5; // water level too high
      double ux = ( lx + 5 )/.2;
      int ix = floor(ux);
      ux = ux-ix;
      double udh = ( ldh + 5 )/.2;
      int idh = floor( udh );
      udh = udh-idh;
      alpha = ( 1-ux )*( ( 1-udh )*fFactor[ix][idh] + udh*fFactor[ix][idh+1] ) + ux*( ( 1-udh )*fFactor[ix+1][idh]+udh*fFactor[ix+1][idh+1]);
    }
  }
  std::cout<<"factorN  fH, fD "<<fHeight<<" "<<fDepth<<" "<<alpha<<std::endl;
  return alpha;
}

void MP::multiPathLight(double *pmt, double *par, vector<double> &results)
{
  using namespace RAT;
  using namespace ROOT;
  results.clear();
 
  double fMaxReflects = 500;
  double fWaterLevel = 0;
  double fSpeedOfLightWater = fSpeedOfLight/fWaterRIeff;
  double fSpeedOfLightScint = fSpeedOfLight/fScintRIeff;

  fFactor.ResizeTo( 50 + 2, 50 + 2 );

  for( int ix = 0; ix < 51; ix++ ) {
    double lx = -5 + ix*0.2;
    double x = exp( lx );
    for( int idh = 0; idh < 51; idh++ ) {
      double ldh = -5 + idh*0.2;
      double dh = exp( ldh );
      double f = Factor( x, dh );
      fFactor[ix][idh] = f;
      //std::cout<<"fFactor["<<ix<<"]["<<idh<<"] = "<<f<<std::endl;
    }
  }

  std::vector<double> fRWater;
  std::vector<double> fRScint;
  /// calculate cosTheta, reflection and transimission probabilities (Fresnel equations)
  double step_ct = 1./fMaxReflects;
  for ( double cosTi = 0; cosTi < 1.; cosTi += step_ct ) { // cosThetaIncident
    double Rp = 0., Rs = 0.; // reflection prob for p-wave and s-wave
    double sinTi = sqrt( 1 - cosTi*cosTi ); // sinThetaIncident square
    // incident in water, transmit in scint
    double cosTt = 1 - (fWaterRI/fScintRI*sinTi)*(fWaterRI/fScintRI*sinTi);
    if( cosTt < 0 ) { // total internal reflection
      Rp = 1;
      Rs = 1;
      fRWater.push_back( ( Rp + Rs )/2. );
      //fTWater.push_back( ( 1-Rp + 1-Rs )/2. );
    } else { // transmission(refraction) and reflection
      cosTt = sqrt( cosTt ); // cosThetaTransmit
      Rs = pow( ( fWaterRI*cosTi - fScintRI*cosTt )/( fWaterRI*cosTi + fScintRI*cosTt ), 2 );
      Rp = pow( ( fWaterRI*cosTt - fScintRI*cosTi )/( fWaterRI*cosTt + fScintRI*cosTi ), 2 );
      fRWater.push_back( ( Rp + Rs )/2. );
      //fTWater.push_back( ( 1-Rp + 1-Rs )/2. );
    }
    // incident in scint, transmit in water
    Rs = 0, Rp = 0, cosTt = 0;
    cosTt = 1 - (fScintRI/fWaterRI*sinTi)*(fScintRI/fWaterRI*sinTi);
    if( cosTt < 0 ) { // total internal reflection
      Rp = 1;
      Rs = 1;
      fRScint.push_back( ( Rp + Rs )/2. );
      //fTScint.push_back( ( 1-Rp + 1-Rs )/2. );
    } else {
    Rs = pow( ( fScintRI*cosTi - fWaterRI*cosTt )/( fScintRI*cosTi + fWaterRI*cosTt ), 2 );
    Rp = pow( ( fScintRI*cosTt - fWaterRI*cosTi )/( fScintRI*cosTt + fWaterRI*cosTi ), 2 );
    fRScint.push_back( ( Rp + Rs )/2. );
    //fTScint.push_back( ( 1-Rp + 1-Rs )/2. );
   }
  } // end of for loop
  fRScint.push_back( 0 );
  fRWater.push_back( 0 );

/// loop pmts 
  double dx = pmt[0] - par[0];
  double dy = pmt[1] - par[1];
  double dz = pmt[2] - par[2];
  double dr = sqrt(dx*dx + dy*dy + dz*dz);

  double tDiff = pmt[3] - par[3];// tDiff = hitTime - trialTime; fRes = tDiff - TOF

  TVector3 fVertex; // trial vertex
  fVertex.SetXYZ( par[0], par[1], par[2] );
  TVector3 fPMTpos;
  fPMTpos.SetXYZ( pmt[0], pmt[1], pmt[2] );
  TVector3 fPosDiff; // pmtPos - vertex, NOT diff! diff is sBeta
  fPosDiff.SetXYZ( dx, dy, dz );
  MP::fDepth = fWaterLevel - par[2];//fVertex.Z();
  MP::fHeight = pmt[2] - fWaterLevel;
  double fAVRadius = 8390;
  double scintpath = 0;
  if( fWaterLevel >= -fAVRadius && fWaterLevel <= fAVRadius ) { // check valid water height
    double tof = 0, tof2 = 0, reflected = 0, alpha = 0;// refracted = 0
    if( fDepth>0 ) {  // vertex is below water level
        if( fHeight<0 ) {  //PMT is below water level
         std::cout<<"MP vbpb"<<std::endl;
        /// geo 1: both PMT and vertex are in water: find direct and reflected light paths 
         tof = dr/fSpeedOfLightWater; // direct time
         alpha = fDepth/(fDepth-fHeight);  //find fractional transverse distance between vertex and water intersection
         TVector3 s = (1-alpha)*fVertex + alpha*fPMTpos;  //find a point directly below the intersection with water
         s.SetZ(fWaterLevel);  //move intersection to water surface.
         TVector3 incident = s-fVertex;
         TVector3 reflectedRay = fPMTpos-s;
         double treflected = (incident.Mag()+reflectedRay.Mag())/fSpeedOfLightWater;  //time of flight for reflected ray
         incident = incident.Unit();
         double ct = incident.Z()*fMaxReflects;  // find cos(theta) between vertical and photon trajectory, multiply by fMaxReflects to look up in table
         int ix=(int)ct; double u=ct-ix;
         reflected = (1-u)*fRWater[ix]+u*fRWater[ix+1];  // find reflection probability for this cos(theta), interpolated
         cout<<"direct tof "<<tof<<" reflected prob "<<reflected<<" reflec tof "<<treflected<<endl;
        }
        else {
        /// geo 2: PMT is above water, vertex below: find refraction path
          std::cout<<"MP vbpt "<<std::endl;
            alpha = FactorN(fPosDiff.Perp(), false);
            TVector3 s = (1-alpha)*fVertex + alpha*fPMTpos;  //find a point directly below the intersection with water
            s.SetZ(fWaterLevel);  //move intersection to water surface.
            TVector3 incident = s - fVertex;
            TVector3 refractedRay = fPMTpos - s;
            double par1[3] = {s.X(),s.Y(),s.Z()};
            scintpath = ScintExWaterPathCalculation(par1,pmt);
            tof = incident.Mag()/fSpeedOfLightWater + scintpath/fSpeedOfLightScint + (refractedRay.Mag()-scintpath)/fSpeedOfLightWater;
            /// no refraction, straight line calculation
            double al = fDepth/(fDepth+fHeight);
            double I = al*dr;
            double T = dr - I;
            TVector3 s2 = (1-al)*fVertex + al*fPMTpos;
            double par2[3] = {s2.X(),s2.Y(),s2.Z()};
            double sp2 = ScintExWaterPathCalculation(par2,pmt);
            tof2 = I/fSpeedOfLightWater + sp2/fSpeedOfLightScint + (T - sp2)/fSpeedOfLightWater;
            std::cout<<"MP tof "<<tof<<" = "<<incident.Mag()/fSpeedOfLightScint<<" + "<<refractedRay.Mag()/fSpeedOfLightScint<<std::endl;
            std::cout<<"MP simple "<<tof2<<" = "<<I/fSpeedOfLightScint<<" + "<<T/fSpeedOfLightScint<<std::endl;
         }
       } else { // vertex is above water level (fDepth<0)
        if(fHeight>0){ // PMT is above water level  
        /// geo 3: both vertex and PMT are in scintillator: find direct and reflected light paths
            std::cout<<"MP vtpt "<<std::endl;
            scintpath = ScintExWaterPathCalculation(par,pmt);
            if(scintpath == 0) tof = dr/fSpeedOfLightWater;
            else tof = scintpath/fSpeedOfLightScint + (dr-scintpath)/fSpeedOfLightWater;
            std::cout<<tof<<std::endl;
        }else {  // PMT is below water level
         /// geo 4: vertex is above, PMT is below: find refracted path
            std::cout<<"MP vtpb "<<std::endl;
            alpha = FactorN(fPosDiff.Perp(), true);
            TVector3 s = alpha*fVertex + (1-alpha)*fPMTpos;  //find a point directly below the intersection with water
            s.SetZ(fWaterLevel);  //move intersection to water surface.
            TVector3 incident = s - fVertex;
            TVector3 refractedRay = fPMTpos - s;
            tof = incident.Mag()/fSpeedOfLightScint + refractedRay.Mag()/fSpeedOfLightWater;  //time of flight for refracted ray
            double al = -fDepth/(-fDepth - fHeight);
            double I = al*dr, T = (1-al)*dr;
            tof2 = I/fSpeedOfLightScint + T/fSpeedOfLightScint; 
            std::cout<<incident.Mag()/fSpeedOfLightScint<<" "<<refractedRay.Mag()/fSpeedOfLightWater<<std::endl;
            std::cout<<"MP tof "<<tof<<" = "<<incident.Mag()/fSpeedOfLightScint<<" + "<<refractedRay.Mag()/fSpeedOfLightScint<<std::endl;
            std::cout<<"MP simple "<<tof2<<" = "<<I/fSpeedOfLightScint<<" + "<<T/fSpeedOfLightScint<<std::endl;
        }
      } //vertex above
        results.push_back(tof);
        results.push_back(tof2);
       //std::cout<<"tof mc "<<tofmc<<std::endl;
      }//valid calculation
}
