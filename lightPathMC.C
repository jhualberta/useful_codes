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
  void multiPathLight(double *pmt, double *par);
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
double fScintRIeff = 1.59;
double fWaterRIeff = 1.38486;
double fScintRI = fScintRIeff;
double fWaterRI = fWaterRIeff;

void lightPathMC()
{
   const char* filename = "FitRat_MC_partialscint_2p5MeVbeta_x0y0z2500_level0_10evts.root";//FitRat_TrackOn_MC_partialscint_2p5MeVbeta_x0y0z100_level0_100evts.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos;
   bool scint = false;
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   // RAT::DU::Utility::Get()->GetLightPathCalculator() must be called *after* the RAT::DU::DSReader constructor.
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
   const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity

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
   bool trigCut = true;//true;
   int nhitCut = 0;//15;
   unsigned int fecdID = 9190;
   TString fitName = "partialFitter";
   //TString specFit = "fit10";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   TVector3 mcPos;
   for( size_t iEntry = 0; iEntry < 2; iEntry++)//dsReader.GetEntryCount(); iEntry++) // dsReader.GetEntryCount(); iEntry++ )
   {
    std:cout << " event ID "<< iEntry <<std::endl ;
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
    Int_t nevC = rDS.GetEVCount();//!! if retriggered events, nevC == 2
    const RAT::DS::MC& rmc= rDS.GetMC();  
    const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
    mcPos = rmcparticle.GetPosition();
    double tmc = rmcparticle.GetTime();
    //cout<<"trigger events "<<nevC<<endl;
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
            for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt+= 10)
            {
              TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
              double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
              cout<<"("<<pmtpos.X()<<", "<<pmtpos.Y()<<", "<<pmtpos.Z()<<", "<<hitTime<<endl;
              lightPath.CalcByPosition( mcPos, pmtpos );
              double distInInnerAV = lightPath.GetDistInInnerAV();
              double distInAV = lightPath.GetDistInAV();
              double distInWater = lightPath.GetDistInWater();
              const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
              cout<<"mc tof "<<transitTime<<endl;
              lightPath.CalcByPositionPartial( pos_fit, pmtpos );
              double distInInnerAV_fit = lightPath.GetDistInInnerAV();
              double distInAV_fit = lightPath.GetDistInAV();
              double distInWater_fit = lightPath.GetDistInWater();
              const double transitTimeFit = groupVelocity.CalcByDistance( distInInnerAV_fit, distInAV_fit, distInWater_fit ); // Assumes a 400nm photon
              cout<<"fit tof "<<transitTimeFit<<endl;
              double pmt[4] = {pmtpos.X(),pmtpos.Y(),pmtpos.Z(),hitTime};
              MP fitMP;
              fitMP.multiPathLight(pmt, par);
              cout<<endl;
             }
           } // valid position 
        } // global valid
      } // fit exists
   }//if triggeretd events
  }//for loop

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

void MP::multiPathLight(double *pmt, double *par)
{
  using namespace RAT;
  using namespace ROOT;
  
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
    double tof = 0, reflected = 0, alpha = 0;// refracted = 0
    if( fDepth>0 ) {  // vertex is below water level
        if( fHeight<0 ) {  //PMT is below water level
         std::cout<<"MP vbpb"<<std::endl;
        /// geo 1: both PMT and vertex are in water: find direct and reflected light paths 
            tof = dr/fSpeedOfLightWater; // direct time
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
            double tof2 = 0;
            double al = fDepth/(fDepth+fHeight);
            double I = al*dr;
            double T = dr - I;
            TVector3 s2 = (1-al)*fVertex + al*fPMTpos;
            double par2[3] = {s2.X(),s2.Y(),s2.Z()};
            double sp2 = ScintExWaterPathCalculation(par2,pmt);
            tof2 = I/fSpeedOfLightWater + sp2/fSpeedOfLightScint + (T - sp2)/fSpeedOfLightWater;
            std::cout<<"MP simple "<<tof2<<std::endl;
         }
       } else { // vertex is above water level (fDepth<0)
        if(fHeight>0){ // PMT is above water level  
        /// geo 3: both vertex and PMT are in scintillator: find direct and reflected light paths
            std::cout<<"MP vtpt "<<std::endl;
            scintpath = ScintExWaterPathCalculation(par,pmt);
            if(scintpath == 0) tof = dr/fSpeedOfLightWater;
            else tof = scintpath/fSpeedOfLightScint + (dr-scintpath)/fSpeedOfLightWater;
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
            double tof2 = I/fSpeedOfLightScint + T/fSpeedOfLightScint; 
            std::cout<<s.X()<<" "<<s.Y()<<" "<<s.Z()<<std::endl;
            std::cout<<"MP simple "<<tof2<<std::endl;
        }
      } //vertex above
       std::cout<<"MP tof "<<tof<<std::endl;
       //std::cout<<"tof mc "<<tofmc<<std::endl;
      }//valid calculation
}
