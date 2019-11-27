////////////////////////////////////////////////////////////////////////
///
/// \class
///
/// \brief
///
/// \author Kalpana Singh  kalpana.singh@ualberta.ca
///
/// REVISION HISTORY:\n
///
/// \detail
////
/////////////////////////////////////////////////////////////////////////*


#include <RAT/MultiPathFitter.hh>
#include <RAT/FitterWLS.hh>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <vector>
#include "TVector3.h"
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

using std::vector;
using std::string;
using namespace RAT;
using namespace RAT::DS;

namespace RAT{

 std::vector<double> FitterWLS::sPDFData_Ch ;
 std::vector<double> FitterWLS::sPDFX_Ch ;
 std::vector<double> FitterWLS::sPDF_Ch ;
 std::vector<double> FitterWLS::sDerivative_Ch ;

 std::vector<double> FitterWLS::sPDFData_WLS ;
 std::vector<double> FitterWLS::sPDFX_WLS ;
 std::vector<double> FitterWLS::sPDF_WLS ;
 std::vector<double> FitterWLS::sDerivative_WLS ;

 std::vector<double> FitterWLS::sPDFData_WLSAngle ;
 std::vector<double> FitterWLS::sPDFX_WLSAngle ;
 std::vector<double> FitterWLS::sPDF_WLSAngle ;
 std::vector<double> FitterWLS::sDerivative_WLSAngle ;

 std::vector<double> FitterWLS::sPDFData_ChAngle ;
 std::vector<double> FitterWLS::sPDF_ChAngle ;
 std::vector<double> FitterWLS::sPDFX_ChAngle ;
 std::vector<double> FitterWLS::sDerivative_ChAngle ;

 TVector3 FitterWLS::sVertex; // trial vertex
 double FitterWLS::sTime0; //trial fit time
 double FitterWLS::sZenith ;
 double FitterWLS::sAzimuth ;

 TVector3 FitterWLS::sStartVertex; // trial vertex
 double FitterWLS::sStartTime0; //trial fit time
 double FitterWLS::sStartZenith ;
 double FitterWLS::sStartAzimuth ;
 int FitterWLS::entriesAngle = 200 ;
 int FitterWLS::entriesTime = 1600 ;

FitterWLS::FitterWLS( double aTime, const TVector3 &aPosition, const int aPMTID,  int nPar): MultiPathFitter(aTime, aPosition, aPMTID, nPar)
 {
  fT = aTime ;   
  sTrialparArray[0]= &sVertex[0] ; 
  sTrialparArray[1]= &sVertex[1] ; 
  sTrialparArray[2]= &sVertex[2] ; 
  sTrialparArray[3] = &sTime0 ; 
  sTrialparArray[4] = &sZenith ; 
  sTrialparArray[5] = &sAzimuth; 
  CLHEP::HepRandom::setTheSeed( sEvent+1234567 );

 }

 void FitterWLS::SetInitParameters()
 {

  //CLHEP::HepRandom::setTheSeed( sEvent+1234567 );
  double ran0 = CLHEP::RandFlat::shoot( 0.,1. );
  double ran1 = CLHEP::RandFlat::shoot( -1.,1. );
  double ran2Pi = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );

  sStartVertex.SetMagThetaPhi( pow( ran0,1.0/3 )*7000,acos( ran1 ),ran2Pi );

  sStartTime0 = CLHEP::RandFlat::shoot(100,300);  //MC
 
  //sStartTime0 = sStartTime0 = CLHEP::RandFlat::shoot(-300,-100);  //Data

  sStartZenith = CLHEP::RandFlat::shoot( 0.,CLHEP::pi ); 
  sStartAzimuth = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );

  //sStartVertex.SetXYZ(0.0,0.0,0.0);
  //sStartZenith = TMath::Pi()/2.0 ;
  //sStartAzimuth = 0.0 ;

  sStartparArray[0]= &sStartVertex[0] ;
  sStartparArray[1]= &sStartVertex[1] ;
  sStartparArray[2]= &sStartVertex[2] ;
  sStartparArray[3]= &sStartTime0 ;
  sStartparArray[4]= &sStartZenith ;
  sStartparArray[5]= &sStartAzimuth ;

 }

 void FitterWLS::BeginOfRun(DS::Run& run) 
 {
    DB* db = DB::Get();
    DBLinkPtr dbLink = db->GetLink( "FIT_WLS_Time" );
    sPDFData_WLS = dbLink->GetDArray( "sPDFTime_WLSFit" );

    DBLinkPtr dbLink1 = db->GetLink( "FIT_WLS_Angle" );
    sPDFData_WLSAngle = dbLink1->GetDArray("sPDFAngle_WLSFit");

    DBLinkPtr effDB = db->GetLink( "EFFECTIVE_VELOCITY" );
    effDB = db->GetLink( "EFFECTIVE_VELOCITY", "lightwater_sno" );
    sSpeedOfLightWater = effDB->GetD( "inner_av_velocity" );    

    DBLinkPtr dbLink2 = db->GetLink( "FIT_Ch_Time" );
    sPDFData_Ch = dbLink2->GetDArray( "sPDFTime_ChFit" );

    DBLinkPtr dbLink3 = db->GetLink("FIT_Ch_Angle");
    sPDFData_ChAngle = dbLink3->GetDArray("sPDFAngle_ChFit");
  
    int j;
    entriesTime = sPDFData_WLS.size(); 
    
    double xtemp = sPDFData_WLS[0];
    for( j = entriesTime-1; j >= 0; j-- ){
	if(sPDFData_WLS[j] == sPDFData_WLS[0]){
	  sPDFData_WLS[j] = xtemp;
	  xtemp -= 0.1;
	}
   }
   
    xtemp = sPDFData_WLS[ entriesTime-1 ];
    for( j = entriesTime-1; sPDFData_WLS[j] == xtemp; j-- );
    for(j++;j<entriesTime ;j++){
	 sPDFData_WLS[j] = xtemp; xtemp -= 0.1;
    }
  
    for ( j = 0; j < entriesTime ; j++ ) {
	sPDFX_WLS.push_back(0.25 * static_cast<Float_t> (j + 1) - 100.25);
	sPDF_WLS.push_back(sPDFData_WLS[j]);
    }
      
    for ( j = 0; j < entriesTime - 1; j++ ) {
	sDerivative_WLS.push_back( ( sPDF_WLS[j + 1] - sPDF_WLS[j] ) * 4.0 );
    }
   
    entriesTime = sPDFData_Ch.size();

    xtemp = sPDFData_Ch[0];
    for( j = entriesTime-1; j >= 0; j-- ){
        if(sPDFData_Ch[j] == sPDFData_Ch[0]){
          sPDFData_Ch[j] = xtemp;
          xtemp -= 0.1;
      }
    }
      
    xtemp = sPDFData_Ch[ entriesTime-1 ];
    for( j = entriesTime-1; sPDFData_Ch[j] == xtemp; j-- ) {
       for(j++;j<entriesTime; j++) sPDFData_Ch[j] = xtemp; xtemp -= 0.1 ; 
     }
    
    for ( j = 0; j < entriesTime; j++ ) {
       sPDFX_Ch.push_back(0.25 * static_cast<Float_t> (j + 1) - 100.25);
       sPDF_Ch.push_back(sPDFData_Ch[j]);
    }
     
    for ( j = 0; j < entriesTime - 1; j++ ) sDerivative_Ch.push_back( ( sPDF_Ch[j + 1] - sPDF_Ch[j] ) * 4.0 ) ;

    entriesAngle = sPDFData_WLSAngle.size();

    for ( j = 0; j < entriesAngle; j++ ) {
       sPDFX_WLSAngle.push_back( 0.01 * static_cast<Float_t> (j + 1) - 1.01);
       sPDF_WLSAngle.push_back( sPDFData_WLSAngle[j] );
    }
     
    for ( j = 0; j < entriesAngle - 1; j++ ) sDerivative_WLSAngle.push_back( ( sPDF_WLSAngle[j + 1] - sPDF_WLSAngle[j] )*100.0 );

    for ( j = 0; j < entriesAngle; j++ ) {
       sPDFX_ChAngle.push_back( 0.01 * static_cast<Float_t> (j + 1) - 1.01);
       sPDF_ChAngle.push_back( sPDFData_ChAngle[j] );
    }
     
    for ( j = 0; j < entriesAngle - 1; j++ ) sDerivative_ChAngle.push_back( ( sPDF_ChAngle[j + 1] - sPDF_ChAngle[j] )*100.0 );
    
 }


  double FitterWLS::LAnddLdCosTheta_Ch( double cosTheta, double &dLdCosTheta_Ch)
  {
   double L_Ch, dAngle, u  ;
   
   long ibin = (int)floor((cosTheta + 1.0)*100) ;
 
  if(ibin > 0 && ibin < entriesAngle -2 )
   {
      dAngle = cosTheta - sPDFX_ChAngle[ibin] ;
      L_Ch = (sPDF_ChAngle[ibin] + sDerivative_ChAngle[ibin]*dAngle) ;

      if(dAngle< 0.005 ){
         
	  u = dAngle*100 +0.02 ;
	  dLdCosTheta_Ch = u*sDerivative_ChAngle[ibin] + (1-u)*sDerivative_ChAngle[ibin-1] ;
      }
      else{
        
	 u = 100*dAngle -0.02 ;
         dLdCosTheta_Ch = (1-u)*sDerivative_ChAngle[ibin] + u*sDerivative_ChAngle[ibin+1] ;
      }
   }
  else if(ibin<=0)
    {
   
      L_Ch = +sPDF_ChAngle[0] ;
      dLdCosTheta_Ch = sDerivative_ChAngle[0] ;
	
   }       
   else 
    {
     
       L_Ch = +sPDF_ChAngle[200-2] ;
       dLdCosTheta_Ch = sDerivative_ChAngle[200-2] ;
     }
      return L_Ch ;
  }


 double FitterWLS::LAnddLdCosTheta_WLS( double cosTheta, double &dLdCosTheta_WLS)
  {
   double L_WLS, dAngle, u  ;
   
   long ibin = (int)floor((cosTheta + 1.0)*100) ;
  
   if(ibin > 0 && ibin < entriesAngle -2 )
   {
          
      dAngle = cosTheta - sPDFX_WLSAngle[ibin] ;
      L_WLS = (sPDF_WLSAngle[ibin] + sDerivative_WLSAngle[ibin]*dAngle) ;
    
      if(dAngle< 0.005 ){
          
	u = dAngle*100 +0.02 ;
	dLdCosTheta_WLS = u*sDerivative_WLSAngle[ibin] + (1-u)*sDerivative_WLSAngle[ibin-1] ;
	
      }
      else{
    
	 u = 100*dAngle -0.02 ;
         dLdCosTheta_WLS = (1-u)*sDerivative_WLSAngle[ibin] + u*sDerivative_WLSAngle[ibin+1] ;
      }
   }
   else if(ibin<=0)
    {
          L_WLS = +sPDF_WLSAngle[0] ;
	  dLdCosTheta_WLS = sDerivative_WLSAngle[0] ;
    }       
   else {

         L_WLS = +sPDF_WLSAngle[200-2] ;
         dLdCosTheta_WLS = sDerivative_WLSAngle[200-2] ;
      
    }
 
   return L_WLS ;
  }


  
  double FitterWLS::LAnddLdt_WLS( double t, double &dLdt_WLS )
  {
    double L_WLS, dt, u;
    fRes=fT-t-sTime0;
   
    if( fRes > 3000 )fRes = 3000;
    if( fRes < -3000 )fRes = -3000;
    long ibin = (int) floor( ( fRes + 100.0 ) * 4 );
    if( ibin > 0 && ibin < entriesTime - 2 )
    {
       dt = fRes-sPDFX_WLS[ibin];
       L_WLS = ( sPDF_WLS[ibin] + sDerivative_WLS[ibin] * dt );

        if( dt < 0.125 ){
          u = 0.5+dt*4;
          dLdt_WLS = -u*sDerivative_WLS[ibin]-(1-u)*sDerivative_WLS[ibin-1];
        }
        else{
           u = 4*dt-0.5;
           dLdt_WLS = -(1-u)*sDerivative_WLS[ibin]-u*sDerivative_WLS[ibin+1];
 
        } // end of the if statement if dt<0.125
     
     } else if (ibin <= 0) {

  	  L_WLS = +sPDF_WLS[0];
          dLdt_WLS = -sDerivative_WLS[0];
     } else {

          L_WLS = +sPDF_WLS[1600-2];
          dLdt_WLS = -sDerivative_WLS[1600-2];
    }
    return L_WLS;
  }

  double FitterWLS::LAnddLdt_Ch( double t, double &dLdt_Ch )
  {
    double L_Ch, dt, u;
    fRes=fT-t-sTime0;
   
    if( fRes > 3000 )fRes = 3000;
    if( fRes < -3000 )fRes = -3000;
    long ibin = (int) floor( ( fRes + 100.0 ) * 4 );

    if( ibin > 0 && ibin < entriesTime - 2 )
     {
         dt = fRes-sPDFX_Ch[ibin];
         L_Ch = ( sPDF_Ch[ibin] + sDerivative_Ch[ibin] * dt );  
         if( dt < 0.125 ){

          u = 0.5+dt*4;
          dLdt_Ch = -u*sDerivative_Ch[ibin]-(1-u)*sDerivative_Ch[ibin-1];
         }
         else{

	   u = 4*dt-0.5;
           dLdt_Ch = -(1-u)*sDerivative_Ch[ibin]-u*sDerivative_Ch[ibin+1];
          } // end of the if statement if dt<0.125
     
     } else if (ibin <= 0) {
	      L_Ch = +sPDF_Ch[0];
	      dLdt_Ch = -sDerivative_Ch[0];
     } else {
            L_Ch = +sPDF_Ch[1600-2];	
            dLdt_Ch = -sDerivative_Ch[1600-2];
     }

    return L_Ch;
  }


 void FitterWLS::CalculatePMTLikelihood()
 {
 
  TVector3 PMTPosition, sVertex2, diffWLS ;
  TVector3 diffCh =  fPosition -sVertex ;

  double sFactor = 1/2.5 ;

  double L, dLdt_Ch, dLdt_WLS, dLdx, dLdy,dLdz, dLdt0, dLdtheta, dLdphi, dLdCosTheta_Ch, dLdCosTheta_WLS ;

  TVector3 fDirection, fDerivativeTheta, fDerivativePhi  ;
 
  fDirection.SetMagThetaPhi(1.0, sZenith, sAzimuth);

  fDerivativeTheta.SetXYZ(cos(sZenith)*cos(sAzimuth), cos(sZenith)*sin(sAzimuth), -sin(sZenith)) ;

  fDerivativePhi.SetXYZ(-sin(sZenith)*sin(sAzimuth), sin(sZenith)*cos(sAzimuth), 0.0 ) ;
  
  double cosTheta = fDirection*(diffCh.Unit()) ;

  sVertex2 = sVertex + 100*fDirection ;//Jie: remove 100 factor

  diffWLS = fPosition -sVertex2 ;

  double tCh = diffCh.Mag()/sSpeedOfLightWater ;

  double tWLS = diffWLS.Mag()/sSpeedOfLightWater ;

  dLdt_Ch =0.0, dLdt_WLS =0.0, dLdCosTheta_Ch = 0.0, dLdCosTheta_WLS=0.0  ;

  double L_Ch =0.0,  L_WLS = 0.0, L_ChAngle =0.0, L_WLSAngle = 0.0 ;
 
  L_Ch = LAnddLdt_Ch( tCh, dLdt_Ch ) ;
  L_WLS = LAnddLdt_WLS( tWLS, dLdt_WLS ) ;

  L_ChAngle = LAnddLdCosTheta_Ch( cosTheta, dLdCosTheta_Ch) ; 
  L_WLSAngle = LAnddLdCosTheta_WLS( cosTheta, dLdCosTheta_WLS) ; 
 
  L = sFactor*L_ChAngle*L_Ch + (1.0-sFactor)*L_WLSAngle*L_WLS ;

  dLdt0 = sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L ;

  dLdx = -(sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L)*(diffCh.X())*1.0/(diffCh.Mag()*sSpeedOfLightWater) ;
  dLdy = -(sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L)*(diffCh.Y())*1.0/(diffCh.Mag()*sSpeedOfLightWater) ; 
  dLdz = -(sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L)*(diffCh.Z())*1.0/(diffCh.Mag()*sSpeedOfLightWater) ;

  dLdtheta = (sFactor*dLdCosTheta_Ch*L_Ch + (1.0-sFactor)*dLdCosTheta_WLS*L_WLS)*(fDerivativeTheta*diffCh.Unit())/L ;
  dLdphi = (sFactor*dLdCosTheta_Ch*L_Ch + (1.0-sFactor)*dLdCosTheta_WLS*L_WLS)*(fDerivativePhi*diffCh.Unit())/L ;

  L = log(L) ;

  MultiPathFitter::sCovariance[0][0] += dLdx*dLdx;
  MultiPathFitter::sCovariance[1][1] += dLdy*dLdy;
  MultiPathFitter::sCovariance[2][2] += dLdz*dLdz;
  MultiPathFitter::sCovariance[3][3] += dLdt0*dLdt0;
  MultiPathFitter::sCovariance[4][4] += dLdtheta*dLdtheta ;
  MultiPathFitter::sCovariance[5][5] += dLdphi*dLdphi ;
  
  MultiPathFitter::sCovariance[1][0] += dLdy*dLdx;
  MultiPathFitter::sCovariance[2][0] += dLdz*dLdx;
  MultiPathFitter::sCovariance[3][0] += dLdt0*dLdx;
  MultiPathFitter::sCovariance[4][0] += dLdtheta*dLdx;
  MultiPathFitter::sCovariance[5][0] += dLdphi*dLdx;

  MultiPathFitter::sCovariance[2][1] += dLdz*dLdy;
  MultiPathFitter::sCovariance[3][1] += dLdt0*dLdy;
  MultiPathFitter::sCovariance[4][1] += dLdtheta*dLdy;
  MultiPathFitter::sCovariance[5][1] += dLdphi*dLdy;

  MultiPathFitter::sCovariance[3][2] += dLdt0*dLdz;
  MultiPathFitter::sCovariance[4][2] += dLdtheta*dLdz;
  MultiPathFitter::sCovariance[5][2] += dLdphi*dLdz;

  MultiPathFitter::sCovariance[4][3] += dLdtheta*dLdt0;
  MultiPathFitter::sCovariance[5][3] += dLdphi*dLdt0;

  MultiPathFitter::sCovariance[5][4] += dLdphi*dLdtheta;

  MultiPathFitter::sBeta[0][0] += dLdx;
  MultiPathFitter::sBeta[1][0] += dLdy;
  MultiPathFitter::sBeta[2][0] += dLdz;
  MultiPathFitter::sBeta[3][0] += dLdt0;
  MultiPathFitter::sBeta[4][0] += dLdtheta;
  MultiPathFitter::sBeta[5][0] += dLdphi;

  MultiPathFitter::sLikelihood += L;

 } 


} /* namespace RAT */
