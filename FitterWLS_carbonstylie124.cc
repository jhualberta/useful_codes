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

  sStartVertex.SetXYZ(0.0,0.0,0.0);
  sStartZenith = TMath::Pi()/2.0 ;
  sStartAzimuth = 0.0 ;

  sStartparArray[0]= &sStartVertex[0] ;
  sStartparArray[1]= &sStartVertex[1] ;
  sStartparArray[2]= &sStartVertex[2] ;
  sStartparArray[3]= &sStartTime0 ;
  sStartparArray[4]= &sStartZenith ;
  sStartparArray[5]= &sStartAzimuth ;

 }

 void FitterWLS::Initialise() 
 {
    DB* db = DB::Get();
    DBLinkPtr dbLink = db->GetLink( "FIT_WLS_Time" );
    sPDFData_WLS = dbLink->GetDArray( "sPDFTime_WLSFit" );

    DBLinkPtr dbLink1 = db->GetLink( "FIT_WLS_Angle" );
    sPDFData_WLSAngle = dbLink1->GetDArray("sPDFAngle_WLSFit");

    DBLinkPtr effDB = db->GetLink( "EFFECTIVE_VELOCITY" );
    effDB = db->GetLink( "EFFECTIVE_VELOCITY", "lightwater_sno" );
    sSpeedOfLightWater = effDB->GetD( "inner_av_velocity" );    

    DBLinkPtr dbLink2 = db->GetLink( "FIT_MULTIPATH" );
    sPDFData_Ch = dbLink2->GetDArray( "sPDF_multipathfit" );

    DBLinkPtr dbLink3 = db->GetLink("FIT_Ch_Angle");
    sPDFData_ChAngle = dbLink3->GetDArray("sPDFAngle_ChFit");

    TFile *f = new TFile("pdfDerivatives.root","RECREATE") ;
    TH1F *h1 = new TH1F("h1","pdfCherenkov",1600,1,1600);
    TH1F *h1x= new TH1F("h1x","pdfXCh",1600,1,1600);
    TH1F *dh1 = new TH1F("dh1","DerivativeCherenkov",1600,1,1600);
    TH1F *h2 = new TH1F("h2","pdfWLS",1600,1,1600);
    TH1F *h2x= new TH1F("h2x","pdfXWLS",1600,1,1600);
    TH1F *dh2 = new TH1F("dh2","DerivativeWLS",1600,1,1600);
    TH1F *h3= new TH1F("h3","pdfAngleWLS",200,1,200);
    TH1F *h3x= new TH1F("h3x","pdfXAngleWLS",200,1,200);
    TH1F *dh3= new TH1F("dh3","DerivativeAngleWLS",200,1,200);
    TH1F *h4= new TH1F("h4","pdfAngleCh",200,1,200);
    TH1F *h4x= new TH1F("h4x","pdfXAngleCh",200,1,200);
    TH1F *dh4= new TH1F("dh4","DerivativeAngleCh",200,1,200);

    int j;
    entriesTime = sPDFData_WLS.size(); 
    //  Add a very slight slope to the original and final bins to get
    //  rid of singular matrices
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
   //sTolerance = 0.001;
   //  Add a very slight slope to the original and final bins to get
   //  rid of singular matrices
    xtemp = sPDFData_Ch[0];
    for( j = entriesTime-1; j >= 0; j-- ){
        if(sPDFData_Ch[j] == sPDFData_Ch[0]){
          sPDFData_Ch[j] = xtemp;
          xtemp -= 0.1;
      }
    }
      
    xtemp = sPDFData_Ch[ entriesTime-1 ];
    for( j = entriesTime-1; sPDFData_Ch[j] == xtemp; j-- );
    for(j++;j<entriesTime; j++){
        sPDFData_Ch[j] = xtemp; xtemp -= 0.1;
    }

    for ( j = 0; j < entriesTime; j++ ) {
       sPDFX_Ch.push_back(0.25 * static_cast<Float_t> (j + 1) - 100.25);
       sPDF_Ch.push_back(sPDFData_Ch[j]);
    }

    for ( j = 0; j < entriesTime - 1; j++ ) {
       sDerivative_Ch.push_back( ( sPDF_Ch[j + 1] - sPDF_Ch[j] ) * 4.0 );
    }

   for(j=0; j< entriesTime -1 ; j++ ){
     h1->SetBinContent(j+1,sPDF_Ch[j]);
     h1x->SetBinContent(j+1,sPDFX_Ch[j]);
     dh1->SetBinContent(j+1,sDerivative_Ch[j]);
     h2->SetBinContent(j+1,sPDF_WLS[j]);
     h2x->SetBinContent(j+1,sPDFX_WLS[j]);
     dh2->SetBinContent(j+1,sDerivative_WLS[j]);
   }

   entriesAngle = sPDFData_WLSAngle.size();
  
   //sTolerance = 0.001;
   //  Add a very slight slope to the original and final bins to get
   //  rid of singular matrices
   /*    xtemp = sPDFData_WLSAngle[0];
    for( j = entriesAngle-1; j >= 0; j-- ){
        if(sPDFData_WLSAngle[j] == sPDFData_WLSAngle[0]){
          sPDFData_WLSAngle[j] = xtemp;
	  xtemp -= 0.01;
        }
   }

    xtemp = sPDFData_WLSAngle[ entriesAngle-1 ];
    for( j = entriesAngle-1; sPDFData_WLSAngle[j] == xtemp; j-- );
    for(j++; j<entriesAngle ;j++){
        sPDFData_WLSAngle[j] = xtemp; xtemp -= 0.01;
   }
   */

    for ( j = 0; j < entriesAngle; j++ ) {
       sPDFX_WLSAngle.push_back( 0.01 * static_cast<Float_t> (j + 1) - 1.01);
       sPDF_WLSAngle.push_back( sPDFData_WLSAngle[j] );
    }

    for ( j = 0; j < entriesAngle - 1; j++ ) {
       sDerivative_WLSAngle.push_back( ( sPDF_WLSAngle[j + 1] - sPDF_WLSAngle[j] )*100.0 );
    }


   for(j=0; j< entriesAngle -1 ; j++ ){
     h3->SetBinContent(j+1, sPDF_WLSAngle[j]);
     h3x->SetBinContent(j+1, sPDFX_WLSAngle[j]);
     dh3->SetBinContent(j+1, sDerivative_WLSAngle[j]);
   }

    for ( j = 0; j < entriesAngle; j++ ) {
       sPDFX_ChAngle.push_back( 0.01 * static_cast<Float_t> (j + 1) - 1.01);
       sPDF_ChAngle.push_back( sPDFData_ChAngle[j] );
    }

    for ( j = 0; j < entriesAngle - 1; j++ ) {
       sDerivative_ChAngle.push_back( ( sPDF_ChAngle[j + 1] - sPDF_ChAngle[j] )*100.0 );
    }


   for(j=0; j< entriesAngle -1 ; j++ ){
     h4->SetBinContent(j+1, sPDF_ChAngle[j]);
     h4x->SetBinContent(j+1, sPDFX_ChAngle[j]);
     dh4->SetBinContent(j+1, sDerivative_ChAngle[j]);
   }

   f->cd();
   h1->Write(); h1x->Write(); dh1->Write(); h2->Write(); h2x->Write(); dh2->Write();h3->Write(); h3x->Write(); dh3->Write(); h4->Write(); h4x->Write(); dh4->Write() ;
   f->Close();	
 }


   // returns Likelihood and derivative for a given time of flight t.

  double FitterWLS::LAnddLdCosTheta( double cosTheta, double &dLdCosTheta, std::string pdfType)
  {
   double L, dAngle, u  ;
   
   long ibin = (int)floor((cosTheta + 1.0)*100) ;
   // std::cout<<" ibin in angle pdf "<< ibin << std::endl ;

   if(ibin > 0 && ibin < entriesAngle -2 )
   {
     if(pdfType =="Cherenkov") {
       dAngle = cosTheta - sPDFX_ChAngle[ibin] ;
       L = (sPDF_ChAngle[ibin] + sDerivative_ChAngle[ibin]*dAngle) ;
     }
     else{
       dAngle = cosTheta - sPDFX_WLSAngle[ibin] ;
       L = (sPDF_WLSAngle[ibin] + sDerivative_WLSAngle[ibin]*dAngle) ;
     }
     // if(dAngle< 0.005) {
     if(dAngle< 0.005 ){

        if(pdfType =="Cherenkov") {
	  u = dAngle*100 +0.02 ;
	  dLdCosTheta = u*sDerivative_ChAngle[ibin] + (1-u)*sDerivative_ChAngle[ibin-1] ;
	}
	else{
	  u = dAngle*100 +0.02 ;
	  dLdCosTheta = u*sDerivative_WLSAngle[ibin] + (1-u)*sDerivative_WLSAngle[ibin-1] ;
	}
     }
     else{
       if(pdfType =="Cherenkov") {
	 u = 100*dAngle -0.02 ;
         dLdCosTheta = (1-u)*sDerivative_ChAngle[ibin] + u*sDerivative_ChAngle[ibin+1] ;
       }
       else{
	 u = 100*dAngle -0.02 ;
         dLdCosTheta = (1-u)*sDerivative_WLSAngle[ibin] + u*sDerivative_WLSAngle[ibin+1] ;
       }
     }
   }
   else if(ibin<=0)
	{
	  if(pdfType =="Cherenkov") {
	    L = +sPDF_ChAngle[0] ;
	    dLdCosTheta = sDerivative_ChAngle[0] ;
	  }
	  else{
	    L = +sPDF_WLSAngle[0] ;
	    dLdCosTheta = sDerivative_WLSAngle[0] ;
	  }
	}       
   else {
     if(pdfType =="Cherenkov") {
       L = +sPDF_ChAngle[200-2] ;
       dLdCosTheta = sDerivative_ChAngle[200-2] ;
     }
     else{
       L = +sPDF_WLSAngle[200-2] ;
       dLdCosTheta = sDerivative_WLSAngle[200-2] ;
     }
   }
   //   std::cout<< " angle L value "<< L << " and its derivative "<< dLdCosTheta <<std::endl ;
   return L ;
  }
  

  double FitterWLS::LAnddLdt( double t, double &dLdt, std::string pdfType )
  {
    double L, dt, u;
    fRes=fT-t-sTime0;
    //std::cout << " hit Time "<< fT<< " TOF "<< t<< " trialTime " << sTrialTime0<< " fResidual time "<< fRes <<std::endl ;
    if( fRes > 3000 )fRes = 3000;
    if( fRes < -3000 )fRes = -3000;
    long ibin = (int) floor( ( fRes + 100.0 ) * 4 );
    if( ibin > 0 && ibin < entriesTime - 2 ) {
      if(pdfType =="Cherenkov") {
           dt = fRes-sPDFX_Ch[ibin];
            L = ( sPDF_Ch[ibin] + sDerivative_Ch[ibin] * dt );
	}
      else{
          dt = fRes-sPDFX_WLS[ibin];
          L = ( sPDF_WLS[ibin] + sDerivative_WLS[ibin] * dt );
      }
      if( dt < 0.125 ){
        if(pdfType =="Cherenkov"){
          u = 0.5+dt*4;
          dLdt = -u*sDerivative_Ch[ibin]-(1-u)*sDerivative_Ch[ibin-1];
         }
        else{
           u = 0.5+dt*4;
          dLdt = -u*sDerivative_WLS[ibin]-(1-u)*sDerivative_WLS[ibin-1];
         }
      }else{
        if(pdfType =="Cherenkov") {
	   u = 4*dt-0.5;
           dLdt = -(1-u)*sDerivative_Ch[ibin]-u*sDerivative_Ch[ibin+1];
         }
	else{
           u = 4*dt-0.5;
           dLdt = -(1-u)*sDerivative_WLS[ibin]-u*sDerivative_WLS[ibin+1];
        }
      }
    } else if (ibin <= 0) {
  	if(pdfType == "Cherenkov") {
	      L = +sPDF_Ch[0];
	      dLdt = -sDerivative_Ch[0];
        }
 	else{
              L = +sPDF_WLS[0];
              dLdt = -sDerivative_WLS[0];
	}	
    } else {
      if(pdfType == "Cherenkov") {
            L = +sPDF_Ch[1600-2];	
            dLdt = -sDerivative_Ch[1600-2];
	}
      else{
            L = +sPDF_WLS[1600-2];
            dLdt = -sDerivative_WLS[1600-2];

	}
    }
    return L;
  }


 void FitterWLS::CalculatePMTLikelihood()
 {
 
  TVector3 PMTPosition ;
  TVector3 diff=  fPosition -sVertex ;
//  std::cout<<" sVertex "<<sVertex.X()<<","<<sVertex.Y()<<","<<sVertex.Z()<<" PMT "<<fPosition.X()<<","<<fPosition.Y()<<","<<fPosition.Z()<<" diff "<<diff.X()<<","<<diff.Y()<<","<<diff.Z()<<std::endl ;

  double sFactor = 1/3.0 ;

  double t = diff.Mag()/sSpeedOfLightWater ;

  double L, dLdt_Ch, dLdt_WLS, dLdx, dLdy,dLdz, dLdt0, dLdtheta, dLdphi, dLdCosTheta_Ch, dLdCosTheta_WLS ;

  TVector3 fDirection, fDerivativeTheta, fDerivativePhi  ;
 
  fDirection.SetMagThetaPhi(1.0, sZenith, sAzimuth);

  fDerivativeTheta.SetXYZ(cos(sZenith)*cos(sAzimuth), cos(sZenith)*sin(sAzimuth), -sin(sZenith)) ;

  fDerivativePhi.SetXYZ(-sin(sZenith)*sin(sAzimuth), sin(sZenith)*cos(sAzimuth), 0.0 ) ;
  
  double cosTheta = fDirection*(diff.Unit()) ;

  dLdt_Ch =0.0, dLdt_WLS =0.0, dLdCosTheta_Ch = 0.0, dLdCosTheta_WLS=0.0  ;

  double L_Ch =0.0,  L_WLS = 0.0, L_ChAngle =0.0, L_WLSAngle = 0.0 ;
 
  L_Ch = LAnddLdt( t, dLdt_Ch, "Cherenkov" ) ;
  L_WLS = LAnddLdt( t, dLdt_WLS, "WLS" ) ;

  L_ChAngle = LAnddLdCosTheta( cosTheta, dLdCosTheta_Ch, "Cherenkov") ; 
  L_WLSAngle = LAnddLdCosTheta( cosTheta, dLdCosTheta_WLS, "WLS") ; 

  //if(L_Angle< 9.0) sFactor =0.0  ;
  //  L = sFactor*L_Angle*L_Ch + (1.0-sFactor)*L_WLS ;
 
  L = sFactor*L_ChAngle*L_Ch + (1.0-sFactor)*L_WLSAngle*L_WLS ;

  //std::cout<<" cosTheta "<< cosTheta<<" L_Ch "<<L_Ch<<" L_WLS "<<L_WLS<<" L_ChAngle "<<L_ChAngle<<" L_WLSAngle "<< L_WLSAngle<<" L " << L << std::endl ;

  dLdt0 = sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L ;

  dLdx = -(sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L)*(diff.X())*1.0/(diff.Mag()*sSpeedOfLightWater) ;
  dLdy = -(sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L)*(diff.Y())*1.0/(diff.Mag()*sSpeedOfLightWater) ; 
  dLdz = -(sFactor*L_ChAngle*dLdt_Ch/L + (1.0-sFactor)*L_WLSAngle*dLdt_WLS/L)*(diff.Z())*1.0/(diff.Mag()*sSpeedOfLightWater) ;

  dLdtheta = (sFactor*dLdCosTheta_Ch*L_Ch + (1.0-sFactor)*dLdCosTheta_WLS*L_WLS)*(fDerivativeTheta*diff.Unit())/L ;
  dLdphi = (sFactor*dLdCosTheta_Ch*L_Ch + (1.0-sFactor)*dLdCosTheta_WLS*L_WLS)*(fDerivativePhi*diff.Unit())/L ;

    // dLdtheta = (sFactor*dLdCosTheta*(fDerivativeTheta*diff.Unit())*L_Ch)/L ;
    //dLdphi =  (sFactor*dLdCosTheta*(fDerivativePhi*diff.Unit())*L_Ch)/L ;

  L = log(L) ;

  //std::cout<<" sVertex "<<sVertex.X()<<","<<sVertex.Y()<<","<<sVertex.Z()<<" PMT "<<fPosition.X()<<","<<fPosition.Y()<<","<<fPosition.Z()<<" diff "<<diff.X()<<","<<diff.Y()<<","<<diff.Z()<< " Direction "<<fDirection.X()<<" ,"<<fDirection.Y()<<", "<<fDirection.Z()<<" cosTheta "<< cosTheta<< " L_Angle "<<L_Angle<<" L_Ch "<<L_Ch<< " L_WLS "<<L_WLS<< " sum "<< sFactor*L_Angle*L_Ch + (1.0-sFactor)*L_WLS  << " and Log L " << L<< std::endl ;


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
