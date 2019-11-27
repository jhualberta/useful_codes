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
//////////////////////////////////////////////////////////////////////////

#include <RAT/MultiPathFitter.hh>
#include <RAT/FitterWater.hh>
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


namespace RAT{

  TVector3 FitterWater::fVertex; // trial vertex
  double FitterWater::fTime0; //trial fit time
  double FitterWater::fZenith ;
  double FitterWater::fAzimuth ;
  double FitterWater::fFOM = 0;

  std::vector<double> FitterWater::fPDFDataCh ;
  std::vector<double> FitterWater::fPDFxCh ;
  std::vector<double> FitterWater::fPDFCh ;
  std::vector<double> FitterWater::fDerivativeCh ;

  std::vector<double> FitterWater::fPDFDataChAngle ;
  std::vector<double> FitterWater::fPDFChAngle ;
  std::vector<double> FitterWater::fPDFxChAngle ;
  std::vector<double> FitterWater::fDerivativeChAngle ;

  TVector3 FitterWater::fStartVertex; // trial vertex
  double FitterWater::fStartTime0; //trial fit time
  double FitterWater::fStartZenith ;
  double FitterWater::fStartAzimuth ;

  double FitterWater::fWaterRI; //index of refraction water
  double FitterWater::fSpeedOfLight;
  double FitterWater::fSpeedOfLightWater ;
  int FitterWater::fEntriesAngle ;
  int FitterWater::fEntriesTimeCh ;
  double FitterWater::fChBinWidth;
  double FitterWater::fAngleBinWidth;

  FitterWater::FitterWater( double aTime, const TVector3 &aPosition, const int aPMTID,  int nPar): MultiPathFitter(aTime, aPosition, aPMTID, nPar)
  {
    fT = aTime ;
    sTrialparArray[0]= &fVertex[0] ;
    sTrialparArray[1]= &fVertex[1] ;
    sTrialparArray[2]= &fVertex[2] ;
    sTrialparArray[3] = &fTime0 ;
    sTrialparArray[4] = &fZenith ;
    sTrialparArray[5] = &fAzimuth;
    CLHEP::HepRandom::setTheSeed( sEvent );
  }

  void FitterWater::SetInitParameters(size_t mcevent)
  {
  
    double ran0 = CLHEP::RandFlat::shoot( 0.,1. );
    double ran1 = CLHEP::RandFlat::shoot( -1.,1. );
    double ran2Pi = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );

    fStartVertex.SetMagThetaPhi( pow( ran0,1.0/3 )*7000,acos( ran1 ),ran2Pi );
    //MC and data have different PMT hit time
    if(mcevent>0){
      fStartTime0 = CLHEP::RandFlat::shoot(100,300);  //MC
    }else{
      fStartTime0 = CLHEP::RandFlat::shoot(-300,-100);  //Data
    }
    fStartZenith = CLHEP::RandFlat::shoot( 0.,CLHEP::pi );
    fStartAzimuth = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );

    sStartparArray[0]= &fStartVertex[0] ;
    sStartparArray[1]= &fStartVertex[1] ;
    sStartparArray[2]= &fStartVertex[2] ;
    sStartparArray[3]= &fStartTime0 ;
    sStartparArray[4]= &fStartZenith ;
    sStartparArray[5]= &fStartAzimuth ;
 }


  void FitterWater::BeginOfRun(DS::Run& run)
  {
    DB* db = DB::Get();
    DBLinkPtr dbLink2 = db->GetLink( "FIT_FITTERAW" );
    fPDFDataCh = dbLink2->GetDArray( "sPDF_multipathfit" );
    fWaterRI = dbLink2->GetD( "water_RI" );

    fChBinWidth = dbLink2->GetD( "time_bin_width" );
    double fChOffset = dbLink2->GetD( "time_offset" );

    DBLinkPtr dbLink3 = db->GetLink("FIT_Ch_Angle" );
    fPDFDataChAngle = dbLink3->GetDArray( "sPDFAngle_ChFit" );
    fAngleBinWidth = dbLink3->GetD( "angle_bin_width" );
    double fAngleOffset = dbLink3->GetD( "angle_offset" ) ;

    fSpeedOfLight = CLHEP::c_light;
    fSpeedOfLightWater = fSpeedOfLight/fWaterRI ; 

   // std::cout<<" fSpeedOfLight "<< fSpeedOfLight<<" refractive index in water "<<fWaterRI<< " Speed of light in water "<<fSpeedOfLightWater<<std::endl ;

    TFile *f = new TFile("pdfDerivatives.root","RECREATE") ;
    TH1F *h1 = new TH1F("h1","pdfCherenkov",1600,1,1600);
    TH1F *h1x= new TH1F("h1x","pdfXCh",1600,1,1600);
    TH1F *dh1 = new TH1F("dh1","DerivativeCherenkov",1600,1,1600);

    TH1F *h3= new TH1F("h3","pdfAngleCh",200,1,200);
    TH1F *h3x= new TH1F("h3x","pdfXAngleCh",200,1,200);
    TH1F *dh3= new TH1F("dh3","DerivativeAngleCh",200,1,200);

    fEntriesTimeCh = fPDFDataCh.size();
    double xtemp = fPDFDataCh[0];
    int j ;
    for(j = fEntriesTimeCh-1; j >= 0; j-- ){
      if(fPDFDataCh[j] == fPDFDataCh[0]){
        fPDFDataCh[j] = xtemp;
        xtemp -= 0.0001;
      }
    }
  
    xtemp = fPDFDataCh[ fEntriesTimeCh-1 ];

    for( j = fEntriesTimeCh-1; fPDFDataCh[j] == xtemp; j-- ) ;
      for(j++;j<fEntriesTimeCh; j++) { 

	fPDFDataCh[j] = xtemp; 
	xtemp -= 0.0001 ; 
	}
    
  
    for (j = 0; j < fEntriesTimeCh; j++ ) {
      fPDFxCh.push_back(fChBinWidth * static_cast<Float_t> (j + 1) - fChOffset);
      fPDFCh.push_back(fPDFDataCh[j]);
    }

    for (j = 0; j < fEntriesTimeCh - 1; j++ ) fDerivativeCh.push_back( ( fPDFCh[j + 1] - fPDFCh[j] ) * 1/fChBinWidth ) ;


   for(j=0; j< fEntriesTimeCh -1 ; j++ ){
    
    h1->SetBinContent(j+1, fPDFCh[j]);
    h1x->SetBinContent(j+1, fPDFxCh[j]);
    dh1->SetBinContent(j+1, fDerivativeCh[j]);
   }

    fEntriesAngle = fPDFDataChAngle.size();

    for (j = 0; j < fEntriesAngle; j++ ) {
      fPDFxChAngle.push_back( fAngleBinWidth * static_cast<Double_t> (j + 1) - fAngleOffset );
      fPDFChAngle.push_back( fPDFDataChAngle[j] );
    }	

  //Low weight to pmts opposite to the direction of incoming particle in LLH calculations
 /*   for (j = 0; j < fEntriesAngle - 1; j++ ) {
	if (fPDFChAngle[j] >0.08) { fPDFChAngle[j] =fPDFChAngle[j]; std::cout <<" fPDFChAngle[j] "<< fPDFChAngle[j]<<std::endl ;}
        else if (fPDFChAngle[j]==0.0) { fPDFChAngle[j] = 0.0 ; std::cout <<" should be only on ebin "<<std::endl ;}
	else  { fPDFChAngle[j] = 0.08 ; std::cout <<" rest of the bins "<<std::endl ;}
    }
*/    

    for (j = 0; j < fEntriesAngle - 1; j++ ) fDerivativeChAngle.push_back( ( fPDFChAngle[j + 1] - fPDFChAngle[j] ) * 1/fAngleBinWidth );
 
   for( j=0; j< fEntriesAngle -1 ; j++ ){
    
    h3->SetBinContent(j+1, fPDFChAngle[j]);
    h3x->SetBinContent(j+1, fPDFxChAngle[j]);
     dh3->SetBinContent(j+1, fDerivativeChAngle[j]);
   }
   f->cd();
   h1->Write(); h1x->Write(); dh1->Write(); h3->Write(); h3x->Write(); dh3->Write();
   f->Close();
  }

  
  double FitterWater::LAnddLdCosTheta_Ch( double cosTheta, double &dLdCosTheta_Ch)
  {
   double L_Ch, dAngle, u  ;
   
   long ibin = (int)floor((cosTheta + 1.0)*100) ;
 
  if(ibin > 0 && ibin < fEntriesAngle -2 )
   {
      dAngle = cosTheta - fPDFxChAngle[ibin] ;
      L_Ch = (fPDFChAngle[ibin] + fDerivativeChAngle[ibin]*dAngle) ;

      if(dAngle< 0.005 ){
         
	  u = dAngle*100 +0.02 ;
	  dLdCosTheta_Ch = u*fDerivativeChAngle[ibin] + (1-u)*fDerivativeChAngle[ibin-1] ;
      }
      else{
        
	 u = 100*dAngle -0.02 ;
         dLdCosTheta_Ch = (1-u)*fDerivativeChAngle[ibin] + u*fDerivativeChAngle[ibin+1] ;
      }
   }
  else if(ibin<=0)
    {
   
      L_Ch = +fPDFChAngle[0] ;
      dLdCosTheta_Ch = fDerivativeChAngle[0] ;
	
   }       
   else 
    {
     
       L_Ch = +fPDFChAngle[200-2] ;
       dLdCosTheta_Ch = fDerivativeChAngle[200-2] ;
     }
      return L_Ch ;
  }

  double FitterWater::LAnddLdt_Ch( double t, double &dLdt_Ch )
  {
    double L_Ch, dt, u;
    fRes=fT-t-fTime0;
   
    if( fRes > 3000 )fRes = 3000;
    if( fRes < -3000 )fRes = -3000;
    long ibin = (int) floor( ( fRes + 100.0 ) * 4 );

    if( ibin > 0 && ibin < fEntriesTimeCh - 2 )
     {
         dt = fRes-fPDFxCh[ibin];
         L_Ch = ( fPDFCh[ibin] + fDerivativeCh[ibin] * dt );  
         if( dt < 0.125 ){

          u = 0.5+dt*4;
          dLdt_Ch = -u*fDerivativeCh[ibin]-(1-u)*fDerivativeCh[ibin-1];
         }
         else{

	   u = 4*dt-0.5;
           dLdt_Ch = -(1-u)*fDerivativeCh[ibin]-u*fDerivativeCh[ibin+1];
          } // end of the if statement if dt<0.125
     
     } else if (ibin <= 0) {
	      L_Ch = +fPDFCh[0];
	      dLdt_Ch = -fDerivativeCh[0];
     } else {
            L_Ch = +fPDFCh[1600-2];	
            dLdt_Ch = -fDerivativeCh[1600-2];
     }

    return L_Ch;
  }


  void FitterWater::CalculatePMTLikelihood()
  {
    TVector3 diffCh =  fPosition -fVertex ;

    TVector3 fDirection, fDerivativeTheta, fDerivativePhi  ;
    fDirection.SetMagThetaPhi(1.0, fZenith, fAzimuth);
    fDerivativeTheta.SetXYZ(cos(fZenith)*cos(fAzimuth), cos(fZenith)*sin(fAzimuth), -sin(fZenith)) ;
    fDerivativePhi.SetXYZ(-sin(fZenith)*sin(fAzimuth), sin(fZenith)*cos(fAzimuth), 0.0 ) ;

    double cosTheta = fDirection*(diffCh.Unit()) ;
   
    double tCh = diffCh.Mag()/fSpeedOfLightWater ;
  
    double dLdtCh = 0.0, dLdCosThetaCh = 0.0 ;
    double lCh =0.0, lChAngle =0.0 ;

    lCh = LAnddLdt_Ch( tCh, dLdtCh ) ;
    
    lChAngle= LAnddLdCosTheta_Ch( cosTheta, dLdCosThetaCh) ; 
   
    double L = lChAngle*lCh  ;

    double dLdt0 = lChAngle*dLdtCh/L  ;
    double dLdx = -(lChAngle*dLdtCh/L)*(diffCh.X())*1.0/(diffCh.Mag()*fSpeedOfLightWater) ;
    double dLdy = -(lChAngle*dLdtCh/L)*(diffCh.Y())*1.0/(diffCh.Mag()*fSpeedOfLightWater) ;
    double dLdz = -(lChAngle*dLdtCh/L)*(diffCh.Z())*1.0/(diffCh.Mag()*fSpeedOfLightWater) ;
    double dLdtheta = (dLdCosThetaCh*lCh)*(fDerivativeTheta*diffCh.Unit())/L ;
    double dLdphi = (dLdCosThetaCh*lCh)*(fDerivativePhi*diffCh.Unit())/L ;

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


  void FitterWater::SummariseOtherOutputs(){
    fFOM=MultiPathFitter::sBestLikelihood/MultiPathFitter::sNHits;
  }


  void FitterWater::MakeParameterHistogram(){

    TString sNameHx = "Hx_";
    sNameHx += sEvent;
    TString sNameHy = "Hy_";
    sNameHy += sEvent;
    TString sNameHz = "Hz_";
    sNameHz += sEvent;
    TString sNameHt = "Ht_";
    sNameHt += sEvent;
    TString sNameHtheta = "Htheta_";
    sNameHtheta += sEvent;
    TString sNameHphi = "Htphi_";
    sNameHphi += sEvent;

    TString sNameDHx = "DHx_";
    sNameDHx += sEvent;
    TString sNameDHy = "DHy_";
    sNameDHy += sEvent;
    TString sNameDHz = "DHz_";
    sNameDHz += sEvent;
    TString sNameDHt = "DHt_";
    sNameDHt += sEvent;
    TString sNameDHtheta = "DHtheta_";
    sNameDHtheta += sEvent;
    TString sNameDHphi = "DHtphi_";
    sNameDHphi += sEvent;

    TString sNameNDHx = "NDHx_";
    sNameNDHx += sEvent;
    TString sNameNDHy = "NDHy_";
    sNameNDHy += sEvent;
    TString sNameNDHz = "NDHz_";
    sNameNDHz += sEvent;
    TString sNameNDHt = "NDHt_";
    sNameNDHt += sEvent;
    TString sNameNDHtheta = "NDHtheta_";
    sNameNDHtheta += sEvent;
    TString sNameNDHphi = "NDHtphi_";
    sNameNDHphi += sEvent;

    //histograms for likihood for position
    TH1F sHx= TH1F(sNameHx,sNameHx,1000,-5000,5000);
    TH1F sHy= TH1F(sNameHy,sNameHy,1000,-5000,5000);
    TH1F sHz= TH1F(sNameHz,sNameHz,1000,-5000,5000);
    TH1F sHt= TH1F(sNameHt,sNameHt,400,-100,300);
    TH1F sHtheta= TH1F(sNameHtheta,sNameHtheta,628,-TMath::Pi(),TMath::Pi());
    TH1F sHphi= TH1F(sNameHphi,sNameHphi,942,-TMath::Pi(),2*TMath::Pi());

    //histograms for liklihod of dirivative
    TH1F sDHx= TH1F(sNameDHx,sNameDHx,1000,-5000,5000);
    TH1F sDHy= TH1F(sNameDHy,sNameDHy,1000,-5000,5000);
    TH1F sDHz= TH1F(sNameDHz,sNameDHz,1000,-5000,5000);
    TH1F sDHt= TH1F(sNameDHt,sNameDHt,400,-100,300);
    TH1F sDHtheta= TH1F(sNameDHtheta,sNameDHtheta,628,-TMath::Pi(),TMath::Pi());
    TH1F sDHphi= TH1F(sNameDHphi,sNameDHphi,942,-TMath::Pi(),2*TMath::Pi());

    //histogram for likelihood of numrical derivative
    TH1F sNDHx= TH1F(sNameNDHx,sNameNDHx,1000,-5000,5000);
    TH1F sNDHy= TH1F(sNameNDHy,sNameNDHy,1000,-5000,5000);
    TH1F sNDHz= TH1F(sNameNDHz,sNameNDHz,1000,-5000,5000);
    TH1F sNDHt= TH1F(sNameNDHt,sNameNDHt,400,-100,300);
    TH1F sNDHtheta = TH1F(sNameNDHtheta,sNameNDHtheta,628,-TMath::Pi(),TMath::Pi());
    TH1F sNDHphi = TH1F(sNameNDHphi,sNameNDHphi,942,-TMath::Pi(),2*TMath::Pi());

    vHists.push_back(sHx);
    vHists.push_back(sHy);
    vHists.push_back(sHz);
    vHists.push_back(sHt);
    vHists.push_back(sHtheta);
    vHists.push_back(sHphi);

    vHists.push_back(sDHx);
    vHists.push_back(sDHy);
    vHists.push_back(sDHz);
    vHists.push_back(sDHt);
    vHists.push_back(sDHtheta);
    vHists.push_back(sDHphi);

    vHists.push_back(sNDHx);
    vHists.push_back(sNDHy);
    vHists.push_back(sNDHz);
    vHists.push_back(sNDHt);
    vHists.push_back(sNDHtheta);
    vHists.push_back(sNDHphi);

  }

} /* namespace RAT */
