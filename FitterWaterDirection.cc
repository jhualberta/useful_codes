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
#include <RAT/FitterWaterDirection.hh>
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

  double FitterWaterDirection::fZenith ;
  double FitterWaterDirection::fAzimuth ;
  double FitterWaterDirection::fFOM = 0;
  TVector3 FitterWaterDirection::sVertex ;

  std::vector<double> FitterWaterDirection::fPDFDataChAngle ;
  std::vector<double> FitterWaterDirection::fPDFChAngle ;
  std::vector<double> FitterWaterDirection::fPDFxChAngle ;
  std::vector<double> FitterWaterDirection::fDerivativeChAngle ;

  double FitterWaterDirection::fStartZenith ;
  double FitterWaterDirection::fStartAzimuth ;

  double FitterWaterDirection::fWaterRI; //index of refraction water
  double FitterWaterDirection::fSpeedOfLight;
  double FitterWaterDirection::fSpeedOfLightWater ;
  int FitterWaterDirection::fEntriesAngle ;
  double FitterWaterDirection::fAngleBinWidth;

  FitterWaterDirection::FitterWaterDirection( double aTime, const TVector3 &aPosition, const int aPMTID,  int nPar): MultiPathFitter(aTime, aPosition, aPMTID, nPar),fPosition(aPosition)
  {
    fT = aTime ;
    sTrialparArray[0] = &fZenith ;
    sTrialparArray[1] = &fAzimuth;
    CLHEP::HepRandom::setTheSeed( sEvent );
  }

  void FitterWaterDirection::SetInitParameters(size_t mcevent)
  { 
   
    fStartZenith = CLHEP::RandFlat::shoot( 0.,CLHEP::pi );
    fStartAzimuth = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );
    
    sStartparArray[0]= &fStartZenith ;
    sStartparArray[1]= &fStartAzimuth ;
 }


  void FitterWaterDirection::BeginOfRun(DS::Run& run)
  {
    DB* db = DB::Get();

    DBLinkPtr dbLink3 = db->GetLink("FIT_Ch_Angle" );
    fPDFDataChAngle = dbLink3->GetDArray( "sPDFAngle_ChFit" );
    fAngleBinWidth = dbLink3->GetD( "angle_bin_width" );
    double fAngleOffset = dbLink3->GetD( "angle_offset" ) ;
 
    TFile *f = new TFile("pdfDerivatives.root","RECREATE") ;
    TH1F *h3= new TH1F("h3","pdfAngleCh",200,1,200);
    TH1F *h3x= new TH1F("h3x","pdfXAngleCh",200,1,200);
    TH1F *dh3= new TH1F("dh3","DerivativeAngleCh",200,1,200);

    fEntriesAngle = fPDFDataChAngle.size();
    int j ;
    for (j = 0; j < fEntriesAngle; j++ ) {
      fPDFxChAngle.push_back( fAngleBinWidth * static_cast<Double_t> (j + 1) - fAngleOffset );
      fPDFChAngle.push_back( fPDFDataChAngle[j] );
    }	

   for (j = 0; j < fEntriesAngle - 1; j++ ) fDerivativeChAngle.push_back( ( fPDFChAngle[j + 1] - fPDFChAngle[j] ) * 1/fAngleBinWidth );
 
   for( j=0; j< fEntriesAngle -1 ; j++ ){
    
    h3->SetBinContent(j+1, fPDFChAngle[j]);
    h3x->SetBinContent(j+1, fPDFxChAngle[j]);
     dh3->SetBinContent(j+1, fDerivativeChAngle[j]);
   }
   f->cd();
   h3->Write(); h3x->Write(); dh3->Write();
   f->Close();
  }

  
  double FitterWaterDirection::LAnddLdCosTheta_Ch( double cosTheta, double &dLdCosTheta_Ch)
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

  void FitterWaterDirection::SetVertex(TVector3 aVertex)
 {
  sVertex = aVertex ;
 }

  void FitterWaterDirection::CalculatePMTLikelihood()
  {
    TVector3 diffCh =  fPosition -sVertex ;

    TVector3 fDirection, fDerivativeTheta, fDerivativePhi  ;
    fDirection.SetMagThetaPhi(1.0, fZenith, fAzimuth);
    fDerivativeTheta.SetXYZ(cos(fZenith)*cos(fAzimuth), cos(fZenith)*sin(fAzimuth), -sin(fZenith)) ;
    fDerivativePhi.SetXYZ(-sin(fZenith)*sin(fAzimuth), sin(fZenith)*cos(fAzimuth), 0.0 ) ;

    double cosTheta = fDirection*(diffCh.Unit()) ;
   
    double dLdCosThetaCh = 0.0 ;
    double lChAngle =0.0 ;

    lChAngle= LAnddLdCosTheta_Ch( cosTheta, dLdCosThetaCh) ; 
   
    double L = lChAngle  ;

    double dLdtheta = (dLdCosThetaCh)*(fDerivativeTheta*diffCh.Unit())/L ;
    double dLdphi = (dLdCosThetaCh)*(fDerivativePhi*diffCh.Unit())/L ;

    L = log(L) ;

    MultiPathFitter::sCovariance[0][0] += dLdtheta*dLdtheta ;
    MultiPathFitter::sCovariance[1][1] += dLdphi*dLdphi ;

    MultiPathFitter::sCovariance[1][0] += dLdtheta*dLdphi;

    MultiPathFitter::sBeta[0][0] += dLdtheta;
    MultiPathFitter::sBeta[1][0] += dLdphi;

    MultiPathFitter::sLikelihood += L;
  }


  void FitterWaterDirection::SummariseOtherOutputs(){
    fFOM=MultiPathFitter::sBestLikelihood/MultiPathFitter::sNHits;
  }


  void FitterWaterDirection::MakeParameterHistogram(){

    TString sNameHtheta = "Htheta_";
    sNameHtheta += sEvent;
    TString sNameHphi = "Htphi_";
    sNameHphi += sEvent;

    TString sNameDHtheta = "DHtheta_";
    sNameDHtheta += sEvent;
    TString sNameDHphi = "DHtphi_";
    sNameDHphi += sEvent;

    TString sNameNDHtheta = "NDHtheta_";
    sNameNDHtheta += sEvent;
    TString sNameNDHphi = "NDHtphi_";
    sNameNDHphi += sEvent;

    //histograms for likihood for position
    TH1F sHtheta= TH1F(sNameHtheta,sNameHtheta,628,-TMath::Pi(),TMath::Pi());
    TH1F sHphi= TH1F(sNameHphi,sNameHphi,942,-TMath::Pi(),2*TMath::Pi());

    //histograms for liklihod of dirivative
    TH1F sDHtheta= TH1F(sNameDHtheta,sNameDHtheta,628,-TMath::Pi(),TMath::Pi());
    TH1F sDHphi= TH1F(sNameDHphi,sNameDHphi,942,-TMath::Pi(),2*TMath::Pi());

    //histogram for likelihood of numrical derivative
    TH1F sNDHtheta = TH1F(sNameNDHtheta,sNameNDHtheta,628,-TMath::Pi(),TMath::Pi());
    TH1F sNDHphi = TH1F(sNameNDHphi,sNameNDHphi,942,-TMath::Pi(),2*TMath::Pi());

    vHists.push_back(sHtheta);
    vHists.push_back(sHphi);

    vHists.push_back(sDHtheta);
    vHists.push_back(sDHphi);

    vHists.push_back(sNDHtheta);
    vHists.push_back(sNDHphi);

  }

} /* namespace RAT */
