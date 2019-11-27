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
#include <RAT/FitterWaterPosition.hh>
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

  TVector3 FitterWaterPosition::fVertex; // trial vertex
  double FitterWaterPosition::fTime0; //trial fit time
  double FitterWaterPosition::fFOM = 0;

  std::vector<double> FitterWaterPosition::fPDFDataCh ;
  std::vector<double> FitterWaterPosition::fPDFxCh ;
  std::vector<double> FitterWaterPosition::fPDFCh ;
  std::vector<double> FitterWaterPosition::fDerivativeCh ;


  TVector3 FitterWaterPosition::fStartVertex; // trial vertex
  double FitterWaterPosition::fStartTime0; //trial fit time

  double FitterWaterPosition::fWaterRI; //index of refraction water
  double FitterWaterPosition::fSpeedOfLight;
  double FitterWaterPosition::fSpeedOfLightWater ;
  int FitterWaterPosition::fEntriesTimeCh ;
  double FitterWaterPosition::fChBinWidth;


  FitterWaterPosition::FitterWaterPosition( double aTime, const TVector3 &aPosition, const int aPMTID,  int nPar): MultiPathFitter(aTime, aPosition, aPMTID, nPar),fPosition(aPosition)
  {
    fT = aTime ;
    sTrialparArray[0]= &fVertex[0] ;
    sTrialparArray[1]= &fVertex[1] ;
    sTrialparArray[2]= &fVertex[2] ;
    sTrialparArray[3] = &fTime0 ;
    CLHEP::HepRandom::setTheSeed( sEvent );
  }

  void FitterWaterPosition::SetInitParameters(size_t mcevent)
  {
 
    double ran0 = CLHEP::RandFlat::shoot( 0.,1. );
    double ran1 = CLHEP::RandFlat::shoot( -1.,1. );
    double ran2Pi = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );
    fStartVertex.SetMagThetaPhi( pow( ran0,1.0/3 )*7000,acos( ran1 ),ran2Pi );
    //MC and data have different PMT hit time
    if(mcevent>0){
      fStartTime0 = CLHEP::RandFlat::shoot(100,300);  //MC
    }else{
      fStartTime0 = CLHEP::RandFlat::shoot(-100,100);  //Data
    }
    sStartparArray[0]= &fStartVertex[0] ;
    sStartparArray[1]= &fStartVertex[1] ;
    sStartparArray[2]= &fStartVertex[2] ;
    sStartparArray[3]= &fStartTime0 ;
 }


  void FitterWaterPosition::BeginOfRun(DS::Run& run)
  {
    DB* db = DB::Get();
    DBLinkPtr dbLink2 = db->GetLink( "FIT_FITTERAW" );
    fPDFDataCh = dbLink2->GetDArray( "sPDF_multipathfit" );
    fWaterRI = dbLink2->GetD( "water_RI" );

    fChBinWidth = dbLink2->GetD( "time_bin_width" );
    double fChOffset = dbLink2->GetD( "time_offset" );

    fSpeedOfLight = CLHEP::c_light;
    fSpeedOfLightWater = fSpeedOfLight/fWaterRI ; 

    TFile *f = new TFile("pdfDerivatives_Position.root","RECREATE") ;
    TH1F *h1 = new TH1F("h1","pdfCherenkov",1600,1,1600);
    TH1F *h1x= new TH1F("h1x","pdfXCh",1600,1,1600);
    TH1F *dh1 = new TH1F("dh1","DerivativeCherenkov",1600,1,1600);

    fEntriesTimeCh = fPDFDataCh.size();
    double xtemp = fPDFDataCh[0];
    int j ;
    for(j = fEntriesTimeCh-1; j >= 0; j-- ){
      if(fPDFDataCh[j] == fPDFDataCh[0]){
        fPDFDataCh[j] = xtemp;
        xtemp -= 0.001;
      }
    }
  
    xtemp = fPDFDataCh[ fEntriesTimeCh-1 ];

    for( j = fEntriesTimeCh-1; fPDFDataCh[j] == xtemp; j-- ) ;
      for(j++;j<fEntriesTimeCh; j++) { 

	fPDFDataCh[j] = xtemp; 
	xtemp -= 0.001; //0.0000001 ; 
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

   f->cd();
   h1->Write(); h1x->Write(); dh1->Write(); f->Close();
  }


  double FitterWaterPosition::LAnddLdt_Ch( double t, double &dLdt_Ch )
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


  void FitterWaterPosition::CalculatePMTLikelihood()
  {
    TVector3 diffCh =  fPosition -fVertex ;

    double tCh = diffCh.Mag()/fSpeedOfLightWater ;
    double dLdtCh = 0.0 ;
    double lCh =0.0 ;

    lCh = LAnddLdt_Ch( tCh, dLdtCh ) ;
  
    double L = lCh  ;

    double dLdt0 = dLdtCh/L  ;
    double dLdx = -(dLdtCh/L)*(diffCh.X())*1.0/(diffCh.Mag()*fSpeedOfLightWater) ;
    double dLdy = -(dLdtCh/L)*(diffCh.Y())*1.0/(diffCh.Mag()*fSpeedOfLightWater) ;
    double dLdz = -(dLdtCh/L)*(diffCh.Z())*1.0/(diffCh.Mag()*fSpeedOfLightWater) ;

    L = log(L) ;

    MultiPathFitter::sCovariance[0][0] += dLdx*dLdx;
    MultiPathFitter::sCovariance[1][1] += dLdy*dLdy;
    MultiPathFitter::sCovariance[2][2] += dLdz*dLdz;
    MultiPathFitter::sCovariance[3][3] += dLdt0*dLdt0;

    MultiPathFitter::sCovariance[1][0] += dLdy*dLdx;
    MultiPathFitter::sCovariance[2][0] += dLdz*dLdx;
    MultiPathFitter::sCovariance[3][0] += dLdt0*dLdx;

    MultiPathFitter::sCovariance[2][1] += dLdz*dLdy;
    MultiPathFitter::sCovariance[3][1] += dLdt0*dLdy;

    MultiPathFitter::sCovariance[3][2] += dLdt0*dLdz;

    MultiPathFitter::sBeta[0][0] += dLdx;
    MultiPathFitter::sBeta[1][0] += dLdy;
    MultiPathFitter::sBeta[2][0] += dLdz;
    MultiPathFitter::sBeta[3][0] += dLdt0;
    MultiPathFitter::sLikelihood += L;
  }


  void FitterWaterPosition::SummariseOtherOutputs(){
    fFOM=MultiPathFitter::sBestLikelihood/MultiPathFitter::sNHits;
  }


  void FitterWaterPosition::MakeParameterHistogram(){

    TString sNameHx = "Hx_";
    sNameHx += sEvent;
    TString sNameHy = "Hy_";
    sNameHy += sEvent;
    TString sNameHz = "Hz_";
    sNameHz += sEvent;
    TString sNameHt = "Ht_";
    sNameHt += sEvent;

    TString sNameDHx = "DHx_";
    sNameDHx += sEvent;
    TString sNameDHy = "DHy_";
    sNameDHy += sEvent;
    TString sNameDHz = "DHz_";
    sNameDHz += sEvent;
    TString sNameDHt = "DHt_";
    sNameDHt += sEvent;

    TString sNameNDHx = "NDHx_";
    sNameNDHx += sEvent;
    TString sNameNDHy = "NDHy_";
    sNameNDHy += sEvent;
    TString sNameNDHz = "NDHz_";
    sNameNDHz += sEvent;
    TString sNameNDHt = "NDHt_";
    sNameNDHt += sEvent;

    //histograms for likihood for position
    TH1F sHx= TH1F(sNameHx,sNameHx,1000,-5000,5000);
    TH1F sHy= TH1F(sNameHy,sNameHy,1000,-5000,5000);
    TH1F sHz= TH1F(sNameHz,sNameHz,1000,-5000,5000);
    TH1F sHt= TH1F(sNameHt,sNameHt,400,-100,300);

    //histograms for liklihod of dirivative
    TH1F sDHx= TH1F(sNameDHx,sNameDHx,1000,-5000,5000);
    TH1F sDHy= TH1F(sNameDHy,sNameDHy,1000,-5000,5000);
    TH1F sDHz= TH1F(sNameDHz,sNameDHz,1000,-5000,5000);
    TH1F sDHt= TH1F(sNameDHt,sNameDHt,400,-100,300);

    //histogram for likelihood of numrical derivative
    TH1F sNDHx= TH1F(sNameNDHx,sNameNDHx,1000,-5000,5000);
    TH1F sNDHy= TH1F(sNameNDHy,sNameNDHy,1000,-5000,5000);
    TH1F sNDHz= TH1F(sNameNDHz,sNameNDHz,1000,-5000,5000);
    TH1F sNDHt= TH1F(sNameNDHt,sNameNDHt,400,-100,300);

    vHists.push_back(sHx);
    vHists.push_back(sHy);
    vHists.push_back(sHz);
    vHists.push_back(sHt);

    vHists.push_back(sDHx);
    vHists.push_back(sDHy);
    vHists.push_back(sDHz);
    vHists.push_back(sDHt);

    vHists.push_back(sNDHx);
    vHists.push_back(sNDHy);
    vHists.push_back(sNDHz);
    vHists.push_back(sNDHt);

  }

} /* namespace RAT */
