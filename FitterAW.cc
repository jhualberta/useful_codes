////////////////////////////////////////////////////////////////////////
///
/// \class
///
/// \brief
///
/// \author David Auty auty@ualberta.ca
///
/// REVISION HISTORY:\n
///
/// \detail
///       To fit data in a partial fill geometry
////
////////////////////////////////////////////////////////////////////////


#include <RAT/MultiPathFitter.hh>
#include <RAT/FitterAW.hh>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/Log.hh>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <vector>
#include <TF1.h>
#include <TH1F.h>
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"

namespace RAT{

  TVector3 FitterAW::fVertex; // trial vertex
  double FitterAW::fTime0; //trial fit time
  TVector3 FitterAW::fBestFitDirection;
  double FitterAW::fFOM = 0;
  TVector3 FitterAW::fStartVertex; // random start vertex
  double FitterAW::fStartTime0; //random start time
  double FitterAW::fWaterLevel; //distance of water level from 0 z
  double FitterAW::fAVRadius;//radius of the av
  double FitterAW::fAirRI; //index of refraction air
  double FitterAW::fWaterRI; //index of refraction water
  double FitterAW::fMaxRes; // maximum allowed resdiual
  int FitterAW::fMaximumReflections;
  int FitterAW::fStepsBelowTolerance;
  std::vector<double> FitterAW::fPDFData; //hold the PDF data
  std::vector<double> FitterAW::fPDFX;
  double FitterAW::fSpeedOfLight;
  double FitterAW::fSpeedOfLightAir; //(mm/s)
  unsigned int FitterAW::fNumFitFactors;
  unsigned int FitterAW::fGoodFitThreshold;
  unsigned int FitterAW::fMaxNumStartPos; //Maximum number of start positions to be tested
  double FitterAW::fRhoCut; //Limit in rho to generate cnadidate vertexes
  double FitterAW::fBinWidth; //time resolution in ns
  double FitterAW::fTimeOffset;
  double FitterAW::fPSupRad; //radius of PSUP
  int FitterAW::fPDFentries; //number of fPDFentries in pdf
  TMatrixT<double> FitterAW::fFactor(52, 52 );
  std::vector<double> FitterAW::fRWater;
  std::vector<double> FitterAW::fRAir;
  std::vector<double> FitterAW::fDerivative;
  double FitterAW::fDepth; //depth of vertex below water
  double FitterAW::fFitTolerance; //chisquare tolerance
  TF1 *FitterAW::fX=new TF1( "Xfun",FitterAW::X, 0,0.850908,1 );

  FitterAW::FitterAW( double aTime, const TVector3 &aPosition, const int aPMTID,  int nPar, MultipathFitterClassData &data):
		  MultiPathFitter(aTime, aPosition, aPMTID, nPar),
		  fShared(data)
  {
    fT = aTime;
    sTrialparArray[0] = &fShared.fVertex[0];
    sTrialparArray[1] = &fVertex[1];
    sTrialparArray[2] = &fVertex[2];
    sTrialparArray[3] = &fTime0;
    CLHEP::HepRandom::setTheSeed( sEvent );
    if( fWaterLevel >= -fAVRadius && fWaterLevel <= fAVRadius ) {
      fHeight = aPosition.Z()-fWaterLevel;
    } else { fHeight = aPosition.Z(); }
  }


  void FitterAW::SetInitParameters(size_t mcevent)
  {
    double ran0 = CLHEP::RandFlat::shoot( 0.,1. );
    double ran1 = CLHEP::RandFlat::shoot( -1.,1. );
    double ran2Pi = CLHEP::RandFlat::shoot( 0.,2*CLHEP::pi );

    fStartVertex.SetMagThetaPhi( pow( ran0,1.0/3. )*7000,acos( ran1 ),ran2Pi );
    //MC and data have different PMT hit time
    if( mcevent>0 ) {
      fStartTime0 = CLHEP::RandFlat::shoot(100,300); //MC
    } else {
      fStartTime0 = CLHEP::RandFlat::shoot(-300,-100); //Data
    }
    sStartparArray[0]= &fStartVertex[0];
    sStartparArray[1]= &fStartVertex[1];
    sStartparArray[2]= &fStartVertex[2];
    sStartparArray[3]= &fStartTime0;
  }


  void FitterAW::BeginOfRun(DS::Run& run)
  {
    DB* db = DB::Get();
    DBLinkPtr dbLink = db->GetLink( "FIT_FITTERAW" );
    fAirRI = dbLink->GetD( "air_RI" );
    fWaterRI = dbLink->GetD( "water_RI" );
    fMaxRes = dbLink->GetD( "max_PMT_res" );
    MultiPathFitter::fMaximumIterations = dbLink->GetI( "maximum_iterations" );
    fMaximumReflections = dbLink->GetI( "maximum_reflections");
    fStepsBelowTolerance = dbLink->GetI( "steps_below_tolerance" );
    fPDFData = dbLink->GetDArray( "sPDF_multipathfit" );
    fMaximumIterations = dbLink->GetI( "maximum_iterations" );
    fNumFitFactors = dbLink->GetI( "num_fit_factors" );
    fGoodFitThreshold = dbLink->GetI( "low_DeltaLogL_threshold" );
    fMaxNumStartPos = dbLink->GetI( "max_start_positions" );
    fRhoCut = dbLink->GetD( "rho_sampling_cut" );
    fBinWidth = dbLink->GetD( "time_bin_width" );
    fTimeOffset = dbLink->GetD( "time_offset" );
    fFitTolerance = dbLink->GetD( "fit_tolerance" );
    fPSupRad = dbLink->GetD( "psup_rad" );

    fSpeedOfLight = CLHEP::c_light;
    fSpeedOfLightWater = fSpeedOfLight/fWaterRI;
    fSpeedOfLightAir = fSpeedOfLight/fAirRI;

    DBLinkPtr innerAvDb = db->GetLink( "GEO", "inner_av" );
    fWaterLevel = innerAvDb->GetD( "split_z" );

    DBLinkPtr avDB = db->GetLink( "GEO", "av" );
    fAVRadius = db->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_sphere" );

    //set entries now so don't have to keep the
    fPDFentries = fPDFData.size();
    //  Add a very slight slope to the original and final bins to get
    //  rid of singular matrices
    double xtemp = fPDFData[0];
    for(int j = fPDFentries-1; j >= 0; j-- ) {
      if(fPDFData[j] == fPDFData[0]) {
        fPDFData[j] = xtemp;
        xtemp -= 0.1;
      }
    }
    xtemp = fPDFData[ fPDFentries-1 ];
    int k;
    for( k = fPDFentries-1; fPDFData[k] == xtemp; k-- );
    for( k++; k<fPDFentries; k++ ) {
      fPDFData[k] = xtemp; xtemp -= 0.1;
    }
    for( int j = 0; j < fPDFentries; j++ ) {
      fPDFX.push_back(fBinWidth * static_cast<Double_t> (j + 1) - fTimeOffset);
    }
    for( int j = 0; j < fPDFentries - 1; j++ ) {
      fDerivative.push_back( ( fPDFData[j + 1] - fPDFData[j] ) * 1/fBinWidth );
    }
    for( int ix = 0; ix < 51; ix++ ) {
      double lx = -5+ix*0.2;
      double x = exp( lx );
      for( int idh = 0; idh < 51; idh++ ) {
        double ldh = -5+idh*0.2;
        double dh = exp( ldh );
        double f = Factor( x,dh );
        fFactor[ix][idh]=f;
      }
    }
    //calculate cos(theta),reflection, and refraction probabilities:
    for ( double ct = 0; ct < 1.; ct+=0.002 ) {
      double R1 = 0.,R2=0.;
      double st = sqrt( 1-ct*ct );
      double nct = fWaterRI*ct;
      double srist2 = (1-fWaterRI*fWaterRI*st*st);
      if( srist2 < 0 ) {
        R1 = 1;  //total internal reflection
        R2 = 1;
      } else { // refraction and reflection
        srist2 = sqrt( srist2 );
        R1 = pow( ( nct-srist2 )/( nct+srist2 ), 2 );
        R2 = pow( ( ct-fWaterRI*srist2 )/( ct+fWaterRI*srist2 ), 2 );
      }
      fRWater.push_back( ( R1+R2 )/2. );
      srist2 = sqrt(fWaterRI*fWaterRI-st*st);
      R1 = pow( ( ct-srist2 )/( ct+srist2 ),2 );
      R2 = pow(( fWaterRI*fWaterRI*ct-srist2)/(fWaterRI*fWaterRI*ct+srist2), 2 );
      fRAir.push_back( ( R1+R2 )/2. );
    }
    fRAir.push_back( 0 );
    fRWater.push_back( 0 );
  }


  void FitterAW::CalculatePMTOther()
  {
    fBestFitDirection += CalculatePMTDirection();
  }


  void FitterAW::SummariseOtherOutputs()
  {
    fBestFitDirection=fBestFitDirection.Unit();
    fFOM=MultiPathFitter::sBestLikelihood/MultiPathFitter::sNHits;
  }


  // returns Likelihood and derivative for a given time of flight t.
  double FitterAW::LAnddLdt( double t, double &dLdt )
  {
    double L = 0.;
    fRes=fT-t-fTime0;
    if( fRes > fMaxRes )fRes = fMaxRes;//longer than valid time (late light, crosstalk?)
    if( fRes < -fMaxRes )fRes = -fMaxRes;//longer than valid time (late light, crosstalk?)
    long ibin = static_cast<int>(floor( ( fRes + 100.0 ) * 1/fBinWidth ));
    //proper derivative won't be calculated for the first and last bin, but
    //will be calculated for the other bins
    if( ibin > 0 && ibin < fPDFentries- 2 ) {
      double dt = fRes-fPDFX[ibin];
      L = ( fPDFData[ibin] + fDerivative[ibin] * dt );
      // find the derivertive between bins for less than and greater than
      // the bin centre
      if( dt < 400 / (fPDFentries*2.) ) {
        double u = 0.5+dt*(1/fBinWidth);
        dLdt = -u*fDerivative[ibin]-(1-u)*fDerivative[ibin-1];
      } else {
          double u = (1/fBinWidth)*dt-0.5;
          dLdt = -(1-u)*fDerivative[ibin]-u*fDerivative[ibin+1];
        }
    } else if (ibin <= 0) {
        L = +fPDFData[0];
        dLdt = -fDerivative[0];
      } else {
          L = +fPDFData[fPDFentries-2];
          dLdt = -fDerivative[fPDFentries-2];
        }
    return L;
  }

  double FitterAW::X( double *theta, double *p )
  {
    // parameter is the vertical distance between PMT and vertex
    return fDepth*tan( *theta )+ *p *tan(asin( fWaterRI*sin( *theta ) ) );
  }

  // calculates fraction of transverse distance from dh to the surface
  // of the water
  double FitterAW::Factor( double transpos, double dh )
  {
    double height = 1;
    fX->SetParameters( &height );
    double saveDepth = fDepth;
    fDepth = dh;
    double theta = fX->GetX(transpos);
    fDepth = saveDepth;
    //starting point is below intersection
    if( abs( transpos ) < 0.001 )return 0;
    // dh tan(theta)/x= fraction of straight line distance to point
    // directly below intersection.
    else return dh*tan( theta )/transpos;
  }


  double FitterAW::FactorN( double transpos, bool negative )
  {
    //calculates fraction of transverse distance from fDepth to the
    //surface by interpolating, based on fDepth and fHeight.
    double ret2 = 0, dh = 0, lx = 0;
    if(fHeight == 0) ret2 = 1;//vertex at water surface
    else{
      if ( negative == false ){// if the PMT is above the water
        dh = fDepth/fHeight;
        //transpos is projection of vector (vertex to PMT) onto detector z-axis
        //lx is the log of the ratio of transpos and fHeight to define the extreme
        //conditions
        lx = log(transpos/fHeight);
      }
      else if (negative == true ){// if the PMT is below the water
        dh = fHeight/fDepth;
        lx = log(transpos/-fDepth);
      }
      double ldh = log(dh);
      if( lx < -5 ) ret2 = ( fWaterRI*dh+1 )/( 1+dh );//sVertex close to water level
      else {
        if( lx>5 ) lx = 5;//vetex too low
        if( ldh<-5 ) ldh = -5; // water level too low
        else if( ldh > 5 )ldh = 5.0;// water level too high
        double ux=( lx+5 )/.2;
        int ix = floor(ux);
        ux = ux-ix;
        double udh = ( ldh+5 )/.2;
        int idh = floor( udh );
        udh = udh-idh;
        ret2 = ( 1-ux )*( ( 1-udh )*fFactor[ix][idh]+udh*fFactor[ix][idh+1] )+
        ux*( ( 1-udh )*fFactor[ix+1][idh]+udh*fFactor[ix+1][idh+1]);
      }
    }
    return ret2;
  }


  void FitterAW::VertexInBottomMat(double &L, double &dLdt, double &dLdx, double &dLdy, double &dLdz, double &dLdt0)
  {
    TVector3 diff=fPosition-fVertex;
    if( fHeight<0 ){  //Pmt and vertex are below water level
      // both PMT and vertex are under water -correct the direct time of flight
      PMTVertexSame(L, dLdx, dLdy, dLdz, dLdt0, false);
    } else {  //PMT is above water, vertex below- refraction
      PMTVertexDiff(L, dLdt, dLdx, dLdy, dLdz, dLdt0, false);
      }
  }

  TVector3 FitterAW::VertexInBottomMat() {
    TVector3 diff=fPosition-fVertex;
    TVector3 direction;
    if( fHeight<0 ){  //Pmt and vertex are below water level
      // both PMT and vertex are under water -correct the direct time of flight
      direction = PMTVertexSame( false );
    } else {  //PMT is above water, vertex below- refraction
      direction = PMTVertexDiff( false );
    }
    return direction;
  }


  void FitterAW::VertexInTopMat(double &L, double &dLdt, double &dLdx, double &dLdy, double &dLdz, double &dLdt0)
  {
    TVector3 diff=fPosition-fVertex;
    if( fHeight>0 ) {
      //pmt is also above water level- find direct and reflected light paths
      PMTVertexSame(L, dLdx, dLdy, dLdz, dLdt0, true);
    } else {  //vertex above water, PMT below- find refracted path
      PMTVertexDiff(L, dLdt, dLdx, dLdy, dLdz, dLdt0, true);
    }
  }


  TVector3 FitterAW::VertexInTopMat() {
    TVector3 diff=fPosition-fVertex;
    TVector3 direction;
    if( fHeight>0 ) {
      //pmt is also above water level- find direct and reflected light paths
      direction = PMTVertexSame( true );
    } else { //vertex above water, PMT below- find refracted path
        direction = PMTVertexDiff( true );
      }
    return direction;
  }


  void FitterAW::PMTVertexDiff(double &L, double &dLdt, double &dLdx, double &dLdy, double &dLdz, double &dLdt0,  bool vertexTopMat){
    double sSpeedMaterial1=0, sSpeedMaterial2=0,alpha = 0;
    TVector3 s;
    TVector3 diff=fPosition-fVertex;
    if(vertexTopMat == false) {
      sSpeedMaterial1 = fSpeedOfLightWater;
      sSpeedMaterial2 = fSpeedOfLightAir;
      alpha = FactorN( diff.Perp(), false );
      //find a point directly below the intersection with water
      s = ( 1-alpha )*fVertex+alpha*fPosition;
    } else if(vertexTopMat == true) {
      sSpeedMaterial1 = fSpeedOfLightAir;
      sSpeedMaterial2 = fSpeedOfLightWater;
      alpha = FactorN( diff.Perp() , true );
      s = alpha*fVertex+( 1-alpha )*fPosition;
    }
    s.SetZ( fWaterLevel );  //move intersection to water surface.
    TVector3 incident = s-fVertex;
    TVector3 refracted = fPosition-s;
    //time of flight for reflected ray
    double t = incident.Mag()/sSpeedMaterial1+refracted.Mag()/sSpeedMaterial2;
    incident = incident.Unit();
    // Refracted light that hits the water surface.
    // 1.0001 means that we don't get an in infinity if reflected=1.
    L = ( LAnddLdt( t, dLdt ) );

    dLdt = dLdt/L;
    L = log( L );
    diff = diff.Unit(); //components are dt/dx,dt/dy,dt/dz; dt/dt0=-1
    dLdx = -dLdt*incident.X()/sSpeedMaterial1;
    dLdy = -dLdt*incident.Y()/sSpeedMaterial1;
    dLdz = -dLdt*incident.Z()/sSpeedMaterial1;
    dLdt0 = dLdt;
  }


  TVector3 FitterAW::PMTVertexDiff( bool vertexTopMat ) {
    double sSpeedMaterial1=0, sSpeedMaterial2=0,alpha = 0;
    TVector3 s;
    TVector3 diff=fPosition-fVertex;
    if(vertexTopMat == false) {
      sSpeedMaterial1 = fSpeedOfLightWater;
      sSpeedMaterial2 = fSpeedOfLightAir;
      alpha = FactorN( diff.Perp(), false );
      //find a point directly below the intersection with water
      s = ( 1-alpha )*fVertex+alpha*fPosition;
    } else if(vertexTopMat == true) {
        sSpeedMaterial1 = fSpeedOfLightAir;
        sSpeedMaterial2 = fSpeedOfLightWater;
        alpha = FactorN( diff.Perp() , true );
        s = alpha*fVertex+( 1-alpha )*fPosition;
    }

    s.SetZ( fWaterLevel ); //move intersection to water surface.
    TVector3 incident = s-fVertex;
    TVector3 refracted = fPosition-s;
    //time of flight for reflected ray
    double t = incident.Mag()/sSpeedMaterial1+refracted.Mag()/sSpeedMaterial2;
    incident = incident.Unit();
    // Refracted light that hits the water surface.
    // 1.0001 means that we don't get an in infinity if reflected=1.
    double dLdt = 0.;
    double L = ( LAnddLdt( t, dLdt ) );
    TVector3 direction = L*incident;

    return direction;
  }


  void FitterAW::PMTVertexSame(double &L, double &dLdx, double &dLdy, double &dLdz, double &dLdt0,  bool vertexTopMat){
    double sSpeedMaterial1=-1, sSpeedMaterial2=-1,alpha = 0;
    std::vector<double> sRInterface;
    short belowLevel=0;
    TVector3 s;
    if(vertexTopMat == false) {
      // find fractional transverse distance between vertex and water intersection
      alpha = fDepth/( fDepth-fHeight );
      // find a point directly below the intersection with water
      s = ( 1-alpha )*fVertex+alpha*fPosition;
      sSpeedMaterial2 = fSpeedOfLightAir/fWaterRI;
      sSpeedMaterial1 = fSpeedOfLightWater;
      sRInterface = fRWater;
      belowLevel = 1;
    }
    if(vertexTopMat == true) {
      //find fractional transverse distance between pmt and water intersection
      alpha = fHeight/( fHeight-fDepth );
      //find a point directly above the intersection with water
      s = alpha*fVertex+( 1-alpha )*fPosition;
      sSpeedMaterial2 = fSpeedOfLightAir;
      sSpeedMaterial1 = fSpeedOfLightAir;
      sRInterface = fRAir;
      belowLevel = -1;
    }
    TVector3 diff=fPosition-fVertex;
    // both PMT and vertex are under water -correct the direct time of flight
    double t = diff.Mag()/sSpeedMaterial2;
    s.SetZ( fWaterLevel );  //move intersection to water surface.
    TVector3 incident = s-fVertex;
    TVector3 reflectedRay = fPosition-s;
    // time of flight for reflected ray
    double treflected = ( incident.Mag()+reflectedRay.Mag() )/sSpeedMaterial1;
    incident = incident.Unit();
    // find cos(theta) between vertical and photon trajectory,multiply
    // by fMaximumReflections to look up in table
    double ct = belowLevel*incident.Z()*fMaximumReflections;
    int ix = static_cast<int>(ct);
    double u = ct-ix;
    // find reflection probability for this cos theta, interpolated
    double reflected = ( 1-u )*sRInterface[ix]+u*sRInterface[ix+1];
    double dLdt1 = 0., dLdt2 = 0.;
    // Second term is for direct light that never hits the water surface;
    // first term for water-water reflections
    L = (reflected*LAnddLdt( treflected, dLdt2 )+LAnddLdt( t, dLdt1 ) );
    dLdt1 = dLdt1/L;
    dLdt2 = reflected*dLdt2/L;
    L = log( L/( 1+reflected ) );
    // components are proportional to dt/dx,dt/dy,dt/dz; dt/dt0=-1;
    // incident is already set to a unit vector
    diff = diff.Unit();
    dLdx = -( dLdt1*diff.X()+dLdt2*incident.X() )/sSpeedMaterial1;
    dLdy = -( dLdt1*diff.Y()+dLdt2*incident.Y() )/sSpeedMaterial1;
    dLdz = -( dLdt1*diff.Z()+dLdt2*incident.Z() )/sSpeedMaterial1;
    dLdt0 = dLdt1+dLdt2;
  }

  TVector3 FitterAW::PMTVertexSame( bool vertexTopMat ) {
    double sSpeedMaterial1 = -1, sSpeedMaterial2 = -1,alpha = 0;
    std::vector<double> sRInterface;
    short belowLevel = 0;
    TVector3 s;
    if(vertexTopMat == false) {
      // find fractional transverse distance between vertex and water intersection
      alpha = fDepth/( fDepth-fHeight );
      // find a point directly below the intersection with water
      s = ( 1-alpha )*fVertex+alpha*fPosition;
      sSpeedMaterial2 = fSpeedOfLightAir/fWaterRI;
      sSpeedMaterial1 = fSpeedOfLightWater;
      sRInterface = fRWater;
      belowLevel = 1;
    }
    if(vertexTopMat == true) {
      //find fractional transverse distance between pmt and water intersection
      alpha = fHeight/( fHeight-fDepth );
      //find a point directly above the intersection with water
      s=alpha*fVertex+( 1-alpha )*fPosition;
      sSpeedMaterial2 = fSpeedOfLightAir;
      sSpeedMaterial1 = fSpeedOfLightAir;
      sRInterface = fRAir;
      belowLevel = -1;
    }
    TVector3 diff=fPosition-fVertex;
    // both PMT and vertex are under water -correct the direct time of flight
    double t = diff.Mag()/sSpeedMaterial2;
    s.SetZ( fWaterLevel );  //move intersection to water surface.
    TVector3 incident = s-fVertex;
    TVector3 reflectedRay = fPosition-s;
    // time of flight for reflected ray
    double treflected = ( incident.Mag()+reflectedRay.Mag() )/sSpeedMaterial1;
    incident = incident.Unit();
    // find cos(theta) between vertical and photon trajectory,multiply
    // by fMaximumReflections to look up in table
    double ct = belowLevel*incident.Z()*fMaximumReflections;
    int ix = static_cast<int>(ct);
    double u = ct-ix;
    // find reflection probability for this cos theta, interpolated
    double reflected = ( 1-u )*sRInterface[ix]+u*sRInterface[ix+1];
    // Second term is for direct light that never hits the water surface;
    // first term for water-water reflections
    double dLdt1 = 0, dLdt2 = 0;
    double L2 = reflected*LAnddLdt( treflected, dLdt2 );
    double L1 = LAnddLdt( t, dLdt1 );
    TVector3 direction = L1*diff.Unit()+L2*incident;
    return direction;
  }

  void FitterAW::CalculatePMTLikelihood()
  {
    double L=0, dLdt=0, dLdx=0, dLdy=0,dLdz=0,dLdt0=0;
    if( fWaterLevel >= -fAVRadius && fWaterLevel <= fAVRadius ){
      fDepth = fWaterLevel-fVertex.Z();
      if( fDepth>0 ){  //vertex is below water level
        VertexInBottomMat(L, dLdt, dLdx, dLdy, dLdz, dLdt0);
      } else { //vertex is above water level (fDepth is negative)
          VertexInTopMat(L, dLdt, dLdx, dLdy, dLdz, dLdt0);
        }
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
      sLikelihood += L;
    } else { Log::Die( "Not a valid water height "); }
  }


  TVector3 FitterAW::CalculatePMTDirection()
  {
    TVector3 direction;
    if( fDepth>0 ) {  //vertex is below water level
      direction = VertexInBottomMat();
    }else { //vertex is above water level (fDepth is negative)
      direction = VertexInTopMat();
    }
    return direction;
  }

  void FitterAW::MakeParameterHistogram()
  {
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
    TH1F sHx= TH1F(sNameHx,sNameHx,2000,-10000,10000);
    TH1F sHy= TH1F(sNameHy,sNameHy,2000,-10000,10000);
    TH1F sHz= TH1F(sNameHz,sNameHz,2000,-10000,10000);
    TH1F sHt= TH1F(sNameHt,sNameHt,400,-100,300);

    //histograms for liklihod of dirivative
    TH1F sDHx= TH1F(sNameDHx,sNameDHx,2000,-10000,10000);
    TH1F sDHy= TH1F(sNameDHy,sNameDHy,2000,-10000,10000);
    TH1F sDHz= TH1F(sNameDHz,sNameDHz,2000,-10000,10000);
    TH1F sDHt= TH1F(sNameDHt,sNameDHt,400,-100,300);

    //histogram for likelihood of numrical derivative
    TH1F sNDHx= TH1F(sNameNDHx,sNameNDHx,2000,-10000,10000);
    TH1F sNDHy= TH1F(sNameNDHy,sNameNDHy,2000,-10000,10000);
    TH1F sNDHz= TH1F(sNameNDHz,sNameNDHz,2000,-10000,10000);
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
