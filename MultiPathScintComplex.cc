////////////////////////////////////////////////////////////////////////
///
/// \class
///
/// \brief
///
/// \author David Auty auty@ualberta.ca
/// 
/// REVISION HISTORY: Jie Hu jhu9@ualberta.ca - modified to fit J Tseng's MultiPathFunction\n
///
/// \detail
///       To fit data in a partial scint-water fill geometry
////
////////////////////////////////////////////////////////////////////////

#include <RAT/MultiPathScint.hh>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

#define DEBUG

using namespace RAT;
using namespace ROOT;

MultiPathScint::MultiPathScint() {
  // parameter definitions
  fPars.push_back(Parameter("x",2000,-10000.0,10000.0));
  fPars.push_back(Parameter("y",2000,-10000.0,10000.0));
  fPars.push_back(Parameter("z",2000,-10000.0,10000.0));
  fPars.push_back(Parameter("t",400,-100.0,300.0));
}

void MultiPathScint::BeginOfRun(DS::Run&) {
#ifdef MULTIPATH_FIXED_SEED
  //step = 0.0001;
  step = 0.5001;
#endif
  DB* db = DB::Get();
  DBLinkPtr dbLink = db->GetLink( "FIT_MULTIPATH" );
  fPDFDataCh = dbLink->GetDArray( "sPDF_scint" );
  double fScintRIeff = dbLink->GetD( "scint_RI" );
  double fWaterRIeff = dbLink->GetD( "water_RI" );
  double fChBinWidth = dbLink->GetD( "time_bin_width" );
  double fChOffset = dbLink->GetD( "time_offset" );
  fMaxPMTtRes = dbLink->GetD( "max_PMT_tRes" );
  double r = dbLink->GetD( "psup_rad" );
  fPSUPRadius2 = r * r;
  fAVRadius = 6050; 

  double fSpeedOfLight = CLHEP::c_light;
  fSpeedOfLightScint = fSpeedOfLight / fScintRIeff;
  fSpeedOfLightWater = fSpeedOfLight / fWaterRIeff;
  fEntriesTimeCh = fPDFDataCh.size();
  double xtemp = fPDFDataCh[0];
  int j;
  for( j = fEntriesTimeCh - 1; j >= 0; j-- ) {
    if( fPDFDataCh[j] == fPDFDataCh[0] ) {
      fPDFDataCh[j] = xtemp;
      xtemp -= 0.001;
    }
  }

  xtemp = fPDFDataCh[ fEntriesTimeCh-1 ];

  for( j = fEntriesTimeCh - 1; fPDFDataCh[j] == xtemp; j-- );

  for( j++; j < fEntriesTimeCh; j++ ) {
    fPDFDataCh[j] = xtemp;
    xtemp -= 0.001;
  }

  for( j = 0; j < fEntriesTimeCh; j++ ) {
    fPDFxCh.push_back(fChBinWidth * static_cast<double>(j + 1) - fChOffset);
    fPDFCh.push_back(fPDFDataCh[j]);
  }

  for( j = 0; j < fEntriesTimeCh - 1; j++ )
    fDerivativeCh.push_back( ( fPDFCh[j + 1] - fPDFCh[j] ) / fChBinWidth );

#ifdef DEBUG
  //Check the fitting pdfs and their derivatives for any event
  TFile f("pdfDerivatives_MultiPathScint.root","RECREATE");
  TH1F h1("h1","pdfCherenkov",fEntriesTimeCh, -100, 300);
  TH1F h1x("h1x","pdfXCh",fEntriesTimeCh, -100, 300);
  TH1F dh1("dh1","DerivativeCherenkov",fEntriesTimeCh,-100, 300);

  for( j = 0; j < fEntriesTimeCh - 1; j++ ) {
    h1.SetBinContent(j+1, fPDFCh[j]);
    h1x.SetBinContent(j+1, fPDFxCh[j]);
    dh1.SetBinContent(j+1, fDerivativeCh[j]);
  }

  f.cd();
  h1.Write();
  h1x.Write();
  dh1.Write();
  f.Close();
#endif // DEBUG
}

void MultiPathScint::Initialise(const std::string&) {
}

void MultiPathScint::SetSeed(const DS::FitResult&) {
}

DS::FitResult MultiPathScint::MakeResult(const std::vector<double>& pars) {
  TVector3 fittedVertex(pars[0], pars[1], pars[2]);
  DS::FitVertex vertex;
  vertex.SetPosition(fittedVertex);
  vertex.SetTime(pars[3]);
  DS::FitResult result;
  result.SetVertex(0, vertex);
  return result;
}

bool MultiPathScint::ValidResult(const std::vector<double>& pars) {
  double r2 = pars[0]*pars[0] + pars[1]*pars[1] + pars[2]*pars[2];
  return r2 < fPSUPRadius2;
}

std::vector<double> MultiPathScint::GetStart() {
  // Could set up a seed here, but I'm not using it for now (from aTime and aVertex)
#ifdef MULTIPATH_FIXED_SEED
  double ran0 = step;
  double ran1 = 2.0*step - 1.0;
  double ran2Pi = step*2.0*CLHEP::pi;
  double t = 100.0 + 200.0*step;
  step += 0.0001;
  if (step > 1.0) step = 0.0001;
#else
  double ran0 = CLHEP::RandFlat::shoot( 0.0, 1.0 );
  double ran1 = CLHEP::RandFlat::shoot( -1.0, 1.0 );
  double ran2Pi = CLHEP::RandFlat::shoot( 0.0, 2.0*CLHEP::pi );
  double t = CLHEP::RandFlat::shoot(100.0, 300.0);
#endif

  double r = pow(ran0, 1.0/3.0) * 10000.0; // mm
  double costheta = ran1;
  double sintheta = sqrt(1.0 - costheta*costheta);
  std::vector<double> v;
  v.push_back(r*sintheta*cos(ran2Pi));
  v.push_back(r*sintheta*sin(ran2Pi));
  v.push_back(r*costheta);
  // in principle, MC and data could have different PMT hit time offset,
  // but in latest incarnation of MultiPathFitter code, this doesn't seem to be the case
  v.push_back(t);
  return v;
}

void MultiPathScint::LAnddLdt_Ch(double& L_Ch, double& dLdt_Ch) {
  if (fRes > fMaxPMTtRes) fRes = fMaxPMTtRes; //max_PMT_tRes
  else if (fRes < -fMaxPMTtRes) fRes = -fMaxPMTtRes;

  int ibin = (int) floor( ( fRes + 100.0 ) * 4.0 );
  if(ibin > 0 && ibin < fEntriesTimeCh - 2) {
    double dt = fRes-fPDFxCh[ibin];
    L_Ch = (fPDFCh[ibin] + fDerivativeCh[ibin] * dt);
    if (dt < 0.125) {
      double u = 0.5 + dt * 4.0;
      dLdt_Ch = -u*fDerivativeCh[ibin] - (1.0-u)*fDerivativeCh[ibin-1];
    } else {
      double u = 4.0 * dt - 0.5;
      dLdt_Ch = -(1.0-u)*fDerivativeCh[ibin] - u*fDerivativeCh[ibin+1];
    }// end of the if statement if dt<0.125
  } else if (ibin <= 0) {
    L_Ch = fPDFCh[0];
    dLdt_Ch = -fDerivativeCh[0];
  } else {
    L_Ch = fPDFCh[fEntriesTimeCh - 2];
    dLdt_Ch = -fDerivativeCh[fEntriesTimeCh - 2];
  }
}

void MultiPathScint::CalculateScintExWaterPath( double *vtxdata, std::vector<double>& pathResults )
{
  double pathInScint = 0, dpdx = 0, dpdy = 0, dpdz = 0;
  TVector3 startpos, pmtpos, fitpos, xpDiff, pathDiff, incidentDirect;
  startpos.SetXYZ(vtxdata[0],vtxdata[1],vtxdata[2]);
  pmtpos.SetXYZ(vtxdata[3], vtxdata[4], vtxdata[5]);
  pathDiff = (pmtpos - startpos);
  xpDiff = (pmtpos - fitpos);
  incidentDirect = (pmtpos - startpos).Unit();

  if(startpos.Mag()<= fAVRadius)
  {
     double sqrVal = pow((startpos * incidentDirect),2) - (startpos).Mag2() + fAVRadius * fAVRadius; 
     if(sqrVal<0) pathInScint = 0;
     else {
       pathInScint = -(startpos * incidentDirect) + sqrt(sqrVal);// Note: this is the one of two solutions from the quadratic equation
       double dsndx = 0, dsndy = 0, dsndz = 0, dsqrValdx = 0, dsqrValdy = 0, dsqrValdz = 0;
       dsndx = (pmtpos.X()-2*startpos.X())/pathDiff.Mag() + (startpos*pmtpos - startpos*startpos)*( (pmtpos-startpos).X() )/pow(pathDiff.Mag(),3);
       dsndy = (pmtpos.Y()-2*startpos.Y())/pathDiff.Mag() + (startpos*pmtpos - startpos*startpos)*( (pmtpos-startpos).Y() )/pow(pathDiff.Mag(),3);
       dsndz = (pmtpos.Z()-2*startpos.Z())/pathDiff.Mag() + (startpos*pmtpos - startpos*startpos)*( (pmtpos-startpos).Z() )/pow(pathDiff.Mag(),3);
       dsqrValdx = 1./sqrt(sqrVal)*( (startpos*incidentDirect)*dsndx - startpos.X() );
       dsqrValdy = 1./sqrt(sqrVal)*( (startpos*incidentDirect)*dsndy - startpos.Y() );
       dsqrValdz = 1./sqrt(sqrVal)*( (startpos*incidentDirect)*dsndz - startpos.Z() );
       dpdx = -1*dsndx + dsqrValdx;
       dpdy = -1*dsndy + dsqrValdy;
       dpdz = -1*dsndz + dsqrValdz;
     }
  }
  else { //for exwater vertex, check whehther ray pass through AV
     if( startpos*pmtpos/pathDiff.Mag() < fAVRadius )
     {
       double sqrVal1 = pow((startpos * incidentDirect), 2) - (startpos).Mag2() + fAVRadius * fAVRadius;                         
       double sqrVal2 = pow((pmtpos * incidentDirect), 2) - (pmtpos).Mag2() + fAVRadius * fAVRadius;
       if( sqrVal1<0 || sqrVal2<0 ) pathInScint = 0;
       else {
         double pathInWater1 = -(startpos * incidentDirect) + sqrt(sqrVal1);
         double pathInWater2 = (pmtpos * incidentDirect) - sqrt(sqrVal2);
         pathInScint = pathDiff.Mag() - pathInWater1 - pathInWater2;
         double dsn1dx = 0, dsn1dy = 0, dsn1dz = 0, dsqrVal1dx = 0, dsqrVal1dy = 0, dsqrVal1dz = 0;
         double dsn2dx = 0, dsn2dy = 0, dsn2dz = 0, dsqrVal2dx = 0, dsqrVal2dy = 0, dsqrVal2dz = 0;
         double dP1dx = 0, dP1dy = 0, dP1dz = 0;
         double dP2dx = 0, dP2dy = 0, dP2dz = 0;
         // dpathInWater1/dx, ...
         dsn1dx = (pmtpos.X()-2*startpos.X())/pathDiff.Mag() + (startpos*pmtpos - startpos*startpos)*( (pmtpos-startpos).X() )/pow(pathDiff.Mag(),3);
         dsn1dy = (pmtpos.Y()-2*startpos.Y())/pathDiff.Mag() + (startpos*pmtpos - startpos*startpos)*( (pmtpos-startpos).Y() )/pow(pathDiff.Mag(),3);
         dsn1dz = (pmtpos.Z()-2*startpos.Z())/pathDiff.Mag() + (startpos*pmtpos - startpos*startpos)*( (pmtpos-startpos).Z() )/pow(pathDiff.Mag(),3);
         dsqrVal1dx = 1./sqrt(sqrVal1)*( (startpos*incidentDirect)*dsn1dx - startpos.X() );
         dsqrVal1dy = 1./sqrt(sqrVal1)*( (startpos*incidentDirect)*dsn1dy - startpos.Y() );
         dsqrVal1dz = 1./sqrt(sqrVal1)*( (startpos*incidentDirect)*dsn1dz - startpos.Z() );
         dP1dx = -1*dsn1dx + dsqrVal1dx;
         dP1dy = -1*dsn1dy + dsqrVal1dy;
         dP1dz = -1*dsn1dz + dsqrVal1dz;
 
         // dpathInWater2/dx, ...    
         dsn2dx = (-pmtpos.X())/pathDiff.Mag() + (pmtpos*pmtpos - pmtpos*startpos)*( (pmtpos-startpos).X() )/pow(pathDiff.Mag(),3);
         dsn2dy = (-pmtpos.Y())/pathDiff.Mag() + (pmtpos*pmtpos - pmtpos*startpos)*( (pmtpos-startpos).Y() )/pow(pathDiff.Mag(),3);
         dsn2dz = (-pmtpos.Z())/pathDiff.Mag() + (pmtpos*pmtpos - pmtpos*startpos)*( (pmtpos-startpos).Z() )/pow(pathDiff.Mag(),3);
         dsqrVal2dx = 1./sqrt(sqrVal2)*( (pmtpos*incidentDirect)*dsn2dx );
         dsqrVal2dy = 1./sqrt(sqrVal2)*( (pmtpos*incidentDirect)*dsn2dy );
         dsqrVal2dz = 1./sqrt(sqrVal2)*( (pmtpos*incidentDirect)*dsn2dz );
         dP2dx = dsn2dx - dsqrVal2dx;
         dP2dy = dsn2dy - dsqrVal2dy;
         dP2dz = dsn2dz - dsqrVal2dz;

         dpdx = -(pmtpos-startpos).X()/pathDiff.Mag() - dP1dx - dP2dx;
         dpdy = -(pmtpos-startpos).Y()/pathDiff.Mag() - dP1dy - dP2dy;
         dpdz = -(pmtpos-startpos).Z()/pathDiff.Mag() - dP1dz - dP2dz;
       }
     }
  }

  pathResults[0] = pathInScint;
  pathResults[1] = dpdx;
  pathResults[2] = dpdy;
  pathResults[3] = dpdz;
}

double MultiPathScint::Calculate(const std::vector<double>& pmt,
                                         const std::vector<double>& par,
                                         std::vector<double>& diff) {

  double dx = pmt[0] - par[0];
  double dy = pmt[1] - par[1];
  double dz = pmt[2] - par[2];
  double dr = sqrt(dx*dx + dy*dy + dz*dz);

  double scintpath = 0;
  std::vector<double> pathResults(4,0.0);
  double vtxdata[6] = {par[0], par[1], par[2], pmt[0], pmt[1], pmt[2]}; 
  CalculateScintExWaterPath( vtxdata, pathResults ); 
  scintpath = pathResults[0];
  double dpdx = pathResults[1];
  double dpdy = pathResults[2];
  double dpdz = pathResults[3];

  double tof = scintpath / fSpeedOfLightScint + (dr - scintpath)/ fSpeedOfLightWater;

  fRes = pmt[3] - tof - par[3]; // residual:  (PMT time - distance) - trial time

  double L = 0.0;
  double dLdtCh = 0.0;
  LAnddLdt_Ch(L, dLdtCh); // uses fRes to get L and dLdtCh

  double dLdt0 = dLdtCh / L;
//  double factor = -dLdt0 / (dr * fSpeedOfLightScint );
  double dLdx = dLdt0*( dpdx/fSpeedOfLightScint + ( -dx/dr - dpdx )/fSpeedOfLightWater );
  double dLdy = dLdt0*( dpdy/fSpeedOfLightScint + ( -dy/dr - dpdy )/fSpeedOfLightWater );
  double dLdz = dLdt0*( dpdz/fSpeedOfLightScint + ( -dz/dr - dpdz )/fSpeedOfLightWater );

  diff.assign(4, 0.0);
  diff[0] = dLdx;
  diff[1] = dLdy;
  diff[2] = dLdz;
  diff[3] = dLdt0;

  L = log(L);

#ifdef DEBUG
//  debug << "AWP::CL pmt=" << pmt[0] << "," << pmt[1] << "," << pmt[2] << ")\n";
//  debug << "        fRes=" << fRes << " dL=" << dLdtCh << "\n";
//  debug << "        L=" << L << " diff=(" << diff[0] << "," << diff[1] << "," << diff[2] << "," << diff[3] << ")\n";
#endif

  return L;
}

