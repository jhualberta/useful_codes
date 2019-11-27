////////////////////////////////////////////////////////////////////////
///
/// \class
///
/// \brief
///
/// \author David Auty auty@ualberta.ca
///
/// REVISION HISTORY: Jie Hu jhu9@ualberta.ca - new strategy and modified
/// to fit J Tseng's MultiPathFunction
///
/// \detail
///       To fit data in a partial scint-water fill geometry
////
////////////////////////////////////////////////////////////////////////

#include <RAT/MultiPathScintWater.hh>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include <TF1.h>
//#define DEBUG

using namespace RAT;
using namespace ROOT;

MultiPathScintWater::MultiPathScintWater() {
  // parameter definitions
  fPars.push_back(Parameter("x",1678,-8390.0,8390.0));
  fPars.push_back(Parameter("y",1678,-8390.0,8390.0));
  fPars.push_back(Parameter("z",1678,-8390.0,8390.0));
  fPars.push_back(Parameter("t",400,-100.0,300.0));
  fPars.push_back(Parameter("wl",2000,4108.0,6108.0));
}

void MultiPathScintWater::BeginOfRun(DS::Run&) {
#ifdef MULTIPATH_FIXED_SEED
  //step = 0.0001;
  step = 0.5001;
#endif
  DB* db = DB::Get();
  DBLinkPtr dbLink = db->GetLink( "FIT_MULTIPATH" );
  fPDFDataCh = dbLink->GetDArray( "sPDF_multipathfit_shift_labppo_0p5" );
  fPDFDataScint = dbLink->GetDArray( "sPDF_scint_labppo_0p5" );
  fScintRI = dbLink->GetD( "scint_RI_labppo_0p5" ); // effective grVelocity in scint
  fWaterRI = dbLink->GetD( "water_RI_tuned_labppo_0p5" ); // effective grVelocity in water, based on partial and scint phase MC
  double fChBinWidth = dbLink->GetD( "time_bin_width" ); // time resolution in ns, set to 0.25 ns
  double fChOffset = dbLink->GetD( "time_offset" );
  fMaxPMTtRes = dbLink->GetD( "max_PMT_tRes" );
  fBoundary = dbLink->GetD( "boundary_tolerance" );
  fNeckpathEnable = dbLink->GetZ( "neckpath_enable" );
  fPSUPRadius = dbLink->GetD( "psup_rad" );
  fPSUPRadius2 = fPSUPRadius * fPSUPRadius;

  fZoff = dbLink->GetD( "av_offset_z" ); //z offset of AV to PSUP center
  fAVRadius = db->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_sphere" ); // 6005.0, AV inner radius
  fAVRadiusOuter = db->GetLink( "SOLID", "acrylic_vessel_outer" )->GetD( "r_sphere" ); // 6060, AV outer radius
  fNeckRadius = db->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_neck" ); // inner radius of neck = 730 mm
  fNeckRadiusOuter = db->GetLink( "SOLID", "acrylic_vessel_outer" )->GetD( "r_neck" ); // outer radius of neck = 785 mm
  fZneckLo = sqrt( fAVRadiusOuter*fAVRadiusOuter - fNeckRadiusOuter*fNeckRadiusOuter ) + fZoff; // bottom of the neck
  double fSpeedOfLight = CLHEP::c_light;
  fSpeedOfLightWater = fSpeedOfLight / fWaterRI;
  fSpeedOfLightScint = fSpeedOfLight / fScintRI;

  fMaxNumStartPos = dbLink->GetI( "max_start_positions" );
  fRhoCut = dbLink->GetD( "rho_sampling_cut" ); // sampling only in AV, for future test
  //fGoodFitThreshold = dbLink->GetI( "low_DeltaLogL_threshold" );// for future test

  DBLinkPtr innerAvDb = db->GetLink( "GEO", "inner_av" );
  //fWaterLevel = innerAvDb->GetD( "split_z" );
  //fWaterLevel = fWaterLevel + fZoff; // correct the water level

  fEntriesTimeCh = fPDFDataCh.size();// 1600
  double xtemp = fPDFDataCh[0];
  int j;
  /// make slow slopes for the flat de-weighted lines in pdfs to avoid zeros in derivatives
  double slope_step = 0.001; // default: 0.001
  for( j = fEntriesTimeCh - 1; j >= 0; j-- ) {
    if( fPDFDataCh[j] == fPDFDataCh[0] ) {
      fPDFDataCh[j] = xtemp;
      xtemp -= slope_step;
    }
  }

  xtemp = fPDFDataCh[ fEntriesTimeCh-1 ];

  for( j = fEntriesTimeCh - 1; fPDFDataCh[j] == xtemp; j-- );

  for( j++; j < fEntriesTimeCh; j++ ) {
    fPDFDataCh[j] = xtemp;
    xtemp -= slope_step;
  }

  for( j = 0; j < fEntriesTimeCh; j++ ) {
    fPDFxCh.push_back(fChBinWidth * static_cast<double>(j + 1) - fChOffset);
    fPDFCh.push_back(fPDFDataCh[j]);
  }

  for( j = 0; j < fEntriesTimeCh - 1; j++ )
    fDerivativeCh.push_back( ( fPDFCh[j + 1] - fPDFCh[j] ) / fChBinWidth );

  /// fill scint pdf
  double xtemp1 = fPDFDataScint[0];
  for( j = fEntriesTimeCh-1; j >= 0; j-- ) {
    if( fPDFDataScint[j] == fPDFDataScint[0] ) {
      fPDFDataScint[j] = xtemp1;
      xtemp1 -= slope_step;
    }
  }

  xtemp1 = fPDFDataScint[ fEntriesTimeCh-1 ];

  for( j = fEntriesTimeCh-1; fPDFDataScint[j] == xtemp1; j-- );

  for( j++; j < fEntriesTimeCh; j++ ) {
    fPDFDataScint[j] = xtemp1;
    xtemp1 -= slope_step;
  }

  for( j = 0; j < fEntriesTimeCh; j++ ) {
    fPDFxScint.push_back(fChBinWidth * static_cast<double>(j + 1) - fChOffset);
    fPDFScint.push_back(fPDFDataScint[j]);
  }

  for( j = 0; j < fEntriesTimeCh - 1; j++ )
    fDerivativeScint.push_back( ( fPDFScint[j + 1] - fPDFScint[j] ) / fChBinWidth );

#ifdef DEBUG
  //Check the fitting pdfs and their derivatives for any event
  TFile f("pdfDerivatives_MultiPathPartial.root","RECREATE");
  TH1F h1("h1","pdfCherenkov",fEntriesTimeCh, -100, 300);
  TH1F h1x("h1x","pdfXCh",fEntriesTimeCh, -100, 300);
  TH1F dh1("dh1","DerivativeCherenkov",fEntriesTimeCh,-100, 300);

  TH1F h2("h2","pdfScint",fEntriesTimeCh, -100, 300);
  TH1F h2x("h2x","pdfXScint",fEntriesTimeCh, -100, 300);
  TH1F dh2("dh2","DerivativeScint",fEntriesTimeCh, -100, 300);

  for( j = 0; j < fEntriesTimeCh; j++ ) {
    h1.SetBinContent(j+1, fPDFCh[j]);
    h1x.SetBinContent(j+1, fPDFxCh[j]);
    dh1.SetBinContent(j+1, fDerivativeCh[j]);
  }
  for( j = 0; j < fEntriesTimeCh; j++ ) {
    h2.SetBinContent(j+1, fPDFScint[j]);
    h2x.SetBinContent(j+1, fPDFxScint[j]);
    dh2.SetBinContent(j+1, fDerivativeScint[j]);
  }

  f.cd();
  h1.Write();h1x.Write();dh1.Write();
  h2.Write();h2x.Write();dh2.Write();
  f.Close();
#endif // DEBUG
}

void MultiPathScintWater::Initialise(const std::string&) {
}

void MultiPathScintWater::SetSeed(const DS::FitResult&) {
}

DS::FitResult MultiPathScintWater::MakeResult(const std::vector<double>& pars) {
  TVector3 fittedVertex(pars[0], pars[1], pars[2]);
  DS::FitVertex vertex;
  vertex.SetPosition(fittedVertex);
  vertex.SetTime(pars[3]);
  DS::FitResult result;
  result.SetVertex(0, vertex);
  return result;
}

bool MultiPathScintWater::ValidResult(const std::vector<double>& pars) {
  double r2 = pars[0]*pars[0] + pars[1]*pars[1] + pars[2]*pars[2];
  double wl = pars[4];
  return ( (r2 < fPSUPRadius2 ) && (wl >= -fAVRadius + fZoff) && (wl <= fAVRadius + fZoff) );
}

std::vector<double> MultiPathScintWater::GetStart() {
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
  double t = CLHEP::RandFlat::shoot( 100.0, 300.0 );
#endif
  double ran2 = CLHEP::RandFlat::shoot( 0.0, 1.0 );
  double wl = 4108 + 2000*ran2;
  //std::cout<<"ran2 wl "<<wl<<std::endl;
  double r = pow(ran0, 1.0/3.0) * fPSUPRadius; // mm
  double costheta = ran1;
  double sintheta = sqrt(1.0 - costheta*costheta);
  std::vector<double> v;
  v.push_back(r*sintheta*cos(ran2Pi));
  v.push_back(r*sintheta*sin(ran2Pi));
  v.push_back(r*costheta+fZoff);
  // in principle, MC and data could have different PMT hit time offset,
  // but in latest incarnation of MultiPathFitter code, this doesn't seem to be the case
  v.push_back(t);
  v.push_back(wl);// water level
  return v;
}

/// returns Likelihood and derivative for a given time of flight t.
double MultiPathScintWater::LAnddLdt( double tDiff, double tof, double &dLdt ) {
  double L = 0.0;
  fRes = tDiff - tof;
  if( fRes > fMaxPMTtRes ) fRes = fMaxPMTtRes; // longer than valid time (late light, crosstalk?)
  if( fRes < -fMaxPMTtRes ) fRes = -fMaxPMTtRes; // eariler than valid time

  int ibin = (int) floor( ( fRes + 100.0 ) * 4.0 );

  if(ibin > 0 && ibin < fEntriesTimeCh - 2) {
    double dt = fRes-fPDFxCh[ibin];
    L = (fPDFCh[ibin] + fDerivativeCh[ibin] * dt);
    if (dt < 0.125) {
      double u = 0.5 + dt * 4.0;
      dLdt = -u*fDerivativeCh[ibin] - (1.0-u)*fDerivativeCh[ibin-1];
    } else {
      double u = 4.0 * dt - 0.5;
      dLdt = -(1.0-u)*fDerivativeCh[ibin] - u*fDerivativeCh[ibin+1];
    }// end of the if statement if dt<0.125
  } else if (ibin <= 0) {
    L = fPDFCh[0];
    dLdt = -fDerivativeCh[0];
  } else {
    L = fPDFCh[fEntriesTimeCh - 2];
    dLdt = -fDerivativeCh[fEntriesTimeCh - 2];
  }
  return L;
}

double MultiPathScintWater::LAnddLdtScint( double tDiff, double tof, double &dLdt ) {
  double L = 0.0;
  fRes = tDiff - tof;
  if( fRes > fMaxPMTtRes ) fRes = fMaxPMTtRes; // longer than valid time (late light, crosstalk?)
  else if( fRes < -fMaxPMTtRes ) fRes = -fMaxPMTtRes; // eariler than valid time

  int ibin = (int) floor( ( fRes + 100.0 ) * 4.0 );
  if(ibin > 0 && ibin < fEntriesTimeCh - 2) {
    double dt = fRes - fPDFxScint[ibin];
    L = (fPDFScint[ibin] + fDerivativeScint[ibin] * dt);
    if (dt < 0.125) {
      double u = 0.5 + dt * 4.0;
      dLdt = -u*fDerivativeScint[ibin] - (1.0-u)*fDerivativeScint[ibin-1];
    } else {
      double u = 4.0 * dt - 0.5;
      dLdt = -(1.0-u)*fDerivativeScint[ibin] - u*fDerivativeScint[ibin+1];
    }// end of the if statement if dt<0.125
  } else if (ibin <= 0) {
    L = fPDFScint[0];
    dLdt = -fDerivativeScint[0];
  } else {
    L = fPDFScint[fEntriesTimeCh- 2];
    dLdt = -fDerivativeScint[fEntriesTimeCh - 2];
  }
  return L;
}

double MultiPathScintWater::Calculate(const std::vector<double>& pmt,
                                      const std::vector<double>& par,
                                      std::vector<double>& diff) {
  // look up position each time - is this efficient enough?
  // const TVector3& pos = DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
  // double tpmt = pmt.GetTime();

  double dx = pmt[0] - par[0];
  double dy = pmt[1] - par[1];
  double dz = pmt[2] - par[2];
  double dr = sqrt(dx*dx + dy*dy + dz*dz); // | pmtpos - X0|

  TVector3 fVertex( par[0], par[1], par[2] ); // X0 = (x0,y0,z0), trial vertex
  TVector3 fPMTpos( pmt[0], pmt[1], pmt[2] );
  /// following is to calculate time residual by: PMTtime - tof - trial time;
  // fRes = pmt[3] - tof - par[3];
  // in partial fitter, tof is calculated separately in different ways in different situations;
  // instead of fRes, (tDiff, tof) is passed to the likelihood calculation
  double tDiff = pmt[3] - par[3];// tDiff = hitTime - trialTime; then fRes = tDiff - tof

  double L = 0.0;
  double dLdt = 0.0, dLdx = 0.0, dLdy = 0.0, dLdz = 0.0, dLdt0 = 0.0, dLdwl = 0.0;
  double fx = par[0]; // trial x0
  double fy = par[1]; // trial y0
  double fz = par[2]; // trial z0
  /// water level set relative to the AV
  double fWaterLevel = par[4]; // trial waterLevel
  fDepth = fWaterLevel - fz; // L - z0
  fHeight = pmt[2] - fWaterLevel; // Zpmt - L
  //std::cout<<"trial water: "<<fWaterLevel<<std::endl;
  if( fWaterLevel >= -fAVRadius + fZoff && fWaterLevel <= fAVRadius + fZoff ) { // valid water level
    double tof = 0.0;

    /// evaluate the distanceInScint = scintpathInNeck + scintpathInAV
    double distInScint = 0;

    /// evaluate line-AV sphere intersection first
    double scintpathInAV = 0;
    double udDotVtx = fx*dx/dr + fy*dy/dr + (fz-fZoff)*dz/dr; // (P-X0)/|P-X0| * (X0 - oAV)
    double vtx2 = fx*fx + fy*fy + (fz-fZoff)*(fz-fZoff); // (X0 - oAV)*(X0 - oAV)
    double rVertex = sqrt( vtx2 ); // Roffset = |X0 - oAV|
    double sqrVal = udDotVtx*udDotVtx - vtx2 + fAVRadiusOuter*fAVRadiusOuter;
    double a1 = 0, a2 = 0, aplus = 0, abig = 0, asmall = 0;
    double da3dwl = 0, dspAVdwl = 0;
    if( sqrVal>0 ) { // line passes AV sphere
      /// find the line-sphere intersect points; a1, a2 are the path lengths btw vertex and intersection points
      a1 = -udDotVtx + sqrt(sqrVal);
      a2 = -udDotVtx - sqrt(sqrVal); // a2<a1
      /// find the line-interface plane intersection; a3 is the path length btw vertex and intersection point
      double a3 = -9999; // a3 is not defined if lightPath is parallel to interface
      if( fabs(dz)>fBoundary ) {
	a3 = ( fWaterLevel-fz )*dr/dz; // well-defined a3
        da3dwl = dr/dz;
      }
      if( rVertex<fAVRadiusOuter ) { // vertex inside the AV, equivalent to a1*a2<0
        if ( a1*a2<0 ) {
          aplus = a1; // always has a1>0>a2
          if( a3>0 ) { // a3>0, hit interface plane
            if( fz<fWaterLevel ) { // vertex below
              if( a3<aplus ) { scintpathInAV = aplus - a3; dspAVdwl = - da3dwl; }
            }
            else { // vertex above
              if( a3<aplus ) { scintpathInAV = a3; dspAVdwl = da3dwl; }
              else scintpathInAV = aplus;
            }
          }
          else { // a3<=0, not hit interface plane
            if( fz>=fWaterLevel ) scintpathInAV = aplus; // vertex must above
          }
        }
      } // vertex in AV

      else { // rVertex>=fAVRadiusOuter, vertex in external
        if( a1>0 && a2>0) { // always has a1>a2>0
          abig = a1; // far intersection point
          asmall = a2; // near intersection point
          double zbig = fz + abig*dz/dr; // z position of the intersection points
          double zsmall = fz + asmall*dz/dr;
          if( zsmall>=fWaterLevel && zbig>=fWaterLevel )
            scintpathInAV = abig - asmall;
          if( zsmall<fWaterLevel && zbig>fWaterLevel && a3>0 )
	  { scintpathInAV = abig - a3; dspAVdwl = -da3dwl; }
          if( zsmall>fWaterLevel && zbig<fWaterLevel && a3>0 )
	  { scintpathInAV = a3 - asmall; dspAVdwl = da3dwl; }
        } // ensure a1 and a2 are positive
      } // vertex in external
    } // pass through AV

    /// evaluate line-neck cylinder intersection
    double scintpathInNeck = 0;
    if( fNeckpathEnable ) { // turn off to ignore high z(neck) events

      double fZneckHi = fPSUPRadius; // top Z of the neck = Zpsup
      double bneck = fx*dx + fy*dy;
      double drPerp2 = dx*dx + dy*dy;
      double sqrValneck = bneck*bneck - ( fx*fx+fy*fy-fNeckRadiusOuter*fNeckRadiusOuter )*drPerp2;

      if( sqrValneck > 0 ) { // check the line passes through the neck
        double aneck1 = dr*( -bneck + sqrt(sqrValneck) )/drPerp2; // dr comes from direction.Unit()
        double aneck2 = dr*( -bneck - sqrt(sqrValneck) )/drPerp2;
        /// the z position of intersection point is checked, to cancel out the duplications in AV

        if( aneck1*aneck2<0 ) { // vertex inside cylinder, one hit point on cylinder
          double aneckplus = aneck1; // aneck1>0>aneck2
          double zplus = fz + aneckplus*dz/dr;
          if( zplus < fZneckHi && zplus > fZneckLo ) { // hit point in neck region
            if( rVertex>=fAVRadiusOuter ) { // vertex outside AV
              scintpathInNeck = aneckplus;
            }
            else { // vertex inside AV
              if( a1*a2<0 ) scintpathInNeck = aneckplus - aplus; // hit AV at one point and then hit neck
            }
          }
          if( zplus < fZneckLo ) { // hit below the neck but on cylinder
            if( rVertex>=fAVRadiusOuter && fz>fZneckLo && fz<fZneckHi ) { // this case forces vertex inside neck and outside AV
              if( a1>0 && a2>0 ) { // must hits the AV at 2 points
                scintpathInNeck = asmall;
              }
            }
          }
        } // vertex inside cylinder

        if( aneck1>0 && aneck2>0 ) { // vertex outside cylinder, two intersection points on cylinder
          double aneckbig = aneck1>aneck2? aneck1:aneck2; // far hit point
          double anecksmall = aneck1<aneck2? aneck1:aneck2; // near hit point
          double zbig = fz + aneckbig*dz/dr; // z of far hit point
          double zsmall = fz + anecksmall*dz/dr; // z of near hit point
          /// if two intersection points are in neck region
          if( zbig<fZneckHi && zbig>fZneckLo && zsmall<fZneckHi && zsmall>fZneckLo )
          {
            scintpathInNeck = aneckbig - anecksmall;
            // an extreme condition: pass thru neck as well as virtual AV sphere inside neck
            if(a1>0 && a2>0) {
              double zAVbig = fz + abig*dz/dr, zAVsmall = fz + asmall*dz/dr; // z of AV hit points
              if(zAVbig>=fZneckLo && zAVsmall>=fZneckLo) scintpathInNeck = scintpathInNeck - (abig-asmall);
            }
          }
          /// if one intersection point below neck while the other in neck
          if( zbig<fZneckLo && zsmall>fZneckLo && zsmall<fZneckHi ) {
            if(a1>0 && a2>0) { // must pass the AV at 2 points and outside the AV
              scintpathInNeck = asmall - anecksmall;
            }
          }
          if( zsmall<fZneckLo && zbig>fZneckLo && zbig<fZneckHi) {
            if(a1*a2<0) {// inside AV
              scintpathInNeck = aneckbig - aplus;
            }
            if(a1>0 && a2>0) { // outside AV
              scintpathInNeck = aneckbig - abig;
            }
          }
        } // vertex outside the neck, aneck1>0 and aneck2>0
      } // line hits the neck
    } // enalbe neck path calculation

    #ifdef DEBUG
    debug<<"scint path: "<<scintpathInAV<<" + "<<scintpathInNeck<<" = "<<scintpathInAV+scintpathInNeck<<"\n";
    #endif

    if( scintpathInAV<0 ) scintpathInAV = 0; // avoid negative calculated value
    if( scintpathInNeck<0 ) scintpathInNeck = 0; // avoid negative calculated value

    distInScint = scintpathInAV + scintpathInNeck;
    tof = ( dr - distInScint)/fSpeedOfLightWater + distInScint/fSpeedOfLightScint;

    /// check whether vertex is in scintillator or water
    bool checkVertexInScint = false;
    if ( (rVertex<fAVRadiusOuter && fDepth<0) || (fz>=fZneckLo && fz<fPSUPRadius && fVertex.Perp()<=fNeckRadiusOuter) )
      checkVertexInScint = true;

    int layout = 0; // flag cases
    if( distInScint < 0 ) distInScint = 0; // remove wrong calculations

    if( distInScint < 0.01 ) { // if light path is always in water; except for geo 1, the other layout go to default
      layout = 0;
    }
    else layout = 1;

    switch ( layout ) { // the layout is related to PMT and vertex; is kept to introduce more complicated situations
      case 0: // default, MPW like: direct light path always in water and never pass AV
      {
        tof = dr/fSpeedOfLightWater;
        L = LAnddLdt(tDiff, tof, dLdt);
        dLdt0 = dLdt/L;
        double factor = -dLdt0/(dr*fSpeedOfLightWater);
        dLdx = factor * dx;
        dLdy = factor * dy;
        dLdz = factor * dz;
        dLdwl = 1e-12;
	L = log(L);
        break;
      }

      case 1: // light path passes through scintillator; vertex in water or scintillator
      {
        // if( checkVertexInScint ) L = LAnddLdtScint(tDiff, tof, dLdt);
        // else L = LAnddLdt(tDiff, tof, dLdt);
        L = LAnddLdtScint(tDiff, tof, dLdt);
	double fSpeedOfLightEff = checkVertexInScint? fSpeedOfLightScint:fSpeedOfLightWater;
        dLdt0 = dLdt/L;
        /// derivatives are simplified by moving the vertex a little bit in the opposite direction to the light path in water or scintillator
        double factor = -dLdt0/(dr*fSpeedOfLightEff);
        dLdx = factor * dx;
        dLdy = factor * dy;
        dLdz = factor * dz;
	dLdwl = dLdt0*( -dspAVdwl/fSpeedOfLightWater + dspAVdwl/fSpeedOfLightScint );
        L = log(L);
        break;
      }
    } // switch layout

    diff.assign(5, 0.0);// sBeta
    diff[0] = dLdx;
    diff[1] = dLdy;
    diff[2] = dLdz;
    diff[3] = dLdt0;
    diff[4] = dLdwl;
  } // valid water level
   else {
    fHeight = pmt[2];
    // Log::Die( "Not a valid water height ");
   }

#ifdef DEBUG
  debug<< "AWP::CL fX= (" << par[0] << "," << par[1] << "," << par[2] << ") , "<< par[3]<<" \n";
  debug<< "AWP::CL pmt= (" << pmt[0] << "," << pmt[1] << "," << pmt[2] << ")" << pmt[3]<<"\n";
  debug<< "        fRes=" << fRes << " ," << dLdx << " "<< dLdy << " " <<dLdz << "\n";
  debug<< "        L=" << L << " diff=(" << diff[0] << "," << diff[1] << "," << diff[2] << "," << diff[3] << ")\n";
#endif

  return L;
}
