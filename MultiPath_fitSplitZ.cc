#include <cmath>
#include <string>
#include <vector>
#include <ostream>
#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <RAT/DS/FitVertex.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/MultiPath.hh>
#include <RAT/MultiPathWaterPosition.hh>
#include <RAT/MultiPathWaterDirection.hh>
#include <RAT/MultiPathScintWater.hh>
#include <RAT/MultiPathScint.hh>
#include <RAT/SDecompQRH.hh>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TError.h" // to set gErrorIgnoreLevel

//#define DEBUG

namespace RAT {

namespace Methods {

MultiPath::MultiPath() {
  fNEvents = 0;
  fFunction = 0;
  //gErrorIgnoreLevel = kSysError; // higher level is kFatal
}

MultiPath::~MultiPath() {
  if (fFunction) delete fFunction;
}

void MultiPath::SetFunction(VertexFunction* func) {
  fFunction = func;
  size_t n = fFunction->GetNumberOfParameters();
  fBeta.ResizeTo(n);
  fCov.ResizeTo(n,n);
}

// initialize with multipath fitter for different phases
// (default is multipath-waterposition)
void MultiPath::Initialise(const std::string& param) {
  fInitString = param;
  if (param == "waterdirection") SetFunction(new RAT::MultiPathWaterDirection());
  else if (param == "waterposition") {
    SetFunction(new RAT::MultiPathWaterPosition());
    fInitString = "waterposition";
  }
  else if (param == "scintwater") {
    SetFunction(new RAT::MultiPathScintWater());
    fInitString = "scintwater";
  }
  else if (param == "scint") {
    SetFunction(new RAT::MultiPathScint());
    fInitString = "scint";
  }

  fFunction->Initialise(param);
}

void MultiPath::BeginOfRun(DS::Run& run) {
  if (fFunction) fFunction->BeginOfRun(run);

  DB* db0 = DB::Get();
  DBLinkPtr dbLinkFit = db0->GetLink( "FIT_MULTIPATH" );
  //number of steps to do after reaching tolerance
  fStepsBelowTolerance = dbLinkFit->GetD( "steps_below_tolerance" );
  //convergence parameter
  fMaximumIterations = dbLinkFit->GetD( "maximum_iterations" );
  fTolerance = dbLinkFit->GetD( "fit_tolerance" ); // chisquare tolerance
  fMaxGoodFits = dbLinkFit->GetD( "maximum_goodfits" ); // max good fits
  fBoundaryCut = dbLinkFit->GetZ( "boundaryCut_valid" );
  fDumpValid = dbLinkFit->GetZ( "dump_valid" );
  //interested events with given gtid
  fGTID = dbLinkFit->GetIArray( "dump_gtid" );
  std::sort(fGTID.begin(),fGTID.end());

  if (fInitString == "scintwater") { // extra parameters needed by partial fitter
    DBLinkPtr innerAvDb = db0->GetLink( "GEO", "inner_av" );
    fFitterPar = innerAvDb->GetD( "split_z" ); // a parameter needed by specific fitter, eg. fWaterLevel for partial fitter
    fMax_nStart = dbLinkFit->GetD( "max_start_positions" );
    fBoundaryTolerance = dbLinkFit->GetD( "boundary_tolerance" );
    fZoff = dbLinkFit->GetD( "av_offset_z" ); //z offset of AV to PSUP center
    fAVRadius = db0->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_sphere" ); // 6005.0, AV inner radius
    fAVRadiusOuter = db0->GetLink( "SOLID", "acrylic_vessel_outer" )->GetD( "r_sphere" ); // 6060, AV outer radius
    fNeckRadius = db0->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_neck" ); // inner radius of neck = 730 mm
    fNeckRadiusOuter = db0->GetLink( "SOLID", "acrylic_vessel_outer" )->GetD( "r_neck" ); // outer radius of neck = 785 mm
    fZneckLo = sqrt( fAVRadiusOuter*fAVRadiusOuter - fNeckRadiusOuter*fNeckRadiusOuter ) + fZoff;
  }
}

void MultiPath::EndOfRun(DS::Run& run) {
  if (fFunction) fFunction->EndOfRun(run);
}

void MultiPath::DefaultSeed() {
  DS::FitVertex vertex;
  // Initialise the SeedResult at the centre and with arbitrary
  // 3000 mm errors in all directions (note units are mm and ns)
  const TVector3 seedPosition( 0.0, 0.0, 0.0 );
  const ROOT::Math::XYZVectorF seedPositionError( 3000.0, 3000.0, 3000.0 );
  // 230 ns is roughly the peak from ET1D and GV1D PDFs
  const double seedTime = 230.0;
  const double seedTimeError = 100.0;
  vertex.SetPosition( seedPosition, false, true );
  vertex.SetPositionErrors( seedPositionError, false, true );
  vertex.SetTime( seedTime, false, true );
  vertex.SetTimeErrors( seedTimeError, false, true );
  fSeedResult.SetVertex( 0, vertex );
}

double MultiPath::SumLikelihoods(const std::vector<double>& par) {
  size_t npars = par.size();
  std::vector<double> beta(npars, 0.0);
  std::vector<double> cov(npars*(npars+1)/2, 0.0);

  std::vector<double> dbeta(npars);
  double L = 0.0;
  for (size_t ipmt = 0; ipmt < fPMTReduced.size(); ipmt++) {
    L += fFunction->Calculate(fPMTReduced[ipmt], par, dbeta);
    for (size_t i = 0; i < npars; i++) beta[i] += dbeta[i];
    size_t k = 0;
    for (size_t i = 0; i < npars; i++) {
      for (size_t j = i; j < npars; j++) {
        cov[k++] += dbeta[i]*dbeta[j];
      }
    }
  }

  for (size_t i = 0; i < npars; i++) fBeta[i] = beta[i];
  size_t k = 0;
  for (size_t i = 0; i < npars; i++) {
    for (size_t j = i; j < npars; j++) {
      fCov[i][j] = cov[k++];
      if (i != j) fCov[j][i] = fCov[i][j];
    }
  }
  return L;
}

DS::FitResult MultiPath::GetBestFit() {
#ifdef DEBUG
  debug << "MultiPath::GetBestFit() called (gtid " << fEvent->GetGTID() << ")\n";
#endif

  // get seed
  fFitResult.Reset();
  if (fPMTData.empty()) return fFitResult;
  CopySeedToResult();
  fFunction->SetSeed(fFitResult);

  SelectPMTData(fFitResult.GetVertex(0));
  if (fSelectedPMTData.empty()) return fFitResult;

  // copy pmt data after PMT selection
  fPMTReduced.clear();
  for (size_t i = 0; i < fSelectedPMTData.size(); i++) {
    const RAT::FitterPMT& pmt = fSelectedPMTData[i];
    const TVector3& pos = DU::Utility::Get()->GetPMTInfo().GetPosition(pmt.GetID());
    std::vector<double> data;
    data.push_back(pos.X());
    data.push_back(pos.Y());
    data.push_back(pos.Z());
    data.push_back(pmt.GetTime());
    fPMTReduced.push_back(data);
  }

#ifdef DEBUG
  debug << "MultiPath::GetBestFit:  " << fPMTReduced.size() << " PMTs\n";
#endif

  int nGoodFits = 0;
  int nStart = 0;
  fLBest = -100000.0;
  size_t npars;
  for (nStart = 0; nGoodFits < fMaxGoodFits && nStart < 250; nStart++) {
    int iStepsBelowTolerance = 0;
    double alambda = 0.001;
    bool oneStep = false; // set to true the first time we find a real step

    // set initial parameters
    std::vector<double> pars = fFunction->GetStart();
    npars = pars.size();

#ifdef DEBUG
    debug << "Seed " << nStart << " (";
    for (size_t i = 0; i < npars; i++) debug << " " << pars[i];
    debug << ")\n";
#endif

    double L = SumLikelihoods(pars);

    std::vector<double> parsSaved = pars;
    double LSaved = L;
    TVectorD betaSaved = fBeta;
    TMatrixT<double> covSaved = fCov;

    for (int iteration = 0; iteration < fMaximumIterations; iteration++) {

#ifdef DEBUG
      debug << "    step (" << iteration << " " << iStepsBelowTolerance << " " << alambda
             << " " << L << ")";
      for (size_t i = 0; i < npars; i++) debug << " " << pars[i];
      debug << ")\n";
#endif

      // augment diagonal elements to implement Marquardt
      for (size_t i = 0; i < npars; i++) fCov[i][i] *= (1.0 + alambda);

      TVectorD dA = fBeta;
      if (npars > 2) {
#ifdef DEBUG
        for (size_t i = 0; i < npars; i++) {
          debug << "      A" << i << ":  ";
          for (size_t j = 0; j < npars; j++) {
            debug << fCov[i][j] << " ";
          }
          debug << "\n";
        }
        debug << "       b=(";
        for (int i = 0; i < npars; i++) debug << " " << dA[i];
        debug << ")\n";
#endif
        SDecompQRH lu(fCov);
        if (!lu.Decompose()) {
#ifdef DEBUG
          debug << "    LU non-invertible\n";
#endif
          break; // break out if not invertible
        }
        if (!lu.Solve(dA)) {
#ifdef DEBUG
          debug << "    cannot solve LU system\n";
#endif
          break; // break out if not invertible
        }
#ifdef DEBUG
        debug << "       x0=" << dA[0] << " x1=" << dA[1] << " x2=" << dA[2] << " x3=" << dA[3] <<  "\n";
#endif
      } else {
        // 2x2 problem:  quicker to do it directly?
        double a = fCov[0][0];
        double b = fCov[0][1];
        double c = fCov[1][1];
        double det = a*c - b*b; // symmetric
#ifdef DEBUG
        debug << "      A00=" << a << " A01=" << b << " A11=" << c << " det=" << det << "\n";
        debug << "       b0=" << dA[0] << " b1=" << dA[1] << "\n";
#endif
        if (std::fabs(det) < 1e-16) break; // machine epsilon for double
        double x = (c*dA[0] - b*dA[1]) / det;
        double y = (a*dA[1] - b*dA[0]) / det;
        dA[0] = x;
        dA[1] = y;
#ifdef DEBUG
        debug << "       x0=" << x << " x1=" << y << "\n";
#endif
      }

      // likelihood still not minimized: proceed
      bool zero = true;
      for (size_t j = 0; j < pars.size(); j++) {
        pars[j] += dA[j];
        if (dA[j] != 0.0) zero = false;
      }
      if (zero) {
#ifdef DEBUG
        debug << "      zero step\n";
#endif
        break; // no step at all - otherwise infinite loop
      }
      L = SumLikelihoods(pars); // overwrites fCov and fBeta

      // test if update is small
#ifdef DEBUG
      debug << "      dL = " << L-LSaved << " oneStep=" << oneStep << "\n";
#endif
      if (oneStep && std::fabs(L - LSaved) < fTolerance) {
        iStepsBelowTolerance++;
        if (iStepsBelowTolerance >= fStepsBelowTolerance) break;
      } else {
        iStepsBelowTolerance = 0;
      }

      // test if update is better
      if (L >= LSaved) {
        oneStep = true;
        // if likelihood greater than previous likelihood, alambda multiplied by 0.1
        alambda *= 0.1;
        parsSaved = pars;
        LSaved = L;
        covSaved = fCov;
        betaSaved = fBeta;
      } else {
        // if likelihood not greater than previous, alambda multiplied by 10
        alambda *= 10.0;
        pars = parsSaved;
        L = LSaved; // reset to previous step data
        fBeta = betaSaved;
        fCov = covSaved;
        // don't test whether L step is small; even if it is, we're not optimized
      }
      if (alambda > 1.0e90) {
#ifdef DEBUG
        debug << "    big alambda break\n";
#endif
        break;
      }

    } // for 500 or 100 iterations

#ifdef DEBUG
    debug << "    done with iterations L=" << L << " LBest=" << fLBest << " (";
    for (size_t jk = 0; jk < npars; jk++) debug << pars[jk] << " ";
    debug << ")\n";
#endif

    if (L > fLBest-0.5 && fFunction->ValidResult(pars)) {
      if (fLBest < L-0.25) nGoodFits = 1; // reset good fit counter (significant change)
      else nGoodFits++; // increments good fit counter (need 6 to have a valid fit)
      if (L > fLBest) {
        // reset best likelihood to current likelihood value
        fParsBest = pars;
        fLBest = L;
      }
#ifdef DEBUG
      debug << "  Good fit (" << nGoodFits << " " << nStart << ") " << fLBest << " (";
      for (size_t jk = 0; jk < npars; jk++) debug << fParsBest[jk] << " ";
      debug << ")\n";
#endif
    }

  } // start position iterations

  if (nGoodFits == 0) fParsBest.assign(npars, -10000.0);
  // at this point, may calculate other quantities (CalculatePMTOther)

#ifdef DEBUG
  debug << "  Best fit (" << nGoodFits << " " << nStart << ") " << fLBest << " (";
  for (size_t jk = 0; jk < npars; jk++) debug << fParsBest[jk] << " ";
  debug << ")\n";
#endif

  if (fDumpValid && std::find(fGTID.begin(), fGTID.end(), fEvent->GetGTID()) != fGTID.end()) {
    DumpLikelihood(fNEvents, fEvent->GetGTID(), fInitString);
  }
  fNEvents++;
  double fwaterlevel = 0;
  fwaterlevel = fParsBest[4];
  // copy into fit result for return
  fFitResult = fFunction->MakeResult(fParsBest);
  fFitResult.SetFOM("multipath_" + fInitString, fLBest);
  fFitResult.SetFOM("multipath_SelectedNHit_" + fInitString, static_cast<double>( fPMTReduced.size() ));
  fFitResult.SetFOM("starts", static_cast<double>(nStart));
  fFitResult.SetFOM("goodfits", static_cast<double>(nGoodFits));
  fFitResult.SetFOM("fitted_splitZ", static_cast<double>(fwaterlevel));
  return fFitResult;
}

void MultiPath::DumpLikelihood(int index, int gtid, std::string dumpVariable) {
  // make histograms
  TH1F LLHFitPosition("LLHFitPosition","LLHFitPosition", 10000, 0.0, 10000.0);
  TH1F LLHoriginalPosition("LLHoriginalPosition","LLHoriginalPosition", 10000, 0.0, 10000.0);
  TH1F TresFitPosition("TresFitPosition","TresFitPosition", 10000, 0.0, 10000.0);
  TH1F TresOriginalPosition("TresOriginalPosition","TresOriginalPosition", 10000, 0.0, 10000.0);
  TH2F pmtGEO("pmtGEO","PMT Geometry", 628, -M_PI, M_PI, 1000, -1.0, 1.0);
  TH2F LZR("LZR","likelihood surface", 2000, -9000, 9000, 1000,0,9000);

  std::vector<TH1F*> vHists;
  const size_t nprefix = 3;
  const char* prefix[] = { "H", "DH", "NDH" };
  const std::vector<VertexFunction::Parameter>& specs = fFunction->GetParameters();
  size_t npars = specs.size();
  for (size_t ip = 0; ip < nprefix; ip++) {
    for (size_t in = 0; in < npars; in++) {
      std::ostringstream s;
      s << prefix[ip] << specs[in].GetName().c_str() << "_" << index;
      const char* sn = s.str().c_str();
      TH1F* h = new TH1F(sn, sn, specs[in].GetNBins(),
                         specs[in].GetLowValue(), specs[in].GetHighValue());
      vHists.push_back(h);
    }
  }

  // fill histograms
  for (size_t ipmt = 0; ipmt < fPMTReduced.size(); ipmt++) {
    std::vector<double> dbeta(npars);
    double L = fFunction->Calculate(fPMTReduced[ipmt], fParsBest, dbeta);
    double tres = fFunction->GetResidual(); // GetResidual
    int PMTid = fSelectedPMTData[ipmt].GetID();
    LLHFitPosition.Fill(PMTid, L); // likelihood contribution of each pmt
    TresFitPosition.Fill(PMTid, tres);
    const TVector3& pos = DU::Utility::Get()->GetPMTInfo().GetPosition(PMTid);
    pmtGEO.Fill(pos.Phi(), pos.CosTheta()); // pmt locations
  }

  // re-evaluate with trial parameters set to 0?
  std::vector<double> origin(npars, 0.0);
  //if (npars >= 3) origin[3] = fParsBest[3];
  for (size_t ipmt = 0; ipmt < fPMTReduced.size(); ipmt++) {
    std::vector<double> dbeta(npars);
    double L = fFunction->Calculate(fPMTReduced[ipmt], origin, dbeta);
    double tres = fFunction->GetResidual(); // GetResidual
    int PMTid = fSelectedPMTData[ipmt].GetID();
    TresOriginalPosition.Fill(PMTid, tres);
    LLHoriginalPosition.Fill(PMTid, L);
  }

  // use best parameters
  std::vector<double> trial(npars, 0.0);
  for (size_t j = 0; j < npars; j++) {
   trial[j] = fParsBest[j];
  }

  // varying one param. while keeping the others fixed
  for (size_t j = 0; j < npars; j++) {
    int count = 0;
    double lowBin = vHists[j]->GetBinLowEdge(1);
    double highBin = vHists[j]->GetBinLowEdge(vHists[j]->GetNbinsX()+1);
    double increment = lowBin - vHists[j]->GetBinLowEdge(0);
    int NIncrement = (int)(highBin-lowBin)/increment;
    double oldL = 0.0;
    for (int k = 0; k < NIncrement; k++) {
      //trial[j] = fParsBest[j] + (lowBin + k*increment);
      trial[j] = (lowBin + k*increment);
      double L = SumLikelihoods(trial);
      vHists[j]->Fill(trial[j], L);
      vHists[j+npars]->Fill(trial[j], fBeta[j]); // fill fitter calculated derivatives (DHx, DHy, ...)
      if (count++ == 0) oldL = L; // for numeric derivative, initialize the first bin
      if (count) {
        double w = (L - oldL) / (int)((increment + 0.005)*100.0);
        vHists[j+2*npars]->Fill(trial[j], w*100); // fill numerical derivatives (NDHx, NDHy, ...)
        oldL = L;
      }
    }
    trial[j] = fParsBest[j];// set varing param. back to the best fit values and varying the next one
  }

  // write histograms
  std::ostringstream s;
  s << "multipath_" << fInitString << "_gtID_" << gtid << "_" << dumpVariable << ".root";
  TFile file(s.str().c_str(), "RECREATE");
  LLHFitPosition.Write();
  LLHoriginalPosition.Write();
  TresFitPosition.Write();
  TresOriginalPosition.Write();
  pmtGEO.Write();
  for (size_t i = 0; i < vHists.size(); i++) {
    TH1F* h = vHists[i];
    h->Write();
  }
  file.Close();

  // delete histograms
  for (size_t i = 0; i < vHists.size(); i++) {
    TH1F* h = vHists[i];
    delete h;
  }
}

} //::Methods

} //::RAT

