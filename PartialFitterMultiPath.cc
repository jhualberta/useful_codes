#include <RAT/PartialFitterMultiPath.hh>
#include <RAT/MethodFactory.hh>
#include <RAT/PDFFactory.hh>
#include <RAT/OptimiserFactory.hh>
#include <RAT/ClassifierFactory.hh>
#include <RAT/PMTSelectorFactory.hh>
#include <RAT/Method.hh>
#include <RAT/Classifier.hh>
#include <RAT/OptimisedMethod.hh>
#include <RAT/PDFMethod.hh>
#include <RAT/SeededMethod.hh>
#include <RAT/SelectorMethod.hh>
#include <RAT/SelectorClassifier.hh>
#include <RAT/Optimiser.hh>
#include <RAT/SeededClassifier.hh>
#include <RAT/OptimisedClassifier.hh>
#include <RAT/FitterPMT.hh>
#include <RAT/MetaInformation.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <RAT/ListHelp.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DS/Entry.hh>
using namespace RAT;
using namespace RAT::Methods;
using namespace RAT::PDFs;
using namespace RAT::Optimisers;
using namespace RAT::PMTSelectors;
using namespace RAT::DS;
using namespace RAT::Classifiers;

#include <sstream>
using namespace std;

#include <TStopwatch.h>
using namespace ROOT;

PartialFitterMultiPath::PartialFitterMultiPath()
  : Processor("partialFitterMP")
{
  fQuadSeed = MethodFactory::Get()->GetMethod( "quad" );
  fPositionTime = MethodFactory::Get()->GetMethod( "multipath-scintwater" );
  fNullSelector = PMTSelectorFactory::Get()->GetPMTSelector( "null" );
  fPMTCalSelector = PMTSelectorFactory::Get()->GetPMTSelector( "PMTCalSelector" );
  fModeCut = PMTSelectorFactory::Get()->GetPMTSelector( "modeCut" );
  fModeCut->SetD( "lowLimit", -100.5 );
  fModeCut->SetD( "highLimit", 100.5 );

  fTimeResidualCut = PMTSelectorFactory::Get()->GetPMTSelector( "timeResidualCut" );
  // fET1D = PDFFactory::Get()->GetPDF( "partialET1D" );
  fEnergySeed = MethodFactory::Get()->GetMethod( "energyRThetaFunctional" );
  fPartialEnergy = MethodFactory::Get()->GetMethod( "partialEnergy" );
  fIsotropy = ClassifierFactory::Get()->GetClassifier( "isotropy" );
  fITR = ClassifierFactory::Get()->GetClassifier( "ITR" );
  fQPDT = ClassifierFactory::Get()->GetClassifier( "QPDT" );
  fTimingPeaks = ClassifierFactory::Get()->GetClassifier( "timingPeaks" );
  fMeanTime = ClassifierFactory::Get()->GetClassifier( "meanTime" );
  fPreTriggerHits = ClassifierFactory::Get()->GetClassifier( "preTriggerHits" );
  fIsoRegions = ClassifierFactory::Get()->GetClassifier( "isoRegions" );
  fEarlyTime = ClassifierFactory::Get()->GetClassifier( "earlyTime" );
  fBiPoCumulTimeResid = ClassifierFactory::Get()->GetClassifier( "BiPoCumulTimeResid" );
  fBiPo212LikelihoodDiff = ClassifierFactory::Get()->GetClassifier( "BiPoLikelihoodDiff-212" );
  fBiPo214LikelihoodDiff = ClassifierFactory::Get()->GetClassifier( "BiPoLikelihoodDiff-214" );
  fBiPoLikelihoodDiffOptimiser = OptimiserFactory::Get()->GetOptimiser( "grid-55" );
  fBerkeleyAlphaBeta = ClassifierFactory::Get()->GetClassifier("BerkeleyAlphaBeta");
  fNearAVAngular = ClassifierFactory::Get()->GetClassifier( "nearAVAngular" );
  // Now combine the components where appropriate
  dynamic_cast< SelectorMethod* >( fQuadSeed )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorMethod* >( fQuadSeed )->AddPMTSelector( fModeCut);
  // dynamic_cast< OptimisedMethod* >( fPositionTime )->SetOptimiser( fPowell );
  dynamic_cast< SelectorMethod* >( fPositionTime )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorMethod* >( fPositionTime )->AddPMTSelector( fModeCut );
  // dynamic_cast< PDFMethod* >( fPositionTime )->SetPDF( fET1D );
  dynamic_cast< SelectorClassifier* >( fIsotropy )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorClassifier* >( fIsotropy )->AddPMTSelector( fTimeResidualCut );
  dynamic_cast< SelectorClassifier* >( fITR )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorClassifier* >( fQPDT )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< OptimisedClassifier* >( fBiPo212LikelihoodDiff )->SetOptimiser( fBiPoLikelihoodDiffOptimiser );
  dynamic_cast< OptimisedClassifier* >( fBiPo214LikelihoodDiff )->SetOptimiser( fBiPoLikelihoodDiffOptimiser );

  fCutOff = 3.0; // Maximum nhit for which quad cannot run
}

PartialFitterMultiPath::~PartialFitterMultiPath()
{
  delete fQuadSeed;
  delete fPositionTime;
  // delete fPowell;
  delete fNullSelector;
  // delete fET1D;
  delete fEnergySeed;
  delete fPartialEnergy;
  delete fIsotropy;
  delete fITR;
  delete fQPDT;
  delete fTimingPeaks;
  delete fMeanTime;
  delete fPreTriggerHits;
  delete fIsoRegions;
  delete fEarlyTime;
  delete fBiPoCumulTimeResid;
  delete fBiPo212LikelihoodDiff;
  delete fBiPo214LikelihoodDiff;
  delete fBiPoLikelihoodDiffOptimiser;
  delete fBerkeleyAlphaBeta;
}

void
PartialFitterMultiPath::BeginOfRun( DS::Run& run )
{
  /// First call the begin of run functions
  fQuadSeed->BeginOfRun( run );
  fPositionTime->BeginOfRun( run );
  // fPowell->BeginOfRun( run );
  fNullSelector->BeginOfRun( run );
  // fET1D->BeginOfRun( run );
  fEnergySeed->BeginOfRun( run );
  fPartialEnergy->BeginOfRun( run );
  fITR->BeginOfRun( run );
  fQPDT->BeginOfRun( run );
  fTimingPeaks->BeginOfRun( run );
  fMeanTime->BeginOfRun( run );
  fPreTriggerHits->BeginOfRun( run );
  fIsoRegions->BeginOfRun( run );
  fEarlyTime->BeginOfRun( run );
  fBiPoCumulTimeResid->BeginOfRun( run );
  fBiPo212LikelihoodDiff->BeginOfRun( run );
  fBiPo214LikelihoodDiff->BeginOfRun( run );
  fBerkeleyAlphaBeta->BeginOfRun( run );
}

Processor::Result
PartialFitterMultiPath::DSEvent( DS::Run& run,
                        DS::Entry& ds )
{
  for( size_t iEV = 0; iEV < ds.GetEVCount(); iEV++ )
    Event( run, ds.GetEV( iEV ) );
  return OK;
}

Processor::Result
PartialFitterMultiPath::Event( DS::Run& run,
                      DS::EV& ev )
{

  /// PartialFitter Logic:
  ///  This is a simple overview of the logic that is expressed in the next 80 lines
  ///
  ///  1. Fit the position & time as the seedResult from  quad. The results are used for energy fitter if fit fails.
  ///     The MP fitter for position & time has its own random seeding and not using quad results as seed.
  ///      if nhits < fNhitCutoff then abort the fit and return no FitResult. \n
  ///  2. Fit the position & time as the partialResult using the MultiPath method and MultiPathScintWater
  ///     specifically.
  ///     If this fails then abort the fit and return no FitResult
  ///  3. Generate energy seed using energyRThetaFunctional and partialResult as seed (or seedResult if the
  ///      positionTimeLikelihood result is invalid)
  ///  4. Fit the energy using partialEnergy as the method, energySeed and partialResult (or seedResult)
  ///      as seeds.
  ///     If this fails then abort the fit and return no FitResult
  ///  5. Run beta14 (isotropy) with partialResult as seed
  ///  6. Run itr with partialResult as seed \n
  ///  7. Run qpdt with partialResult as seed. \n
  ///  8. Run timingPeaks \n
  ///  9. Run meanTime with partialResult as seed \n
  ///  10. Run preTriggerHits \n
  ///  11. Run isoRegions with partialResult as seed \n
  ///  12. Run earlyTime with partialResult as seed \n
  ///  13. Run BiPoCumulTimeResid with partialResult as seed \n
  ///  14. Run BiPoLikelihoodDiff with partialResult as seed and grid-55 as optimiser
  ///      for both 212BiPo and 214BiPo PDFs \n
  ///  15. Run BerkeleyAlphaBeta with partialResult as seed \n

  TStopwatch timer;
  timer.Start( true );

  const size_t currentPass = MetaInformation::Get()->GetCurrentPass();
  ev.AddFitterPass( currentPass ); // Ensure the EV knows fitters were run for this pass
  ev.AddClassifierPass( currentPass ); // Ensure the EV knows classifiers were run for this pass

  vector<FitterPMT> pmtData;
  for( size_t iPMTCal =0; iPMTCal < ev.GetCalPMTs().GetCount(); iPMTCal++ )
    pmtData.push_back( FitterPMT( ev.GetCalPMTs().GetPMT( iPMTCal ) ) );

  /// Initialize the seed
  fQuadSeed->SetEventData( pmtData, &ev, &run );
  // Run the seed method
  FitResult seedResult;
  try
    {
      if( ev.GetCalPMTs().GetCount() > fCutOff )
        seedResult = fQuadSeed->GetBestFit();
      else
        {
          warn << "PartialFitterMultiPath::Event: insufficient points for quad to run, exiting" << newline;
          return Processor::FAIL;
        }
    }
  catch( Method::MethodFailError& error ) { warn << "PartialFitterMultiPath::Event: Seed failed " << error.what() << ", continuing." << newline; }

  // Now initialise the position time method
  fPositionTime->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededMethod* >( fPositionTime )->DefaultSeed();
  // dynamic_cast< SeededMethod* >( fPositionTime )->SetSeed( seedResult );
  /// Run the position time method
  FitResult partialResult;
  try
    {
      FitResult positionResult = fPositionTime->GetBestFit();
      partialResult.SetVertex(0, positionResult.GetVertex(0));
      SetFOMs(partialResult, positionResult, "Position");
    }
  catch( Method::MethodFailError& error )
    {
      warn << "PartialFitter::Event: Main method failed " << error.what() << ", exiting" << newline;
      return Processor::FAIL;
    }

  // Now initialise the energy seed
  FitResult energySeed;
  fEnergySeed->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededMethod* >( fEnergySeed )->DefaultSeed();
  if( partialResult.GetValid() )
    dynamic_cast< SeededMethod* >( fEnergySeed )->SetSeed( partialResult );
  else if( seedResult.GetValid() )
    dynamic_cast< SeededMethod* >( fEnergySeed )->SetSeed( seedResult );
  else
    warn << "PartialFitter::Event: No seed for the energy lookup method." << newline;
  // Run the energy functional method
  try
    {
      energySeed = fEnergySeed->GetBestFit();
    }
  catch( Method::MethodFailError& error )
    {
      warn << "PartialFitter::Event: Energy seed method failed " << error.what() << ", energy fit will use default seed" << newline;
    }

  // Now initialise the energy method
  fPartialEnergy->SetEventData( pmtData, &ev, &run );
  // Seed position
  dynamic_cast< SeededMethod* >( fPartialEnergy )->DefaultSeed();
  if( partialResult.GetValid() )
    dynamic_cast< SeededMethod* >( fPartialEnergy )->SetSeed( partialResult );
  else if( seedResult.GetValid() )
    dynamic_cast< SeededMethod* >( fPartialEnergy )->SetSeed( seedResult );
  else
    warn << "PartialFitter::Event: No seed for the energy lookup method." << newline;
  // Seed energy
  dynamic_cast< SeededMethod* >( fPartialEnergy )->SetSeed( energySeed );
  // Run the energy functional method
  try
    {
      FitVertex vertex = partialResult.GetVertex(0);
      FitResult energyResult = fPartialEnergy->GetBestFit();
      FitVertex energyVertex = energyResult.GetVertex(0);
      vertex.SetEnergy( energyVertex.GetEnergy(), energyVertex.ValidEnergy() );
      vertex.SetPositiveEnergyError( energyVertex.GetPositiveEnergyError(), energyVertex.ValidEnergy() );
      vertex.SetNegativeEnergyError( energyVertex.GetNegativeEnergyError(), energyVertex.ValidEnergy() );
      partialResult.SetVertex( 0, vertex );
      SetFOMs(partialResult, energyResult, "Energy");
    }
  catch( Method::MethodFailError& error )
    {
      warn << "PartialFitter::Event: Energy fit method failed " << error.what() << newline;
    }

  // Now initialise and run the Isotropy classifier
  fIsotropy->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fIsotropy )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fIsotropy )->SetSeed( partialResult );
  try { ev.SetClassifierResult( currentPass, fIsotropy->GetName() + ":partialFitter", fIsotropy->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the ITR classifier
  fITR->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fITR )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fITR )->SetSeed( partialResult );
  try { ev.SetClassifierResult( currentPass, fITR->GetName() + ":partialFitter", fITR->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the QPDT classifier
  fQPDT->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fQPDT )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fQPDT )->SetSeed( partialResult );
  try { ev.SetClassifierResult( currentPass, fQPDT->GetName() + ":partialFitter", fQPDT->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the TimingPeaks classifier
  fTimingPeaks->SetEventData( pmtData, &ev, &run );
  try { ev.SetClassifierResult( currentPass, fTimingPeaks->GetName() + ":partialFitter", fTimingPeaks->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the MeanTime classifier
  fMeanTime->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fMeanTime )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fMeanTime )->SetSeed( partialResult );
  try { ev.SetClassifierResult( currentPass, fMeanTime->GetName() + ":partialFitter", fMeanTime->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the PreTriggerHits classifier
  fPreTriggerHits->SetEventData( pmtData, &ev, &run );
  try { ev.SetClassifierResult( currentPass, fPreTriggerHits->GetName() + ":partialFitter", fPreTriggerHits->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the IsoRegions classifier
  fIsoRegions->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fIsoRegions )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fIsoRegions )->SetSeed( partialResult );
  try { ev.SetClassifierResult( currentPass, fIsoRegions->GetName() + ":partialFitter", fIsoRegions->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the EarlyTime classifier
  fEarlyTime->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fEarlyTime )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fEarlyTime )->SetSeed( partialResult );
  try { ev.SetClassifierResult( currentPass, fEarlyTime->GetName() + ":partialFitter", fEarlyTime->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the BiPoCumulTimeResid classifier
  fBiPoCumulTimeResid->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fBiPoCumulTimeResid )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fBiPoCumulTimeResid )->SetSeed( partialResult );
  try { ev.SetClassifierResult ( currentPass, fBiPoCumulTimeResid->GetName() + ":partialFitter", fBiPoCumulTimeResid->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the BiPoLikelihoodDiff classifiers
  fBiPo212LikelihoodDiff->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fBiPo212LikelihoodDiff )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fBiPo212LikelihoodDiff )->SetSeed( partialResult );
  try { ev.SetClassifierResult ( currentPass, fBiPo212LikelihoodDiff->GetName() + ":212:partialFitter", fBiPo212LikelihoodDiff->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }
  fBiPo214LikelihoodDiff->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fBiPo214LikelihoodDiff )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fBiPo214LikelihoodDiff )->SetSeed( partialResult );
  try { ev.SetClassifierResult ( currentPass, fBiPo214LikelihoodDiff->GetName() + ":214:partialFitter", fBiPo214LikelihoodDiff->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the BerkeleyAlphaBeta classifier
  fBerkeleyAlphaBeta->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fBerkeleyAlphaBeta )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fBerkeleyAlphaBeta )->SetSeed( partialResult );
  try { ev.SetClassifierResult ( currentPass, fBerkeleyAlphaBeta->GetName() + ":partialFitter", fBerkeleyAlphaBeta->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  timer.Stop();
  partialResult.SetExecutionTime( timer.RealTime() );

  return SetResult( partialResult, ev, currentPass);
}


Processor::Result
PartialFitterMultiPath::SetResult( DS::FitResult& partialResult,
                          DS::EV& ev,
                          size_t currentPass)
{

  ev.SetFitResult( currentPass, "partialFitter", partialResult );
  // Now set the 0 vertex as the default fit vertex, if the fit is valid
  DS::FitVertex defaultVertex = partialResult.GetVertex(0);
  // Pos, energy and time are calculated by the partial fitter so the DS::INVALID flag must be set
  defaultVertex.SetDirection(TVector3(DS::INVALID, DS::INVALID, DS::INVALID),false);
  defaultVertex.SetDirectionErrors(TVector3(DS::INVALID, DS::INVALID, DS::INVALID),false);
  ev.SetDefaultFitVertex( "partialFitter", defaultVertex );
  return Processor::OK;
}

void
PartialFitterMultiPath::EndOfRun( DS::Run& run )
{
  // First call the begin of run functions
  fQuadSeed->EndOfRun( run );
  fPositionTime->EndOfRun( run );
  // fPowell->EndOfRun( run );
  fNullSelector->EndOfRun( run );
  // fET1D->EndOfRun( run );
  fEnergySeed->EndOfRun( run );
  fPartialEnergy->EndOfRun( run );
  fITR->EndOfRun( run );
  fQPDT->EndOfRun( run );
  fTimingPeaks->EndOfRun( run );
  fMeanTime->EndOfRun( run );
  fPreTriggerHits->EndOfRun( run );
  fIsoRegions->EndOfRun( run );
  fEarlyTime->EndOfRun( run );
  fBiPoCumulTimeResid->EndOfRun( run );
  fBiPo212LikelihoodDiff->EndOfRun( run );
  fBiPo214LikelihoodDiff->EndOfRun( run );
  fBerkeleyAlphaBeta->EndOfRun( run );
}


void
PartialFitterMultiPath::SetFOMs( DS::FitResult& fitResult,
                      const DS::FitResult& partialResult,
                      const string& prefix )
{
  vector<string> fomNames = partialResult.GetFOMNames();
  for(vector<string>::const_iterator iter = fomNames.begin(); iter != fomNames.end(); ++iter)
    fitResult.SetFOM( (prefix+*iter), partialResult.GetFOM(*iter) );
}
