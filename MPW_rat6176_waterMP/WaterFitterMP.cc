#include <RAT/WaterFitterMP.hh>
#include <RAT/MethodFactory.hh>
#include <RAT/PDFFactory.hh>
#include <RAT/OptimiserFactory.hh>
#include <RAT/ClassifierFactory.hh>
#include <RAT/PMTSelectorFactory.hh>
#include <RAT/Method.hh>
#include <RAT/SmearResult.hh>
#include <RAT/Classifier.hh>
#include <RAT/OptimisedMethod.hh>
#include <RAT/PDFMethod.hh>
#include <RAT/SeededMethod.hh>
#include <RAT/SelectorMethod.hh>
#include <RAT/SelectorClassifier.hh>
#include <RAT/Optimiser.hh>
#include <RAT/SeededClassifier.hh>
#include <RAT/FitterPMT.hh>
#include <RAT/MetaInformation.hh>
#include <RAT/TimeResidualCut.hh>
#include <RAT/StraightTimeResidualCut.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <RAT/ListHelp.hh>
#include <RAT/MetaInformation.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DS/FitVertex.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/DataQCFlags.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/DataCleaningBits.hh>
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

#include <sstream>
#include <vector>
#include <string>
using std::vector;
using std::string;


WaterFitterMP::WaterFitterMP()
  : Processor("waterFitterMP")
{
  fQuadSeed = MethodFactory::Get()->GetMethod( "quad" ); // keep the quad for energy seed
  fQuadNhitCut = 3.0;
  /// set the multipath as default position time fitter
  fPositionTime = MethodFactory::Get()->GetMethod( "multipath-waterposition" );
  // As discussed in the reconstruction and analysis calls a nhit threshold, after selectors
  // is applied to decide whether the event should or not be fitted
  // By default the threshold is 0 (ie. all events), but for water phase it was agreed
  // to set a 6 nhit threshold. The reasoning is that under this number of hits
  // the fitter almost invariably fails, making the process needlessly slow
  fPositionTime->SetI("nhit_cut",4);//!!!nhit>4, lower down from 6 for MPW

  fPMTCalSelector = PMTSelectorFactory::Get()->GetPMTSelector( "PMTCalSelector" );
  fModeCut = PMTSelectorFactory::Get()->GetPMTSelector( "modeCut" );
  // optimized modeCut window for multipath fitter
  fModeCut->SetD( "lowLimit", -50.5 );
  fModeCut->SetD( "highLimit", 100.5 );
  fEnergySeed = MethodFactory::Get()->GetMethod( "energyPromptLookup" );
  fEnergyMethod = MethodFactory::Get()->GetMethod( "energyRSP" );
  fIsotropy = ClassifierFactory::Get()->GetClassifier( "isotropy" );
  fTimeResidualCut0 = PMTSelectorFactory::Get()->GetPMTSelector( "timeResidualCut" );
  fTimeResidualCut = PMTSelectorFactory::Get()->GetPMTSelector( "straightTimeResidualCut" ); // change from timeResidualCut
  fTimeResidualCut->SetD( "lowLimit", -10.5 );
  fTimeResidualCut->SetD( "highLimit", 120.5 );
  fDirectionModeCut = PMTSelectorFactory::Get()->GetPMTSelector( "modeCut" );// same to the position for MP
  fDirection = MethodFactory::Get()->GetMethod( "multipath-waterdirection" );
  fSmearResult = MethodFactory::Get()->GetMethod("smearResult");
  fPMTCalSelector = PMTSelectorFactory::Get()->GetPMTSelector( "PMTCalSelector" );
  fITR = ClassifierFactory::Get()->GetClassifier( "ITR" );
  fQPDT = ClassifierFactory::Get()->GetClassifier( "QPDT" );
  /// Now combine the components where appropriate
  dynamic_cast< SelectorMethod* >( fQuadSeed )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorMethod* >( fQuadSeed )->AddPMTSelector( fModeCut);
  dynamic_cast< SelectorMethod* >( fPositionTime )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorMethod* >( fPositionTime )->AddPMTSelector( fModeCut );
  dynamic_cast< SelectorClassifier* >( fIsotropy )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorClassifier* >( fIsotropy )->AddPMTSelector( fTimeResidualCut0 );
  dynamic_cast< SelectorClassifier* >( fITR )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorClassifier* >( fQPDT )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorMethod* >( fDirection )->AddPMTSelector( fPMTCalSelector );
  dynamic_cast< SelectorMethod* >( fDirection )->AddPMTSelector( fModeCut );
//  dynamic_cast< SelectorMethod* >( fDirection )->AddPMTSelector( fTimeResidualCut );

  fMuonWater = MethodFactory::Get()->GetMethod( "MuonWater" );

  // Get drive correction parameters from DB
  DBLinkPtr fDriveCorrection = DB::Get()->GetLink("FIT_MULTIPATH","snoplus_water");
  fDrivePositionScale = fDriveCorrection->GetD("driveCor_p0");
  fDriveDirectionScale = fDriveCorrection->GetD("driveCor_p1");
  fDriveEnable = fDriveCorrection->GetZ("driveCor_enable");
}

void WaterFitterMP::SetI(const std::string& param, const int value)
{
  // Use this to pass processor calls to other processors using . syntax
  vector<string> parts = split( param, "."  );
  // For now only pass paramets and use none in WaterFitterMP
  if( parts.size() != 2 )
    throw Processor::ParamUnknown( param );
  if( parts[0] == string( "smear" ) )
    fSmearResult->SetI(parts[1], value);
  else if( parts[0] == string( "posTime" ) )
    fPositionTime->SetI(parts[1], value);
  else
    throw Processor::ParamUnknown( param );
}

void WaterFitterMP::SetD(const std::string& param, const double value)
{
  // Use this to pass processor calls to other processors using . syntax
  vector<string> parts = split( param, "."  );
  // For now only pass paramets and use none in WaterFitterMP
  if( parts.size() != 2 )
    throw Processor::ParamUnknown( param );
  if( parts[0] == string( "smear" ) )
    fSmearResult->SetD(parts[1], value);
  else
    throw Processor::ParamUnknown( param );
}

void WaterFitterMP::SetS(const std::string& param, const std::string& value)
{
  // Use this to pass processor calls to other processors using . syntax
  vector<string> parts = split( param, "."  );
  // For now only pass paramets and use none in WaterFitterMP
  if( parts.size() != 2 )
    throw Processor::ParamUnknown( param );
  if( parts[0] == string( "smear" ) )
    fSmearResult->SetS(parts[1], value);
  else
    throw Processor::ParamUnknown( param );
}

WaterFitterMP::~WaterFitterMP()
{
  delete fIsotropy;
  delete fTimeResidualCut0;
  delete fTimeResidualCut;
  delete fQuadSeed;
  delete fPositionTime;
  delete fModeCut;
  delete fEnergySeed;
  delete fEnergyMethod;
  delete fDirection;
  delete fDirectionModeCut;
  delete fITR;
  delete fQPDT;
  delete fMuonWater;
  delete fSmearResult;
  delete fPMTCalSelector;
}

void
  WaterFitterMP::BeginOfRun( DS::Run& run )
{
  /// First call the begin of run functions
  fTimeResidualCut0->BeginOfRun( run );
  fTimeResidualCut->BeginOfRun( run );
  fQuadSeed->BeginOfRun( run );
  fPositionTime->BeginOfRun( run );
  fModeCut->BeginOfRun( run );
  fEnergySeed->BeginOfRun( run );
  fEnergyMethod->BeginOfRun( run );
  fDirectionModeCut->BeginOfRun( run );
  fDirection->BeginOfRun( run );
  fITR->BeginOfRun( run );
  fQPDT->BeginOfRun( run );
  fMuonWater->BeginOfRun( run );
  fSmearResult->BeginOfRun( run );
  fPMTCalSelector->BeginOfRun( run );
}

Processor::Result
WaterFitterMP::DSEvent( DS::Run& run,
                      DS::Entry& ds )
{
  for( size_t iEV = 0; iEV < ds.GetEVCount(); iEV++ )
    Event( run, ds.GetEV( iEV ) );
  return OK;
}

void
  WaterFitterMP::EndOfRun( DS::Run& run )
{
  // Now call the end of run functions
  fTimeResidualCut->EndOfRun( run );
  fTimeResidualCut0->EndOfRun( run );
  fQuadSeed->EndOfRun( run );
  fPositionTime->EndOfRun( run );
  fModeCut->EndOfRun( run );
  fEnergySeed->EndOfRun( run );
  fEnergyMethod->EndOfRun( run );
  fDirectionModeCut->EndOfRun( run );
  fDirection->EndOfRun( run );
  fITR->EndOfRun( run );
  fQPDT->EndOfRun( run );
  fMuonWater->EndOfRun( run );
  fSmearResult->EndOfRun( run );
  fPMTCalSelector->EndOfRun( run );
}

Processor::Result
WaterFitterMP::Event( DS::Run& run,
                    DS::EV& ev )
{
  /// WaterFitter Logic:
  ///  This is a simple overview of the logic that is expressed in the next 90 lines
  ///
  ///  0. See if event was caused by a muon (muon datacleaning tag
  ///     check), if so, run the muon fitter routines, else continue. \n
  ///  1. Fit the position & time as the seedResult from quad,
  ///      if nhits < fNhitCutoff then abort the fit and return no FitResult. \n
  ///  2. Fit the direction as directionSeedResult using simpleDirection with the quad as seed.
  ///     If this fails then abort the fit and return no FitResult. \n
  ///  3. Fit the position & time as waterResult using the positionTimeLikelihood, metaDriveCorrectSeed-powell,
  ///      modeCut, gv1d-lightwater-sno with the seedResult and directionSeedResult as the seeds
  ///     If this fails then abort the fit and return no FitResult. \n
  ///     If smearResult position is set to "on", then apply the given position smear
  ///  4. Fit the direction as directionResult using positionTimeDirectionLikelihood, simulatedAnnealing
  ///      modeCut, positionTimeDirectionPDF with waterResult and dirSeed as seeds
  ///     If this fails then abort the fit and return no FitResult. \n
  ///     If smearResult direction is set to "on", then apply the given direction smear
  ///  5. Fit the energy seed as energySeedResult using energyPromptLookup and the waterResult as seed
  ///     If this fails then abort the fit and return no FitResult. \n
  ///  6. Apply drive correction using fitted position and direction. See docdb #4592 for details. \n
  ///  7. Fit the energy as waterResult using energyRSP and energySeedResult as seed
  ///     If this fails then abort the fit and return no FitResult. \n
  ///     If smearResult energy is set to "on", then apply the given energy smear function
  ///  8. Run Isotropy with waterResult as seed and TimeResidualCut as PMT selector. \n
  ///  9. Run itr with waterResult as seed. \n
  ///  10. Run qpdt with waterResult as seed. \n

  /// NFB: FIXME : Need to reorganize the fitter to catch fit failures
  ///              that look like good fits.
  TStopwatch timer;
  timer.Start( true );

  const size_t currentPass = MetaInformation::Get()->GetCurrentPass();
  ev.AddFitterPass( currentPass ); // Ensure the EV knows fitters were run for this pass
  ev.AddClassifierPass( currentPass ); // Ensure the EV knows classifiers were run for this pass

  vector<FitterPMT> pmtData;
  for( size_t iPMTCal =0; iPMTCal < ev.GetCalPMTs().GetCount(); iPMTCal++ )
    pmtData.push_back( FitterPMT( ev.GetCalPMTs().GetPMT( iPMTCal ) ) );
  /// Additionally select pmtData so it *only* contains PMTs which pass the PMTCal groups recommended selector cuts.
  DS::FitVertex dummyVertex;
  pmtData = fPMTCalSelector->GetSelectedPMTs( pmtData, dummyVertex );

  // if event was tagged as muon, run muon fitting routines
  const DS::DataQCFlags& dcFlags = ev.GetDataCleaningFlags();
  const DU::DataCleaningBits& dcBits = DU::Utility::Get()->GetDataCleaningBits();
  // If the data cleaning tags have processed and the muon tag is set for this pass then run the muon fitter

  bool is_muon = false;
  Int_t dcPassNumber = dcFlags.GetLatestPass();
  if (dcFlags.ExistFlags( dcPassNumber )) {
    is_muon = (dcFlags.GetApplied(dcPassNumber).Get(dcBits.GetBitIndex( "muontag" )) && !dcFlags.GetFlags(dcPassNumber).Get(dcBits.GetBitIndex( "muontag" )));
  }
  if( is_muon )
    {
      FitResult muonResult;
      fMuonWater->SetEventData( pmtData, &ev, &run );
      try
        {
          muonResult = fMuonWater->GetBestFit();
          timer.Stop();
          muonResult.SetExecutionTime( timer.RealTime() );

          return SetResult( muonResult, ev, currentPass );
        }
      catch( Method::MethodFailError& error )
        {
          debug << "WaterFitterMP::Event: Muon fitting failed: "<< error.what() << ", quiting." << newline;
          // -- we still have to return a result (including a vertex), or the output
          //    processors will fail
          return SetResult(muonResult,ev,currentPass);
        }
    }

  // Initialize the seeds
  fQuadSeed->SetEventData( pmtData, &ev, &run );
  // Run the seed method : Quad
  FitResult seedResult;
  try
    {
      if( ev.GetNhitsCleaned() > fQuadNhitCut )
        seedResult = fQuadSeed->GetBestFit();
      else
        {
          detail << "WaterFitterMP::Event: insufficient points for quad to run, exiting" << newline;
          // Again, this will blow on our faces if there are not enough hits to fit
          // as there is no vertex for the output processor to grab
          return SetResult( seedResult, ev, currentPass );
        }
    }
  catch( Method::MethodFailError& error ) {
    debug << "WaterFitterMP::Event: Seed failed " << error.what() << ", exiting." << newline;
    // Quad only fails if no points were found inside the PSUP
    // But in that case we still have to prevent problems of the vertex not having been assigned
    return SetResult( seedResult, ev, currentPass );
  }
  /// Now initialise the position time method
  fPositionTime->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededMethod* >( fPositionTime )->DefaultSeed();
  /// Run the position time method
  FitResult waterResult;
  // Set a dummy (invalid) vertex so from here on we no longer have to worry about
  // unexisting vertices
  SetResult( waterResult, ev, currentPass );
  try
    {
      FitResult positionResult = fPositionTime->GetBestFit();
      // If position smear result is set then that will occur here
      // otherwise the positionResult is unchanged. The various types
      // of smearing are applied in order. Isotropic->xyz->rthetaphi.
      if( dynamic_cast< SmearResult* >( fSmearResult )->IsAppliedPosition() )
      {
        dynamic_cast< SeededMethod* >( fSmearResult )->SetSeed( positionResult );
        dynamic_cast< SmearResult* >( fSmearResult )->SetActivePosition();
        positionResult = fSmearResult->GetBestFit();
      }

      waterResult.SetVertex(0, positionResult.GetVertex(0));
      SetFOMs(waterResult, positionResult, "Position");
    }
  catch( Method::MethodFailError& error )
    {
      debug << "WaterFitterMP::Event: Main method failed " << error.what() << ", exiting" << newline;
      return Processor::FAIL;
    }

  /// Now initialise the direction method
  fDirection->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededMethod* >( fDirection )->DefaultSeed(); // no default seed for MPW direction
  dynamic_cast< SeededMethod* >( fDirection )->SetSeed( waterResult );
  /// Run the direction method
  try
  {
    FitVertex vertex = waterResult.GetVertex(0);
    FitResult directionResult = fDirection->GetBestFit();
    // If direction smear result is set then that will occur here
    // otherwise the directionResult is unchanged
    if( dynamic_cast< SmearResult* >( fSmearResult )->IsAppliedDirection() )
    {
      dynamic_cast< SeededMethod* >( fSmearResult )->SetSeed( directionResult );
      dynamic_cast< SmearResult* >( fSmearResult )->SetActiveDirection();
      directionResult = fSmearResult->GetBestFit();
    }
    FitVertex directionVertex = directionResult.GetVertex(0);//// problem here!!!!!
    vertex.SetDirection( directionVertex.GetDirection(), directionVertex.ValidDirection() );
    // vertex.SetPositiveDirectionError( directionVertex.GetPositiveDirectionError(), directionVertex.ValidDirection() );
    // vertex.SetNegativeDirectionError( directionVertex.GetNegativeDirectionError(), directionVertex.ValidDirection() );
    waterResult.SetVertex( 0, vertex );
    SetFOMs(waterResult, directionResult, "Direction");
    //cout<<"!!! direction OK"<<endl;
  }
  catch( Method::MethodFailError& error )
  {
    cout << "WaterFitterMP::Event: Direction method failed " << error.what() << ", exiting" << newline;
    //debug << "WaterFitterMP::Event: Direction method failed " << error.what() << ", exiting" << newline;
    return Processor::FAIL;
  }
  catch( Optimiser::OptimiserFailError& error )
  {
    cout << "WaterFitterMP::Event: Direction method failed " << error.what() << ", exiting" << newline;
    //debug << "WaterFitterMP::Event: Direction method failed " << error.what() << ", exiting" << newline;
    return Processor::FAIL;
  }
  // Apply dirve correction to the fitted position
  if(fDriveEnable == 1)
    DriveCorrection(waterResult);

  SetResult( waterResult, ev, currentPass );
  // NFB: Check if the current waterResult vertex is valid.
  // if not, just return. No point in doing a energy fit without a position and time
  if (!waterResult.GetValid()) {
    detail << "WaterFitterMP::Event: Invalid vertex. Bypassing energy fit." << newline;
    return Processor::FAIL;
  }

  // Now initialise and run the ITR classifier
  fITR->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fITR )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fITR )->SetSeed( waterResult );
  try { ev.SetClassifierResult( currentPass, fITR->GetName() + ":waterFitterMP", fITR->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the Isotropy classifier
  fIsotropy->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fIsotropy )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fIsotropy )->SetSeed( waterResult );
  try { ev.SetClassifierResult( currentPass, fIsotropy->GetName() + ":waterFitterMP", fIsotropy->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise and run the QPDT classifier
  fQPDT->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededClassifier* >( fQPDT )->DefaultSeed();
  dynamic_cast< SeededClassifier* >( fQPDT )->SetSeed( waterResult );
  try { ev.SetClassifierResult( currentPass, fQPDT->GetName() + ":waterFitterMP", fQPDT->GetClassification() ); }
  catch( Classifier::ClassifierFailError& error ) { warn << error.what() << newline; }

  // Now initialise the energy seed
  //cout<<"are you OK? energy"<<endl; 
  fEnergySeed->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededMethod* >( fEnergySeed )->DefaultSeed();
  if( waterResult.GetVertexCount() && waterResult.GetVertex(0).ValidPosition() )
    dynamic_cast< SeededMethod* >( fEnergySeed )->SetSeed( waterResult );
  else if( seedResult.GetVertexCount() && seedResult.GetVertex(0).ValidPosition() )
    dynamic_cast< SeededMethod* >( fEnergySeed )->SetSeed( seedResult );
  else
    warn << "WaterFitterMP::Event: No seed for the energyPromptLookup method." << newline;
  // Run the energy seed
  FitResult energySeedResult;
  try
    {
      energySeedResult = fEnergySeed->GetBestFit();
    }
  catch( Method::MethodFailError& error )
    {
      warn << "WaterFitterMP::Event: Energy seed method failed " << error.what() << ", exiting" << newline;
      return Processor::FAIL;
    }

  // Now initialise the energy method
  fEnergyMethod->SetEventData( pmtData, &ev, &run );
  dynamic_cast< SeededMethod* >( fEnergyMethod )->DefaultSeed();
  dynamic_cast< SeededMethod* >( fEnergyMethod )->SetSeed( energySeedResult );

  // Run the energyRSP method
  try
    {
      FitVertex vertex = waterResult.GetVertex(0);
      FitResult energyResult = fEnergyMethod->GetBestFit();
      FitVertex energyVertex;
      if( dynamic_cast< SmearResult* >( fSmearResult )->IsAppliedEnergy() )
      {
        dynamic_cast< SeededMethod* >( fSmearResult )->SetSeed( energySeedResult );
        dynamic_cast< SmearResult* >( fSmearResult )->SetActiveEnergy();
        energyResult = fSmearResult->GetBestFit();
      }
      energyVertex = energyResult.GetVertex(0);
      vertex.SetEnergy( energyVertex.GetEnergy(), energyVertex.ValidEnergy() );
      vertex.SetPositiveEnergyError( energyVertex.GetPositiveEnergyError(), energyVertex.ValidEnergy() );
      vertex.SetNegativeEnergyError( energyVertex.GetNegativeEnergyError(), energyVertex.ValidEnergy() );
      waterResult.SetVertex( 0, vertex );
      SetFOMs(waterResult, energyResult, "Energy");
      waterResult.SetFOM("EnergyRSP", 1.0);
    }
  catch( Method::MethodFailError& error )
  {
      warn << "WaterFitter::Event: Energy method failed " << error.what() << ", exiting" << newline;
      return Processor::FAIL;
  }

  timer.Stop();
  waterResult.SetExecutionTime( timer.RealTime() );
  return SetResult(waterResult, ev, currentPass);
}

Processor::Result
WaterFitterMP::SetResult( DS::FitResult& waterResult,
                        DS::EV& ev,
                        size_t currentPass)
{
  // NFB: If the resultset does not contain a vertex the job fails
  ev.SetFitResult( currentPass, "waterFitterMP", waterResult );
  try {
    // Now set the 0 vertex as the default fit vertex, if the fit is valid
    ev.SetDefaultFitVertex( "waterFitterMP", waterResult.GetVertex(0) );
  }
  catch (FitResult::NoVertexError&) {
     // Add a dummy/invalid vertex so that other processors down the line still work with it
     waterResult.SetVertex(0,FitVertex());
     ev.SetFitResult( currentPass, "waterFitterMP", waterResult );
     ev.SetDefaultFitVertex( "waterFitterMP", waterResult.GetVertex(0) );
     return Processor::FAIL;
   }
  return Processor::OK;

}

void
WaterFitterMP::SetFOMs( DS::FitResult& fitResult,
                      const DS::FitResult& partialResult,
                      const string& prefix )
{
  vector<string> fomNames = partialResult.GetFOMNames();
  for(vector<string>::const_iterator iter = fomNames.begin(); iter != fomNames.end(); ++iter)
    fitResult.SetFOM( (prefix+*iter), partialResult.GetFOM(*iter) );
}

void
WaterFitterMP::DriveCorrection( DS::FitResult& fitResult )
{
  // Apply drive correction based on docdb #4592
  if( fitResult.GetVertexCount() && fitResult.GetVertex(0).ValidPosition() && fitResult.GetVertex(0).ValidDirection() ){
    // Get current position and direction
    TVector3 fitPos = fitResult.GetVertex(0).GetPosition();
    TVector3 fitDir = fitResult.GetVertex(0).GetDirection();
    // scale vectors by appropriate parameters
    fitPos.SetMag(fitPos.Mag()*fDrivePositionScale);
    fitDir.SetMag(fitDir.Mag()*fDriveDirectionScale);
    // Set the corrected position
    fitResult.GetVertex(0).SetPosition(fitPos + fitDir);
  } else {
    debug << "WaterFitterMP::DriveCorrection: Couldn't find valid position and or direction in fit" << newline;
  }
}
