#include <RAT/MultiPathProcessor.hh>
#include <RAT/PMTSelectorFactory.hh>
#include <RAT/FitterPMT.hh>
#include <RAT/MetaInformation.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/MultiPathFitter.hh>
#include <RAT/FitterWLS.hh>
#include <RAT/FitterAW.hh>
#include <RAT/FitterWaterPosition.hh>
#include <RAT/FitterWaterDirection.hh>
#include <string>

using namespace RAT;
using namespace RAT::PMTSelectors;
using namespace RAT::DS;

#include <sstream>
using namespace std;

#include <TStopwatch.h>
using namespace ROOT;

int MultiPathProcessor::sNumberOfEvents = 0;


MultiPathProcessor::MultiPathProcessor()
  : Processor("multiPathProcessor")
{

}

MultiPathProcessor::~MultiPathProcessor()
{
  return;
}

void
MultiPathProcessor::BeginOfRun( DS::Run& run )
{
    // First call the begin of run functions
  if (sValue=="WLSFitter"){FitterWLS::BeginOfRun(run);}
  if (sValue=="AWFitter"){FitterAW::BeginOfRun(run);}
  if (sValue=="WATERFitter"){FitterWaterPosition::BeginOfRun(run);FitterWaterDirection::BeginOfRun(run);}

}


void MultiPathProcessor::EndOfRun( DS::Run& run  )
{
  return ;
}


void
MultiPathProcessor::SetI( const std::string & param, const int value)
{
  if (param == "find_likelihood"){
    sNumberOfEvents = value;
  }
  else {
    throw ParamUnknown( param ) ;
  }
}


void
MultiPathProcessor::SetS( const std::string& param, const std::string& value )
{
  if ( param == "sFitter" ){
    if( value == "WLSFitter") {
      sValue = "WLSFitter" ;
      nFitParameter = 6 ;
    }
    if( value == "AWFitter") {
      sValue = "AWFitter" ;
      nFitParameter = 4 ;
    }
    if( value == "WATERFitter") {
      sValue = "WATERFitter" ;
      nFitParameter = 4 ;
    }
   }
   else {
    throw ParamUnknown( param ) ;
  }
}

Processor::Result
MultiPathProcessor::DSEvent( DS::Run& run,
                             DS::Entry& ds )
{
  size_t mcevent = run.GetNMCEvents();
  
  int count=1;
  for( size_t iEV = 0; iEV < ds.GetEVCount(); iEV++ ){
    size_t nPMT=ds.GetEV( iEV ).GetCalPMTs().GetCount();
    if(nPMT==0) continue;
    Event( run, ds.GetEV( iEV ), mcevent  );
   
  }
  return OK;
}


Processor::Result
MultiPathProcessor::Event( DS::Run& run,
                           DS::EV& ev, size_t mcevent )
{
  TStopwatch timer;
  timer.Start( true );

  const size_t currentPass = MetaInformation::Get()->GetCurrentPass();
  // Ensure the EV knows fitters were run for this pass
  ev.AddFitterPass( currentPass );
  // Ensure the EV knows classifiers were run for this pass
  ev.AddClassifierPass( currentPass );
  MultiPathFitter::Clear();
  size_t nPMT=ev.GetCalPMTs().GetCount();
  for( size_t iPMTCal = 0; iPMTCal < nPMT; iPMTCal++ )
    {
     //if(ev.GetGTID()==1286306 && iPMTCal == nPMT-1) {
	//	std::cout<<"number of pmts "<< nPMT <<" and we are at "<< iPMTCal <<std::endl ; (time<-15000)
     // }
      PMTCal &cal=ev.GetCalPMTs().GetPMT(iPMTCal);
      double time= cal.GetTime();
      //std::cout<< " pmt ID "<< cal.GetID()<<" pmt Hit Time "<< time << std::endl ;
      //if(time > 192.5 && time < 92.5) continue ;
      if(time<-9000) continue ;
      const TVector3 &position=DU::Utility::Get()->GetPMTInfo().GetPosition(cal.GetID() );
      
      if(sValue == "WLSFitter") new FitterWLS(time, position ,cal.GetID(), nFitParameter );
      if(sValue == "AWFitter") new FitterAW(time, position ,cal.GetID(), nFitParameter );
      if(sValue == "WATERFitter") new FitterWaterPosition(time, position ,cal.GetID(), nFitParameter );
    }

  MultiPathFitter::Fit(mcevent, ev.GetGTID(),1);

  RAT::DS::FitResult multiPathResult;
  RAT::DS::FitVertex vertex;

  //FitVertex vertex;
  TVector3 fittedDirection ;
  //FitResult multiPathResult;
  if(sValue == "WLSFitter") {
    vertex.SetPosition(FitterWLS::fVertex);
    fittedDirection.SetMagThetaPhi(1.0, FitterWLS::fZenith, FitterWLS::fAzimuth);
    vertex.SetDirection(fittedDirection);
    vertex.SetTime(FitterWLS::fTime0);
    multiPathResult.SetFOM("FitterWLS",FitterWLS::fFOM);
  }
  if(sValue == "AWFitter") {
    vertex.SetPosition(FitterAW::fVertex);
    vertex.SetDirection(FitterAW::fBestFitDirection);
    vertex.SetTime(FitterAW::fTime0);
    multiPathResult.SetFOM("FitterAW",FitterAW::fFOM);
  }
  if(sValue == "WATERFitter") 
   {
    TVector3 FittedVertex ;double fittedTime ;
    if(nPMT==0) { FittedVertex.SetXYZ(-10000.0,0.0,0.0) ; fittedTime=-9000 ; }
    else{
    FittedVertex= FitterWaterPosition::fVertex ;
    fittedTime = FitterWaterPosition::fTime0 ;
    }
    MultiPathFitter::Clear();
    nFitParameter=2 ;
    //size_t nPMT=ev.GetCalPMTs().GetCount();
  
    for( size_t iPMTCal = 0; iPMTCal < nPMT; iPMTCal++ )
    {
      PMTCal &cal=ev.GetCalPMTs().GetPMT(iPMTCal);
      double time= cal.GetTime();
      if(time > 192.5 && time < 92.5 )continue;
      const TVector3 &position=DU::Utility::Get()->GetPMTInfo().GetPosition(cal.GetID() );
      new FitterWaterDirection(time, position ,cal.GetID(), nFitParameter );
      if(iPMTCal==0) FitterWaterDirection::SetVertex(FittedVertex);
    }
    MultiPathFitter::Fit(mcevent, ev.GetGTID(),0);
    
    vertex.SetPosition(FittedVertex);
    vertex.SetTime(fittedTime);
    fittedDirection.SetMagThetaPhi(1.0,FitterWaterDirection::fZenith, FitterWaterDirection::fAzimuth);
    vertex.SetDirection(fittedDirection);
    nFitParameter=4 ;    
    multiPathResult.SetFOM("FitterWater",FitterWaterPosition::fFOM);
  }
  multiPathResult.SetVertex(0,vertex);
  if(sValue != "WATERFitter")MultiPathFitter::DumpLikelihood(sNumberOfEvents);
  
  return SetResult( multiPathResult, ev, currentPass);
}


Processor::Result
MultiPathProcessor::SetResult( DS::FitResult& multiPathResult,
                               DS::EV& ev,
                               size_t currentPass)
{
  ev.SetFitResult( currentPass, "MultiPathProcessor", multiPathResult );
  // Now set the 0 vertex as the default fit vertex, if the fit is valid
  ev.SetDefaultFitVertex( "MultiPathProcessor", multiPathResult.GetVertex(0) );

  return Processor::OK;
}
