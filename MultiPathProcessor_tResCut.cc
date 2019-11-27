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
#include <algorithm> //Jie add
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
 //fModeCut = PMTSelectorFactory::Get()->GetPMTSelector( "modeCut" );//Jie, this not compatible with Rat631 PMTselector!!
}

MultiPathProcessor::~MultiPathProcessor()
{
  //delete fModeCut;//Jie
  return;
}

void
MultiPathProcessor::BeginOfRun( DS::Run& run )
{

  //fModeCut->BeginOfRun( run );//Jie
    // First call the begin of run functions
  if (sValue=="WLSFitter"){FitterWLS::BeginOfRun(run);}
  if (sValue=="AWFitter"){FitterAW::BeginOfRun(run);}
  if (sValue=="WATERFitter"){FitterWaterPosition::BeginOfRun(run);FitterWaterDirection::BeginOfRun(run);}

}


void MultiPathProcessor::EndOfRun( DS::Run& run  )
{
  //fModeCut->EndOfRun( run );//Jie
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
 
 //Jie: test Mode cuts, loop calPMTs and calculate Mode Time before fitting
/*
  double modeCutVal = 100;
  double fLowCut  = -1*modeCutVal;
  double fHighCut =  modeCutVal;
  vector<Int_t> iPMTselect;
  vector<double> hitTimeData;
  TH1D *htempMode = new TH1D("Temp", "Temp", 500, 0.0, 500.0);
  for( size_t iPMTCal = 0; iPMTCal < nPMT; iPMTCal++ )
  {
      PMTCal &cal=ev.GetCalPMTs().GetPMT(iPMTCal);
      double time= cal.GetTime();
      //std::cout<< " pmt ID "<< cal.GetID()<<" pmt Hit Time "<< time << std::endl ;
      //if(time > 192.5 && time < 92.5) continue ;
      htempMode->Fill(time);
      hitTimeData.push_back(time);
  }
   
  double modeTime = htempMode->GetXaxis()->GetBinCenter( htempMode->GetMaximumBin() );
  double medianTime;
  if(hitTimeData.size()%2==0 )
   medianTime =   ( hitTimeData[hitTimeData.size()/2] + hitTimeData[hitTimeData.size()/2 + 1] ) / 2;
  else medianTime = hitTimeData[hitTimeData.size()/2];
 
  for( size_t iPMTCal = 0; iPMTCal < nPMT; iPMTCal++ )
  {
      PMTCal &cal=ev.GetCalPMTs().GetPMT(iPMTCal);
      double time= cal.GetTime();
      if(time > (modeTime+fLowCut) && time < (modeTime+fHighCut)) 
       {iPMTselect.push_back(iPMTCal);}
      else {cout<<"GTID "<<ev.GetGTID()<<" "<<iPMTCal<<" modeTime "<<modeTime<<" hitTime "<<time<<endl;}
  }
  
 delete htempMode;
*/
 //Jie

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
    double tResCutVal = 20;//Jie
    double grVelocity = 2.17554021555098529e+02;//Jie
    if(nPMT==0) { FittedVertex.SetXYZ(-10000.0,0.0,0.0) ; fittedTime=-9000 ; }
    else{//Jie: Get MPW fit position, calculate tRes
       FittedVertex= FitterWaterPosition::fVertex ;//original two lines, save fit
       fittedTime = FitterWaterPosition::fTime0 ;//original two lines, save fit
        //Jie: test time residual cuts, using position, time results from first fit 
       vector<Int_t> iPMTselect;//select PMTs to throw away  
       vector<double> tResData;
       for( size_t iPMTCal = 0; iPMTCal < nPMT; iPMTCal++ )
       {
         PMTCal &cal=ev.GetCalPMTs().GetPMT(iPMTCal);
         double time= cal.GetTime();
         const TVector3 &PMTpos = DU::Utility::Get()->GetPMTInfo().GetPosition(cal.GetID() );
         double tRes = time-fittedTime-(PMTpos-FittedVertex).Mag()/grVelocity;
         tResData.push_back(tRes);
         cout<<"first fit "<<ev.GetGTID()<<" tRes "<<tRes<<endl;
         if(tRes<tResCutVal && tRes>-50) iPMTselect.push_back(iPMTCal);//throw away PMTs
       }//Jie

      int oldCount = nPMT;
      int newCount = iPMTselect.size();
      vector<int> tempNew;
      vector<int> tempOld;
      for(std::vector<Int_t>::iterator it = iPMTselect.begin(); it < iPMTselect.end(); it++)
      {tempOld.push_back(*it);}

      for(int loop=0;loop<3;loop++)
      { 
		cout<<"new calPMTs "<<newCount<<" old calPMTs "<<oldCount<<endl;
		MultiPathFitter::Clear();
		if(newCount!=oldCount)
		{ 
			tempNew.clear();
      	  for(std::vector<Int_t>::iterator it1 = tempOld.begin(); it1 < tempOld.end(); it1++)
     	  {
		   PMTCal &cal=ev.GetCalPMTs().GetPMT(*it1);
           double time= cal.GetTime();
           const TVector3 &PMTpos = DU::Utility::Get()->GetPMTInfo().GetPosition(cal.GetID() );
           double tRes = time-fittedTime-(PMTpos-FittedVertex).Mag()/grVelocity;
           cout<<"fit again, loop "<<loop<<" "<<ev.GetGTID()<<" tRes "<<tRes<<endl;
           if(tRes<tResCutVal && tRes>-50) tempNew.push_back(*it1);//throw away PMTs
	   new FitterWaterPosition(time, PMTpos ,cal.GetID(), nFitParameter ); 
	  }
	 oldCount = tempOld.size();
	 newCount = tempNew.size();
   	 MultiPathFitter::Fit(mcevent, ev.GetGTID(),1); 
         tempOld.clear();

         for(std::vector<Int_t>::iterator it2 = tempNew.begin(); it2 < tempNew.end(); it2++)
         {tempOld.push_back(*it2);}

	}	
	else{break;}
      }//fitting loop
       MultiPathFitter::Clear();
       FittedVertex= FitterWaterPosition::fVertex ;//original two lines, save fit
       fittedTime = FitterWaterPosition::fTime0 ;//original two lines, save fit
       iPMTselect.clear();
   }//for position
  
    MultiPathFitter::Clear();
    nFitParameter=2 ;
    //size_t nPMT=ev.GetCalPMTs().GetCount();

//Fit Direction
  for( size_t iPMTCal = 0; iPMTCal < nPMT; iPMTCal++ )   
//    for(std::vector<Int_t>::iterator it2 = iPMTselect.begin(); it2 < iPMTselect.end(); it2++)
    {
      PMTCal &cal=ev.GetCalPMTs().GetPMT(iPMTCal);
//      PMTCal &cal=ev.GetCalPMTs().GetPMT(*it2);//Jie
      double time= cal.GetTime();
      if(time > 192.5 && time < 92.5 )continue;//Jie: currently turned off 
      const TVector3 &position=DU::Utility::Get()->GetPMTInfo().GetPosition(cal.GetID() );
      new FitterWaterDirection(time, position ,cal.GetID(), nFitParameter );
      if(iPMTCal==0) FitterWaterDirection::SetVertex(FittedVertex);
//      if(*it2==0) FitterWaterDirection::SetVertex(FittedVertex);//Jie
    }
    MultiPathFitter::Fit(mcevent, ev.GetGTID(),0);
    
    vertex.SetPosition(FittedVertex);
    vertex.SetTime(fittedTime);
    fittedDirection.SetMagThetaPhi(1.0,FitterWaterDirection::fZenith, FitterWaterDirection::fAzimuth);
    vertex.SetDirection(fittedDirection);
    
   // Do the drive correction here:
    double p0 = 1, p1 = 0;//currently turned off. For n=1.40, p0=0.997143, p1= -64.2801
    vertex.SetPosition(p0*FittedVertex+p1*fittedDirection);//Jie
    nFitParameter=4 ; 
    multiPathResult.SetFOM("FitterWater",FitterWaterPosition::fFOM);
  }
  multiPathResult.SetVertex(0,vertex);
  if(sValue == "WATERFitter" && sNumberOfEvents ==3125)//Jie: check one event with gtid 
   MultiPathFitter::DumpLikelihood(sNumberOfEvents);
  
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
