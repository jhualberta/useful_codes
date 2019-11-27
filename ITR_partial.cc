#include <RAT/ITR.hh>
#include <RAT/DB.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DS/FitVertex.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/EffectiveVelocity.hh>
using namespace RAT;
using namespace RAT::Classifiers;
using namespace RAT::DS;

#include <math.h>
using namespace std;

void
ITR::Initialise( const std::string& param)
{
  fIndex = param;
}

void
ITR::BeginOfRun( DS::Run& )
{
  DB *db = DB::Get();
  string index;
  if( fIndex.empty() )
  {
    try{
      index = db->GetLink("GEO","inner_av")->GetS("material");
    }
    catch(DBNotFoundError & ){
      try{
        // If we are in a partial-fill geometry, use the upper
        // material
        index = db->GetLink("GEO","inner_av")->GetS("material_top");
      }
      catch(DBNotFoundError & ) {
        Log::Die("ITR: inner_av material and material_top "
        "undefined. Try something like:\n"
        "/rat/db/set GEO[inner_av] material "
        "\"te_labppo_scintillator\" or "
        "/rat/db/set GEO[inner_av] material_top "
        "\"te_labppo_scintillator\"");
      }
    }
  }
  else
  {
    index = fIndex;
  }
  string defaultIndex = "labppo_scintillator";
  DBLinkPtr dbITRLink = db->GetLink("CLASSIFIER_ITR",defaultIndex);
  DBLinkGroup grp = db->GetLinkGroup("CLASSIFIER_ITR");
  DBLinkGroup::iterator it;

  for(it = grp.begin(); it != grp.end(); ++it)
  {
    if(index == it->first)
    {
      dbITRLink = db->GetLink("CLASSIFIER_ITR",index);
    }
  }

  fLowerTimeLimit = dbITRLink->GetD("t1");
  fUpperTimeLimit = dbITRLink->GetD("t2");

  fLightPath = DU::Utility::Get()->GetLightPathCalculator();

}

void
ITR::DefaultSeed()
{
  fEventPos = TVector3();
  fEventTime = 0.0;
  fSeedVertex.SetPosition( fEventPos, false, true );
  fSeedVertex.SetTime( fEventTime, false, true );
}

void
ITR::SetSeed( const DS::FitResult& seed )
{
  try
    {
      fSeedVertex = seed.GetVertex(0);
    }
  catch( FitResult::NoVertexError& error )
    {
      // Nothing to seed from, strange request
      return;
    }
  try { fEventPos = fSeedVertex.GetPosition(); }
  catch( FitVertex::NoValueError& error ) { /* No data to seed from */ }
  try { fEventTime = fSeedVertex.GetTime(); }
  catch( FitVertex::NoValueError& error ) { /* No data to seed from */ }
}

DS::ClassifierResult
ITR::GetClassification()
{
  fClassifierResult.Reset();
  if( fPMTData.empty() )
    return fClassifierResult;

  SelectPMTData( fSeedVertex );

  const DU::PMTInfo& pmtInfo = DU::Utility::Get()->GetPMTInfo();
  int numInTime = 0;
  for( vector<FitterPMT>::const_iterator iPMT = fSelectedPMTData.begin(); iPMT != fSelectedPMTData.end(); ++iPMT )
    {
      const TVector3 pmtPos = pmtInfo.GetPosition( iPMT->GetID() );
      fLightPath.CalcByPosition(fEventPos, pmtPos);
      double distInAV = fLightPath.GetDistInAV();
      double distInWater = fLightPath.GetDistInWater();
      double distInUpperTarget = fLightPath.GetDistInUpperTarget();
      double distInLowerTarget = fLightPath.GetDistInLowerTarget();
      const double transitTime = RAT::DU::Utility::Get()->GetEffectiveVelocity().CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );
      const double corTime = iPMT->GetTime() - transitTime - fEventTime;

      if( corTime < fUpperTimeLimit && corTime > fLowerTimeLimit )
        {
          numInTime++;
        }
    }

  fClassifierResult.SetClassification( "ITR", static_cast<double>( numInTime ) / static_cast<double>( fSelectedPMTData.size() ) );
  fClassifierResult.SetValid( true );
  return fClassifierResult;
}
