////////////////////////////////////////////////////////////////////
/// \file PlotHitTimeResiduals.cc
///
/// \brief Functions to plot hit time residuals.
///
/// \author P G Jones <p.g.jones@qmul.ac.uk>
///
/// REVISION HISTORY:\n
///     2014-03-27 : P G Jones - First Revision.\n
///
/// \details EV Calibrated hit times are plotted minus transit times
/// based on the MC position or the fitted position.
///
////////////////////////////////////////////////////////////////////

#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <string>
using namespace std;
//TH2D* h_tres_cos = new TH2D("h_tres_cos","cos theta vs t_res for water",300,-500,800,100,-1.2,1.2);
TH2D* h_tres_cos = new TH2D("h_tres_cos","cos theta vs t_res for scint",900,-100,800,200,-1.2,1.2);
  
/// Plot the hit time residuals for the MC position
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the histogram plot
TH1D* PlotHitTimeResidualsMCPosition( const string& fileName )
{ 
  double anglez; 
  TH1D* hHitTimeResiduals = new TH1D( "hHitTimeResidualsMC events in scint (z>0)", "Hit time residuals using the MC position events in scint (z>0)", 400, -100.0, 300.0 );

  // If this is being done on data that does not require remote database connection
  // eg.: a simple simulation with default run number (0)
  // We can disable the remote connections:
  //
  // NOTE: Don't do this if you are using real data!!!
  //RAT::DB::Get()->SetAirplaneModeStatus(true);

  RAT::DU::DSReader dsReader( fileName );

  // RAT::DU::Utility::Get()->GetLightPathCalculator() must be called *after* the RAT::DU::DSReader constructor.
  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo(); // The PMT positions etc..
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      const TVector3 eventPosition = rDS.GetMC().GetMCParticle(0).GetPosition(); // At least 1 is somewhat guaranteed
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
        {
          const RAT::DS::EV& rEV = rDS.GetEV( iEV );
          const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
          for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ )
            {
              const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );
	      //lightPath.CalcByPositionPartial( eventPosition, pmtInfo.GetPosition( pmtCal.GetID() ) );		
              //double distInUpperTarget = lightPath.GetDistInUpperTarget();
              //double distInLowerTarget = lightPath.GetDistInLowerTarget();
              //double distInAV = lightPath.GetDistInAV();
              //double distInWater = lightPath.GetDistInWater();
              //const double transitTime =groupVelocity.CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );

              lightPath.CalcByPosition( eventPosition, pmtInfo.GetPosition( pmtCal.GetID() ) );
	      double distInInnerAV = lightPath.GetDistInInnerAV();
              double distInAV = lightPath.GetDistInAV();
              double distInWater = lightPath.GetDistInWater();
              const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon

              //if(eventPosition.Z()>0)
	      {
                  hHitTimeResiduals->Fill( pmtCal.GetTime() - transitTime - 390 + rDS.GetMCEV(iEV).GetGTTime());
                  TVector3 evePos = rDS.GetMC().GetMCParticle(0).GetPosition();
                  TVector3 pmtpos=pmtInfo.GetPosition( pmtCal.GetID() );
                  TVector3 evemom = rDS.GetMC().GetMCParticle(0).GetMomentum();
                  anglez = evemom*(evePos-pmtpos).Unit();
                  h_tres_cos->Fill(pmtCal.GetTime() - transitTime - 390 + rDS.GetMCEV(iEV).GetGTTime(),anglez);
               }
             }
          }
    }

   TCanvas *C1 = new TCanvas("C1");
   C1->cd();
   gPad->SetLogy();
   hHitTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
   hHitTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
   hHitTimeResiduals->Draw();
   TCanvas *C2 = new TCanvas("C2");
   C2->cd();
   h_tres_cos->GetYaxis()->SetTitle( "Cos theta (angle betwwen hit pmt and momentum)" );
   h_tres_cos->GetXaxis()->SetTitle( "Hit time residuals" );
   h_tres_cos->Draw("COLZ");
   return hHitTimeResiduals;
}

/// Plot the hit time residuals for the fit position
///
/// @param[in] fileName of the RAT::DS root file to analyse
/// @return the histogram plot
TH1D* PlotHitTimeResidualsFitPosition( const string& fileName, std::string fitName = "partialFitter")
{
  TH1D* hHitTimeResiduals = new TH1D( "hHitTimeResidualsFit", "Hit time residuals using the Fit position", 1000, -500.0, 500.0 );
  // If this is being done on data that does not require remote database connection
  // eg.: a simple simulation with default run number (0)
  // We can disable the remote connections:
  //
  // NOTE: Don't do this if you are using real data!!!
  RAT::DB::Get()->SetAirplaneModeStatus(true);

  RAT::DU::DSReader dsReader( fileName );

  // RAT::DU::Utility::Get()->GetLightPathCalculator() must be called *after* the RAT::DU::DSReader constructor.
  RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
  const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo(); // The PMT positions etc..
  for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
    {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      for( size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++ )
        {
          const RAT::DS::EV& rEV = rDS.GetEV( iEV );

          // grab the fit information
          if(fitName == "")
              fitName = rEV.GetDefaultFitName();

          TVector3 eventPosition;
          double   eventTime;

          try{
              const RAT::DS::FitVertex& rVertex = rEV.GetFitResult(fitName).GetVertex(0);
              if(!(rVertex.ValidPosition() && rVertex.ValidTime()))
                  continue; // fit invalid

              eventPosition = rVertex.GetPosition();
              eventTime = rVertex.GetTime();
          }
          catch(const RAT::DS::FitCollection::NoResultError&){
              // no fit result by the name of fitName
              continue;
          }
          catch (const RAT::DS::FitResult::NoVertexError&){
              // no fit vertex
              continue;
          }
          catch(const RAT::DS::FitVertex::NoValueError&){
              // position or time missing
              continue;
          }
          // DataNotFound --> implies no fit results are present, don't catch.

          // calculate time residuals
          const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
          for( size_t iPMT = 0; iPMT < calibratedPMTs.GetCount(); iPMT++ )
            {
              const RAT::DS::PMTCal& pmtCal = calibratedPMTs.GetPMT( iPMT );

              lightPath.CalcByPosition( eventPosition, pmtInfo.GetPosition( pmtCal.GetID() ) );
              double distInInnerAV = lightPath.GetDistInInnerAV();
              double distInAV = lightPath.GetDistInAV();
              double distInWater = lightPath.GetDistInWater();
              const double transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater ); // Assumes a 400nm photon
              // Time residuals estimate the photon emission time relative to the event start so subtract off the transit time and eventTime
              hHitTimeResiduals->Fill( pmtCal.GetTime() - transitTime - eventTime);
            }
        }
   }
  hHitTimeResiduals->GetYaxis()->SetTitle( "Count per 1 ns bin" );
  hHitTimeResiduals->GetXaxis()->SetTitle( "Hit time residuals [ns]" );
  hHitTimeResiduals->Draw();
  return hHitTimeResiduals;
}

/// Plot both the MC and Fitted position residuals
///
/// @param[in] fileName of the RAT::DS root file to analyse
void PlotHitTimeResiduals( const string& fileName )
{
  gStyle->SetFillColor( kWhite );
  TCanvas* c1 = new TCanvas();
  TH1D* mc = PlotHitTimeResidualsMCPosition( fileName );
  TH1D* fit = PlotHitTimeResidualsFitPosition( fileName );
  mc->Draw();
  fit->SetLineColor( kGreen + 2 );
  fit->Draw("SAME");
  TLegend* t1 = new TLegend( 0.7, 0.7, 0.9, 0.9 );
  t1->AddEntry( mc, "MC Position", "l" );
  t1->AddEntry( fit, "Fit Position", "l" );
  t1->Draw();
  c1->Update();
}
