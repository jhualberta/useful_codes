#include <RAT/DS/Meta.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include <vector>
#include "TH2.h"
#include "TH1.h"
#include <TVector3.h>
#include "TF1.h"
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"

void lightPathCalculator()
{
   const char* filename = "FitRat_opticalphoton_shortWL_x1760_y1760_z5464_WL4400.root";

//   const char* filename = "FitRat_opticalphoton_shortWL_z1000.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   bool scint = false;
   TVector3 mcPos, pmtpos;
// z = 5464 
   int pmtid[]={1824};//9405};//,1824};
// z = 5450
 //  int pmtid[]= {734, 1832, 1321, 1760, 9159, 1390, 6074, 397, 6805, 8938, 1751, 4898, 8186, 7330, 781, 249, 8501, 1360, 1834, 9617, 2739, 1763};

// z = -1000
//   int pmtid[]={5946 ,2575 ,4490 ,27,1781 ,8865 ,712,8657 ,2731 ,1344 ,1894 ,5600 ,2166 ,9540 ,2854 ,7378 ,4733 ,1770 ,4787,7938 ,2012 ,582,594,460,4529,486,706,6904 ,1595 ,2298 ,3680 ,8021 ,4106};

// z = +1000
//   int pmtid[] = {3318, 5472, 7131, 9079, 3611, 9084, 9609, 5301, 5525, 7078, 4254, 6725, 3052, 845 , 5473, 1153, 9574, 2528, 4700, 5155, 7813, 7662, 5798, 8710, 3511, 2819, 4707, 2161, 4185, 5096, 643 , 3711, 3115, 6829, 7969, 8090, 5802};
   size_t evt = sizeof(pmtid)/sizeof(pmtid[0]);

   mcPos.SetXYZ(1760,1760,5000);

   for(size_t i = 0;i<1;i++) {
     //std:cout << " event ID "<< iEntry <<std::endl ;
     //const RAT::DS::Entry& rDS = dsReader.GetEntry( evt );
     //const RAT::DS::MC& rmc= rDS.GetMC();
     //const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
     //mcPos = rmcparticle.GetPosition();
     
     pmtpos = pmtInfo.GetPosition( pmtid[i] );
     TVector3 u = (pmtpos - mcPos).Unit();
     for(double j = 0;j<800;j = j+100) {     
      mcPos = mcPos+j*u;
     //cout<<mcPos.X()<<", "<<mcPos.Y()<<", "<<mcPos.Z()<<endl;
     //cout<<mcPos.Mag()<<endl;
     const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
     // RAT::DU::Utility::Get()->GetLightPathCalculator() must be called *after* the RAT::DU::DSReader constructor.
     RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
     lightPath.CalcByPositionPartial( mcPos, pmtpos );
     double distInInnerAV = lightPath.GetDistInInnerAV();
     double distInAV = lightPath.GetDistInAV();
     double distInWater = lightPath.GetDistInWater();
     double distInUpperTarget = lightPath.GetDistInUpperTarget();
     double distInLowerTarget = lightPath.GetDistInLowerTarget();
     RAT::DU::EffectiveVelocity effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();
     effectiveVelocity.BeginOfRun();
     const double transitTime = effectiveVelocity.CalcByDistance( distInUpperTarget, distInAV, distInWater+distInLowerTarget );
     cout<<"pmtid "<<pmtid[i]<<" ("<<pmtpos.X()<<", "<<pmtpos.Y()<<", "<<pmtpos.Z()<<") "<<"tof "<<transitTime<<endl;
     //cout<<(mcPos-pmtpos).Mag()<<", "<<endl;//<<
     cout<<transitTime<<", ";
    }
   }
   cout<<endl;
}
