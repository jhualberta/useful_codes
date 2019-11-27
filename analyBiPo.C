//For Rat-v6.2.8 SNOP data
//2017-6-25
//#include <RAT/DataCleaningUtility.hh>
#include <RAT/DS/Meta.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DataCleaningUtility.hh>
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
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = 
fit9 = FTP
*/
using namespace std ;
const double ITRval = 0.55;
double fv1 = 5700;
double fv2 = 5800;

void analyBiPo()
{
   const char* filename = "Analysis10_r0000250652_s000_p001.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos; 
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   RAT::DU::LightPathCalculator lightPath = RAT::DU::Utility::Get()->GetLightPathCalculator(); // To calculate the light's path
//const RAT::DU::GroupVelocity& groupVelocity = RAT::DU::Utility::Get()->GetGroupVelocity(); // To get the group velocity
   const RAT::DU::EffectiveVelocity &effectiveVelocity = RAT::DU::Utility::Get()->GetEffectiveVelocity();
   TH2F* hPhiTheta = new TH2F("hPhiTheta","theta vs phi",1000,-TMath::Pi(),TMath::Pi(),1000,-1.0,1.0);
   // use lightPathCalculator
   TH1F* htRes = new TH1F("htRes", "Time Residue, using lightPath, trig", 400,-100,300);
   TH1F* htRes_fvcut = new TH1F("htRes_fvcut", "Time Residue, using lightPath, fvcut", 400,-100,300);
   TH1F* htRes_nhit100 = new TH1F("htRes_nhit100", "Time Residue, using lightPath, trig", 400,-100,300);
   TH1F* htRes_nhit300 = new TH1F("htRes_nhit300", "Time Residue, using lightPath, trig", 400,-100,300);
   TH2F* htResVsPMTz_nhit100 = new TH2F("htResVsPMTz_nhit100", "Time Residue vs PMTz, using lightPath, nhit>100", 400,-100,300, 2000,-9000,9000);

   //in micro seconds !!
   TH1F* hDeltaT_bi212 = new TH1F("hDeltaT_bi212","deltaT, Bi212",1000,0,5000);
   TH1F* hDeltaT_po212 = new TH1F("hDeltaT_po212","deltaT, Po212",1000,0,5000);
   TH1F* hDeltaT_bi214 = new TH1F("hDeltaT_bi214","deltaT, Bi214",1000,0,5000);
   TH1F* hDeltaT_po214 = new TH1F("hDeltaT_po214","deltaT, Po214",1000,0,5000);
//
   TH2F* hRZ_bi212 = new TH2F("hRZ_bi212","fit r vs z, bi212", 1000, 0, 6000, 2000, 4000, 6000);
   TH2F* hRZ_po212 = new TH2F("hRZ_po212","fit r vs z, po212", 1000, 0, 6000, 2000, 4000, 6000);
   TH2F* hRZ_bi214 = new TH2F("hRZ_bi214","fit r vs z, bi214", 1000, 0, 6000, 2000, 4000, 6000);
   TH2F* hRZ_po214 = new TH2F("hRZ_po214","fit r vs z, po214", 1000, 0, 6000, 2000, 4000, 6000);

   TH2F* hYZ_bi212 = new TH2F("hYZ_bi212","fit y vs z, bi212", 2000, -6000, 6000, 2000, 4000, 6000);
   TH2F* hYZ_po212 = new TH2F("hYZ_po212","fit y vs z, po212", 2000, -6000, 6000, 2000, 4000, 6000);
   TH2F* hYZ_bi214 = new TH2F("hYZ_bi214","fit y vs z, bi214", 2000, -6000, 6000, 2000, 4000, 6000);
   TH2F* hYZ_po214 = new TH2F("hYZ_po214","fit y vs z, po214", 2000, -6000, 6000, 2000, 4000, 6000);

   TH1F* hnhits_bi212 = new TH1F("hnhits_bi212","nhits, bi212", 1000, 0, 5000);
   TH1F* hnhits_po212 = new TH1F("hnhits_po212","nhits, po212", 1000, 0, 5000);
   TH1F* hnhits_bi214 = new TH1F("hnhits_bi214","nhits, bi214", 1000, 0, 5000);
   TH1F* hnhits_po214 = new TH1F("hnhits_po214","nhits, po214", 1000, 0, 5000);

   int nhitCut = 175;//15;
   int nhitBi212 = 1200; //(175,1200)
   int nhitPo212 = 700;  //(175,700)
   int nhitBi214 = 1700; //(175,1700)
   int nhitPo214 = 320;  //(175, 320)
   // time window
   double dtBipo212 = 0.400; //400 to 800 ns
   double dtBipo214 = 3.690; //3690 to 1000000 ns

   unsigned int fecdID = 9207;
   TString fitName = "partialFitter";
   ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   size_t evtNum = dsReader.GetEntryCount(); 
   TVector3 pos_fit;
   for( size_t iEntry = 0; iEntry < evtNum; iEntry++ )
   {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      Int_t nevC = rDS.GetEVCount();
      if(nevC == 1) {
      //for(Int_t iev=0;iev<nevC; iev++) {
       /// Get the Event information    
       //const RAT::DS::EV& rev = rDS.GetEV(iev);
       const RAT::DS::EV& rev = rDS.GetEV(0);
       std::vector<std::string> fitname = rev.GetFitNames();
       std::vector<std::string>::iterator it;
       int trig = rev.GetTrigType();//get triggerType
       const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
       RAT::DS::DataQCFlags dcflags = rev.GetDataCleaningFlags();
       bool tagBi212 = false, tagPo212 = false, tagBi214 = false, tagPo214 = false;
       if( RAT::EventIsClean( rev, dcAnalysisWord ) )
       {
         double nhits = rev.GetNhitsCleaned();//NOTE: use nhits_cleaned
         if(nhits>nhitCut) {
            if(rev.FitResultExists("partialFitter")){//find partialFitter exists
            try{
                 RAT::DS::FitVertex fitVertex = rev.GetFitResult("partialFitter").GetVertex(0);
                 if(fitVertex.ValidPosition())
                 {
                   pos_fit=fitVertex.GetPosition();
                   double tfit = fitVertex.GetTime();
                   double posX=(pos_fit.X()); double posY=(pos_fit.Y()); double posZ=(pos_fit.Z());
                   Double_t day = rev.GetUniversalTime().GetDays();
                   Double_t second = rev.GetUniversalTime().GetSeconds();
                   Double_t nanosecond= rev.GetUniversalTime().GetNanoSeconds();
                   //double timePri = (day*24*3600*1e9 + second*1e9 + nanosecond)/1000;
                   double timePri = rev.GetClockCount50()/1000;
		   //std::cout<<day<<" "<<second<<" "<<nanosecond<<" "<<timePri<<std::endl;
		   if(sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108))<fv1 && posZ-108>5100)
                   {
		     // event pair
                     double deltaT = 0;
      	             int count = 0;//how many events checking next
		     for( size_t jEntry = iEntry+1; jEntry < evtNum; jEntry++ ) // delayed events
                     {  
                        const RAT::DS::Entry& rDely = dsReader.GetEntry( jEntry );
                        Int_t nevC2 = rDely.GetEVCount();
                        if(nevC2 == 1)//for(Int_t iev2=0;iev2<nevC2; iev2++) 
			{
                         //std::cout<<iEntry<<" "<<jEntry<<std::endl;
	     		/// Get the Event information    
                        const RAT::DS::EV& revDely = rDely.GetEV(0);
                        int trig = revDely.GetTrigType();//get triggerType
                        //const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
                        RAT::DS::DataQCFlags dcflags = revDely.GetDataCleaningFlags();
                        if( RAT::EventIsClean( revDely, dcAnalysisWord ) )
                        {
                          double nhits2 = revDely.GetNhitsCleaned();//NOTE: use nhits_cleaned
                          if(nhits2>nhitCut) {
                             if(revDely.FitResultExists("partialFitter")){//find partialFitter exists
                             try{
                                  RAT::DS::FitVertex fitVertex2 = revDely.GetFitResult("partialFitter").GetVertex(0);
                                  if(fitVertex2.ValidPosition())
                                  {
				     TVector3 pos_fit2 = fitVertex2.GetPosition();
                                     //Double_t day2 = revDely.GetUniversalTime().GetDays();
                                     //Double_t second2 = revDely.GetUniversalTime().GetSeconds();
                                     //Double_t nanosecond2 = revDely.GetUniversalTime().GetNanoSeconds();
				     //double timeDely = (day2*24*3600*1e9 + second2*1e9 + nanosecond2)/1000;
                                     double timeDely = rev.GetClockCount50()/1000;
				     deltaT = ( long(timeDely - timePri) & 0x7FFFFFFFFFF )*20/1000; //micro seconds
                                     std::cout<<deltaT<<" "<<timeDely<<" "<<timePri<<std::endl;
				     if(deltaT>dtBipo212 && (pos_fit-pos_fit2).Mag()<500 ) // dR cut
			             {  
					//count++;
					std::cout<<"deltaT = " <<deltaT<<std::endl;
                                        if(deltaT>dtBipo214 && deltaT<1000) // tag bipo214
                                        {
				          if(nhits<nhitBi214) { // tag bi214
                                            tagBi214 = true;
                                            size_t gtid_bi214 = rev.GetGTID();
                                            cout<<"bi214 "<<gtid_bi214<<endl;
					  }
                                          if(nhits<nhitPo214) { // tag Po214
                                            tagPo214 = true;
					    size_t gtid_po214 = rev.GetGTID();
					    cout<<"Po214 "<<gtid_po214<<endl;
					  }
					  if(deltaT<dtBipo212+400./1000) // tag bipo212
					  {     
					   if(nhits<nhitBi212) 
					   {
				             tagBi212 = true;
                                           }
					   if(nhits<nhitPo212)
                                           {
					     tagPo214 = true; 
					   }
					   std::cout<<"delay time "<<deltaT<<std::endl;
					   break;
					  }
				        }
			                else if(deltaT>1000) //beyond time
				        {break;}
				     }
                                    }//valid position
                                  }//try catch
                             catch(exception& e)
                             {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
                           } //partialFitter
			 } //nhit cut
			}// event2 is clean
                      } //2nd triggered
                       if(tagBi212) {
		         hRZ_bi212->Fill(pos_fit.Perp(),pos_fit.Z());
                         hYZ_bi212->Fill(pos_fit.Y(),pos_fit.Z());   
                         hDeltaT_bi212->Fill(deltaT);
                         hnhits_bi212->Fill(nhits);
		       }
                       if(tagPo212) {
                         hRZ_po212->Fill(pos_fit.Perp(),pos_fit.Z());
                         hYZ_po212->Fill(pos_fit.Y(),pos_fit.Z());   
                         hDeltaT_po212->Fill(deltaT);
			 hnhits_po212->Fill(nhits);
                       }
                       if(tagPo214) {
                         hRZ_po214->Fill(pos_fit.Perp(),pos_fit.Z());
                         hYZ_po214->Fill(pos_fit.Y(),pos_fit.Z());
                         hDeltaT_po214->Fill(deltaT);
			 hnhits_po214->Fill(nhits);
                       }
		       if(tagBi214) {
                         hRZ_bi214->Fill(pos_fit.Perp(),pos_fit.Z());
                         hYZ_bi214->Fill(pos_fit.Y(),pos_fit.Z());   
                         hDeltaT_bi214->Fill(deltaT);
			 hnhits_bi214->Fill(nhits);
                       }
                       //if(count>10){ std::cout<<"break"<<std::endl; break;}
		     }//loop 2nd event
                   }//FV cut
  		   }//valid position
             }//try catch
             catch(exception& e)
             {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
            } // partialFitter exists
         }//nhit>175
       } // dataClean 
      } // triggered events
   }//loop events

   TString newfilename = "ResolBiPo_"+TString(filename);
   TFile *fp = new TFile(newfilename,"recreate");
   fp->cd();
		   
   hDeltaT_bi212->Write(); hDeltaT_po212->Write(); hDeltaT_bi214->Write(); hDeltaT_po214->Write();
   hRZ_bi212->Write(); hRZ_po212->Write(); hRZ_bi214->Write(); hRZ_po214->Write();
   hYZ_bi212->Write(); hYZ_po212->Write(); hYZ_bi214->Write(); hYZ_po214->Write();
   hnhits_bi212->Write();hnhits_po212->Write();hnhits_bi214->Write();hnhits_po214->Write();

}
