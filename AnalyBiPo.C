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
 
using namespace std ;
const double ITRval = 0.55;
double fv1 = 5700;
double fv2 = 5800;
double deltaRcut = 500;
long int outTimeWin = 2000; // 2000 micro s
void AnalyBiPo()
{
   //gSystem->cd("/scratch/djauty/data/full/Jan/6_16_11/");
   const char* filename = "FitNormal_bipo214_0p5PPO_IsoTop_splitZ4500.root";
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
   TH1F* hNhitClean = new TH1F("hNhitClean","nhitCleaned, all triggered",2000,0,2000); //all fitted, nhitCleaned
   TH1F* hDeltaR = new TH1F("hDeltaR","deltaR",1000,0,9000);
   //in micro seconds !!
   TH1F* hDeltaT = new TH1F("hDeltaT","deltaT",1000,0,2000); //us
   TH1F* hDeltaT_bi212 = new TH1F("hDeltaT_bi212","deltaT, Bi212",1000,0,5);
   TH1F* hDeltaT_po212 = new TH1F("hDeltaT_po212","deltaT, Po212",1000,0,5);
   TH1F* hDeltaT_bi214 = new TH1F("hDeltaT_bi214","deltaT, Bi214",2000,0,2000);
   TH1F* hDeltaT_po214 = new TH1F("hDeltaT_po214","deltaT, Po214",2000,0,2000);
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

   int nhitCut = 100;//15;
   int nhitBi212 = 1200; //(175,1200)
   int nhitPo212 = 700;  //(175,700)
   int nhitBi214 = 1700; //(175,1700)
   int nhitPo214 = 320;  //(175, 320)
   // time window
//   double dtBipo212 = 0.400; //400 to 800 ns
//   double dtBipo214 = 3.690; //3.69 to 1000 us

   unsigned int fecdID = 9207;
   TString fitName = "partialFitter";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   size_t evtNum = dsReader.GetEntryCount(); 
   TVector3 pos_fit, pos_fit2;
   TTree *Tpair = new TTree("T","event pairs");

   TTree *TpairBipo212 = new TTree("TBipo212","event pairs, bipo212");
   TTree *TpairBipo214 = new TTree("TBipo214","event pairs, bipo214");

   long int evtPriID = 0;
   long int nhitsPri = 0;
   double posxPri = 0, posyPri = 0, poszPri = 0, tPri = 0;

   long int evtDelyID = 0; 
   long int nhitsDely = 0;
   double posxDely = 0, posyDely = 0, poszDely = 0, tDely = 0;

   Tpair->Branch("evtPriID",&evtPriID,"evtPriID/l");
   Tpair->Branch("nhitPri",&nhitsPri,"nhitsPri/l");
   Tpair->Branch("posxPri",&posxPri,"posxPri/D");
   Tpair->Branch("posyPri",&posyPri,"posyPri/D");
   Tpair->Branch("poszPri",&poszPri,"poszPri/D");
   Tpair->Branch("tPri",&tPri,"tPri/D");

   Tpair->Branch("evtDelyID",&evtDelyID,"evtDelyID/l");
   Tpair->Branch("nhitDely",&nhitsDely,"nhitsDely/l");
   Tpair->Branch("posxDely",&posxDely,"posxDely/D");
   Tpair->Branch("posyDely",&posyDely,"posyDely/D");
   Tpair->Branch("poszDely",&poszDely,"poszDely/D");
   Tpair->Branch("tDely",&tDely,"tDely/D");
   //
   //
   TpairBipo212->Branch("evtPriID",&evtPriID,"evtPriID/l");
   TpairBipo212->Branch("nhitPri",&nhitsPri,"nhitsPri/l");
   TpairBipo212->Branch("posxPri",&posxPri,"posxPri/D");
   TpairBipo212->Branch("posyPri",&posyPri,"posyPri/D");
   TpairBipo212->Branch("poszPri",&poszPri,"poszPri/D");
   TpairBipo212->Branch("tPri",&tPri,"tPri/D");
   ///
   TpairBipo212->Branch("evtDelyID",&evtDelyID,"evtDelyID/l");
   TpairBipo212->Branch("nhitDely",&nhitsDely,"nhitsDely/l");
   TpairBipo212->Branch("posxDely",&posxDely,"posxDely/D");
   TpairBipo212->Branch("posyDely",&posyDely,"posyDely/D");
   TpairBipo212->Branch("poszDely",&poszDely,"poszDely/D");
   TpairBipo212->Branch("tDely",&tDely,"tDely/D");
   //
   //
   TpairBipo214->Branch("evtPriID",&evtPriID,"evtPriID/l");
   TpairBipo214->Branch("nhitPri",&nhitsPri,"nhitsPri/l");
   TpairBipo214->Branch("posxPri",&posxPri,"posxPri/D");
   TpairBipo214->Branch("posyPri",&posyPri,"posyPri/D");
   TpairBipo214->Branch("poszPri",&poszPri,"poszPri/D");
   TpairBipo214->Branch("tPri",&tPri,"tPri/D");
   ///
   TpairBipo214->Branch("evtDelyID",&evtDelyID,"evtDelyID/l");
   TpairBipo214->Branch("nhitDely",&nhitsDely,"nhitsDely/l");
   TpairBipo214->Branch("posxDely",&posxDely,"posxDely/D");
   TpairBipo214->Branch("posyDely",&posyDely,"posyDely/D");
   TpairBipo214->Branch("poszDely",&poszDely,"poszDely/D");
   TpairBipo214->Branch("tDely",&tDely,"tDely/D");

   for( size_t iEntry = 0; iEntry < evtNum; iEntry++ )
   {
      const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
      Int_t nevC = rDS.GetEVCount();
      for(Int_t iev=0;iev<nevC; iev++) {
       /// Get the Event information    
       const RAT::DS::EV& rev = rDS.GetEV(0);//iev);
       //std::vector<std::string> fitname = rev.GetFitNames();
       //std::vector<std::string>::iterator it;
       int trig = rev.GetTrigType();//get triggerType
       const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
       bool tagBi212 = false, tagPo212 = false, tagBi214 = false, tagPo214 = false;
       //RAT::DS::DataQCFlags dcflags = rev.GetDataCleaningFlags();
       //if( RAT::EventIsClean( rev, dcAnalysisWord ) )
       {
         double nhits = rev.GetNhitsCleaned();//NOTE: use nhits_cleaned
         hNhitClean->Fill(nhits);
	 if(nhits>nhitCut) {
            if(rev.FitResultExists("partialFitter")) {//find partialFitter exists
            try{
                 RAT::DS::FitVertex fitVertex = rev.GetFitResult("partialFitter").GetVertex(0);
                 if(fitVertex.ValidPosition())
                 {
                   pos_fit = fitVertex.GetPosition();
                   double tfit = fitVertex.GetTime();
                   double posX=(pos_fit.X()); double posY=(pos_fit.Y()); double posZ=(pos_fit.Z());
                   Double_t day = rev.GetUniversalTime().GetDays();
                   Double_t second = rev.GetUniversalTime().GetSeconds();
                   Double_t nanosecond= rev.GetUniversalTime().GetNanoSeconds();
                   Double_t timePri0 = (day*24*3600*1e9 + second*1e9 + nanosecond);
                   long int timePri = double(rev.GetClockCount50()); //!!! go to seconds
		   //std::cout<<day<<" "<<second<<" "<<nanosecond<<" "<<timePri<<std::endl;
		   if(sqrt(posX*posX+posY*posY+(posZ-108)*(posZ-108))<fv1 && posZ-108>5100)
                   {
		     // event pair
		     size_t evtPri = rev.GetGTID();
                     double deltaT0 = 0;
		     double deltaT = 0;
                     int count = 0;//how many events checking next
                     bool timeOut = false;
                     /// loop 2nd events
		     for( size_t jEntry = iEntry+1; jEntry < evtNum; jEntry++ ) // delayed events
                     {
                        const RAT::DS::Entry& rDely = dsReader.GetEntry( jEntry );
                        Int_t nevC2 = rDely.GetEVCount();
                        //if(nevC2==1)
			for(Int_t iev2=0;iev2<nevC2; iev2++)
			{
                         //std::cout<<iEntry<<" "<<jEntry<<std::endl;
	     		/// Get the Event information    
                        const RAT::DS::EV& revDely = rDely.GetEV(iev2);
                        int trig = revDely.GetTrigType();//get triggerType
                        double nhits2 = revDely.GetNhitsCleaned();//NOTE: use nhits_cleaned
                        //const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
                        //RAT::DS::DataQCFlags dcflags = revDely.GetDataCleaningFlags();
                        //if( RAT::EventIsClean( revDely, dcAnalysisWord ) )
                        {
                          if(nhits2>nhitCut) {
                             if(revDely.FitResultExists("partialFitter")){//find partialFitter exists
                             try{
                                  RAT::DS::FitVertex fitVertex2 = revDely.GetFitResult("partialFitter").GetVertex(0);
                                  if(fitVertex2.ValidPosition())
                                  {
				     TVector3 pos_fit2 = fitVertex2.GetPosition();
                                     Double_t day2 = revDely.GetUniversalTime().GetDays();
                                     Double_t second2 = revDely.GetUniversalTime().GetSeconds();
                                     Double_t nanosecond2 = revDely.GetUniversalTime().GetNanoSeconds();
				     Double_t timeDely0 = (day2*24*3600*1e9 + second2*1e9 + nanosecond2);
                                     long int timeDely = revDely.GetClockCount50();
                                     size_t evtDely = revDely.GetGTID();
				     //deltaT = ( (timeDely - timePri) & 0x7FFFFFFFFFF )*20/1e3; //micro seconds
				     deltaT = ( (timeDely - timePri) )*20/1e3; //micro seconds
                                     deltaT0 = (timeDely0 - timePri0)/1e3; //universal time for comparison
				     //std::cout<<"universal "<<deltaT0<<" 50MHz clock "<<deltaT<<" "<<timeDely<<" "<<timePri<<std::endl;
                                     hDeltaT->Fill(deltaT); hDeltaR->Fill( (pos_fit-pos_fit2).Mag() );
				     if(deltaT>outTimeWin) timeOut = true;
				     else if( deltaT>0.4 && deltaT<outTimeWin && (pos_fit-pos_fit2).Mag()<deltaRcut) // dR cut
			             {
					//count++;
                                       evtPriID = evtPri; 
                                       nhitsPri = nhits;
                                       posxPri = pos_fit.X(), posyPri = pos_fit.Y(), poszPri = pos_fit.Z();
                                       tPri = timePri*20/1e3; // us

                                       evtDelyID = evtDely;
                                       nhitsDely = nhits2;
                                       posxDely = pos_fit2.X(), posyDely = pos_fit2.Y(), poszDely = pos_fit2.Z();
                                       tDely = timeDely*20/1e3; // us

				       Tpair->Fill();
                                       //std::cout<<"event pairID "<<evtPri<<", "<<evtDely<<std::endl;
                                       //std::cout<<"event pairNhits "<<nhits<<", "<<nhits2<<std::endl;
				       //std::cout<<"deltaT = " <<deltaT<<std::endl;
                                       if(deltaT>3.69 && deltaT<1800) // tag bipo214 (3.69 to 1800 us)
                                       {
				         if(nhits<nhitBi214) { // tag bi214
                                           tagBi214 = true;
                                           size_t gtid_bi214 = rev.GetGTID();
                                           hRZ_bi214->Fill(pos_fit.Perp(),pos_fit.Z());
                                           hYZ_bi214->Fill(pos_fit.Y(),pos_fit.Z());
                                           hDeltaT_bi214->Fill(deltaT);
                                           hnhits_bi214->Fill(nhits);
					   // cout<<"bi214,id "<<gtid_bi214<<endl;
				         }
                                         if(nhits2<nhitPo214) { // tag Po214
                                           tagPo214 = true;
				           size_t gtid_po214 = rev.GetGTID();
                                           hRZ_po214->Fill(pos_fit2.Perp(),pos_fit2.Z());
                                           hYZ_po214->Fill(pos_fit2.Y(),pos_fit2.Z());
                                           hDeltaT_po214->Fill(deltaT);
                                           hnhits_po214->Fill(nhits2);
					   // cout<<"Po214,id "<<gtid_po214<<endl;
				         }
                                         TpairBipo214->Fill();
                                       } //214Bipo (3.69, 1800) us
				       else if(deltaT<3.69) // tag bipo212
				       {     
				          if(nhits<nhitBi212) 
				          {
				            tagBi212 = true;
                                            hRZ_bi212->Fill(pos_fit.Perp(),pos_fit.Z());
                                            hYZ_bi212->Fill(pos_fit.Y(),pos_fit.Z());
                                            hDeltaT_bi212->Fill(deltaT);
                                            hnhits_bi212->Fill(nhits);
                                          }
				          if(nhits2<nhitPo212)
                                          {
				            tagPo212 = true;
                                            hRZ_po212->Fill(pos_fit2.Perp(),pos_fit2.Z());
                                            hYZ_po212->Fill(pos_fit2.Y(),pos_fit2.Z());
                                            hDeltaT_po212->Fill(deltaT);
                                            hnhits_po212->Fill(nhits2);
				          }
					  TpairBipo212->Fill();
				          //std::cout<<"delay time "<<deltaT<<std::endl;
				       } //212Bipo (0.4, 3.69) us
				     else {break;}
			            } //deltaT>0.4
                                    }//valid position
                                  }//try catch
                             catch(exception& e)
                             {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
                           } //partialFitter
			 } //nhit cut
			}// event2 is clean
                      } //2nd triggered
                      if(timeOut){ break; }
		      //if(!tagBi214 && !tagPo214 && !tagBi212 && !tagPo212) {break;}
		     }//loop 2nd event
                   }//1st evt FV cut
                 }//evt1, fitValid position
             }//try catch
             catch(exception& e)
             {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
            } // partialFitter exists
         }//nhit>175
       } // dataClean 
      } // triggered events
   }//loop events

   TString newfilename = "ResolBiPo_"+TString(filename);
   //gSystem->cd("/home/jhu9/scratch/partialAnalysis/");
   TFile *fp = new TFile(newfilename,"recreate");
   fp->cd();
   hNhitClean->Write();
   hDeltaT->Write();hDeltaR->Write();Tpair->Write();TpairBipo212->Write();TpairBipo214->Write(); 
   hDeltaT_bi212->Write(); hDeltaT_po212->Write(); hDeltaT_bi214->Write(); hDeltaT_po214->Write();
   hRZ_bi212->Write(); hRZ_po212->Write(); hRZ_bi214->Write(); hRZ_po214->Write();
   hYZ_bi212->Write(); hYZ_po212->Write(); hYZ_bi214->Write(); hYZ_po214->Write();
   hnhits_bi212->Write();hnhits_po212->Write();hnhits_bi214->Write();hnhits_po214->Write();
   fp->Close();
}
