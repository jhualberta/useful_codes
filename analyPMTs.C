//For Rat-v6.2.8 SNOP data
//2017-6-25
//#include <RAT/DataCleaningUtility.hh>
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
//NOTE: no MC data for real data!! otherwise overflow error
//for SNO+ data
//Unit: mm
/*
fit10 = 
fit9 = FTP
*/
using namespace std ;
const double ITRval = 0.55;
void analyPMTs()
{
   const char* filename = "FitRat_MC_partialscint_2p5MeVbeta_x0_y0_z1000_level0.root";
   TFile *f1 = new TFile(filename);
   RAT::DU::DSReader dsReader(filename);
   TVector3 sourcePos;
   int panel = 16; 
   bool scint = true;
   sourcePos.SetXYZ(0,0,0);
   const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
   TVector3 u_fit,pos_fit, pos_fit_ItrBeta14, u_fit_pos_ItrBeta14;
   Double_t theta_e;
   //u_e.SetXYZ(1,0,0);//initial electron direction
   Double_t grVelocity = 2.17554021555098529e+02 ;//light water:2.17554021555098529e+02; heavy water:2.18254558686904687e+02
   if(scint) grVelocity = 183.503698029;//labppo 
   Double_t rPSUP = 8900;
   string process1 = "Scintillation" ;
   string process2 = "Cerenkov" ;
   string process3 = "OpAbsorption" ;
   Double_t energy, wavelength ;
 
   Double_t countFitValid = 0;//total fitted events
   Double_t countFECDtotal =0;//total FECD==9188 events
   Double_t countSuccess = 0;
   Double_t countTimeWindow = 0;//fitValid && FECD && tFit cuts
   Double_t countClean = 0;
   Double_t countTrig = 0;
   Double_t countTrigFECD = 0;
   Double_t countITRbeta14 = 0;
   Double_t countNoCleanTotal=0; Double_t countNoCleanTrig=0;Double_t countNoCleanFECD=0;
   int trigWord = 0x11;//trig number,1<<6, http://www.snoplus.ca/docs/rat/user_manual/html/node47.html#t:trigword
   bool trigCut = true;//true;
   int nhitCut = 0;//15;
   unsigned int fecdID = 9190;
   TString fitName = "partialFitter";
   //TString specFit = "fit10";
   //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );
   double countAll;
   double countUp;
   TH1F *hPMTdist = new TH1F("hPMTdist","pmt distribution ratio",100,0,1);
   TH1F *htRes = new TH1F("htRes","time res",400,-100,300);
   TH1F *htRes_early = new TH1F("htRes_early","time res, early",400,-100,300);

   vector<TH1F*> hlist;
   vector<TH2F*> hlist2;

   TVector3 pos_mc;
   for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
   {
    //std:cout << " event ID "<< iEntry <<std::endl ;
    const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
    Int_t nevC = rDS.GetEVCount();//!! if retriggered events, nevC == 2
    const RAT::DS::MC& rmc= rDS.GetMC();  
    const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
    pos_mc =rmcparticle.GetPosition();
    double tmc = rmcparticle.GetTime();
    //cout<<"trigger events "<<nevC<<endl;
    if(nevC>1) nevC = 1;//!!! a retrigger cut
    for(Int_t iev=0;iev<nevC; iev++){
     /// Get the Event information    
     const RAT::DS::EV& rev = rDS.GetEV(iev);
     std::vector<std::string> fitname = rev.GetFitNames();
     std::vector<std::string>::iterator it;
     double nhits = rev.GetNhits();
     int trig = rev.GetTrigType();//get triggerType
     const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
     countNoCleanTotal++;
     countUp = 0; countAll = 0;
     countAll = calpmts.GetCount();

     size_t totalbin = panel*panel;
     int index1, index2, index3;//costheta, phi, 2-D table indices
     double step1 = 2./panel;
     double step2 = 2*TMath::Pi()/panel;
     vector<vector<double> > group_hitTime(totalbin+1);
     vector<vector<int> > group_index(totalbin+1);
     group_hitTime.clear();
     vector<int> selectId;
     selectId.clear();
     if(rev.FitResultExists("partialFitter")) {//find partialFitter exists
        RAT::DS::FitVertex fitVertex = rev.GetFitResult("partialFitter").GetVertex(0);
        if(rev.GetFitResult("partialFitter").GetValid())//Global Validity of the fitter !!! NOTE: this is different
        {
          pos_fit=fitVertex.GetPosition(); 
	  //u_fit=fitVertex.GetDirection().Unit();
 	  double tfit = fitVertex.GetTime();
	  if(fitVertex.ValidPosition())//NOTE: check fitPosition Valid !!
 	  { 
	    double posX=(pos_fit.X());double posY=(pos_fit.Y());double posZ=(pos_fit.Z());

            for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
            {
             TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
             if(pmtpos.Z()>0) countUp++;
             double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
             double costheta = pmtpos.CosTheta();
             double phi = pmtpos.Phi();
             index1 = floor((costheta-(-1))/step1);
             index2 = floor((phi-(-TMath::Pi()))/step2);
             index3 = index1*16 + index2;
             //cout<<"Loop PMTs "<<endl;
             //cout<<costheta<<" "<<phi<<" "<<index1<<" "<<index2<<" "<<index3<<endl;
             group_index[index3].push_back(ipmt);
             group_hitTime[index3].push_back(hitTime);
             double tres = hitTime - (pos_fit - pmtpos).Mag()/grVelocity - tfit;
             htRes->Fill(tres);
            }
            hPMTdist->Fill(countUp/countAll);
 
            vector<double> temp_time;
            for (size_t ind1 = 0; ind1 < totalbin; ind1 ++) {//group pmt and select
             //TH1F *htemp = new TH1F("htemp","",500,0,500);
             //TH2F *htemp2 = new TH2F("htemp2","",16,-1,1,16,-TMath::Pi(),TMath::Pi());
             if(!group_hitTime[ind1].empty()) {
              //cout<<"panel: "<<ind1<<" hitTime ";
              for (size_t ind2 = 0;ind2 < group_hitTime[ind1].size(); ind2 ++) {// ind2
                temp_time.push_back(group_hitTime[ind1][ind2]);
                //cout<<group_hitTime[ind1][ind2]<<", ";
                //htemp->Fill(group_hitTime[ind1][ind2]);
                TVector3 tempV = pmtInfo.GetPosition(calpmts.GetPMT(group_index[ind1][ind2]).GetID());
                //htemp2->Fill( tempV.CosTheta(), tempV.Theta(), group_hitTime[ind1][ind2]);
              }
              vector<double>::iterator iterTime = min_element(temp_time.begin(), temp_time.end());
              //std::cout<<"min "<<*iterTime<<endl;
              int IndexMin = iterTime - temp_time.begin();
  
              // calculate median
              vector<double> copy_temp = temp_time;
              vector<double>::iterator middle = copy_temp.begin() + (copy_temp.size() / 2);
              // this sets iterator middle to the median element
              std::nth_element(copy_temp.begin(), middle, copy_temp.end());
              double median = *middle;
              vector<double>::iterator itOfMedian = std::find(temp_time.begin(), temp_time.end(), median);
              int IndexMedian = itOfMedian - temp_time.begin();
  
              /// an accurate way to calculate median
              //vector<double>::iterator median_it1 = temp_time.begin()+temp_time.size()/2 - 1;
              //vector<double>::iterator median_it2 = temp_time.begin()+temp_time.size()/2;
              //std::nth_element(temp_time.begin(), median_it1 , temp_time.end()); // e1
              //std::nth_element(temp_time.begin(), median_it2 , temp_time.end()); // e2
              //double median = (temp_time.size()%2 == 0) ? (*median_it1 + *median_it2)/2 : *median_it2;
              //std::cout<<"median "<<median<<endl;
  
              int selectIndex = group_index[ind1][IndexMin];
              selectId.push_back(selectIndex);
              //std::cout<<"index min "<<IndexMin<<" index median "<<IndexMedian<<std::endl; 
              //hlist.push_back((TH1F*)htemp->Clone(Form("hpanel_%u",ind1)));//%u for size_t
              //hlist2.push_back((TH2F*)htemp2->Clone(Form("hpanel2_%lu",ind1)));//%u for size_t
              //delete htemp; delete htemp2;
              temp_time.clear();
            }
            else{   
              //delete htemp;delete htemp2;
              temp_time.clear();
              continue;
             }
         } //group pmt and select 
         int i = 0; 
         for(vector<int>::iterator it = selectId.begin(); it!= selectId.end(); it++)
         {
           TVector3 pmtpos = pmtInfo.GetPosition(*it);
           if(pmtpos.Z()>0) countUp++;
           double hitTime =(calpmts.GetPMT(i)).GetTime();
           i++;
           double tres_early = hitTime - (pos_fit - pmtpos).Mag()/grVelocity - tfit; 
           htRes_early->Fill(tres_early);
         }
      }//valid position
     }//valid global 
   }//fit exits
   }//if triggered events
  }//for loop



   TString newfile(filename);
   TFile *fnew = new TFile("PMT_selector.root","recreate");
   fnew->cd();
   hPMTdist->Write();htRes->Write();htRes_early->Write();
   // for (size_t j=0;j<hlist2.size();j++)
   // { 
   //   //hlist[j]->Write();
   //   hlist2[j]->Write();
   // }
   fnew->Close();
}
