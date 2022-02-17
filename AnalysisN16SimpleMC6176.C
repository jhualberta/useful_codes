// This code loops over all of the events in a given root file and makes plots
// It is run as a function analysisSolar.cc+(path/inputFileName.root, path/outputFileName.root)
// Run it as root 'AnalysisN16SimpleMC6176.C+("AntiNeutrinoMC_WaterN16sourceRun_r107055_s0_p10.root")'
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
#include <RAT/DataCleaningUtility.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <vector>
//DATA
void AnalysisN16SimpleMC6176(const char* filename)
{
  const double c_light = 299.792458;
  double waterRI = 1.38486;
  double grVelocity = c_light/waterRI;
  double nhitCut = 25;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  //!!!Load source manipulate position here
  //Run 107055 (-5.283	-0.209	-1.057)
  TVector3 srcPos(-5.283,-0.209,-1.057);
  //NOTE!!!! energy correction for Rat-6.17.x
  const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
  // Define all of the histograms and ntuples that are to be made  
  //MPW no drive correction, with drive correction and ITR+beta14 cuts
  TH1D *hITR = new TH1D("hITR","trig+FECD, ITR",100,0,1.0);
  TH1D *hBeta14 = new TH1D("hBeta14","trig+FECD, beta14",200,-1.0,4.0);
  TH1D *hThetaij = new TH1D("hThetaij","trig+FECD, thetaij",200,-1.0,4.0);
  TH2D *hITR_vs_Beta14 = new TH2D("hITR_vs_Beta14","trig+FECD, itr vs beta14",100,0,1.0,200,-1.0,4.0);
  TH1D *hTrig = new TH1D("hTrig","trigger type",50,0,50);  
 
  TH1D* radius = new TH1D( "radius", "Radius of all triggered events", 100, 0., 15000.);
  
  TH1D* nHit = new TH1D( "nHit", "nHit of all triggered events", 1000, 0.5, 1000.5);   
  TH1D* nHitAV = new TH1D( "nHitAV", "nHit of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5);

  TTree *tree= new TTree("T","Simulation");
  
  double posx, posy, posz, posRad, energyOrigin, energy, time, dirx, diry, dirz, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, thij;
  double posxcor, posycor, poszcor, posRadCor;
  double posTheta, posPhi;//for berkeley blob  
  double prompt_cut; 
  UInt_t nhits, day, sec, runNumber, subRun, eventGTID, fecdID, triggerWord;
  double nsecs;
  double Gtest, Utest, posFOM, posFOM2, scaleLogL, dirFOM, dirFOM2, dirscaleLogL;

  //std::vector<double> vtRes;
  std::vector<UInt_t> vfecd;   

  //std::vector<double>* pvtRes = &vtRes;
  std::vector<UInt_t>* pvfecd = &vfecd; 
  tree->Branch("posx",&posx,"posx/D");
  tree->Branch("posy",&posy,"posy/D");
  tree->Branch("posz",&posz,"posz/D");
  tree->Branch("posTheta",&posTheta,"posTheta/D");
  tree->Branch("posPhi",&posPhi,"posPhi/D");
  tree->Branch("time",&time,"time/D");
  //tree->Branch("timeRes",&vtRes); // time residual
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("dirx",&dirx,"dirx/D");
  tree->Branch("diry",&diry,"diry/D");
  tree->Branch("dirz",&dirz,"dirz/D");
  tree->Branch("nhits",&nhits,"nhits/i");
  tree->Branch("energyOrigin",&energyOrigin,"energyOrigin/D");
  tree->Branch("energy",&energy,"energy/D");
  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("subRun", &subRun, "subRun/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("triggerWord",&triggerWord,"triggerWord/I");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("beta14",&beta14,"beta14/D");
  tree->Branch("thetaij",&thij,"thij/D");
  tree->Branch("Gtest",&Gtest,"Gtest/D");
  tree->Branch("Utest",&Utest,"Utest/D");
  tree->Branch("posFOM",&posFOM,"posFOM/D");
  tree->Branch("posFOM2",&posFOM2,"posFOM2/D");
  tree->Branch("scaleLogL",&scaleLogL,"scaleLogL/D");
  tree->Branch("dirFOM",&dirFOM,"dirFOM/D");
  tree->Branch("dirFOM2",&dirFOM2,"dirFOM2/D");
  tree->Branch("dirscaleLogL",&dirscaleLogL,"dirscaleLogL/D");
  tree->Branch("prompt_cut",&prompt_cut,"prompt_cut/D");
  //char filename[1500];
  Int_t countValid = 0; 
  Int_t i=0;
  UInt_t TriggerType;
  Int_t nentries=0;
  ULong64_t bitword;
  cout<<" filename "<<filename<<endl;
    // Load the RAT file
  RAT::DU::DSReader dsReader(filename);
  
  // To later get the run info
  const RAT::DS::Run& run = dsReader.GetRun();
  
  size_t numEvents = dsReader.GetEntryCount();

  cout<<numEvents<<endl; 
  // Loop over all of the events
  for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
  {
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
     // const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
     
     size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
     ///!!! this is dataClean cuts, only for data, not MC     
     //ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

     if(iEntry%10000 == 0) {cout<<iEntry<<endl;}
      
  	// Looping over triggered events in each ds
      if (iEVCount) {
       for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
       {
           // Get the event
           const RAT::DS::EV& ev = rDS.GetEV(iEV);
           const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
           runNumber = run.GetRunID();
           subRun = run.GetSubRunID();
           eventGTID = ev.GetGTID();
           //RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
           //if( RAT::EventIsClean( ev, dcAnalysisWord ) )
           {	  
      	  // Fill the nHit histogram, all events should have nHit
      	  Float_t nHitOfEvent = ev.GetNhitsCleaned();
      	  nHit->Fill(nHitOfEvent);
      	  nhits= nHitOfEvent;
      	  TriggerType = ev.GetTrigType();
          //pvtRes->clear();
      	  // Now get the radius, if applicable
      	  Double_t radiusOfEvent = -1., radiusOfEventCor = -1.;
      	  TVector3 pos_fit, pos_cor;
      	  TVector3 u_fit;
      	  const string fitName = "waterFitter"; // reconstruction fitter name
      	  // this FECD is for N16 "source tag", after RAT-6.17.6, this is automatically applied on data and MC
          pvfecd->clear();
	  for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
          {
                  fecdID = calpmts.GetFECDPMT(ipmt).GetID();pvfecd->push_back(fecdID);
          }
          bool checkFitExist = ev.FitResultExists(fitName);
          // trigger cuts, not used here
      	  //if(((TriggerType&2)==2) || ((TriggerType&18)==18) )
      	  {
      	    if(checkFitExist) // reconstruction valid
      	    {  // reconstruction fitter name needs to exist
                   try {//!!! fVertex somehow broken here 
                   RAT::DS::FitVertex fVertex = ev.GetFitResult("waterFitter").GetVertex(0);
                   
      	     if( fVertex.ValidPosition() ) // valid reconstruction position 
      	     {
                    //std::cout<<"eventGTID "<<eventGTID<<" "<<iEntry<<std::endl;
                    pos_fit = fVertex.GetPosition(); 
        
                    RAT::DS::FitVertex fVertex1 =ev.GetFitResult("waterFitter").GetVertex(0);
                    if(fVertex1.ValidDirection()) { // valid reconstruction direction
                     try{
                       u_fit = fVertex1.GetDirection();
                     }
                     catch(exception& e)
                     { std::cout<<e.what()<<" problems in dir fit"<<std::endl;}
                      //NOTE: drive correction here
      	              //pos_fit = driveP0*pos_fit + driveP1*u_fit;
                     countValid++;
                     time = fVertex.GetTime();
      	             //!!! z - 108
                     posx = pos_fit.X();posy = pos_fit.Y();posz = pos_fit.Z();
      	             // reconstruction performance figure of merit: position and direction
                     posFOM = ev.GetFitResult(fitName).GetFOM("PositionLogL");
                     posFOM2 = ev.GetFitResult(fitName).GetFOM("PositionSelectedNHit");
                     scaleLogL = posFOM/posFOM2;
                     dirFOM = ev.GetFitResult(fitName).GetFOM("DirectionLogL");
                     dirFOM2 = ev.GetFitResult(fitName).GetFOM("DirectionSelectedNHit");
                     dirscaleLogL = dirFOM/dirFOM2;

                     radiusOfEvent = sqrt(posx*posx+posy*posy+(posz-108)*(posz-108));
      	             posTheta = pos_fit.CosTheta(); posPhi = pos_fit.Phi();
                     posRad = radiusOfEvent;
                     radius->Fill(radiusOfEvent);
                     dirx = u_fit.X();diry = u_fit.Y();dirz = u_fit.Z();		
                     if( !ev.GetFitResult("waterFitter").GetVertex(0).ContainsEnergy() ) continue;
      	             RAT::DS::FitResult fResult = ev.GetFitResult("waterFitter");
      	             if(!fResult.GetValid()) continue;
      	             //cout<< " Fit is Valid "<<endl ;
      	             energyOrigin = fResult.GetVertex(0).GetEnergy();
                     //Note !!! correct energy for rat-6.17.x
		     energy = eCorr.CorrectEnergyRSP(energyOrigin,2);
		     // reconstruction performance figure of metir: energy
                     Gtest =  fResult.GetFOM("EnergyGtest");
                     Utest =  fResult.GetFOM("EnergyUtest");
                     //!!!!!!!!!!!!!!!
                     double countPMT_ITR = 0;
                     vector<TVector3> pmtDir;
                     const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
                     double countSelectPMT = 0;

                     if( !ev.ClassifierResultExists("ITR:waterFitter") ) continue;
                     if ( !ev.GetClassifierResult( "ITR:waterFitter" ).GetValid() ) continue;
                     itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );

                     if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
                     if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
                     beta14 = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );
                     if(nhits>5 && itr>0.55 && beta14<0.95 && beta14>-0.12 && energy>3.5 && posRad<5850 && scaleLogL>10)
                     {
                        hBeta14->Fill(beta14);
     			// LETA cuts (a small effect): the LETA proximity+direction cut is applied. This
                        // cut keeps all events that reconstruct more than 70 cm away from the source, and events that reconstruct
                        // within 70 cm are kept if the event direction is within 45 degree of the vector from the source to the event vertex.
      		       //if( posRad>700 ) {
      		       // hBeta14->Fill(beta14);}
                       // else {
                       //   if( (srcPos - pos_fit).Unit()*u_fit>sqrt(2)/2 ) //costheta >sqrt2/2
                       //      hBeta14->Fill(beta14);
		      //}
      	             }
      	             if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
                     if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
                     thij = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "thetaij" );

                     // Fill the nHit histogram if the event was inside the AV
                     if (radiusOfEventCor >= 0. && posRad<= 6000.) {
                        nHitAV->Fill(nHitOfEvent);
                     }
                     
                    //------------------------------------  
      	            cosThetaToSun = -2.;  
      	            // Now look the direction to the Sun
                    RAT::DS::UniversalTime eventTime = ev.GetUniversalTime();
      	            day = eventTime.GetDays() ;
      	            sec = eventTime.GetSeconds() ;
      	            nsecs = eventTime.GetNanoSeconds() ;
      	            TVector3 directionSun = RAT::SunDirection(eventTime.GetDays(),eventTime.GetSeconds(),eventTime.GetNanoSeconds());
      	            sunDirX = directionSun.X();sunDirY = directionSun.Y();sunDirZ = directionSun.Z();
      	            TVector3 eventDir; eventDir.SetXYZ(dirx,diry,dirz);
      	            cosThetaToSun = eventDir*(-1.0*directionSun);
                    tree->Fill();
             }//fit Valid Direction
            }//fit Valid Position
           }//try catch
           catch(exception& e)
           {std::cout<<e.what()<<" problems in pos fit"<<std::endl;}
            }//fit exist
          } // trig cuts
      	   hTrig->Fill(TriggerType);
         } //data cleaning
       }//loop triggered events
  }//triggered events
   }//loop ds events 
  // Give everything an axis label, as should always be
  radius->SetXTitle( "Radius [mm]" );radius->SetYTitle( "Counts" );
  nHit->SetXTitle( "nHit" );nHit->SetYTitle( "Counts" );  
  nHitAV->SetXTitle( "nHit" );nHitAV->SetYTitle( "Counts" );  
  cout<<countValid<<endl;
  // Write the histograms to fileName
  TString newname = "ExtractValues_";
  TString fileName(filename);
  TString runID(fileName(fileName.Index("r00")+5,6));
//  TString processname = newname+runID+"_test.root";
  TString processname = newname+filename;

  TFile *file=new TFile(processname,"RECREATE");
  file->cd();
  tree->Write();
  hITR->Write();hBeta14->Write();hITR_vs_Beta14->Write();hTrig->Write();hThetaij->Write();
  radius->Write();nHit->Write();nHitAV->Write();
  // candidateEventInfo->Write();
}
