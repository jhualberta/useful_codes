#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include "Event.h"
#include "TNtuple.h"
#include <TMVA/Reader.h>
using namespace std;
int main() {
// ONly for DATA!!!!! MC needs posMC, energyMC info, use files_solar
  bool checkMC = 0;
  TString dir = "../../";

// ------ data ------
  TString fileName = "GetSolar_promptCut_E5to15_fullData_Merged_MP_p1mask_ExtractTres_WaterMP6176general_nhit15_Analysis_r200004to207718_p004.root";

  TFile* inputFile = new TFile(dir+fileName);
  //// Set up an output file and book some histograms
  TFile* histFile = new TFile("tmva8parsFinal_E5to15_fullData_Merged_MP_p1mask_whole.root", "RECREATE");

  TH1D* hBDTsig = new TH1D("hBDTsig", "BDT, signal", 100, -10.0, 10.0);
  TH1D* hBDTbkg = new TH1D("hBDTbkg", "BDT, background", 100, -10.0, 10.0);
  TH1D* hMLPsig = new TH1D("hMLPsig", "MLP, signal", 100, -0.5, 1.5);
  TH1D* hMLPbkg = new TH1D("hMLPbkg", "MLP, background", 100, -0.5, 1.5);
  TH1D* hFISHERsig = new TH1D("hFISHERsig", "FISHER, signal", 100, -0.5, 1.5);
  TH1D* hFISHERbkg = new TH1D("hFISHERbkg", "FISHER, background", 100, -0.5, 1.5);

  TH1D* hBDTsigCosThetaToSun = new TH1D("hBDTsigCosThetaToSun", "BDT, signal", 40, -1.0, 1.0);
  TH1D* hBDTbkgCosThetaToSun = new TH1D("hBDTbkgCosThetaToSun", "BDT, background", 40, -1.0, 1.0);
  TH1D* hMLPsigCosThetaToSun = new TH1D("hMLPsigCosThetaToSun", "MLP, signal", 40, -1.0, 1.0);
  TH1D* hMLPbkgCosThetaToSun = new TH1D("hMLPbkgCosThetaToSun", "MLP, background", 40, -1.0, 1.0);
  TH1D* hFISHERsigCosThetaToSun = new TH1D("hFISHERsigCosThetaToSun", "FISHER, signal", 40, -1.0, 1.0);
  TH1D* hFISHERbkgCosThetaToSun = new TH1D("hFISHERbkgCosThetaToSun", "FISHER, background", 40, -1.0, 1.0);

  TH2D* hBDTsigCosThetaToSunVsE = new TH2D("hBDTsigCosThetaToSunVsE","BDT, signal, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hBDTbkgCosThetaToSunVsE = new TH2D("hBDTbkgCosThetaToSunVsE","BDT, background, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hMLPsigCosThetaToSunVsE = new TH2D("hMLPsigCosThetaToSunVsE","MLP, signal, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hMLPbkgCosThetaToSunVsE = new TH2D("hMLPbkgCosThetaToSunVsE","MLP, background, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hFISHERsigCosThetaToSunVsE = new TH2D("hFISHERsigCosThetaToSunVsE","FISHER, signal, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hFISHERbkgCosThetaToSunVsE = new TH2D("hFISHERbkgCosThetaToSunVsE","FISHER, background, cos vs E", 40, -1.0, 1.0, 100, 5, 15);

  TH2D* hBDTsigPosRhoZ = new TH2D("hBDTsigPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hBDTbkgPosRhoZ = new TH2D("hBDTbkgPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPsigPosRhoZ = new TH2D("hMLPsigPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPbkgPosRhoZ = new TH2D("hMLPbkgPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hFISHERsigPosRhoZ = new TH2D("hFISHERsigPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hFISHERbkgPosRhoZ = new TH2D("hFISHERbkgPosRhoZ","",1000,0,6000,2000,-6000,6000);

  TH1D* hBDTsigEnergy = new TH1D("hBDTsigEnergy","BDT, signal, Energy", 10, 5, 15);
  TH1D* hBDTbkgEnergy = new TH1D("hBDTbkgEnergy","BDT, background, Energy", 10, 5, 15);
  TH1D* hMLPsigEnergy = new TH1D("hMLPsigEnergy","MLP, signal, Energy", 10, 5, 15);
  TH1D* hMLPbkgEnergy = new TH1D("hMLPbkgEnergy","MLP, background, Energy", 10, 5, 15);
  TH1D* hFISHERsigEnergy = new TH1D("hFISHERsigEnergy","FISHER, signal, Energy", 10, 5, 15);
  TH1D* hFISHERbkgEnergy = new TH1D("hFISHERbkgEnergy","FISHER, background, Energy", 10, 5, 15);

  // cosThetaToSun, energy, energymc, klDiv, posx, posy, posz, beta14, Gtest, Utest, scaleLogL,zfactor,udotR
  TNtuple *ntupleBDTsig = new TNtuple("ntupleBDTsig","BDT signal","cosThetaToSun:energy:itr:beta14:posx:posy:posz:Utest:Gtest:scaleLogL:zfactor:udotR:klDiv");
  TNtuple *ntupleBDTbkg = new TNtuple("ntupleBDTbkg","BDT bkg","cosThetaToSun:energy:itr:beta14:posx:posy:posz:Utest:Gtest:scaleLogL:zfactor:udotR:klDiv");
  TNtuple *ntupleMLPsig = new TNtuple("ntupleMLPsig","MLP signal","cosThetaToSun:energy:itr:beta14:posx:posy:posz:Utest:Gtest:scaleLogL:zfactor:udotR:klDiv");
  TNtuple *ntupleMLPbkg = new TNtuple("ntupleMLPbkg","MLP bkg","cosThetaToSun:energy:itr:beta14:posx:posy:posz:Utest:Gtest:scaleLogL:zfactor:udotR:klDiv");
  TNtuple *ntupleFISHERsig = new TNtuple("ntupleFISHERsig","FISHER signal","cosThetaToSun:energy:itr:beta14:posx:posy:posz:Utest:Gtest:scaleLogL:zfactor:udotR:klDiv");
  TNtuple *ntupleFISHERbkg = new TNtuple("ntupleFISHERbkg","FISHER bkg","cosThetaToSun:energy:itr:beta14:posx:posy:posz:Utest:Gtest:scaleLogL:zfactor:udotR:klDiv");

  /// for MC
  TH1D* hBDTsigEnergyMC = new TH1D("hBDTsigEnergyMC","BDT, signal, MC Energy", 10, 5, 15);
  TH1D* hBDTbkgEnergyMC = new TH1D("hBDTbkgEnergyMC","BDT, background, MC Energy", 10, 5, 15);
  TH2D* hBDTsigPosRhoZMC = new TH2D("hBDTsigPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hBDTbkgPosRhoZMC = new TH2D("hBDTbkgPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hBDTsigCosThetaToSunVsEmc = new TH2D("hBDTsigCosThetaToSunVsEmc","BDT, signal, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);
  TH2D* hBDTbkgCosThetaToSunVsEmc = new TH2D("hBDTbkgCosThetaToSunVsEmc","BDT, background, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);

  TH1D* hMLPsigEnergyMC = new TH1D("hMLPsigEnergyMC","MLP, signal, MC Energy", 10, 5, 15);
  TH1D* hMLPbkgEnergyMC = new TH1D("hMLPbkgEnergyMC","MLP, background, MC Energy", 10, 5, 15);
  TH2D* hMLPsigPosRhoZMC = new TH2D("hMLPsigPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPbkgPosRhoZMC = new TH2D("hMLPbkgPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPsigCosThetaToSunVsEmc = new TH2D("hMLPsigCosThetaToSunVsEmc","MLP, signal, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);
  TH2D* hMLPbkgCosThetaToSunVsEmc = new TH2D("hMLPbkgCosThetaToSunVsEmc","MLP, background, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);

  TH1D* hFISHERsigEnergyMC = new TH1D("hFISHERsigEnergyMC","FISHER, signal, MC Energy", 10, 5, 15);
  TH1D* hFISHERbkgEnergyMC = new TH1D("hFISHERbkgEnergyMC","FISHER, background, MC Energy", 10, 5, 15);
  TH2D* hFISHERsigPosRhoZMC = new TH2D("hFISHERsigPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hFISHERbkgPosRhoZMC = new TH2D("hFISHERbkgPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hFISHERsigCosThetaToSunVsEmc = new TH2D("hFISHERsigCosThetaToSunVsEmc","FISHER, signal, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);
  TH2D* hFISHERbkgCosThetaToSunVsEmc = new TH2D("hFISHERbkgCosThetaToSunVsEmc","FISHER, background, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);

  // store true energy for oscillation
  TNtuple *ntupleBDTsigMC = new TNtuple("ntupleBDTsigMC","BDT signal, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleBDTbkgMC = new TNtuple("ntupleBDTbkgMC","BDT bkg, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleMLPsigMC = new TNtuple("ntupleMLPsigMC","MLP signal, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleMLPbkgMC = new TNtuple("ntupleMLPbkgMC","MLP bkg, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleFISHERsigMC = new TNtuple("ntupleFISHERsigMC","FISHER signal, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleFISHERbkgMC = new TNtuple("ntupleFISHERbkgMC","FISHER bkg, Emc","cosThetaToSun:energymc");

// Set up the TMVA Reader object.
// The names in AddVariable must be same as in the input (weight) file.
  TTree *T2 = (TTree*)inputFile->Get("T2");
  TMVA::Reader* reader = new TMVA::Reader();
  float posx, posy, posz, posRad, energy, itr, beta14, Gtest, Utest, zfactor, scaleLogL, udotR, klDiv;// TMVA needs float, not double
  float cosThetaToSun;
  // MC info
  float posxmc, posymc, poszmc, energymc;

  // exactly in the same order of the training
  reader->AddVariable("energy",& energy);
  reader->AddVariable("itr",&itr); 
  reader->AddVariable("beta14", &beta14);
  reader->AddVariable("Gtest", &Gtest);
  reader->AddVariable("Utest", &Utest);
  reader->AddVariable("zfactor", &zfactor);
  reader->AddVariable("scaleLogL", &scaleLogL);
  reader->AddVariable("udotR", &udotR);
//  reader->AddVariable("klDiv", &klDiv);

//===================================================
//  reader->AddVariable("posxmc", &posxmc);
//  reader->AddVariable("posymc", &posymc);
//  reader->AddVariable("poszmc", &poszmc);
//  reader->AddVariable("energymc", &energymc);
//  reader->AddVariable("posRad",& posRad);

////  T2->SetBranchAddress("posxmc", &posxmc);
////  T2->SetBranchAddress("posymc", &posymc);
////  T2->SetBranchAddress("poszmc", &poszmc);
////  T2->SetBranchAddress("energymc", &energymc);
////  T2->SetBranchAddress("dirxmc", &dirxmc);
////  T2->SetBranchAddress("dirymc", &dirymc);
////  T2->SetBranchAddress("dirzmc", &dirzmc);

//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_external/train/weights_E5to15_half/";
//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_external/train/weightsNumu_E5to15_half/";

//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_allbkg/train/weights_9pars_E5to15_whole/";
//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_allbkg/train/weightsNumu_9pars_E5to15_whole/";

  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_allbkg/train/weights_8pars_E5to15_whole/";

//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_allbkg/train/weights_7pars_E5to15_whole/";
//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_allbkg/train/weightsNumu_9pars_E5to15_whole/";
//  std::string dirWeight  = "../../files_MCsolar/tmvaCowan_external/train/weights_E5to15_external/";

  std::string prefix = "tmvaTest"; // weights!!
  reader->BookMVA("BDT", dirWeight + prefix + "_BDT.weights.xml");
  reader->BookMVA("MLP", dirWeight + prefix + "_MLP.weights.xml");
  reader->BookMVA("Fisher", dirWeight + prefix + "_Fisher.weights.xml");

// Open input file, get the TTrees, put into a vector

  inputFile->ls();
  TTree* sig = dynamic_cast<TTree*>(inputFile->Get("T2"));
//  TTree* bkg = dynamic_cast<TTree*>(inputFile->Get("bkg"));
  std::vector<TTree*> treeVec;
  treeVec.push_back(sig);
//  treeVec.push_back(bkg);

// Loop over TTrees

  int nSigAccBDT = 0;
  int nbkgAccBDT = 0;
  int nSigAccMLP = 0;
  int nbkgAccMLP = 0;
  int nSigAccFISHER = 0;
  int nbkgAccFISHER = 0;

  int nSig, nbkg;
  const double tCutBDT = -0.0578;
  const double tCutMLP = 0.3476;
  const double tCutFISHER = -0.12;

  /* convert float to double !!!*/
  double posRad_d, energy_d, itr_d, beta14_d, Gtest_d, Utest_d, zfactor_d, scaleLogL_d, udotR_d, klDiv_d, cosThetaToSun_d;
  double posx_d, posy_d, posz_d, dirx_d, diry_d, dirz_d;
  double posxmc_d, posymc_d, poszmc_d, energymc_d;

  for (int i=0; i<treeVec.size(); i++){

    // treeVec[i]->Print();
    treeVec[i]->SetBranchAddress("posRad", &posRad_d);
    treeVec[i]->SetBranchAddress("energy", &energy_d);
    treeVec[i]->SetBranchAddress("itr", &itr_d);
    treeVec[i]->SetBranchAddress("beta14", &beta14_d);
    treeVec[i]->SetBranchAddress("Gtest", &Gtest_d);
    treeVec[i]->SetBranchAddress("Utest", &Utest_d);
    treeVec[i]->SetBranchAddress("zfactor", &zfactor_d);
    treeVec[i]->SetBranchAddress("scaleLogL", &scaleLogL_d);
    treeVec[i]->SetBranchAddress("udotR", &udotR_d);
    treeVec[i]->SetBranchAddress("klDiv", &klDiv_d);
    treeVec[i]->SetBranchAddress("cosThetaToSun", &cosThetaToSun_d);
    treeVec[i]->SetBranchAddress("posx", &posx_d);
    treeVec[i]->SetBranchAddress("posy", &posy_d);
    treeVec[i]->SetBranchAddress("posz", &posz_d);
    treeVec[i]->SetBranchAddress("dirx", &dirx_d);
    treeVec[i]->SetBranchAddress("diry", &diry_d);
    treeVec[i]->SetBranchAddress("dirz", &dirz_d);
 
    //if(checkMC) {
    //  treeVec[i]->SetBranchAddress("energymc", &energymc_d);//for MC!!!
    //  treeVec[i]->SetBranchAddress("posxmc", &posxmc_d);//for MC!!!
    //  treeVec[i]->SetBranchAddress("posymc", &posymc_d);//for MC!!!
    //  treeVec[i]->SetBranchAddress("poszmc", &poszmc_d);//for MC!!!
    //}
    int numEntries = treeVec[i]->GetEntries();
    if ( i == 0 ) { nSig = numEntries; }
    // if ( i == 1 ) { nbkg = numEntries; }
    std::cout << "number of entries = " << numEntries << std::endl;

// Loop over events.  The test statistic is identified by the first 
// argument used above in BookMVA (below, e.g., "BDT", "MLP").

    for (int j=0; j<numEntries; j++) {
      treeVec[i]->GetEntry(j);// sets evt
      energy = energy_d, itr = itr_d, beta14 = beta14_d, Gtest = Gtest_d;
      Utest = Utest_d, zfactor = zfactor_d, scaleLogL = scaleLogL_d, udotR = udotR_d, klDiv = klDiv_d;
      posRad = posRad_d;
      posx = posx_d; posy = posy_d; posz = posz_d; cosThetaToSun = cosThetaToSun_d;

      if(checkMC) { // for solarMC nue numu
       energy = energymc_d;
      }//
      double tBDT = reader->EvaluateMVA("BDT");
      double tMLP = reader->EvaluateMVA("MLP");
      double tFISHER = reader->EvaluateMVA("Fisher");
      /*
      //if ( i == 0 ){                       // signal
        hBDTsig->Fill(tBDT);
        hMLPsig->Fill(tMLP);
        if ( tBDT >= tCutBDT ) { nSigAccBDT++; hBDTsigCosThetaToSun->Fill(cosThetaToSun);}
	if ( tMLP >= tCutMLP ) { nSigAccMLP++; hMLPsigCosThetaToSun->Fill(cosThetaToSun);}
      //}
      //else if ( i == 1 ){                  // background
         hBDTbkg->Fill(tBDT);
         hMLPbkg->Fill(tMLP);
         if ( tBDT >= tCutBDT ) { nbkgAccBDT++; hBDTbkgCosThetaToSun->Fill(cosThetaToSun);}
	 if ( tMLP >= tCutMLP ) { nbkgAccMLP++; hMLPbkgCosThetaToSun->Fill(cosThetaToSun);}
      //}
      */
        // cout<<tBDT<<" "<<tMLP<<endl;
       if ( tBDT >= tCutBDT ) {
           hBDTsig->Fill(tBDT); nSigAccBDT++; hBDTsigCosThetaToSun->Fill(cosThetaToSun);
           hBDTsigCosThetaToSunVsE->Fill(cosThetaToSun, energy);
           hBDTsigEnergy->Fill(energy);
           hBDTsigPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
	   ntupleBDTsig->Fill(cosThetaToSun, energy, itr, beta14, posx, posy, posz, Utest, Gtest,scaleLogL,zfactor,udotR,klDiv);
           //if(checkMC) {
           //  hBDTsigEnergyMC->Fill(energymc);
           //  hBDTsigPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
           //}  
       }
       if ( tMLP >= tCutMLP ) { 
           hMLPsig->Fill(tMLP); nSigAccMLP++; hMLPsigCosThetaToSun->Fill(cosThetaToSun);
           hMLPsigCosThetaToSunVsE->Fill(cosThetaToSun, energy);
           hMLPsigEnergy->Fill(energy);
           hMLPsigPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
           ntupleMLPsig->Fill(cosThetaToSun, energy, klDiv, itr, posx, posy, posz, beta14, Gtest, Utest, scaleLogL,zfactor,udotR);
           //if(checkMC) {
           //  hMLPsigEnergyMC->Fill(energymc);
           //  hMLPsigPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
           //}
       }
       if ( tFISHER >= tCutFISHER ) {
           hFISHERsig->Fill(tFISHER); nSigAccFISHER++; hFISHERsigCosThetaToSun->Fill(cosThetaToSun);
           hFISHERsigCosThetaToSunVsE->Fill(cosThetaToSun, energy);
           hFISHERsigEnergy->Fill(energy);
           hFISHERsigPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
           ntupleFISHERsig->Fill(cosThetaToSun, energy, energymc, klDiv, itr, posx, posy, posz, beta14, Gtest, Utest, scaleLogL,zfactor,udotR);

           hFISHERsigEnergyMC->Fill(energymc);
           hFISHERsigPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
           hFISHERsigCosThetaToSunVsEmc->Fill(cosThetaToSun, energymc);
       }

       if ( tBDT < tCutBDT ) {
            hBDTbkg->Fill(tBDT); nbkgAccBDT++; hBDTbkgCosThetaToSun->Fill(cosThetaToSun);
	    hBDTbkgCosThetaToSunVsE->Fill(cosThetaToSun, energy);
	    hBDTbkgEnergy->Fill(energy);
	    hBDTbkgPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
	    ntupleBDTbkg->Fill(cosThetaToSun, energy, klDiv, itr, posx, posy, posz, beta14, Gtest, Utest, scaleLogL,zfactor,udotR);
	    //if(checkMC) {
            //  hBDTbkgEnergyMC->Fill(energymc);
            //  hBDTbkgPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
            //}

	 }

       if ( tMLP < tCutMLP ) { 
            hMLPbkg->Fill(tMLP); nbkgAccMLP++; hMLPbkgCosThetaToSun->Fill(cosThetaToSun);
	    hMLPbkgCosThetaToSunVsE->Fill(cosThetaToSun, energy);
            hMLPbkgEnergy->Fill(energy);
	    hMLPbkgPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
	    ntupleMLPbkg->Fill(cosThetaToSun, energy, klDiv, itr, posx, posy, posz, beta14, Gtest, Utest, scaleLogL,zfactor,udotR);
	    //if(checkMC) {
            //  hMLPbkgEnergyMC->Fill(energymc);
            //  hMLPbkgPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
            //}
	 }
       if ( tFISHER < tCutFISHER ) {
            hFISHERbkg->Fill(tFISHER); nbkgAccFISHER++; hFISHERbkgCosThetaToSun->Fill(cosThetaToSun);
            hFISHERbkgCosThetaToSunVsE->Fill(cosThetaToSun, energy);
            hFISHERbkgEnergy->Fill(energy);
            hFISHERbkgPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
            ntupleFISHERbkg->Fill(cosThetaToSun, energy, energymc, klDiv, itr, posx, posy, posz, beta14, Gtest, Utest, scaleLogL,zfactor,udotR);
       }

    }
  }
  double epsSigBDT = static_cast<double>(nSigAccBDT)/
                      static_cast<double>(nSig);
  double epsbkgBDT = static_cast<double>(nbkgAccBDT)/
                      static_cast<double>(nbkg);
  double epsSigMLP = static_cast<double>(nSigAccMLP)/
                      static_cast<double>(nSig);
  double epsbkgMLP = static_cast<double>(nbkgAccMLP)/
                      static_cast<double>(nbkg);

  double epsSigFISHER = static_cast<double>(nSigAccFISHER)/
                      static_cast<double>(nSig);
  double epsbkgFISHER = static_cast<double>(nbkgAccFISHER)/
                      static_cast<double>(nbkg);

  std::cout << "FISHER signal efficiency        = " << epsSigFISHER << std::endl;
//  std::cout << "FISHER background efficiency    = " << epsbkgFISHER << std::endl;
  std::cout << "BDT signal efficiency     = " << epsSigBDT << std::endl;
//  std::cout << "BDT background efficiency = " << epsbkgBDT << std::endl;
  std::cout << "MLP signal efficiency        = " << epsSigMLP << std::endl;
//  std::cout << "MLP background efficiency    = " << epsbkgMLP << std::endl;

  histFile->Write();
  histFile->Close();
  return 0;

}