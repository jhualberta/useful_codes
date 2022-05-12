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
//!!!!! ONLY for MC
int main() {

  bool checkMC = 1;
  TString dir = "../../";
  TString fileName = "GetSolarMC_containMC_E5to15_Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to207718.root";
  // TString fileName = "GetSolarMC_containMC_E5to15_Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NumuRun_r200004to207718.root";
  //----------------------------------------------------------------------------------------------------------------------------------------------------
  TFile* inputFile = new TFile(dir+fileName);
  // Set up an output file and book some histograms
  TFile* histFile = new TFile("TMVA_"+fileName, "RECREATE");

  TH1D* hBDTsig = new TH1D("hBDTsig", "BDT, signal", 100, -10.0, 10.0);
  TH1D* hBDTbkg = new TH1D("hBDTbkg", "BDT, background", 100, -10.0, 10.0);
  TH1D* hMLPsig = new TH1D("hMLPsig", "MLP, signal", 100, -0.5, 1.5);
  TH1D* hMLPbkg = new TH1D("hMLPbkg", "MLP, background", 100, -0.5, 1.5);

  TH1D* hBDTsigCosThetaToSun = new TH1D("hBDTsigCosThetaToSun", "BDT, signal", 40, -1.0, 1.0);
  TH1D* hBDTbkgCosThetaToSun = new TH1D("hBDTbkgCosThetaToSun", "BDT, background", 40, -1.0, 1.0);
  TH1D* hMLPsigCosThetaToSun = new TH1D("hMLPsigCosThetaToSun", "MLP, signal", 40, -1.0, 1.0);
  TH1D* hMLPbkgCosThetaToSun = new TH1D("hMLPbkgCosThetaToSun", "MLP, background", 40, -1.0, 1.0);

  TH2D* hBDTsigCosThetaToSunVsE = new TH2D("hBDTsigCosThetaToSunVsE","BDT, signal, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hBDTbkgCosThetaToSunVsE = new TH2D("hBDTbkgCosThetaToSunVsE","BDT, background, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hMLPsigCosThetaToSunVsE = new TH2D("hMLPsigCosThetaToSunVsE","MLP, signal, cos vs E", 40, -1.0, 1.0, 100, 5, 15);
  TH2D* hMLPbkgCosThetaToSunVsE = new TH2D("hMLPbkgCosThetaToSunVsE","MLP, background, cos vs E", 40, -1.0, 1.0, 100, 5, 15);

  TH2D* hBDTsigPosRhoZ = new TH2D("hBDTsigPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPsigPosRhoZ = new TH2D("hMLPsigPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hBDTbkgPosRhoZ = new TH2D("hBDTbkgPosRhoZ","",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPbkgPosRhoZ = new TH2D("hMLPbkgPosRhoZ","",1000,0,6000,2000,-6000,6000);

  TH1D* hBDTsigEnergy = new TH1D("hBDTsigEnergy","BDT, signal, Energy", 10, 5, 15);
  TH1D* hBDTbkgEnergy = new TH1D("hBDTbkgEnergy","BDT, background, Energy", 10, 5, 15);
  TH1D* hMLPsigEnergy = new TH1D("hMLPsigEnergy","MLP, signal, Energy", 10, 5, 15);
  TH1D* hMLPbkgEnergy = new TH1D("hMLPbkgEnergy","MLP, background, Energy", 10, 5, 15);

  TNtuple *ntupleBDTsig = new TNtuple("ntupleBDTsig","BDT signal","cosThetaToSun:energy:energymc:posx:posy:posz:Utest:Gtest:scaleLogL:beta14:udotR:klDiv");
  TNtuple *ntupleBDTbkg = new TNtuple("ntupleBDTbkg","BDT bkg","cosThetaToSun:energy:energymc:posx:posy:posz:Utest:Gtest:scaleLogL:beta14:udotR:klDiv");
  TNtuple *ntupleMLPsig = new TNtuple("ntupleMLPsig","MLP signal","cosThetaToSun:energy:energymc:posx:posy:posz:Utest:Gtest:scaleLogL:beta14:udotR:klDiv");
  TNtuple *ntupleMLPbkg = new TNtuple("ntupleMLPbkg","MLP bkg","cosThetaToSun:energy:energymc:posx:posy:posz:Utest:Gtest:scaleLogL:beta14:udotR:klDiv");

  /// for MC
  TH1D* hBDTsigEnergyMC = new TH1D("hBDTsigEnergyMC","BDT, signal, MC Energy", 10, 5, 15);
  TH1D* hBDTbkgEnergyMC = new TH1D("hBDTbkgEnergyMC","BDT, background, MC Energy", 10, 5, 15);
  TH1D* hMLPsigEnergyMC = new TH1D("hMLPsigEnergyMC","MLP, signal, MC Energy", 10, 5, 15);
  TH1D* hMLPbkgEnergyMC = new TH1D("hMLPbkgEnergyMC","MLP, background, MC Energy", 10, 5, 15);
  TH2D* hBDTsigPosRhoZMC = new TH2D("hBDTsigPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPsigPosRhoZMC = new TH2D("hMLPsigPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hBDTbkgPosRhoZMC = new TH2D("hBDTbkgPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);
  TH2D* hMLPbkgPosRhoZMC = new TH2D("hMLPbkgPosRhoZMC","MC position",1000,0,6000,2000,-6000,6000);

  TH2D* hBDTsigCosThetaToSunVsEmc = new TH2D("hBDTsigCosThetaToSunVsEmc","BDT, signal, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);
  TH2D* hBDTbkgCosThetaToSunVsEmc = new TH2D("hBDTbkgCosThetaToSunVsEmc","BDT, background, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);
  TH2D* hMLPsigCosThetaToSunVsEmc = new TH2D("hMLPsigCosThetaToSunVsEmc","MLP, signal, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);
  TH2D* hMLPbkgCosThetaToSunVsEmc = new TH2D("hMLPbkgCosThetaToSunVsEmc","MLP, background, cos vs Emc", 40, -1.0, 1.0, 10, 5, 15);

  // store true energy for oscillation
  TNtuple *ntupleBDTsigMC = new TNtuple("ntupleBDTsigMC","BDT signal, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleBDTbkgMC = new TNtuple("ntupleBDTbkgMC","BDT bkg, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleMLPsigMC = new TNtuple("ntupleMLPsigMC","MLP signal, Emc","cosThetaToSun:energymc");
  TNtuple *ntupleMLPbkgMC = new TNtuple("ntupleMLPbkgMC","MLP bkg, Emc","cosThetaToSun:energymc");

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
  reader->AddVariable("klDiv", &klDiv);
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

//  std::string dirWeight    = "../train/testFinal_E5to15_Nue_whole/";
//  std::string dirWeight    = "../train/testFinal_E5to15_Numu_whole/";
  std::string dirWeight = "../train/weights_9pars_E5to15_whole/";
//  std::string dirWeight  = "../train/weightsNumu_9pars_E5to15_whole/";

  std::string prefix = "tmvaTest";
  reader->BookMVA("BDT", dirWeight + prefix + "_BDT.weights.xml");
  reader->BookMVA("MLP", dirWeight + prefix + "_MLP.weights.xml");

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
  int nSig, nbkg;
  const double tCutBDT = -0.08;//0.0;
  const double tCutMLP = 0.65;

  /* convert float to double !!!*/
  double posRad_d, energy_d, itr_d, beta14_d, Gtest_d, Utest_d, zfactor_d, scaleLogL_d, udotR_d, klDiv_d, cosThetaToSun_d;
  double posx_d, posy_d, posz_d, dirx_d, diry_d, dirz_d;
  double posxmc_d, posymc_d, poszmc_d, energymc_d;

  for (int i=0; i<treeVec.size(); i++){

    // treeVec[i]->Print();
    // treeVec[i]->SetBranchAddress("posRad", &posRad_d);
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
    // treeVec[i]->SetBranchAddress("dirx", &dirx_d);
    // treeVec[i]->SetBranchAddress("diry", &diry_d);
    // treeVec[i]->SetBranchAddress("dirz", &dirz_d);
 
    treeVec[i]->SetBranchAddress("energymc", &energymc_d);//for MC!!!
    // treeVec[i]->SetBranchAddress("posxmc", &posxmc_d);//for MC!!!
    // treeVec[i]->SetBranchAddress("posymc", &posymc_d);//for MC!!!
    // treeVec[i]->SetBranchAddress("poszmc", &poszmc_d);//for MC!!!
    
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
      // !!!!!
      energymc = energymc_d;
      //if(checkMC) { // for solarMC nue numu
      // energy = energymc_d;
      //}//
      double tBDT = reader->EvaluateMVA("BDT");
      double tMLP = reader->EvaluateMVA("MLP");
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
           ntupleBDTsig->Fill(cosThetaToSun, energy, energymc, itr, beta14, posx, posy, posz, Utest, Gtest,scaleLogL,beta14_d,udotR,klDiv);

	   hBDTsigEnergyMC->Fill(energymc);
           hBDTsigPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
	   hBDTsigCosThetaToSunVsEmc->Fill(cosThetaToSun, energymc);
       }
       if ( tMLP >= tCutMLP ) { 
           hMLPsig->Fill(tMLP); nSigAccMLP++; hMLPsigCosThetaToSun->Fill(cosThetaToSun);
           hMLPsigCosThetaToSunVsE->Fill(cosThetaToSun, energy);
           hMLPsigEnergy->Fill(energy);
           hMLPsigPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
           ntupleMLPsig->Fill(cosThetaToSun, energy, energymc, itr, beta14, posx, posy, posz, Utest, Gtest,scaleLogL,beta14_d,udotR,klDiv);

	   hMLPsigEnergyMC->Fill(energymc);
           hMLPsigPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
           hMLPsigCosThetaToSunVsEmc->Fill(cosThetaToSun, energymc);
       }
       if ( tBDT < tCutBDT ) {
          hBDTbkg->Fill(tBDT); nbkgAccBDT++; hBDTbkgCosThetaToSun->Fill(cosThetaToSun);
          hBDTbkgCosThetaToSunVsE->Fill(cosThetaToSun, energy);
          hBDTbkgEnergy->Fill(energy);
          hBDTbkgPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
          ntupleBDTbkg->Fill(cosThetaToSun, energy, energymc, itr, beta14, posx, posy, posz, Utest, Gtest,scaleLogL,udotR,klDiv);

          hBDTbkgEnergyMC->Fill(energymc);
          hBDTbkgPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
          hBDTbkgCosThetaToSunVsEmc->Fill(cosThetaToSun, energymc);
       }
       if ( tMLP < tCutMLP ) { 
          hMLPbkg->Fill(tMLP); nbkgAccMLP++; hMLPbkgCosThetaToSun->Fill(cosThetaToSun);
          hMLPbkgCosThetaToSunVsE->Fill(cosThetaToSun, energy);
          hMLPbkgEnergy->Fill(energy);
          hMLPbkgPosRhoZ->Fill(sqrt(posx*posx+posy*posy),posz-108);
          ntupleMLPbkg->Fill(cosThetaToSun, energy, energymc, itr, beta14, posx, posy, posz, Utest, Gtest,scaleLogL,udotR,klDiv);

          hMLPbkgEnergyMC->Fill(energymc);
          hMLPbkgPosRhoZMC->Fill(sqrt(posxmc*posxmc+posymc*posymc),poszmc-108);
          hMLPbkgCosThetaToSunVsEmc->Fill(cosThetaToSun, energymc);
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

  std::cout << "BDT signal efficiency     = " << epsSigBDT << std::endl;
//  std::cout << "BDT background efficiency = " << epsbkgBDT << std::endl;
  std::cout << "MLP signal efficiency        = " << epsSigMLP << std::endl;
//  std::cout << "MLP background efficiency    = " << epsbkgMLP << std::endl;

  histFile->Write();
  histFile->Close();
  return 0;

}
