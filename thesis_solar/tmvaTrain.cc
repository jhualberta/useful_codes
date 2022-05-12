#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TString.h>
#include <TMVA/Factory.h>
// 5<E<15 !!!!!!!!!!!!!!
// half dataset 70% sig: 222051 70% bkg: 4451, Only test on Nue for half
//
// whole dataset:      70% sig: 459721 bkg: 9325
// whole dataset Numu: 70% sig: 408167 bkg: 9325
//
int main() {  

// Create ouput file, factory object and open the input file

  TFile* outputFile = TFile::Open( "TMVA.root", "RECREATE" );
  TMVA::Factory* factory = new TMVA::Factory("tmvaTest", outputFile, "");
// ---- 5<E<15
//  TFile* trainingFile = new TFile("../../trainRandomDatasetNueAllBkg_FV5p5_E5to15_whole.root");//trainAllDatasetNueAllBkg_FV5p5_E5to15_whole.root");
  TFile* trainingFile = new TFile("../../trainRandomDatasetNumuAllBkg_FV5p5_E5to15_whole.root");

// ---- 4<E<15  
//  TFile* trainingFile = new TFile("../../trainRandomDatasetNueAllBkg_FV5p5_E4to15_whole.root");
//  TFile* trainingFile = new TFile("../../trainRandomDatasetNumuAllBkg_FV5p5_E4to15_whole.root");
// ---- 4<E<5
//  TFile* trainingFile = new TFile("../../trainRandomDatasetNueAllBkg_FV5p5_E4to5_whole.root");
//  TFile* trainingFile = new TFile("../../trainRandomDatasetNumuAllBkg_FV5p5_E4to5_whole.root");
  
// --- half data test ---
//  TFile* trainingFile = new TFile("../../trainRandomDatasetNueAllBkg_FV5p5_E5to15_Half.root");//trainAllDatasetNueAllBkg_FV5p5_E5to15_whole.root");
  
// get the TTree objects from the input files

  TTree* sigTrain = (TTree*)trainingFile->Get("sig");
  TTree* bkgTrain = (TTree*)trainingFile->Get("bkg");
  int nSigTrain = sigTrain->GetEntries();
  int nBkgTrain = bkgTrain->GetEntries();
//
//  TTree* sigTest = (TTree*)testFile->Get("sig");
//  TTree* bkgTest = (TTree*)testFile->Get("bkg");
//  int nSigTest = sigTest->GetEntries();
//  int nBkgTest = bkgTest->GetEntries();

// global event weights (see below for setting event-wise weights)

  double sigWeight = 1.0;
  double bkgWeight = 1.0;

  factory->AddSignalTree ( sigTrain, 1.0, "Training and Testing" );
  factory->AddBackgroundTree ( bkgTrain, 1.0,"Training and Testing" );

// half--test
//  factory->PrepareTrainingAndTestTree("","nTrain_Signal=222051:nTrain_Background=4451:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );

// whole dataset
// -- Nue
//  factory->PrepareTrainingAndTestTree("","nTrain_Signal=459721:nTrain_Background=9325:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );
// -- Numu
  factory->PrepareTrainingAndTestTree("","nTrain_Signal=408167:nTrain_Background=9325:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );
  
// for 4<E<15
// -- Nue
//  factory->PrepareTrainingAndTestTree("","nTrain_Signal=628463:nTrain_Background=241177:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );
// -- Numu
//  factory->PrepareTrainingAndTestTree("","nTrain_Signal=565025:nTrain_Background=241177:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );

// for 4<E<5
// // -- Nue   70% sig: 168742 70% bkg: 231853
//  factory->PrepareTrainingAndTestTree("","nTrain_Signal=168742:nTrain_Background=231853:nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V" );

//  factory->AddSignalTree(sigTrain, sigWeight, TMVA::Types::kTraining);
//  factory->AddBackgroundTree(bkgTrain, bkgWeight, TMVA::Types::kTraining);
//  factory->AddSignalTree(sigTest, sigWeight, TMVA::Types::kTesting);
//  factory->AddBackgroundTree(bkgTest, bkgWeight, TMVA::Types::kTesting);
   
// Define the input variables that shall be used for the MVA training
// (the variables used in the expression must exist in the original TTree).
/* 9/8/6 params are used !!!! */
  factory->AddVariable("energy", 'F');
  factory->AddVariable("itr", 'F');
  factory->AddVariable("beta14", 'F');
  factory->AddVariable("Gtest", 'F');
  factory->AddVariable("Utest", 'F');
  factory->AddVariable("zfactor", 'F');
  factory->AddVariable("scaleLogL", 'F');
  factory->AddVariable("udotR", 'F');
//  factory->AddVariable("klDiv", 'F');

// Book MVA methods (see TMVA manual).  E.g. for Multilayer perceptron
// (MLP) below, HiddenLayers=3 means 1 hidden layer with 3 nodes.

  factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher");  
  factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=400:MinNodeSize=10%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
  factory->BookMethod( TMVA::Types::kMLP, "MLP","H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=200:HiddenLayers=N+6:TestRate=5");
//		  H:!V:NeuronType=linear:VarTransform=N:NCycles=200:HiddenLayers=N+2:TestRate=2:EstimatorType=CE:BPMode=batch");//TrainingMethod=BFGS");
//  factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN","!H:!V:NCycles=200:HiddenLayers=4+1,4:LearningMethod=BFGS:ValidationFraction=0.1" );
//  factory->BookMethod( TMVA::Types::kBDT, "BDTB","!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
//  factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None" );

  // Train, test and evaluate all methods

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

// Save the output and finish up

  outputFile->Close();
  std::cout << "==> wrote root file TMVA.root" << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl; 

  delete factory;
  return 0;

}
