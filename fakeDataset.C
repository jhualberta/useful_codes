#include "TRandom3.h"
#include "TRandom1.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <iostream>
using namespace std;
void fakeDataset()
{
  const int Nseed = 5000;
  TFile *fsolar = new TFile("GetSolar_promptCut_E5to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to203602_s0_p0.root");
  double energy, cosThetaToSun;
  const double Ethreshold = 5;

  // T0 for signal, T1 for bkg
  TTree* T0 = (TTree*)fsolar->Get("T2");

  int Nfsig = 71;
  int Nfbkg = 38;

  int nFile0 = T0->GetEntries();
  T0->SetBranchAddress("energy", &energy);
  T0->SetBranchAddress("cosThetaToSun", &cosThetaToSun);
  vector<double> venergy_sig;//[NsigSelect];
  vector<double> vcosThetaToSun_sig;//[NsigSelect];
  vector<int> entry_sig; //!!! also save the bkg entry #
  int countSig = 0;
  for(int i = 0; i<nFile0; i++)
  {
    T0->GetEntry(i);
    //!!!!!!!!!!! energy Cut >5 MeV !!!!!!!!!!!
    if(energy>Ethreshold)
    {
    //Venergy_sig[i] = energy;
    //VcosThetaToSun_sig[i] = cosThetaToSun;
      venergy_sig.push_back(energy);
      vcosThetaToSun_sig.push_back(cosThetaToSun);
      entry_sig.push_back(i);
      countSig++; //!! true Entry of bkg
    }
  }

  const int NsigSelect = countSig;
  double Venergy_sig[NsigSelect];
  double VcosThetaToSun_sig[NsigSelect];
  int Ventry_sig[NsigSelect];

  for(int i = 0; i<NsigSelect; i++)
  {
    Venergy_sig[i] = venergy_sig[i];
    VcosThetaToSun_sig[i] = vcosThetaToSun_sig[i];
    Ventry_sig[i] = entry_sig[i];
  }

  TFile *fbkg = new TFile("Merged_MP_GetSolarMC_promptCut_allBkgsTl208andBi214_r200004to203602_s0_p0.root");
  double energy1, cosThetaToSun1;
  TTree* T1 = (TTree*)fbkg->Get("T2");

  int nFile1 = T1->GetEntries();
  T1->SetBranchAddress("energy", &energy1);
  T1->SetBranchAddress("cosThetaToSun", &cosThetaToSun1);
  const int NbkgSelect1 = nFile1;
  vector<double> venergy_bkg;//[NbkgSelect1];
  vector<double> vcosThetaToSun_bkg;//[NbkgSelect1];
  vector<int> entry_bkg; //!!! also save the bkg entry #
  int countBkg = 0;
  for(int i = 0; i<NbkgSelect1; i++)
  {
    T1->GetEntry(i);
    //!!!!!!!!!!! energy Cut >5 MeV !!!!!!!!!!!
    if(energy1>Ethreshold)
    {	    
      //venergy_bkg[i] = energy1;
      //vcosThetaToSun_bkg[i] = cosThetaToSun1;
      venergy_bkg.push_back(energy1);
      vcosThetaToSun_bkg.push_back(cosThetaToSun1);
      entry_bkg.push_back(i);
      countBkg++; //!! true Entry of bkg
    }
  }
  
  const int NbkgSelect = countBkg;
  double Venergy_bkg[NbkgSelect];
  double VcosThetaToSun_bkg[NbkgSelect];
  int Ventry_bkg[NbkgSelect];
  /// cout<<countBkg<<" "<<venergy_bkg.size()<<endl;

  for(int i = 0; i<NbkgSelect; i++)
  {
    Venergy_bkg[i] = venergy_bkg[i];
    VcosThetaToSun_bkg[i] = vcosThetaToSun_bkg[i];
    Ventry_bkg[i] = entry_bkg[i];
  }

  // cout<<"sig = "<<NsigSelect<<" bkg = "<<NbkgSelect<<endl;
  TFile *foutput = new TFile("ensemble_output5000evts.root","RECREATE");  
  int nBkg = 0;
  int nSig = 0;
  int entrySig = 0;
  int entryBkg = 0;
  int countBkgSmall = 0;
  int countBkgBig = 0;
  TRandom *rs1 = new TRandom();
  TRandom *rs2 = new TRandom();

  for(int iseed = 0; iseed<Nseed; iseed++)
  {
     TString tName, tName1, tName2;
     tName.Form("hcosThetaSunVsE_test%d",iseed);
     tName1.Form("hcosTheta_test%d",iseed);
     tName2.Form("hE_test%d",iseed);

     // TNtuple *tn = new TNtuple(tName,"test ","cosThetaToSun:energy");
     TH2F *hcosThetaSunVsE = new TH2F(tName,"",40,-1,1,10,5,15); 
     TH1F *hcosTheta = new TH1F(tName1,"",40,-1,1);
     TH1F *hE = new TH1F(tName2,"",10,5,15);

     TString tNameSig, tName1Sig, tName2Sig, tNameBkg, tName1Bkg, tName2Bkg;
     tNameSig.Form("hcosThetaSunVsE_sig%d",iseed);
     tName1Sig.Form("hcosTheta_sig%d",iseed);
     tName2Sig.Form("hE_sig%d",iseed);

     tNameBkg.Form("hcosThetaSunVsE_bkg%d",iseed);
     tName1Bkg.Form("hcosTheta_bkg%d",iseed);
     tName2Bkg.Form("hE_bkg%d",iseed);

     TH2F *hcosThetaSunVsE_sig = new TH2F(tNameSig,"",40,-1,1,10,5,15);
     TH1F *hcosTheta_sig = new TH1F(tName1Sig,"",40,-1,1);
     TH1F *hE_sig = new TH1F(tName2Sig,"",10,5,15);

     TH2F *hcosThetaSunVsE_bkg = new TH2F(tNameBkg,"",40,-1,1,10,5,15);
     TH1F *hcosTheta_bkg = new TH1F(tName1Bkg,"",40,-1,1);
     TH1F *hE_bkg = new TH1F(tName2Bkg,"",10,5,15);

     if (iseed%50 == 0) 
	     cout<<"seed:   "<<iseed<<endl;
     for(int j = 0;j<Nfsig;j++)
     {
       entrySig = rs1->Uniform(0,NsigSelect);
       double cosVal = VcosThetaToSun_sig[entrySig];
       double energyVal = Venergy_sig[entrySig]; 
       // cout<<" sig "<<cosVal<<" "<<energyVal<<endl;
       // T0->GetEntry(entrySig);
       // tn->Fill(cosThetaToSun, energy);
       hcosThetaSunVsE->Fill(cosVal, energyVal);
       hcosTheta->Fill(cosVal);
       hE->Fill(energyVal);
       hcosThetaSunVsE_sig->Fill(cosVal, energyVal);
       hcosTheta_sig->Fill(cosVal);
       hE_sig->Fill(energyVal);
     } 
     for(int j = 0;j<Nfbkg;j++)
     {
       entryBkg = rs2->Uniform(0,NbkgSelect);
       // double cosVal1 = rs2->Uniform(-1,1);
       double cosVal1 = VcosThetaToSun_bkg[entryBkg];
       double energyVal1 = Venergy_bkg[entryBkg];
       // cout<<"bkg "<<entryBkg<<" "<<cosVal1<<" "<<energyVal1<<endl;
       // tn->Fill(cosThetaToSun1, energy1);
       //cout<<" bkg "<<cosVal1<<" "<<energyVal1<<endl;
       //if(cosVal1<0) countBkgSmall++;
       //if(cosVal1>0) countBkgBig++;
       hcosThetaSunVsE->Fill(cosVal1, energyVal1);
       hcosTheta->Fill(cosVal1);
       hE->Fill(energyVal1);
       hcosThetaSunVsE_bkg->Fill(cosVal1, energyVal1);
       hcosTheta_bkg->Fill(cosVal1);
       hE_bkg->Fill(energyVal1);
     }
     foutput->cd();
     // tn->Write();
     hcosThetaSunVsE->Write();hcosTheta->Write();hE->Write();
     hcosThetaSunVsE_bkg->Write();hcosTheta_bkg->Write();hE_bkg->Write();
     hcosThetaSunVsE_sig->Write();hcosTheta_sig->Write();hE_sig->Write();

     // foutput->Close();
     delete hcosThetaSunVsE;
     delete hcosTheta;
     delete hE; 

     delete hcosThetaSunVsE_bkg;delete hcosTheta_bkg;delete hE_bkg;
     delete hcosThetaSunVsE_sig;delete hcosTheta_sig;delete hE_sig;

  }

//  cout<<"small "<<countBkgSmall<<" big "<<countBkgBig<<endl; 


}
