//Do chi2 Minimize: ((p0*xfit+p1*ufit)-Xmc)^2
//Find p0, p1, p2
//Minimize |Xcor-Xmc|
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DU/PMTInfo.hh>
#include "TVector3.h"
#include <vector>
#include <TMath.h>
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/Minimizer.h"
#include "Math/GSLMinimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <TMinuit.h>
#include "TString.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/GSLMinimizer1D.h"
#include "Math/Functor.h"
using namespace std;

//global data
 
TString specFit = "alberta";
//Spatial Resolution
TString filename = "FitAP_Water_6MeV_10000evt_ISOFill.root";
TFile *f1 = new TFile(filename);
TVector3 u_mc,u_fit,pos_mc,pos_fit;

vector<double> posfitX,posfitY,posfitZ;
vector<double> posmcX, posmcY,posmcZ;
vector<double> ufitX,ufitY, ufitZ;

double gaussianMean = 0;
double gaussianSigma = 100;//mm

void RetrieveData(vector<double> &posfitX,vector<double> &posfitY,vector<double> &posfitZ,
vector<double> &posmcX, vector<double> &posmcY,vector<double> &posmcZ,vector<double> &ufitX, vector<double> &ufitY, vector<double> &ufitZ)
{
 
 RAT::DU::DSReader dsReader("FitAP_Water_6MeV_10000evt_ISOFill.root");
 for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
 {
     //std:cout << " event ID "<< iEntry <<std::endl ;
     const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
     const RAT::DS::MC& rmc= rDS.GetMC(); 
    int particleNum = rmc.GetMCParticleCount();
    //if(particleNum>1)
    {
     const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
     
     u_mc=(rmcparticle.GetMomentum()).Unit();
     pos_mc =rmcparticle.GetPosition();  
     //cout<<pos_mc.X()<<" "<<pos_mc.Y()<<" "<<pos_mc.Z()<<endl;
     double time_mc = rmcparticle.GetTime();
     Int_t nevC = rDS.GetEVCount();
     TVector3 sourcePos;     
//     sourcePos.SetXYZ(0,0,0);
//     pos_mc = sourcePos;
     for(Int_t iev=0;iev<nevC; iev++){
     /// Get the Event information    
      const RAT::DS::EV& rev = rDS.GetEV(iev);  
      std::vector<std::string> fitname = rev.GetFitNames();
      std::vector<std::string>::iterator it;
 
      const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
      const string fitName = "alberta"; 
      if(rev.FitResultExists(fitName)){
		   //std::cout<<"get fitted "<<calpmts.GetFECDCount()<<std::endl;
       RAT::DS::FitVertex fVertex = rev.GetFitResult("alberta").GetVertex(0);
       pos_fit = fVertex.GetPosition();
       RAT::DS::FitVertex fVertex1 = rev.GetFitResult("albertadirection").GetVertex(0);
       u_fit = fVertex1.GetDirection();
       if(fVertex.ValidPosition()){ 
	  ufitX.push_back(u_fit.X());ufitY.push_back(u_fit.Y());ufitZ.push_back(u_fit.Z());
	  double tfit = fVertex.GetTime();
	  Double_t posX=(pos_fit.X());posfitX.push_back(posX);
	  Double_t posY=(pos_fit.Y());posfitY.push_back(posY);
 	  Double_t posZ=(pos_fit.Z());posfitZ.push_back(posZ);
 	  Double_t pos_mcX=(pos_mc.X());posmcX.push_back(pos_mcX);
          Double_t pos_mcY=(pos_mc.Y());posmcY.push_back(pos_mcY);
          Double_t pos_mcZ=(pos_mc.Z());posmcZ.push_back(pos_mcZ);         
	  Double_t tMC = time_mc;
         }
     }
   }
  }
 }
}

double CalculateData(const double *par)  
{
 double chisq = 0;
 const double p0 = 0;
 const double p1 = par[0];
 const double p2 = par[1];
 
 //catch (posfitX.size()!= posmcX.size()) 
 //{
   //std::cout<<"size error: posfitX "<<posfitX.size()<<" posmcX "<<posmcX.size()<<std::endl;
 //}
 
 double sum = 0; 
 for(unsigned int i =0;i<posfitX.size();i++)
 {
	TVector3 POSmc, POSfit,Ufit;
	POSmc.SetXYZ(posmcX[i],posmcY[i],posmcZ[i]);
	POSfit.SetXYZ(posfitX[i],posfitY[i],posfitZ[i]);
	Ufit.SetXYZ(ufitX[i],ufitY[i],ufitZ[i]);
 	double correctedbias = (POSmc-(p1*POSfit+p2*Ufit)).Mag();	 
	cout<<POSfit.X()<<" "<<POSfit.Y()<<" "<<POSfit.Z()<<endl;
 	sum = sum + correctedbias*correctedbias/(1000*1000);
 	//std::cout<<sum<<std::endl;
 }

  chisq = sum;
  return chisq;
}

void MinimizeDrive4()
{
   void RetrieveData(vector<double> &posfitX,vector<double> &posfitY,vector<double> &posfitZ,
   vector<double> &posmcX, vector<double> &posmcY,vector<double> &posmcZ,
   vector<double> &ufitX,vector<double> &ufitY, vector<double> &ufitZ);
   RetrieveData(posfitX,posfitY,posfitZ,posmcX,posmcY,posmcZ,ufitX,ufitY,ufitZ);//global settings
 
   	
  // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined, 
   // kScan, kFumili
   ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
 
   min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
   const int Npar = 2; 
   ROOT::Math::Functor f(&CalculateData,Npar); 
   double step[Npar] = {0.01,0.01};
   double variable[Npar] = {1.0,100.0};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"p1",variable[0], step[0]);
   min.SetVariable(1,"p2",variable[1], step[1]);
//   min.SetVariable(2,"p2",variable[2], step[2]);
   min.Minimize(); 
 
   const double *xs = min.X();
   cout << "Minimum: f(" << xs[0] << "," << xs[1] <<"): " 
        << CalculateData(xs) << endl;
 
}
