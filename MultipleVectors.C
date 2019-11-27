#include <vector>
#include <TMath.h>
#include <TROOT.h>
#include "TVector3.h"
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
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <TMinuit.h>
#include "TString.h"
using namespace ROOT;

int function(vector<int> &a, vector<ROOT::Math::XYZVectorF> &pos) 
{

  for(int i=0;i<5;i++)
  {x.push_back(i);}
  
  for(int i=0;i<5;i++) 
  { ROOT::Math::XYZVectorF temp(i,i*3,i*2); 
    pos.push_back(temp);
  }
  return 0;
}

void MultipleVectors()
{

  vector<int> a;
  vector<ROOT::Math::XYZVectorF> pos;
  int function(vector<int> &a, vector<ROOT::Math::XYZVectorF> &pos);
  
  function(a,pos);
  std::cout<<"data size "<<pos.size()<<std::endl;
  std::cout<<"now print out"<<std::endl;
  for(int i=0;i<5;i++) std::cout<<"a: "<<a[i]<<std::endl;
  for(int i=0;i<5;i++) std::cout<<"pos: "<<pos[i].X()<<" "<<pos[i].Y()<<" "<<pos[i].Z()<<std::endl;
 
}
