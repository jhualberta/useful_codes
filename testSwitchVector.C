#include <iostream>
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <TMinuit.h>
#include <algorithm>
////global data


void test(int Type, vector<double>& hist)
{
  switch(Type) {
   case 0:// 
        for(int i =0;i<5;i++)
        {
          hist.push_back(i);
        }
	break;
   case 1: // smearing sigmaX 
        for(int i =0;i<3;i++)
        {
          hist.push_back(i*2);
        }
        break;
   case 2: // smearing sigmaY 
        for(int i =0;i<2;i++)
        {
	  hist.push_back(i*3);
        }
	break;
   default:
  }
}

void testSwitchVector()
{

  vector<double> xx1;
  vector<double> xx2;
  vector<double> xx3;
  vector<double> xx4;

  test(0,xx1);
  test(1,xx2);
  test(2,xx3);
  test(0,xx4);

  for(vector<double>::iterator it=xx1.begin();it!=xx1.end();it++)
        cout<<*it<<", ";
  cout<<endl; 

  for(vector<double>::iterator it=xx2.begin();it!=xx2.end();it++)
        cout<<*it<<", ";
  cout<<endl;

  for(vector<double>::iterator it=xx3.begin();it!=xx3.end();it++)
        cout<<*it<<", ";
  cout<<endl;
  for(vector<double>::iterator it=xx4.begin();it!=xx4.end();it++)
        cout<<*it<<", ";

  cout<<endl;

}
