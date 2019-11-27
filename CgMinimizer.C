//Conjugate gradient method
#include <vector>
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"

void CgMinimizer()
{
  //double x[10];
  //double x0[10];//trial values
  //for(int i=1;i<11;i++)
  //{
	//x0[i-1]=i/10;
   //}
  
  TMatrixD A(2,2);
  A[0][0]=2,A[0][1]=0,A[1][0]=0,A[1][1]=20;
  
  TMatrixD x0(2,1); x0[0]=2; x0[1]=1;
  TMatrixD x(2,1);
  TMatrixD g0(2,1);
  TMatrixD g(2,1);
  
  
  
  
	
	
	
	
	
	
}

