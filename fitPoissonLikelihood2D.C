#include <iostream>
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
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
#include "TGraph.h"
#include "TError.h"
#include <TMinuit.h>
#include <algorithm>
const int Nfit = 40;
double Ethresh = 5;
double data[Nfit][10];
double mc[Nfit][10];
double Norm;
vector<double> *fcnPtr = new vector<double>();
vector<double> *sigPtr = new vector<double>();
vector<double> *bkgPtr = new vector<double>();
double p_bkg = (1-(-1))/Nfit;

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  int n = Nfit; 

  double lnL = 0.0;
  for (int j=0; j<15-Ethresh; j++)
  {
    for (int i=0; i<n; i++)
    {	    
      double x = data[i][j];
      double pi = mc[i][j];

      double lambda = par[0]*mc[i][j]+par[1]*p_bkg;// lambda = S*pdf+B
      double pdf = 0;
      if ( x >0  && lambda>0 ) {
        pdf = x*log(lambda)-lambda;
      }

      if ( pdf > 0.0 ) {
       lnL += pdf;    // need positive f
       cout << "data, mc, S, B, pdf = "<<x<<" "<<pi<<" "<<par[0]<<" "<<par[1]<<" "<<pdf<<endl;
       fcnPtr->push_back(f);
       sigPtr->push_back(par[0]);
       bkgPtr->push_back(par[1]);
      }
      else {
        // cout << "WARNING -- pdf is negative!!!" << endl;
        // cout << "x, S, B, pdf = " << x << "  " << par[0] << "  " << par[1]<< " "<< pdf << endl;
      }
    }
  }
  f = -2.0 * lnL;         // factor of -2 so minuit gets the errors right
}                         // end of fcn

void fitPoissonLikelihood2D()
{
  int binTot = Nfit;
  TFile *ff = new TFile("analysis_mpw_E5to15.root","read");
  TFile *fpdf = new TFile("PDF_cosThetaSun_solarNueMC_FV5p5.root","read");//!!! always 200 bins

  TH1F *hpdf1 = (TH1F*)fpdf->Get("hcosSun4to5");//pdf
  TH1F *hpdf2 = (TH1F*)fpdf->Get("hcosSun5to6");//pdf
  TH1F *hpdf3 = (TH1F*)fpdf->Get("hcosSun6to7");//pdf
  TH1F *hpdf4 = (TH1F*)fpdf->Get("hcosSun7to8");//pdf
  TH1F *hpdf5 = (TH1F*)fpdf->Get("hcosSun8to9");//pdf
  TH1F *hpdf6 = (TH1F*)fpdf->Get("hcosSun9to10");//pdf
  TH1F *hpdf7 = (TH1F*)fpdf->Get("hcosSun10to15");//pdf

  TH2F* hdata = (TH2F*)ff->Get("hMLPSigCosThetaToSunVsE");
  TH1F* hfitAngle= new TH1F("hfitAngle","fitted",binTot,-1,1);

  hpdf1->SetLineColor(kBlue);
  hdata->RebinX(40/binTot);
  hpdf1->RebinX(40/binTot);
  hpdf2->RebinX(40/binTot);
  hpdf3->RebinX(40/binTot);
  hpdf4->RebinX(40/binTot);
  hpdf5->RebinX(40/binTot);
  hpdf6->RebinX(40/binTot);
  hpdf7->RebinX(40/binTot);
//  cout<<"Rebin: "<<40/binTot<<" "<<40/binTot<<endl;
  Norm = hdata->ProjectionX()->Integral();
  
///check drawing
//   hpdf->Scale(Norm/hpdf->Integral());
//  hpdf->Scale(hdata->GetMaximum()/hpdf->GetMaximum());
//  hpdf->Draw(); hdata->Draw("sames");

  hpdf1->Scale(1./hpdf1->Integral());
  hpdf2->Scale(1./hpdf2->Integral());
  hpdf3->Scale(1./hpdf3->Integral());
  hpdf4->Scale(1./hpdf4->Integral());
  hpdf5->Scale(1./hpdf5->Integral());
  hpdf6->Scale(1./hpdf6->Integral());
  hpdf7->Scale(1./hpdf7->Integral());

  ///!!!!!!!!!!!! cheating 
///  int kk = hdata->GetBinContent(9);
///  hdata->SetBinContent(8,kk);
 
  for(int j = 0;j<15-Ethresh;j++) //loop cosThetaSun
  {
    for(int i = 0;i<binTot; i++) //loop E
    {  
     data[i][j] = hdata->GetBinContent(i+1,j+1);
     double Edata = Ethresh + j;
     // cout<<"E = "<<j+4<<"-"<<j+5<<" MeV, "<<data[i][j]<<endl;
     //if(data[i][j]!=0)
     {
       // mc[i][j] = hpdf1->GetBinContent(i+1);	     
       if(4<=Edata && Edata<5) mc[i][j] = hpdf1->GetBinContent(i+1);
       if(5<=Edata && Edata<6) mc[i][j] = hpdf2->GetBinContent(i+1);
       if(6<=Edata && Edata<7) mc[i][j] = hpdf3->GetBinContent(i+1);
       if(7<=Edata && Edata<8) mc[i][j] = hpdf4->GetBinContent(i+1);
       if(8<=Edata && Edata<9) mc[i][j] = hpdf5->GetBinContent(i+1);
       if(9<=Edata && Edata<10) mc[i][j] = hpdf6->GetBinContent(i+1);
       if(10<=Edata) mc[i][j] = hpdf7->GetBinContent(i+1); // merge into 1 bin
     }
    }
     //cout<<data[i]<<endl;
  }

  const int npar = 2;              // the number of parameters
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];

  par[0] = 190.0;            // a guess
  stepSize[0] = 0.1;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 200;
  parName[0] = "signal";

  par[1] = 50.0;            // a guess
  stepSize[1] = 0.1;       // take e.g. 0.1 of start value
  minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.
  maxVal[1] = 100;
  parName[1] = "background";

  for (int i=0; i<npar; i++){
    minuit.DefineParameter(i, parName[i].c_str(), 
      par[i], stepSize[i], minVal[i], maxVal[i]);
  }

// Do the minimization!

  minuit.Migrad();       // Minuit's best minimization algorithm
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
    minuit.GetParameter(i,outpar[i],err[i]);
  }

// Plot the result.  For this example plot x values as tick marks.

  double xmin = -1.0;
  double xmax = 1.0;
//  TF1* func = new TF1("funcplot", expPdf, xmin, xmax, npar);
//  func->SetParameters(outpar);
//  func->Draw();
//
//  func->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
//  func->SetLineColor(1);             //  black (default)
//  func->SetLineWidth(1);
//
//  func->GetXaxis()->SetTitle("x");
//  func->GetYaxis()->SetTitle("f(x;#xi)");

  TH1F *hFitted = new TH1F("hFitted","fitted", Nfit, xmin, xmax);
  for(int i = 0;i<Nfit;i++) 
  { 
    for(int j = 0;j<15-Ethresh;j++)
    {	    
     int counts = outpar[0]*mc[i][j]+outpar[1]*p_bkg;
     hFitted->SetBinContent(i+1, counts);
    }

  }
  hFitted->SetLineColor(kRed);
  
  double chi2 = 0;
  int countBins = 0;
  for(int i = 0; i<Nfit; i++)
  {
    double fitCount = hFitted->GetBinContent(i+1);
    double dataCount = hdata->GetBinContent(i+1);
          
    chi2 += (fitCount-dataCount)*(fitCount-dataCount);
    countBins++;
  }

  cout<<chi2<<" chi2/ndf "<<chi2/countBins<<endl;

  TCanvas *c = new TCanvas("c","",800,600);
  hdata->ProjectionX()->Draw();
  hFitted->Draw("same"); 
  cout<<fcnPtr->size()<<endl;

  int nsize = fcnPtr->size();
  vector<double> fcnVal;
  vector<double> bkgVal;
  vector<double> sigVal;
  for(int i = 0;i<nsize;i++)
  { 
    fcnVal.push_back(fcnPtr->at(i));
    sigVal.push_back(sigPtr->at(i));
    bkgVal.push_back(bkgPtr->at(i));
  }

//  TGraph *gSig = new TGraph(nsize,(&fcnVal)[0],(&bkgVal)[0]); 
//  TGraph *gBkg = new TGraph(nsize,(&fcnVal)[0],(&sigVal)[0]); 

//  TGraph *gSig = new TGraph(fcnVal->size(), (fcnVal), (&sigVal)[0]);
//  TGraph *gBkg = new TGraph(fcnVal->size(), (fcnVal), (&bkgVal)[0]);

//  TCanvas *c1 = new TCanvas("c1","",800,600);
//  c1->Divide(1,2);
//  c1->cd(1);gSig->Draw("AP*");
//  c1->cd(2);gBkg->Draw("AP*");

//  const double tickHeight = 0.1;
//  TLine* tick = new TLine();
//  for (int i=0; i<Nfit; i++){
//    tick->DrawLine(data[i], 0, data[i], tickHeight);
//  }
  double radius = 5.5;
  double volume = TMath::Pi()*pow(radius,3)*4./3;
  double mass = volume*0.997/1000;
  double day = 92.54;
  cout<<"signal = "<<outpar[0]<<"+-"<<err[0]<<" rate = "<<outpar[0]/(mass*day)<<" +- "<<err[0]/(mass*day)<<endl;
  cout<<"bkg = "<<outpar[1]<<"+-"<<err[1]<<" rate = "<<outpar[1]/(mass*day)<<" +- "<<err[1]/(mass*day)<<endl;

  cout << "To exit, quit ROOT from the File menu of the plot" << endl;
  //theApp.Run(true);
  //canvas->Close();

  //delete canvas, tick, xVecPtr;


}
