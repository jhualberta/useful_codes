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
#include "TError.h"
#include <TMinuit.h>
#include <algorithm>
const int Nfit = 40; //bins
double p_bkg = 1.0/double(Nfit); // uniform background prob.
double NsigTrue = 71;
double NbkgTrue = 38;

double data[Nfit];
double mc[Nfit];
double Norm = 1; //not used !!
vector<double> *fcnPtr = new vector<double>();
vector<double> *sigPtr = new vector<double>();
vector<double> *bkgPtr = new vector<double>();

double pdfBkg[40] = {0.0258900649507,0.025781813808,0.0252946836661,0.0261787346644,0.0260223719028,0.0252345441424,0.0260163579504,0.0255292278085,0.0256976184749,0.0260945393312,0.0253969208564,0.0254751022372,0.0248195814289,0.0248255953813,0.0254751022372,0.0246511907626,0.0254209766659,0.0253187394756,0.0249819581429,0.0245850372865,0.0254570603801,0.025276641809,0.0248195814289,0.0248316093337,0.0248436372384,0.0248616790955,0.0248556651431,0.0243986047631,0.0247414000481,0.0242061582872,0.0244226605725,0.0245369256675,0.0240618234304,0.0250841953332,0.0237611258119,0.0234724560981,0.024837623286,0.024170074573,0.0243264373346,0.0243444791917};
// double pdfBkg[40]={0.518576,0.517647,0.516718,0.515789,0.51486,0.513932,0.513003,0.512074,0.511145,0.510217,0.509288,0.508359,0.50743,0.506501,0.505573,0.504644,0.503715,0.502786,0.501858,0.500929,0.5,0.499071,0.498142,0.497214,0.496285,0.495356,0.494427,0.493499,0.49257,0.491641,0.490712,0.489783,0.488855,0.487926,0.486997,0.486068,0.48514,0.484211,0.483282,0.482353}; /// 20.018576;

double bestLikelihood( double par[] )
{
  int n = Nfit;
  double lnL = 0.0;
  for (int i=0; i<n; i++) {
    double x = data[i];// x = n_i
    double p_i = mc[i];
/// use pdg2020 (40.16)
    double mu_i = par[0]*p_i+par[1]*p_bkg;
    // double mu_i = par[0]*p_i+par[1]*pdfBkg[i];
    double pdf = 0;

    if( x/mu_i>0 )
    {
       pdf = -1*(mu_i - x + x*log(x/mu_i));
    }
    else if (x == 0)
    {
      pdf = -1*(mu_i); //cout<<"zero bin "<<pdf<<endl;
    }

    // cout<<mu_i<<" "<<x/mu_i<<" "<<pdf<<endl;

    if ( pdf != 0.0 ) {
     lnL += pdf;    // need positive f
    }
  } // looping data
  return -2.0*lnL; 
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag) {

  int n = Nfit; 

  double lnL = 0.0;
  for (int i=0; i<n; i++) {
    double x = data[i];// x = n_i
    double p_i = mc[i];
    // if (x == 0) cout<<" zero bins !!"<<endl;
    //// double lambda = par[0]*p_i+par[1]*p_bkg;// lambda = S*pdf+B
    //// double pdf = 0;
    //// if ( x >0  && lambda>0 ) {
    ////   pdf = x*log(lambda)-lambda;
    //// }

/// use pdg2020 (40.16)
    // double mu_i = par[0]*p_i+par[1]*p_bkg;
    double mu_i = par[0]*p_i+par[1]*pdfBkg[i];
    double pdf = 0;

    if( x/mu_i>0 )
    {
       pdf = -1*(mu_i - x + x*log(x/mu_i));
    }
    else if (x == 0)
    {
      pdf = -1*(mu_i); //cout<<"zero bin "<<pdf<<endl;
    }

    // cout<<mu_i<<" "<<x/mu_i<<" "<<pdf<<endl;

    if ( pdf != 0.0 ) {
     lnL += pdf;    // need positive f
     // cout<<"Good fit -- pdf is fine" << endl;
     // cout << "data, mc, S, B, pdf = "<<x<<" "<<p_i<<" "<<par[0]<<" "<<par[1]<<" "<<pdf<<endl;
     fcnPtr->push_back(f);
     sigPtr->push_back(par[0]);
     bkgPtr->push_back(par[1]);
    }
    else {
    //  cout<<"WARNING -- pdf is negative!!!" << endl;
    //  cout<<"data, mc, S, B, pdf = "<<x<<" "<<p_i<<" "<<par[0]<<" "<<par[1]<<" "<<pdf<<endl;
    }
  }
  f = -2.0 * lnL;         // factor of -2 so minuit gets the errors right
}                         // end of fcn

void loopFitPoissonLikelihood1() // for ensemble test
{
  int binTot = Nfit;
  TFile *ff = new TFile("ensembleTMVA_output500000evts_poisson.root","read");
  TFile *fpdf = new TFile("PDF_MPW_cosThetaSun_solarNueMC_FV5p5.root","read");//!!! always 200 bins

  TH1F *hpdf = (TH1F*)fpdf->Get("hcosSun5to15");//pdf
//  TH1F* hdata = (TH1F*)ff->Get("hcos_5to15");
//  TH1F* hdata = (TH1F*)ff->Get("hMLPSigCosThetaToSun");
  int iseed = 0;
  const int Nseed = 500000;
  double vSignal[Nseed];
  double vSignal_sigma[Nseed];
  double vBkg[Nseed];
  double vBkg_sigma[Nseed];
  double vChi2[Nseed];
  double vLbest[Nseed];
  for(; iseed<Nseed; iseed++)
  {	  
    TString tName;
    tName.Form("T2_%d",iseed);//"hcosTheta_test%d",iseed);
    if (iseed%1000 == 0)
             cout<<"seed:   "<<iseed<<endl;

    TString tNameSig;
    TString tNameBkg;

    tNameSig.Form("hcosTheta_sig%d",iseed);
    tNameBkg.Form("hcosTheta_bkg%d",iseed);

    TH1F* hdata = new TH1F("hdata","",40,-1,1);
    TTree *ttt = (TTree*)(ff->Get(tName));
    (TH1F*)ttt->Project("hdata","cosThetaToSun");//hbdtE5to15_sig, hmlpE5to15_sig, hcos_default_E5to15
//    TH1F* hfitAngle= new TH1F("hfitAngle","fitted",binTot,-1,1);

    TH1F* hsig = (TH1F*)ff->Get(tNameSig);
    TH1F* hbkg = (TH1F*)ff->Get(tNameBkg);

    double numberSig = hsig->GetEntries();
    double numberBkg = hbkg->GetEntries();
 

    hpdf->SetLineColor(kBlue);
    hdata->Rebin(40/binTot);hpdf->Rebin(40/binTot);
    //cout<<"Rebin: "<<40/binTot<<" "<<40/binTot<<endl;
    Norm = hdata->Integral();
    
// //check drawing
//    hpdf->Scale(Norm/hpdf->Integral());
//    hpdf->Scale(hdata->GetMaximum()/hpdf->GetMaximum());
//    hpdf->Draw(); hdata->Draw("sames");

    hpdf->Scale(1.0/hpdf->Integral());

    ///!!!!!!!!!!!! cheating 
    //int kk = hdata->GetBinContent(9);
    //hdata->SetBinContent(8,kk);
 
    for(int i = 0;i<Nfit;i++)
    {
       data[i] = hdata->GetBinContent(i+1);
       mc[i] = hpdf->GetBinContent(i+1);
       //cout<<data[i]<<endl;
    }

    const int npar = 2;              // the number of parameters
    TMinuit minuit(npar);
    minuit.SetFCN(fcn);
    minuit.SetPrintLevel(-1); //quiet print nothing
    double par[npar];               // the start values
    double stepSize[npar];          // step sizes 
    double minVal[npar];            // minimum bound on parameter 
    double maxVal[npar];            // maximum bound on parameter
    string parName[npar];

    par[0] = numberSig;            // a guess
    stepSize[0] = 0.01;       // take e.g. 0.1 of start value
    minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.
    maxVal[0] = 109;
    parName[0] = "signal";

    par[1] = numberBkg;            // a guess
    stepSize[1] = 0.01;       // take e.g. 0.1 of start value
    minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.
    maxVal[1] = 109;
    parName[1] = "background";

    for (int i=0; i<npar; i++){
      minuit.DefineParameter(i, parName[i].c_str(), 
        par[i], stepSize[i], minVal[i], maxVal[i]);
    }

//   Do the minimization!

    minuit.Migrad();       // Minuit's best minimization algorithm
    double outpar[npar], err[npar];
    for (int i=0; i<npar; i++){
      minuit.GetParameter(i,outpar[i],err[i]);
    }

//   Plot the result.  For this example plot x values as tick marks.

    double xmin = -1.0;
    double xmax = 1.0;
//    TF1* func = new TF1("funcplot", expPdf, xmin, xmax, npar);
//    func->SetParameters(outpar);
//    func->Draw();
//  
//    func->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
//    func->SetLineColor(1);             //  black (default)
//    func->SetLineWidth(1);
//  
//    func->GetXaxis()->SetTitle("x");
//    func->GetYaxis()->SetTitle("f(x;#xi)");
    // cout<<NsigTrue<<" "<<NbkgTrue<<endl;

    TH1F *hFitted = new TH1F("hFitted","fitted", Nfit, xmin, xmax);
    for(int i = 0;i<Nfit;i++) 
    {
        double counts = outpar[0]*mc[i]+outpar[1]*p_bkg;
	// cout<<"bin "<<i<<": sig= "<<outpar[0]*mc[i]<<" + bkg = "<<outpar[1]*p_bkg<<" counts = "<<counts<<endl;
        hFitted->SetBinContent(i+1, counts);
    }
    hFitted->SetLineColor(kRed);
    hFitted->SetLineWidth(2);
    hFitted->GetXaxis()->SetTitle("cos#theta_{sun}");
    hFitted->GetXaxis()->SetTitle("Counts/92.54 days/0.05 bin");
    double chi2_40bin = 0;
    double chi2 = 0;
    double countChi2Bins = 0;
    hdata->GetSumw2()->Set(1); 
    hdata->Sumw2(kTRUE);
    for(int i = 0; i<Nfit; i++)
    {
      double fitCount = hFitted->GetBinContent(i+1);
      double dataCount = hdata->GetBinContent(i+1);

      double delta = (fitCount-dataCount);
      double error = hdata->GetBinError(i);
      if(error == 0) error = 1;
         //cout<<"delta = "<<delta<<", delta2 = "<<delta*delta<<" error= "<<error<<endl;
      chi2_40bin += delta*delta/(error*error);
      if(dataCount>0) 
      {
         if(error == 0) error = 1;
         //cout<<"delta = "<<delta<<", delta2 = "<<delta*delta<<" error= "<<error<<endl;
         chi2 += delta*delta/(error*error);
         countChi2Bins++;
      }
    }

    //*----- don't draw -----*//
    double Lbest = bestLikelihood(outpar);
    vLbest[iseed] = Lbest;
    // cout<<"L_best = "<<Lbest<<endl;
    // cout<<"non-zero bins = "<<countChi2Bins<<endl;
    // cout<<"chi2/ndf "<<chi2<<"/"<<countChi2Bins<<" "<<chi2/countChi2Bins<<endl;
    // cout<<"chi2/ndf (40 bin) "<<chi2_40bin<<"/40 "<<chi2_40bin/40<<endl;
    vChi2[iseed] = chi2_40bin/40;

    // TCanvas *c = new TCanvas("c","",800,600);
    // hpdf->Scale(Norm);
    // hpdf->Draw();
    // hdata->Draw("same");
    // hFitted->Draw("same");

    // cout<<"fcn calcu "<<fcnPtr->size()<<endl;

    //int nsize = fcnPtr->size();
    //vector<double> fcnVal;
    //vector<double> bkgVal;
    //vector<double> sigVal;
    //for(int i = 0;i<nsize;i++)
    //{ 
    //  fcnVal.push_back(fcnPtr->at(i));
    //  sigVal.push_back(sigPtr->at(i));
    //  bkgVal.push_back(bkgPtr->at(i));
    //}

    //TGraph *gSig = new TGraph(nsize,(&fcnVal)[0],(&bkgVal)[0]); 
    //TGraph *gBkg = new TGraph(nsize,(&fcnVal)[0],(&sigVal)[0]); 

//    TGraph *gSig = new TGraph(fcnVal->size(), (fcnVal), (&sigVal)[0]);
//    TGraph *gBkg = new TGraph(fcnVal->size(), (fcnVal), (&bkgVal)[0]);

    // TCanvas *c1 = new TCanvas("c1","",800,600);
    // c1->Divide(1,2);
    // c1->cd(1);gSig->Draw("AP*");
    // c1->cd(2);gBkg->Draw("AP*");

//    const double tickHeight = 0.1;
//    TLine* tick = new TLine();
//    for (int i=0; i<Nfit; i++){
//      tick->DrawLine(data[i], 0, data[i], tickHeight);
//    }
    double radius = 5.5;
    double volume = TMath::Pi()*pow(radius,3)*4./3;
    double mass = volume*0.997/1000;
    double day = 92.54;
    /// cout<<"actual number: "<<numberSig<<", "<<numberBkg<<endl;
    /// cout<<"signal = "<<outpar[0]<<"+-"<<err[0]<<" rate = "<<outpar[0]/(mass*day)<<" +- "<<err[0]/(mass*day)<<endl;
    /// cout<<"bkg = "<<outpar[1]<<"+-"<<err[1]<<" rate = "<<outpar[1]/(mass*day)<<" +- "<<err[1]/(mass*day)<<endl;
    vSignal[iseed] = outpar[0];
    vSignal_sigma[iseed] = err[0];
    vBkg[iseed] = outpar[1];
    vBkg_sigma[iseed] = err[1];
    //cout<<vSignal[iseed]<<" "<<vSignal_sigma[iseed]<<endl;
    //
    delete hFitted;
    delete hdata;
   } // loop of the seed

   TH1F *hpull = new TH1F("hpull","fit pull, sig", 200, -10,10);
   TH1F *hbias = new TH1F("hbias","bias, sig", 200, -1,1);

   TH1F *hpullBkg = new TH1F("hpullBkg","fit pull, bkg", 200, -10,10);
   TH1F *hbiasBkg = new TH1F("hbiasBkg","bias, bkg", 200, -1,1);

   TH1F *hchi2 = new TH1F("hchi2", "chi2/40",200,0,2);
   
   TH1F *hLbest = new TH1F("hLbest", "best -2logL",1000,0,100);
   
   for(int iseed = 0; iseed<Nseed; iseed++)
   {

       TString tNameSig;
       TString tNameBkg;

       tNameSig.Form("hcosTheta_sig%d",iseed);
       tNameBkg.Form("hcosTheta_bkg%d",iseed);

       TH1F* hsig = (TH1F*)ff->Get(tNameSig);
       TH1F* hbkg = (TH1F*)ff->Get(tNameBkg);

       NsigTrue = hsig->GetEntries();
       NbkgTrue = hbkg->GetEntries();
	   
       double pull = (vSignal[iseed] - NsigTrue)/vSignal_sigma[iseed];
       double bias = (vSignal[iseed] - NsigTrue)/vSignal[iseed];
       // cout<<vSignal[iseed]<<" "<<vBkg[iseed]<<" "<<vBkg_sigma[iseed]<<" "<<pull<<" "<<bias<<endl;
       hpull->Fill(pull); 
       hbias->Fill(bias);
       double pullBkg = (vBkg[iseed] - NbkgTrue)/vBkg_sigma[iseed];
       double biasBkg = (vBkg[iseed] - NbkgTrue)/vBkg[iseed];
       hpullBkg->Fill(pullBkg);
       hbiasBkg->Fill(biasBkg);
       hchi2->Fill(vChi2[iseed]);
       hLbest->Fill(vLbest[iseed]);
   }
   
   TCanvas *c = new TCanvas("c","",800,600);
   hpull->Draw();//hpullBkg->SetLineColor(kRed);hpullBkg->Draw("sames");
   TCanvas *c1 = new TCanvas("c1","",800,600);
   hbias->Draw();
  
   TCanvas *c2 = new TCanvas("c2","",800,600);
   hpullBkg->Draw();
   TCanvas *c3 = new TCanvas("c3","",800,600);
   hbiasBkg->Draw();

   TCanvas *c4 = new TCanvas("c4","best -2lnL",800,600);
   hLbest->Draw();

   TFile *fnew = new TFile("record_ensemble500000Test.root","recreate");
   hbias->Write();
   hpull->Write();
   hbiasBkg->Write();
   hpullBkg->Write();
   hchi2->Write();
   hLbest->Write();
    // cout << "To exit, quit ROOT from the File menu of the plot" << endl;
  //theApp.Run(true);
  //canvas->Close();

  //delete canvas, tick, xVecPtr;
}
