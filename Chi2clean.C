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
//global data
TString filename = "SaveHisto_Rat_100934.root";//SaveHisto_Jeff_107055.root";
//TString fileSx = "fSpatial_1e5evts_100934.root";
//fSpatial_1e5evts_Etot.root

double lowBin = -10000, hiBin = 10000;
const int iNum = 2000;
int binTot = iNum;
//fit region for X mm
double hiCut = 2000; 
double loCut = -hiCut;
double Nfit[iNum];//data
double NfitX[iNum];//data
double NSx[iNum];//Spatial Resolution data
double Nerror[iNum];//errors in data
double Nx[iNum];//x values
double mu_P0;//initial param. mu_P for R(x) 
double scalepeak;
TH1F *hchi_sigmaP = new TH1F("hchi_sigmaP","chi2 vs sigmaP",400,100,500);
TH1F *hchi_muP = new TH1F("hchi_muP","chi2 vs muP",200,-100,100);
TH1F *hchi_tauP = new TH1F("hchi_tauP","chi2 vs tauP",200,200,400);
 
double posResol(double *x, double *par)
{   //par[0]=aE,par[1]=sigma_P, par[2]=mu_P, par[3]=tau_P, par[4]=scale

/*Fit 3: for fixed aE, scale */ //fix a_e = 0.55 and fit for 3 params only, otherwise hard to converge
   double funcPosResol = (1.0-0.55)/(TMath::Sqrt(2*TMath::Pi())*par[0])*TMath::Exp(-0.5*pow((x[0]-par[1])/(par[0]),2))+0.55/(2*par[2])*TMath::Exp(-TMath::Abs(x[0]-par[1])/(par[2]));
   return funcPosResol*par[3];
}

void cal_chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //Position Resolution

  //Fit 4 
  int funcpar = 4;
  TF1 *fun1 = new TF1("fun1",posResol,lowBin,hiBin,funcpar); 
  fun1->SetParLimits(0,100,400);
  fun1->SetParLimits(1,-500,500);
  fun1->SetParLimits(2,100,400);
  fun1->SetParameters(par);

  double beginX = (lowBin + lowBin+ (hiBin-lowBin)/binTot)/2;//begin center
  double NR[iNum];//need to calculate S(x)*R(xfit^i-x) bin by bin 
  TH1F *hRx_temp = new TH1F("hRx_temp","tempt Rx",binTot,lowBin,hiBin);	

  for(int ii=0;ii<binTot;ii++) 
  {
   double tempX = beginX+ii*(hiBin-lowBin)/binTot;
   double tempRx = 0;
   if(tempX>=loCut && tempX<=hiCut)//cuts
     tempRx = fun1->Eval(tempX);
   else tempRx = 0;
   //if(tempRx<1e-6) tempRx = 0;
   hRx_temp->SetBinContent(ii+1,tempRx);
  }
  //hRx_temp->Scale(scalepeak/hRx_temp->GetMaximum());
// use histogram for convolution

  double convl=0;   
  TH1F *htempG = new TH1F("htempG","htempG",binTot,lowBin,hiBin);
  TH1F *hConvl= new TH1F("hConvl","tempt convl",binTot,lowBin,hiBin);

  for(int j=0; j<binTot; j++){
    double val = NSx[j];
    for(int k=0; k<binTot; k++)
    {
  double convVal = hRx_temp->GetBinContent(k);
  double t = htempG->GetBinContent(j+k-1);
  htempG->Fill(NfitX[j]+hRx_temp->GetBinCenter(k),val*convVal);
    }
  }

  for(int m=0; m<binTot; m++) {
    hConvl->SetBinContent(m+1,htempG->GetBinContent(m));
  }
  //hConvl->Scale(scalepeak/hConvl->GetMaximum());

  for(int ii=0;ii<iNum;ii++)
  {NR[ii]=hConvl->GetBinContent(ii+1);}

  //htemp->Draw();
  delete hRx_temp;
  delete htempG; delete hConvl;
  //delete htemp;

//delete fun1;	 
//calculate chi2
  double chiVal = 0;
  double chi2 = 0;
  for(int ibin=0;ibin<iNum;ibin++)
  {
    double temp = (NR[ibin]-Nfit[ibin])/((hiBin-lowBin)/binTot*100);//fix error bin
//    double temp = (NR[ibin]-Nfit[ibin])/Nerror[ibin];
    if(Nfit[ibin]>10)//cout on Counts
    chi2+=temp*temp;
    else continue;
  }
	 		  
  //std::cout<<"cal chi2  "<<chi2<<std::endl; 
  f = chi2;
  std::cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<f<<std::endl; 
  hchi_sigmaP->Fill(par[0],chi2);
  hchi_muP->Fill(par[1],chi2);
  hchi_tauP->Fill(par[2],chi2);
  return;
}

void Chi2clean()
{  
  TFile *ff = new TFile("fSpatial_1e5evts_100934.root");
  TFile *ff2 = new TFile(filename);
  //original histo, Spatial Resolution
  TH1F* hSx_temp = (TH1F*)ff->Get("hSx");
  //original histo, Data
  TH1F* hNfit_temp = (TH1F*)ff2->Get("hfitX");//copy directly	
  ///shift, scale and rebin
  int rebinX = 2000/binTot;
  hSx_temp->Rebin(rebinX);hNfit_temp->Rebin(rebinX);
  TH1F* hSx = new TH1F("hSx","rebin hSx",binTot,lowBin,hiBin);
  TH1F* hNfit_temp2 = new TH1F("hNfit_temp2","rebin hfitX",binTot,lowBin,hiBin);//rebin
  TH1F* hNfit = new TH1F("hNfit","rebin hfitX",binTot,lowBin,hiBin);//shifted
  double mu_fit = hNfit_temp2->GetMean();
  for(int i = 0;i<binTot;i++)
  {
    //cut hSx within [-2000,2000] mm, set others 0
    double xx = hSx_temp->GetXaxis()->GetBinCenter(i+1);
    if(xx<loCut+mu_fit || xx>hiCut+mu_fit) 
    {//2000 mm
      hSx->SetBinContent(i+1,0);
    } 
    else{
      hSx->SetBinContent(i+1,hSx_temp->GetBinContent(i));
    }
    double xxfit = hNfit_temp2->GetXaxis()->GetBinCenter(i+1);
    if(xxfit<-2000+mu_fit|| xx>2000+mu_fit)
      hNfit_temp2->SetBinContent(i+1,0);
    else
    hNfit_temp2->SetBinContent(i,hNfit_temp->GetBinContent(i));
  }
  /// Find rough centers for hSx and hNfit 
  //int max1Bin = hSx->GetMaximumBin();
  hSx_temp->GetXaxis()->SetRangeUser(-20,40);
  int centerBinSx = hSx_temp->GetMinimumBin();//center "dip"
  //std::cout<<"center X Sx "<<hSx->GetBinCenter(centerBinSx)<<std::endl;
  hSx_temp->GetXaxis()->UnZoom();//resume hSx

  int meanBinFitX = hNfit_temp2->GetXaxis()->FindBin(hNfit_temp2->GetMean());  

  //hSx->Scale(1e4/hSx->Integral());
  //hNfit->Scale(1e4/hNfit->Integral());

  int shift = centerBinSx - meanBinFitX;

  //std::cout<<"shift "<<leftBin<<" "<<meanBin<<" "<<rightBin<<" "<<centerBin<<" "<<shift<<std::endl;
  int ishift;
  if(shift<0) ishift = -1*shift;
  else ishift = 0;
//  for (int i=ishift;i<=binTot;i++) 
  for(int i = 0;i<binTot;i++)
   hNfit->SetBinContent(i,hNfit_temp2->GetBinContent(i+0));
//  for(int i=0;i<ishift;i++)
//   hNfit->SetBinContent(0);
 
//  mu_P0 = hNfit->GetBinCenter(centerBin);//after shift, hNfit mean bin is at hSx center
  mu_P0 = hNfit->GetBinCenter(meanBinFitX) - hSx->GetBinCenter(centerBinSx);

  std::cout<<"mu_P0 "<<mu_P0<<std::endl;
  double max1 = hSx->GetMaximum();//GetBinContent(centerBinSx);
  double max2 = hNfit->GetBinContent(meanBinFitX);
  hNfit->Scale(max1/max2);//scale to hSx

  hNfit->SetLineColor(kRed);
  //hSx->GetXaxis()->UnZoom();
  hSx->GetXaxis()->SetRangeUser(-2000,2000);
  hSx->Draw();//hNfit_temp2->Draw("same");
  hNfit->Draw("sames");
  //input data	
  for(int iter=0;iter<binTot;iter++)
  {
    if(hNfit->GetBinContent(iter+1)>5)
    Nfit[iter]=hNfit->GetBinContent(iter+1);
    else Nfit[iter]=0;
  }	

  for(int iter=0;iter<binTot;iter++)
  {NfitX[iter]=hNfit->GetXaxis()->GetBinCenter(iter+1);}
  
  //hNfit->Sumw2();//get sigma_i
  for(int iter=0;iter<binTot;iter++)
  {
    if(hNfit->GetBinError(iter+1) == 0) Nerror[iter]=0.1;
    else Nerror[iter]=hNfit->GetBinError(iter+1);
  }

  for(int iter=0;iter<binTot;iter++)
  {NSx[iter]=hSx->GetBinContent(iter+1);}

  for(int iter=0;iter<binTot;iter++)
  {Nx[iter]= hSx->GetXaxis()->GetBinCenter(iter+1);}//std::cout<<Nx[iter]<<", ";}

  Int_t npar = 4;//fit sigmaP, tauP, muP
  TMinuit *ptMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
  //
  //  select verbose level:
  //    default :     (58 lines in this test)
  //    -1 : minimum  ( 4 lines in this test)
  //     0 : low      (31 lines)
  //     1 : medium   (61 lines)
  //     2 : high     (89 lines)
  //     3 : maximum (199 lines in this test)
  //

  //-1 no output; 1 standard output
  ptMinuit->SetPrintLevel(1);
  //1 for chi2; 0.5 for negative log-likelihood
  ptMinuit->SetErrorDef(1);
  // set the user function that calculates chi_square (the value to minimize)
  ptMinuit->SetFCN(cal_chi2);

  Double_t arglist[10];
  Int_t ierflg = 0;
  //1 standard; 2 improve but slower
  arglist[0] = 2;
  ptMinuit->mnexcm("SET STR", arglist, 1, ierflg);

  //fix and releasing parameters
  //arglist[0] = 0.55;
  //minuit.mnexcm("FIX ", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  //Fit 3
  const int Npar = 4;

  static Double_t vstart[] = {170, mu_P0, 320, 20};
  static Double_t step[] = {0.1,0.5,0.1,0.1};

  ptMinuit->mnparm(0, "sigmaP", vstart[0], step[0], 0,0,ierflg);
  ptMinuit->mnparm(1, "muP", vstart[1], step[1], 0,0,ierflg);
  ptMinuit->mnparm(2, "tauP", vstart[2], step[2], 0,0,ierflg);
  ptMinuit->mnparm(3, "scale", vstart[3], step[3], 0,0,ierflg);
  // Now ready for minimization step
  //call Migrad 500 iterations max
  arglist[0] = 5000;//5000;
  arglist[1] = 1;//tolerance
  ptMinuit->mnexcm("MIGRAD", arglist ,Npar,ierflg);
  std::cout << "\nPrint results from minuit\n";
  double fParamVal;
  double fParamErr;

  ptMinuit->GetParameter(0,fParamVal,fParamErr);
  double sigmaP_r = fParamVal;

  ptMinuit->GetParameter(1,fParamVal,fParamErr);
  double muP_r = fParamVal;

  ptMinuit->GetParameter(2,fParamVal,fParamErr);
  double tauP_r = fParamVal;

  ptMinuit->GetParameter(3,fParamVal,fParamErr);
  double scale = fParamVal;
  // Set starting values and step sizes for parameters

  TF1 *funMain = new TF1("fitPosResol",posResol,-3000,3000,Npar);

  funMain->SetParameters(sigmaP_r,muP_r,tauP_r,scale);

  double peakfunc = funMain->GetMaximum();
  //funMain->SetParameter(5,hNfit->GetMaximum()/peakfunc);
  //std::cout<<"peak resol "<<peakfunc<<" peak fitX "<<hNfit->GetMaximum()<<std::endl;
  funMain->Draw("sames");

  TH1F *htempG = new TH1F("htempG","htempG",binTot,lowBin,hiBin);
  TH1F *hConvlOutput= new TH1F("hConvlOutput","output convl",binTot,lowBin,hiBin);
  TH1F *hRxOutput= new TH1F("hRxOutput","output hRx",binTot,lowBin,hiBin);

  for(int ii=0;ii<binTot;ii++)
  {
   double tempRx = funMain->Eval(lowBin+ii*(hiBin-lowBin)/binTot);
   if(tempRx<1e-6) tempRx = 0;
   hRxOutput->SetBinContent(ii+1,tempRx);
  }

  for(int j=0; j<binTot; j++){
    double val = NSx[j];
    for(int k=0; k<binTot; k++)
    {
      double convVal = hRxOutput->GetBinContent(k);
      double t = htempG->GetBinContent(j+k-1);
      htempG->Fill(NfitX[j]+hRxOutput->GetBinCenter(k),val*convVal);
    }
  }

  for(int m=0; m<binTot; m++){
    hConvlOutput->SetBinContent(m+1,htempG->GetBinContent(m));
  }

  hConvlOutput->SetLineColor(kBlue);  
  hConvlOutput->SetLineWidth(2);  
  //hConvlOutput->Scale(scalepeak/hConvlOutput->GetMaximum());
  hConvlOutput->Draw("sames");
  
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  std::cout << "\n";
  std::cout << " Minimum chi square = " << amin << "\n";
  std::cout << " Estimated vert. distance to min. = " << edm << "\n";
  std::cout << " Number of variable parameters = " << nvpar << "\n";
  std::cout << " Highest number of parameters defined by user = " << nparx << "\n";
  std::cout << " Status of covariance matrix = " << icstat << "\n";
  std::cout << "\n";

  //TCanvas c2("c2","",500,600);
  //c2.cd();
  ////Plot chi2 for a param.
  //arglist[0] = 0;
  //ptMinuit->mnexcm("SCAN", arglist,1,ierflg);

  TString runID(filename(filename.Index("r0"),11));  
  TString newfile = "Convolved_"+ runID +".root";
  TFile *fconv = new TFile(newfile,"recreate");
  fconv->cd();
  hSx->Write();hNfit->Write();  
  funMain->Write();
  hConvlOutput->Write(); 
  hchi_muP->Write();
  hchi_sigmaP->Write();
  hchi_tauP->Write();

}
