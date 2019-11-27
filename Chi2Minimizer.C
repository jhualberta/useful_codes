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
TString filename = "SaveHisto_mpw_r0000106958.root";

int binTot = 2000;
double lowBin = -9000, hiBin = 9000;
double binStep = (hiBin-lowBin)/binTot;
const int iNum = 2000;
double Nfit[iNum];//data
double NfitX[iNum];//data
double NSx[iNum];//Spatial Resolution data
double Nerror[iNum];//errors in data
double Nx[iNum];//x values
double mu_P0;//initial param. mu_P for R(x) 
double scalepeak;
 
double posResol(double *x, double *par)//fix a_e = 0.55 and fit for 3 params only, otherwise hard to converge
{   //par[0]=aE,par[1]=sigma_P, par[2]=mu_P, par[3]=tau_P, par[4]=scale
   double funcPosResol = (1.0-0.55)/(TMath::Sqrt(2*TMath::Pi())*par[0])*TMath::Exp(-0.5*pow((x[0]-par[1])/(par[0]),2))+0.55/(2*par[2])*TMath::Exp(-TMath::Abs(x[0]-par[1])/(par[2]));
   return funcPosResol*1e6;
   cout<<"func "<<x<<" "<<funcPosResol<<endl;
}

double posResolOutput(double *x, double *par)//output, 5 parameters
{   //par[0]=aE,par[1]=sigma_P, par[2]=mu_P, par[3]=tau_P, par[4]=scale
   double funcPosResol = (1.0-par[0])/(TMath::Sqrt(2*TMath::Pi())*par[1])*TMath::Exp(-0.5*pow((x[0]-par[2])/(par[1]),2))+0.55/(2*par[3])*TMath::Exp(-TMath::Abs(x[0]-par[2])/(par[3]));
   return funcPosResol*par[4];
   cout<<"func "<<x<<" "<<funcPosResol<<endl;
}

TH1F *hchi_sigmaP = new TH1F("hchi_sigmaP","chi2 vs sigmaP",400,100,500);
TH1F *hchi_muP = new TH1F("hchi_muP","chi2 vs muP",200,-100,100);
TH1F *hchi_tauP = new TH1F("hchi_tauP","chi2 vs tauP",200,200,400);

void cal_chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //Position Resolution
  int funcpar = 3; 
  TF1 *fun1 = new TF1("fun1",posResol,lowBin,hiBin,funcpar);
  //a_e, sigmaP, muP, tauP, scale
  //fix a_e=0.55  muP=40 //mm
//  fun1->FixParameter(0,0.55);
  fun1->SetParameters(par);
  fun1->SetParLimits(0,10,400);
  fun1->SetParLimits(1,-500,500);
  fun1->SetParLimits(2,10,400);
  //fun1->SetParameters(0.55,par[0],mu_P0,par[1],1e6);
//
  double beginX = (lowBin + lowBin+ (hiBin-lowBin)/binTot)/2;//begin center
  double NR[iNum];//need to calculate S(x)*R(xfit^i-x) bin by bin 

  TH1F *hRx_temp = new TH1F("hRx_temp","tempt Rx",binTot,lowBin,hiBin);	

  for(int ii=0;ii<binTot;ii++) 
  {
   double tempX = beginX+ii*(hiBin-lowBin)/binTot;
   double tempRx = 0;
   if(tempX>=-2000 && tempX<=2000)//cuts
     tempRx = fun1->Eval(tempX);
   else tempRx = 0;
   //if(tempRx<1e-6) tempRx = 0;
   hRx_temp->SetBinContent(ii+1,tempRx);
  }
  hRx_temp->Scale(scalepeak/hRx_temp->GetMaximum());
// use histogram for convolution

  double convl=0;   
  TH1F *htempG = new TH1F("htempG","htempG",binTot,lowBin,hiBin);
  TH1F *hConvl= new TH1F("hConvl","tempt convl",binTot,lowBin,hiBin);

  for(int j=0; j<binTot; j++){
    double val = NSx[j];
    for(int k=0; k<binTot; k++)
    {
	double convVal = hRx_temp->GetBinContent(k) ;
	double t = htempG->GetBinContent(j+k-1);
	htempG->Fill(NfitX[j]+hRx_temp->GetBinCenter(k),val*convVal);
    }
  }

  for(int m=0; m<binTot; m++) {
    hConvl->SetBinContent(m+1,htempG->GetBinContent(m));
  }
  hConvl->Scale(scalepeak/hConvl->GetMaximum());

//use function integral for convolution
/*
  ///for each xfit, loop all Sx to calculate convol
  TH1F *htemp = new TH1F("htemp","htemp",binTot,lowBin,hiBin);
  for(int ii=0;ii<binTot;ii++)
  {
    double xfit = NfitX[ii]; 
    double convl=0;
    for(int jj=0;jj<binTot;jj++)
    {
      double Sx = NSx[jj];	
      double getX = Nx[jj];
      if(getX<=-2000 || getX>=2000) convl += 0;//Sx = 0
      else{
        double Rx = fun1->Eval(xfit-getX); 
        if(Rx<1e-4) Rx = 0;
        convl += Sx*Rx;
      }		 
      //std::cout<<xfit<<" "<<getX<<"  "<<Sx<<"  "<<Rx<<std::endl;	   
    }
    htemp->SetBinContent(ii+1,convl);
  }
  htemp->Scale(scalepeak/htemp->GetMaximum());
*/
  for(int ii=0;ii<iNum;ii++)
  {NR[ii]=hConvl->GetBinContent(ii+1);}//,convl);}//std::cout<<NR[ii]<<", ";}
//  {NR[ii]=htemp->GetBinContent(ii+1);}//std::cout<<NR[ii]<<", ";}

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
    double temp = (NR[ibin]-Nfit[ibin])/((hiBin-lowBin)/binTot);//fix error bin
//    double temp = (NR[ibin]-Nfit[ibin])/Nerror[ibin];
    chi2+=temp*temp;
  }
	 		  
  //std::cout<<"cal chi2  "<<chi2<<std::endl; 
  f = chi2;
  std::cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<f<<std::endl; 
  hchi_sigmaP->Fill(par[0],chi2);
  hchi_muP->Fill(par[1],chi2);
  hchi_tauP->Fill(par[2],chi2);
  return;
}

void Chi2Minimizer()
{  
  TFile *ff = new TFile("fSpatial_NewRat.root");
  TFile *ff2 = new TFile(filename);
  //Spatial Resolution
  TH1F* hSx_temp = (TH1F*)ff->Get("hSx");//!!! 1200 bins
  //Data
  TH1F* hNfit_temp = (TH1F*)ff2->Get("hfitX");//copy directly	
  ///shift, scale and rebin
  int rebinX = 2000/binTot;
  hSx_temp->Rebin(rebinX);hNfit_temp->Rebin(rebinX);
  TH1F* hSx = new TH1F("hSx","rebin hSx",binTot,lowBin,hiBin);
  TH1F* hNfit_temp2 = new TH1F("hNfit_temp2","rebin hfitX",binTot,lowBin,hiBin);//rebin
  TH1F* hNfit = new TH1F("hNfit","rebin hfitX",binTot,lowBin,hiBin);//shifted
  for(int i = 0;i<binTot;i++)
  {
    //cut hSx within [-2000,2000] mm, set others 0
    double xx = hSx_temp->GetXaxis()->GetBinCenter(i+1);
    if(xx<-2000 || xx>2000)//2000 mm
      hSx->SetBinContent(i+1,0);
    else
      hSx->SetBinContent(i+1,hSx_temp->GetBinContent(i));

    hNfit_temp2->SetBinContent(i,hNfit_temp->GetBinContent(i));
  }
  ///now shift the fitX to center
  int max1Bin = hSx->GetMaximumBin();
  hSx->GetXaxis()->SetRangeUser(-20,40);
  int centerBin = hSx->GetMinimumBin();//center "dip"
  hSx->GetXaxis()->UnZoom();//resume hSx
  int max2Bin = hNfit_temp2->GetMaximumBin();//NOTE: peak of fitX
  double xMean = hNfit_temp2->GetBinCenter(max2Bin);
  int meanBin = max2Bin;//default
  //search left and right around peak(+-130 mm) to find a second peak, define the real peak as the mean 
  std::vector<double> search2nd;
  double searchbin = TMath::Nint(130/((hiBin-lowBin)/binTot));
  int maxbinIndex = 0, countTempBin = 0;
  for(int i = -1*searchbin + max2Bin;i<searchbin+max2Bin+1; i++)
  {
    if(i!=max2Bin)//remove max peak to find 2nd peak
      search2nd.push_back(hNfit_temp->GetBinContent(i));
    else { search2nd.push_back(-1); maxbinIndex = countTempBin;}
    countTempBin++;
  }
  //now 2nd peak should be max in search2nd
  std::vector<double>::iterator secondbig= std::max_element(search2nd.begin(), search2nd.end());
  int secondbig_tempBin = std::distance(search2nd.begin(), secondbig);
  int secondbigBin = max2Bin + (secondbig_tempBin - maxbinIndex);//second max peak in histogram
  meanBin = ((secondbigBin+max2Bin)/2);//mean peak
  std::cout<<"peak max,second,mean(bin): "<<max2Bin<<" "<<secondbigBin<<" "<<meanBin<<std::endl;  
  std::cout<<"peak max,second,mean(pos): "<<lowBin+binStep*max2Bin<<" "<<lowBin+binStep*secondbigBin<<" "<<lowBin+binStep*meanBin<<" mm"<<std::endl;
  search2nd.clear();

  int shift = meanBin - centerBin;
  //std::cout<<"shift "<<leftBin<<" "<<meanBin<<" "<<rightBin<<" "<<centerBin<<" "<<shift<<std::endl;
  int ishift;
  if(shift<0) ishift = -1*shift;
  else ishift = 0;
//  for (int i=ishift;i<=binTot;i++) 
  for(int i = 0;i<binTot;i++)
   hNfit->SetBinContent(i,hNfit_temp2->GetBinContent(i+shift));
//  for(int i=0;i<ishift;i++)
//   hNfit->SetBinContent(0);
 
  mu_P0 = hNfit->GetBinCenter(centerBin);//after shift, hNfit mean bin is at hSx center
  std::cout<<"mu_P0 "<<mu_P0<<std::endl;
  hSx->Scale(1000.0/hSx->Integral());
  hNfit->Scale(1000.0/hNfit->Integral());
  hNfit_temp2->Scale(1000.0/hNfit_temp2->Integral());
  double max1 = hSx->GetMaximum();
  double max2 = hNfit->GetMaximum();
  scalepeak = max1;
  hNfit->Scale(scalepeak/max2);//scale to hSx
  hNfit_temp2->Scale(scalepeak/max2);
  hNfit->SetLineColor(kRed);
  hSx->GetXaxis()->UnZoom();
  hSx->GetXaxis()->SetRangeUser(-3000,3000);
  hSx->Draw();//hNfit_temp2->Draw("same");
  hNfit->Draw("sames");
  //input data	
  for(int iter=0;iter<binTot;iter++)
  {Nfit[iter]=hNfit->GetBinContent(iter+1);}	

  for(int iter=0;iter<binTot;iter++)
  {NfitX[iter]=hNfit->GetXaxis()->GetBinCenter(iter+1);}
  
  hNfit->Sumw2();//get sigma_i
  for(int iter=0;iter<binTot;iter++)
  {
    if(hNfit->GetBinError(iter+1) == 0) Nerror[iter]=0.1;
    else Nerror[iter]=hNfit->GetBinError(iter+1);
  }
  
  for(int iter=0;iter<binTot;iter++)
  {NSx[iter]=hSx->GetBinContent(iter+1);}
  
  for(int iter=0;iter<binTot;iter++)
  {Nx[iter]= hSx->GetXaxis()->GetBinCenter(iter+1);}//std::cout<<Nx[iter]<<", ";}

  Int_t npar = 3;//fit sigmaP, tauP, muP
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
  const int Npar = 3; 
  static Double_t vstart[] = {175,mu_P0,310};
  static Double_t step[] = {0.1,0.5,0.1};

  //ptMinuit->mnparm(0, "a_e", vstart[0], step[0], 0,0,ierflg);
  ptMinuit->mnparm(0, "sigmaP", vstart[0], step[0], 0,0,ierflg);
  ptMinuit->mnparm(1, "muP", vstart[1], step[1], 0,0,ierflg);
  ptMinuit->mnparm(2, "tauP", vstart[2], step[2], 0,0,ierflg);
  //ptMinuit->mnparm(4, "scale", vstart[3], step[3], 0,0,ierflg);

  // Now ready for minimization step
  //call Migrad 500 iterations max
  arglist[0] = 5000;//5000;
  arglist[1] = 0.1;//tolerance
  ptMinuit->mnexcm("MIGRAD", arglist ,Npar,ierflg);
  std::cout << "\nPrint results from minuit\n";
  double fParamVal;
  double fParamErr;
  //ptMinuit->GetParameter(0,fParamVal,fParamErr);
  //double a_e = fParamVal;

  ptMinuit->GetParameter(0,fParamVal,fParamErr);
  double sigmaP_r = fParamVal;

  ptMinuit->GetParameter(1,fParamVal,fParamErr);
  double muP_r = fParamVal;

  ptMinuit->GetParameter(2,fParamVal,fParamErr);
  double tauP_r = fParamVal;

//  ptMinuit->GetParameter(3,fParamVal,fParamErr);
//  double scale = fParamVal;

  double a_e = 0.55;
  TF1 *funMain = new TF1("fitPosResol",posResolOutput,-3000,3000,5);
  funMain->SetParameters(a_e,sigmaP_r,muP_r,tauP_r,1);
  double peakfunc = funMain->GetMaximum();
  funMain->SetParameter(4,hNfit->GetMaximum()/peakfunc);
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
      double convVal = hRxOutput->GetBinContent(k) ;
      double t = htempG->GetBinContent(j+k-1);
      htempG->Fill(NfitX[j]+hRxOutput->GetBinCenter(k),val*convVal);
    }
  }

  for(int m=0; m<binTot; m++){
    hConvlOutput->SetBinContent(m+1,htempG->GetBinContent(m));
  }

  hConvlOutput->SetLineColor(kBlue);  
  hConvlOutput->SetLineWidth(2);  
  hConvlOutput->Scale(scalepeak/hConvlOutput->GetMaximum());
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
