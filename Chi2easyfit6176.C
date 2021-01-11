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
#include <iomanip>      // std::setprecision
//global data
TString filename = "Merged_MPWn16_WaterMP6176_nhit6_Calibration_r0000107055_p009.root";//Merged_MPW_nhit6_Calibration_r106925_s0.root";//Merged_RatWater_Calibration_r0000107055_p009.root";//Merged_MPW_nhit6_Calibration_r106925_s0.root";
TString fileSx = "Merged_fSpatial_New6176_test_6and7MeVgamma_1e5evts.root";
int choice = 1;// 1: posx, 2:posy, 3:posz
TString hSxName;
bool defaultFit = false;

double rangeDefault = 1500;
double fitRange = rangeDefault;
double lowBin = -10000, hiBin = 10000;// the bin width must be 10 mm/bin to match hSx !!!
const int iNum = 2000;
int binTot = iNum;
//fit region for X mm
double hiCut = 3000;
double loCut = -hiCut;
double Nfit[iNum];//data
double NfitX[iNum];//data
double NSx[iNum];//Spatial Resolution data
double Nerror[iNum];//errors in data
double Nx[iNum];//x values
double mu_P0;//initial param. mu_P for R(x) 
double scalepeak;

double posResol(double *x, double *par)
{   // par[0]=mu_P, par[1]=sigma_P, par[2]=tau_P, par[4]=scale

//Fit 3: for fixed aE, scale //fix a_e = 0.45 and fit for 3 params only, otherwise hard to converge
   double funcPosResol = (1.0-par[3])/(TMath::Sqrt(2*TMath::Pi())*par[1])*TMath::Exp(-0.5*pow((x[0]-par[0])/(par[1]),2))+par[3]/(2*par[2])*TMath::Exp(-TMath::Abs(x[0]-par[0])/(par[2]));
   return funcPosResol*par[4];
}

double cal_convl(double *x, double *par)
{
  //Position Resolution
  double xx = x[0];
  //Fit 4 
  int funcpar = 5;
  TF1 *fun1 = new TF1("fun1",posResol,lowBin,hiBin,funcpar); 
  fun1->SetParLimits(0,100,400);
  fun1->SetParLimits(1,-500,500);
  fun1->SetParLimits(2,100,400);
  fun1->SetParLimits(4,0.3,0.8);
  fun1->SetParameters(par);

  //double beginX = (lowBin + lowBin+ (hiBin-lowBin)/binTot)/2;//begin center
  double convl = 0, xi=0, sxi=0, rxi=0; 
 
  for(int ii=0;ii<binTot;ii++) 
  {
     xi = Nx[ii];//beginX+ii*(hiBin-lowBin)/binTot;
     sxi = NSx[ii];
     rxi = fun1->Eval(xx-xi);//Wrong!!
     convl = convl+ sxi*rxi;
     //cout<<convl<<" "<<sxi<<" "<<rxi<<endl;
  }
  return convl;
}

void Chi2easyfit6176()
{
  TFile *ff = new TFile(fileSx);//track simulation
  TFile *ff2 = new TFile(filename);
  //original histo, Spatial Resolution
  //
  if(choice == 1) hSxName = "hSx";
  if(choice == 2) hSxName = "hSy";
  if(choice == 3) hSxName = "hSz";
  
  TH1F* hSx_temp = (TH1F*)ff->Get(hSxName);
  //original histo, Data
  TTree *Tdata = (TTree*)ff2->Get("T");
  TH1F* hNfit_temp = new TH1F("hNfit_temp","rebin "+hSxName,binTot,lowBin,hiBin);
  TH1F* hNfit= new TH1F("hNfit","rebin hfitX",binTot,lowBin,hiBin);//fitted energy
  double posx, posy, posz, energy, posRad, itr;
  UInt_t nhits;
  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("posRad", &posRad);
  Tdata->SetBranchAddress("itr", &itr);
  for(int i =0;i<Tdata->GetEntries();i++)
  {
    Tdata->GetEntry(i);
    //cout<<energy<<" "<<nhits<<endl;
    if(nhits>5 )
    {
      if(choice == 1) hNfit_temp->Fill(posx);	    
      if(choice == 2) hNfit_temp->Fill(posy);
      if(choice == 3) hNfit_temp->Fill(posz);
    }
  }

  ///shift, scale and rebin
  int rebinX = ( (hiBin-lowBin)/binTot )/10; // origin 10 mm/bin to x mm/bin
  if(rebinX == 0) rebinX = 1;
  //cout<<"rebinX "<<rebinX<<endl;
  // hSx_temp->Rebin(2);hNfit_temp->Rebin(2);
  //hSx_temp->Rebin(rebinX);//hNfit_temp->Rebin(rebinX);
  TH1F* hSx = new TH1F(hSxName,"rebin hSx",binTot,lowBin,hiBin);
  TH1F* hNfit_temp2 = new TH1F("hNfit_temp2","rebin hfitZ",binTot,lowBin,hiBin);//rebin

  // double muSx = hSx_temp->GetMean();

  /// calculate shifts between hNfit and hSx
  TF1 *fgaus = new TF1("fgaus","gaus",-9000,9000);
  int shift = 0;
  double mu_P0 = 0;
  int centerBinSx = 0;
  hNfit_temp->Fit("fgaus");
  double mean_hNfit = fgaus->GetParameter(1);
  int meanBinFitX = hNfit_temp->GetXaxis()->FindBin(mean_hNfit);
  cout<<"mean of hNfit "<<mean_hNfit<<" bin: "<<meanBinFitX<<endl;
  //if(hSxName == "hSz")
  //{
  hSx_temp->Fit("fgaus","q");
  double muSx = fgaus->GetParameter(1);
  int meanBin_hSx = hNfit_temp2->GetXaxis()->FindBin(hNfit_temp2->GetMean());
  // shift = meanBin_hSx - meanBinFitX;
  // mu_P0 = mean_hNfit - mean_hSx; 
  //}
  //else 
  for(int i = 0;i<binTot;i++)
  {
     //cut hSx within [-2000,2000] mm, set others 0
     double xx = hSx_temp->GetXaxis()->GetBinCenter(i);
     if(xx<loCut+muSx || xx>hiCut+muSx) 
     {//2000 mm
       hSx->SetBinContent(i,0);
     } 
     else{
       hSx->SetBinContent(i,hSx_temp->GetBinContent(i));
     }
     //cut hNfit to lowBin, hiBin
     double xxfit = hNfit_temp2->GetXaxis()->GetBinCenter(i);//(i+1);
     if(xxfit<loCut+mean_hNfit|| xx>hiCut+mean_hNfit)
       hNfit_temp2->SetBinContent(i+1,0);
     else
     hNfit_temp2->SetBinContent(i,hNfit_temp->GetBinContent(i));
  }
  /// Find rough centers for hSx and hNfit 
  //int max1Bin = hSx->GetMaximumBin();
  //hSx_temp->GetXaxis()->SetRangeUser(-20,40);
  // centerBinSx = hSx_temp->GetMinimumBin();//center "dip"
  if(choice == 1) centerBinSx = 1000;
  if(choice == 2) centerBinSx = 1001;
  if(choice == 3) centerBinSx = 1011;
  //std::cout<<"center X Sx "<<hSx->GetBinCenter(centerBinSx)<<std::endl;
  hSx_temp->GetXaxis()->UnZoom();//resume hSx
  shift = centerBinSx - meanBinFitX;
  /// mu_P0 = hNfit->GetBinCenter(meanBinFitX) - hSx->GetBinCenter(centerBinSx);
  //hSx->Scale(1e4/hSx->Integral());
  //hNfit->Scale(1e4/hNfit->Integral());
  std::cout<<"hSx center bin "<<centerBinSx<<" shift "<<shift<<" "<<endl;//leftBin<<" "<<meanBin<<" "<<rightBin<<" "<<centerBin<<" "<<shift<<std::endl;
  int ishift = 0;
  /// shift>0, left shift; else right shift
  ///for (int i=0;i<=binTot;i++) 
  ///{
  ///  hSx->SetBinContent(i+1,hSx_temp->GetBinContent(i+0));
  ///}
  if(shift<0) ishift = -1*shift;//
  else ishift = shift;
  for (int i=0;i<=binTot;i++) 
  {
    if(i<=binTot-shift) hSx->SetBinContent(i,hSx_temp->GetBinContent(i+shift));
    else hSx->SetBinContent(i,0);
  }
  for (int i=0;i<=binTot;i++)
  {
    hNfit->SetBinContent(i,hNfit_temp2->GetBinContent(i));
  }
  // hSx->Draw();hNfit->Draw("same");

//for(int i=0;i<ishift;i++)
//hNfit->SetBinContent(0);
//mu_P0 = mean_hNfit - hSx->GetBinCenter(centerBinSx); 
//  mu_P0 = hNfit->GetBinCenter(centerBin);//after shift, hNfit mean bin is at hSx center
  for(int iter=0;iter<binTot;iter++)
  {NSx[iter]=hSx->GetBinContent(iter+1);}

  for(int iter=0;iter<binTot;iter++)
  {Nx[iter]= hSx->GetXaxis()->GetBinCenter(iter+1);}

  double max1 = hSx->GetMaximum();//GetBinContent(centerBinSx);
  double max2 = hNfit->GetBinContent(meanBinFitX);
  hSx->Scale(1./hSx->Integral());
  //hNfit->Scale(max1/max2);//scale to hSx
  //hNfit->Draw("same");

  hNfit->Rebin(4); hSx->Rebin(4);
  //hSx->GetXaxis()->UnZoom();
  hNfit->GetXaxis()->SetRangeUser(loCut+mean_hNfit,hiCut+mean_hNfit);
  const int Npar = 5;
  double range1 = mean_hNfit - fitRange;
  double range2 = mean_hNfit + fitRange;
  if(range1<-6000) range1 = -6000;
  if(range2>6000) range2 = 6000;

  if(choice == 1 || choice == 2) defaultFit = true;
  if(defaultFit)
  {
    range1 = -rangeDefault;
    range2 = rangeDefault;
  }
  mean_hNfit = mean_hNfit - hSx->GetMean();
  cout<<"fit range: "<<range1<<", "<<range2<<endl;
  TF1 *funMain = new TF1("fitConvl",cal_convl,range1, range2, Npar); 
  // funMain->FixParameter(4,0.45);
  funMain->SetParameters(mu_P0,180, 320,0.45,20);
  funMain->SetParName(0,"#mu_{P}");
  funMain->SetParName(1,"#sigma_{P}");
  funMain->SetParName(2,"#tau_{P}");
  funMain->SetParName(3,"#alpha_{e}");
  funMain->SetParName(4,"scale");
  funMain->SetLineColor(kRed);
  TCanvas *c = new TCanvas("c","",800,600);
  c->cd();
  hNfit->SetMarkerStyle(4);
  hNfit->Fit("fitConvl","rb");
  hNfit->GetXaxis()->SetTitle("x_{fit} [mm]");
  hNfit->GetYaxis()->SetTitle("Counts (40 mm/bin)");
  // hSx->Draw("sames");
  double aE = funMain->GetParameter(4),  aEerr= funMain->GetParError(4);
  double sigma = funMain->GetParameter(0), sigmaerr= funMain->GetParError(0);
  double tau = funMain->GetParameter(2), tauerr= funMain->GetParError(2);
  double mu = funMain->GetParameter(1), muerr= funMain->GetParError(1);
  double ndf = funMain->GetNDF();
  double chi2 = funMain->GetChisquare();
  
  hNfit->SetMarkerStyle(4);   
  std::cout<<setprecision(4);
  std::cout<<" & "<<aE<<"$\\pm$"<<aEerr<<" & "<<sigma<<"$\\pm$"<<sigmaerr<<" & "<<tau<<"$\\pm$"<<tauerr<<" & "<<mu<<"$\\pm$"<<muerr<<" & "<<chi2<<"/"<<ndf<<std::endl;
  std::cout<<aE<<", "<<aEerr<<", "<<mu<<", "<<muerr<<", "<<sigma<<", "<<sigmaerr<<", "<<tau<<", "<<tauerr<<", "<<chi2<<"/"<<ndf<<std::endl;

  TString newfile = "Convolved_"+ filename;
  TFile *fconv = new TFile(newfile,"recreate");
  TH1F *hConvl= new TH1F("hConvl","tempt convl",binTot,lowBin,hiBin);

 // for(int m=0; m<binTot; m++) {
 //   hConvl->SetBinContent(m+1,htempG->GetBinContent(m));
 // }
 // hConvl->Scale(scalepeak/hConvl->GetMaximum());

  fconv->cd();
  hSx->Write();hNfit->Write();
  funMain->Write();
//  hConvlOutput->Write();

}
