#include "TH1D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TPad.h"
#include <iomanip>
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "Math/MinimizerOptions.h"
#include "TTree.h"
#include "TLegend.h"
double angularResol(double *x, double *par)
{
    //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale
	double jillings = par[0]*par[1]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*par[2]*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    return jillings*par[3];	
}

/// for calculating cosTheta
double inte_angularResol(double *x, double *par)//an integral of the angularResol, to calculate cosTheta
{
    //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale
	double inte_jillings = par[0]*par[1]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*par[2]*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    return inte_jillings*par[3];	
}
double fSolve(double *x, double *par)//an integral of the angularResol, to calculate cosTheta
{
   //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale, par[4]=F1, par[5]=percent1, par[6] = integral
    double solve_jillings = par[3]*par[0]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    double result = solve_jillings-par[4]+par[5]*par[6];
    return result;	
}

double boulayAngularResol(double *x, double *par)
{
    //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale

	double jillings = par[0]*par[1]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*par[2]*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    return jillings*par[3];	
}

void SaveAngular6176_simpleBeta(TString infile_name)//, double xsrc, double ysrc, double zsrc)
{
  double ae,sig_p,mu_p,tau_p;//par[0],par[1],par[2],par[3]
  int bin = 200;
  //double range1 = -0.5;
  double range1 = 0.3;//Ed Leming range
  //TVector3 srcPos(xsrc,ysrc,zsrc);
  TFile *ff2 = new TFile(infile_name,"read");
  TTree *Tdata = (TTree*)ff2->Get("T");

  double posx, posy, posz, energy, posRad, itr, dirx, diry, dirz, beta14, scaleLogL, mcDirx, mcDiry, mcDirz;
  UInt_t nhits;

  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("posRad", &posRad);
  Tdata->SetBranchAddress("itr", &itr);
  Tdata->SetBranchAddress("dirx", &dirx);
  Tdata->SetBranchAddress("diry", &diry);
  Tdata->SetBranchAddress("dirz", &dirz);

  Tdata->SetBranchAddress("mcDirx", &mcDirx);
  Tdata->SetBranchAddress("mcDiry", &mcDiry);
  Tdata->SetBranchAddress("mcDirz", &mcDirz);

  Tdata->SetBranchAddress("beta14", &beta14);
  Tdata->SetBranchAddress("energy", &energy);
  Tdata->SetBranchAddress("scaleLogL", &scaleLogL);

  TH1F *hAngularCut0 = new TH1F("hAngularCut0","",bin,-1,1);
  TH1F *hAngularCut1 = new TH1F("hAngularCut1","",bin,-1,1);
  TH1F *hAngularCut2 = new TH1F("hAngularCut2","",bin,-1,1);
  TH1F *hAngularCut3 = new TH1F("hAngularCut3","",bin,-1,1);

  for(int i =0;i<Tdata->GetEntries();i++)
  {
    Tdata->GetEntry(i);
    //cout<<energy<<" "<<nhits<<endl;
    double posRad = sqrt(posx*posx+posy*posy+(posz-108)*(posz-108));
    // if(nhits>5 && itr>0.55 && beta14<0.95 && beta14>-0.12 && energy>3.5 && posRad<5850 && scaleLogL>10)
    //if(nhits>5 && itr>0.55 && beta14<0.95 && beta14>-0.12 && energy>3.5)
    {
      TVector3 evtPos(posx,posy,posz-108);//source position is corrected for posz-108
      TVector3 u_fit(dirx,diry,dirz);
      TVector3 u_mc(mcDirx, mcDiry, mcDirz);
      //TVector3 Xdiff = evtPos-srcPos;
      //double diff = Xdiff.Mag();
      // if(diff>1000 && diff<2300) 
      hAngularCut0->Fill(u_fit*u_mc);
      //if(diff>1000) hAngularCut1->Fill(Xdiff.Unit()*Dir);
      //if(diff>1200 && diff<2300) hAngularCut2->Fill(Xdiff.Unit()*Dir);
      //if(diff>1500) hAngularCut3->Fill(Xdiff.Unit()*Dir);
    }
  }

  TString runID(infile_name(infile_name.Index("10"),6));
  TFile *f1 = new TFile(infile_name);    
  //From SNO+ data
  double inteCount = hAngularCut0->Integral(range1,1);

  //while(inteCount>5000)
  //{	  
  //  if(inteCount<5000) 
  //  {
  //    range1 = range1 - 0.1;	    
  //    inteCount = hAngularCut0->Integral(range1,1);
  //  }
  //  if(inteCount>5000 || range1<-0.5) break;
  //}
  TF1 *funcAngularResol = new TF1("fitAngularResol",angularResol,range1,1.0,4);
///fit Angular Resolution 
  // funcAngularResol->FixParameter(0,0.5);
  funcAngularResol->SetParameters(0.5,4.0,18.8,200);
  /// funcAngularResol->SetParLimits(0,0.1,0.9);
  funcAngularResol->SetParNames("#alpha_{M}","#beta_{M}","#beta_{S}","scale");
///Rebin histograms
   //hAngular->Rebin(10.0);hAngularCut1->Rebin(10.0);hAngularCut2->Rebin(10.0);hAngularCut3->Rebin(10.0);hAngularCut4->Rebin(10.0);hAngularCut5->Rebin(10.0);    
   //hAngular->Sumw2();hAngularCut1->Sumw2();hAngularCut2->Sumw2();hAngularCut3->Sumw2();hAngularCut4->Sumw2();hAngularCut5->Sumw2();
     
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiFit");
   hAngularCut0->Fit(funcAngularResol,"PRq");
   //   hAngularCut2->Fit(funcAngularResol,"PRq");//q for quiet
  
   TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(hAngularCut0,"data","lep");
   //if(choice == 1) legend->AddEntry(hAngularCut1,"data","lep");
   legend->AddEntry(funcAngularResol,"resolution fit","l");
   legend->Draw();
 
   double ndf = funcAngularResol->GetNDF();
   double chi2 = funcAngularResol->GetChisquare();    
  
   TF1 *inte_funcAngularResol = new TF1("inte_fitAngularResol",inte_angularResol,-1.0,1.0,4);
   //fitting parameters for Integral: par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale
   double aM1, bM1, bS1, scale1;
   aM1 = funcAngularResol->GetParameter(0);
   bM1 = funcAngularResol->GetParameter(1);
   bS1 = funcAngularResol->GetParameter(2);
   scale1 = funcAngularResol->GetParameter(3);
   inte_funcAngularResol->SetParameters(aM1,bM1,bS1,scale1);
  //check the integral function is right
  // cout<<inte_funcAngularResol->Eval(1.0)-inte_funcAngularResol->Eval(-1.0)<<" "<<funcAngularResol->Integral(-1.0,1.0)<<endl;
  
   double integral =  inte_funcAngularResol->Eval(1.0)-inte_funcAngularResol->Eval(-1.0); 
   double percent1 = 0.5;double percent2 = 0.8;double percent3 = 0.95;
   double F1 = inte_funcAngularResol->Eval(1.0);
   //Solve function: F(x)-F(1) + percent1*integral = 0
   //par[0]=aM, par[1]=bM, par[2]=bS,apar[3]=scale, par[4]=F1, par[5]=percent1, par[6] = integral
   TF1 SolvFunc("solve_func",fSolve,-1.0,1.0,7);
   SolvFunc.SetParameters(aM1,bM1,bS1,scale1,F1,percent1,integral);
   ROOT::Math::WrappedTF1 wf1(SolvFunc); 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf;
   // Set parameters of the method
   brf.SetFunction( wf1, -1.0, 1.0);
   brf.Solve();
   double cosTheta50 = brf.Root();
  // cout << "cosTheta50% "<<brf.Root() << endl;
   
   SolvFunc.SetParameters(aM1,bM1,bS1,scale1,F1,percent2,integral);
   ROOT::Math::WrappedTF1 wf2(SolvFunc); 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf2;
   // Set parameters of the method
   brf2.SetFunction( wf2, -1.0, 1.0);
   brf2.Solve();
   double cosTheta80 = brf2.Root();
   //cout << "cosTheta80% "<<brf2.Root() << endl;

   SolvFunc.SetParameters(aM1,bM1,bS1,scale1,F1,percent3,integral);
   ROOT::Math::WrappedTF1 wf3(SolvFunc); 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf3;
   // Set parameters of the method
   brf3.SetFunction( wf3, -1.0, 1.0);
   brf3.Solve();
   double cosTheta95 = brf3.Root();
  // cout << "cosTheta95% "<<brf3.Root() << endl;
   TString savefilename = "SaveAngle_"+infile_name;
   TFile *savefile = new TFile(savefilename,"recreate");

   double aM2 = funcAngularResol->GetParameter(0), aM2err= funcAngularResol->GetParError(0);
   double bM2 = funcAngularResol->GetParameter(1), bM2err= funcAngularResol->GetParError(1);
   double bS2 = funcAngularResol->GetParameter(2), bS2err= funcAngularResol->GetParError(2);
   cout<<setprecision(4);
   //std::cout<<runID<<" & "<<bM2<<"$\\pm$"<<bM2err<<" & "<<bS2<<"$\\pm$"<<bS2err<<" & "<<aM2<<"$\\pm$"<<aM2err<<" & "<<chi2<<"/"<<ndf<<" & "<<cosTheta50<<" & "<<cosTheta80<<" & "<<cosTheta95<<std::endl;
   // std::cout<<infile_name<<", "
   std::cout<<"range1 "<<range1<<","<<bM2<<","<<bM2err<<","<<bS2<<","<<bS2err<<","<<aM2<<","<<aM2err<<","<<chi2<<"/"<<ndf<<","<<chi2/ndf<<std::endl;
   savefile->cd();
   hAngularCut0->Write();
   //if(choice == 1) hAngularCut1->Write();
   funcAngularResol->Write();

///The Exact Boulay Function    
   //TF1 *funcBoulayPosResol = new TF1("fBoulayPosResol",boulayPosResol,-200,200,5);
   //funcBoulayPosResol->SetParameters(0.55,17.23,-0.5414,32.25,1.0e5);
   //funcBoulayPosResol->SetParNames("#alpha_{e}","#sigma_{P}","#mu_{P}","#tau_{P}","scale");
   //cout<<funcBoulayPosResol->Integral(-200,200)<<endl;
   
   //TF1 *funcBoulayAngularResol = new TF1("fBoulayAngularResol",boulayAngularResol,-0.5,1.0,4);
   //funcBoulayAngularResol->SetParameters(0.5,4.0,18.8,200);
   //cout<<funcBoulayAngularResol->Integral(-1,1)<<endl;
   //funcBoulayAngularResol->SetParNames("#alpha_{M}","#beta_{M}","#beta_{S}","scale");
   //gPad->Update();
/// 
   //TPaveStats *st = (TPaveStats*)hpx->FindObject("stats");
   //st->SetOptStat(1111);
   //hpx->SetStats(1).
   
   //gStyle->SetOptStat(1111);
//   gStyle->SetOptFit(1);
   //hpx->SetStats(1);
   //hpx->Write();
   //hConvol->Write();
   //hAngular->SetStats(1);
 //  hAngular->Write();hAngularCut1->Write();hAngularCut2->Write();hAngularCut3->Write();hAngularCut4->Write();     
 //  funcBoulayPosResol->Write(); funcBoulayAngularResol->Write();

}
