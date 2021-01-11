#include <algorithm>
#include "TGraph.h"
#include "TCanvas.h"

// Use ROOT6
void plotCosChDumpSinusoidal()
{

   // Plot pmt angular distributions based on MPW dump direction 
   //T->Scan("sunDirX:sunDirY:sunDirZ","eventGTID==7882360")
   TVector3 sunDir(-0.538142, 0.7354762, 0.4116766);
   TVector3 evtPos(5124.0131, 4292.3745, -3525.042);
   TVector3 evtDir(-0.601272, -0.755969, 0.2588064);

// Rat water
   TVector3 evtPosRat(3185.7912, 2047.4488, -2613.959);
   TVector3 evtDirRat(-0.428308, -0.897953, 0.1011471);

   TVector3 xplus(1,0,0);
   TVector3 assumeDir;
   assumeDir = xplus;
   //assumeDir = -sunDir;

   TVector3 diff; 
   /// rotation to the sun direction
   double cosPhi = evtDir.Unit()*assumeDir.Unit();
   double phiRot = acos(cosPhi);
   double sinPhi = sin(phiRot);
   TVector3 rotN = evtDir.Cross(assumeDir);
   double nx = rotN.X(), ny = rotN.Y(), nz = rotN.Z();
   TMatrixD Rot(3,3);
   Rot[0][0] = nx*nx*(1-cosPhi)+cosPhi; Rot[0][1] = nx*ny*(1-cosPhi)-nz*sinPhi; Rot[0][2] = nx*nz*(1-cosPhi)+ny*sinPhi; 
   Rot[1][0] = nx*ny*(1-cosPhi)+nz*sinPhi; Rot[1][1] = ny*ny*(1-cosPhi)+cosPhi; Rot[1][2] = ny*nz*(1-cosPhi)-nx*sinPhi; 
   Rot[2][0] = nx*nz*(1-cosPhi)-ny*sinPhi; Rot[2][1] = ny*nz*(1-cosPhi)+nx*sinPhi; Rot[2][2] = nz*nz*(1-cosPhi)+cosPhi; 
   //cout<<R[0][0]<<" "<<R[0][1]<<" "<<R[0][2]<<endl;
   //cout<<R[1][0]<<" "<<R[1][1]<<" "<<R[1][2]<<endl;
   //cout<<R[2][0]<<" "<<R[2][1]<<" "<<R[2][2]<<endl;


   TVector3 dirRot = Rot*evtDir;
   TVector3 evtRot = Rot*evtPos;
   TMatrixD MevtPos(3,1);
   MevtPos[0][0] = evtPos.X(), MevtPos[1][0] = evtPos.Y(), MevtPos[2][0] = evtPos.Z();
   TVector3 pmtpos;
   double rpsup = 8390;

   double band = 200; //20 cm
   double rad = band/rpsup;

   //plot Cherenkov cone band
   //
   //
   const int nbin = 1000;
   double thetaCh = 41./180*TMath::Pi();

   double ycone = rpsup*cos(thetaCh);
   double xcone[nbin], zcone[nbin];

   double conePhiMin[nbin], coneThetaMin[nbin];
   double conePhiMax[nbin], coneThetaMax[nbin];

   TGraph *grmin = new TGraph(nbin,conePhiMin, coneThetaMin);
   TGraph *grmax = new TGraph(nbin,conePhiMax, coneThetaMax);
   TGraph *grshade = new TGraph(2*nbin);
   TVector3 coneTemp;

   thetaCh = thetaCh - rad;

   double step = 2*TMath::Pi()/double(nbin);
   for (int i=0;i<nbin;i++) {
    zcone[i] = rpsup*sin(thetaCh)*sin(0+i*step);
    xcone[i] = rpsup*sin(thetaCh)*cos(0+i*step);
    // cout<<xcone[i]<<" "<<zcone[i]<<" "<<(sqrt(xcone[i]*xcone[i]+zcone[i]*zcone[i])/rpsup+band/rpsup)*180/3.14<<endl;
    coneTemp.SetXYZ(xcone[i],ycone,zcone[i]);
    //conePhiMin[i] = coneTemp.Phi();
    //coneThetaMin[i] = coneTemp.CosTheta();
    conePhiMin[i] = coneTemp.Phi()*sin(coneTemp.Theta());//sinu
    coneThetaMin[i] = coneTemp.Theta();//sinu
    //cout<<coneTemp.Phi()<<" "<<coneTemp.CosTheta()<<endl;
   }

   thetaCh = thetaCh + rad;
   for (int i=0;i<nbin;i++) {
    zcone[i] = rpsup*sin(thetaCh)*sin(0+i*step);
    xcone[i] = rpsup*sin(thetaCh)*cos(0+i*step);
    coneTemp.SetXYZ(xcone[i],ycone,zcone[i]);
    conePhiMax[i] = coneTemp.Phi()*sin(coneTemp.Theta());//sinu
    coneThetaMax[i] = coneTemp.Theta();//sinu
   }
 
   for (int i=0;i<nbin;i++) {
     grshade->SetPoint(i,conePhiMax[i], coneThetaMax[i]);
     grshade->SetPoint(nbin +i,conePhiMin[nbin-i-1], coneThetaMin[nbin-i-1]);
   }
   grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
//   grmin->Draw("al");
//   grmax->Draw("l");
//   grshade->Draw("f");

   
   //gr->SetLineWidth(4);
//gr->SetMarkerColor(4);
//gr->SetMarkerStyle(21);
//gr->Draw("CP");
//

   TFile *f0 = new TFile("allPMTs_100915_s001_gtID_169700_waterdirection.root");
   TFile *f = new TFile("mpw_100915_s001_gtID_169700_waterdirection.root");

   TH2F *hpmt = (TH2F*)f->Get("pmtGEO");

   TH2F *hpmtBeforeSelection = (TH2F*)f0->Get("pmtGEO");

   int totx = hpmt->GetXaxis()->GetNbins();
   int toty = hpmt->GetYaxis()->GetNbins();
    
   double theta = 0;
   double phi = 0;
   double stepx = (TMath::Pi()+TMath::Pi())/totx;
   double stepy = 2./toty;
 
   double cosCh = 0;
   TH1F *hCh = new TH1F("hCh","",200,-1,1);
   
   //sinusodial map
   int fNBinsPhi = 45, fNBinsTheta = 45;
   TH2F* hpmtPhiTheta = new TH2F("hpmtPhiTheta", "PMT distributions, sinusoidal ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi());
   TH2F* hpmtPhiThetaRot = new TH2F("hpmtPhiThetaRot", "after rotate PMT distributions, sinusoidal ", fNBinsPhi, -TMath::Pi(), TMath::Pi(), fNBinsTheta, 0, TMath::Pi());
   //// sinusodial map
   //// pmtpos.Phi() * sin(pmtpos.Theta()), pmtpos.Theta()

   //TH2F *hpmtPhiTheta = new TH2F("hpmtPhiTheta","original", 628, -TMath::Pi(), TMath::Pi(), 200, -1, 1);
   //TH2F *hpmtPhiThetaRot = new TH2F("hpmtPhiThetaRot","after rotate", 628, -TMath::Pi(), TMath::Pi(), 200, -1, 1);
   int npoints = 0;
   vector<double> vPhi, vTheta;
   vector<double> vRotPhi, vRotTheta;

   for(int i = 0;i<totx;i++) {
     for(int j = 0;j<toty;j++) {
         double costheta = -1 + j*stepy;
	 theta = acos(costheta);
         phi = -TMath::Pi() + stepx*i;
	 if(hpmt->GetBinContent(i+1,j+1)!=0)
	 {
           //cout<<theta<<" "<<phi<<endl;	 
           //hpmtPhiTheta->Fill(phi, costheta);
	   hpmtPhiTheta->Fill(phi*sin(theta),theta);
           pmtpos.SetMagThetaPhi(rpsup, theta, phi);
	   //vPhi.push_back(phi); vTheta.push_back(costheta);
           vPhi.push_back(phi*sin(theta)); vTheta.push_back(theta);
	   TMatrixD Mpmtpos(3,1);
           Mpmtpos[0][0] = pmtpos.X(); Mpmtpos[1][0] = pmtpos.Y();Mpmtpos[2][0] = pmtpos.Z();
           TMatrixD MpmtposRot(3,1);
	   MpmtposRot = Rot*Mpmtpos - Rot*MevtPos;
	   //hpmtPhiThetaRot->Fill(pmtposRot.Phi(), pmtposRot.CosTheta());
	   TVector3 pmtposRot;
	   pmtposRot.SetXYZ(MpmtposRot[0][0],MpmtposRot[1][0],MpmtposRot[2][0]);
	   hpmtPhiThetaRot->Fill(pmtposRot.Phi()*sin(pmtposRot.Theta()), pmtposRot.Theta());
           // cout<<pmtposRot.Phi()<<" "<<sin(pmtposRot.Theta())<<" "<<pmtposRot.Theta()<<endl;
	   //vRotPhi.push_back(pmtposRot.Phi()); vRotTheta.push_back(pmtposRot.CosTheta()); 
           vRotPhi.push_back(pmtposRot.Phi()*sin(pmtposRot.Theta()));vRotTheta.push_back(pmtposRot.Theta());
	   diff = (pmtpos - evtPos).Unit();
           cosCh = diff*evtDir;
           hCh->Fill(cosCh);
	   npoints++;
	 }
     }
   }

   vector<double> vPhiAll, vThetaAll;
   int nAll = 0;
   for(int i = 0;i<totx;i++) {
     for(int j = 0;j<toty;j++) {
         double costheta = -1 + j*stepy;
         theta = acos(costheta);
         phi = -TMath::Pi() + stepx*i;
         if(hpmtBeforeSelection->GetBinContent(i+1,j+1)!=0)
         {
           //vPhiAll.push_back(phi); vThetaAll.push_back(costheta);
	   vPhiAll.push_back(phi*sin(theta));vThetaAll.push_back(theta);
           nAll++;
	 }
     }
   }

   TGraph *gpmtpos = new TGraph(npoints);
   TGraph *gpmtposAll = new TGraph(nAll-npoints);

   TGraph *gpmtposRot = new TGraph(npoints);


   TGraph *gSunDir = new TGraph(1); // Direction Marker
   TGraph *gDir = new TGraph(1); // Direction Marker
   TGraph *gDirMPW = new TGraph(1); // Direction Marker

   TGraph *gPos = new TGraph(1); // Direction Marker
   TGraph *gPosMPW = new TGraph(1); // Direction Marker

   for (int i=0;i<npoints;i++) {
     gpmtpos->SetPoint(i,vPhi[i], vTheta[i]);
     gpmtposRot->SetPoint(i,vRotPhi[i], vRotTheta[i]);
   }
  
   vector<double> cutPMTphi, cutPMTtheta;
   //not need to sort since it already sorted
   std::set_difference(vPhiAll.begin(), vPhiAll.end(), vPhi.begin(), vPhi.end(),std::inserter(cutPMTphi, cutPMTphi.begin()));
   std::set_difference(vThetaAll.begin(), vThetaAll.end(), vTheta.begin(), vTheta.end(),std::inserter(cutPMTtheta, cutPMTtheta.begin()));

   if(nAll - npoints == 0) cout<<"no pmt cuts"<<endl;
   if(nAll - npoints != 0) {
     for (int i=0;i<nAll-npoints;i++) {
       gpmtposAll->SetPoint(i,cutPMTphi[i], cutPMTtheta[i]);
     }
   }

   gSunDir->SetMarkerStyle(43);gSunDir->SetMarkerSize(3);
   TVector3 sunPoint = -sunDir;
   //gSunDir->SetPoint(0,sunPoint.Phi(), sunPoint.CosTheta());
   gSunDir->SetPoint(0,sunPoint.Phi()*sin(sunPoint.Theta()), sunPoint.Theta());

   gpmtposAll->SetMarkerStyle(5);gpmtposAll->SetMarkerSize(3);

   gDir->SetMarkerStyle(29);gDir->SetMarkerSize(3);
   gDir->SetMarkerColor(kBlue);
   //gDir->SetPoint(0,evtDirRat.Phi(),evtDirRat.CosTheta());
   gDir->SetPoint(0,evtDirRat.Phi()*sin(evtDirRat.Theta()),evtDirRat.Theta());

   gPos->SetMarkerStyle(21);gPos->SetMarkerSize(2);
   gPos->SetMarkerColor(kBlue);
   //gPos->SetPoint(0,evtPosRat.Phi(),evtPosRat.CosTheta());
   gPos->SetPoint(0,evtPosRat.Phi()*sin(evtPosRat.Theta()),evtPosRat.Theta());

   gDirMPW->SetMarkerStyle(30);gDirMPW->SetMarkerSize(3);
   gDirMPW->SetMarkerColor(kRed);
   //gDirMPW->SetPoint(0,evtDir.Phi(),evtDir.CosTheta());
   gDirMPW->SetPoint(0,evtDir.Phi()*sin(evtDir.Theta()),evtDir.Theta());

   gPosMPW->SetMarkerStyle(25);gPosMPW->SetMarkerSize(2);
   gPosMPW->SetMarkerColor(kRed);
   //gPosMPW->SetPoint(0,evtPos.Phi(),evtPos.CosTheta());
   gPosMPW->SetPoint(0,evtPos.Phi()*sin(evtPos.Theta()),evtPos.Theta());

   gpmtpos->SetMarkerStyle(4);gpmtpos->SetMarkerSize(2);
//   gpmtpos->GetXaxis()->SetTitle("#phi");
//   gpmtpos->GetYaxis()->SetTitle("cos#theta");
   gpmtpos->GetXaxis()->SetTitle("#phi sin#theta");
   gpmtpos->GetYaxis()->SetTitle("#theta");

   gpmtpos->GetXaxis()->SetTitleOffset(0.7);
   gpmtpos->GetXaxis()->SetTitleSize(0.05);
   gpmtpos->GetYaxis()->SetTitleOffset(0.7);
   gpmtpos->GetYaxis()->SetTitleSize(0.05);

   TCanvas *c1 = new TCanvas("c1","",800,600);
   c1->cd();
   auto *legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(gpmtpos,"PMTs","p");
   legend->AddEntry(gSunDir,"-1*sunDir","p");
   legend->AddEntry(gDir,"eventDir, waterFitter","p");
   legend->AddEntry(gDirMPW,"eventDir, MPW","p");
   legend->AddEntry(gPos,"eventPos, waterFitter","p");
   legend->AddEntry(gPosMPW,"eventPos, MPW","p");
   gpmtpos->Draw("AP");
   gDir->Draw("P");
   gDirMPW->Draw("P");
   gPos->Draw("P");
   gPosMPW->Draw("P");
   gSunDir->Draw("P");
   if(nAll - npoints != 0) 
   { gpmtposAll->Draw("P");
     legend->AddEntry("gpmtposAll","PMTs removed by selectors","p");
   }
   legend->Draw("sames");
   c1->Modified();
   c1->Update();

//   //hCh->Draw();
   TCanvas *c3 = new TCanvas("c3","",800,600);
   c3->Divide(2,1);
   c3->cd(1);
   hpmt->Draw();
   hpmtPhiTheta->Draw();
   //gpmtpos->Draw("AP");
   //grmin->Draw("l");
   //grmax->Draw("l");
   grshade->Draw("f");

   c3->cd(2);
   hpmtPhiThetaRot->Draw();


}
