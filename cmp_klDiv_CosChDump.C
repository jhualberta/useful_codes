//#include <algorithm>
//#include "TGraph.h"
//#include "TCanvas.h"

//!!!!! Use ROOT6
// calculate KL-div
{
//100207
   // Plot pmt angular distributions based on MPW dump direction 
   //T->Scan("sunDirX:sunDirY:sunDirZ","eventGTID==7882360")
   TVector3 sunDir(-0.538142, 0.7354762, 0.4116766);
   TVector3 evtPos(5124.0131, 4292.3745, -3525.042);
   TVector3 evtDir(-0.601272, -0.755969, 0.2588064);

// Rat water
   TVector3 evtPosRat(3185.7912, 2047.4488, -2613.959);
   TVector3 evtDirRat(-0.428308, -0.897953, 0.1011471);

   TFile *f0 = new TFile("allPMTs_100207_s001_gtID_5079885_waterdirection.root");
   TFile *f = new TFile("mpw_100207_s001_gtID_5079885_waterdirection.root");

   TH2F *hpmt = (TH2F*)f->Get("pmtGEO");
   TH2F *hpmtBeforeSelection = (TH2F*)f0->Get("pmtGEO");

   int totx = hpmt->GetXaxis()->GetNbins();
   int toty = hpmt->GetYaxis()->GetNbins();
   
   double rpsup = 8390;
   double band = 200; //20 cm
   double rad = band/rpsup;

   double theta = 0;
   double phi = 0;
   double stepx = (TMath::Pi()+TMath::Pi())/totx;
   double stepy = 2./toty;
 
   double cosCh = 0;
   TH1F *hCh = new TH1F("hCh","",200,-1,1);
   TH1F *hChRat = new TH1F("hChRat","",200,-1,1);

   TVector3 pmtpos, diff;
   TH2F *hpmtPhiTheta = new TH2F("hpmtPhiTheta","original", 628, -TMath::Pi(), TMath::Pi(), 200, -1, 1);
   TH2F *hpmtPhiThetaRot = new TH2F("hpmtPhiThetaRot","after rotate", 628, -TMath::Pi(), TMath::Pi(), 200, -1, 1);
   int npoints = 0;
   vector<double> vPhi, vTheta;
   vector<double> vRotPhi, vRotTheta;
   double noiseAve = 0, countNoise = 0;
   for(int i = 0;i<totx;i++) {
     for(int j = 0;j<toty;j++) {
         double costheta = -1 + j*stepy;
	 theta = acos(costheta);
         phi = -TMath::Pi() + stepx*i;
	 if(hpmt->GetBinContent(i+1,j+1)!=0)
	 {
           //cout<<theta<<" "<<phi<<endl;	 
           hpmtPhiTheta->Fill(phi, costheta);
           pmtpos.SetMagThetaPhi(rpsup, theta, phi);
	   vPhi.push_back(phi); vTheta.push_back(costheta);
           TMatrixD Mpmtpos(3,1);
           Mpmtpos[0][0] = pmtpos.X(); Mpmtpos[1][0] = pmtpos.Y();Mpmtpos[2][0] = pmtpos.Z();
           TMatrixD MpmtposRot(3,1);
	   diff = (pmtpos - evtPos).Unit();
           cosCh = diff*evtDir;
           hCh->Fill(cosCh);
	   hChRat->Fill((pmtpos-evtPosRat).Unit()*evtDirRat);
	   npoints++;
	 }
     }
   }
  noiseAve = noiseAve/countNoise;
  //water angular pdf, 200, -1,1
  double angularPDF[200]={0.0140342, 0.0141184, 0.0143289, 0.0145342, 0.0148263, 0.0151684, 0.0152658, 0.0150974, 0.0156105, 0.0155079, 0.0157974, 0.0162526, 0.0168789, 0.0167263, 0.0168289, 0.0175921, 0.0176895, 0.0181237, 0.0182763, 0.0186289, 0.0186, 0.0193868, 0.0194711, 0.0200395, 0.0198184, 0.0208079, 0.0209, 0.0214684, 0.0218921, 0.0224632, 0.0230737, 0.0233053, 0.0238079, 0.0241263, 0.0246211, 0.0248105, 0.0259684, 0.0259658, 0.0266605, 0.0274895, 0.0280605, 0.0284184, 0.0289421, 0.0294237, 0.0300816, 0.0305421, 0.0313974, 0.0316184, 0.0328868, 0.0332158, 0.0333158, 0.0349079, 0.0346316, 0.0361368, 0.0371474, 0.0369553, 0.0383737, 0.0394105, 0.0402605, 0.0404658, 0.0409132, 0.0419579, 0.0429605, 0.0440158, 0.0450132, 0.0461395, 0.0470763, 0.0479711, 0.0487868, 0.0500579, 0.0512711, 0.0515763, 0.0533842, 0.0542184, 0.0553368, 0.0565789, 0.0570026, 0.0590368, 0.0604763, 0.0614895, 0.0627368, 0.0642579, 0.0665026, 0.0664237, 0.0681132, 0.0696211, 0.0720237, 0.0737974, 0.0745342, 0.0763921, 0.0784211, 0.0802868, 0.0812921, 0.0829921, 0.0848526, 0.0870789, 0.0892105, 0.0913632, 0.0935263, 0.0953184, 0.0971263, 0.0994947, 0.101568, 0.10505, 0.106989, 0.108747, 0.112716, 0.115024, 0.117761, 0.120861, 0.122789, 0.126653, 0.129711, 0.133795, 0.137197, 0.140158, 0.143882, 0.147392, 0.150787, 0.153826, 0.1595, 0.162261, 0.166605, 0.171576, 0.174908, 0.180482, 0.186182, 0.191463, 0.196187, 0.201682, 0.206774, 0.213574, 0.218197, 0.224616, 0.231679, 0.239163, 0.245353, 0.253068, 0.260413, 0.268808, 0.275805, 0.284239, 0.293247, 0.302471, 0.310608, 0.321092, 0.332032, 0.344274, 0.355942, 0.367397, 0.380934, 0.392561, 0.407068, 0.421218, 0.438945, 0.454755, 0.473458, 0.490255, 0.508824, 0.529003, 0.552076, 0.572608, 0.599029, 0.625458, 0.654121, 0.684313, 0.717371, 0.74935, 0.792632, 0.829653, 0.877292, 0.92525, 0.970547, 0.996903, 1.00237, 0.987555, 0.947292, 0.907053, 0.869334, 0.841003, 0.809482, 0.776424, 0.750021, 0.724474, 0.700124, 0.674739, 0.650079, 0.629787, 0.611989, 0.589055, 0.569929, 0.550379, 0.532632, 0.515161, 0.496537, 0.482984, 0.466082, 0.450729, 0.437766, 0.425111}; 
  // calculate KL-div

  cout<<hCh->Integral()<<" "<<hChRat->Integral()<<endl;

//  hCh->Scale(1./hCh->Integral());
//  hChRat->Scale(1./hChRat->Integral());
  TH1F *hpdf = new TH1F("hpdf","",200,-1,1);
  for(int i = 0;i<200;i++)
  {
   double q = angularPDF[i];
   hpdf->SetBinContent(i+1,q);
  }

  double lowVal = 0;
  int countlowVal = 0;
  for(int i = 0;i<200;i++)
  {
    double val = -1 + 0.01*i;
    double y = hCh->GetBinContent(i);
    if(val<0.4 && y!=0 )
    {	    
      lowVal += y;
      countlowVal++;
    }  
  }
  cout<<"scale "<<countlowVal<<" "<<lowVal<<" "<<lowVal/countlowVal<<endl;
  double scale = lowVal/countlowVal;
  hpdf->Scale(scale/angularPDF[199]);
  //hpdf->Scale(hCh->Integral()/hpdf->Integral());
  double kl1 = 0;
  double kl2 = 0;
  for(int i = 0;i<200;i++)
  {
    double p1 = hCh->GetBinContent(i+1);
    double p2 = hChRat->GetBinContent(i+1);
    double q = hpdf->GetBinContent(i+1);
    if(p1!=0)
    {
      kl1 = kl1+p1*log(p1/q);  
    }
  
    if(p2!=0)
    {
      kl2 = kl2+p2*log(p2/q);
    }
  }
  cout<<"kL-div: MPW= "<<kl1<<" rat= "<<kl2<<endl;
  hpdf->SetLineColor(kGreen+2);
TH1::SetDefaultSumw2();
  TFile *ff = new TFile("save_compare_ChAngle.root","recreate");
   hCh->SetLineColor(kRed);
   hCh->Draw();
   hpdf->Draw("sames");
   hChRat->Draw("sames");

}
