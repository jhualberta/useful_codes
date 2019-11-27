double asymGaus(double *x, double *par)
{
  //0: mu, 1:sigma, 2: shaping, 3: scale
  double arg = -par[2]*(x[0]-par[0])/(sqrt(2)*par[1]);// erfc(arg)
  double asym = par[3]*TMath::Gaus(x[0],par[0],par[1],true)*TMath::Erfc(arg);
  return asym;
}

double expo(double *x, double *par)
{
  //0:scale, 1:slope, 2:shift to y 
  double exp = par[1]*TMath::Exp(par[0]*x[0])+par[2];
  return exp;
}

double gausexp(double *x, double *par)
{
  double gaus = par[3]*TMath::Gaus(x[0],par[0],par[1],true)*TMath::Erfc(-par[2]*(x[0]-par[0])/(sqrt(2)*par[1]));
  double expo = par[5]*TMath::Exp(par[4]*x[0])+par[6];
  double gausexp = gaus+expo;
  return par[7]*gausexp+par[8];
}

void fitAsymGaus()
{
  TFile *f1 = new TFile("Merged_SaveAmbe_r252665to252684.root");
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiFit");
  TTree *gout = (TTree*)f1->Get("sno2p2");
  TH1F *hNhitDely = new TH1F("hNhitDely","delay",2500,0,2500);
  TH1F *hNhitPrompt = new TH1F("hNhitPrompt","delay",2500,0,2500);
  gout->Draw("DNhits>>hNhitDely");
  gout->Draw("PNhits>>hNhitPrompt");
  double dRange1 = 650, dRange2 = 850;
  double PRange1 = 1100, PRange2 = 1500;

  hNhitDely->Sumw2();
  hNhitDely->Rebin(2);
  //hNhitDely->Sumw2();
  hNhitDely->GetXaxis()->SetRangeUser(dRange1, dRange2);

  hNhitPrompt->Sumw2();
  hNhitPrompt->Rebin(2);
  //hNhitPrompt->Sumw2();
  hNhitPrompt->GetXaxis()->SetRangeUser(PRange1, PRange2);

  TF1 *fasymgaus = new TF1("fasymgaus",asymGaus,dRange1,dRange2,4);
  double scale = hNhitDely->GetMaximum();
  fasymgaus->SetParameters(hNhitDely->GetXaxis()->GetBinCenter(hNhitDely->GetMaximumBin()),100,1,scale);
  fasymgaus->SetParNames("#mu","#sigma","scale1","scale2");
  hNhitDely->Fit(fasymgaus,"QR");

  TF1 *fasymgaus1 = new TF1("fasymgaus1",asymGaus,PRange1,PRange2,4);
  double scale2 = hNhitPrompt->GetMaximum();
  fasymgaus1->SetParameters(hNhitPrompt->GetXaxis()->GetBinCenter(hNhitPrompt->GetMaximumBin()),100,1,scale);
  fasymgaus1->SetParNames("#mu","#sigma","scale1","scale2");
  hNhitPrompt->Fit(fasymgaus1,"QR");

}
