//pars 1 = mean gaus
//pars 2 = sig gaus
//pars 3 = +fractional exponential component
//pars 4 = +slope exp
//pars 5 = is a manually added mean as found by fitting a ordanary Gausian
void davidFit()
{
  TFile *f = new TFile("Merged_gamma2_0.ntuple.root");
  TH1F *h = new TH1F("h","",100,0,10);
  TTree *t = (TTree*)f->Get("output");
  t->Project("h","energy","fitValid & scintFit");
  h->Scale(1./h->GetMaximum());
  TF1 *fgausexp = new TF1("fgausexp",fgausexp,1.5,5,6);
  fgausexp->SetParameter(0,1.797);
  fgausexp->SetParameter(1,0.08);
//  fgausexp->SetParameter(4,-2);
  h->Fit(fgausexp,"rq"); 
  //TF1 *g1    = new TF1("g1","gaus",1.5,2.2);
  //TF1 *exp = new TF1("exp","exp",2.2,6);
  //TF1 *total = new TF1("total","gaus(0)+exp(3)",1.5,6);
  //h->Fit(total,"R+");


}

double fgausexp(double *x,double *par)
{
  double gaus = par[2]*TMath::Gaus(x[0],par[0],par[1],true);
  double expo = par[3]*TMath::Exp(par[4]*x[0])+par[5];
  double gausexp = gaus+expo;
  return gausexp;
  
//  if(pars[2]!=0){
//    if(x[0]>=pars[5]){
//      f = pars[0]*(1-pars[3])*exp(-0.5*pow(((x[0]-pars[1])/pars[2]),2))+pars[3]/(2*pars[4])*exp(-TMath::Abs(x[0]-pars[1])/pars[4]);
//    }
//    else{
//      f = pars[0]*exp(-0.5*pow(((x[0]-pars[1])/pars[2]),2));
//    }
//  }
//  return f;
}
