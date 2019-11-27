//pars 1 = mean gaus
//pars 2 = sig gaus
//pars 3 = +fractional exponential component
//pars 4 = +slope exp
//pars 5 = is a manually added mean as found by fitting a ordanary Gausian

double t1 = -7.19, t2 = -24.81, t3 = -269.87;
double a1 = 0.553, a2 = 0.331, a3 = 0.116;
double rt = 0.8;

void timingFit()
{
  TFile *f = new TFile("ResolLeon_Analysis10_r0000250815_s000_p001.root");
  TH1F *htResInScint_nhit100 = (TH1F*)f->Get("htResInScint_nhit100");
    htResInScint_nhit100->Sumw2();
  htResInScint_nhit100->Scale(1./htResInScint_nhit100->Integral());
  TF1 *fLabppo = new TF1("fLabppo",fLabppo,3.5,50,6);
  fLabppo->SetParameter(0,t1);
  fLabppo->SetParameter(1,t2);
  fLabppo->SetParameter(2,t3);
  fLabppo->SetParameter(3,a1);
  fLabppo->SetParameter(4,a2);
  fLabppo->SetParameter(5,a3);

//  fLabppo->SetParameter(4,-2);
  htResInScint_nhit100->Fit(fLabppo,"rq"); 
  //TF1 *g1    = new TF1("g1","gaus",1.5,2.2);
  //TF1 *exp = new TF1("exp","exp",2.2,6);
  //TF1 *total = new TF1("total","gaus(0)+exp(3)",1.5,6);
  //h->Fit(total,"R+");


}

double fLabppo(double *x,double *par)
{

  double y1 = par[3]*(exp(-x[0]/-par[0])-exp(-x[0]/rt))/(-par[0]-rt);
  double y2 = par[4]*(exp(-x[0]/-par[1])-exp(-x[0]/rt))/(-par[1]-rt);
  double y3 = par[5]*(exp(-x[0]/-par[2])-exp(-x[0]/rt))/(-par[2]-rt);
  //double y4 = a4*(exp(-x/-t4)-exp(-x/rt))/(-t4-rt);
  return y1+y2+y3;//+y4
}
