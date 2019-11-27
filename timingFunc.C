//pars 1 = mean gaus
//pars 2 = sig gaus
//pars 3 = +fractional exponential component
//pars 4 = +slope exp
//pars 5 = is a manually added mean as found by fitting a ordanary Gausian
#include "TF1.h"
#include "TLegend.h"

double rt = 0.8;

double calpar4(double *x,double *par)
{
  double y1 = par[4]*(exp(-x[0]/par[0])-exp(-x[0]/rt))/(par[0]-rt);
  double y2 = par[5]*(exp(-x[0]/par[1])-exp(-x[0]/rt))/(par[1]-rt);
  double y3 = par[6]*(exp(-x[0]/par[2])-exp(-x[0]/rt))/(par[2]-rt);
  double y4 = par[7]*(exp(-x[0]/par[3])-exp(-x[0]/rt))/(par[3]-rt);
  return y1+y2+y3+y4;
}

double calpar3(double *x,double *par)
{
  double y1 = par[3]*(exp(-x[0]/par[0])-exp(-x[0]/rt))/(par[0]-rt);
  double y2 = par[4]*(exp(-x[0]/par[1])-exp(-x[0]/rt))/(par[1]-rt);
  double y3 = par[5]*(exp(-x[0]/par[2])-exp(-x[0]/rt))/(par[2]-rt);
  return y1+y2+y3;
}
double range2 = 500;

void timingFunc()
{
  double calpar3(double *x, double *par);
  double calpar4(double *x, double *par);

  TF1 *fLabppo= new TF1("fLabppo",calpar4,0,range2,8);
  TF1 *fLabppo_0p5= new TF1("fLabppo_0p5",calpar3,0,range2,6);
  TF1 *fLabppo_dda = new TF1("fLabppo_dda",calpar4,0,range2,8);
  TF1 *fLabppo_te = new TF1("fLabppo_te",calpar4,0,range2,8);

  TF1 *fAlphaLabppo = new TF1("fAlphaLabppo",calpar4,0,range2,8);
  TF1 *fAlphaLabppo_0p5 = new TF1("fAlphaLabppo_0p5",calpar3,0,range2,6);
  TF1 *fAlphaLabppo_dda = new TF1("fAlphaLabppo_dda",calpar4,0,range2,8);
  TF1 *fAlphaLabppo_te = new TF1("fAlphaLabppo_te",calpar4,0,range2,8);
  fLabppo_0p5->SetLineColor(kGreen+2);
  fLabppo_dda->SetLineColor(kOrange);
  fLabppo_te->SetLineColor(kBlue);

  fAlphaLabppo_0p5->SetLineColor(kGreen+2);
  fAlphaLabppo_dda->SetLineColor(kOrange);
  fAlphaLabppo_te->SetLineColor(kBlue);

  fAlphaLabppo->SetLineStyle(9);
  fAlphaLabppo_0p5->SetLineStyle(9);
  fAlphaLabppo_dda->SetLineStyle(9);
  fAlphaLabppo_te->SetLineStyle(9);

  fLabppo->SetParameters(4.88,15.4,66.0,400,0.665,0.218,0.083,0.0346);
  fLabppo_dda->SetParameters(5,12.1,33.3,499,0.68,0.21,0.07,0.04);
  fLabppo_0p5->SetParameters(7.19,24.81,269.87,0.553,0.331,0.116);
  fLabppo_te->SetParameters(3.7,10,52,500,0.72,0.23,0.02,0.03);

  fAlphaLabppo->SetParameters(4.79,18.4,92,900,0.427,0.313,0.157,0.1027);
  fAlphaLabppo_0p5->SetParameters(6.56,23.82,224.19,0.574,0.311,0.115);
  fAlphaLabppo_dda->SetParameters(3.8,11.3,65.3,758,0.48,0.32,0.14,0.06);
  fAlphaLabppo_te->SetParameters(3.69,15.5,79.3,489,0.63,0.23,0.07,0.07);


  fLabppo->Draw();fLabppo_dda->Draw("same");fLabppo_0p5->Draw("same");fLabppo_te->Draw("same");
  fAlphaLabppo->Draw("same");fAlphaLabppo_dda->Draw("same");fAlphaLabppo_0p5->Draw("same");fAlphaLabppo_te->Draw("same");
  fLabppo->GetXaxis()->SetTitle("ns");
  TLegend * l = new TLegend(0.78, 0.25, 0.97 ,0.45);
  l->AddEntry(fLabppo, "LAB + 2g/L PPO (e-)");
  l->AddEntry(fLabppo_0p5, "LAB + 0.5g/L PPO (e-)");
  l->AddEntry(fLabppo_dda, "LAB + 2g/L PPO + 0.5% DDA (e-)");
  l->AddEntry(fLabppo_te, "LAB + 2g/L PPO + te + 0.5% DDA (e-)");

  l->AddEntry(fAlphaLabppo, "LAB + 2g/L PPO (#alpha)");
  l->AddEntry(fAlphaLabppo_0p5, "#LAB + 0.5g/L PPO (#alpha)");
  l->AddEntry(fAlphaLabppo_dda, "LAB + 2g/L PPO + 0.5% DDA (#alpha)");
  l->AddEntry(fAlphaLabppo_te, "LAB + 2g/L PPO + te + 0.5% DDA (#alpha)");

  l->Draw();
}
