#include <TF1.h>

double SIN2(double *, double *);
double COS2(double *, double *);

double SIN2PlusCos2(double *, double*);

TF1 *f_sin2;
TF1 *f_cos2;
void add_TF1s(){

  f_sin2 = new TF1("f_sin2", SIN2, -10., 10., 0);
  f_sin2->SetLineColor(2);
  f_sin2->SetNpx(4500);
  f_cos2 = new TF1("f_cos2", COS2, -10., 10., 0);
  f_cos2->SetNpx(4500);
  f_cos2->SetLineColor(4);
  
  TF1 *f_sin2_plus_cos2 = new TF1("f_sin2_plus_cos2", SIN2PlusCos2, -10., 10., 0.);
  f_sin2_plus_cos2->SetLineColor(1);

    f_sin2->Draw();
  f_cos2->Draw("Same");
  f_sin2_plus_cos2->Draw("Same");
  
}


double SIN2( double *x, double *par ){
  return sin(x[0])*sin(x[0]);
}
double COS2( double *x, double *par ){
  return cos(x[0])*cos(x[0]);
}

double SIN2PlusCos2( double *x, double *par){
  return f_sin2->Eval(x[0]) + f_cos2->Eval(x[0]);
}
