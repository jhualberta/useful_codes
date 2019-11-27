double R(double *x, double *p)
{
  double n1 = 1.50;
  double n2 = 1.33;
  double sinTi = x[0];
  double cosTi = sqrt(1-x[0]*x[0]);
  double rs = pow( (n1*cosTi - n2*sqrt(1- pow( n1/n2*sinTi, 2)))/(n1*cosTi + n2*sqrt(1- pow( n1/n2*sinTi, 2))) ,2);
  double rp = pow( (n1*sqrt(1- pow( n1/n2*sinTi, 2)) - n2*cosTi)/(n1*sqrt(1- pow( n1/n2*sinTi, 2))+n2*cosTi), 2);
  return (rs+rp)/2;
}

void fresnelPlot()
{

  TF1 *f1 = new TF1("f1",R,0,1,0);    

  f1->Draw();

}
