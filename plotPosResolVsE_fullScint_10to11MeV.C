{
  double E[] = {10,11,12,13,14,15,16,17,18,19,20};

  //// fill inner_av iso
  double xresol[] = {7.94138e+01, 8.48104e+01, 9.11535e+01, 9.73609e+01, 9.67649e+01, 1.04722e+02, 1.08080e+02, 1.07176e+02, 1.13414e+02, 1.20769e+02, 1.26842e+02};
  double yresol[] = {8.55495e+01, 8.62108e+01, 8.76837e+01, 9.29311e+01, 9.95121e+01, 9.91880e+01, 1.06960e+02, 1.12131e+02, 1.15150e+02, 1.20439e+02, 1.25178e+02};
  double zresol[] = {8.04833e+01, 9.03408e+01, 9.08782e+01, 9.06763e+01, 1.00712e+02, 9.87777e+01, 9.98450e+01, 1.13144e+02, 1.13341e+02, 1.22786e+02, 1.25530e+02};
//  double Rresol[] = {};

  double xresolErr[] = {1.28791e+00,1.29403e+00,1.53798e+00,1.51344e+00,1.62599e+00,1.90998e+00,1.73658e+00,1.74991e+00,1.79851e+00,2.18370e+00,2.36537e+00};
  double yresolErr[] = {1.50049e+00,1.32714e+00,1.42953e+00,1.46949e+00,1.62754e+00,1.65675e+00,1.79896e+00,1.85377e+00,1.91759e+00,2.18174e+00,2.21771e+00};
  double zresolErr[] = {1.30133e+00,1.55287e+00,1.51599e+00,1.54840e+00,1.74720e+00,1.51576e+00,1.53641e+00,2.17387e+00,1.97189e+00,2.15407e+00,2.23129e+00};
//  double RresolErr[] = {};

  double xbias[] = {1.00415e+00,-1.94555e+00,-3.27215e+00,3.59758e+00,-4.72519e+00,4.01377e+00,-8.94133e-01,1.40061e+00,-1.84042e+00,8.03495e-01,7.92856e-01}; // mu
  double ybias[] = {3.32821e+00,-2.44416e+00,-1.12114e+00,-1.67161e+00,4.83944e+00,-5.54428e+00,-3.21380e+00,5.57990e-01,4.98040e+00,6.57235e+00,-2.54451e+00};
  double zbias[] = {3.79079e+00, 1.92557e+00, 1.45339e+00, 3.55728e+00,3.16212e+00, 9.92673e-01,-5.84700e-01,2.70534e+00,8.04093e-01,7.12232e+00, 7.17366e+00};
//  double rbias[] = {};

  double xbiasErr[] = {1.87746e+00,2.00054e+00,2.17680e+00,2.31506e+00,2.30415e+00,2.45448e+00,2.52505e+00,2.50458e+00,2.60471e+00,2.81840e+00,2.97472e+00};
  double ybiasErr[] = {2.04598e+00,2.04830e+00,2.07540e+00,2.23231e+00,2.38458e+00,2.42228e+00,2.52267e+00,2.61866e+00,2.66695e+00,2.79775e+00,2.88215e+00};
  double zbiasErr[] = {1.89282e+00,2.17988e+00,2.15613e+00,2.14249e+00,2.40025e+00,2.34663e+00,2.33077e+00,2.66591e+00,2.60100e+00,2.85683e+00,2.91161e+00};
//  double rbiasErr[] = {};

  int N = sizeof(E)/sizeof(double);
  TGraphErrors *gr1 = new TGraphErrors(N,E,xresol);
  TGraphErrors *gr2 = new TGraphErrors(N,E,yresol);
  TGraphErrors *gr3 = new TGraphErrors(N,E,zresol);

  TGraphErrors *gr4 = new TGraphErrors(N,E,rbias);
  gr1->SetMarkerStyle(24);
  gr2->SetMarkerStyle(25);
  gr3->SetMarkerStyle(26);
  gr1->SetMarkerColor(kRed);
  gr2->SetMarkerColor(kGreen+2);
  gr1->SetMarkerSize(2);
  gr2->SetMarkerSize(2);
  gr3->SetMarkerSize(2);

  TCanvas *c = new TCanvas("c","",1000,800);
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(gr1,"#sigma_{x}","p");
  legend->AddEntry(gr2,"#sigma_{y}","p");
  legend->AddEntry(gr3,"#sigma_{z}","p");
  // gpad->SetGridx();
  gr1->GetXaxis()->SetTitle("E [MeV]");
  gr1->GetYaxis()->SetTitle("resolution [mm]");
  gr1->Draw("AP");
  gr2->Draw("P");
  gr3->Draw("P");
  legend->Draw();

  TCanvas *c1 = new TCanvas("c1","",1000,800);
  gr4->SetMarkerStyle(25);
  gr4->Draw("AP");

}
