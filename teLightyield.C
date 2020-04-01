{

  double concen[] = {0.5, 1, 2, 3, 4};
  double biller[] = {0.7,0.6,0.45,0.4,0.32};
  double ua1900[] = {0.927,0.593,0.221,0.0367,0.0031};
  double lu1900[] = {0.224,0.077,0.016,0.006,0.003};
  double ua2000[] = {0.929, 0.621, 0.294, 0.0288, 0.019}; 
  double lu2000[] = {0.295,0.146, 0.053, 0.029, 0.016};
 
  TGraph *g = new TGraph(5, concen, biller);
  g->SetMarkerStyle(21);
  g->SetMarkerColor(kOrange+1);
  g->Draw("AP");

  TGraph *gua1900 = new TGraph(5, concen, ua1900);
  gua1900->SetMarkerStyle(21);
  gua1900->SetMarkerColor(kBlue);
  gua1900->Draw("P");

  TGraph *gua2000 = new TGraph(5, concen, ua2000);
  gua2000->SetMarkerStyle(21);
  gua2000->SetMarkerColor(kGreen);
  gua2000->Draw("P");

  TGraph *glu1900 = new TGraph(5, concen, lu1900);
  glu1900->SetMarkerStyle(21);
  glu1900->SetMarkerColor(kRed);
  glu1900->Draw("P");

  TGraph *glu2000 = new TGraph(5, concen, lu2000);
  glu2000->SetMarkerStyle(21);
  glu2000->SetMarkerColor(kBlack);
  glu2000->Draw("P");




}
