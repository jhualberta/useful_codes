{


  double z[] ={-4999.899,-4500.2,-4001,-3501.399,-2999.7,-2500.399,-1998.499,-1000.099,-499.6,500.8,1000.3,1500.9,2000.001,2500.3,3000.3,3500.9,4000.5,4500.4,4973.567};
  int ndataz = int(sizeof(z)/sizeof(double));
  const int nz = 19;//ndataz;
  double y[]={-0.201,-998.068,-2000.578,-2998.017,-3992.167,-4999.882,5002.057,3973.021,2980.035,1986.669,994.183,496.858,1494.71,2487.539,3496.338,4505.371,-501.268,-1494.927,-2487.912,-3475.769,-4498.62};
  int ndatay = int(sizeof(y)/sizeof(double));
  const int ny = 21;//ndatay;
  double x[]={-4999.043,-4002.525,-3004.229,-2000.155,-992.994,998.133,2011.103,4003.323,5004.868,4503.445,-4489.301,3511.147,2502.805,-3476.709};
  int ndatax = int(sizeof(x)/sizeof(double));
  const int nx = 15;//ndatax; 
  /// for z-scan
  double z_betaM[] = {4.269,4.144,4.051,4.044,4.158,4.231,4.249,4.198,4.6,3.919,4.289,4.134,3.904,4.132,3.95,4.026,3.895,3.756,3.944};
  double z_betaM_rms[]={0.09955,0.1162,0.1063,0.1107,0.1107,0.117,0.1188,0.1197,0.1216,0.1125,0.1277,0.09398,0.09414,0.1092,0.1355,0.1315,0.09658,0.08698,0.0996};
  double z_betaS[] = {16.36,18.19,16.89,17.78,17.81,18.15,18.14,18.14,18.09,18.18,18.19,18.1,17.83,18.04,18.23,17.41,16.59,17.71,19.02};
  double z_betaS_rms[]={0.7386,0.9655,0.5571,0.6232,0.6285,0.6666,0.673,0.6826,0.6958,0.6392,0.7342,0.5364,0.7828,0.7493,0.7497,0.7139,0.7022,0.7402,0.866};
  double z_chi[]={1.172,1.07,1.029,1.076,0.9029,1.111,1.555,0.8367,0.8617,0.9837,1.003,1.374,1.083,1.317,1.479,0.96,1.079,1.026,1.296};

  double zMC_betaM[] = {4.353,4.298,4.419,4.37,4.375,4.394,4.45,4.244,4.33,4.568,4.381,4.396,4.521,4.339,4.412,4.314,4.46,4.431,4.586};
  double zMC_betaM_rms[]={0.08035,0.07536,0.07571,0.07492,0.07497,0.07412,0.07498,0.07467,0.07424,0.07618,0.07427,0.07458,0.07659,0.07534,0.07483,0.07435,0.07672,0.07488,0.07961};
  double zMC_betaS[] = {17.04,17.31,19.4,18.63,19.29,19.75,19.13,18.47,19.25,18.62,18.72,19.05,18.32,18.92,19.17,17.97,17.43,17.23,17.84};
  double zMC_betaS_rms[]={0.3891,0.37,0.4089,0.395,0.4057,0.4136,0.4019,0.3881,0.4013,0.4004,0.3954,0.3988,0.3915,0.3961,0.4037,0.3794,0.3766,0.3677,0.4036};
  double zMC_chi[] = {0.8944,1.218,1.111,1.089,1.065,1.286,1.032,1.14,0.8703,1.725,1,1.217,1.246,1.04,0.9364,0.8872,1.082,1.097,1.104};

  double ez[nz] ={0};

  // y scan
  double y_betaM[] = {3.961,4.163,4.353,4.5,4.029,4.132,4.083,4.086,4.194,4.151,4.138,4.147,4.163,4.133,4.196,4.018,4.092,4.003,4.079,4.127,4.202};
  double y_betaM_rms[]={0.1072,0.1048,0.1203,0.1205,0.12,0.1284,0.1192,0.1177,0.1078,0.1084,0.1184,0.1191,0.1135,0.12,0.1217,0.1024,0.09031,0.1072,0.1042,0.1116,0.1115};
  double y_betaS[] = {18.81,17.9,18.01,17.62,18.01,16.41,16.73,17.51,18.25,18.69,18.95,18.02,18.69,18.64,17.94,17.61,17.91,18.37,18.57,16.95,17.69};
  double y_betaS_rms[]={0.6217,0.58,0.6787,0.6611,0.6715,0.6601,0.6201,0.6418,0.6137,0.6341,0.6957,0.6663,0.6603,0.6946,0.6822,0.5565,0.5033,0.6038,0.5995,0.603,0.6188};
  double y_chi[]={1.188,0.8863,1.227,1.581,1.052,1.36,1.416,1.045,0.9532,0.9704,1.019,0.9075,1.185,1.144,0.8667,0.8173,1.143,0.9856,1.011,1.102,0.8898};

  double yMC_betaM[] = {4.291,4.326,4.445,4.493,4.357,4.393,4.336,4.471,4.356,4.427,4.384,4.422,4.291,4.478,4.288,4.3,4.337,4.357,4.338,4.451,4.249};
  double yMC_betaM_rms[]={0.07522,0.07413,0.07551,0.0754,0.07656,0.08164,0.08208,0.07524,0.07495,0.07514,0.07483,0.07452,0.0737,0.07505,0.07554,0.0751,0.07471,0.07481,0.07604,0.0757,0.0748};
  double yMC_betaS[] = {18.83,18.28,18.9,18.4,18.3,17.17,17.75,18.55,18.4,18.36,18.37,18.94,18.73,18.82,18.63,18.64,18.7,18.38,17.98,18.69,18.12};
  double yMC_betaS_rms[]={0.3965,0.3842,0.4004,0.3879,0.3931,0.3977,0.4103,0.3936,0.3911,0.3891,0.3877,0.3982,0.3928,0.3987,0.3935,0.3923,0.392,0.3872,0.3825,0.3995,0.3808};
  double yMC_chi[] = {1.01,1.117,1.11,1.366,1.468,0.9478,1.43,1.273,1.507,1.155,1.301,1.152,1.274,0.9234,1.062,1.255,1.027,1.112,1.254,0.9995,1.088};

  double ey[ny];// ={0};

  // x scan
  double x_betaM[] = {4.13,4.401,4.223,4.059,3.731,4.406,4.231,3.936,4.16,4.209,4.116,4.057,4.093,3.99,4.06};
  double x_betaM_rms[]={0.06573,0.1357,0.1332,0.1311,0.1236,0.1304,0.1303,0.1278,0.1281,0.09329,0.1287,0.08778,0.08606,0.08979,0.1501};
  double x_betaS[] = {18.05,16.2,18.12,16.45,18.19,20.57,17.99,18.26,17.47,17.38,16.87,17.14,17.62,18.69,16.03};
  double x_betaS_rms[]={0.3682,0.6901,0.7543,0.6717,0.6924,0.8343,0.7429,0.7208,0.7089,0.7688,0.6825,0.469,0.474,0.5187,1.04};
  double x_chi[]={1.314,1.15,1.085,0.8207,1.412,0.8948,1.157,1.118,1.139,1.27,1.322,1.189,1.137,1.161,1.104};

  double xMC_betaM[] = {4.361,4.477,4.336,4.396,4.322,4.317,4.395,4.436,4.423,4.519,4.383,4.454,4.494,4.311,4.519};
  double xMC_betaM_rms[]={0.07481,0.08087,0.07531,0.07526,0.07515,0.07462,0.07515,0.07482,0.0757,0.08114,0.07503,0.07564,0.07575,0.07437,0.07644};
  double xMC_betaS[] = {19.55,17.64,18.3,18.21,19.17,19.13,18.57,19.12,18.36,17.48,18.74,17.89,18.47,19.21,18.08};
  double xMC_betaS_rms[]={0.4081,0.4064,0.387,0.385,0.401,0.3983,0.3906,0.3996,0.3921,0.4082,0.3977,0.3811,0.3936,0.4018,0.3884};
  double xMC_chi[] = {1.123,0.9458,0.8097,1.137,1.135,1.365,1.324,1.14,0.793,1.654,0.9989,1.062,1.065,1.394,1.054};

  double ex[nx];// ={0};

  for(int i = 0;i<nx;i++)
  {ex[i] =0;}
  for(int i = 0;i<ny;i++)
  {ey[i] =0;}
  for(int i = 0;i<nz;i++)
  {ez[i] =0;}

//============ bias data - MC=============================
  double XdeltaBetaM[nx];
  double YdeltaBetaM[ny];
  double ZdeltaBetaM[nz];

  double XdeltaBetaS[nx];
  double YdeltaBetaS[ny];
  double ZdeltaBetaS[nz];

  for(int i = 0;i<nx;i++)
  {
    double deltaBetaM = (x_betaM[i]-xMC_betaM[i])/xMC_betaM[i]*100;
    double deltaBetaS = (x_betaS[i]-xMC_betaS[i])/xMC_betaS[i]*100;
    XdeltaBetaM[i] = deltaBetaM;
    XdeltaBetaS[i] = deltaBetaS;
  }

  for(int i = 0;i<ny;i++)
  {
    double deltaBetaM = (y_betaM[i]-yMC_betaM[i])/yMC_betaM[i]*100;
    double deltaBetaS = (y_betaS[i]-yMC_betaS[i])/yMC_betaS[i]*100;
    YdeltaBetaM[i] = deltaBetaM;
    YdeltaBetaS[i] = deltaBetaS;
  }

  for(int i = 0;i<nz;i++)
  {
    double deltaBetaM = (z_betaM[i]-zMC_betaM[i])/zMC_betaM[i]*100;
    double deltaBetaS = (z_betaS[i]-zMC_betaS[i])/zMC_betaS[i]*100;
    ZdeltaBetaM[i] = deltaBetaM;
    ZdeltaBetaS[i] = deltaBetaS;
  }

  TGraphErrors *gz_betaM = new TGraphErrors(nz,z,z_betaM, ez, z_betaM_rms);
  gz_betaM->SetMarkerStyle(21);
  gz_betaM->SetMarkerSize(1.5);
  gz_betaM->GetXaxis()->SetTitle("mm");
  gz_betaM->SetTitle("#beta_{M}");
  gz_betaM->GetYaxis()->SetRangeUser(3,5);

  TGraphErrors *gz_betaS = new TGraphErrors(nz,z,z_betaS,ez, z_betaS_rms);
  gz_betaS->SetMarkerStyle(21);
  gz_betaS->SetMarkerSize(1.5);
  gz_betaS->GetXaxis()->SetTitle("mm");
  gz_betaS->SetTitle("#beta_{S}");
  gz_betaS->GetYaxis()->SetRangeUser(12,22);

  TGraphErrors *MCgz_betaM = new TGraphErrors(nz,z,zMC_betaM, ez, zMC_betaM_rms);
  MCgz_betaM->SetMarkerStyle(24);
  MCgz_betaM->SetMarkerSize(1.5);
  MCgz_betaM->GetXaxis()->SetTitle("mm");

  TGraphErrors *MCgz_betaS = new TGraphErrors(nz,z,zMC_betaS, ez, zMC_betaS_rms);
  MCgz_betaS->SetMarkerStyle(24);
  MCgz_betaS->SetMarkerSize(1.5);
  MCgz_betaS->GetXaxis()->SetTitle("mm");

  TGraph *gz_chi = new TGraph(nz,z,z_chi);
  gz_chi->SetMarkerSize(1.5);gz_chi->SetMarkerStyle(21);gz_chi->GetYaxis()->SetRangeUser(0,2);

  TGraph *MCgz_chi = new TGraph(nz,z,zMC_chi);
  MCgz_chi->SetMarkerSize(1.5);MCgz_chi->SetMarkerStyle(24);MCgz_chi->GetYaxis()->SetRangeUser(0,2);
  gz_chi->GetXaxis()->SetTitle("mm");
  gz_chi->SetTitle("#chi^{2}/ndf");

/// y scan
//
  TGraphErrors *gy_betaM = new TGraphErrors(ny,y,y_betaM, ey, y_betaM_rms);
  gy_betaM->SetMarkerStyle(21);
  gy_betaM->SetMarkerSize(1.5);
  gy_betaM->GetXaxis()->SetTitle("mm");
  gy_betaM->SetTitle("#beta_{M}");
  gy_betaM->GetYaxis()->SetRangeUser(3,5);

  TGraphErrors *gy_betaS = new TGraphErrors(ny,y,y_betaS, ey, y_betaS_rms);
  gy_betaS->SetMarkerStyle(21);
  gy_betaS->SetMarkerSize(1.5);
  gy_betaS->GetXaxis()->SetTitle("mm");
  gy_betaS->SetTitle("#beta_{S}");
  gy_betaS->GetYaxis()->SetRangeUser(12,22);

  TGraphErrors *MCgy_betaM = new TGraphErrors(ny,y,yMC_betaM, ey, yMC_betaM_rms);
  MCgy_betaM->SetMarkerStyle(24);
  MCgy_betaM->SetMarkerSize(1.5);
  MCgy_betaM->GetXaxis()->SetTitle("mm");

  TGraphErrors *MCgy_betaS = new TGraphErrors(ny,y,yMC_betaS, ey, yMC_betaS_rms);
  MCgy_betaS->SetMarkerStyle(24);
  MCgy_betaS->SetMarkerSize(1.5);
  MCgy_betaS->GetXaxis()->SetTitle("mm");

  TGraph *gy_chi = new TGraph(ny,y,y_chi);
  gy_chi->SetMarkerSize(1.5);gy_chi->SetMarkerStyle(21);gy_chi->GetYaxis()->SetRangeUser(0,2);
  gy_chi->GetXaxis()->SetTitle("mm");
  gy_chi->SetTitle("#chi^{2}/ndf");

  TGraph *MCgy_chi = new TGraph(ny,y,yMC_chi);
  MCgy_chi->SetMarkerSize(1.5);MCgy_chi->SetMarkerStyle(24);MCgy_chi->GetYaxis()->SetRangeUser(0,2);
  //x scan 
//
  TGraphErrors *gx_betaM = new TGraphErrors(nx,x,x_betaM, ex, x_betaM_rms);
  gx_betaM->SetMarkerStyle(21);
  gx_betaM->SetMarkerSize(1.5);
  gx_betaM->GetXaxis()->SetTitle("mm");
  gx_betaM->SetTitle("#beta_{M}");
  gx_betaM->GetYaxis()->SetRangeUser(3,5);

  TGraphErrors *gx_betaS = new TGraphErrors(nx,x,x_betaS, ex, x_betaS_rms);
  gx_betaS->SetMarkerStyle(21);
  gx_betaS->SetMarkerSize(1.5);
  gx_betaS->GetXaxis()->SetTitle("mm");
  gx_betaS->SetTitle("#beta_{S}");
  gx_betaS->GetYaxis()->SetRangeUser(12,22);

  TGraphErrors *MCgx_betaM = new TGraphErrors(nx,x,xMC_betaM, ex, xMC_betaM_rms);
  MCgx_betaM->SetMarkerStyle(24);
  MCgx_betaM->SetMarkerSize(1.5);
  MCgx_betaM->GetXaxis()->SetTitle("mm");

  TGraphErrors *MCgx_betaS = new TGraphErrors(nx,x,xMC_betaS, ex, xMC_betaS_rms);
  MCgx_betaS->SetMarkerStyle(24);
  MCgx_betaS->SetMarkerSize(1.5);
  MCgx_betaS->GetXaxis()->SetTitle("mm");

  TGraph *gx_chi = new TGraph(nx,x,x_chi);
  gx_chi->GetXaxis()->SetTitle("mm");
  gx_chi->SetTitle("#chi^{2}/ndf");
  gx_chi->SetMarkerSize(1.5);gx_chi->SetMarkerStyle(21);gx_chi->GetYaxis()->SetRangeUser(0,2);
  TGraph *MCgx_chi = new TGraph(nx,x,xMC_chi);
  MCgx_chi->SetMarkerSize(1.5);MCgx_chi->SetMarkerStyle(24);MCgx_chi->GetYaxis()->SetRangeUser(0,2);
///////////////////////////////////////////////////////////
  TCanvas *cc = new TCanvas("cc","x",1500,500);
  cc->Divide(3,1);//2,2);
  cc->cd(1);gPad->SetGridx();gPad->SetGridy();
  gx_betaM->Draw("AP"); MCgx_betaM->Draw("P");
  cc->cd(2);gPad->SetGridx();gPad->SetGridy();
  gx_betaS->Draw("AP"); MCgx_betaS->Draw("P");
  cc->cd(3);gPad->SetGridx();gPad->SetGridy();
  gx_chi->Draw("AP"); MCgx_chi->Draw("P");

  TCanvas *cc1 = new TCanvas("cc1","y",1500,500);
  cc1->Divide(3,1);
  cc1->cd(1);gPad->SetGridx();gPad->SetGridy();
  gy_betaM->Draw("AP"); MCgy_betaM->Draw("P");
  cc1->cd(2);gPad->SetGridx();gPad->SetGridy();
  gy_betaS->Draw("AP"); MCgy_betaS->Draw("P");
  cc1->cd(3);gPad->SetGridx();gPad->SetGridy();
  gy_chi->Draw("AP"); MCgy_chi->Draw("P");

  TCanvas *cc2 = new TCanvas("cc2","z",1500,500);
  cc2->Divide(3,1);
  cc2->cd(1);gPad->SetGridx();gPad->SetGridy();
  gz_betaM->Draw("AP"); MCgz_betaM->Draw("P");
  cc2->cd(2);gPad->SetGridx();gPad->SetGridy();
  gz_betaS->Draw("AP"); MCgz_betaS->Draw("P");
  cc2->cd(3);gPad->SetGridx();gPad->SetGridy();
  gz_chi->Draw("AP"); MCgz_chi->Draw("P");

////////data - MC /////////////////
  TCanvas *cc = new TCanvas("cc","",800,800);
  cc->Divide(1,2);
  cc->cd(1);
  TGraph *gbetaMbiasX = new TGraph(nx,x, XdeltaBetaM);
  gbetaMbiasX->GetXaxis()->SetTitle("source position [mm]");
  gbetaMbiasX->GetYaxis()->SetRangeUser(-20,20);
  gbetaMbiasX->SetMarkerStyle(25);
  gbetaMbiasX->SetMarkerColor(kRed);
  gbetaMbiasX->SetMarkerSize(2.5);
//  gbetaMbiasX->SetTitle("X scan");

  TGraph *gbetaMbiasY = new TGraph(ny,y, YdeltaBetaM);
  gbetaMbiasY->GetXaxis()->SetTitle("source position [mm]");
  gbetaMbiasY->SetMarkerStyle(24);
  gbetaMbiasY->SetMarkerColor(kGreen+2);
  gbetaMbiasY->SetMarkerSize(2.5);
//  gbetaMbiasY->SetTitle("Y scan");

  TGraph *gbetaMbiasZ = new TGraph(nz,z, ZdeltaBetaM);
  gbetaMbiasZ->GetXaxis()->SetTitle("source position [mm]");
  gbetaMbiasZ->SetMarkerStyle(22);
  gbetaMbiasZ->SetMarkerColor(kBlue);
  gbetaMbiasZ->SetMarkerSize(2.5);
//  gbetaMbiasZ->SetTitle("Z scan");

  TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9);
  legend1->AddEntry(gbetaMbiasX,"along X-axis","p");
  legend1->AddEntry(gbetaMbiasY,"along Y-axis","p");
  legend1->AddEntry(gbetaMbiasZ,"along Z-axis","p");

  gbetaMbiasX->Draw("AP");
  gbetaMbiasY->Draw("P");
  gbetaMbiasZ->Draw("P");
  legend1->Draw("same");

  cc->cd(2);
  TGraph *gbetaSbiasX = new TGraph(nx,x, XdeltaBetaM);
  gbetaSbiasX->GetXaxis()->SetTitle("source position [mm]");
  gbetaSbiasX->GetYaxis()->SetRangeUser(-20,20);
  gbetaSbiasX->SetMarkerStyle(25);
  gbetaSbiasX->SetMarkerColor(kRed);
  gbetaSbiasX->SetMarkerSize(2.5);
//  gbetaSbiasX->SetTitle("X scan");

  TGraph *gbetaSbiasY = new TGraph(ny,y, YdeltaBetaM);
  gbetaSbiasY->GetXaxis()->SetTitle("source position [mm]");
  gbetaSbiasY->SetMarkerStyle(24);
  gbetaSbiasY->SetMarkerColor(kGreen+2);
  gbetaSbiasY->SetMarkerSize(2.5);
//  gbetaSbiasY->SetTitle("Y scan");

  TGraph *gbetaSbiasZ = new TGraph(nz,z, ZdeltaBetaM);
  gbetaSbiasZ->GetXaxis()->SetTitle("source position [mm]");
  gbetaSbiasZ->SetMarkerStyle(22);
  gbetaSbiasZ->SetMarkerColor(kBlue);
  gbetaSbiasZ->SetMarkerSize(2.5);
//  gbetaSbiasZ->SetTitle("Z scan");

  TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
  legend2->AddEntry(gbetaSbiasX,"along X-axis","p");
  legend2->AddEntry(gbetaSbiasY,"along Y-axis","p");
  legend2->AddEntry(gbetaSbiasZ,"along Z-axis","p");

  gbetaSbiasX->Draw("AP");
  gbetaSbiasY->Draw("P");
  gbetaSbiasZ->Draw("P");
  legend1->Draw("same");




}
