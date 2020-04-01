//root -l 'scanNtuple.C("Analysis_r0000200003_s000_p000.ntuple.root")'

//fitValid==1
//skyShine>1
//DC mask Bitmask: 
//0x210000000242
//AV coordinates
//3705<z_AV<4500
/*
-5861.0, -2524.0, -998.923
-5861.0, -2524.0, -1.62
-5861.0, -2524.0, -5000.525
-5861.0, -2524.0, -4000.021
-5861.0, -2524.0, -3000.151
-5861.0, -2524.0, -1999.248
-5861.0, -2524.0, 1000.798
-5861.0, -2524.0, 2000.597
-5861.0, -2524.0, 3000.734
-5861.0, -2524.0, 4000.977
-5861.0, -2524.0, 5000.86
-5861.0, -2524.0, 4498.167
-5861.0, -2524.0, 3498.838
-5861.0, -2524.0, 2498.641
-5861.0, -2524.0, 1498.713
-5861.0, -2524.0, -1501.717
-5861.0, -2524.0, -2500.89
-5861.0, -2524.0, -3500.764
-5861.0, -2524.0, -4500.812
*/

void scanWaterN16(const char*fname)
{
  TFile *_file0=TFile::Open(fname);
  TTree *output = (TTree*)_file0->Get("T");

  double sourcePos[3] = {-5861.0, -2524.0, -4000.021};//-5861.0, -2524.0, -5000.525};
  double manipX = sourcePos[0], manipY = sourcePos[1], manipZ = sourcePos[2];
  double sourceR = sqrt(manipX*manipX+manipY*manipY+manipZ*manipZ);

  TH1F *hx = new TH1F("hx","fitted x, fecd==9188",2000,-9000,9000);
  TH1F *hy = new TH1F("hy","fitted y, fecd==9188",2000,-9000,9000);
  TH1F *hz = new TH1F("hz","fitted z, fecd==9188",2000,-9000,9000);
  TH1F *hR = new TH1F("hR","fitted R, fecd==9188",1000,0,9000);

  TH1F *hXbias = new TH1F("hXbias","fitted x bias, fecd==9188",2000,-9000,9000);
  TH1F *hYbias = new TH1F("hYbias","fitted y bias, fecd==9188",2000,-9000,9000);
  TH1F *hZbias = new TH1F("hZbias","fitted z bias, fecd==9188",2000,-9000,9000);

  TH1F *hRbias = new TH1F("hRbias","fitted Rbias, fecd==9188",2000,-9000,9000);

  TH2F *hRZ = new TH2F("hRZ","",1000,0,9000,2000,-9000,9000);
  TH1F *hnhits = new TH1F("hnhits","nhits at interface",1000,0,2000);
  TH1F *hnhits_interface = new TH1F("hnhits_interface","nhits at interface",1000,0,2000);	  
  TH1F *hnhits_av = new TH1F("hnhits_av","nhits at AV",1000,0,2000);
  TH1F *hnhits_boundary = new TH1F("hnhits_boundary","nhits at boundary",1000,0,2000);	  

  TH1F *hmcEdep = new TH1F("hmcEdep","mcEdep",1000,0,20);
  TH1F *hmcEdep_interface = new TH1F("hmcEdep_interface","mcEdep at interface",1000,0,20);
  TH1F *hmcEdep_av = new TH1F("hmcEdep_av","mcEdep at AV",1000,0,20);
  TH1F *hmcEdep_boundary = new TH1F("hmcEdep_boundary","mcEdep at boundary",1000,0,20);

  TH1F *htrigWord = new TH1F("htrigWord","trigWord",100,0,100);
  TH1F *htrigWord_interface = new TH1F("htrigWord_interface","trigWord at interface",100,0,100); 
  TH1F *htrigWord_av = new TH1F("htrigWord_av","trigWord at AV",100,0,10);
  TH1F *htrigWord_boundary = new TH1F("htrigWord_boundary","trigWord at boundary",100,0,100);

  output->Project("hRZ","(posZcor-108):sqrt(posXcor**2+posYcor**2)","");
  gStyle->SetPalette(81);
  TH2F* hZrho = new TH2F("hZrho","",1000,0,9000,2000,-9000,9000);
  output->Project("hZrho","(posZcor-108):sqrt(posXcor**2+posYcor**2)","fecd==9188");

  output->Project("hx","posXcor","fecd==9188");
  output->Project("hy","posYcor","fecd==9188");
  output->Project("hz","posZcor-108","fecd==9188");
  output->Project("hR","sqrt( (posZcor-108)**2+posXcor**2+posYcor**2 )","fecd==9188");

  output->Project("hXbias",Form("posXcor-%f",manipX),"fecd==9188");
  output->Project("hYbias",Form("posYcor-%f",manipY),"fecd==9188");
  output->Project("hZbias",Form("posZcor-108-%f",manipZ), "fecd==9188");

  output->Project("hRbias",Form("( (posXcor-%f)*(%f)+(posYcor-%f)*(%f)+(posZcor-%f)*(%f) )/%f", sourcePos[0], sourcePos[0], sourcePos[1], sourcePos[1], sourcePos[2], sourcePos[2], sourceR), "fecd==9188");
  output->Project("hnhits","nhits","nhits>0");
  output->Project("hnhits_interface","nhits","nhits>0 && posZcor<4450 && posZcor>4425");
  output->Project("hnhits_av","nhits","nhits>0 && sqrt(posXcor**2+posYcor**2+posZcor**2)<6005 && sqrt(posXcor**2+posYcor**2+posZcor**2)>5999");
  output->Project("hnhits_boundary","nhits","nhits>0 && sqrt(posXcor**2+posYcor**2+posZcor**2)<6005 && sqrt(posXcor**2+posYcor**2+posZcor**2)>5999 || (posZcor<4450 && posZcor>4425)");

  //output->Project("hmcEdep","mcEdep","mcEdep>0");
  //output->Project("hmcEdep_interface","mcEdep","posZcor<4450 && posZcor>4425");
  //output->Project("hmcEdep_av","mcEdep","sqrt(posXcor**2+posYcor**2+posZcor**2)<6005 && sqrt(posXcor**2+posYcor**2+posZcor**2)>5999");
  //output->Project("hmcEdep_boundary","mcEdep","sqrt(posXcor**2+posYcor**2+posZcor**2)<6005 && sqrt(posXcor**2+posYcor**2+posZcor**2)>5999 || (posZcor<4450 && posZcor>4425)");

  //output->Project("htrigWord","triggerWord","");
  //output->Project("htrigWord_interface","triggerWord","posZcor<4450 && posZcor>4425");
  //output->Project("htrigWord_av","triggerWord","sqrt(posXcor**2+posYcor**2+posZcor**2)<6005 && sqrt(posXcor**2+posYcor**2+posZcor**2)>5999");
  //output->Project("htrigWord_boundary","triggerWord","sqrt(posXcor**2+posYcor**2+posZcor**2)<6005 && sqrt(posXcor**2+posYcor**2+posZcor**2)>5999 || (posZcor<4450 && posZcor>4425)");
 // hZrho->Draw("colz");
  TString fileName(fname);
  TFile *fresol = new TFile("Resol"+fileName,"recreate");
  fresol->cd();
  hRZ->Write();hZrho->Write();
  hx->Write();hy->Write();hz->Write();hR->Write();hXbias->Write();hYbias->Write(); hZbias->Write(); hRbias->Write();
  hnhits->Write();hnhits_interface->Write();hnhits_av->Write();hnhits_boundary->Write();
  hmcEdep->Write();hmcEdep_interface->Write();hmcEdep_av->Write();hmcEdep_boundary->Write();
  htrigWord->Write();htrigWord_interface->Write();htrigWord_av->Write();htrigWord_boundary->Write();
  fresol->Close();
}
