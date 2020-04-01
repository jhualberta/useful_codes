//root -l 'scanNtuple.C("Analysis_r0000200003_s000_p000.ntuple.root")'

//fitValid==1
//skyShine>1
//DC mask Bitmask: 
//0x210000000242
//AV coordinates
//3705<z_AV<4500


void scanWaterNtuple(const char*fname)
{
  TFile *_file0=TFile::Open(fname);
  TTree *output = (TTree*)_file0->Get("T");
  TH2F *hRZ = new TH2F("hRZ","",1000,0,9000,1000,3000,9000);
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

  hZrho->Draw("colz");
  TString fileName(fname);
  TFile *fresol = new TFile("Resol"+fileName,"recreate");
  fresol->cd();
  hRZ->Write();hZrho->Write();
  hnhits->Write();hnhits_interface->Write();hnhits_av->Write();hnhits_boundary->Write();

  hmcEdep->Write();hmcEdep_interface->Write();hmcEdep_av->Write();hmcEdep_boundary->Write();

  htrigWord->Write();htrigWord_interface->Write();htrigWord_av->Write();htrigWord_boundary->Write();

  fresol->Close();
}
