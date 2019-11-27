{
  int pmtID = 44;
  TString pmtid; pmtid.Form("%d",pmtID);
  TString filename = "multipath_scintwater_gtID_0_2scintpdf.root";
  TFile *ff = new TFile(filename,"read");
  // get LLH vs PMTid and LLH derivative vs PMTid
  TH2D *Hx_2_0 = (TH2D*)ff->Get("Hx_2_0");
  TH2D *DHx_2_0 = (TH2D*)ff->Get("DHx_2_0");

  TH2D *Hy_2_0 = (TH2D*)ff->Get("Hy_2_0");
  TH2D *DHy_2_0 = (TH2D*)ff->Get("DHy_2_0");

  TH2D *Hz_2_0 = (TH2D*)ff->Get("Hz_2_0");
  TH2D *DHx_2_0 = (TH2D*)ff->Get("DHx_2_0");

  // project x, y, z and project for certain PMT
  TH1D *h1x = Hx_2_0->ProjectionX("Hx project to PMT",pmtID, pmtID+1);
  TH1D *hxall = Hx_2_0->ProjectionX();
  TH1D *h1dx = DHx_2_0->ProjectionX("DHx project to PMT",pmtID, pmtID+1);
  TH1D *h1dxall = DHx_2_0->ProjectionX();
 
  TH1D *h1y = Hy_2_0->ProjectionX("Hy project to PMT",pmtID, pmtID+1);
  TH1D *hyall = Hy_2_0->ProjectionX();
  TH1D *h1dy = DHy_2_0->ProjectionX("DHy project to PMT",pmtID, pmtID+1);
  TH1D *h1dyall = DHy_2_0->ProjectionX();

  TH1D *h1z = Hz_2_0->ProjectionX("Hz project to PMT",pmtID, pmtID+1);
  TH1D *hzall = Hz_2_0->ProjectionX();
  TH1D *h1dz = DHz_2_0->ProjectionX("DHz project to PMT",pmtID, pmtID+1);
  TH1D *h1dzall = DHz_2_0->ProjectionX();
 
  TString newfile(filename(0,filename.Index(".")));
  newfile = newfile+"_pmtID"+pmtid+".root";
  TFile *f = new TFile(newfile,"recreate");
  f->cd();
  h1x->Write();h1y->Write();h1z->Write();
  hxall->Write();hyall->Write();hzall->Write();
  h1dx->Write();h1dy->Write();h1dz->Write();
  h1dxall->Write();h1dyall->Write();h1dzall->Write();
  f->Close();
}
