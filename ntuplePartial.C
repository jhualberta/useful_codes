//root -l 'scanNtuple.C("Analysis_r0000200003_s000_p000.ntuple.root")'

//fitValid==1
//skyShine>1
//DC mask Bitmask: 
//0x210000000242
//AV coordinates
//3705<z_AV<4500


void ntuplePartial(const char*fname)
{
  TFile *_file0=TFile::Open(fname);
  TTree *output = (TTree*)_file0->Get("output");
  TH2F *hRhoZ = new TH2F("hRhoZ","",1000,0,9000,1000,3000,9000);

  output->Project("hRhoZ","(posz-108):sqrt(posx**2+posy**2)","fitValid && partialFit && (dcFlagged & 0x210000000242) == 0x210000000242 && sqrt(posx**2+posy**2+(posz-108)**2)<5500 && skyShine>1");
  _file0->Write();
  _file0->Close();
}
