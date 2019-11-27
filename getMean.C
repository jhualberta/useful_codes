
//root -l 'getMean.C("sim_2p5MeVbeta_z0_splitZ0_1.root")'
void getMean(char *fname)
{
 TFile *_file0 = TFile::Open(fname);
 TTree *T = (TTree*)_file0->Get("T");
 TH1F *hNhits = new TH1F("hNhits","",2000,0,2000); 
 T->Draw("ds.evs.nhits_cleaned>>hNhits","ds.evs.trigType");
 hNhits->Fit("gaus","q");
 double mean = hNhits->GetFunction("gaus")->GetParameter(1);
 double sigma = hNhits->GetFunction("gaus")->GetParameter(2);
 cout<<mean<<" "<<sigma<<endl;
}
