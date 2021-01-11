{
  TString fname;	
  cout<<"put filename"<<endl;
  cin>>fname;
  TFile *file = new TFile(fname);

  TH2F *hRhoZ = new TH2F("hRhoZ","",1000,0,9000,2000,-9000,9000);
  T->Draw("posz:sqrt(posx**2+posy**2)>>hRhoZ","","colz");
  hRhoZ->RebinX(4);
  hRhoZ->RebinY(8);
  hRhoZ->GetXaxis()->SetTitle("#rho [mm]");
  hRhoZ->GetYaxis()->SetTitle("z [mm]");


}
