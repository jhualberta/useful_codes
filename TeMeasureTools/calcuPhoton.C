{

 TFile f("plots_caen_coin.root");
 h2d = (TH2F*)f->Get("e2d_cor");

 //h2d->Draw();
 double single = 13.89;
 double nPhotons = 0;
 for( int i = 100;i<700;i++)
 {
  for(int j = 100;j<700;j++)
  {
    double count = h2d->GetBinContent(i,j); 
    nPhotons + = count/single*9.36742934051144061e-01;
  } 
 }
 cout<<nPhotons<<endl;
}
