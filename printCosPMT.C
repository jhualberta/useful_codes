{

  TFile *ff = new TFile("PDF_cosPMT_solarNueMC.root");
  hcosPMT4to15 = (TH1F*)ff->Get("hcosPMT4to15");
  hcosPMT4to5 = (TH1F*)ff->Get("hcosPMT4to5");
  hcosPMT5to6 = (TH1F*)ff->Get("hcosPMT5to6");
  hcosPMT6to7 = (TH1F*)ff->Get("hcosPMT6to7");
  hcosPMT7to8 = (TH1F*)ff->Get("hcosPMT7to8");
  hcosPMT8to9 = (TH1F*)ff->Get("hcosPMT8to9");
  hcosPMT9to10 = (TH1F*)ff->Get("hcosPMT9to10");
  hcosPMT10to15 = (TH1F*)ff->Get("hcosPMT10to15");
  // hcosPMT10to11 = (TH1F*)ff->Get("hcosPMT10to11");
  // hcosPMT11to12 = (TH1F*)ff->Get("hcosPMT11to12");
  // hcosPMT12to13 = (TH1F*)ff->Get("hcosPMT12to13");
  // hcosPMT13to14 = (TH1F*)ff->Get("hcosPMT13to14");
  // hcosPMT14to15 = (TH1F*)ff->Get("hcosPMT14to15");
  int binTot = 40;
  cout<<"4 - 15 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT4to15->GetBinContent(i+1)<<",";
  }
  cout<<endl; 
  cout<<"4 - 5 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT4to5->GetBinContent(i+1)<<",";
  }
  cout<<endl;
  cout<<"5 - 6 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT5to6->GetBinContent(i+1)<<",";
  }
  cout<<endl;
  cout<<"6 - 7 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT6to7->GetBinContent(i+1)<<",";
  }
  cout<<endl;
  cout<<"7 - 8 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT7to8->GetBinContent(i+1)<<",";
  }
  cout<<endl;
  cout<<"8 - 9 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT8to9->GetBinContent(i+1)<<",";
  }
  cout<<endl;
  cout<<"9 - 10 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT9to10->GetBinContent(i+1)<<",";
  }
  cout<<endl;
  cout<<"10 - 15 MeV"<<endl;
  for(int i = 0;i<binTot;i++)
  {
    cout<<hcosPMT10to15->GetBinContent(i+1)<<",";
  }

//  cout<<"10 - 11 MeV"<<endl;
//  for(int i = 0;i<binTot;i++)
//  {
//    cout<<hcosPMT10to11->GetBinContent(i+1)<<",";
//  }
//
//  cout<<"11 - 12 MeV"<<endl;
//  for(int i = 0;i<binTot;i++)
//  {
//    cout<<hcosPMT11to12->GetBinContent(i+1)<<",";
//  }
//
//  cout<<"12 - 13 MeV"<<endl;
//  for(int i = 0;i<binTot;i++)
//  {
//    cout<<hcosPMT12to13->GetBinContent(i+1)<<",";
//  }
//
//  cout<<"13 - 14 MeV"<<endl;
//  for(int i = 0;i<binTot;i++)
//  {
//    cout<<hcosPMT13to14->GetBinContent(i+1)<<",";
//  }
//
//  cout<<"14 - 15 MeV"<<endl;
//  for(int i = 0;i<binTot;i++)
//  {
//    cout<<hcosPMT14to15->GetBinContent(i+1)<<",";
//  }

}
