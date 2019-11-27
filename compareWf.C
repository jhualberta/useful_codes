{
  TFile *f1 = new TFile("dumpWaveform_ch0_Humidity_TeDDA_noCoin3_28Feb_0.root"); 
  TFile *f2 = new TFile("dumpWaveform_ch1_Humidity_TeDDA_noCoin3_28Feb_0.root"); 
  int num1, num2;
  cout<<"put event in ch0"<<endl;
  cin>>num1; 
  cout<<"put event in ch1"<<endl;
  cin>>num2; 
  TH1F *h1 = (TH1F*)f1->Get(Form("hwf%u",num1));
  TH1F *h2 = (TH1F*)f2->Get(Form("hwf%u",num1));
  double s1 = h1->GetMaximum();
  double s2 = h2->GetMaximum();
  h2->SetLineColor(kRed);
  
  double coarseGain = 80;
  Int_t npeaks = 4;
  Int_t recordLength = 252;
  Int_t adcChannel = 8192;
  double attenu = 0.25;//25%
  int delayW = 1; // 1 ns for DT5751
  TH1F *hp1 = new TH1F("hp1","", recordLength, 0, recordLength);
  TH1F *hp2 = new TH1F("hp2","", recordLength, 0, recordLength);

  int shortgate[2] = {32,47};
  int longgate[2] = {32,70};
  
  int preTrigSum1 = 0;
  for(int i = 0;i<recordLength;i++) {
    if(i<32) preTrigSum1+=h1->GetBinContent(i+1);
    double val = h1->GetBinContent(i+1);
    val = val*attenu;
    double val_c = h1->GetBinContent(i+1+delayW)*(-1);
    val_c = val_c+val;
    if(i+delayW<recordLength) {
      hp1->SetBinContent(i+1,val_c);
    }
    else hp1->SetBinContent(i+1,hp1->GetBinContent(i-1));
  }

  int preTrigSum2 = 0;
  for(int i = 0;i<recordLength;i++) {
    if(i<32) preTrigSum2+=h2->GetBinContent(i+1);
    double val = h2->GetBinContent(i+1);
    val = val*attenu;
    double val_c = h2->GetBinContent(i+1+delayW)*(-1);
    val_c = val_c+val;
    if(i+delayW<recordLength) {
      hp2->SetBinContent(i+1,val_c);
    }
    else hp2->SetBinContent(i+1,hp2->GetBinContent(i-1));
  }

  double deltaGate = shortgate[1]-shortgate[0]; 
  double baseline1 = preTrigSum1/32;
  double baseline2 = preTrigSum2/32;
  double e1 = h1->Integral(shortgate[0],shortgate[1])-baseline1*deltaGate;
  double e2 = h2->Integral(shortgate[0],shortgate[1])-baseline2*deltaGate;
 
  double baselineL = 315*(longgate[1]-longgate[0]);
  double e1L = h1->Integral(longgate[0],longgate[1])-baselineL;
  double e2L = h2->Integral(longgate[0],longgate[1])-baselineL;
 
  cout<<"baseline1,2 "<<baseline1<<", "<<baseline2<<" "<<e1<<" "<<e2<<endl; 
  cout<<e1L<<" "<<e2L<<endl; 

  double s3 = hp1->GetMaximum();
  double s4 = hp2->GetMaximum();
  hp2->SetLineColor(kRed);
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
  c1->Divide(1,2);
  c1->cd(1);
  if(s1>s2) {
  h1->Draw();
  h2->Draw("same");  
  }
  else {
  h2->Draw();
  h1->Draw("same");
  }

  c1->cd(2);
  if(s3>s4) {
  hp1->Draw();
  hp2->Draw("same");
  }
  else {
  hp2->Draw();
  hp1->Draw("same");
  }


}
