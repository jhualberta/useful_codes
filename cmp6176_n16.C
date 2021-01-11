{
  int number;
  cout<<"run number:  "<<endl;
  cin>>number;
  TString runN;
  runN.Form("%d",number);
  TString ratFile = "Merged_Rat6176_Calibration_r";
  TString ratMCFile = "Merged_Filtered_AntiNeutrinoMC_WaterN16source_r";

  TString mpwFile = "Merged_MPW_nhit6_Calibration_r";
  TString mpwMCFile = "Merged_MPWsolarTest_WaterMP6176_AntiNeutrinoMC_WaterN16sourceRun_r";

  TFile *fRat = new TFile(ratFile+runN+"_s0.root");
  TFile *fRatmc = new TFile(ratMCFile+runN+"_s0.root");
  TFile *fMP = new TFile(mpwFile+runN+"_s0.root");
  TFile *fMPmc = new TFile(mpwMCFile+runN+"_s0.root");

  TTree *tMP = (TTree*)fMP->Get("T");
  TTree *tMPmc = (TTree*)fMPmc->Get("T");
  TTree *tRat = (TTree*)fRat->Get("T");
  TTree *tRatmc = (TTree*)fRatmc->Get("T");
  // energy
  //
  //
  double hi = 12, lo = 0;
  int bin = 80; // 0.15
  TH1F *hMP = new TH1F("hMP","",bin, lo, hi);
  TH1F *hMPmc = new TH1F("hMPmc","",bin,lo,hi);
  TH1F *hRat = new TH1F("hRat","",bin,lo,hi);
  TH1F *hRatmc = new TH1F("hRatmc","",bin,lo,hi);

//  TH1F *hR1 = new TH1F("hR1","",1000,0,9000);
//  TH1F *hR2 = new TH1F("hR2","",1000,0,9000);
//
//  TH2F *hRhoZ1 = new TH2F("hRhoZ1","",1000,0,9000, 2000,-9000,9000);
//  TH2F *hRhoZ2 = new TH2F("hRhoZ2","",1000,0,9000, 2000,-9000,9000);

 // posRad<5300 && beta14<0.95 && beta14>-0.12 && itr>0.55 && nhits>=20
//  t1->Draw("energyOrigin>>h1","posRad<5300 && nhits>=20 && energyOrigin>9.9");// && energyOrigin<10.1");
//  t2->Draw("energy>>h2","posRad<5300 && nhits>=20 && energyOrigin>9.9","sames");// && energyOrigin<10.1","sames");

 tMP->Project("hMP","energy","nhits>5");//"posRad<5300");// && energyOrigin<10.1");
 tMPmc->Project("hMPmc","energy","nhits>5");//"posRad<5300");// && energyOrigin<10.1","sames");
 tRat->Project("hRat","energy","nhits>5");//posRad<5300");// && energyOrigin<10.1");
 tRatmc->Project("hRatmc","energy","nhits>5");//posRad<5300");// && energyOrigin<10.1","sames");


// nhits
// TH1F *hMP = new TH1F("hMP","",100,0,100);
// TH1F *hMPmc = new TH1F("hMPmc","",100,0,100);
// TH1F *hRat = new TH1F("hRat","",100,0,100);
// TH1F *hRatmc = new TH1F("hRatmc","",100,0,100);
//
// tMP->Project("hMP","nhits","");// && nhitsOrigin<10.1");
// tMPmc->Project("hMPmc","nhits","");// && nhitsOrigin<10.1","sames");
// tRat->Project("hRat","nhits","");// && nhitsOrigin<10.1");
// tRatmc->Project("hRatmc","nhits","");// && nhitsOrigin<10.1","sames");

//// beta14
// TH1F *hMP = new TH1F("hMPW_data","",200,-2,2);
// TH1F *hMPmc = new TH1F("hMPW_mc","",200,-2,2);
// TH1F *hRat = new TH1F("hRat","",200,-2,2);
// TH1F *hRatmc = new TH1F("hRatmc","",200,-2,2);
//
// tMP->Project("hMPW_data","beta14","");
// tMPmc->Project("hMPW_mc","beta14","");
// tRat->Project("hRat","beta14","");
// tRatmc->Project("hRatmc","beta14","");
// hMP->GetXaxis()->SetTitle("#beta_{14}");
// hMP->GetYaxis()->SetTitle("Scaled");
// hMP->SetLineColor(kRed);
// hMPmc->SetLineColor(kRed);
//
// hMPmc->SetLineStyle(2);
// hRatmc->SetLineStyle(2);
// hMP->Sumw2();
// hMPmc->Sumw2();
//// hMP->Draw();
//// hMPmc->Draw("sames");
// hRat->Draw();
// hRatmc->Draw("sames");
// hMP->Scale(1./hMP->Integral());
// hMPmc->Scale(1./hMPmc->Integral());
// hRat->Scale(1./hRat->Integral());
// hRatmc->Scale(1./hRatmc->Integral());



TCanvas c("c","",800,600);
c.cd();
hMP->GetYaxis()->SetTitle("Normalized counts/0.15 MeV");
hMP->GetXaxis()->SetTitle("E [MeV]");

hMP->SetLineWidth(2);
hMP->SetMarkerStyle(21);
hMPmc->SetLineColor(kBlue);
hMP->Scale(hMPmc->Integral()/hMP->Integral());
hMP->Draw("E2");
hMPmc->Draw("same");


// TH2F *hMP = new TH2F("hMPW_data","",125,0,15,125,-1,2);
// TH2F *hMPmc = new TH2F("hMPW_mc","",125,0,15,125,-1,2);
// TH2F *hRat = new TH2F("hRat","",125,0,15,125,-1,2);
// TH2F *hRatmc = new TH2F("hRatmc","",125,0,15,125,-1,2);
//
// tMP->Project("hMPW_data","beta14:energy","");
// tMPmc->Project("hMPW_mc","beta14:energy","");
// tRat->Project("hRat","beta14:energy","");
// tRatmc->Project("hRatmc","beta14:energy","");
//
// TCanvas c1("c1","MPW data",800,600);
// hMP->GetXaxis()->SetTitle("E [MeV]");
// hMP->GetYaxis()->SetTitle("#beta_{14}");
// c1.SetLogz();c1.SetGridx();c1.SetGridy();
// hMP->Draw("colz");

// TCanvas c2("c2","MPW mc",800,600);
// hMPmc->GetXaxis()->SetTitle("E [MeV]");
// hMPmc->GetYaxis()->SetTitle("#beta_{14}");
// c2.SetLogz();c2.SetGridx();c2.SetGridy();
// hMPmc->Draw("colz");
//
// TCanvas c3("c3","Rat data",800,600);
// hRat->GetXaxis()->SetTitle("E [MeV]");
// hRat->GetYaxis()->SetTitle("#beta_{14}");
// c3.SetLogz();c3.SetGridx();c3.SetGridy();
// hRat->Draw("colz");
//
// TCanvas c4("c4","Rat mc",800,600);
// hRatmc->GetXaxis()->SetTitle("E [MeV]");
// hRatmc->GetYaxis()->SetTitle("#beta_{14}");
// c4.SetLogz();c4.SetGridx();c4.SetGridy();
// hRatmc->Draw("colz");





//  TCanvas c1("c1","MPW ",400,300);
//  t1->Draw("posz:sqrt(posx*posx+posy*posy)>>hRhoZ1","energyOrigin>=10 && nhits>20","colz");
//
//
//  TCanvas c2("c2","Rat ",400,300);
//  t2->Draw("posz:sqrt(posx*posx+posy*posy)>>hRhoZ2","energy>=10 && nhits>20","colz");

//  TH1F *hdirx1 = new TH1F("hdirx1","",1000,-100,100);
//  TH1F *hdirx2 = new TH1F("hdirx2","",1000,-100,100);
//  hdirx1->SetLineColor(kRed);
//  t1->Draw("diry/dirx>>hdirx1","nhits>=20 && energyOrigin<10.1 && energyOrigin>9.9");
//  t2->Draw("diry/dirx>>hdirx2","nhits>=20 && energyOrigin<10.1 && energyOrigin>9.9","sames");

//  h1->Scale(1./h1->Integral());
//  h2->Scale(1./h2->Integral());

//  t1->Draw("posRad>>hR1","posRad<5300 && nhits>=20");
//  t2->Draw("posRad>>hR2","posRad<5300 && nhits>=20","sames");


}
