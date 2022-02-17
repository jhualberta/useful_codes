{
  TFile *f1 = new TFile("GetSolarMC_E4to15_symKlDiv_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to203602_s0_p0.root");
  TFile *f2 = new TFile("GetSolarMC_E4to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214_AvRun_r200004to203602_s0_p0.root");
  TFile *f3 = new TFile("GetSolarMC_E4to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214_ExwaterRun_r200004to203602_s0_p0.root");
  TFile *f4 = new TFile("GetSolarMC_E4to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214Run_r200004to203602_s0_p0.root");
  TFile *f5 = new TFile("GetSolarMC_E4to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208_AvRun_r200004to203602_s0_p0.root");
  TFile *f6 = new TFile("GetSolarMC_E4to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208_ExwaterRun_r200004to203602_s0_p0.root");
  TFile *f7 = new TFile("GetSolarMC_E4to15_Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208Run_r200004to203602_s0_p0.root");

  TFile *f8 = new TFile("GetSolarMC_E4to15_MergedHalf_symKldiv_allbkgsTl208andBi214_r200004to203602_s0_p0.root");

  // ("GetSolarMC_E4to15_MergedHalf_symKldiv_allbkgsTl208andBi214_r200004to203602_s0_p0.root");
  TTree *t1 = (TTree*)f1->Get("T2");
  TTree *t2 = (TTree*)f2->Get("T2");
  TTree *t3 = (TTree*)f3->Get("T2");
  TTree *t4 = (TTree*)f4->Get("T2");
  TTree *t5 = (TTree*)f5->Get("T2");
  TTree *t6 = (TTree*)f6->Get("T2");
  TTree *t7 = (TTree*)f7->Get("T2");
  TTree *t8 = (TTree*)f8->Get("T2");

  TH1F *h1 = new TH1F("h1","klDiv 1",200,0,20);
  TH1F *h2 = new TH1F("h2","klDiv 2",200,0,20);
  TH1F *h3 = new TH1F("h3","klDiv 3",200,0,20);
  TH1F *h4 = new TH1F("h4","klDiv 4",200,0,20);
  TH1F *h5 = new TH1F("h5","klDiv 5",200,0,20);
  TH1F *h6 = new TH1F("h6","klDiv 6",200,0,20);
  TH1F *h7 = new TH1F("h7","klDiv 7",200,0,20);
  TH1F *h8 = new TH1F("h8","klDiv 8",200,0,20);

  t1->Project("h1","klDiv");
  t2->Project("h2","klDiv");
  t3->Project("h3","klDiv");
  t4->Project("h4","klDiv");
  t5->Project("h5","klDiv");
  t6->Project("h6","klDiv");
  t7->Project("h7","klDiv");
  t8->Project("h8","klDiv");

  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
  h4->Scale(1./h4->Integral());
  h5->Scale(1./h5->Integral());
  h6->Scale(1./h6->Integral());
  h7->Scale(1./h7->Integral());
  h8->Scale(1./h8->Integral());

  h2->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h4->SetLineColor(kRed);
  h5->SetLineColor(kRed);
  h6->SetLineColor(kRed);
  h7->SetLineColor(kRed);
  h8->SetLineColor(kRed);

  TCanvas c1("c1","Bi214 AV",800,600);c1->cd();h1->Draw(); h2->Draw("same");
  TCanvas c2("c2","Bi214 ex",800,600);  c2->cd();h1->Draw(); h3->Draw("same");
  TCanvas c3("c3","Bi214 run",800,600); c3->cd();h1->Draw(); h4->Draw("same");
  TCanvas c4("c4","Tl208 AV",800,600);  c4->cd();h1->Draw(); h5->Draw("same");
  TCanvas c5("c5","Tl208 ex",800,600);  c5->cd();h1->Draw(); h6->Draw("same");
  TCanvas c6("c6","Tl208 run",800,600); c6->cd();h1->Draw(); h7->Draw("same");
  TCanvas c7("c7","all",800,600);       c7->cd();h1->Draw(); h8->Draw("same");


}
