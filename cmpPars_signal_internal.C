{
  TFile *f1 = new TFile("GetSolarMC_containMC_E4to15_Merged_ExtractTres_WaterMP6176_nhit20_GeneralPhysicsMC_WaterSolar_NueRun_r200004to207718.root");
  TFile *f2 = new TFile("fullMerged_InternalOnly_Bi214Tl208_WaterMP6176_nhit20_r200004to207718_s0_p0.root");

  TTree *t1 = (TTree*)f1->Get("T2");
  TTree *t2 = (TTree*)f2->Get("T2");

  TH1F *h1Kldiv = new TH1F("h1Kldiv","klDiv 1",200,0,20);
  TH1F *h2Kldiv = new TH1F("h2Kldiv","klDiv 2",200,0,20);

  TH1F *h1Beta14 = new TH1F("h1Beta14","Beta14 1",400,-2,2);
  TH1F *h2Beta14 = new TH1F("h2Beta14","Beta14 2",400,-2,2);

  TH1F *h1Gtest = new TH1F("h1Gtest","Gtest 1",200,0,2);
  TH1F *h2Gtest = new TH1F("h2Gtest","Gtest 2",200,0,2);

  TH1F *h1Utest = new TH1F("h1Utest","Utest 1",200,0,2);
  TH1F *h2Utest = new TH1F("h2Utest","Utest 2",200,0,2);

  TH1F *h1UdotR = new TH1F("h1UdotR","UdotR 1",400,-2,2);
  TH1F *h2UdotR = new TH1F("h2UdotR","UdotR 2",400,2,2);

  TH1F *h1scaleLogL = new TH1F("h1scaleLogL","scaleLogL 1",200,0,20);
  TH1F *h2scaleLogL = new TH1F("h2scaleLogL","scaleLogL 2",200,0,20);

  TH1F *h1energy = new TH1F("h1energy","energy 1",200,0,20);
  TH1F *h2energy = new TH1F("h2energy","energy 2",200,0,20);

  TH1F *h1Itr = new TH1F("h1Itr","Itr 1",200,0,2);
  TH1F *h2Itr = new TH1F("h2Itr","Itr 2",200,0,2);

  TH1F *h1zfactor = new TH1F("h1zfactor","zfactor 1",200,0,2);
  TH1F *h2zfactor = new TH1F("h2zfactor","zfactor 2",200,0,2);

  t1->Project("h1Kldiv","klDiv");
  t2->Project("h2Kldiv","klDiv");

  t1->Project("h1Beta14","beta14");
  t2->Project("h2Beta14","beta14");

  t1->Project("h1Gtest","Gtest");
  t2->Project("h2Gtest","Gtest");

  t1->Project("h1Utest","Utest");
  t2->Project("h2Utest","Utest");

  t1->Project("h1UdotR","udotR");
  t2->Project("h2UdotR","udotR");

  t1->Project("h1scaleLogL","scaleLogL");
  t2->Project("h2scaleLogL","scaleLogL");

  t1->Project("h1energy","energy");
  t2->Project("h2energy","energy");

  t1->Project("h1zfactor","zfactor");
  t2->Project("h2zfactor","zfactor");

  t1->Project("h1Itr","itr");
  t2->Project("h2Itr","itr");

  h1Kldiv->Scale(1./h1Kldiv->GetMaximum());
  h2Kldiv->Scale(1./h2Kldiv->GetMaximum());

  h1Beta14->Scale(1./h1Beta14->GetMaximum());
  h2Beta14->Scale(1./h2Beta14->GetMaximum());

  h1Gtest->Scale(1./h1Gtest->GetMaximum());
  h2Gtest->Scale(1./h2Gtest->GetMaximum());

  h1Utest->Scale(1./h1Utest->GetMaximum());
  h2Utest->Scale(1./h2Utest->GetMaximum());

  h1UdotR->Scale(1./h1UdotR->GetMaximum());
  h2UdotR->Scale(1./h2UdotR->GetMaximum());

  h1scaleLogL->Scale(1./h1scaleLogL->GetMaximum());
  h2scaleLogL->Scale(1./h2scaleLogL->GetMaximum());

  h1zfactor->Scale(1./h1zfactor->GetMaximum());
  h2zfactor->Scale(1./h2zfactor->GetMaximum());

  h1energy->Scale(1./h1energy->GetMaximum());
  h2energy->Scale(1./h2energy->GetMaximum());

  h1Itr->Scale(1./h1Itr->GetMaximum());
  h2Itr->Scale(1./h2Itr->GetMaximum());

  h2Itr->SetLineColor(kRed);
  h2Kldiv->SetLineColor(kRed);
  h2Beta14->SetLineColor(kRed);
  h2Gtest->SetLineColor(kRed);
  h2Utest->SetLineColor(kRed);
  h2UdotR->SetLineColor(kRed);
  h2scaleLogL->SetLineColor(kRed);
  h2zfactor->SetLineColor(kRed);
  h2energy->SetLineColor(kRed);

  TCanvas c1("c1","klDiv",800,600);     c1->cd();h1Kldiv->Draw(); h2Kldiv->Draw("sames");
  TCanvas c2("c2","beta14",800,600);    c2->cd();h1Beta14->Draw(); h2Beta14->Draw("sames");
  TCanvas c3("c3","Gtest",800,600);     c3->cd();h1Gtest->Draw(); h2Gtest->Draw("sames");
  TCanvas c4("c4","Utest",800,600);     c4->cd();h1Utest->Draw(); h2Utest->Draw("sames");
  TCanvas c5("c5","zfactor",800,600);   c5->cd();h1UdotR->Draw(); h2UdotR->Draw("sames");
  TCanvas c6("c6","energy",800,600);    c6->cd();h1scaleLogL->Draw(); h2scaleLogL->Draw("sames");
  TCanvas c7("c7","itr",800,600);       c7->cd();h1zfactor->Draw(); h2zfactor->Draw("sames");
  TCanvas c8("c8","scaleLogL",800,600); c8->cd();h1energy->Draw(); h2energy->Draw("sames");
  TCanvas c9("c9","udotR",800,600);     c9->cd();h1Itr->Draw(); h2Itr->Draw("sames");

}
