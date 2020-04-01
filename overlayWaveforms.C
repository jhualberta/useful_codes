{

   TFile *f1 = new TFile("sumup_ch0_rmSaturate_DistilledMarLabPPO_26Aug2019_attenu_NoCoin.root");
   TFile *f2 = new TFile("sumup_ch1_rmSaturate_DistilledMarLabPPO_26Aug2019_attenu_NoCoin.root");
//   TFile *f3 = new TFile("sumup_ch0_testAttenu_LABPPO_csv_db5_19July.root");

   TH1F *h1 = (TH1F*)f1->Get("hwf0");
   TH1F *h2 = (TH1F*)f2->Get("hwf0");
//   TH1F *h3 = (TH1F*)f3->Get("hwf0");

   h1->SetLineColor(kBlue); 
   h2->SetLineColor(kRed);
//   h3->SetLineColor(kGreen+2);

   h1->Scale(1./h1->Integral());
   h2->Scale(1./h2->Integral());
//   h3->Scale(1./h3->Integral());

   double baseline0 = h1->Integral(0,20)/20;
   double baseline1 = h2->Integral(0,20)/20;
//   double baseline2 = h3->Integral(0,20)/20;

   h2->Scale(baseline0/baseline1);
//   h3->Scale(baseline0/baseline2);

   h2->Draw();
   h1->Draw("sames");
//   h3->Draw("sames");

}
