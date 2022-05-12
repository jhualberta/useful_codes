{
//  TFile *f1 = new TFile("saveMCsum_nue_numu_pdfs.root");
//  TFile *f2 = new TFile("saveMCsum_smearEUp_nue_numu_pdfs.root");//DirDown_nue_numu_pdfs.root");
//  TFile *f3 = new TFile("saveMCsum_smearEresolDown_nue_numu_pdfs.root");//DirUp_nue_numu_pdfs.root");

  TFile *f1 = new TFile("saveFlux8parsFinal_MLP_AllBkgsE5to15_ntupleTMVA_nue_numu_pdfs.root");
//  TFile *f2 = new TFile("saveMCsum_smearEUp_nue_numu_pdfs.root");//DirDown_nue_numu_pdfs.root");
  TFile *f3 = new TFile("saveSmear8pars_MLP_SmearedMC_EresolDown_5to15MeV.root");//DirUp_nue_numu_pdfs.root");

  TH1F *h1 = (TH1F*)f1->Get("hcosSunSumNueNumu_peeEtrue");
//  TH1F *h2 = (TH1F*)f2->Get("hcosSunSumNueNumu_peeEtrue");
  TH1F *h3 = (TH1F*)f3->Get("hcosSunSumNueNumu_peeEtrue");
  TCanvas c("c","",800,600);
  h3->Draw();h3->SetLineStyle(2); h3->SetLineColor(kRed);h3->GetXaxis()->SetTitle("cos#theta_{sun}");h3->GetYaxis()->SetTitle("Counts/190.33 days/0.05 bin");
//  h2->Draw();h2->SetLineStyle(2); h2->SetLineColor(kRed);h2->GetXaxis()->SetTitle("cos#theta_{sun}");h2->GetYaxis()->SetTitle("counts");
  h1->Draw("same");

  //h2->Draw("same");h2->SetLineStyle(3);h2->SetLineColor(kBlue);








}
