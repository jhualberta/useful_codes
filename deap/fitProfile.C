{
  TFile *ff = new TFile("results_saveOptMC_Opt2023MC_vac_run32149_2p6Gamma_-XYdir_1E5.root");
  hb = (TH2D*)ff->Get("H2D_distVsTdiff_qw");
  // hb->Draw();
  // hb->GetYaxis()->SetTitle("mm");
  TProfile *pp = hb->ProfileX();
  pp->Draw();
  TF1 *g0 = new TF1("g0","pol1",20,120);
  pp->Fit(g0,"R");
  pp->GetXaxis()->SetTitle("ns");
  pp->GetYaxis()->SetTitle("mm");

}
