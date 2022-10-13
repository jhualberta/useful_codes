{
    TFile ff1("SaveTrack_trackMCusher_run30695_Po210_vacuum_avBulk_thin0p00000001_1000evts_rmStepLimiter.root");
    TFile ff2("SaveTrack_trackMCusher_run30695_Po210_vacuum_avBulk_thin0p00001_1000evts_rmStepLimiter.root");
    TFile ff3("SaveTrack_trackMCusher_run30695_Po210_vacuum_avBulk_thin0p001_1000evts_rmStepLimiter.root");
    TFile ff4("SaveTrack_trackMCusher_run30695_Po210_vacuum_avBulk_thin0p05_1e4evts_rmStepLimiter.root");
    TFile ff5("SaveTrack_trackMCusher_run30695_Po210_vacuum_avBulk_thin0p1_1e4evts_rmStepLimiter.root");

    ff1->cd();
    h_tpb2av = (TH1F*)ff1->Get("H_E_tpb2av");

    h_tpb2vac_thin0p00000001 = (TH1F*)ff1->Get("H_E_tpb2vac");
    h_tpb2vac_thin0p00001    = (TH1F*)ff2->Get("H_E_tpb2vac");
    h_tpb2vac_thin0p001      = (TH1F*)ff3->Get("H_E_tpb2vac");
    h_tpb2vac_thin0p05       = (TH1F*)ff4->Get("H_E_tpb2vac");
    h_tpb2vac_thin0p1        = (TH1F*)ff5->Get("H_E_tpb2vac");

    h_tpb2vac_thin0p00000001->SetLineColor(kRed);
    h_tpb2vac_thin0p00001->SetLineColor(kBlue);
    h_tpb2vac_thin0p001->SetLineColor(kBlack);
    h_tpb2vac_thin0p05->SetLineColor(kGreen+2);
    h_tpb2vac_thin0p1->SetLineColor(kOrange+1);

    double total = h_tpb2vac_thin0p00000001->Integral();
    double int1 = h_tpb2vac_thin0p00000001->Integral(h_tpb2vac_thin0p00000001->GetXaxis()->FindBin(4.6-0.6), h_tpb2vac_thin0p00000001->GetXaxis()->FindBin(5.2));
    double int2 = h_tpb2vac_thin0p00001->Integral(h_tpb2vac_thin0p00001->GetXaxis()->FindBin(4.6-0.6), h_tpb2vac_thin0p00001->GetXaxis()->FindBin(5.2));
    double int3 = h_tpb2vac_thin0p001->Integral(h_tpb2vac_thin0p001->GetXaxis()->FindBin(4.6-0.6), h_tpb2vac_thin0p001->GetXaxis()->FindBin(5.2));
    double int4 = h_tpb2vac_thin0p05->Integral(h_tpb2vac_thin0p05->GetXaxis()->FindBin(4.6-0.6), h_tpb2vac_thin0p05->GetXaxis()->FindBin(5.2));
    double int5 = h_tpb2vac_thin0p1->Integral(h_tpb2vac_thin0p1->GetXaxis()->FindBin(4.6-0.6), h_tpb2vac_thin0p1->GetXaxis()->FindBin(4.6+0.6));

    cout<<int1<<", "<<int2<<", "<<int3<<", "<<int4<<","<<int5<<" "<<total<<endl;

    double scale0 = h_tpb2vac_thin0p00000001->Integral();
    h_tpb2vac_thin0p00001->Scale(scale0/h_tpb2vac_thin0p00001->Integral());
    h_tpb2vac_thin0p001->Scale(scale0/h_tpb2vac_thin0p001->Integral());
    h_tpb2vac_thin0p05->Scale(scale0/h_tpb2vac_thin0p05->Integral());
    h_tpb2vac_thin0p1->Scale(scale0/h_tpb2vac_thin0p1->Integral());
///// vac2tpb
    h_vac2tpb_thin0p00000001 = (TH1F*)ff1->Get("H_E_vac2tpb");
    h_vac2tpb_thin0p00001    = (TH1F*)ff2->Get("H_E_vac2tpb");
    h_vac2tpb_thin0p001      = (TH1F*)ff3->Get("H_E_vac2tpb");
    h_vac2tpb_thin0p05       = (TH1F*)ff4->Get("H_E_vac2tpb");
    h_vac2tpb_thin0p1        = (TH1F*)ff5->Get("H_E_vac2tpb");

    h_vac2tpb_thin0p00000001->SetLineColor(kRed);
    h_vac2tpb_thin0p00001->SetLineColor(kBlue);
    h_vac2tpb_thin0p001->SetLineColor(kBlack);
    h_vac2tpb_thin0p05->SetLineColor(kGreen+2);
    h_vac2tpb_thin0p1->SetLineColor(kOrange+1);

    double scale0 = h_vac2tpb_thin0p00000001->Integral();
    h_vac2tpb_thin0p00001->Scale(scale0/h_vac2tpb_thin0p00001->Integral());
    h_vac2tpb_thin0p001->Scale(scale0/h_vac2tpb_thin0p001->Integral());
    h_vac2tpb_thin0p05->Scale(scale0/h_vac2tpb_thin0p05->Integral());
    h_vac2tpb_thin0p1->Scale(scale0/h_vac2tpb_thin0p1->Integral());

    // delta E loss from E0 to E(at av2tpb), loss in AV
    h_deltaEp_av_thin0p00000001 = (TH1F*)ff1->Get("H_deltaEp_av");
    h_deltaEp_av_thin0p00001= (TH1F*)ff2->Get("H_deltaEp_av");
    h_deltaEp_av_thin0p001 = (TH1F*)ff3->Get("H_deltaEp_av");
    h_deltaEp_av_thin0p05 = (TH1F*)ff4->Get("H_deltaEp_av");
    h_deltaEp_av_thin0p1 = (TH1F*)ff5->Get("H_deltaEp_av");

    h_deltaEp_av_thin0p00000001->SetLineColor(kRed);
    h_deltaEp_av_thin0p00001->SetLineColor(kBlue);
    h_deltaEp_av_thin0p001->SetLineColor(kBlack);
    h_deltaEp_av_thin0p05->SetLineColor(kGreen+2);
    h_deltaEp_av_thin0p1->SetLineColor(kOrange+1);

    double scale = h_deltaEp_av_thin0p00000001->Integral();
    h_deltaEp_av_thin0p00001->Scale(scale/h_deltaEp_av_thin0p00001->Integral());
    h_deltaEp_av_thin0p001->Scale(scale/h_deltaEp_av_thin0p001->Integral());
    h_deltaEp_av_thin0p05->Scale(scale/h_deltaEp_av_thin0p05->Integral());
    h_deltaEp_av_thin0p1->Scale(scale/h_deltaEp_av_thin0p1->Integral());

    // E at av2tpb
    h_av2tpb_thin0p00000001 = (TH1F*)ff1->Get("H_E_av2tpb");
    h_av2tpb_thin0p00001= (TH1F*)ff2->Get("H_E_av2tpb");
    h_av2tpb_thin0p001 = (TH1F*)ff3->Get("H_E_av2tpb");
    h_av2tpb_thin0p05 = (TH1F*)ff4->Get("H_E_av2tpb");
    h_av2tpb_thin0p1 = (TH1F*)ff5->Get("H_E_av2tpb");

    h_av2tpb_thin0p00000001->SetLineColor(kRed);
    h_av2tpb_thin0p00001->SetLineColor(kBlue);
    h_av2tpb_thin0p001->SetLineColor(kBlack);
    h_av2tpb_thin0p05->SetLineColor(kGreen+2);
    h_av2tpb_thin0p1->SetLineColor(kOrange+1);

    double scale1 = h_av2tpb_thin0p00000001->Integral();
    h_av2tpb_thin0p00001->Scale(scale1/h_av2tpb_thin0p00001->Integral());
    h_av2tpb_thin0p001->Scale(scale1/h_av2tpb_thin0p001->Integral());
    h_av2tpb_thin0p05->Scale(scale1/h_av2tpb_thin0p05->Integral());
    h_av2tpb_thin0p1->Scale(scale1/h_av2tpb_thin0p1->Integral());

    //////
    //
    // delta E tpb to vac, E loss in tpb layer 1
    h_deltaEp_tpb1_thin0p00000001 = (TH1F*)ff1->Get("H_deltaEp_tpb1");
    h_deltaEp_tpb1_thin0p00001    = (TH1F*)ff2->Get("H_deltaEp_tpb1");
    h_deltaEp_tpb1_thin0p001      = (TH1F*)ff3->Get("H_deltaEp_tpb1");
    h_deltaEp_tpb1_thin0p05       = (TH1F*)ff4->Get("H_deltaEp_tpb1");
    h_deltaEp_tpb1_thin0p1        = (TH1F*)ff5->Get("H_deltaEp_tpb1");

    h_deltaEp_tpb1_thin0p00000001->SetLineColor(kRed);
    h_deltaEp_tpb1_thin0p00001->SetLineColor(kBlue);
    h_deltaEp_tpb1_thin0p001->SetLineColor(kBlack);
    h_deltaEp_tpb1_thin0p05->SetLineColor(kGreen+2);
    h_deltaEp_tpb1_thin0p1->SetLineColor(kOrange+1);

    double scale2 = h_deltaEp_tpb1_thin0p00000001->Integral();
    h_deltaEp_tpb1_thin0p00001->Scale(scale2/h_deltaEp_tpb1_thin0p00001->Integral());
    h_deltaEp_tpb1_thin0p001->Scale(scale2/h_deltaEp_tpb1_thin0p001->Integral());
    h_deltaEp_tpb1_thin0p05->Scale(scale2/h_deltaEp_tpb1_thin0p05->Integral());
    h_deltaEp_tpb1_thin0p1->Scale(scale2/h_deltaEp_tpb1_thin0p1->Integral());

    //////
    //
    // delta E from vac to tpb2, or E loss in tpb layer 2
    h_deltaEp_tpb2_thin0p00000001 = (TH1F*)ff1->Get("H_deltaEp_tpb2");
    h_deltaEp_tpb2_thin0p00001    = (TH1F*)ff2->Get("H_deltaEp_tpb2");
    h_deltaEp_tpb2_thin0p001      = (TH1F*)ff3->Get("H_deltaEp_tpb2");
    h_deltaEp_tpb2_thin0p05       = (TH1F*)ff4->Get("H_deltaEp_tpb2");
    h_deltaEp_tpb2_thin0p1        = (TH1F*)ff5->Get("H_deltaEp_tpb2");

    h_deltaEp_tpb2_thin0p00000001->SetLineColor(kRed);
    h_deltaEp_tpb2_thin0p00001->SetLineColor(kBlue);
    h_deltaEp_tpb2_thin0p001->SetLineColor(kBlack);
    h_deltaEp_tpb2_thin0p05->SetLineColor(kGreen+2);
    h_deltaEp_tpb2_thin0p1->SetLineColor(kOrange+1);

    double scale3 = h_deltaEp_tpb2_thin0p00000001->Integral();
    h_deltaEp_tpb2_thin0p00001->Scale(scale3/h_deltaEp_tpb2_thin0p00001->Integral());
    h_deltaEp_tpb2_thin0p001->Scale(scale3/h_deltaEp_tpb2_thin0p001->Integral());
    h_deltaEp_tpb2_thin0p05->Scale(scale3/h_deltaEp_tpb2_thin0p05->Integral());
    h_deltaEp_tpb2_thin0p1->Scale(scale3/h_deltaEp_tpb2_thin0p1->Integral());

    //////
    // delta E from tpb2vac to vac2tpb, or E loss in vac
    h_deltaEp_vac_thin0p00000001 = (TH1F*)ff1->Get("H_deltaEp_vac");
    h_deltaEp_vac_thin0p00001    = (TH1F*)ff2->Get("H_deltaEp_vac");
    h_deltaEp_vac_thin0p001      = (TH1F*)ff3->Get("H_deltaEp_vac");
    h_deltaEp_vac_thin0p05       = (TH1F*)ff4->Get("H_deltaEp_vac");
    h_deltaEp_vac_thin0p1        = (TH1F*)ff5->Get("H_deltaEp_vac");

    h_deltaEp_vac_thin0p00000001->SetLineColor(kRed);
    h_deltaEp_vac_thin0p00001->SetLineColor(kBlue);
    h_deltaEp_vac_thin0p001->SetLineColor(kBlack);
    h_deltaEp_vac_thin0p05->SetLineColor(kGreen+2);
    h_deltaEp_vac_thin0p1->SetLineColor(kOrange+1);

    double scale4 = h_deltaEp_vac_thin0p00000001->Integral();
    h_deltaEp_vac_thin0p00001->Scale(scale4/h_deltaEp_vac_thin0p00001->Integral());
    h_deltaEp_vac_thin0p001->Scale(scale4/h_deltaEp_vac_thin0p001->Integral());
    h_deltaEp_vac_thin0p05->Scale(scale4/h_deltaEp_vac_thin0p05->Integral());
    h_deltaEp_vac_thin0p1->Scale(scale4/h_deltaEp_vac_thin0p1->Integral());

    ///////
    // delta E from tpb2av to final, or E loss in 2nd layer AV
    h_deltaEp_final_thin0p00000001 = (TH1F*)ff1->Get("H_deltaEp_final");
    h_deltaEp_final_thin0p00001    = (TH1F*)ff2->Get("H_deltaEp_final");
    h_deltaEp_final_thin0p001      = (TH1F*)ff3->Get("H_deltaEp_final");
    h_deltaEp_final_thin0p05       = (TH1F*)ff4->Get("H_deltaEp_final");
    h_deltaEp_final_thin0p1        = (TH1F*)ff5->Get("H_deltaEp_final");

    h_deltaEp_final_thin0p00000001->SetLineColor(kRed);
    h_deltaEp_final_thin0p00001->SetLineColor(kBlue);
    h_deltaEp_final_thin0p001->SetLineColor(kBlack);
    h_deltaEp_final_thin0p05->SetLineColor(kGreen+2);
    h_deltaEp_final_thin0p1->SetLineColor(kOrange+1);

    double scale5 = h_deltaEp_final_thin0p00000001->Integral();
    h_deltaEp_final_thin0p00001->Scale(scale5/h_deltaEp_final_thin0p00001->Integral());
    h_deltaEp_final_thin0p001->Scale(scale5/h_deltaEp_final_thin0p001->Integral());
    h_deltaEp_final_thin0p05->Scale(scale5/h_deltaEp_final_thin0p05->Integral());
    h_deltaEp_final_thin0p1->Scale(scale5/h_deltaEp_final_thin0p1->Integral());


///
//
//  all E 
    h_allStepE_thin0p00000001 = (TH1F*)ff1->Get("H_energy_all");
    h_allStepE_thin0p00001    = (TH1F*)ff2->Get("H_energy_all");
    h_allStepE_thin0p001      = (TH1F*)ff3->Get("H_energy_all");
    h_allStepE_thin0p05       = (TH1F*)ff4->Get("H_energy_all");
    h_allStepE_thin0p1        = (TH1F*)ff5->Get("H_energy_all");

    h_allStepE_thin0p00000001->SetLineColor(kRed);
    h_allStepE_thin0p00001->SetLineColor(kBlue);
    h_allStepE_thin0p001->SetLineColor(kBlack);
    h_allStepE_thin0p05->SetLineColor(kGreen+2);
    h_allStepE_thin0p1->SetLineColor(kOrange+1);

    double scaleEtot = h_allStepE_thin0p00000001->Integral();
    h_allStepE_thin0p00001->Scale(scaleEtot/h_allStepE_thin0p00001->Integral());
    h_allStepE_thin0p001->Scale(scaleEtot/h_allStepE_thin0p001->Integral());
    h_allStepE_thin0p05->Scale(scaleEtot/h_allStepE_thin0p05->Integral());
    h_allStepE_thin0p1->Scale(scaleEtot/h_allStepE_thin0p1->Integral());

    //h_deltaEp_av->Draw();
    // h_deltaEp_tpb1->Draw("same");

    TCanvas *c1 = new TCanvas("c1","delta E0 - Eav2tpb",800,600);
    h_deltaEp_av_thin0p00000001->Draw();
    h_deltaEp_av_thin0p00001->Draw("same");
    h_deltaEp_av_thin0p001->Draw("same");
    h_deltaEp_av_thin0p05->Draw("same");
    h_deltaEp_av_thin0p1->Draw("same");

    TCanvas *c2 = new TCanvas("c2","Eav2tpb",800,600);
    h_av2tpb_thin0p00000001->Draw();
    h_av2tpb_thin0p00001->Draw("same");
    h_av2tpb_thin0p001->Draw("same");
    h_av2tpb_thin0p05->Draw("same");
    h_av2tpb_thin0p1->Draw("same");

    TCanvas *c3 = new TCanvas("c3","delta Etpb - Evac (1st layer)",800,600);
    h_deltaEp_tpb1_thin0p00000001->GetXaxis()->SetTitle("MeV");h_deltaEp_tpb1_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_deltaEp_tpb1_thin0p00000001->Draw();
    h_deltaEp_tpb1_thin0p00001->Draw("same");
    h_deltaEp_tpb1_thin0p001->Draw("same");
    h_deltaEp_tpb1_thin0p05->Draw("same");
    h_deltaEp_tpb1_thin0p1->Draw("same");

    TCanvas *c4 = new TCanvas("c4","delta Evac - Etpb (2nd layer)",800,600);
    h_deltaEp_tpb2_thin0p00000001->GetXaxis()->SetTitle("MeV");h_deltaEp_tpb2_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_deltaEp_tpb2_thin0p00000001->Draw();
    h_deltaEp_tpb2_thin0p00001->Draw("same");
    h_deltaEp_tpb2_thin0p001->Draw("same");
    h_deltaEp_tpb2_thin0p05->Draw("same");
    h_deltaEp_tpb2_thin0p1->Draw("same");

    TCanvas *c5 = new TCanvas("c5","deltaE in vac",800,600);
    h_deltaEp_vac_thin0p00000001->GetXaxis()->SetTitle("MeV");h_deltaEp_vac_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_deltaEp_vac_thin0p00000001->Draw();
    h_deltaEp_vac_thin0p00001->Draw("same");
    h_deltaEp_vac_thin0p001->Draw("same");
    h_deltaEp_vac_thin0p05->Draw("same");
    h_deltaEp_vac_thin0p1->Draw("same");

    TCanvas *c6 = new TCanvas("c6","delta Etpb2av - Ef, or E loss in 2nd AV layer",800,600);
    h_deltaEp_final_thin0p00000001->GetXaxis()->SetTitle("MeV");h_deltaEp_final_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_deltaEp_final_thin0p00000001->Draw();
    h_deltaEp_final_thin0p00001->Draw("same");
    h_deltaEp_final_thin0p001->Draw("same");
    h_deltaEp_final_thin0p05->Draw("same");
    h_deltaEp_final_thin0p1->Draw("same");

    TCanvas *c7 = new TCanvas("c7","E at tpb2vac",800,600);

    h_tpb2vac_thin0p00000001->GetXaxis()->SetTitle("MeV");h_tpb2vac_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_tpb2vac_thin0p00000001->Draw();
    h_tpb2vac_thin0p00001->Draw("same");
    h_tpb2vac_thin0p001->Draw("same");
    h_tpb2vac_thin0p05->Draw("same");
    h_tpb2vac_thin0p1->Draw("same");

    TCanvas *c8 = new TCanvas("c8","E at tpb2vac",800,600);
    h_vac2tpb_thin0p00000001->GetXaxis()->SetTitle("MeV");h_vac2tpb_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_vac2tpb_thin0p00000001->Draw();
    h_vac2tpb_thin0p00001->Draw("same");
    h_vac2tpb_thin0p001->Draw("same");
    h_vac2tpb_thin0p05->Draw("same");
    h_vac2tpb_thin0p1->Draw("same");

    TCanvas *c9 = new TCanvas("c9","E at tpb2vac",800,600);
    h_allStepE_thin0p00000001->GetXaxis()->SetTitle("MeV");h_allStepE_thin0p00000001->GetYaxis()->SetTitle("Scaled");
    h_allStepE_thin0p00000001->Draw();
    h_allStepE_thin0p00001->Draw("same");
    h_allStepE_thin0p001->Draw("same");
    h_allStepE_thin0p05->Draw("same");
    h_allStepE_thin0p1->Draw("same");


}
