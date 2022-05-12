{

  TFile *inputfile = new TFile("ensembleTMVA_output5000evts_poisson.root");
  TH1F *hRatio_sig = new TH1F("hRatio_sig", "signal ratios", 200,0,2);
  TH1F *hRatio_bkg = new TH1F("hRatio_bkg", "background ratios", 200,0,2);
 
  for(int i = 0; i<5000;i++)
  {
   TString number;
   number = number.Itoa(i,10);

   TString hSigName = "hcosTheta_sig"+number;
   TString hBkgName = "hcosTheta_bkg"+number;

   TH1F *hSig = (TH1F*)inputfile->Get(hSigName);
   TH1F *hBkg = (TH1F*)inputfile->Get(hBkgName);

   double Nsig_true = hSig->Integral();
   double Nbkg_true = hBkg->Integral();
   double Nsig_cut, Nbkg_cut;
    
   TTree *Tdata = inputfile->Get("T2_"+number);
   float scaleLogL, Gtest, Utest, zfactor, itr, beta14, cosThetaToSun;
   float signal;
   Tdata->SetBranchAddress("scaleLogL", &scaleLogL);
   Tdata->SetBranchAddress("Gtest", &Gtest);
   Tdata->SetBranchAddress("Utest", &Utest);
   Tdata->SetBranchAddress("zfactor", &zfactor);
   Tdata->SetBranchAddress("itr", &itr);
   Tdata->SetBranchAddress("beta14", &beta14);
   Tdata->SetBranchAddress("cosThetaToSun", &cosThetaToSun);
   Tdata->SetBranchAddress("signal", &signal);
   TH1F *hsig = new TH1F("hsig"+number,"",40,-1,1);
   TH1F *htot = new TH1F("htot"+number,"",40,-1,1);
// default cuts
   Tdata->Project("hsig"+number, "cosThetaToSun", "scaleLogL>10.0 && -11<zfactor && zfactor<1 && 0<Gtest && Gtest<1.9 && Utest<0.95 && itr>0.55 && -0.12<beta14 && beta14<0.95");

/// suggested by Brian (FOM10)
//
//   Tdata->Project("hsig"+number, "cosThetaToSun", "scaleLogL>9.9 && itr>0.55 && -0.12<beta14 && beta14<0.95 && 0<Gtest && Gtest<1.9 && Utest>0.5");

   Tdata->Project("htot"+number, "cosThetaToSun", "");

   double Nsig_cut = hsig->Integral();
   double Nbkg_cut = htot->Integral() - Nsig_cut;

   // cout<<"#"<<i<<" cuts "<<Nsig_cut<<" "<<Nbkg_cut<<" "<<Nsig_true<<" "<<Nbkg_true<<endl;
   hRatio_sig->Fill(Nsig_cut/Nsig_true);
   hRatio_bkg->Fill(Nbkg_cut/Nbkg_true);
   // delete hSig; delete hBkg; delete Tdata;
  }
  TFile *fnew = new TFile("outputTMVAratiosSimpleCuts_test1.root","recreate");
  fnew->cd();
  hRatio_sig->Write();
  hRatio_bkg->Write();

  hRatio_sig->Draw();
  hRatio_bkg->Draw("sames");

}
