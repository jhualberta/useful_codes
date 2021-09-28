double scaleNue = 1./1700;
double scaleNumu = 1./9600;
double scaleTime = 0.96;

TH1F *hPee(TH2F* hcosVsE, bool checkOsc, bool checkNumu)
{
 //PSelmaa
  double Etot[155] ={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15,15.1,15.2,15.3,15.4,15.5};

  double pee_bs05op_b8[155] = {0.555112,0.552053,0.54893,0.545744,0.542496,0.539188,0.535822,0.532399,0.528922,0.525395,0.521819,0.518197,0.514534,0.510832,0.507095,0.503328,0.499535,0.495719,0.491885,0.488038,0.484183,0.480323,0.476464,0.47261,0.468767,0.464937,0.461127,0.457339,0.453579,0.44985,0.446157,0.442502,0.438889,0.435322,0.431802,0.428334,0.424919,0.42156,0.418258,0.415016,0.411834,0.408714,0.405657,0.402663,0.399734,0.39687,0.394071,0.391336,0.388666,0.386061,0.38352,0.381043,0.378629,0.376276,0.373986,0.371755,0.369585,0.367473,0.365418,0.363419,0.361476,0.359587,0.35775,0.355965,0.35423,0.352544,0.350906,0.349314,0.347768,0.346266,0.344807,0.343389,0.342011,0.340673,0.339374,0.338111,0.336884,0.335692,0.334534,0.333408,0.332315,0.331252,0.330219,0.329215,0.328238,0.327289,0.326367,0.32547,0.324597,0.323749,0.322923,0.32212,0.321339,0.320578,0.319839,0.319118,0.318417,0.317735,0.31707,0.316423,0.315792,0.315178,0.31458,0.313997,0.313428,0.312874,0.312334,0.311808,0.311295,0.310794,0.310306,0.309829,0.309364,0.308911,0.308468,0.308036,0.307614,0.307202,0.306799,0.306406,0.306023,0.305648,0.305281,0.304923,0.304573,0.304231,0.303896,0.303569,0.303249,0.302936,0.30263,0.30233,0.302037,0.30175,0.301469,0.301194,0.300925,0.300662,0.300404,0.300151,0.299903,0.299661,0.299423,0.29919,0.298961,0.298738,0.298518,0.298303,0.298092,0.297885,0.297682,0.297483,0.297288,0.297096,0.296908};

  TH1F *hcosSun_pee = new TH1F("hcosSun_pee","",40,-1,1);
  TH1F *hcosSun_peeEtrue = new TH1F("hcosSun_peeEtrue","use Etrue to get Pee",40,-1,1);
  TH1F *hcosSun_noOsci= new TH1F("hcosSun_noOsci","",40,-1,1);

  if(!checkOsc)
  {
    for(int i = 0; i<40; i++)
     {
       for(int j = 0; j<10; j++)
       {
         int counts = hcosVsE->GetBinContent(i+1, j+1);
         double cosThetaToSun = -1+i*(2./40);
         hcosSun_noOsci->Fill(cosThetaToSun, counts*scaleNue*scaleTime);
       }
     }
    cout<<"No oscillation "<<hcosSun_noOsci->Integral()<<endl;
    return hcosSun_noOsci;
  }
  else {
    for(int i = 0; i<40; i++)
    {
      for(int j = 0; j<10; j++) 
      {
        int counts = hcosVsE->GetBinContent(i+1, j+1);         
        int index = 0;
        double weight = 0;
        double energyfit = 5+j*1;
        double cosThetaToSun = -1+i*(2./40);
        index = TMath::Nint((energyfit/0.1));
        double peeTrue = pee_bs05op_b8[index-1];
        weight = counts*scaleNue*scaleTime*peeTrue;
        if(checkNumu) weight = counts*scaleNumu*scaleTime*(1 - peeTrue);
        //cout<<counts<<" "<<cosThetaToSun<<" "<<energyfit<<" "<<peeTrue<<" "<<weight<<endl;
        hcosSun_pee->Fill(cosThetaToSun, weight);
      }
    }
    cout<<"Erecon: "<<hcosSun_pee->Integral()<<", Etrue: "<<hcosSun_pee->Integral()<<endl;
    hcosSun_pee->SetLineColor(kRed);
    return hcosSun_pee;
  }

}
  
void scalePeeMC_tmva()
{

  bool checkNumu = 0;
  TFile *fMC = new TFile("tmva_mpw_MCsolarNue_E5to15_whole_9params.root");
  TFile *fMC_numu = new TFile("tmva_mpw_MCsolarNumu_E5to15_whole_9params.root");

  TString ss = "BDT"; //  "MLP"  
//----- smearing files -------
/// analysis_MPWdata_smearDirDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearDirUp_E5to15_FV5p5m.root
/// analysis_MPWdata_smearEresolDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearEresolUp_E5to15_FV5p5m.root
/// analysis_MPWdata_smearEscaleDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearEscaleUp_E5to15_FV5p5m.root
/// analysis_MPWdata_smearPosScaleDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearPosScaleUp_E5to15_FV5p5m.root
/// analysis_MPWdata_smearXshiftDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearXshiftUp.root
/// analysis_MPWdata_smearYshiftDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearYshiftUp_E5to15_FV5p5m.root
/// analysis_MPWdata_smearZshiftDown_E5to15_FV5p5m.root
/// analysis_MPWdata_smearZshiftUp_E5to15_FV5p5m.root


  TString histName1 = "h"+ss+"SigCosThetaToSunVsE";
  TString histName2 = "h"+ss+"SigCosThetaToSunVsE";
  TH2F *hcosNue = (TH2F*)fMC->Get(histName1);
  TH2F *hcosNumu = (TH2F*)fMC_numu->Get(histName2);

//  TH2F *hcosNue = (TH2F*)fMC->Get("hMLPSigCosThetaToSunVsE");
//  TH2F *hcosNumu = (TH2F*)fMC_numu->Get("hMLPSigCosThetaToSunVsE");

  TH1F *hcosSunNue_pee;// = new TH1F("hcosSunNue_pee","",40,-1,1);
  TH1F *hcosSunNumu_pee;// = new TH1F("hcosSunNumu_pee","",40,-1,1);
  TH1F *hcosSun_noOsci;

  TString fileName = "saveMCtmva"+ss+"_nue_numu_pdfs_9params.root";
  TFile *fnew = new TFile(fileName,"recreate");
  // TFile *fnew = new TFile("saveMCtmvaMLP_nue_numu_pdfs.root","recreate");

  hcosSunNue_pee = hPee(hcosNue, 1, 0);
  hcosSunNumu_pee = hPee(hcosNumu, 1, 1);
  hcosSun_noOsci = hPee(hcosNue, 0, 0);

  hcosSunNue_pee->SetTitle("hcosSunNue_pee");
  hcosSunNumu_pee->SetTitle("hcosSunNumu_pee");
  hcosSunNue_pee->SetName("hcosSunNue_pee");
  hcosSunNumu_pee->SetName("hcosSunNumu_pee");

  TH1F *hcosSun_peeEtrue_sum = (TH1F*)hcosSunNue_pee->Clone();
  hcosSun_peeEtrue_sum->Add(hcosSunNumu_pee);
  hcosSun_peeEtrue_sum->SetName("hcosSunSumNueNumu_peeEtrue");

  fnew->cd();
  hcosSun_noOsci->Write();
  hcosSunNue_pee->Write();
  hcosSunNumu_pee->Write();
  hcosSun_peeEtrue_sum->Write();
}
