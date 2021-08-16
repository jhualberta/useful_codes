double scaleNue = 1./1700;
double scaleNumu = 1./9600;
double scaleTime = 0.96;
double scaleDataClean = (1-1.7/100);

// scale Nue, Numu MC to detection counts, 5-15 MeV
TH1F *hPee(TTree* T, bool checkNumu)
{
  UInt_t nhits;
  double energy = 0, klDiv=-1000, itr = 0, posx = 0, posy = 0, posz = 0, posRad = 0, dirx = 0, diry = 0, dirz = 0, beta14 =0, Gtest= 0, Utest = 0, scaleLogL = 0, cosThetaToSun = 0;
  double medianProbHit = 0, medianProb = 0, medianDevHit=0, medianDev=0;
  double udotR = 0, zfactor = 0;
  double posxmc = 0, posymc = 0, poszmc = 0, dirxmc = 0, dirymc = 0, dirzmc = 0, energymc = 0;
  UInt_t eventGTID = 0, runNumber = 0;
  //  std::vector<double> *gcosPMT0 = new std::vector<double>;
  //std::vector<double> *gCosPMTfit = new std::vector<double>;
  //std::vector<double> *gtRes = new std::vector<double>;
  T->SetBranchAddress("energy", &energy);////!!!! original energy!!!
  T->SetBranchAddress("nhits", &nhits);
  T->SetBranchAddress("itr", &itr);
  T->SetBranchAddress("posx", &posx);
  T->SetBranchAddress("posy", &posy);
  T->SetBranchAddress("posz", &posz);
  T->SetBranchAddress("dirx", &dirx);
  T->SetBranchAddress("diry", &diry);
  T->SetBranchAddress("dirz", &dirz);
  T->SetBranchAddress("beta14", &beta14);
  T->SetBranchAddress("cosThetaToSun", &cosThetaToSun);
  T->SetBranchAddress("Gtest", &Gtest);
  T->SetBranchAddress("Utest", &Utest);
  T->SetBranchAddress("scaleLogL", &scaleLogL);
  T->SetBranchAddress("runNumber",&runNumber);
  T->SetBranchAddress("eventGTID",&eventGTID);
  T->SetBranchAddress("zfactor", &zfactor);

  T->SetBranchAddress("posxmc", &posxmc);
  T->SetBranchAddress("posymc", &posymc);
  T->SetBranchAddress("posymc", &poszmc);
  T->SetBranchAddress("dirxmc", &dirxmc);
  T->SetBranchAddress("dirymc", &dirymc);
  T->SetBranchAddress("dirzmc", &dirzmc);
  T->SetBranchAddress("energymc", &energymc);

 //PSelmaa
  double Etot[155] ={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15,15.1,15.2,15.3,15.4,15.5};

  double pee_bs05op_b8[155] = {0.555112,0.552053,0.54893,0.545744,0.542496,0.539188,0.535822,0.532399,0.528922,0.525395,0.521819,0.518197,0.514534,0.510832,0.507095,0.503328,0.499535,0.495719,0.491885,0.488038,0.484183,0.480323,0.476464,0.47261,0.468767,0.464937,0.461127,0.457339,0.453579,0.44985,0.446157,0.442502,0.438889,0.435322,0.431802,0.428334,0.424919,0.42156,0.418258,0.415016,0.411834,0.408714,0.405657,0.402663,0.399734,0.39687,0.394071,0.391336,0.388666,0.386061,0.38352,0.381043,0.378629,0.376276,0.373986,0.371755,0.369585,0.367473,0.365418,0.363419,0.361476,0.359587,0.35775,0.355965,0.35423,0.352544,0.350906,0.349314,0.347768,0.346266,0.344807,0.343389,0.342011,0.340673,0.339374,0.338111,0.336884,0.335692,0.334534,0.333408,0.332315,0.331252,0.330219,0.329215,0.328238,0.327289,0.326367,0.32547,0.324597,0.323749,0.322923,0.32212,0.321339,0.320578,0.319839,0.319118,0.318417,0.317735,0.31707,0.316423,0.315792,0.315178,0.31458,0.313997,0.313428,0.312874,0.312334,0.311808,0.311295,0.310794,0.310306,0.309829,0.309364,0.308911,0.308468,0.308036,0.307614,0.307202,0.306799,0.306406,0.306023,0.305648,0.305281,0.304923,0.304573,0.304231,0.303896,0.303569,0.303249,0.302936,0.30263,0.30233,0.302037,0.30175,0.301469,0.301194,0.300925,0.300662,0.300404,0.300151,0.299903,0.299661,0.299423,0.29919,0.298961,0.298738,0.298518,0.298303,0.298092,0.297885,0.297682,0.297483,0.297288,0.297096,0.296908};

  TH1F *hcosSun_pee = new TH1F("hcosSun_pee","",40,-1,1);
  TH1F *hcosSun_peeEtrue = new TH1F("hcosSun_peeEtrue","use Etrue to get Pee",40,-1,1);
  TH1F *hcosSun_noOsci= new TH1F("hcosSunNue_noOsci","",40,-1,1);

  for(int i =0;i<T->GetEntries();i++)
  {
    //gCosPMTfit->clear();
    T->GetEntry(i);
    if(sqrt(posx**2+posy**2+(posz-108)**2)<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && energy<=15 && energy>=5 && zfactor>-11 && zfactor<1)
    {
// Pee Erecon
       // if( energy<15.5 && energy>0.1 )
       // {
       //    int index = 0;
       //    double weight = 0;
       //    index = TMath::Nint((energy/0.1));
       //    double peeTrue = pee_bs05op_b8[index-1];
       //    weight = scaleNue*scaleTime*peeTrue;
       //    if(checkNumu) weight = scaleNumu*scaleTime*(1 - peeTrue);
       //    hcosSun_pee->Fill(cosThetaToSun, weight);
       //    //cout<<energy<<", "<<index<<", "<<Etot[index-1]<<", "<<pee_bs05op_b8[index-1]<<endl;
       // }
// Pee Etrue
       if( energymc<15.5 && energymc>0.1 )
       {
          int index = 0;
          double weight = 0;
          index = TMath::Nint((energymc/0.1));
          double peeTrue = pee_bs05op_b8[index-1];
          weight = scaleNue*scaleTime*scaleDataClean*peeTrue;
	  if(checkNumu) { weight = scaleNumu*scaleTime*scaleDataClean*(1 - peeTrue);}
          hcosSun_peeEtrue->Fill(cosThetaToSun, weight);
          //cout<<energy<<", "<<index<<", "<<Etot[index-1]<<", "<<pee_bs05op_b8[index-1]<<endl;
        }
     }
   }

  cout<<"Erecon: "<<hcosSun_pee->Integral()<<", Etrue: "<<hcosSun_peeEtrue->Integral()<<endl;
//  cout<<"No oscillation "<<hcosSun_noOsci->Integral()<<endl;

  return hcosSun_peeEtrue;
}
 
// for no oscillation
TH1F *hPee(TTree* T)
{
  UInt_t nhits;
  double energy = 0, klDiv=-1000, itr = 0, posx = 0, posy = 0, posz = 0, posRad = 0, dirx = 0, diry = 0, dirz = 0, beta14 =0, Gtest= 0, Utest = 0, scaleLogL = 0, cosThetaToSun = 0;
  double medianProbHit = 0, medianProb = 0, medianDevHit=0, medianDev=0;
  double udotR = 0, zfactor = 0;
  double posxmc = 0, posymc = 0, poszmc = 0, dirxmc = 0, dirymc = 0, dirzmc = 0, energymc = 0;
  UInt_t eventGTID = 0, runNumber = 0;
  //  std::vector<double> *gcosPMT0 = new std::vector<double>;
  //std::vector<double> *gCosPMTfit = new std::vector<double>;
  //std::vector<double> *gtRes = new std::vector<double>;
  T->SetBranchAddress("energy", &energy);////!!!! original energy!!!
  T->SetBranchAddress("nhits", &nhits);
  T->SetBranchAddress("itr", &itr);
  T->SetBranchAddress("posx", &posx);
  T->SetBranchAddress("posy", &posy);
  T->SetBranchAddress("posz", &posz);
  T->SetBranchAddress("dirx", &dirx);
  T->SetBranchAddress("diry", &diry);
  T->SetBranchAddress("dirz", &dirz);
  T->SetBranchAddress("beta14", &beta14);
  T->SetBranchAddress("cosThetaToSun", &cosThetaToSun);
  T->SetBranchAddress("Gtest", &Gtest);
  T->SetBranchAddress("Utest", &Utest);
  T->SetBranchAddress("scaleLogL", &scaleLogL);
  T->SetBranchAddress("runNumber",&runNumber);
  T->SetBranchAddress("eventGTID",&eventGTID);
  T->SetBranchAddress("zfactor", &zfactor);

  T->SetBranchAddress("posxmc", &posxmc);
  T->SetBranchAddress("posymc", &posymc);
  T->SetBranchAddress("posymc", &poszmc);
  T->SetBranchAddress("dirxmc", &dirxmc);
  T->SetBranchAddress("dirymc", &dirymc);
  T->SetBranchAddress("dirzmc", &dirzmc);
  T->SetBranchAddress("energymc", &energymc);

 //PSelmaa
  double Etot[155] ={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15,15.1,15.2,15.3,15.4,15.5};

  double pee_bs05op_b8[155] = {0.555112,0.552053,0.54893,0.545744,0.542496,0.539188,0.535822,0.532399,0.528922,0.525395,0.521819,0.518197,0.514534,0.510832,0.507095,0.503328,0.499535,0.495719,0.491885,0.488038,0.484183,0.480323,0.476464,0.47261,0.468767,0.464937,0.461127,0.457339,0.453579,0.44985,0.446157,0.442502,0.438889,0.435322,0.431802,0.428334,0.424919,0.42156,0.418258,0.415016,0.411834,0.408714,0.405657,0.402663,0.399734,0.39687,0.394071,0.391336,0.388666,0.386061,0.38352,0.381043,0.378629,0.376276,0.373986,0.371755,0.369585,0.367473,0.365418,0.363419,0.361476,0.359587,0.35775,0.355965,0.35423,0.352544,0.350906,0.349314,0.347768,0.346266,0.344807,0.343389,0.342011,0.340673,0.339374,0.338111,0.336884,0.335692,0.334534,0.333408,0.332315,0.331252,0.330219,0.329215,0.328238,0.327289,0.326367,0.32547,0.324597,0.323749,0.322923,0.32212,0.321339,0.320578,0.319839,0.319118,0.318417,0.317735,0.31707,0.316423,0.315792,0.315178,0.31458,0.313997,0.313428,0.312874,0.312334,0.311808,0.311295,0.310794,0.310306,0.309829,0.309364,0.308911,0.308468,0.308036,0.307614,0.307202,0.306799,0.306406,0.306023,0.305648,0.305281,0.304923,0.304573,0.304231,0.303896,0.303569,0.303249,0.302936,0.30263,0.30233,0.302037,0.30175,0.301469,0.301194,0.300925,0.300662,0.300404,0.300151,0.299903,0.299661,0.299423,0.29919,0.298961,0.298738,0.298518,0.298303,0.298092,0.297885,0.297682,0.297483,0.297288,0.297096,0.296908};

  TH1F *hcosSun_pee = new TH1F("hcosSun_pee","",40,-1,1);
  TH1F *hcosSun_peeEtrue = new TH1F("hcosSun_peeEtrue","use Etrue to get Pee",40,-1,1);
  TH1F *hcosSun_noOsci= new TH1F("hcosSun_noOsci","",40,-1,1);

  for(int i =0;i<T->GetEntries();i++)
  {
    //gCosPMTfit->clear();
    T->GetEntry(i);
    if(sqrt(posx**2+posy**2+(posz-108)**2)<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 && scaleLogL>10 && Utest<0.95 && 0<Gtest && Gtest<1.9 && energy<=15 && energy>=5 && zfactor>-11 && zfactor<1)
    {
       if( energymc<15.5 && energymc>0.1 )
       {
          double weight = 0;
          weight = scaleNue*scaleTime*scaleDataClean;
          hcosSun_noOsci->Fill(cosThetaToSun, weight);
          //cout<<energy<<", "<<index<<", "<<Etot[index-1]<<", "<<pee_bs05op_b8[index-1]<<endl;
        }
     }
   }

  cout<<"Erecon, no oscillation: "<<hcosSun_noOsci->Integral()<<endl;
//  cout<<"No oscillation "<<hcosSun_noOsci->Integral()<<endl;
  return hcosSun_noOsci;
}

void scalePeeMCsmear()
{
  double scaleNue = 1./1700;
  double scaleNumu = 1./9600;
  double scaleTime = 0.96;
  double scaleDataClean = (1-1.7/100);

  bool checkNumu = 0;
  TString path1 = "";//SmearNue/";
  TString path2 = "";//SmearNumu/";
  TFile *fMC = new TFile(path1+"GetSolarMC_smearDirDown_5to15MeV_Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to207718.root");
  TFile *fMC_numu = new TFile(path2+"GetSolarMC_smearDirDown_5to15MeV_Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NumuRun_r200004to207718.root");
  TFile *fnew = new TFile("saveMCsum_smearDirDown_nue_numu_pdfs.root","recreate");

  TTree *Tnue = (TTree*)fMC->Get("T2");
  TTree *Tnumu = (TTree*)fMC_numu->Get("T2");
  TH1F *hcosSunNue_peeEtrue = new TH1F("hcosSunNue_peeEtrue","",40,-1,1);
  TH1F *hcosSunNumu_peeEtrue = new TH1F("hcosSunNumu_peeEtrue","",40,-1,1);
  TH1F *hcosSun_noOsci = new TH1F("hcosSun_noOsci","",40,-1,1);

  hcosSunNue_peeEtrue->SetName("hcosSunNue_peeEtrue");
  hcosSunNumu_peeEtrue->SetName("hcosSunNumu_peeEtrue");
  hcosSun_noOsci->SetName("hcosSun_noOsci");

  hcosSunNue_peeEtrue = hPee(Tnue, 0);
  double countNue = hcosSunNue_peeEtrue->Integral();
  hcosSunNumu_peeEtrue = hPee(Tnumu, 1);
  double countNumu = hcosSunNumu_peeEtrue->Integral();
  hcosSun_noOsci = hPee(Tnue);

  cout<<"oscillated Nue+Numu "<<countNue+countNumu<<endl;

//  hcosSunNue_peeEtrue->SetName("hcosSunNue_peeEtrue");
//  hcosSunNumu_peeEtrue->SetName("hcosSunNumu_peeEtrue");
//  hcosSun_noOsci->SetName("hcosSun_noOsci");

  TH1F *hcosSun_peeEtrue_sum = (TH1F*)hcosSunNue_peeEtrue->Clone();;
  hcosSun_peeEtrue_sum->Add(hcosSunNumu_peeEtrue);
  hcosSun_peeEtrue_sum->SetName("hcosSunSumNueNumu_peeEtrue");

  fnew->cd();
  hcosSunNue_peeEtrue->Write();
  hcosSunNumu_peeEtrue->Write();
  hcosSun_peeEtrue_sum->Write();
  hcosSun_noOsci->Write();

}
