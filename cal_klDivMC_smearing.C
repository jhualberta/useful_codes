//#include <algorithm>
//#include "TGraph.h"
//#include "TCanvas.h"
//!!!!! Use ROOT6
// calculate KL-div
{
  //water angular pdf, 200, -1,1
//  double angularPDF[200]={0.0140342, 0.0141184, 0.0143289, 0.0145342, 0.0148263, 0.0151684, 0.0152658, 0.0150974, 0.0156105, 0.0155079, 0.0157974, 0.0162526, 0.0168789, 0.0167263, 0.0168289, 0.0175921, 0.0176895, 0.0181237, 0.0182763, 0.0186289, 0.0186, 0.0193868, 0.0194711, 0.0200395, 0.0198184, 0.0208079, 0.0209, 0.0214684, 0.0218921, 0.0224632, 0.0230737, 0.0233053, 0.0238079, 0.0241263, 0.0246211, 0.0248105, 0.0259684, 0.0259658, 0.0266605, 0.0274895, 0.0280605, 0.0284184, 0.0289421, 0.0294237, 0.0300816, 0.0305421, 0.0313974, 0.0316184, 0.0328868, 0.0332158, 0.0333158, 0.0349079, 0.0346316, 0.0361368, 0.0371474, 0.0369553, 0.0383737, 0.0394105, 0.0402605, 0.0404658, 0.0409132, 0.0419579, 0.0429605, 0.0440158, 0.0450132, 0.0461395, 0.0470763, 0.0479711, 0.0487868, 0.0500579, 0.0512711, 0.0515763, 0.0533842, 0.0542184, 0.0553368, 0.0565789, 0.0570026, 0.0590368, 0.0604763, 0.0614895, 0.0627368, 0.0642579, 0.0665026, 0.0664237, 0.0681132, 0.0696211, 0.0720237, 0.0737974, 0.0745342, 0.0763921, 0.0784211, 0.0802868, 0.0812921, 0.0829921, 0.0848526, 0.0870789, 0.0892105, 0.0913632, 0.0935263, 0.0953184, 0.0971263, 0.0994947, 0.101568, 0.10505, 0.106989, 0.108747, 0.112716, 0.115024, 0.117761, 0.120861, 0.122789, 0.126653, 0.129711, 0.133795, 0.137197, 0.140158, 0.143882, 0.147392, 0.150787, 0.153826, 0.1595, 0.162261, 0.166605, 0.171576, 0.174908, 0.180482, 0.186182, 0.191463, 0.196187, 0.201682, 0.206774, 0.213574, 0.218197, 0.224616, 0.231679, 0.239163, 0.245353, 0.253068, 0.260413, 0.268808, 0.275805, 0.284239, 0.293247, 0.302471, 0.310608, 0.321092, 0.332032, 0.344274, 0.355942, 0.367397, 0.380934, 0.392561, 0.407068, 0.421218, 0.438945, 0.454755, 0.473458, 0.490255, 0.508824, 0.529003, 0.552076, 0.572608, 0.599029, 0.625458, 0.654121, 0.684313, 0.717371, 0.74935, 0.792632, 0.829653, 0.877292, 0.92525, 0.970547, 0.996903, 1.00237, 0.987555, 0.947292, 0.907053, 0.869334, 0.841003, 0.809482, 0.776424, 0.750021, 0.724474, 0.700124, 0.674739, 0.650079, 0.629787, 0.611989, 0.589055, 0.569929, 0.550379, 0.532632, 0.515161, 0.496537, 0.482984, 0.466082, 0.450729, 0.437766, 0.425111}; 

  const int binTot = 40;	
  // 4 - 5 MeV
  double angularPDF1[binTot] = {0.00496146,0.00504869,0.00500962,0.0055209,0.00502723,0.00547247,0.00580516,0.00573003,0.00612657,0.0063335,0.0063024,0.00681038,0.00622783,0.00717169,0.00767004,0.00786569,0.00808914,0.00872315,0.00932469,0.0103924,0.0105462,0.011569,0.0124757,0.0136538,0.0155071,0.0176554,0.0188984,0.0232261,0.0257,0.0301125,0.0362866,0.0446785,0.055029,0.0694291,0.149666,0.104073,0.0750265,0.0618837,0.050365,0.0406056};
  //5 - 6 MeV
  double angularPDF2[binTot] = {0.00484151,0.0048969,0.00496314,0.00531255,0.00498975,0.00543391,0.00561689,0.0054016,0.00584033,0.00615363,0.00615309,0.00652395,0.00598612,0.00708376,0.00735227,0.00752982,0.00788276,0.00827995,0.00898909,0.0100512,0.0102309,0.0112028,0.0120624,0.0131464,0.0151362,0.0171504,0.0186862,0.0228615,0.0255778,0.0305216,0.0367958,0.0462551,0.057739,0.0732452,0.14346,0.106413,0.0781052,0.0632,0.0496911,0.0392371};
	//6 - 7 MeV
  double angularPDF3[binTot] = {0.00465902,0.00476542,0.00481659,0.00527483,0.00483133,0.00530461,0.00549831,0.00540059,0.00571052,0.0059421,0.00600021,0.00630031,0.00587271,0.00685309,0.00715232,0.00736163,0.0077343,0.00815177,0.00871901,0.00968724,0.0099726,0.0108538,0.0116512,0.0129872,0.0147259,0.0167583,0.0181547,0.0223379,0.0252943,0.0303578,0.0370741,0.0471326,0.0595739,0.0765364,0.139737,0.108814,0.0810711,0.0640314,0.0492024,0.0376975};
  //7 - 8 MeV
   double angularPDF4[binTot] = {0.00450628,0.00462101,0.00469356,0.00513191,0.00479615,0.00518961,0.00538432,0.00526723,0.00558274,0.00590399,0.00581153,0.00629982,0.00593132,0.00669193,0.0070422,0.00716908,0.00759393,0.00795196,0.0084578,0.00952582,0.00972559,0.0106785,0.0115876,0.0127454,0.0144198,0.0162876,0.0181301,0.0221167,0.0251233,0.0306126,0.0373194,0.0475242,0.0610815,0.0787698,0.137481,0.110518,0.0827311,0.0641354,0.0485868,0.0368736};
//8 - 9 MeV
   double angularPDF5[binTot] = {0.00437522,0.00448784,0.00466786,0.00502833,0.00466525,0.00515356,0.00529618,0.00518704,0.00548186,0.00574232,0.00577275,0.0060593,0.00574014,0.00657544,0.00677416,0.00699418,0.00731334,0.00781557,0.00847911,0.00917223,0.00957705,0.010401,0.0112246,0.012376,0.0141832,0.0161216,0.017757,0.0217839,0.0248803,0.030133,0.0374312,0.0482118,0.0622097,0.0814682,0.136735,0.112162,0.0842819,0.0646721,0.0477922,0.0358166};
   //9 - 10 MeV
   double angularPDF6[binTot] = {0.0044098,0.00453877,0.00449558,0.00487262,0.00463134,0.00504048,0.00519599,0.00502258,0.0054132,0.00564153,0.00569584,0.00593712,0.00578532,0.00632096,0.00668875,0.00697755,0.00716083,0.00761872,0.00828025,0.00902756,0.00946447,0.0101229,0.0108258,0.0121371,0.0137391,0.0156848,0.0174658,0.0212387,0.0245242,0.0301423,0.0376721,0.0488176,0.0636798,0.0837936,0.136244,0.113645,0.0856597,0.0645308,0.0471563,0.0347014};
//10 - 15 MeV
   double angularPDF7[binTot] = {0.00423936,0.00436983,0.00445795,0.00469844,0.0045558,0.00477634,0.00506308,0.0048922,0.00525246,0.00545499,0.00557231,0.00576169,0.00558205,0.0062495,0.00647685,0.00671686,0.00708102,0.0074403,0.00801525,0.00880052,0.00912767,0.00988859,0.0108559,0.0118588,0.0134814,0.0153635,0.0170933,0.0208107,0.0243252,0.0298722,0.0377171,0.0488807,0.0650095,0.0860534,0.135095,0.116815,0.0874833,0.0646726,0.0463817,0.0337581};
   //   //10 - 11 MeV
//   double angularPDF7[binTot] = {0.00429573,0.00434097,0.00455658,0.00469422,0.0046971,0.00474331,0.00504939,0.00490116,0.00527943,0.005573,0.00563749,0.00578187,0.0056885,0.00629393,0.00636612,0.00675305,0.00710918,0.0075375,0.00816411,0.00896781,0.0092402,0.00997846,0.0109275,0.0119978,0.0136033,0.0154215,0.0172792,0.0209425,0.0247763,0.0299902,0.0379368,0.0488711,0.0643725,0.0845411,0.1351,0.116068,0.0868213,0.0644293,0.0468902,0.0343832};
//   //11 - 12 MeV
//   double angularPDF8[binTot] = {0.00430733,0.00456722,0.00442688,0.00483058,0.004465,0.00483232,0.00516845,0.00506795,0.00523429,0.00535211,0.00552364,0.00588229,0.00548032,0.00622016,0.00664638,0.00680925,0.00699811,0.00744166,0.00791121,0.00865277,0.00912925,0.0101376,0.0107042,0.0119153,0.0135041,0.0154464,0.0169608,0.0206097,0.0239831,0.0296869,0.0372863,0.0484999,0.0652874,0.0869661,0.135047,0.116733,0.0876089,0.0647988,0.0462076,0.0336703};
//   //12 - 13 MeV
//   double angularPDF9[binTot] = {0.00394153,0.00417759,0.00417401,0.00449592,0.0043457,0.0047284,0.00479994,0.00460322,0.00523987,0.00545447,0.00547236,0.00559397,0.00549739,0.00632719,0.00639515,0.00660617,0.00714983,0.00717129,0.00790452,0.00856978,0.00874147,0.00935308,0.0110413,0.0115456,0.0132767,0.015226,0.0169285,0.0209273,0.0238745,0.0299156,0.0378558,0.0490116,0.0654072,0.0883196,0.13531,0.117817,0.0890134,0.0653678,0.0460465,0.0323728};
//   //13 - 14 MeV
//   double angularPDF10[binTot] = {0.00413432,0.00408209,0.00443895,0.00444766,0.00438673,0.00494377,0.00521359,0.0045434,0.00527452,0.004996,0.0056923,0.00543119,0.00551823,0.00594471,0.00698917,0.00657139,0.0070414,0.00741566,0.00759844,0.00855586,0.00896494,0.00919994,0.0104794,0.0108624,0.0129774,0.0147269,0.016546,0.0201058,0.0228388,0.0292623,0.0380792,0.0502994,0.0674547,0.0884831,0.135667,0.1202,0.0889618,0.0649218,0.0445462,0.0322042};
//   //14 - 15 MeV
//   double angularPDF11[binTot] = {0.00418285,0.00443566,0.00438969,0.00501023,0.00418285,0.00468847,0.00528602,0.00512514,0.00487233,0.00521707,0.00498724,0.00563076,0.00510216,0.00588357,0.00604445,0.00572269,0.0071706,0.00689481,0.00765324,0.0088943,0.00933097,0.0096987,0.0109627,0.0124336,0.0129163,0.0154444,0.0169153,0.0214429,0.0249132,0.0308428,0.0363356,0.0495737,0.0675232,0.0890809,0.132725,0.12036,0.0878858,0.063685,0.0435522,0.0330031};

  TH1F *hpdf1 = new TH1F("hpdf1","",binTot,-1,1);
  TH1F *hpdf2 = new TH1F("hpdf2","",binTot,-1,1);
  TH1F *hpdf3 = new TH1F("hpdf3","",binTot,-1,1);
  TH1F *hpdf4 = new TH1F("hpdf4","",binTot,-1,1);
  TH1F *hpdf5 = new TH1F("hpdf5","",binTot,-1,1);
  TH1F *hpdf6 = new TH1F("hpdf6","",binTot,-1,1);
  TH1F *hpdf7 = new TH1F("hpdf7","",binTot,-1,1);
//  TH1F *hpdf8 = new TH1F("hpdf8","",binTot,-1,1);
//  TH1F *hpdf9 = new TH1F("hpdf9","",binTot,-1,1);
//  TH1F *hpdf10 = new TH1F("hpdf10","",binTot,-1,1);
//  TH1F *hpdf11 = new TH1F("hpdf11","",binTot,-1,1);

  for(int i = 0;i<binTot;i++)
  {
   double q1 = angularPDF1[i];
   hpdf1->SetBinContent(i+1,q1);
   double q2 = angularPDF2[i];
   hpdf2->SetBinContent(i+1,q2);
   double q3 = angularPDF3[i];
   hpdf3->SetBinContent(i+1,q3);
   double q4 = angularPDF4[i];
   hpdf4->SetBinContent(i+1,q4);
   double q5 = angularPDF5[i];
   hpdf5->SetBinContent(i+1,q5);
   double q6 = angularPDF6[i];
   hpdf6->SetBinContent(i+1,q6);
   double q7 = angularPDF7[i];
   hpdf7->SetBinContent(i+1,q7);
   //double q8 = angularPDF8[i];
   //hpdf8->SetBinContent(i+1,q8);
   //double q9 = angularPDF9[i];
   //hpdf9->SetBinContent(i+1,q9);
   //double q10 = angularPDF10[i];
   //hpdf10->SetBinContent(i+1,q10);
   //double q11 = angularPDF11[i];
   //hpdf11->SetBinContent(i+1,q11);
  }
//// calculate KL-div
//  hpdf1->Draw();
  const char* filename = "Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to207718.root";
//Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NumuRun_r200004to207718.root";
//Merged_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to207718.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214_AvRun_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214_ExwaterRun_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterBi214Run_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_Nue_ExwaterRun_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterSolar_NueRun_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208_AvRun_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208_ExwaterRun_r200004to203602_s0_p0.root
// Merged_MP_ExtractTres_WaterMP6176_nhit15_GeneralPhysicsMC_WaterTl208Run_r200004to203602_s0_p0.root
  TFile *fMC = new TFile(filename);
  TTree *Tdata = (TTree*)fMC->Get("T");
  UInt_t nhits;
  double energy0, cosThetaToSun0, dirx0, diry0, dirz0;
  double energy = 0, klDiv=-1000, itr = 0, posx = 0, posy = 0, posz = 0, posRad = 0, dirx = 0, diry = 0, dirz = 0, beta14 =0, Gtest= 0, Utest = 0, scaleLogL = 0, cosThetaToSun = 0;
  double medianProbHit = 0, medianProb = 0, medianDevHit=0, medianDev=0;
  double udotR = 0, zfactor = 0;
  double posxmc = 0, posymc = 0, poszmc = 0, dirxmc = 0, dirymc = 0, dirzmc = 0, energymc = 0;

  UInt_t eventGTID = 0, runNumber = 0;
  //  std::vector<double> *gcosPMT0 = new std::vector<double>;
  std::vector<double> *gCosPMTfit = new std::vector<double>;
  std::vector<double> *gtRes = new std::vector<double>;
  Tdata->SetBranchAddress("energy", &energy0);////!!!! original energy!!!
  Tdata->SetBranchAddress("nhits", &nhits);
  Tdata->SetBranchAddress("gCosPMTfit", &gCosPMTfit);
  Tdata->SetBranchAddress("gtRes",&gtRes);
  Tdata->SetBranchAddress("itr", &itr);
  Tdata->SetBranchAddress("posx", &posx);
  Tdata->SetBranchAddress("posy", &posy);
  Tdata->SetBranchAddress("posz", &posz);
  Tdata->SetBranchAddress("dirx", &dirx);
  Tdata->SetBranchAddress("diry", &diry);
  Tdata->SetBranchAddress("dirz", &dirz);
  Tdata->SetBranchAddress("beta14", &beta14);
  Tdata->SetBranchAddress("cosThetaToSun", &cosThetaToSun);
  Tdata->SetBranchAddress("Gtest", &Gtest);
  Tdata->SetBranchAddress("Utest", &Utest);
  Tdata->SetBranchAddress("scaleLogL", &scaleLogL);
  Tdata->SetBranchAddress("runNumber",&runNumber);
  Tdata->SetBranchAddress("eventGTID",&eventGTID);
  Tdata->SetBranchAddress("medianProbHit", &medianProbHit);
  Tdata->SetBranchAddress("medianProb", &medianProb);
  Tdata->SetBranchAddress("medianDevHit", &medianDevHit);
  Tdata->SetBranchAddress("medianDev", &medianDev);
  Tdata->SetBranchAddress("posxmc", &posxmc);
  Tdata->SetBranchAddress("posymc", &posymc);
  Tdata->SetBranchAddress("posymc", &poszmc);
  Tdata->SetBranchAddress("dirxmc", &dirxmc);
  Tdata->SetBranchAddress("dirymc", &dirymc);
  Tdata->SetBranchAddress("dirzmc", &dirzmc);
  Tdata->SetBranchAddress("energymc", &energymc);

  TString ffname0(filename);
  ffname0 = "GetSolarMC_smearShiftYup_5to15MeV_"+ffname0;
  // GetSolarMC_smearBeta14up_5to15MeV_
  // GetSolarMC_smearBeta14down_5to15MeV_
  // GetSolarMC_smearDirScaleUp_5to15MeV_

  TFile *fnew = new TFile(ffname0,"recreate");

  TTree *tree= new TTree("T2","solar");
  std::vector<double> *gCosPMT = new std::vector<double>;
  //std::vector<double> *gcosPMTMC = new std::vector<double>;
  tree->Branch("nhits",&nhits);
  tree->Branch("energy",&energy,"energy/D");
  tree->Branch("klDiv", &klDiv);
  // tree->Branch("gtRes". &gtRes);
  tree->Branch("gCosPMTfit",&gCosPMT);
  tree->Branch("itr", &itr);
  tree->Branch("posx", &posx);
  tree->Branch("posy", &posy);
  tree->Branch("posz", &posz);
  tree->Branch("posRad", &posRad);
  tree->Branch("dirx", &dirx);
  tree->Branch("diry", &diry);
  tree->Branch("dirz", &dirz);
  tree->Branch("beta14", &beta14);
  tree->Branch("cosThetaToSun", &cosThetaToSun);
  tree->Branch("Gtest", &Gtest);
  tree->Branch("Utest", &Utest);
  tree->Branch("scaleLogL", &scaleLogL);
  tree->Branch("runNumber",&runNumber);
  tree->Branch("eventGTID",&eventGTID);  
  tree->Branch("medianProbHit", &medianProbHit);
  tree->Branch("medianProb", &medianProb);
  tree->Branch("medianDevHit", &medianDevHit);
  tree->Branch("medianDev", &medianDev);
  tree->Branch("udotR", &udotR);
  tree->Branch("zfactor", &zfactor);

  tree->Branch("posxmc", &posxmc);
  tree->Branch("posymc", &posymc);
  tree->Branch("posymc", &poszmc);
  tree->Branch("dirxmc", &dirxmc);
  tree->Branch("dirymc", &dirymc);
  tree->Branch("dirzmc", &dirzmc);
  tree->Branch("energymc", &energymc);

  //  tree->Branch("gcosPMT",&gCosPMTfit);
  TRandom3 rtheta;
  double timeCut1 = -5, timeCut2 = 1;//time res cuts
  for(int i =0;i<Tdata->GetEntries();i++)
  {
    gCosPMT->clear();
    Tdata->GetEntry(i);
    // hpdf->Scale(scaleOrigin/hpdf->Integral());
    energy = energy0;//eCorr->CorrectEnergyRSP(energy0,2);
    vector<double>& vecRef = *gCosPMTfit; // vector is not copied here
    vector<double>& vecTres = *gtRes;
    //cout<<i<<" "<<energy0<<endl;
    //apply cuts!!!

    /**!!!!! smearing beta14: +0.001/-0.036 **/
    // beta14 = beta14+0.01;

    /*!!!!! smearing direction: +0.013/-0.101  original, no changing signs! --------------------------
    // change sign in the transforamtion eq.:
    //    -costhetaToSun' = 1+(-costhetaToSun-1)*(1+delta) => costhetaToSun' = -1+(costhetaToSun+1)*(1+delta)
    // original transformation, not used: costheta' = 1+(costheta-1)*(1+delta),  ///
    // if doing so, must be -0.013/+0.101!!!
    */

    double deltaTheta = 0.101;
    /// double deltaTheta = -0.013;

    double deltaTheta = 0.013; // up, not change sign
    //double deltaTheta = -0.101; // down, not change sign

    // /* ---------- Print out ---------*/
    // // cout<<"before = "<<cosThetaToSun<<endl;
    /// cosThetaToSun = -1+(cosThetaToSun+1)*(1+deltaTheta); // change sign
    cosThetaToSun = 1+(cosThetaToSun-1)*(1+deltaTheta); // don't change sign !!!!
    /// // // cout<<"after = "<<cosThetaToSun<<endl;

    if(cosThetaToSun>1) {
       cosThetaToSun = rtheta.Uniform(-1,1);
       cout<<">+1, randomize "<<cosThetaToSun<<endl;
    }
    if(cosThetaToSun<-1) {
       cosThetaToSun = rtheta.Uniform(-1,1);
       cout<<"<-1, randomize "<<cosThetaToSun<<endl;
    }

    ///!!! smearing position scale
    // up: (0.07/100,0.02/100,0.08/100) 
    // TVector3 scaleP(0.07/100,0.02/100,0.08/100); 
    /// TVector3 scaleP(-0.06/100,-0.07/100,-0.01/100);
    // down: (-0.06/100,-0.07/100,-0.01/100)
    // posx = posx*(1+scaleP.X()); 
    // posy = posy*(1+scaleP.Y()); 
    // posz = posz*(1+scaleP.Z()); 

    /**!!! smearing position shifts, one by one !!!**/
//    posx = posx + 6.48;
//    posy = posy + 6.13;
//    posz = posz + 6.71;
//
//    posx = posx - 5.98;
//    posy = posy - 4.11;
//    posz = posz - 4.82;

    /// !!!!! smearing energy scale
    //double eScale = +2.0/100;
    // double eScale = -2.0/100;
    //energy = (1+eScale)*energy;

    /// !!!!! smearing energy resolution
    // double eResol = 0.011;
    // double sigma = sqrt(energy*((1+eResol)*(1+eResol)-1));
    // energy = energy + TMath::Gaus(0,sigma);
    // energy = energy - TMath::Gaus(0,sigma);

    if(nhits>20 && energy>=5 && energy<=15 && sqrt(posx**2+posy**2+(posz-108)**2)<5500 && itr>0.55 && beta14>-0.12 && beta14<0.95 )
    {
       TH1F *hChtemp = new TH1F("hChtemp","",binTot,-1,1);
       // cout<<gCosPMTfit->size()<<endl;
       for(int j = 0; j<gCosPMTfit->size();j++)
       { 
         double valCh = vecRef[j];
	 double tres = vecTres[j];
         // cout<<valCh<<endl;
         if(tres>timeCut1 && tres<timeCut2)
         {		 
	   hChtemp->Fill(valCh);
	   gCosPMT->push_back(valCh);
	 }
	 //gCosPMTfit.push_back(valCh);
       }
       posRad = sqrt(posx*posx+posy*posy+(posz-108)*(posz-108));
       udotR = (posx*dirx+posy*diry+posz*dirz)/sqrt(posx*posx+posy*posy+posz*posz);
       zfactor = 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb);
       // calculate scale
       // double lowVal = 0;
       // int countlowVal = 0;
       // for(int j1 = 0;j1<200;j1++)
       // {
       //   double val = -1 + 0.01*i;
       //   double y = hChtemp->GetBinContent(i);
       //   if(val<0.4 && y!=0 )
       //   {
       //     lowVal += y;
       //     countlowVal++;
       //   }
       // }
       double scale = hChtemp->Integral();
       // cout<<scale<<" "<<energy<<endl;
       // if(countlowVal!=0) scale = lowVal/countlowVal;
       // cout<<scale<<" "<<angularPDF[binTot-1]<<" "<<scale/angularPDF[binTot-1]<<endl; 
       // hpdf->Scale(scale/angularPDF[binTot-1]);
       hpdf1->Scale(scale/hpdf1->Integral());hpdf2->Scale(scale/hpdf2->Integral());hpdf3->Scale(scale/hpdf3->Integral());hpdf4->Scale(scale/hpdf4->Integral());hpdf5->Scale(scale/hpdf5->Integral());hpdf6->Scale(scale/hpdf6->Integral());hpdf7->Scale(scale/hpdf7->Integral());
       // hpdf8->Scale(scale/hpdf8->Integral());hpdf9->Scale(scale/hpdf9->Integral());hpdf10->Scale(scale/hpdf10->Integral());hpdf11->Scale(scale/hpdf11->Integral());
       // hChtemp->Draw();
       // hpdf10->SetLineColor(kRed);hpdf10->Draw();hChtemp->Draw("same");

       int binCut = 26;
       double kl1 = 0;
       for(int ii = binCut;ii<binTot;ii++) //binCut = 26
       { 
         double p1 = hChtemp->GetBinContent(ii+1);
         double q = 0; //hpdf->GetBinContent(ii+1);
         if(energy>=4 && energy<5) q = hpdf1->GetBinContent(ii+1);
         if(energy>=5 && energy<6) q = hpdf2->GetBinContent(ii+1);
         if(energy>=6 && energy<7) q = hpdf3->GetBinContent(ii+1);
         if(energy>=7 && energy<8) q = hpdf4->GetBinContent(ii+1);
         if(energy>=8 && energy<9) q = hpdf5->GetBinContent(ii+1);
         if(energy>=9 && energy<10) q = hpdf6->GetBinContent(ii+1);
         if(energy>=10 && energy<15) q = hpdf7->GetBinContent(ii+1);
         //if(energy>=11 && energy<12) q = hpdf8->GetBinContent(ii+1);
         //if(energy>=12 && energy<13) q = hpdf9->GetBinContent(ii+1);
         //if(energy>=13 && energy<14) q = hpdf10->GetBinContent(ii+1);
         //if(energy>=14 && energy<15) q = hpdf11->GetBinContent(ii+1);
         if(p1!=0 && q!=0)
         {
            //if() cout<<p1<<" q "<<q<<endl;
            kl1 = kl1+p1*log(p1/q);
         }
       }
       klDiv = kl1;
       if(kl1 == 0) klDiv = -1000;
       //if(klDiv == -1000)cout<<i<<" "<<runNumber<<" "<<eventGTID<<", klDiv = "<<klDiv<<endl;
       // if(klDiv<-1000) cout<<i<<", "<<energy<<" "<<klDiv<<endl;
       delete hChtemp;
       tree->Fill();
    } // cut conditions
  }
//  TFile *fnew = new TFile("results.root","recreate");
  fnew->cd();
  tree->Write();
  fnew->Close();

}
