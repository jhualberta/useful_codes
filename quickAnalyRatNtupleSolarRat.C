{
  TFile *f1 = new TFile("Merged_fullNtuple_Processed_Analysis_r200004to207718.ntuple.root");

  TH1F *hcosE4to15 = new TH1F("hcosE4to15","",40,-1,1);
  TH1F *hcosE5to15 = new TH1F("hcosE5to15","",40,-1,1);
  TH1F *hcosE6to15 = new TH1F("hcosE6to15","",40,-1,1);

  TH2F *hRhoVsZ1 = new TH2F("hRhoVsZ1","",1000,0,9000,2000,-9000,9000);

  TTree* T1 = (TTree*)f1->Get("T");
 /// z-factor: 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb)<1 && 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb)>-11 
  T1->Project("hcosE4to15","cosThetaSun","!(triggerWord & 0x3F) || (triggerWord & 0xBEF9400) && ((dcApplied & 0xFB0000017FFE) & dcFlagged ) != (dcApplied & 0xFB0000017FFE) && zfactor<1 && zfactor>-11 && scaledLogL>10 && energyFOMGtest>0 && energyFOMGtest<1.9 && energyFOMUtest<0.95 && itr>0.55 && nhits>20 && beta14<0.95 && beta14>-0.12 && sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))<5500 && energy>4 && energy<15","colz");
  T1->Project("hcosE5to15","cosThetaSun","!(triggerWord & 0x3F) || (triggerWord & 0xBEF9400) && ((dcApplied & 0xFB0000017FFE) & dcFlagged ) != (dcApplied & 0xFB0000017FFE) && zfactor<1 && zfactor>-11 && scaledLogL>10 && energyFOMGtest>0 && energyFOMGtest<1.9 && energyFOMUtest<0.95 && itr>0.55 && nhits>20 && beta14<0.95 && beta14>-0.12 && sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))<5500 && energy>5 && energy<15","colz");
  T1->Project("hcosE6to15","cosThetaSun","!(triggerWord & 0x3F) || (triggerWord & 0xBEF9400) && ((dcApplied & 0xFB0000017FFE) & dcFlagged ) != (dcApplied & 0xFB0000017FFE) && zfactor<1 && zfactor>-11 && scaledLogL>10 && energyFOMGtest>0 && energyFOMGtest<1.9 && energyFOMUtest<0.95 && itr>0.55 && nhits>20 && beta14<0.95 && beta14>-0.12 && sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))<5500 && energy>6 && energy<15","colz");

//  T1->Project("hRhoVsZ1","(posz-108):sqrt(posx*posx+posy*posy)","1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb)<1 && 1-3*(medianDevHit+medianDev)/(medianProbHit-medianProb)>-11 && scaledLogL>10 && energyFOMGtest>0 && energyFOMGtest<1.9 && energyFOMUtest<0.95 && itr>0.55 && nhits>15 && beta14<0.95 && beta14>-0.12 && (posx*posx+posy*posy+(posz-108)*(posz-108)<5500 && energy>4 && energy<15","colz");

  TFile *ff = new TFile("savedNtupule_rat.root","recreate");
  ff->cd();
  hcosE4to15->Write();
  hcosE5to15->Write();
  hcosE6to15->Write();

}
