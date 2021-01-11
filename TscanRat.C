{
  TFile *fname = new TFile("Merged_Analysis_r100002to100076.root");
  TTree *T = (TTree*)fname->Get("output");
  ((TTreePlayer*)(T->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(T->GetPlayer()))->SetScanFileName("toto.txt");
  T->Scan("runID:eventID:triggerWord:energy:nhits:itr:beta14:posx:posy:posz-108:dirx:diry:dirz:sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))","nhits>30 && itr>=0.55 && beta14>=-0.12 && beta14<=0.95 && sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))<=5300 && !( !(triggerWord & 0x3F) || (triggerWord & 0xBFF9400) ) && ( ((dcApplied & 0x7FFE) & dcFlagged) == (dcApplied & 0x7FFE) ) && energy>5");
  fname->Close();
}
