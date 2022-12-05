{

  TFile *ff = new TFile("DanielRecalcTimeWin_L6Cuts_deap_cal_030681_0000.root");
  TH1F* hNpmt = new TH1F("hNpmt","", 100, 0, 100);
  TH1F* hNsubpulse = new TH1F("hNsubpulse", "", 100, 0, 100);
  TTree *tree = (TTree*)ff->Get("T");

  float qpe, fprompt, nSCBayes, fmaxpe, rprompt60Bayes;
  int nPMTs;
  int nSubpeaks[255];

  Long64_t eventID;
  tree->SetBranchAddress("nPMTs", &nPMTs);
  tree->SetBranchAddress("nSubpeaks", &nSubpeaks);
  tree->SetBranchAddress("eventID", &eventID);
  tree->SetBranchAddress("qpe",&qpe);
  tree->SetBranchAddress("fprompt",&fprompt);
  tree->SetBranchAddress("nSCBayes",&nSCBayes);
  tree->SetBranchAddress("rprompt60Bayes",&rprompt60Bayes);

  tree->Project("hNpmt","nPMTs","");
 
  // hNpmt->Draw();
  int sum = 0;
  for(int i = 0; i<tree.GetEntries(); i++) // loop events
  {
    tree->GetEntry(i);
    sumPMTpulse = 0;
    for(int j = 0; j<nPMTs; j++) // loop pmts
    {	
       npulse = nSubpeaks[j];
       sumPMTpulse += npulse;
       cout<<"event "<<eventID<<" nPMTs "<<nPMTs<<", has subpeaks "<<npulse<<endl; 
    }
    sum += sumPMTpulse;
    cout<<"total subpeaks in the event "<<sumPMTpulse<<endl; 
  }
  cout<<"summed "<<sum<<endl;
}
