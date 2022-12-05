{
  TFile *ff = new TFile("Combined1.root");
  TTree *tt = (TTree*)ff->Get("data_satCorr");  
  TH1F *h = new TH1F("h","",300,28000,28300);
  tt->Project("h","runID");
  // h->Draw();
  for(int i = 0;i<h->GetEntries();i++)
  {
    int run = h->GetXaxis()->GetBinCenter(i+1);
    int counts = h->GetBinContent(i+1);
    if(counts>0)
	    cout<<run<<endl;
  
  }
  h->Draw();
}
