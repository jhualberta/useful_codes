void mergeHistoBins() {
   const Int_t nbins = 100;
   TH1F *h = new TH1F("h","merging high bins",nbins,0,3);
   h->FillRandom("gaus",5000);
   Double_t xbins[nbins];
   Int_t i;
   TAxis *axis = h->GetXaxis();
   for (i=0;i<50;i++)  xbins[i] = axis->GetBinLowEdge(1+i);
   for (i=50;i<=60;i++) xbins[i] = axis->GetBinLowEdge(51+5*(i-50));
   for (i=0;i<=60;i++) printf("xbins[%d]=%g\n",i,xbins[i]);
   TH1F *hnew = new TH1F("hnew","merged",60,xbins);
   for (i=0;i<nbins;i++) {
      Double_t y = h->GetBinContent(i+1);
      Double_t x = axis->GetBinCenter(i+1);
      if (i < 50) hnew->Fill(x,y);
      else        hnew->Fill(x,y/5);
   }
   TCanvas *c1 = new TCanvas("c1","c1",10,10,700,900);
   c1->Divide(1,2);
   c1->cd(1);
   h->Draw();
   c1->cd(2);
   hnew->Draw();
}
