{
  TH1F *h = new TH1F("h","sin x",628,-TMath::Pi(),TMath::Pi());
  TH1F *hd = new TH1F("hd","derivative",628,-TMath::Pi(), TMath::Pi());
  double k = -1;
  for(int i = 0;i<628;i++)
  {
   h->SetBinContent(i+1,sin(k));
   k = k+2*TMath::Pi()/628;
  }

  for(int i = 0;i<628;i++)
  {
   if(i<628-1) hd->SetBinContent(i+1,h->GetBinContent(i+1)-h->GetBinContent(i));
   else hd->SetBinContent(i+1,0);
  }

  hd->Scale(h->GetMaximum()/hd->GetMaximum());
  h->Draw();
  hd->Draw("same");
}
