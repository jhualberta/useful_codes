#include "TH1F.h"
#include "TMath.h"

TH1F *hist(bool check)
{
   TF1 *mygaus = new TF1("mygaus","TMath::Gaus(x,3,.5)",0,6);
   TF1 *mypoisson = new TF1("mypoisson","TMath::Poisson(x,6)",0,6);

   TH1F *h = new TH1F("h","",100,0,6);
   if(check) { h->FillRandom("mygaus",100); h->SetTitle("hgaus"); }
   else { h->FillRandom("mypoisson",100); h->SetTitle("hpoisson"); }
   return h;
}

void functionTH1F()
{
   TH1F *h;
   TH1F *h2;
   h = hist(1);
   h2 = hist(0);
   h2->SetLineColor(kRed);
   h->Draw();
   h2->Draw("same");


}
