//----file scott.C
Double_t g2(Double_t *x, Double_t *par)
{
   Double_t e0 = (x[0]-par[1])/par[2];
   Double_t g0 = par[0]*TMath::Exp(-0.5*e0*e0);
   Double_t e1 = (x[0]-par[4])/par[5];
   Double_t g1 = par[3]*TMath::Exp(-0.5*e1*e1);
   return g0+g1;
}
void scott()
{
   // Create the original function
   TF1 *g2 = new TF1("g2",g2,0,1,6);
   g2->SetParameters(2,0.3,0.1, 1,0.7,0.15);

   //Get random numbers from g2 function
   TH1F *h1 = new TH1F("h1","Scott test",100,0,1);
   h1->FillRandom("g2",10000);

   //Draw h1 full
   h1->SetFillColor(17);
   h1->Draw();

   //Copy h1 in a clone h1c. Set range and color for h1c
   TH1F *h1c = (TH1F*)h1->Clone();
   h1c->SetFillColor(42);
   h1c->GetXaxis()->SetRange(50,90);
   h1c->Draw("same");

   //to draw a shaded area above and below an histogram range, we create
   //a TGraph object (here we shade bins 60 to 80).
   Int_t i;
   Int_t n = 2*(80-60);
   TGraph *gr = new TGraph(2*n);
   for (i=0;i<20;i++) {
      Float_t xlow = h1->GetBinLowEdge(60+i);
      Float_t xup  = h1->GetBinLowEdge(60+i+1);
      Float_t y    = h1->GetBinContent(60+i);
      Float_t yup  = 1.1*y;
      Float_t ydown= 0.9*y;
      gr->SetPoint(2*i,  xlow,yup);
      gr->SetPoint(2*i+1,xup, yup);
      gr->SetPoint(2*n-2*i-1,xlow, ydown);
      gr->SetPoint(2*n-2*i-2,xup, ydown);
   }
   gr->SetFillColor(2);
   gr->Draw("lf");
}
