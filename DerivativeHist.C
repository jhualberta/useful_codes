{
   TF1 f1("f1","x*x",-10,10);
   f1->Draw(); 
   double histBin = 200;
 
   TH1F *h1 = new TH1F("h1","",200,-10,10);
   double xx= -10;
   for(xx = -10;xx<10;xx+=0.1)
   {
    h1->Fill(xx,xx*xx);
   }
   h1->Draw();

   TH1F *hDeriv = new TH1F("hDeriv","",200-1,-10,10);
   for (j = 0; j < fEntriesTimeCh - 1; j++ ) fDerivativeCh.push_back( ( fPDFCh[j + 1] - fPDFCh[j] ) * 1/fChBinWidth ) ;
for(j=0; j< fEntriesTimeCh -1 ; j++ ){

   dh1->SetBinContent(j+1, fDerivativeCh[j]);
}
