#include "TH1F.h"
void printHistBin()
{
  TH1F *h1 = new TH1F("h1","",200,-10,10);
  TF1 *mygaus = new TF1("mygaus","TMath::Gaus(x,0,5)",-10,10); 
  h1->FillRandom("mygaus",1000); 

  double lowBin = h1->GetBinLowEdge(1);
  double highBin = h1->GetBinLowEdge(h1->GetNbinsX()+1);
  double cBin = h1->GetBinLowEdge(0);
  cout<<lowBin<<" "<<highBin<<endl;

  double increment = lowBin - h1->GetBinLowEdge(0);
  int NIncrement = (int)(highBin-lowBin)/increment;

  cout<<increment<<" "<<NIncrement<<endl;

  h1->Draw();



}
