{

  TFile *ff = new TFile("testSolar_nueTl208_lowE0p8to5_Likelihood.root");
  TDirectory *dir = ff->GetDirectory("/dataset/Method_Likelihood/Likelihood");
  TH1D* hRoc; //NOT TH1F!! must be exactly
  dir->GetObject("MVA_Likelihood_rejBvsS",hRoc);
  hRoc->Draw();
  double auc = hRoc->Integral()/100;

  cout<<"AUC: "<<auc<<endl;





}
