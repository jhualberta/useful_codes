{
  TFile *ff = new TFile("ExtractTres_FitterMP6176_Water_2MeV_1E4evt_iso_center.root","read");
  TTree* tree = (TTree*)ff->Get("T");
   
  TH1F* hDeltaX = new TH1F("hDeltaX", "multipath X - mc X, FECD cut", 2000,-9000,9000);
  TH1F* hDeltaY = new TH1F("hDeltaY", "multipath Y - mc Y, FECD cut", 2000,-9000,9000);
  TH1F* hDeltaZ = new TH1F("hDeltaZ", "multipath Z - mc Z, FECD cut", 2000,-9000,9000);
  TH1F* hRadialBias = new TH1F("hRadialBias", "(Xfit-Xmc)*Xmc.Unit()", 2000,-9000,9000);

  TH1F* hRadialBias = new TH1F("hRadialBias", "(Xfit-Xmc)*Xmc.Unit()", 2000,-9000,9000);
  TH1F* hcostheta = new TH1F("hcostheta", "ufit*umc", 200,-1,1);

  double posx, posy, posz, posRad, energyOrigin, energy, time, dirx, diry, dirz, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, iso, thij;
  double mcPosx, mcPosy, mcPosz, mcDirx, mcDiry, mcDirz;

  tree->SetBranchAddress("posx",&posx);
  tree->SetBranchAddress("posy",&posy);
  tree->SetBranchAddress("posz",&posz);
  tree->SetBranchAddress("dirx",&dirx);
  tree->SetBranchAddress("diry",&diry);
  tree->SetBranchAddress("dirz",&dirz);

  tree->SetBranchAddress("mcPosx",&mcPosx);
  tree->SetBranchAddress("mcPosy",&mcPosy);
  tree->SetBranchAddress("mcPosz",&mcPosz);

  tree->SetBranchAddress("mcDirx",&mcDirx);
  tree->SetBranchAddress("mcDiry",&mcDiry);
  tree->SetBranchAddress("mcDirz",&mcDirz);

  for(int i = 0; i<T->GetEntries();i++)
  {
    tree->GetEntry(i);
    hDeltaX->Fill(posx-mcPosx);
    hDeltaY->Fill(posy-mcPosy);
    hDeltaZ->Fill(posz-mcPosz);
    double mcposMag = sqrt(mcPosx*mcPosx+mcPosy*mcPosy+mcPosz*mcPosz);
    if(mcposMag != 0)
    {
      double rbias = ((posx-mcPosx)*mcPosx+(posy-mcPosy)*mcPosx+(posz-mcPosz)*mcPosz)/mcposMag;
      hRadialBias->Fill(rbias);   
    }
    TVector3 ufit(dirx, diry, dirz);
    TVector3 umc(mcDirx, mcDiry, mcDirz);
    hcostheta->Fill(ufit*umc);
  }

  TF1 *gx = new TF1("gx","gaus",hDeltaX->GetMean()-2000,hDeltaX->GetMean()+2000);
  TF1 *gy = new TF1("gy","gaus",hDeltaY->GetMean()-2000,hDeltaY->GetMean()+2000);
  TF1 *gz = new TF1("gz","gaus",hDeltaZ->GetMean()-2000,hDeltaZ->GetMean()+2000);
  TF1 *gR = new TF1("gR","gaus",hRadialBias->GetMean()-2000,hRadialBias->GetMean()+2000);

  gx->SetLineColor(kRed);gy->SetLineColor(kRed);gz->SetLineColor(kRed);
  TCanvas *c = new TCanvas("c","",800,600);
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TCanvas *c3 = new TCanvas("c3","",800,600);

  c->cd(); hDeltaX->GetXaxis()->SetTitle("mm");hDeltaX->GetYaxis()->SetTitle("counts");hDeltaX->Fit(gx);
  c1->cd();hDeltaY->GetXaxis()->SetTitle("mm");hDeltaY->GetYaxis()->SetTitle("counts");hDeltaY->Fit(gy);
  c2->cd();hDeltaZ->GetXaxis()->SetTitle("mm");hDeltaZ->GetYaxis()->SetTitle("counts");hDeltaZ->Fit(gz);
  c3->cd(); hcostheta->Draw();

  hRadialBias->Fit(gR);

  double mux = gx->GetParameter(1), muy = gy->GetParameter(1), muz = gz->GetParameter(1);
  double sigmax = gx->GetParameter(2), sigmay = gy->GetParameter(2), sigmaz = gz->GetParameter(2);
//  cout<<"x,y,z bias "<<mux<<" "<<muy<<" "<<muz<<endl;
//  cout<<"x,y,z resol "<<sigmax<<" "<<sigmay<<" "<<sigmaz<<endl;
  cout<<sigmax<<" "<<sigmay<<" "<<sigmaz<<" "<<mux<<" "<<muy<<" "<<muz<<endl;
  cout<<"total bias: "<<sqrt(mux*mux+muy*muy+muz*muz)<<endl;
  cout<<"radial bias "<<muR<<" "<<sigmaR<<endl;
  cout<<"total resolution: "<<sqrt(sigmax*sigmax+sigmay*sigmay+sigmaz*sigmaz)<<endl;
}
