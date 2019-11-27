{

  TFile *ff = new TFile("Resol_FitMPW_Tl208_BBpointOnAV6000mm.root");
  TFile *f_rat = new TFile("Resol_FitRat_Tl208_BBpointOnAV6000mm.root");
  //phi
  TH1D *hPhi_rat = (TH1D*)f_rat->Get("hPhi");
  TH1D *hPhi_ratCuts = (TH1D*)f_rat->Get("hPhi_bbCuts");
  TH2D *hPhiTheta_ratCuts = (TH2D*)f_rat->Get("hPhiTheta_bbCuts");
  TH1D *hPhiCor = (TH1D*)ff->Get("hPhiCor");
  TH1D *hPhiCor_bbCuts = (TH1D*)ff->Get("hPhiCor_bbCuts");
  TH2D *hPhiThetaCor_bbCuts = (TH2D*)ff->Get("hPhiThetaCor_bbCuts");
  //cosTheta
  TH1D *hTheta_rat = (TH1D*)f_rat->Get("hTheta");
  TH1D *hTheta_ratCuts = (TH1D*)f_rat->Get("hTheta_bbCuts");
  TH2D *hPhiTheta_rat= (TH2D*)f_rat->Get("hPhiTheta");
  TH2D *hPhiTheta_ratCuts = (TH2D*)f_rat->Get("hPhiTheta_bbCuts");
  TH1D *hThetaCor = (TH1D*)ff->Get("hThetaCor");
  TH1D *hThetaCor_bbCuts = (TH1D*)ff->Get("hThetaCor_bbCuts");
  TH2D *hPhiThetaCor = (TH2D*)ff->Get("hPhiThetaCor");
  TH2D *hPhiThetaCor_bbCuts = (TH2D*)ff->Get("hPhiThetaCor_bbCuts");

  hPhi_rat->SetLineColor(kBlack);
  hPhi_ratCuts->SetLineColor(kBlue);
  hPhiCor->SetLineColor(kBlack);
  hPhiCor_bbCuts->SetLineColor(kRed);
  hPhi_rat->SetTitle("Rat");
  hPhiCor->SetTitle("MPW");

  hPhi_rat->GetXaxis()->SetTitle("#phi");
  hPhi_ratCuts->GetXaxis()->SetTitle("#phi");
  hPhiCor->GetXaxis()->SetTitle("#phi");
  hPhiCor_bbCuts->GetXaxis()->SetTitle("#phi");

  hTheta_rat->SetLineColor(kBlack);
  hTheta_ratCuts->SetLineColor(kBlue);
  hThetaCor->SetLineColor(kBlack);
  hThetaCor_bbCuts->SetLineColor(kRed);

  hTheta_rat->SetTitle("Rat");
  hThetaCor->SetTitle("MPW");

  hTheta_rat->GetXaxis()->SetTitle("cos#theta");
  hTheta_ratCuts->GetXaxis()->SetTitle("cos#theta");
  hThetaCor->GetXaxis()->SetTitle("cos#theta");
  hThetaCor_bbCuts->GetXaxis()->SetTitle("cos#theta");
 
//  TExec *ex1 = new TExec("ex1","gStyle->SetStatTextColor(kBlue);"); 
//  TExec *ex2 = new TExec("ex2","gStyle->SetStatTextColor(kBlack);"); 
//  TExec *ex3 = new TExec("ex3","gStyle->SetStatTextColor(kRed);"); 

//Theta
  TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9);
  legend1->AddEntry(hTheta_rat,"Rat WaterFitter fitted #theta","l");
  legend1->AddEntry(hTheta_ratCuts,"Rat WaterFitter fitted #theta after cuts","l");
 
  TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
  legend2->AddEntry(hThetaCor,"MPW fitted cos#theta (after drive correction)","l");
  legend2->AddEntry(hThetaCor_bbCuts,"MPW fitted #theta after cuts","l");

  TLegend *legend3 = new TLegend(0.1,0.7,0.48,0.9);
  legend3->AddEntry(hTheta_ratCuts,"Rat WaterFitter fitted #theta after cuts","l");
  legend3->AddEntry(hThetaCor_bbCuts,"MPW fitted cos#theta after cuts","l");

//Phi
  TLegend *legend11 = new TLegend(0.1,0.7,0.48,0.9);
  legend11->AddEntry(hPhi_rat,"Rat WaterFitter fitted #phi","l");
  legend11->AddEntry(hPhi_ratCuts,"Rat WaterFitter fitted #phi after cuts","l");
 
  TLegend *legend22 = new TLegend(0.1,0.7,0.48,0.9);
  legend22->AddEntry(hPhiCor,"MPW fitted #phi (after drive correction)","l");
  legend22->AddEntry(hPhiCor_bbCuts,"MPW fitted #phi after cuts","l");

  TLegend *legend33 = new TLegend(0.1,0.7,0.48,0.9);
  legend33->AddEntry(hPhi_ratCuts,"Rat WaterFitter fitted #phi after cuts","l");
  legend33->AddEntry(hPhiCor_bbCuts,"MPW fitted #phi after cuts","l");

  hPhiTheta_rat->RebinX(4);hPhiTheta_rat->RebinY(4);
  hPhiThetaCor->RebinX(4);hPhiThetaCor->RebinY(4);
  hPhiTheta_ratCuts->RebinX(4);hPhiTheta_ratCuts->RebinY(4);
  hPhiThetaCor_bbCuts->RebinX(4);hPhiThetaCor_bbCuts->RebinY(4);

  hPhiTheta_ratCuts->GetXaxis()->SetTitle("cos#theta");hPhiTheta_ratCuts->GetYaxis()->SetTitle("cos#theta");
  hPhiThetaCor_bbCuts->GetXaxis()->SetTitle("cos#theta");hPhiThetaCor_bbCuts->GetYaxis()->SetTitle("cos#theta");

  hPhiTheta_ratCuts->GetXaxis()->SetRangeUser(0,1);hPhiTheta_ratCuts->GetYaxis()->SetRangeUser(0.5,1);
  hPhiThetaCor_bbCuts->GetXaxis()->SetRangeUser(0,1);hPhiThetaCor_bbCuts->GetYaxis()->SetRangeUser(0.5,1);

  hPhiTheta_ratCuts->SetTitle("Rat WaterFitter");
  hPhiThetaCor_bbCuts->SetTitle("MPW");
  TCanvas *c = new TCanvas("c","",800,600);

  c->Divide(2,3);
  c->cd(1);gPad->SetLogy();hTheta_rat->Draw();hTheta_ratCuts->Draw("sames");legend1->Draw();
  c->cd(2);gPad->SetLogy();hPhi_rat->Draw();hPhi_ratCuts->Draw("sames");legend1->Draw();

  c->cd(3);gPad->SetLogy();hThetaCor->Draw();hThetaCor_bbCuts->Draw("sames");legend2->Draw();
  c->cd(4);gPad->SetLogy();hPhiCor->Draw();hPhiCor_bbCuts->Draw("sames");legend22->Draw();

  c->cd(5);gPad->SetLogy();hThetaCor_bbCuts->Draw();hTheta_ratCuts->Draw("sames");legend3->Draw();
  c->cd(6);gPad->SetLogy();hPhiCor_bbCuts->Draw();hPhi_ratCuts->Draw("sames");legend33->Draw();

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,2);
  c1->cd(1);gPad->SetLogz();hPhiTheta_rat->Draw("colz");
  c1->cd(2);gPad->SetLogz();hPhiThetaCor->Draw("colz");
  c1->cd(3);gPad->SetLogz();hPhiTheta_ratCuts->Draw("colz");
  c1->cd(4);gPad->SetLogz();hPhiThetaCor_bbCuts->Draw("colz");
}
