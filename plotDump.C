{

TFile *fmpw = new TFile("dumpfit_3391.root");
hHx = (TH1F*)fmpw->Get("Hx");
hHy = (TH1F*)fmpw->Get("Hy");
hHz = (TH1F*)fmpw->Get("Hz");
hHt = (TH1F*)fmpw->Get("Ht");

hDHx = (TH1F*)fmpw->Get("sDHx");
hDHy = (TH1F*)fmpw->Get("sDHy");
hDHz = (TH1F*)fmpw->Get("sDHz");
hDHt = (TH1F*)fmpw->Get("sDHt");

hNDHx = (TH1F*)fmpw->Get("sNDHx");
hNDHy = (TH1F*)fmpw->Get("sNDHy");
hNDHz = (TH1F*)fmpw->Get("sNDHz");
hNDHt = (TH1F*)fmpw->Get("sNDHt");

TCanvas *cDeltaRvsTime = new TCanvas("c1","",800,600);
cDeltaRvsTime->Divide(2,2);
cDeltaRvsTime->cd(1);
//gpad->SetLogy();
//gpad->SetGrid();
TAxis *axis1 = hHx->GetXaxis();
//double bmin1 = axis1->FindBin(hHx->Get); 
//double bmax1 = axis1->FindBin(5000); 

TAxis *axis2 = hDHx->GetXaxis();
//double bmin2 = axis2->FindBin(); 
//double bmax2 = axis2->FindBin(5000); 

TAxis *axis3 = hNDHx->GetXaxis();
//double bmin3 = axis3.FindBin(-5000); 
//double bmax3 = axis3.FindBin(5000); 

double scale1 = hHx->GetBinContent(hHx->GetMaximumBin());
//hHx->Scale(1.0/scale1);
double scale2 = hDHx->GetBinContent(hDHx->GetMaximumBin());
hDHx->Scale(scale1/scale2/25);
double scale3 = hNDHx->GetBinContent(hNDHx->GetMaximumBin());
hNDHx->Scale(scale1/scale3);

hHx->Draw();
hHx->SetLineColor(4);
hDHx->Draw("sames");
hNDHx->SetLineColor(2);
hNDHx->Draw("sames");

cDeltaRvsTime->cd(2);

double scale1 = hHy->GetBinContent(hHy->GetMaximumBin());
//hHx->Scale(1.0/scale1);
double scale2 = hDHy->GetBinContent(hDHy->GetMaximumBin());
hDHy->Scale(scale1/scale2);
double scale3 = hNDHy->GetBinContent(hNDHy->GetMaximumBin());

hNDHy->Scale(scale1/scale3);
hHy->Draw();
hHy->SetLineColor(4);
hDHy->Draw("sames");
hNDHy->SetLineColor(2);
hNDHy->Draw("sames");

cDeltaRvsTime->cd(3);
//gpad->SetLogy();
//gpad->SetGrid();

double scale1 = hHz->GetBinContent(hHz->GetMaximumBin());
//hHx->Scale(1.0/scale1);
double scale2 = hDHz->GetBinContent(hDHz->GetMaximumBin());
hDHz->Scale(scale1/scale2);
double scale3 = hNDHz->GetBinContent(hNDHz->GetMaximumBin());

hNDHz->Scale(scale1/scale3);
hHz->Draw();
hHz->SetLineColor(4);
hDHz->Draw("sames");
hNDHz->SetLineColor(2);
hNDHz->Draw("sames");

cDeltaRvsTime->cd(4);

double scale1 = hHt->GetBinContent(hHt->GetMaximumBin());
//hHx->Scale(1.0/scale1);
double scale2 = hDHt->GetBinContent(hDHt->GetMaximumBin());
hDHx->Scale(scale1/scale2);
double scale3 = hNDHt->GetBinContent(hNDHt->GetMaximumBin());

hNDHt->Scale(scale1/scale3);
hHt->Draw();
hHt->SetLineColor(4);
hDHt->Draw("sames");
hNDHt->SetLineColor(2);
hNDHt->Draw("sames");

}
