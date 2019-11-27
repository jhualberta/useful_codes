{
  TStyle *snoStyle= new TStyle("snoplus","SNO+ plots style for publications");

  // use plain black on white colors
  snoStyle->SetFrameBorderMode(0);
  snoStyle->SetCanvasBorderMode(0);
  snoStyle->SetPadBorderMode(0);
  snoStyle->SetPadBorderSize(0);
  snoStyle->SetPadColor(0);
  snoStyle->SetCanvasColor(0);
  snoStyle->SetTitleColor(0);
  snoStyle->SetStatColor(0);
  snoStyle->SetFillColor(0);

  // use bold lines 
  snoStyle->SetHistLineWidth(2);
  snoStyle->SetLineWidth(2);

  // no title, stats box or fit as default
  snoStyle->SetOptTitle(0);
  //snoStyle->SetOptStat(0);
  //snoStyle->SetOptFit(0);
  
  // postscript dashes
  snoStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // text style and size
  //snoStyle->SetTextFont(132);
  //snoStyle->SetTextSize(0.24);
  snoStyle->SetLabelOffset(0.01,"x");
  snoStyle->SetTickLength(0.015,"x");
  snoStyle->SetTitleOffset(1.5,"x");
  snoStyle->SetLabelOffset(0.01,"y");
  snoStyle->SetTickLength(0.015,"y");
  snoStyle->SetTitleOffset(1.5,"y");
  snoStyle->SetLabelOffset(0.01,"z");
  snoStyle->SetTickLength(0.015,"z");
  snoStyle->SetTitleOffset(1.5,"z");
  snoStyle->SetLabelFont(132,"x");
  snoStyle->SetLabelFont(132,"y");
  snoStyle->SetLabelFont(132,"z");
  snoStyle->SetTitleFont(132,"x");
  snoStyle->SetTitleFont(132,"y");
  snoStyle->SetTitleFont(132,"z");
  snoStyle->SetLabelSize(0.04,"x");
  snoStyle->SetTitleSize(0.05,"x");
  snoStyle->SetTitleColor(1,"x");
  snoStyle->SetLabelSize(0.04,"y");
  snoStyle->SetTitleSize(0.05,"y");
  snoStyle->SetTitleColor(1,"y");
  snoStyle->SetLabelSize(0.04,"z");
  snoStyle->SetTitleSize(0.05,"z");
  snoStyle->SetTitleColor(1,"z");
  snoStyle->SetPadTickX(1);
  snoStyle->SetPadTickY(1);

  // AXIS OFFSETS
  snoStyle->SetTitleOffset(0.8,"x");
  snoStyle->SetTitleOffset(0.8,"y");
  snoStyle->SetTitleOffset(0.8,"z");

  // Legends
  snoStyle->SetLegendBorderSize(0);
  snoStyle->SetLegendFont(132);
  snoStyle->SetLegendFillColor(0);
    
  // graphs - set default marker to cross, rather than .
  snoStyle->SetMarkerStyle(21);  // filled square not .
    
  // SNO+ Preliminary label
  snoStyle->SetTextFont(132);
  snoStyle->SetTextSize(0.06);


  gROOT->SetStyle("snoplus");

  TF1 *f = new TF1("f","gaus",0,1);f->SetParameters(1,0,0.5);f->GetXaxis()->SetTitle("#beta_{14}");f->GetYaxis()->SetTitle("R^{3} / R_{AV}^{3}");f->SetLineColor(kRed);f->Draw();TH1F h("h","h",10,0,1);h.FillRandom("f",100);h.SetLineColor(kBlue);h.Scale(0.05);h.Draw("same");float a[4]={0.1,0.4,0.6,0.8}; float b[4]={0.6,0.5,0.4,0.3}; float berr[4]={0.05,0.05,0.05,0.05};TGraphErrors g(4,a,b,0,berr);g.SetName("g");g.SetMarkerColor(kBlack);g.Draw("P");TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);leg->AddEntry("f","Gaussian Fit","L");leg->AddEntry("h","MC histo","L"); leg->AddEntry("g","Data points","PL"); leg->Draw();TText *tPrelim=new TText(0.6,0.7,"SNO+ Preliminary");tPrelim->Draw();

}

// Cut and paste this line to test:
//TF1 *f = new TF1("f","gaus",0,1);f->SetParameters(1,0,0.5);f->GetXaxis()->SetTitle("#beta_{14}");f->GetYaxis()->SetTitle("R^{3} / R_{AV}^{3}");f->SetLineColor(kRed);f->Draw();TH1F h("h","h",10,0,1);h.FillRandom("f",100);h.SetLineColor(kBlue);h.Scale(0.05);h.Draw("same");float a[4]={0.1,0.4,0.6,0.8}; float b[4]={0.6,0.5,0.4,0.3}; float berr[4]={0.05,0.05,0.05,0.05};TGraphErrors g(4,a,b,0,berr);g.SetName("g");g.SetMarkerColor(kBlack);g.Draw("P");TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);leg->AddEntry("f","Gaussian Fit","L");leg->AddEntry("h","MC histo","L"); leg->AddEntry("g","Data points","PL"); leg->Draw();TText *tPrelim=new TText(0.6,0.7,"SNO+ Preliminary");tPrelim->Draw();
