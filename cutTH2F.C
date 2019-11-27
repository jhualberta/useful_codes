{
   gROOT->Reset();
   TCutG *mycut = new TCutG("mycut",8);
   mycut->SetVarX("");
   mycut->SetVarY("");
   mycut->SetTitle("Graph");
   mycut->SetFillColor(1);
   mycut->SetPoint(0,-0.646552,0.932203);
   mycut->SetPoint(1,-1.26437,0.105932);
   mycut->SetPoint(2,-0.574713,-1.10169);
   mycut->SetPoint(3,0.948276,-0.338983);
   mycut->SetPoint(4,1.07759,0.720339);
   mycut->SetPoint(5,-0.316092,-0.0847458);
   mycut->SetPoint(6,-0.45977,0.402542);
   mycut->SetPoint(7,-0.646552,0.932203);
   
   TFile f("$ROOTSYS/tutorials/hsimple.root");
   TH2F *hpxpy = (TH2F*)f.Get("hpxpy");
   hpxpy->SetMinimum(0);
   hpxpy->SetMaximum(180);
   hpxpy->SetFillColor(kBlue);
   hpxpy->DrawCopy("lego1");
   hpxpy->SetFillColor(kRed);
   hpxpy->Draw("same,surf1,[mycut]");
}
