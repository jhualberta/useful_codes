{
  // sno light water	
  double wave0[] = {200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0};
  double nw0[] = {1.41615, 1.39727, 1.38395, 1.37414, 1.36667, 1.36082, 1.35615, 1.35233, 1.34916, 1.3465, 1.34423, 1.34227, 1.34058, 1.33909, 1.33778, 1.33661, 1.33557, 1.33463, 1.33378, 1.33301, 1.33231, 1.33167, 1.33108, 1.33054, 1.33004, 1.32957, 1.32914, 1.32874, 1.32836, 1.32801, 1.32768};

// snoplus 2018 laser scan
  double wave[] = {200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0};
  double nw[] = {1.41615, 1.39727, 1.38395, 1.37414, 1.36667, 1.36082, 1.35615, 1.35233, 1.34916, 1.3465, 1.34423, 1.34227, 1.34058, 1.33909, 1.33778, 1.33661, 1.33557, 1.33463, 1.33378, 1.33301, 1.33231, 1.33167, 1.33108, 1.33054, 1.33004, 1.32957, 1.32914, 1.32874, 1.32836, 1.32801, 1.32768};
 
 int n0 = sizeof(wave0)/sizeof(double);
 int n = sizeof(wave)/sizeof(double);

 TGraph *gr0 = new TGraph(n0, wave0,nw0);
 gr0->Draw("AP");
 gr0->GetYaxis()->SetTitle("refractive index (n)");
 gr0->GetXaxis()->SetTitle("wavelength [nm]");
// gr0->SetMarkerColor(kBlue);

// TGraph *gr = new TGraph(n, wave,nw);
// gr->Draw("P");
// gr->GetYaxis()->SetTitle("refractive index (n)");
// gr->GetXaxis()->SetTitle("wavelength [nm]");
//



}

