{
   int n = 5;
   double concen[] = {0.25, 0.5, 1, 2, 6};
   double fractionalLightYield[] = {0.57, 0.65, 0.9, 1.0, 0.93};
   double rms[]={125.8, 92.48, 67.08, 53.27, 52.11125.8, 92.48, 67.08, 53.27, 52.11};
   double rmsFV[]={120.6, 90.02, 67.39, 55.81, 58.62};

   double radialBias[] ={-1.46, 8.34, 13.4, 15.3, -7.54};
   double radialBiasFV[] ={-7.22,3.3,10.67, 13.37,-10.97};

//   TGraph *gr = new TGraph(n, concen, rms);
//   gr->GetXaxis()->SetTitle("PPO concentration [g/L]");
//   gr->GetYaxis()->SetTitle("z_{fit} resolution [mm]");
//   // gr->SetLineColor(2);
//   gr->SetMarkerStyle(21);
//   gr->Draw("APL");
//
//   TGraph *grFV = new TGraph(n, concen, rmsFV);
//   grFV->GetXaxis()->SetTitle("PPO concentration [g/L]");
//   grFV->GetYaxis()->SetTitle("z_{fit} resolution [mm]");
//   grFV->SetLineColor(2);
//   grFV->SetMarkerStyle(21);
//   grFV->SetMarkerColor(2);
//   grFV->Draw("PL");

    TGraph *gr = new TGraph(n, concen, fractionalLightYield);
    gr->GetXaxis()->SetTitle("PPO concentration [g/L]");
    gr->GetYaxis()->SetTitle("Fractional light yield");
    // gr->SetLineColor(2);
    gr->SetMarkerStyle(21);
    gr->Draw("APL");


//   TGraph *gr1 = new TGraph(n, concen, radialBias);
//   gr1->GetXaxis()->SetTitle("PPO concentration [g/L]");
//   gr1->GetYaxis()->SetTitle("radial bias [mm]");
//   // gr1->SetLineColor(2);
//   gr1->GetYaxis()->SetRangeUser(-20,20);
//   gr1->SetMarkerStyle(21);
//   gr1->Draw("APL");
//
//   TGraph *gr1FV = new TGraph(n, concen, radialBiasFV);
//   gr1FV->GetXaxis()->SetTitle("PPO concentration [g/L]");
//   gr1FV->GetYaxis()->SetTitle("radial bias [mm]");
//   gr1FV->SetLineColor(2);
//   gr1FV->SetMarkerStyle(21);
//   gr1FV->SetMarkerColor(2);
//   gr1FV->Draw("PL");


}
