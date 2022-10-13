void WaterTankTemp() {
  
  
  //The curl command does all the magic. The dates and sampling frequency can be 
  //varied if you want to be more clever.
  //https://deap-deltav.physics.carleton.ca
  //TString main_command = "curl 'https://deap:3600kgLAr@deap-deltav.physics.carleton.ca:/sample/?tag=DET_OUTER_TEMP/LOWER_RTD_K.CV&start=20161101T120000Z&end=20171031T120000Z&resample_period=60m&timestamp_format=unix&sections=samples' > ";
  //TString main_command = "curl 'https://deap:3600kgLAr@deap-deltav.physics.carleton.ca:/sample/?tag=DET_OUTER_TEMP/LOWER_RTD_K.CV&start=20160101T120000Z&end=20191231T120000Z&resample_period=120m&timestamp_format=unix&sections=samples' > ";
  TString main_command = "curl 'https://deap:3600kgLAr@deap-deltav.physics.carleton.ca:/sample/?tag=DET_OUTER_TEMP/LOWER_RTD_K.CV&start=20201101T120000Z&end=20211118T120000Z&resample_period=120m&timestamp_format=unix&sections=samples' > ";
  TString filename = "TestWaterTemp.txt";
  gSystem->Exec(main_command+filename);
  
  TGraph * gTemp = new TGraph(filename, "%lg %lg"," ");
  gTemp->SetName("WaterTemperature");
  
  RAT::DEAPStyle *fStyle = RAT::DEAPStyle::GetStyle();
  TCanvas * c = new TCanvas("c","c",1000,500);
  c->cd();  
  
  gTemp->SetMarkerStyle(24);
  gTemp->SetMarkerColor(fStyle->Color(1));
  //gTemp->GetXaxis()->SetTimeDisplay(1);
  //gTemp->GetXaxis()->SetTimeFormat("%d\/%m\/%y");  
  //gTemp->GetXaxis()->SetTimeOffset(0,"gmt");
  //gTemp->GetXaxis()->SetNdivisions(507);  
  //gTemp->GetXaxis()->SetTitle("Date (d/m/y)");
  gTemp->GetXaxis()->CenterTitle();
  gTemp->GetYaxis()->CenterTitle();
  gTemp->GetYaxis()->SetTitle("Veto Water Temp. (K)");
  gTemp->SetTitle("");
  gTemp->Draw("ap");  
  
  TFile * fOUT = new TFile("rootfiles/2021Nov_VetoWaterTemp.root", "RECREATE");
  gTemp->Write();
  fOUT->Close();
  /*
  TString main_commandLower = "curl 'https://deap:3600kgLAr@deap-deltav.physics.carleton.ca:/sample/?tag=DET_OUTER_TEMP/LOWER_RTD_K.CV&start=20170801T120000Z&end=20180219T120000Z&resample_period=60m&timestamp_format=unix&sections=samples' > ";
  TString filenameLower = "TestWaterTemp.txt";
  gSystem->Exec(main_commandLower+filenameLower);

  TGraph * gLower = new TGraph(filenameLower,"%lg %lg"," ");
  
  TString main_commandUpper = "curl 'https://deap:3600kgLAr@deap-deltav.physics.carleton.ca:/sample/?tag=DET_OUTER_TEMP/UPPER_RTD_K.CV&start=20170801T120000Z&end=20180219T120000Z&resample_period=60m&timestamp_format=unix&sections=samples' > ";
  TString filenameUpper = "TestWaterTemp.txt";
  gSystem->Exec(main_commandUpper+filenameUpper);

  TGraph * gUpper = new TGraph(filenameUpper,"%lg %lg"," ");  
  
  
  RAT::DEAPStyle *fStyle = RAT::DEAPStyle::GetStyle();
  TCanvas * c = new TCanvas("","",1000,500);
  c->cd();
 
    
  gLower->SetMarkerStyle(24);
  gLower->SetMarkerColor(fStyle->Color(1));
  gLower->GetXaxis()->SetTimeDisplay(1);
  gLower->GetXaxis()->SetTimeFormat("%d\/%m\/%y");  
  gLower->GetXaxis()->SetTimeOffset(0,"gmt");
  gLower->GetXaxis()->SetNdivisions(507);  
  gLower->GetXaxis()->SetTitle("Date (d/m/y)");
  gLower->GetXaxis()->CenterTitle();
  gLower->GetYaxis()->CenterTitle();
  gLower->GetYaxis()->SetTitle("Veto Water Temp. (K)");
  gLower->SetTitle("");
  gLower->Draw("ap");
  
  gUpper->SetMarkerColor(fStyle->Color(2));
  gUpper->SetMarkerStyle(24);
  gUpper->Draw("p same");
  
  TLegend *legend = new TLegend(0.2,0.65,0.46,0.85);
  legend->AddEntry(gUpper,"Upper Sensor","p");
  legend->AddEntry(gLower,"Lower Sensor","p");
  legend->SetFillStyle(0);
  legend->Draw("same");  
  
  c->SaveAs("plots/WaterTemp_Both.pdf");
  
  */
  
}