{

   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("plotDustHist.C","");
   dir.ReplaceAll("/./","/");
   ifstream in, in1;
   in.open("combined_sample5.dat");

   in1.open("combinedClean_sample5.dat");

   TH1F *hradius = new TH1F("hradius","",100,0,50);
   TH1F *hmass = new TH1F("hmass","Total mass contribution",100,0,50);

   TH1F *hradius1 = new TH1F("hradius1","after clean",100,0,50);
   TH1F *hmass1 = new TH1F("hmass1","Total mass contribution",100,0,50);

   // 2720 x 1824 pixels
   // 0.22 um/px
   // 1 mm = 1000 um, 1 cm = 10000 um, 1 cm^3 = 1e12 um^3
   // 1 g/cm^3 = 0.001 ng/um^3
   // 4560.0987 pixels/mm = 4.56 pixels/um -> 20.8 pixels/ um^2
   double num, area;// thresh1, thresh2, thresh3;
   int counts = 0, counts1 = 0;
   while (1) {
      in >> num >> area; //thresh1>> thresh2>>thresh3;
      if (!in.good()){
	      cout<<"end of the file"<<endl;
	      break;
      }
      area = area/20.8;// or *0.22*0.22 
      double radius = sqrt(area/TMath::Pi()); // px^2 to um^2
      hradius->Fill(radius);
      double volume = 4./3*TMath::Pi()*pow(radius,3);
      double mass = volume*0.001;
      hmass->Fill(radius,mass);
      counts++;
   }

   double num1, area1;
   while (1) {
      in1 >> num1 >> area1;
      if (!in1.good()){
              cout<<"end of the file"<<endl;
              break;
      }
      area1 = area1/20.8;
      double radius = sqrt(area1/TMath::Pi()); // px^2 to um^2
      // cout<<area1<<endl;
      hradius1->Fill(radius);
      double volume = 4./3*TMath::Pi()*pow(radius,3);
      double mass = volume*0.001;
      hmass1->Fill(radius,mass);
      counts1++;
   }

   TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry(hradius,"before cleanning","l");
   legend->AddEntry(hradius1,"after cleanning with UPW and Kimwipe","l");

   TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9);
   legend1->AddEntry(hradius,"before cleanning","l");
   legend1->AddEntry(hradius1,"after cleanning with UPW and Kimwipe","l");

   printf(" found %d dusts\n",counts);
   printf(" after clean, %d dusts\n",counts1);

   TCanvas *c = new TCanvas("c","",800,600);
   hradius->SetTitle("For scanned area = 22.07 mm^{2}");
   hradius->GetYaxis()->SetTitle("Number of particles (/0.5 #mum)");
   hradius->GetXaxis()->SetTitle("radius of dust [#mum]");
   hradius->Draw();
   hradius1->SetLineColor(kRed);
   hradius1->Draw("same");
   legend->Draw("same");

   TCanvas *c1 = new TCanvas("c1","",800,600);
   hmass->SetTitle("For scanned area = 22.07 mm^{2}");
   hmass->GetYaxis()->SetTitle("Total mass contribution (ng/0.5 #mum)");
   hmass->GetXaxis()->SetTitle("Radius of dust [#mum]");
   hmass->Draw();
   hmass1->SetLineColor(kRed);
   hmass1->Draw("same");
   legend1->Draw("same");

   in.close();
   in1.close();
}
