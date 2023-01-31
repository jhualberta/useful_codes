{

   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("plotDustHist.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open("Sample4_center1.dat");//Results_sample4_8Nov_center_10times_Retuned_mmUnits.csv");

   double pyData[] = {1070.2692, 402.688, 166.0604, 96.8,73.4712, 42.108, 33.88,31.5084, 27.6848, 23.3772, 22.6512, 21.6832, 21.296, 20.328, 18.1016, 17.424, 16.456, 16.456, 15.6816, 14.8104, 12.1968, 11.0352, 10.6964, 10.164, 9.4864, 9.4864, 8.8088, 8.8088, 8.8088, 8.1796, 8.1312, 7.5504, 6.9696, 6.9212, 6.3888, 6.3888, 5.8564, 5.8564};

   double pyData_sobel[] = {3.6378529112455853, 2.076392583388519, 0.8789219346619663, 0.8109698676000323, 0.7120735950346119, 0.6688766733765055, 0.6564994166285977, 0.6415639670775182, 0.5596223752952112, 0.47767348904401463, 0.46341135291430413, 0.448696113070137, 0.43434085896402985, 0.42303475232525306, 0.408690398178412, 0.408690398178412, 0.3410615209733829, 0.313730969018343};
   double pyData_canny[] = {7.1588210818049065, 1.5552405483182943, 1.5252330878215354, 1.2287422236364127, 1.187294050331903, 0.5338676540950633, 0.5192383591354462, 0.48866765456195255, 0.4726413842690717, 0.41166513503605107, 0.4022003035287298, 0.3925073055536097, 0.3925073055536097, 0.37236512514151915, 0.37236512514151915, 0.3618738553363748, 0.3164495065463376, 0.2775445774222176, 0.23221045351046138, 0.21498510523728476, 0.21498510523728476, 0.19625365277680484, 0.12412170838050639, 0.08776730168831519, 0.0};

   TH1F *hradius = new TH1F("hradius","",100,0,50);

   TH1F *hradius1 = new TH1F("hradius1","",100,0,50);

   TH1F *hradius2 = new TH1F("hradius2","",100,0,50);

   TH1F *hmass = new TH1F("hmass","",100,0,50);

   double totalArea = 0.8*1.2;
   // 2720 x 1824 pixels
   // 0.22 um/px

   double num, area, thresh1, thresh2, thresh3;
   int counts = 0;
   while (1) {
      in >> num >> area>> thresh1>> thresh2>>thresh3;
      if (!in.good()){
	      cout<<"end of the file"<<endl;
	      break;
      }
      // if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
      //
      // area = area*(1000*1000); // mm^2 to um^2
      area = area*0.22*0.22;// px^2 to um^2
      double radius = sqrt(area/TMath::Pi());
      cout<<radius<<",";
      hradius->Fill(radius);
      double volume = 4./3*TMath::Pi()*pow(radius,3);
      double mass = volume*0.001;
      hmass->Fill(radius,mass);
      counts++;
   }

   for(int i = 0; i<sizeof(pyData_sobel)/sizeof(double); i++)
   {
     double xx = pyData_sobel[i];
     hradius1->Fill(xx);
   } 

   for(int i = 0; i<sizeof(pyData_canny)/sizeof(double); i++)
   {  
     double xx = pyData_canny[i];
     hradius2->Fill(xx);     
   }
   hradius1->SetLineColor(kRed);
   hradius2->SetLineColor(kBlue);

   printf(" found %d dusts\n",counts);
   cout<<"Counts "<<hradius->Integral()<<endl;
   cout<<"Density "<<hmass->Integral()/(totalArea/100)<<" ng/cm^2"<<endl;

   hradius->SetTitle("For area = 1.2 mm x 0.8 mm");
   hradius->GetYaxis()->SetTitle("Number of particles (/0.5 #mum)");
   hradius->GetXaxis()->SetTitle("Radius of dust particles [#mum]");

   hmass->GetYaxis()->SetTitle("Total mass contribution (ng/0.5 #mum)");
   hmass->GetXaxis()->SetTitle("Radius of dust particles [#mum]");

   TCanvas *c1 = new TCanvas("c1","",1000,600);
   c1->Divide(2,1);
   c1->cd(1);
   hradius->Draw();
   hradius1->Draw("same");
   hradius2->Draw("same");
   c1->cd(2);
   hmass->Draw();
  

   in.close();

}
