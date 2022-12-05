{

   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("plotDustHistSingleRun.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open("Combined_14Nov.dat");

   TH1F *hradius = new TH1F("hradius","",100,0,50);

   TH1F *hmass = new TH1F("hmass","",100,0,50);

   double totalArea = 0.8*1.2;
   // 2720 x 1824 pixels
   // 0.22 um/px

   double num, area;// thresh1, thresh2, thresh3;
   string fileName;
   int counts = 0;
   while (1) {
      in >> num >> fileName >> area; //>> thresh1>> thresh2>>thresh3;
      if (!in.good()){
	      cout<<"end of the file"<<endl;
	      break;
      }
      // if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
      //
      // area = area*(1000*1000); // mm^2 to um^2
      area = area*0.22*0.22;// px^2 to um^2
      double radius = sqrt(area/TMath::Pi());
      // cout<<radius<<endl;
      hradius->Fill(radius);
      double volume = 4./3*TMath::Pi()*pow(radius,3);
      double mass = volume*0.001;
      hmass->Fill(radius,mass);
      counts++;
   }

   printf(" found %d dusts\n",counts);
   cout<<"Counts "<<hradius->Integral()<<endl;
   cout<<"Density "<<hmass->Integral()/totalArea/100<<" ng/cm^2"<<endl;

   hradius->SetTitle("For area = 1.2 mm x 0.8 mm");
   hradius->GetYaxis()->SetTitle("Number of particles (/0.5 #mum)");
   hradius->GetXaxis()->SetTitle("Radius of dust particles [#mum]");

   hmass->GetYaxis()->SetTitle("Total mass contribution (ng/0.5 #mum)");
   hmass->GetXaxis()->SetTitle("Radius of dust particles [#mum]");

   TCanvas *c1 = new TCanvas("c1","",1000,600);
   c1->Divide(2,1);
   c1->cd(1);
   hradius->Draw();
   c1->cd(2);
   hmass->Draw();
  

   in.close();

}
