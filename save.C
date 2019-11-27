{ 
   int files=30;//20 files
   char filename[100];
   TKey *key;
   TObjArray* hcObj = new TObjArray();
   for(int j=0; j<files+1; j++) 
   {
	     sprintf(filename,"dumpfit-NOset-10MeV-%d.root",j);
	     TFile *f = new TFile(filename,"read"); 
	     if (!f->IsOpen()) 
	     {printf("dumpfit-NOset-10MeV-%d.root, neglect this file\n",j);}
		 printf("dumpfit-NOset-10MeV-%d.root is loaded\n",j);

		 TH1F *hhc = (TH1F*)f.Get("htResVsX");
		 
		 hcObj->AddLast(hhc);  	                    	     
    }   
   TFile *ftot=new TFile("tResVsX.root","recreate");     
   TH1F *hPMT = (TH1F*)f.Get("hPmtX");
   ftot->cd();
   hPMT->Write();//add an extra branch
 
   while(object=next())
   {  
     TH1F* hc = (TH1F*) object; 
     ftot->cd();
     hc->Write();
   }
  
   delete hcObj;	
   ftot->Close();
}
