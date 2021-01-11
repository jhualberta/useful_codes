{
  TString fname = "Merged_Analysis10_r251230-251253_p000_1sep2019.root";
  TFile *ff = new TFile(fname,"read");

  TH1F *hskyShine = new TH1F("hskyShine","sky shine", 1000, 0, 100);
  
  TH2F *hxy = new TH2F("hxy","x vs y, nhitClean>10",2000,-9000,9000,2000,-9000,9000); 
  TH2F *hyz = new TH2F("hyz","y vs z, nhitClean>10",2000,-9000,9000,2000,-9000,9000); 

  TH2F *hyz_leon = new TH2F("hyz_leon","y vs z, nhitClean>10, Leon",2000,-6000,6000,2000,5100,6000); 

  TH2F *hxy_sky1 = new TH2F("hxy_sky1","x vs y, skyShine>1",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_sky1 = new TH2F("hyz_sky1","y vs z, skyShine>1",2000,-9000,9000,2000,-9000,9000);

  TH2F *hxy_sky1p6 = new TH2F("hxy_sky1p6","x vs y, skyShine>1.6",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_sky1p6 = new TH2F("hyz_sky1p6","y vs z, skyShine>1.6",2000,-9000,9000,2000,-9000,9000);

  TH2F *hskyVsNhits = new TH2F("hskyVsNhits","skyShine vs nhitsCleaned", 1000, 0 , 1000, 1000, 0, 100);

  TH2F *hskyVsZ = new TH2F("hskyVsZ","skyShine vs posZ", 2000, -9000, 9000, 1000, 0 , 100);
  TH2F *hskyVsZ_Nhit50 = new TH2F("hskyVsZ_Nhit50","skyShine vs posZ, nhitsClean>50", 2000, -9000, 9000, 1000, 0 , 100);
  TH2F *hskyVsZ_Nhit100 = new TH2F("hskyVsZ_Nhit100","skyShine vs posZ, nhitsClean>100", 2000, -9000, 9000, 1000, 0 , 100);

  TH2F *hskyVsR = new TH2F("hskyVsR","skyShine vs posR", 1000, 0, 9000, 1000, 0 , 100);
  TH2F *hskyVsR_Nhit50 = new TH2F("hskyVsR_Nhit50","skyShine vs posR, nhitsClean>50", 1000, 0, 9000, 1000, 0 , 100);
  TH2F *hskyVsR_Nhit100 = new TH2F("hskyVsR_Nhit100","skyShine vs posR, nhitsClean>100", 1000, 0, 9000, 1000, 0 , 100);

  TH2F *hxy_sky1p6Nhit50 = new TH2F("hxy_sky1p6Nhit50","x vs y, skyShine>1.6, nhitClean>50",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_sky1p6Nhit50 = new TH2F("hyz_sky1p6Nhit50","y vs z, skyShine>1.6, nhitClean>50",2000,-9000,9000,2000,-9000,9000);

  TH2F *hxy_sky1p6Nhit100 = new TH2F("hxy_sky1p6Nhit100","x vs y, skyShine>1.6, nhitClean>100",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_sky1p6Nhit100 = new TH2F("hyz_sky1p6Nhit100","y vs z, skyShine>1.6, nhitClean>100",2000,-9000,9000,2000,-9000,9000);

  TH2F *hxy_sky2 = new TH2F("hxy_sky2","x vs y, skyShine>2",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_sky2 = new TH2F("hyz_sky2","y vs z, skyShine>2",2000,-9000,9000,2000,-9000,9000);

  TH2F *hxy_nhit50 = new TH2F("hxy_nhit50","x vs y, nhitClean>50",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_nhit50 = new TH2F("hyz_nhit50","y vs z, nhitClean>50",2000,-9000,9000,2000,-9000,9000);

  TH2F *hRhoZ = new TH2F("hRhoZ","rho vs z",1000,0,9000,2000,-9000,9000);
  TH2F *hRhoZ_nhit50 = new TH2F("hRhoZ_nhit50","rho vs z, nhitClean>50",1000,0,9000,2000,-9000,9000);
  TH2F *hRhoZ_nhit100 = new TH2F("hRhoZ_nhit100","rho vs z, nhitClean>100",1000,0,9000,2000,-9000,9000);

  TH2F *hxy_nhit100 = new TH2F("hxy_nhit100","x vs y, nhitClean>100",2000,-9000,9000,2000,-9000,9000);
  TH2F *hyz_nhit100 = new TH2F("hyz_nhit100","y vs z, nhitClean>100",2000,-9000,9000,2000,-9000,9000);

  /// nhitCleaned vs positions
  TH2F *hnhitCleanVsX = new TH2F("hnhitCleanVsX","nhitCleaned vs X",1000,0,1000,2000,-9000,9000);
  TH2F *hnhitCleanVsY = new TH2F("hnhitCleanVsY","nhitCleaned vs Y",1000,0,1000,2000,-9000,9000);
  TH2F *hnhitCleanVsZ = new TH2F("hnhitCleanVsZ","nhitCleaned vs Z",1000,0,1000,2000,-9000,9000);
  TH2F *hnhitCleanVsR = new TH2F("hnhitCleanVsR","nhitCleaned vs R",1000,0,1000,1000,0,9000);

  TH2F *hnhitCleanVsX_sky1p6 = new TH2F("hnhitCleanVsX_sky1p6","nhitCleaned vs X, skyShine>1.6",1000,0,1000,2000,-9000,9000);
  TH2F *hnhitCleanVsY_sky1p6 = new TH2F("hnhitCleanVsY_sky1p6","nhitCleaned vs Y, skyShine>1.6",1000,0,1000,2000,-9000,9000);
  TH2F *hnhitCleanVsZ_sky1p6 = new TH2F("hnhitCleanVsZ_sky1p6","nhitCleaned vs Z, skyShine>1.6",1000,0,1000,2000,-9000,9000);
  TH2F *hnhitCleanVsR_sky1p6 = new TH2F("hnhitCleanVsR_sky1p6","nhitCleaned vs R, skyShine>1.6",1000,0,1000,1000,0,9000);

  TTree* output = (TTree*)ff->Get("output");
  output->Project("hskyShine","skyShine","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10");

  output->Project("hskyVsNhits","skyShine:nhitsCleaned", "partialFit && ( ((dcApplied & 0x210000000242) & dcFlagged) == (dcApplied & 0x210000000242) )");

  output->Project("hskyVsZ","skyShine:posz", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10");
  output->Project("hskyVsZ_Nhit50","skyShine:posz", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50");
  output->Project("hskyVsZ_Nhit100","skyShine:posz", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100");

  output->Project("hskyVsR","skyShine:posr", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10");
  output->Project("hskyVsR_Nhit50","skyShine:posr", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50");
  output->Project("hskyVsR_Nhit100","skyShine:posr", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100");

  output->Project("hRhoZ","posz:sqrt(posx*posx+posy*posy)", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10");
  output->Project("hRhoZ_Nhit50","posz:sqrt(posx*posx+posy*posy)", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50");
  output->Project("hRhoZ_Nhit100","posz:sqrt(posx*posx+posy*posy)", "partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100");

  output->Project("hyz","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10"); 
  output->Project("hxy","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10"); 

  output->Project("hyz_sky1","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>1");
  output->Project("hxy_sky1","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>1"); 

  output->Project("hyz_sky1p6","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>1.6");
  output->Project("hxy_sky1p6","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>1.6");

  output->Project("hyz_sky1p6Nhit50","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50 && skyShine>1.6");
  output->Project("hxy_sky1p6Nhit50","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50 && skyShine>1.6");

  output->Project("hyz_sky1p6Nhit100","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100 && skyShine>1.6");
  output->Project("hxy_sky1p6Nhit100","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100 && skyShine>1.6");

  output->Project("hyz_sky2","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>2");
  output->Project("hxy_sky2","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>2");
  ///
  output->Project("hyz_nhit50","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50");
  output->Project("hxy_nhit50","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>50");

  output->Project("hyz_nhit100","posz:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100");
  output->Project("hxy_nhit100","posx:posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>100");

  /// nhitClean vs position
  output->Project("hnhitCleanVsX","posx:nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242))");
  output->Project("hnhitCleanVsY","posy:nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242))");
  output->Project("hnhitCleanVsZ","posz:nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242))");
  output->Project("hnhitCleanVsR","sqrt(posx*posx+posy*posy+posz*posz):nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242))");

  output->Project("hnhitCleanVsX_sky1p6","posx:nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && skyShine>1.6");
  output->Project("hnhitCleanVsY_sky1p6","posy:nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && skyShine>1.6");
  output->Project("hnhitCleanVsZ_sky1p6","posz:nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && skyShine>1.6");
  output->Project("hnhitCleanVsR_sky1p6","sqrt(posx*posx+posy*posy+posz*posz):nhitsCleaned","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && skyShine>1.6");

  /// Leon analysis
  output->Project("hyz_leon","(posz-108):posy","partialFit && (((dcApplied & 0x210000000242) & dcFlagged ) == (dcApplied & 0x210000000242)) && nhitsCleaned>10 && skyShine>1 && sqrt(posx*posx+posy*posy+(posz-108)*(posz-108))<6000 && (posz-108)>5100 && posx<100 && posx>-100");

  TFile *fnew = new TFile("saveHisto_"+fname, "recreate");
  fnew->cd();
  //  hyz->Draw("colz");
  hskyShine->Write(); hskyVsNhits->Write();hskyVsZ->Write();hskyVsZ_Nhit50->Write(); hskyVsZ_Nhit100->Write();
  hskyVsR->Write();hskyVsR_Nhit50->Write(); hskyVsR_Nhit100->Write();
  hxy->Write();hyz->Write();hxy_sky1->Write();hyz_sky1->Write();hxy_sky1p6->Write();hyz_sky1p6->Write();hxy_sky2->Write();hyz_sky2->Write();
  hxy_sky1p6Nhit50->Write();hyz_sky1p6Nhit50->Write();hxy_sky1p6Nhit100->Write();hyz_sky1p6Nhit100->Write();
  hxy_nhit50->Write();hyz_nhit50->Write();hxy_nhit100->Write();hyz_nhit100->Write();
  hRhoZ->Write();hRhoZ_nhit50->Write();hRhoZ_nhit100->Write();
 
  hyz_leon->Write(); 
  hnhitCleanVsX->Write();hnhitCleanVsY->Write();hnhitCleanVsZ->Write();hnhitCleanVsR->Write();
  hnhitCleanVsX_sky1p6->Write();hnhitCleanVsY_sky1p6->Write();hnhitCleanVsZ_sky1p6->Write();hnhitCleanVsR_sky1p6->Write();
  fnew->Close();

}
