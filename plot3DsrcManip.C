{
  gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);
      	// yz-xz scans, x, y, z
  const int N1 = 83;//internal: 21 corner +3 neck + 15 x+ 22 y+ 22 z (+ ex: 19);
  double data_total[N1][3]={{-5.399,1991.039,1988.966},{-9.92,2776.468,-2839.407},{-10.44,3227.323,-3324.9},{-8.377,-2832.328,-2839.69},{-9.284,-3137.206,-3212.182},{-4.719,-2099.246,2110.708},{-121.599,81.021,6496.044},{-123.099,126.177,6997.521},{-124.4,153.386,7499.787},{2824.557,13.507,-2840.127},{3245.741,13.681,-3264.759},{-3250.287,13.438,-3264.698},{-2834.05,13.893,-2839.429},{-2109.283,6.613,2098.574},{-2771.903,16.196,-4773.961},{2771.103,15.596,-4773.273}, {-186.0,254.0,5500.9}, {-186.0,254.0,5999.8}, {-8.765,996.891,0.172}, {-7.518,1990.663,-0.445}, {-7.794,-997.865,0.155}, {-7.083,-1990.651,-0.744}, {3003.718,7.258,-0.69}, {2091.016,6.312,2099.444},\
  {-5.283,-0.209,-1.057},{-4999.043,2.46,-9.899},{-4002.525,5.269,-7.364},{-3004.229,8.101,-2.54},{-2000.155,10.637,-0.361},{-992.994,11.641,0.024},{998.133,10.897,-0.044},{2011.103,9.874,0.057},{4003.323, 4.378, -5.974},{5004.868,2.262,-7.547},{4503.445,3.278,-7.0},{-4489.301,3.918,-10.545},{3511.147,5.742,-2.578},{2502.805,8.681,-2.803},{-3476.709,6.643,-4.886},\
  {-5.995,-0.201,-1.107},{-7.761,-998.068,0.159},{-7.084,-2000.578,-0.716},{-5.491,-2998.017,-4.196},{-3.774,-3992.167,-7.374},{-1.624,-4999.882,-12.012},{-1.745,5002.057,-9.897},{-3.967,3973.021,-7.359},{-5.984,2980.035,-3.441},{-7.952,1986.669,-1.714},{-9.242,994.183,0.553},{-9.867,496.858,0.269},{-8.414,1494.71,0.634},{-6.835,2487.539,-1.126},{-4.949,3496.338,-5.971},{-2.539,4505.371,-7.453},{-7.711,-501.268,0.126},{-7.534,-1494.927,-0.096},{-6.349,-2487.912,-1.434},{-4.366,-3475.769,-6.019},{-2.799,-4498.62,-10.213},{-1.898,-4874.53,-11.077},\
{-186,254,-5501.2},{-186.0,254.0,-4999.899},{-186.0,254.0,-4500.2},{-186.0,254.0,-4001.0},{-186.0,254.0,-3501.399},{-186.0,254.0,-2999.7},{-186.0,254.0,-2500.399},{-186.0,254.0,-1998.499},{-186,254,-1499.5},{-186.0,254.0,-1000.099},{-186.0,254.0,-499.6},{-186,254,0.401},{-186.0,254.0,500.8},{-186.0,254.0,1000.3},{-186.0,254.0,1500.9},{-186.0,254.0,2000.001},{-186.0,254.0,2500.3},{-186.0,254.0,3000.3},{-186.0,254.0,3500.9},{-186.0,254.0,4000.5},{-186.0,254.0,4500.4},{-185.037,247.24,4973.567}};
  const int N2 = 19;
  double data_total2[N2][3] = {{-5861.0,-2524.0,-1.62},{-5861.0,-2524.0,-5000.525},{-5861.0,-2524.0,-4000.021},{-5861.0,-2524.0,-3000.151},{-5861.0,-2524.0,-1999.248},{-5861.0,-2524.0,-998.923},{-5861.0,-2524.0,1000.798},{-5861.0,-2524.0,2000.597},{-5861.0,-2524.0,3000.734},{-5861.0,-2524.0,4000.977},{-5861.0,-2524.0,5000.86},{-5861.0,-2524.0,4498.167},{-5861.0,-2524.0,3498.838},{-5861.0,-2524.0,2498.641},{-5861.0,-2524.0,1498.713},{-5861.0,-2524.0,-1501.717},{-5861.0,-2524.0,-2500.89},{-5861.0,-2524.0,-3500.764},{-5861.0,-2524.0,-4500.812}};

  TH3D *h3D = new TH3D("h3D","",200,-8000,8000,200,-8000,8000,200,-8000,8000);
  TH3D *h3D_ex = new TH3D("h3D_ex","",800,-8000,8000,800,-8000,8000,800,-8000,8000);

  TH2D *h2D = new TH2D("h2D","",160,-8000,8000,160,-8000,8000);

  TH2D *hxy = new TH2D("hxy","",1600,-8000,8000,1600,-8000,8000);
  TH2D *hxz = new TH2D("hxz","",1600,-8000,8000,1600,-8000,8000);
  TH2D *hyz = new TH2D("hyz","",1600,-8000,8000,1600,-8000,8000);

  TH2D *hxy_ex = new TH2D("hxy_ex","",1600,-8000,8000,1600,-8000,8000);
  TH2D *hxz_ex = new TH2D("hxz_ex","",1600,-8000,8000,1600,-8000,8000);
  TH2D *hyz_ex = new TH2D("hyz_ex","",1600,-8000,8000,1600,-8000,8000);

  TNtuple *nt = new TNtuple("nt","","x:y:z");
  TNtuple *nt2 = new TNtuple("nt2","","x:y:z");
  for(int i = 0;i<N1;i++)
  {
      h3D->Fill(data_total[i][0], data_total[i][1], data_total[i][2]);
      //cout<<data_total[i][0]<<", "<<data_total[i][1]<<", "<<data_total[i][2]<<endl;
      nt->Fill(data_total[i][0], data_total[i][1], data_total[i][2]);
  }
  
  for(int i = 0;i<N2;i++)
  {
    h3D_ex->Fill(data_total2[i][0], data_total2[i][1], data_total2[i][2]);  
    nt2->Fill(data_total2[i][0], data_total2[i][1], data_total2[i][2]);
  }

//  h3D->Draw();
  //TGeoManager *geoManager = new TGeoManager("toto","toto");
  //TGeoMaterial *material = new TGeoMaterial("Vaccum",0,0,0);
  //TGeoMedium *medium = new TGeoMedium("Vaccum",1,material);
  //TGeoSphere *worldShape = new TGeoSphere(0.,6000.);
  //// worldShape = ROOT.TGeoCone(10.,0.,0.,0.,10.);
  //TGeoVolume *world = new TGeoVolume("top",worldShape,medium);
  //geoManager->SetTopVolume(world);
  //geoManager->CloseGeometry();
  //world->SetLineColor(8);
  //geoManager->SetTopVisible();
  // world->Draw();
  // world->Draw("same");

  TCanvas c("c","",800,600);
  c.Divide(2,2);
  c.cd(1);
  h3D->Draw();h3D->SetMarkerStyle(8);h3D_ex->SetMarkerSize(0.8);
  //h3D_ex->Draw("same");
  h3D_ex->SetMarkerStyle(21);
  h3D_ex->SetMarkerColor(kRed);

  h3D->GetXaxis()->SetTitle("x [mm]");
  h3D->GetYaxis()->SetTitle("y [mm]");
  h3D->GetZaxis()->SetTitle("z [mm]");
  h3D->GetXaxis()->SetTitleSize(.05);
  h3D->GetYaxis()->SetTitleSize(.05);
  h3D->GetZaxis()->SetTitleSize(.05);
  h3D->GetZaxis()->SetTitleOffset(1.0);

  h3D_ex->Draw("same");

  c.cd(2);
  nt->Draw("y:x>>hxy");hxy->GetXaxis()->SetTitle("x [mm]");hxy->GetYaxis()->SetTitle("y [mm]");hxy->SetMarkerStyle(8);
  hxy->GetXaxis()->SetTitleSize(.06);hxy->GetXaxis()->SetTitleOffset(.7);
  hxy->GetYaxis()->SetTitleSize(.06);hxy->GetYaxis()->SetTitleOffset(.7);
  nt2->Draw("y:x>>hxy_ex","","same");hxy_ex->SetMarkerColor(kRed);hxy_ex->SetMarkerStyle(21);
  
  c.cd(3);
  nt->Draw("z:x>>hxz");hxz->GetXaxis()->SetTitle("x [mm]");hxz->GetYaxis()->SetTitle("z [mm]");hxz->SetMarkerStyle(8);
  hxz->GetXaxis()->SetTitleSize(.06);hxz->GetXaxis()->SetTitleOffset(.7);
  hxz->GetYaxis()->SetTitleSize(.06);hxz->GetYaxis()->SetTitleOffset(.7);

  nt2->Draw("z:x>>hxz_ex","","same");hxz_ex->SetMarkerColor(kRed);hxz_ex->SetMarkerStyle(21);
  c.cd(4);
  nt->Draw("z:y>>hyz");hyz->GetXaxis()->SetTitle("y [mm]");hyz->GetYaxis()->SetTitle("z [mm]");hyz->SetMarkerStyle(8);
  hyz->GetXaxis()->SetTitleSize(.06);hyz->GetXaxis()->SetTitleOffset(.7);
  hyz->GetYaxis()->SetTitleSize(.06);hyz->GetYaxis()->SetTitleOffset(.7);
  nt2->Draw("z:y>>hyz_ex","","same");hyz_ex->SetMarkerColor(kRed);hyz_ex->SetMarkerStyle(21);


}
