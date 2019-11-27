{
  double r = 8390;
  double theta = 0.73;
  double phi = 0.73;
  double theta_f = 2.411865, phi_f = 2.411865;
 
  int total = int((phi_f-phi)/0.01);
  cout<<total<<endl;
  TVector3 cone;
  double cr = (r*r-5621.3*5621.3)/(r*r);

  TH2F *h = new TH2F("h","",total,-TMath::Pi(),TMath::Pi(), 200,-1,1);
  TH2F *h1 = new TH2F("h1","",200,-10,10, 200,-10,10);

//  cone.SetMagThetaPhi(r,theta,phi);
  for(theta=0.73;theta<theta_f;theta+=0.01)
  { 
    phi = acos(sqrt( (cr-cos(theta)*cos(theta)))/sin(theta) ); 
    cone.SetMagThetaPhi(r,theta,phi);
    h->Fill((cone.Unit()).Phi(),(cone.Unit()).CosTheta());
    double chart1_y1 = 2*(cone.Unit()).X()/(1-(cone.Unit()).Z());
    double chart1_y2 = 2*(cone.Unit().Y())/(1-(cone.Unit()).Z());
    h1->Fill(chart1_y1,chart1_y2);
    phi = TMath::Pi()-acos(sqrt( (cr-cos(theta)*cos(theta)))/sin(theta) );
    cone.SetMagThetaPhi(r,theta,phi);
    h->Fill((cone.Unit()).Phi(),(cone.Unit()).CosTheta()); 
    chart1_y1 = 2*(cone.Unit()).X()/(1-(cone.Unit()).Z());
    chart1_y2 = 2*(cone.Unit().Y())/(1-(cone.Unit()).Z());
    h1->Fill(chart1_y1,chart1_y2);
  }
  TCanvas c1("c1","",600,500);
  c1.cd();h->Draw("colz");
  TCanvas c2("c2","",600,500);
  c2.cd();h1->Draw("colz");
}
