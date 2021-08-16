{
  double radius = 5.5;
  double volume = TMath::Pi()*pow(radius,3)*4./3;
  double mass = volume*0.997/1000;
  double counts = 0, days = 0;
  cout<<"put counts:\n";
  cin>>counts;
  cout<<"put days:\n";
  cin>>days;
  cout<<counts/(mass*days)<<" events/(kt*day)"<<endl;
}
