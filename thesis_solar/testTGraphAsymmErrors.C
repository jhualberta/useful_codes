{

  double x[] = {2.1};
  double y[] = {1};
  double errXleft[] = {0.5};
  double errXright[] = {1.4};

  double errYleft[] = {1.5};
  double errYright[] = {3}; 

  TGraphAsymmErrors *gMy  =     new TGraphAsymmErrors(1,x, y, errXleft,errXright, errYleft, errYright);
  
  gMy->SetMarkerStyle(33);

  gMy->Draw();





}
