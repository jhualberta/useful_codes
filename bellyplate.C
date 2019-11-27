{
  const double avRadius = 6005.3;
  const double bellyRadius = 5978.4; 
  const double halfSideMax = 761.2;
  const double halfSideMin = 650.8;


  double cutMax = sqrt( avRadius * avRadius - halfSideMax * halfSideMax );
  double cutMin = sqrt( bellyRadius * bellyRadius - halfSideMin * halfSideMin );
  double rPolyInner = cutMin - ( cutMax - cutMin ) * halfSideMin / ( halfSideMax - halfSideMin );
  double halfHeight = halfSideMin * ( avRadius - rPolyInner ) / ( cutMin - rPolyInner );
 
  cout<<halfHeight<<endl;





}
