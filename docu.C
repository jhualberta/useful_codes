{
  vector<double> b;
  for(int i=0;i<10;i++)
  {
   b.push_back(i);
  }

  ofstream ff;
  ff.open("values.txt");
  

  for(vector<double>::iterator it = b.begin(); it != b.end();++it)
  {
     ff<<*it;
     ff<<" ";
   }
   
   ff.close();

}
