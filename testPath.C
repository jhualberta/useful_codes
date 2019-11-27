// random PMT pos, random vertex
{
 TRandom *r1 = new TRandom3();
 TRandom *r2 = new TRandom3();
 TRandom *r3 = new TRandom3();
 TRandom *r4 = new TRandom3();
 TRandom *r5 = new TRandom3();
 double a1 = 0, a2 = 0;
 double fAVRadius =6050; double P = 8390;
 double fWaterLevel = 0;
 double fSpeedOfLightWater = 299./1.38486, fSpeedOfLightScint = 299./1.64307;
 double fZoff = 0, rpsup = 8390;

 for(int i = 0;i<1;i++)
 {
   double ran0 = r1.Uniform(0,1);
   double rr = pow(ran0, 1.0/3.0) * 8390.0;
   double ran2Pi = r2.Uniform( 0.0, 2.0*TMath::Pi());

   double theta = r3.Uniform( 0.0, TMath::Pi());
   double phi = r4.Uniform( 0.0, 2*TMath::Pi());
 
   TVector3 pmt, pos, ud;
   double costheta = r1.Uniform(-1,1);
   double sintheta = sqrt(1.0 - costheta*costheta);

   //pos.SetXYZ(rr*sintheta*cos(ran2Pi),rr*sintheta*sin(ran2Pi),rr*costheta);
   pos.SetXYZ(0,0,7000);
   //pmt.SetMagThetaPhi(P,theta,phi);
   pmt.SetXYZ(-5497.21, -5470.37, -3201.23);

   double par[3] = {pos.X(), pos.Y(), pos.Z() };
   double vtxmag = pos.Mag();

   double distInScint = 0, tof = 0;
   double dr = (pmt-pos).Mag();
   double dx = (pmt-pos).X();
   double dy = (pmt-pos).Y();
   double dz = (pmt-pos).Z();
   ud = (pmt-pos).Unit();
   cout<<"pmt ("<<pmt.X()<<", "<<pmt.Y()<<", "<<pmt.Z()<<") "<<endl;
   cout<<"vtx ("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<") "<<endl;
   double sqrVal = pow((pos * ud),2) - (pos).Mag2() + fAVRadius*fAVRadius;
   double udDotVtx = pos*ud;
   double fz = pos.Z();
   double vtx2 = pos*pos; // trialVertex*trialVertex
   double sqrVal = udDotVtx*udDotVtx - vtx2 + fAVRadius*fAVRadius;

 /// evaluate line-neck cylinder intersection
   double scintpathInNeck = 0;
   double fNeckRadius = 785; //outer radius of neck
   double neckZlo = 6008.04 + fZoff; //sqrt(6060**2 - 785**2)
   double neckZhi = rpsup; // neck region inside PSUP
   double bneck = dr*(par[0]*dx+par[1]*dy);
   double sqrValneck = bneck*bneck - ( par[0]*par[0]+par[1]*par[1]-fNeckRadius )*vtx2*( dx*dx + dy*dy );

   if( sqrValneck > 1e-3 ) // check the line pass through the neck
   {
       double aneck1 = ( -bneck + sqrt(sqrValneck) )/( dx*dx+dy*dy );
       double aneck2 = ( -bneck - sqrt(sqrValneck) )/( dx*dx+dy*dy );

       if( aneck1*aneck2<0 ) {
         double aneckplus = ( aneck1>0 || aneck2>0 )? aneck1:aneck2;
         double zplus = fz + aneckplus*dz/dr;
         if( zplus < neckZhi && zplus > neckZlo ) {
             scintpathInNeck = aneckplus;
           }
         if( zplus < neckZlo ) {
           if( vtxmag > fAVRadius ) {
             if( sqrVal > 1e-3 ) {
               // cout<<" eh !!!!!!!"<<endl;
               double a1 = -udDotVtx + sqrt(sqrVal);
               double a2 = -udDotVtx - sqrt(sqrVal);
               double asmall = a1<a2? a1:a2;
               scintpathInNeck = asmall;
             }
           }
         }
       }

       if( aneck1>0 && aneck2>0 ) {
         double aneckbig = aneck1>aneck2? aneck1:aneck2;
         double anecksmall = aneck1<aneck2? aneck1:aneck2;
         double zbig = fz + aneckbig*dz/dr;
         double zsmall = fz + anecksmall*dz/dr;
         if( zbig<neckZhi && zbig>neckZlo && zsmall<neckZhi && zsmall>neckZlo )
         {
           scintpathInNeck = aneckbig - anecksmall;
         }

         if( zbig<neckZlo && zsmall>neckZlo && zsmall<neckZhi ) {
           if( sqrVal > 1e-3 ) { 
             double a1 = -udDotVtx + sqrt(sqrVal);
             double a2 = -udDotVtx - sqrt(sqrVal);
             double asmall = a1<a2? a1:a2;
             scintpathInNeck = asmall - anecksmall;
           }
         }
         if( zsmall<neckZlo && zbig>neckZlo && zbig<neckZhi) {
           if( sqrVal > 1e-3 ) { 
             double a1 = -udDotVtx + sqrt(sqrVal);
             double a2 = -udDotVtx - sqrt(sqrVal);
             if(a1*a2<0) {
               double aplus = (a1>0 || a2>0)? a1:a2;
               scintpathInNeck = aneckbig - aplus;
             }
             if(a1>0 && a2>0) {
               double abig = a1>a2? a1:a2;
               scintpathInNeck = aneckbig - abig;
             }
           }
         }
     }
    } // line pass the neck

   cout<<"scint in neck "<<scintpathInNeck<<endl;
   if( sqrVal > 1e-3 ) // check the line pass through the detector, also avoid sqrVal close to 0
   {
     /// find the line-sphere intersect points; a1, a2 are the path lengths btw vertex and intersection points
     double a1 = -udDotVtx + sqrt(sqrVal);
     double a2 = -udDotVtx - sqrt(sqrVal);
     /// find the line-plane intersect point; a3 is the path length btw vertex and intersection point
     double a3 = -9999; //( fWaterLevel-fz )*dr/dz; // not defined and not hit interface
     if( fabs( dz )>1e-3 )
     { a3 = ( fWaterLevel-fz )*dr/dz; } // a3 is defined
     double aplus = 0;
     cout<<"a1 a2 a3 "<<a1<<" "<<a2<<" "<<a3<<endl;

     if( a1*a2<0 ) { // vertex inside the AV
       aplus = ( a1>0 || a2>0 )? a1:a2;
       if( a3<=0 ) { // not hit interface plane
         if( fz>=fWaterLevel ) distInScint = aplus; // vertex must above
       }
       else { // a3>0, hit interface plane
         if( fz>=fWaterLevel ) { // vertex above
           if( a3>=aplus ) distInScint = aplus;
           else { distInScint = aplus - a3; } // vertex below
         }
         else {
           if( a3<aplus ) distInScint = aplus - a3;
         }
       }
     } // a1*a2<0

     if( a1>0 && a2>0) { // vertex in external
       double abig = a1>a2? a1:a2;
       double asmall = a1<a2? a1:a2;
       if( (a3>abig && fz>=fWaterLevel) || (a3!= -9999 && a3<asmall && fz<fWaterLevel) || (a3 == -9999 && fz>fWaterLevel) )
         distInScint = abig - asmall;
       if( a3<abig && a3>asmall ) {
         if( fz<fWaterLevel ) { cout<<" en "<<endl; distInScint = abig - a3; }
         else distInScint = a3 - asmall;
       }
     } // a1>0 and a2>0
   } // pass through AV

   tof = ( dr - distInScint)/fSpeedOfLightWater + distInScint/fSpeedOfLightScint;
   cout<<" tof = ("<<dr<<" - "<<distInScint<<")/"<<fSpeedOfLightWater<<" + "<<distInScint<<"/"<<fSpeedOfLightScint<<" = "<<tof<<" sp/dr "<<distInScint/dr<<endl;

 }// for loop
}  
