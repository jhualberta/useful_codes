////////////////////////////////////////////////////////////////////////
///
/// \class
///
/// \brief
///
/// \author David Auty auty@ualberta.ca
/// 
/// REVISION HISTORY: Jie Hu jhu9@ualberta.ca - modified to fit J Tseng's MultiPathFunction\n
///
/// \detail
///       To fit data in a partial scint-water fill geometry
////
////////////////////////////////////////////////////////////////////////

//#include <RAT/MultiPathScintWaterSimple.hh>
//#include <CLHEP/Random/RandFlat.h>
//#include <CLHEP/Units/SystemOfUnits.h>
//#include <CLHEP/Units/PhysicalConstants.h>
//#include <RAT/DB.hh>
//#include <RAT/DU/Utility.hh>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include <TF1.h>
#include "TVector3.h"
double fDepth = 0, fHeight = 0;
double fWaterRI = 1.33, fScintRI = 1.50;
TMatrixT<double> fFactor;
/// 64 pmts
double datapmt[] = {-4702.56, -6472.42, -2582.24, 222.229,3923.81, -260.85, -7430.4, 262.717,1009.98, -1390.03, -8232.57, 260.289,255.73, -5620.57, -6222.77, 258.198,-201.72, -4846.85, -6872.54, 257.384,4385.98, -3265.07, -6406.34, 261.329,6742.59, 1497.64, -4783.5, 265.838,3473.94, -6149.72, -4568.94, 259.593,-7735.95, -3049.68, -1272.46, 266.185,3513.22, -7356.18, -2054.53, 266.069,5103.06, -6487.5, 1640.88, 262.486,4853.79, -6429.6, 2388.87, 264.335,-3293.48, 7157.67, -2955.59, 271.705,6098.11, -4197.66, 4001.81, 265.954,-3682.5, -1662.4, 7378.3, 270.341,-4805.67, -6614.67, -1561.54, 292.955,-1803.11, 7958.15, -2059.09, 295.114,5651.59, 4472.71, 4289.4, 310.686,-5866.02, -2222.3, -5594.72, 358.588,-201.72, -4846.85, -6872.54, 257.384};

int total = sizeof(datapmt)/sizeof(*datapmt)/4;

//double datapmt[]={-4231.94, -358.47, -7253.2, 305.172, -2731.22, -4261.57, -6614.63, 305.632, 1217.42, -3669.07, -7482.95, 304.138, 5222.62, -687.16, -6558.71, 304.943, 4261.5, 2020.19, -6978.37, 305.057, 896.88, 4653.15, -6943.92, 306.782, -2833.32, 4658.62, -6419.47, 302.759, -3684.61, 3430.4, -6732.86, 302.529, -6248.41, -3581.45, -4353.23, 302.988, -1219.86, -6985.48, -4527.16, 302.874, 4097.97, -6006.91, -4250.23, 302.414, 5458.19, -3553.72, -5340.4, 301.609, 7052.97, 400.17, -4522.56, 302.069, 2406.86, 5241.22, -6109.44, 303.218, -1736.35, 6499.08, -5042.19, 303.908, -6898.1, 270.66, -4825.13, 302.414, -5841.29, -5474.19, -2561.43, 302.414, -147.47, -7383.74, -3977.78, 299.659, 2391.01, -7616.28, -2686.29, 303.333, 5728.35, -5165.1, -3331.17, 301.149, 7589.93, 640.41, -3558.01, 301.494, 1508.06, 7862.13, -2561.43, 299.318, -1798.9, 7100.22, -4087.66, 307.586, -6957.04, 2801.1, -3787.44, 304.023, -8248.37, -1011.26, -1149.76, 295.682, -5177.24, -6623.61, -36.57, 295.455, 5177.4, -6623.47, -36.51, 295, 5979.53, -5915.94, -265.15, 297.5, 7700.16, 3289.55, -805.4, 295.909, 2875.3, 7916.3, -75.09, 293.182, -4928.14, 6783.24, -510.4, 296.023, -7807.98, 3157.92, -36.42, 296.25, -8229.58, -289.67, 1718.82, 289.651, -4593.05, -6858.03, 1640.88, 275.909, 5181.09, -6343.84, 1886.78, 273.182, 7684.58, -2496.74, 2063.86, 290.455, 7400.79, 3663.57, 1621.16, 292.5, 462.84, 8177.63, 1886.77, 292.614, -4069.56, 7368.35, 290.55, 292.614, -8188.81, 1163.2, 1395.85, 293.295, -7516.21, -60.97, 3761.85, 287.674, -4459.52, -5886.75, 3977.78, 282.791, 3707.88, -7078.19, 2640.25, 286.744, 5789.72, -5188.77, 3206.77, 287.791, 6976.78, 2421.94, 3977.78, 289.302, 2016.2, 7372.15, 3497.75, 286.977, -4598.71, 6329.68, 3062.23, 286.744, -7022.88, 3194.47, 3344.33, 287.907, -3989.14, -3977.85, 6247.3, 282.442, -2123.9, -5295.74, 6174.12, 287.674, 3302.21, -5467.24, 5455.66, 281.047, 5478.58, -1624.86, 6187.36, 280.349, 6575.03, 1044.13, 5167.46, 279.318, 1677.82, 5856.42, 5791.01, 283.023, -3427.55, 5220.23, 5628.35, 269.422, -5266.46, 1980.06, 6222.8, 281.279, -2642.39, -1708.78, 7790.18, 277.045, -364.87, -2817.92, 7912.2, 274.545, 2967.78, -3304.35, 7148.06, 275.909, 3040.6, -1298.33, 7727.57, 265.029, 1766.56, 1143.91, 8132.05, 272.5, 1659.04, 2283.54, 7928.06, 274.432, -1693.17, 1963.92, 8013.56, 274.886, -3259.67, 1532.51, 7596.23, 276.705};

//double par[4] = {59.0571, 16.7535, 3834.9, 248.026};
double par[4] = {0,0,-1000,0};//0,0,3400,248.026};


double ScintExWaterPathCalculation(double *incidentVtx, double *pmtvtx )
{
  TVector3 startpos, pmtpos, ndiff;
  startpos.SetXYZ(incidentVtx[0],incidentVtx[1],incidentVtx[2]);
  pmtpos.SetXYZ(pmtvtx[0],pmtvtx[1],pmtvtx[2]);
  ndiff = (pmtvtx - startpos).Unit();
  double pathInScint = 0;
  double sqrVal = pow((startpos * ndiff),2) - (startpos).Mag2() + 6050*6050; 
  if(sqrVal<0) pathInScint = 0;
  else {
    pathInScint = -(startpos * ndiff) + sqrt(sqrVal);
    if(pathInScint>6050)
    pathInScint = -(startpos * ndiff) - sqrt(sqrVal);
 
 }
  return pathInScint;
}

double Factor( double transver, double dh ) {

  /// numerically calculate transverse distance
  /// if PMT above water, vertex below, input (transver,dh) = (transver/h, fDepth/fHeight), then:
  /// TransverseDist = fDepth*tan( theta ) + fHeight*tan(asin( fWaterRI/fScintRI*sin( theta ) ))
  /// if PMT below water, vertex above, input (transver,dh) = (transver/-fDepth, fHeight/fDepth), then:
  /// TransverseDist = (-fHeight)*tan( theta ) + (-fDepth)*tan(asin( fWaterRI/fScintRI*sin( theta ) ))
  /// param: fDepth -- [0]; fWaterRI/fScintRI -- [1];
  /// Find theta (always the angle in water)
  TF1 fTransverseDist( "fTransverseDist", "[0]*tan( x )+ 1*tan( asin( [1]*sin( x ) ) )", 0, 0.850908);
 
  fTransverseDist.SetParameter( 0, dh);
  fTransverseDist.SetParameter( 1, fWaterRI/fScintRI );
  
  /// parameter is the vertical distance between PMT and vertex
  double theta = fTransverseDist.GetX(transver);
  /// starting point is below intersection
  if( fabs( transver ) < 0.001 ) return 0;
  /// fraction of straight line distance to point directly below intersection.
  else return dh*tan( theta )/transver;
}

double FactorN( double transver, bool pmtBot ) {
  /// calculates fraction of transverse distance from fDepth to the surface by interpolating, based on fDepth and fHeight
  double alpha = 0, dh = 0, lx = 0;
  if( fabs(fHeight) < 0.1 || fabs(fDepth) < 0.1 ) {
    if (fabs(fHeight) < 0.1)
    { 
      if( pmtBot == false) alpha = 1;
      else alpha = 0;
      //std::cout<<"case 1 "<<std::endl;
    }
    else if (fabs(fDepth) < 0.1) {
      if (pmtBot == false) alpha = 0;
      else alpha = 1;
      //std::cout<<"case 2 "<<std::endl;
    }
  } 
  else {
    if ( pmtBot == false ) { // if the PMT is above the water
      dh = fDepth/fHeight;
      // transver is projection of vector (vertex to PMT) onto detector z-axis
      // lx is the log of the ratio of transver and fHeight to define the extreme conditions
      lx = log(transver/fHeight);
    }
    else if ( pmtBot == true ) { // if the PMT is below the water
      dh = fHeight/fDepth;
      lx = log(transver/(-fDepth));
    }
    //std::cout<<"case 3 "<<std::endl;
    double ldh = log(dh);
    if( lx < -5 ) alpha = ( fWaterRI*dh + fScintRI )/( 1 + dh ); // sVertex close to water level
    else {
      if( lx>5 ) lx = 5; // vetex too low
      if( ldh<-5 ) ldh = -5; // water level too low
      else if( ldh > 5 ) ldh = 5; // water level too high
      double ux = ( lx + 5 )/.2;
      int ix = floor(ux);
      ux = ux-ix;
      double udh = ( ldh + 5 )/.2;
      int idh = floor( udh );
      udh = udh-idh;
      alpha = ( 1-ux )*( ( 1-udh )*fFactor[ix][idh] + udh*fFactor[ix][idh+1] ) + ux*( ( 1-udh )*fFactor[ix+1][idh]+udh*fFactor[ix+1][idh+1]);
    }
  }
  // std::cout<<"factorN  fH, fD "<<fHeight<<" "<<fDepth<<" "<<alpha<<std::endl;
  return alpha;
}


void multiPathLight()
{
  using namespace RAT;
  using namespace ROOT;
  
  double fMaxReflects = 500;
  
  double fSpeedOfLight = 299;
  double fWaterRIeff = 1.53;
  double fScintRIeff = 1.60;
  double fWaterLevel = 0;
  double fWaterRI = 1.33, fScintRI = 1.50;
  std::vector<double> fRWater;
  std::vector<double> fRScint;
  double fSpeedOfLightWater = fSpeedOfLight/fWaterRIeff;
  double fSpeedOfLightScint = fSpeedOfLight/fScintRIeff;

  fFactor.ResizeTo( 50 + 2, 50 + 2 );

  for( int ix = 0; ix < 51; ix++ ) {
    double lx = -5 + ix*0.2;
    double x = exp( lx );
    for( int idh = 0; idh < 51; idh++ ) {
      double ldh = -5 + idh*0.2;
      double dh = exp( ldh );
      double f = Factor( x, dh );
      fFactor[ix][idh] = f;
      //std::cout<<"fFactor["<<ix<<"]["<<idh<<"] = "<<f<<std::endl;
    }
  }


  /// calculate cosTheta, reflection and transimission probabilities (Fresnel equations)
  double step_ct = 1./fMaxReflects;
  for ( double cosTi = 0; cosTi < 1.; cosTi += step_ct ) { // cosThetaIncident
    double Rp = 0., Rs = 0.; // reflection prob for p-wave and s-wave
    double sinTi = sqrt( 1 - cosTi*cosTi ); // sinThetaIncident square
    // incident in water, transmit in scint
    double cosTt = 1 - (fWaterRI/fScintRI*sinTi)*(fWaterRI/fScintRI*sinTi);
    if( cosTt < 0 ) { // total internal reflection
      Rp = 1;
      Rs = 1;
      fRWater.push_back( ( Rp + Rs )/2. );
      //fTWater.push_back( ( 1-Rp + 1-Rs )/2. );
    } else { // transmission(refraction) and reflection
      cosTt = sqrt( cosTt ); // cosThetaTransmit
      Rs = pow( ( fWaterRI*cosTi - fScintRI*cosTt )/( fWaterRI*cosTi + fScintRI*cosTt ), 2 );
      Rp = pow( ( fWaterRI*cosTt - fScintRI*cosTi )/( fWaterRI*cosTt + fScintRI*cosTi ), 2 );
      fRWater.push_back( ( Rp + Rs )/2. );
      //fTWater.push_back( ( 1-Rp + 1-Rs )/2. );
    }
    // incident in scint, transmit in water
    Rs = 0, Rp = 0, cosTt = 0;
    cosTt = 1 - (fScintRI/fWaterRI*sinTi)*(fScintRI/fWaterRI*sinTi);
    if( cosTt < 0 ) { // total internal reflection
      Rp = 1;
      Rs = 1;
      fRScint.push_back( ( Rp + Rs )/2. );
      //fTScint.push_back( ( 1-Rp + 1-Rs )/2. );
    } else {
    Rs = pow( ( fScintRI*cosTi - fWaterRI*cosTt )/( fScintRI*cosTi + fWaterRI*cosTt ), 2 );
    Rp = pow( ( fScintRI*cosTt - fWaterRI*cosTi )/( fScintRI*cosTt + fWaterRI*cosTi ), 2 );
    fRScint.push_back( ( Rp + Rs )/2. );
    //fTScint.push_back( ( 1-Rp + 1-Rs )/2. );
   }
  } // end of for loop
  fRScint.push_back( 0 );
  fRWater.push_back( 0 );

/// loop pmts 
 for(int i = 0;i <total*4; i=i+4)
 {
  double pmt[4] = {datapmt[0+i],datapmt[1+i],datapmt[2+i],datapmt[3+i]};
  std::cout<<pmt[0]<<" "<<pmt[1]<<" "<<pmt[2]<<" "<<pmt[3]<<std::endl; 
  double dx = pmt[0] - par[0];
  double dy = pmt[1] - par[1];
  double dz = pmt[2] - par[2];
  double dr = sqrt(dx*dx + dy*dy + dz*dz);

  double tDiff = pmt[3] - par[3];// tDiff = hitTime - trialTime; fRes = tDiff - TOF

  TVector3 fVertex; // trial vertex
  fVertex.SetXYZ( par[0], par[1], par[2] );
  TVector3 fPMTpos;
  fPMTpos.SetXYZ( pmt[0], pmt[1], pmt[2] );
  TVector3 fPosDiff; // pmtPos - vertex, NOT diff! diff is sBeta
  fPosDiff.SetXYZ( dx, dy, dz );
  fDepth = fWaterLevel - par[2];//fVertex.Z();
  fHeight = pmt[2] - fWaterLevel;
  double fAVRadius = 8390;
  double scintpath = 0;
  if( fWaterLevel >= -fAVRadius && fWaterLevel <= fAVRadius ) { // check valid water height
    double tof = 0, reflected = 0, alpha = 0;// refracted = 0
    if( fDepth>0 ) {  // vertex is below water level
        if( fHeight<0 ) {  //PMT is below water level
         std::cout<<"vbpb"<<std::endl;
        /// geo 1: both PMT and vertex are in water: find direct and reflected light paths 
            tof = dr/fSpeedOfLightWater; // direct time
        }
        else {
        /// geo 2: PMT is above water, vertex below: find refraction path
          std::cout<<"vbpt "<<std::endl;
            alpha = FactorN(fPosDiff.Perp(), false);
            TVector3 s = (1-alpha)*fVertex + alpha*fPMTpos;  //find a point directly below the intersection with water
            s.SetZ(fWaterLevel);  //move intersection to water surface.
            TVector3 incident = s - fVertex;
            TVector3 refractedRay = fPMTpos - s;
            double par1[3] = {s.X(),s.Y(),s.Z()};
            scintpath = ScintExWaterPathCalculation(par1,pmt);
            tof = incident.Mag()/fSpeedOfLightWater + scintpath/fSpeedOfLightScint + (refractedRay.Mag()-scintpath)/fSpeedOfLightWater;
         }
       } else { // vertex is above water level (fDepth<0)
        if(fHeight>0){ // PMT is above water level  
        /// geo 3: both vertex and PMT are in scintillator: find direct and reflected light paths
            std::cout<<"vtpt "<<std::endl;
            scintpath = ScintExWaterPathCalculation(par,pmt);
            if(scintpath == 0) tof = dr/fSpeedOfLightWater;
            else tof = scintpath/fSpeedOfLightScint + (dr-scintpath)/fSpeedOfLightWater;
        }else {  // PMT is below water level
         /// geo 4: vertex is above, PMT is below: find refracted path
            std::cout<<"vtpb "<<std::endl;
            alpha = FactorN(fPosDiff.Perp(), true);
            TVector3 s = alpha*fVertex + (1-alpha)*fPMTpos;  //find a point directly below the intersection with water
            s.SetZ(fWaterLevel);  //move intersection to water surface.
            TVector3 incident = s - fVertex;
            TVector3 refractedRay = fPMTpos - s;
            tof = incident.Mag()/fSpeedOfLightScint + refractedRay.Mag()/fSpeedOfLightWater;  //time of flight for refracted ray
        }
      } //vertex above
       std::cout<<"tof "<<tof<<std::endl;
       //std::cout<<"tof mc "<<tofmc<<std::endl;
      }//valid calculation
   }
  
}
