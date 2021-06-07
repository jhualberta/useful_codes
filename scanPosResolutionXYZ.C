{

  ///!!! calculate for all three axes    
  // x scan x, y scan y, z scan z
  double x[]={-0.209,-4999.043,-4002.525,-3004.229,-2000.155,-992.994,998.133,2011.103,4003.323,5004.868,4503.445,-4489.301,3511.147,2502.805,-3476.709};
  double y[]={-0.201,-998.068,-2000.578,-2998.017,-3992.167,-4999.882,5002.057,3973.021,2980.035,1986.669,994.183,496.858,1494.71,2487.539,3496.338,4505.371,-501.268,-1494.927,-2487.912,-3475.769,-4498.62, -4874.53};
  double z[] ={-5501.2,-4999.899,-4500.2,-4001,-3501.399,-2999.7,-2500.399,-1998.499,-1499.5,-1000.099,-499.6,0.401,500.8,1000.3,1500.9,2000.001,2500.3,3000.3,3500.9,4000.5,4500.4,4973.567};
  int nx = int(sizeof(x)/sizeof(double));
  int ny = int(sizeof(y)/sizeof(double));
  int nz = int(sizeof(z)/sizeof(double));
//=================================================================================================================================
  /// for x-scan: NOTE: not Z, Y, X !!!
  //x-scan X MPW data !!!!
  double X_muP_x[] ={-5.737,-24.15,1.618,-3.955,-8.659,-7.877,-4.001,0.3639,-10.36,10.08,-1.954,-6.386,-9.3,-10.73,-5.199};

  double X_muPerr_x[]={0.9496,2.009,1.848,1.837,1.719,1.875,1.886,1.856,1.823,2.12,1.843,1.276,1.234,1.297,2.869};
  double X_sigmaP_x[] ={175.7,236,159.8,172.3,168.7,169.1,180.9,178.9,169.9,249,176.1,169.5,177.4,180.6,172.5};
  double X_sigmaPerr_x[]={2.197,5.138,4.144,4.236,3.843,4.255,4.368,4.373,4.111,5.253,4.345,2.998,2.831,3.061,6.438};
//  double X_chi2_xscan[]={
  /// x-scan X MC
  double X_MCmuP_x[] ={-7.462,-8.83,-8.549,-2.958,-9.334,-3.848,-3.102,-1.409,-11.13,0.2388,-4.242,-11.78,-7.557,-5.177,-8.648};
  double X_MCmuPerr_x[]={1.046,1.119,1.025,1.021,1.033,1.042,1.045,1.037,1.034,1.13,1.041,1.041,1.039,1.036,1.027};
  double X_MCsigmaP_x[] ={174.3,209.6,162.8,158.6,168.1,173.2,174.5,170.8,167.1,213.7,169,169.7,172.2,170,163.6};
  double X_MCsigmaPerr_x[]={2.345,2.746,2.285,2.308,2.311,2.377,2.394,2.32,2.311,2.765,2.39,2.367,2.318,2.328,2.292};
//  double X_MCchi2_xscan[]={

  //x-scan Y data -------------------------------------------------------------
  double X_muP_y[] ={4.476,4.033,1.388,-1.11,0.07625,4.396,-0.4823,9.188,4.792,1.095,0.2904,4.565,4.681,7.794,2.629};
  double X_muPerr_y[]={0.9483,1.91,1.904,1.839,1.716,1.908,1.867,1.826,1.802,1.958,1.854,1.299,1.21,1.269,2.908};
  double X_sigmaP_y[] ={175,223.3,180.5,173.9,167.2,180.2,173.1,166.8,163,210.6,183.7,184.3,163,165,182.2};
  double X_sigmaPerr_y[]={2.24,4.757,4.52,4.228,4.083,4.368,4.261,4.142,4.135,4.85,4.357,3.113,2.748,2.926,6.549};
//  double X_chi2_yscan[]={
  /// x-scan Y MC
  double X_MCmuP_y[] ={7.097,0.4619,2.751,7.851,8.141,-1.061,0.01685,7.895,1.742,0.4749,1.536,2.454,4.653,6.606,5.904};
  double X_MCmuPerr_y[]={1.044,1.068,1.036,1.039,1.034,1.049,1.051,1.036,1.039,1.087,1.055,1.042,1.04,1.039,1.041};
  double X_MCsigmaP_y[] ={172.3,189,169.4,169.4,168.1,177.5,178,169.6,169.9,197.4,179.9,172.5,172.5,171.4,172.2};
  double X_MCsigmaPerr_y[]={2.337,2.505,2.354,2.302,2.344,2.367,2.419,2.342,2.314,2.547,2.431,2.414,2.304,2.331,2.325};
//  double X_MCchi2_yscan[]={

  //x-scan Z data -----------------------------------------------------------
  double X_muP_z[] ={0.1965,-2.039,-9.236,0.6057,-1.73,-3.347,-3.732,-2.31,-9.367,-8.726,-2.148,-3.95,-4.384,-2.6,-3.025};
  double X_muPerr_z[]={0.9039,1.741,1.788,1.735,1.652,1.805,1.766,1.75,1.715,1.826,1.759,1.235,1.166,1.209,2.71};
  double X_sigmaP_z[] ={149.7,162.5,144.2,141,148.8,150.6,143.1,145.5,137.3,167.1,154.5,157.6,146.1,140,141.3};
  double X_sigmaPerr_z[]={2.105,4.061,4.19,3.896,3.67,4.043,4.039,4.029,4.005,4.471,4.003,2.855,2.719,2.788,6.587};
//  double X_chi2_zscan[]={
  /// x-scan Z MC
  double X_MCmuP_z[] ={-2.048,-2.077,-0.1985,-5.263,-2.266,-1.603,-2.094,-2.021,-8.126,1.151,0.003679,-2.391,-2.971,-6.223,-6.728};
  double X_MCmuPerr_z[]={0.9866,1.005,0.9772,0.9845,0.9745,0.9858,0.9822,0.969,0.9801,1.017,0.9883,0.9841,0.9756,0.9805,0.9802};
  double X_MCsigmaP_z[] ={142.7,154.3,138.3,141.4,136.9,143.5,140.4,133.6,139,158.9,143.3,141.6,137.9,140.7,140};
  double X_MCsigmaPerr_z[]={2.229,2.3,2.156,2.245,2.201,2.205,2.208,2.218,2.174,2.304,2.269,2.214,2.159,2.204,2.191};
//  double X_MCchi2_zscan[]={
//=================================================================================================================================
//
 /// for y-scan
  //y scan Y MPW data !!!!
  double Y_muP_y[] ={5.024,3.524,5.065,8.754,0.4924,-1.255,17.6,0.7364,6.41,-0.5424,-0.9185,5.718,8.597,5,2.997,2.017,8.41,0.5909,2.986,8.991,0.8327,-3.009};
  double Y_muPerr_y[]={1.533,1.496,1.666,1.666,1.674,1.806,1.812,1.644,1.528,1.564,1.691,1.723,1.68,1.706,1.736,1.507,1.311,1.537,1.488,1.559,1.583,1.091};
  double Y_sigmaP_y[] ={172,173.3,166.2,162.4,160.5,204.9,248.9,162.1,168,174,173.9,178.6,182.2,175,173,171.9,175.7,175.3,169,157.8,160.2,197.5};
  double Y_sigmaPerr_y[]={3.571,3.431,3.877,3.897,3.727,4.501,4.768,3.803,3.499,3.576,3.949,3.991,3.896,4.086,3.996,3.619,3.021,3.562,3.371,3.568,3.608,2.78};
//  double Y_chi2_yscan[]={
  /// y-scan Y MC
  double Y_MCmuP_y[] ={0.1569,0.9977,-0.001831,1.676,0.8258,-3.902,10.83,-2.006,6.258,3.329,1.256,-0.5335,7.662,0.5286,6.091,6.89,2.28,-0.7119,7.815,-0.005193,-3.817,-2.118};
  double Y_MCmuPerr_y[]={1.041,1.042,1.038,1.039,1.031,1.125,1.123,1.027,1.022,1.043,1.045,1.051,1.045,1.031,1.031,1.046,1.045,1.039,1.039,1.028,1.043,1.096};
  double Y_MCsigmaP_y[] ={171.3,172.3,171.4,171.7,164.9,215.6,212.5,162,159.8,172.9,174.9,176.3,175.2,166.1,165.1,173.2,174.1,170.4,172.3,163.8,170.9,196.8};
  double Y_MCsigmaPerr_y[]={2.337,2.373,2.358,2.361,2.285,2.721,2.758,2.321,2.356,2.349,2.347,2.339,2.402,2.327,2.329,2.375,2.309,2.347,2.312,2.29,2.377,2.624};
//  double Y_MCchi2_yscan[]={
  //y scan X ---------------------------------------------
  // MPW data !!!!
  double Y_muP_x[] ={-6.478,-8.748,-7.534,-9.984,-11.44,-1.608,-10.38,-6.879,-2.712,-1.804,-9.927,-4.986,-8.453,-3.34,-5.974,-12.64,-7.355,-6.4,-8.557,-9.069,-3.857,-2.26};
  double Y_muPerr_x[]={1.534,1.503,1.683,1.698,1.71,1.769,1.656,1.632,1.514,1.541,1.681,1.696,1.653,1.683,1.718,1.505,1.306,1.536,1.506,1.574,1.628,1.076};
  double Y_sigmaP_x[] ={172.8,177.4,173,174.9,176.6,200.8,202.5,158.3,162.4,164,170.3,168.1,171.4,162.3,165.2,173.7,172.8,175,177.1,164.6,181.9,197.0};
  double Y_sigmaPerr_x[]={3.671,3.389,3.744,3.926,3.924,4.194,4.182,3.67,3.447,3.535,3.832,4.08,3.864,3.862,3.955,3.646,3.037,3.54,3.461,3.607,3.774,2.555};
//  double Y_chi2_xscan[]={
  // MC
  double Y_MCmuP_x[] ={-7.725,-6.422,-8.527,-4.335,-6.146,-12.13,-3.668,-5.866,-6.692,-9.627,-8.841,-1.157,-10.52,-6.76,-5.35,-4.962,-8.93,-7.953,-8.353,-6.315,-5.113,-4.075};
  double Y_MCmuPerr_x[]={1.043,1.045,1.036,1.039,1.034,1.063,1.066,1.042,1.035,1.036,1.037,1.053,1.045,1.038,1.033,1.046,1.04,1.045,1.036,1.032,1.048,1.064};
  double Y_MCsigmaP_x[] ={172.8,175.3,169.9,171.8,168.2,186,187.2,172.9,169,168.8,169.7,177.4,175.2,170.6,166.9,175.3,170.8,174.6,170.4,167.1,176.3,185.0};
  double Y_MCsigmaPerr_x[]={2.355,2.369,2.322,2.346,2.328,2.441,2.495,2.336,2.31,2.321,2.318,2.391,2.438,2.314,2.333,2.377,2.319,2.339,2.315,2.324,2.397,2.424};
//  double Y_MCchi2_xscan[]={

  //y scan Z ---------------------------------------------
  // MPW data !!!!
  double Y_muP_z[] ={-8.613,0.5589,-6.261,-3.452,-6.425,1.538,-0.7923,-2.828,-7.448,-3.448,-1.913,-1.963,-8.723,-4.277,-2.972,-5.419,-7.395,-6.377,0.8009,0.6093,-2.375,-6.172};
  double Y_muPerr_z[]={1.461,1.416,1.592,1.594,1.625,1.675,1.588,1.588,1.45,1.493,1.605,1.639,1.594,1.625,1.652,1.448,1.246,1.461,1.417,1.501,1.552,1.017};
  double Y_sigmaP_z[] ={147,143.7,142.6,139.6,149.4,170.2,180.2,146.3,139.8,149.6,146.1,151.8,153.6,145.7,146.6,154.5,148.1,147.5,142.8,139.3,156.9,165.2};
  double Y_sigmaPerr_z[]={3.416,3.292,3.663,3.808,3.667,3.893,3.834,3.704,3.347,3.513,3.79,3.743,3.602,3.809,3.815,3.402,2.84,3.396,3.216,3.5,3.636,2.405};
//  double Y_chi2_zscan[]={
  // MC
  double Y_MCmuP_z[] ={-1.005,-1.019,-1.886,-5.995,-8.836,-3.362,-2.199,-1.633,-5.254,-5.508,-1.109,0.2486,-1.69,-2.47,-8.323,1.126,-0.8448,-1.851,-1.956,-6.394,-3.192,-4.264};
  double Y_MCmuPerr_z[]={0.987,0.9855,0.9781,0.9793,0.9786,1.002,0.9964,0.9725,0.9757,0.9853,0.9851,0.9874,0.983,0.9753,0.9721,0.9852,0.9833,0.9845,0.973,0.9783,0.9892,0.9958};
  double Y_MCsigmaP_z[] ={143.2,143.2,139.3,139.2,138.5,153.1,148.1,134.4,136.6,142.3,143,142.2,142,137,133.9,142.5,141.2,142.1,136.3,138.7,144.8,147.2};
  double Y_MCsigmaPerr_z[]={2.219,2.24,2.194,2.18,2.195,2.259,2.3,2.172,2.169,2.177,2.207,2.178,2.23,2.213,2.193,2.211,2.206,2.229,2.135,2.16,2.226,2.298};
//  double Y_MCchi2_zscan[]={

//=================================================================================================================================
//
 /// for z-scan
  //z scan z ---------------------------------------------
  // MPW data !!!!
  //
  double Z_muP_z[] ={-2.516,-10.07,-8.424,-4.118,-0.4166,1.828,-2.9,2.367,-4.469,-1.271,-6.881,-1.412,-3.52,-1.768,-6.612,-2.343,-1.167,-8.968,-3.52,-3.503,-2.256,11.73};
  double Z_muPerr_z[]={2.089,1.931,2.258,1.442,1.532,1.532,1.593,1.599,1.651,1.681,1.682,1.729,1.575,1.738,1.306,1.94,1.795,1.814,1.767,1.869,2.046,2.449};
  double Z_sigmaP_z[] ={208,193.4,144.4,152.4,152.7,144.1,151.2,147.2,158.9,160.6,147.5,147.9,149,144.3,146.8,140.5,137.7,143.5,139.4,141,158.3,218.2};
  double Z_sigmaPerr_z[]={4.798,4.89,5.422,3.313,3.58,3.662,3.735,3.864,3.79,3.888,3.886,4.003,3.692,4.142,2.983,4.542,4.143,4.225,3.901,4.169,4.647,6.197};
//  double Z_chi2_zscan[]={
 /// z-scan z MC
  double Z_MCmuP_z[] ={-0.5407,-4.36,-6.709,-1.22,-6.656,-0.06913,-6.219,0.8981,-2.705,-5.081,-8.99,-6.273,-0.8623,-6.937,-1.377,-7.334,-3.509,-10.34,-8.813,-4.978,-3.691};
  double Z_MCmuPerr_z[]={1.068,1.028,0.9756,0.9671,0.9656,0.9652,0.9718,0.9736,0.9819,0.9889,0.9901,0.9858,0.9875,0.977,0.9811,0.979,0.9848,0.9832,0.9908,1.009,1.091};
  double Z_MCsigmaP_z[] ={170.1,167.1,139,134.9,133.1,133.8,137.8,137.6,142.1,145.2,145.1,142.2,142.5,136.3,139.2,139.1,139.5,139,142,145,178.9};
  double Z_MCsigmaPerr_z[]={2.304,2.414,2.208,2.166,2.153,2.15,2.134,2.183,2.182,2.173,2.211,2.183,2.184,2.2,2.218,2.215,2.184,2.188,2.171,2.248,2.518};
//  double Z_MCchi2_zscan[]={

  //z scan X ---------------------------------------------
  // MPW data !!!!

  double Z_muP_x[] ={-4.976,-4.219,-4.551,-3.438,-5.512,-6.829,-4.814,-4.069,-4.693,-4.604,-3.7,-1.843,-5.21,-3.499,-3.453,-10.02,-11.15,-11.13,-2.539,-6.864,-5.276,-13.37};
  double Z_muPerr_x[]={1.958,1.85,2.347,1.476,1.566,1.578,1.628,1.658,1.681,1.718,1.763,1.82,1.649,1.829,1.385,2.067,1.91,1.934,1.887,2.028,2.221,2.616};
  double Z_sigmaP_x[] ={186.9,162.7,165.9,161.9,160.4,157.4,159.1,164.9,163.9,168.4,172.5,176.7,172.8,172,179.8,176.5,171.6,180,176.5,191.3,214.2,293.6};
  double Z_sigmaPerr_x[]={4.346,4.257,5.539,3.297,3.554,3.532,3.829,3.872,3.847,3.891,4.033,4.289,3.809,4.345,3.212,5.006,4.402,4.454,4.437,4.796,5.515,7.856};
  // double Z_chi2_xscan[]={
  // MC, z scan X
  double Z_MCmuP_x[] ={-6.45,-6.481,-6.202,-5.676,-5.14,-4.748,-4,-4.821,-3.613,-4.269,-3.764,-5.796,-6.615,-4.353,-4.329,-6.786,-8.143,-5.472,-5.722,-4.043,-6.001,-6.495};
  double Z_MCmuPerr_x[]={1.065,1.047,1.035,1.024,1.029,1.022,1.028,1.033,1.041,1.036,1.041,1.05,1.041,1.052,1.044,1.043,1.036,1.049,1.047,1.057,1.09,1.147};
  double Z_MCsigmaP_x[] ={189.7,179.3,171.7,165,167.3,163.6,167.4,169.3,173.6,170.2,172.1,177.0,170.8,177.6,173.1,171.9,169.1,173.9,173.6,177.6,190.6,216.1};
  double Z_MCsigmaPerr_x[]={2.378,2.39,2.272,2.329,2.302,2.278,2.338,2.322,2.334,2.325,2.331,2.375,2.392,2.376,2.373,2.348,2.334,2.373,2.377,2.396,2.531,2.773};
//  double Z_MCchi2_xscan[]={

  //z scan Y ---------------------------------------------
  // MPW data !!!!

  double Z_muP_y[] ={4.725,0.007045,2.348,2.97,4,2.769,6.574,2.489,-0.6907,1.361,-3.314,5.901,2.351,-0.1499,0.5201,2.027,7.413,7.26,4.269,1.461,1.208,14.79};
  double Z_muPerr_y[]={1.945,1.881,2.341,1.483,1.574,1.582,1.641,1.675,1.695,1.73,1.765,1.835,1.667,1.85,1.382,2.074,1.93,1.941,1.899,2.06,2.221,2.658};
  double Z_sigmaP_y[] ={183,174.4,163.2,164.5,163.7,158.6,164.4,172.2,169.1,173,172.9,181.9,180.7,179.9,177.8,178.1,179.8,182.1,181.1,202.6,213.5,308.7};
  double Z_sigmaPerr_y[]={4.543,4.468,5.353,3.393,3.554,3.675,3.733,3.799,3.856,3.959,4.151,4.302,3.929,4.255,3.209,4.75,4.462,4.598,4.503,4.879,5.735,7.921};
//  double Z_chi2_yscan[]={
  // MC, z scan Y
  double Z_MCmuP_y[] ={-1.314,8.723,-1.138,-1.003,0.8264,0.3692,8.259,-1.452,-0.83,-1.403,-1.967,-0.1117,0.1763,-0.2159,8.747,-1.911,0.5352,-1.374,-1.073,8.365,-1.141,2.376};
  double Z_MCmuPerr_y[]={1.062,1.046,1.035,1.024,1.025,1.021,1.033,1.035,1.041,1.035,1.041,1.046,1.04,1.048,1.053,1.049,1.04,1.048,1.048,1.053,1.084,1.153};
  double Z_MCsigmaP_y[] ={187.8,178.6,171,164.8,164.2,163.1,170.7,170.5,172.6,169.4,171.8,174,170.2,174.8,178.5,176.2,171.6,173.2,173.7,175.1,186.8,220};
  double Z_MCsigmaPerr_y[]={2.391,2.371,2.367,2.302,2.336,2.28,2.316,2.315,2.325,2.343,2.357,2.365,2.342,2.377,2.388,2.344,2.365,2.338,2.357,2.384,2.485,2.713};
//  double Z_MCchi2_yscan[]={

  double X_ex[nx],Y_ex[ny],Z_ex[nz];
  for(int i = 0;i<nx;i++) X_ex[i] = 0;
  for(int i = 0;i<ny;i++) Y_ex[i] = 0;
  for(int i = 0;i<nz;i++) Z_ex[i] = 0;

  //////////////////////////////////////////////////////////////////////////
  //vertex shifts
  double Xshift[nx];
  double Yshift[ny];
  double Zshift[nz];

  double XshiftErr[nx];
  double YshiftErr[ny];
  double ZshiftErr[nz];

//  vector<double> Xshift;
//  vector<double> Yshift;
//  vector<double> Zshift;
//
//  vector<double> XshiftErr;
//  vector<double> YshiftErr;
//  vector<double> ZshiftErr;

  TH1F *hXshift = new TH1F("hXshift","",400,-20,20);
  TH1F *hYshift = new TH1F("hYshift","",400,-20,20);
  TH1F *hZshift = new TH1F("hZshift","",400,-20,20);

  cout<<"\n x----"<<endl;
  for(int i = 0;i<nx;i++)
  {
    double deltaX = X_muP_x[i] - X_MCmuP_x[i];
    double deltaXerr =  X_muPerr_x[i] +  X_MCmuPerr_x[i]; 
    Xshift[i] = deltaX;
    XshiftErr[i] = deltaXerr;
    hXshift->Fill(deltaX);
    cout<<deltaX<<", "<<deltaXerr<<" ";
    //Xshift.push_back(deltaX);
    //XshiftErr.push_back(deltaXerr);
    //double deltaY = X_muP_y[i] - X_MCmuP_y[i];
    //double deltaYerr = X_muPerr_y[i] + X_MCmuPerr_y[i];
    //Yshift.push_back(deltaY);
    //YshiftErr.push_back(deltaYerr);
    //double deltaZ = X_muP_z[i] - X_MCmuP_z[i];
    //double deltaZerr =  X_muPerr_z[i] + X_MCmuPerr_z[i];
    //Zshift.push_back(deltaZ);
    //ZshiftErr.push_back(deltaZerr);
  }

  cout<<"\n y----"<<endl;
  for(int i = 0;i<ny;i++)
  { 
    double deltaY = Y_muP_y[i] - Y_MCmuP_y[i];
    double deltaYerr =  Y_muPerr_y[i] + Y_MCmuPerr_y[i];
    Yshift[i] = deltaY;
    YshiftErr[i] = deltaYerr;
    hYshift->Fill(deltaY);
    cout<<deltaY<<","<<deltaYerr<<" ";
//    double deltaY = Y_muP_y[i] - Y_MCmuP_y[i];
//    double deltaYerr =  Y_muPerr_y[i] + Y_MCmuPerr_y[i];
//    Yshift.push_back(deltaY);
//    YshiftErr.push_back(deltaYerr);
//    double deltaZ = Y_muP_z[i] - Y_MCmuP_z[i];
//    double deltaZerr = Y_muPerr_z[i] + Y_MCmuPerr_z[i];
//    Zshift.push_back(deltaZ);
//    ZshiftErr.push_back(deltaZerr);
   }
//
  cout<<"\n z----"<<endl;
  for(int i = 0;i<nz;i++)
  { 
    double deltaZ = Z_muP_z[i] - Z_MCmuP_z[i];
    double deltaZerr =  Z_muPerr_z[i] +  Z_MCmuPerr_z[i];
    Zshift[i] = deltaZ;
    ZshiftErr[i] = deltaZerr;
    hZshift->Fill(deltaZ);
    cout<<deltaZ<<", "<<deltaZerr<<" ";
//    double deltaY = Z_muP_y[i] - Z_MCmuP_y[i];
//    double deltaZerr =  Z_muPerr_y[i] +  Z_MCmuPerr_y[i];
//    Yshift.push_back(deltaY);
//    YshiftErr.push_back(deltaYerr);
//    double deltaZ = Z_muP_z[i] - Z_MCmuP_z[i];
//    double deltaZerr =  Z_muPerr_z[i] +  Z_MCmuPerr_z[i];
//    Zshift.push_back(deltaZ);
//    ZshiftErr.push_back(deltaZerr);
  }
  cout<<endl;

//  double combineX[55] = {-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-186,-185.037,-5.995,-7.761,-7.084,-5.491,-3.774,-1.624,-1.745,-3.967,-5.984,-7.952,-9.242,-9.867,-8.414,-6.835,-4.949,-2.539,-7.711,-7.534,-6.349,-4.366,-2.799,-5.283,-4999.043,-4002.525,-3004.229,-2000.155,-992.994,998.133,2011.103,4003.323,5004.868,4503.445,-4489.301,3511.147,2502.805,-3476.709};
//
//  double combineY[55] = {254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,254,247.24,-0.201,-998.068,-2000.578,-2998.017,-3992.167,-4999.882,5002.057,3973.021,2980.035,1986.669,994.183,496.858,1494.71,2487.539,3496.338,4505.371,-501.268,-1494.927,-2487.912,-3475.769,-4498.62,-0.209,2.46,5.269,8.101,10.637,11.641,10.897,9.874,4.378,2.262,3.278,3.918,5.742,8.681,6.643};
//  double combineZ[55] = {};

//  for(int i = 0;i<ny;i++)
//  {
//    double deltaX = Y_muP_x[i] - Y_MCmuP_x[i];
//    double deltaXerr =  Y_muPerr_x[i] +  Y_MCmuPerr_x[i];
//    Yshift[i] = deltaX;
//    hYshift->Fill(deltaX);
//    YshiftErr[i] = deltaXerr;
//    cout<<deltaX<<",";
//  }
//  cout<<"\n z----"<<endl;
//  for(int i = 0;i<nz;i++)
//  {
//    double deltaX = Z_muP_x[i] - Z_MCmuP_x[i];
//    double deltaXerr =  Z_muPerr_x[i] +  Z_MCmuPerr_x[i];
//    Zshift[i] = deltaX;
//    hZshift->Fill(deltaX);
//    ZshiftErr[i] = deltaXerr;
//    cout<<deltaX<<",";
//  }

  /////////////////////////////////////////////////////////////////////////
  //Resolution shift: calculate x, y, z respectively
  double sum_sigmaP2_scanX[nx];
  double sum_sigmaP2_scanY[ny];
  double sum_sigmaP2_scanZ[nz];
  for(int i = 0;i<nx;i++)
  {
    double delta2X = X_sigmaP_x[i]*X_sigmaP_x[i] - X_MCsigmaP_x[i]*X_MCsigmaP_x[i];
    double delta2Y = X_sigmaP_y[i]*X_sigmaP_y[i] - X_MCsigmaP_y[i]*X_MCsigmaP_y[i];
    double delta2Z = X_sigmaP_z[i]*X_sigmaP_z[i] - X_MCsigmaP_z[i]*X_MCsigmaP_z[i];
    ///sum_sigmaP2_scanX[i]=delta2X+delta2Y+delta2Z;//sqrt(fabs(delta2X+delta2Y+delta2Z));
    sum_sigmaP2_scanX[i] =sqrt(fabs(delta2X+delta2Y+delta2Z));
    // cout<<sum_sigmaP2_scanX[i]<<",";
  }
  for(int i = 0;i<ny;i++)
  { 
    double delta2X = Y_sigmaP_x[i]*Y_sigmaP_x[i] - Y_MCsigmaP_x[i]*Y_MCsigmaP_x[i];
    double delta2Y = Y_sigmaP_y[i]*Y_sigmaP_y[i] - Y_MCsigmaP_y[i]*Y_MCsigmaP_y[i];
    double delta2Z = Y_sigmaP_z[i]*Y_sigmaP_z[i] - Y_MCsigmaP_z[i]*Y_MCsigmaP_z[i];
    //sum_sigmaP2_scanY[i]=delta2X+delta2Y+delta2Z;//(sqrt(fabs(delta2X+delta2Y+delta2Z)));
    sum_sigmaP2_scanY[i] =sqrt(fabs(delta2X+delta2Y+delta2Z));
    // cout<<sum_sigmaP2_scanY[i]<<",";
  }
  for(int i = 0;i<nz;i++)
  { 
    double delta2X = Z_sigmaP_x[i]*Z_sigmaP_x[i] - Z_MCsigmaP_x[i]*Z_MCsigmaP_x[i];
    double delta2Y = Z_sigmaP_y[i]*Z_sigmaP_y[i] - Z_MCsigmaP_y[i]*Z_MCsigmaP_y[i];
    double delta2Z = Z_sigmaP_z[i]*Z_sigmaP_z[i] - Z_MCsigmaP_z[i]*Z_MCsigmaP_z[i];
    //sum_sigmaP2_scanZ[i]=delta2X+delta2Y+delta2Z;//
    sum_sigmaP2_scanZ[i] =sqrt(fabs(delta2X+delta2Y+delta2Z));
    // cout<<sum_sigmaP2_scanZ[i]<<",";
  }
  TCanvas *c1 = new TCanvas("c1","resolution",800,600);
  c1->cd();
  TGraph *gx_bias = new TGraph(nx,x,sum_sigmaP2_scanX);
  gx_bias->GetXaxis()->SetTitle("source position [mm]");
  //gx_bias->GetYaxis()->SetTitle("Data^{2}-MC^{2}");
  gx_bias->SetMarkerStyle(25);
  gx_bias->SetMarkerColor(kRed);
  gx_bias->SetMarkerSize(2.5);
  gx_bias->SetTitle("X scan");

  TGraph *gy_bias = new TGraph(ny,y,sum_sigmaP2_scanY);
  gy_bias->GetXaxis()->SetTitle("mm");
  //gy_bias->GetYaxis()->SetTitle("Data^{2}-MC^{2}");
  gy_bias->SetMarkerStyle(24);
  gy_bias->SetMarkerColor(kGreen+2);
  gy_bias->SetMarkerSize(2.5);
  gy_bias->SetTitle("Y scan");

  TGraph *gz_bias = new TGraph(nz,z,sum_sigmaP2_scanZ);
  gz_bias->GetXaxis()->SetTitle("mm");
  //gz_bias->GetYaxis()->SetTitle("Data^{2}-MC^{2}");
  gz_bias->SetMarkerStyle(22);
  gz_bias->SetMarkerColor(kBlue);
  gz_bias->SetMarkerSize(2.5);
  gz_bias->SetTitle("Z scan");

  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(gx_bias,"along X-axis","p");
  legend->AddEntry(gy_bias,"along Y-axis","p");
  legend->AddEntry(gz_bias,"along Z-axis","p");

  gx_bias->Draw("AP");
  gy_bias->Draw("P");
  gz_bias->Draw("P");
  legend->Draw("same");
//

  TCanvas *c2 = new TCanvas("c2","resolution",800,600);
  c2->cd();
  TGraphErrors *gx_shift = new TGraphErrors(nx,x,Xshift, X_ex,XshiftErr);
  gx_shift->GetXaxis()->SetTitle("source position [mm]");
  gx_shift->GetYaxis()->SetTitle("#mu_{Data}-#mu_{MC} [mm]");
  gx_shift->SetMarkerStyle(25);
  gx_shift->SetMarkerColor(kRed);
  gx_shift->SetMarkerSize(2.5);
  gx_shift->SetTitle("X scan");

  TGraphErrors *gy_shift = new TGraphErrors(ny,y,Yshift,Y_ex,YshiftErr);
//  gy_shift->GetXaxis()->SetTitle("source position [mm]");
//  gy_shift->GetYaxis()->SetTitle("Data-MC [mm]");
  gy_shift->SetMarkerStyle(24);
  gy_shift->SetMarkerColor(kGreen+2);
  gy_shift->SetMarkerSize(2.5);
  gy_shift->SetTitle("Y scan");

  TGraphErrors *gz_shift = new TGraphErrors(nz,z,Zshift,Z_ex,ZshiftErr);
//  gz_shift->GetXaxis()->SetTitle("source position [mm]");
//  gz_shift->GetYaxis()->SetTitle("Data-MC [mm]");
  gz_shift->SetMarkerStyle(22);
  gz_shift->SetMarkerColor(kBlue);
  gz_shift->SetMarkerSize(2.5);
  gz_shift->SetTitle("Z scan");

  TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9);
  legend1->AddEntry(gx_bias,"along X-axis","p");
  legend1->AddEntry(gy_bias,"along Y-axis","p");
  legend1->AddEntry(gz_bias,"along Z-axis","p");

  TCanvas *c001 = new TCanvas("c001","shifts",900,600);
  gx_shift->Draw("AP");
//  gy_shift->Draw("P");
//  gz_shift->Draw("P");
//  TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9);
//  legend1->AddEntry(gx_bias,"along X-axis","p");
//  legend1->AddEntry(gy_bias,"along Y-axis","p");
//  legend1->AddEntry(gz_bias,"along Z-axis","p");
//  legend1->Draw("same");


  TCanvas *c00 = new TCanvas("c00", "muP", 900, 600); 
  c00->Divide(3,2);
  c00->cd(1);
  TGraphErrors *gx_muP = new TGraphErrors(nx,x,X_muP_x,X_ex,X_muPerr_x); //z scan, z
  TGraphErrors *gxMC_muP = new TGraphErrors(nx,x,X_MCmuP_x,X_ex,X_MCmuPerr_x);
  gx_muP->SetMarkerStyle(21);
  gxMC_muP->SetMarkerStyle(24);
  gxMC_muP->SetMarkerColor(kRed);
  gx_muP->SetMarkerSize(2);
  gxMC_muP->SetMarkerSize(2);
  gx_muP->Draw("AP");
  gxMC_muP->Draw("P");
  gx_muP->SetTitle("X scan");
  gx_muP->GetYaxis()->SetRangeUser(-30,30);
  gx_muP->GetYaxis()->SetTitle("#mu_{P}");

  c00->cd(2);
  TGraphErrors *gy_muP = new TGraphErrors(ny,y,Y_muP_y,Y_ex,Y_muPerr_y); //z scan, z
  TGraphErrors *gyMC_muP = new TGraphErrors(ny,y,Y_MCmuP_y,Y_ex,Y_MCmuPerr_y);
  gy_muP->SetMarkerStyle(21);
  gyMC_muP->SetMarkerStyle(24);
  gyMC_muP->SetMarkerColor(kRed);
  gy_muP->SetMarkerSize(2);
  gyMC_muP->SetMarkerSize(2);
  gy_muP->Draw("AP");
  gyMC_muP->Draw("P");
  gy_muP->SetTitle("Y scan");
  gy_muP->GetYaxis()->SetRangeUser(-30,30);
  gy_muP->GetYaxis()->SetTitle("#mu_{P}");

  c00->cd(3);
  TGraphErrors *gz_muP = new TGraphErrors(nz,z,Z_muP_z,Z_ex,Z_muPerr_z); //z scan, z
  TGraphErrors *gzMC_muP = new TGraphErrors(nz,z,Z_MCmuP_z,Z_ex,Z_MCmuPerr_z); 
  gz_muP->SetMarkerStyle(21);
  gzMC_muP->SetMarkerStyle(24);
  gzMC_muP->SetMarkerColor(kRed);
  gz_muP->SetMarkerSize(2);
  gzMC_muP->SetMarkerSize(2);
  gz_muP->Draw("AP");
  gzMC_muP->Draw("P");
  gz_muP->SetTitle("Z scan");
  gz_muP->GetYaxis()->SetRangeUser(-30,30);
  gz_muP->GetYaxis()->SetTitle("#mu_{P}");

  c00->cd(4);
  TGraphErrors *gx_sigmaP = new TGraphErrors(nx,x,X_sigmaP_x,X_ex,X_sigmaPerr_x); //z scan, z
  TGraphErrors *gxMC_sigmaP = new TGraphErrors(nx,x,X_MCsigmaP_x,X_ex,X_MCsigmaPerr_x);
  gx_sigmaP->SetMarkerStyle(21);
  gxMC_sigmaP->SetMarkerStyle(24);
  gxMC_sigmaP->SetMarkerColor(kRed);
  gx_sigmaP->SetMarkerSize(2);
  gxMC_sigmaP->SetMarkerSize(2);
  gx_sigmaP->Draw("AP");
  gxMC_sigmaP->Draw("P");
  gx_sigmaP->GetXaxis()->SetTitle("X scan, source position [mm]");
  gx_sigmaP->GetYaxis()->SetRangeUser(100,330);
  gx_sigmaP->GetYaxis()->SetTitle("#sigma_{P}");

  c00->cd(5);
  TGraphErrors *gy_sigmaP = new TGraphErrors(ny,y,Y_sigmaP_y,Y_ex,Y_sigmaPerr_y); //z scan, z
  TGraphErrors *gyMC_sigmaP = new TGraphErrors(ny,y,Y_MCsigmaP_y,Y_ex,Y_MCsigmaPerr_y);
  gy_sigmaP->SetMarkerStyle(21);
  gyMC_sigmaP->SetMarkerStyle(24);
  gyMC_sigmaP->SetMarkerColor(kRed);
  gy_sigmaP->SetMarkerSize(2);
  gyMC_sigmaP->SetMarkerSize(2);
  gy_sigmaP->Draw("AP");
  gyMC_sigmaP->Draw("P");
  gy_sigmaP->GetXaxis()->SetTitle("Y scan, source position [mm]");
  gy_sigmaP->GetYaxis()->SetRangeUser(100,330);
  gy_sigmaP->GetYaxis()->SetTitle("#sigma_{P}");

  c00->cd(6);
  TGraphErrors *gz_sigmaP = new TGraphErrors(nz,z,Z_sigmaP_z,Z_ex,Z_sigmaPerr_z); //z scan, z
  TGraphErrors *gzMC_sigmaP = new TGraphErrors(nz,z,Z_MCsigmaP_z,Z_ex,Z_MCsigmaPerr_z);
  gz_sigmaP->SetMarkerStyle(21);
  gzMC_sigmaP->SetMarkerStyle(24);
  gzMC_sigmaP->SetMarkerColor(kRed);
  gz_sigmaP->SetMarkerSize(2);
  gzMC_sigmaP->SetMarkerSize(2);
  gz_sigmaP->Draw("AP");
  gzMC_sigmaP->Draw("P");
  gz_sigmaP->GetXaxis()->SetTitle("Z scan, source position [mm]");
  gz_sigmaP->GetYaxis()->SetRangeUser(100,330);
  gz_sigmaP->GetYaxis()->SetTitle("#sigma_{P}");


  TCanvas *c3 = new TCanvas("c3","resolution",800,600);
  hXshift->SetLineColor(kRed);
  hYshift->SetLineColor(kGreen+2);
  hZshift->SetLineColor(kBlue);
  hXshift->Draw();//hYshift->Draw("sames");hZshift->Draw("sames");
  hYshift->Draw();
  hZshift->Draw();


}
