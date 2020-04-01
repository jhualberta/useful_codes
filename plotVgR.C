{
 double vg[]={1.35, 1.36,1.37,1.38,1.383,1.384,1.385,1.39,1.40,1.41,1.42,1.43,1.45};
 //before correction
 //x,xrms,y,yrms,z,zrms,pos,posrms
//Note R = |Rfit| - |Rmc|
 double fitInfo_1p35[]={5.64077, 321.524, 3.22485, 332.926, 0.904442, 317.019, 163.345, 306.324};
 double fitInfo_1p36[]={5.94111, 307.432, 2.06887, 319.65, 0.265798, 306.63, 128.606, 299.137};
 double fitInfo_1p37[]={6.16812, 305.299, 2.3687, 311.701, 0.119632, 305.081, 93.6923, 304.952};
 double fitInfo_1p38[]={6.95801, 297.273, 1.8761, 311.325, 0.294108, 294.332, 61.4275, 295.041};
 double fitInfo_1p383[]={6.25337, 297.642, 2.63445, 308.14, 0.242115, 298.541, 51.7188, 300.989};
 double fitInfo_1p384[]={5.81795, 296.399, 2.47519, 307.007, 0.713736, 297.523, 47.6318, 299.427};
 double fitInfo_1p385[]={6.88338, 296.07, 1.85846, 309.301, 0.772469, 294.033, 45.3939, 293.415};
 double fitInfo_1p39[]={6.3869, 295.319, 2.32596, 304.765, 0.199697, 293.438, 29.0415, 297.357};
 double fitInfo_1p40[]={5.39453, 295.501, 2.57475, 303.568, 1.78653, 296.226, -3.90418, 297.666};
 double fitInfo_1p41[]={6.55162, 294.214, 2.27693, 303.855, 0.965511, 296.371, -34.9067, 296.489};
 double fitInfo_1p42[]={6.95274, 292.343, 2.22646, 307.016, 1.91301, 292.589, -66.0278, 291.049};
 double fitInfo_1p43[]={7.09668, 295.422, 1.61536, 306.074, 1.75428, 297.431, -97.0888, 291.024};
 double fitInfo_1p45[]={7.82997, 301.138, 1.69691, 312.083, 0.528011, 304.146, -156.047, 294.49};

//Note R = (Rfit-Rmc)*Rmc.Unit()
 double fitInfo1_1p35[]={5.64077, 321.524, 3.22485, 332.926, 0.904442, 317.019, 140.308,318.853};
 double fitInfo1_1p36[]={5.94111, 307.432, 2.06887, 319.65, 0.265798, 306.63, 107.295,309.491};
 double fitInfo1_1p37[]={6.16812, 305.299, 2.3687, 311.701, 0.119632, 305.081, 72.845,315.388};
 double fitInfo1_1p38[]={6.95801, 297.273, 1.8761, 311.325, 0.294108, 294.332, 39.5382,310.801};
 double fitInfo1_1p383[]={6.25337, 297.642, 2.63445, 308.14, 0.242115, 298.541, 30.2119,311.913};
 double fitInfo1_1p384[]={5.81795, 296.399, 2.47519, 307.007, 0.713736, 297.523, 26.1897,310.717};
 double fitInfo1_1p385[]={6.88338, 296.07, 1.85846, 309.301, 0.772469, 294.033, 23.8328,305.759};
 double fitInfo1_1p39[]={6.3869, 295.319, 2.32596, 304.765, 0.199697, 293.438, 7.8707,308.309};
 double fitInfo1_1p40[]={5.39453, 295.501, 2.57475, 303.568, 1.78653, 296.226, -25.1402,309.067};
 double fitInfo1_1p41[]={6.55162, 294.214, 2.27693, 303.855, 0.965511, 296.371, -56.1317,306.774};
 double fitInfo1_1p42[]={6.95274, 292.343, 2.22646, 307.016, 1.91301, 292.589, -87.4152,300.61};
 double fitInfo1_1p43[]={7.09668, 295.422, 1.61536, 306.074, 1.75428, 297.431, -117.949,302.322};
 double fitInfo1_1p45[]={7.82997, 301.138, 1.69691, 312.083, 0.528011, 304.146, -176.717,302.765};

//after correction
 double fitInfoCor_1p35[]=={5.84304, 321.492, 2.68501, 333.43, 1.6432, 317.227, 166.16, 307.587};
 double fitInfoCor_1p36[]={6.16079, 307.109, 1.51769, 319.82, 0.996309, 306.678, 131.199, 299.476};
 double fitInfoCor_1p37[]={6.37953, 304.884, 2.4457, 317.588, 0.862657, 304.994, 96.1, 304.596};
 double fitInfoCor_1p38[]={7.17575, 296.751, 1.34298, 311.418, 1.03443, 293.912, 63.6663, 293.903};
 double fitInfoCor_1p383[]={6.4585, 297.149, 2.08764, 308.132, 0.975201, 298.305, 53.8823, 299.826};
 double fitInfoCor_1p384[]={6.02063, 295.875, 1.94926, 307.048, 1.44134, 297.229, 49.7792, 298.111};
 double fitInfoCor_1p385[]={7.08902, 295.484, 1.34693, 309.297, 1.50083, 293.635, 47.5232, 291.922};
 double fitInfoCor_1p39[]={6.57957, 294.748, 1.80261, 304.716, 0.930854, 292.914, 31.0645, 295.671};
 double fitInfoCor_1p40[]={5.57565, 294.968, 2.05852, 303.492, 2.52157, 295.877, -2.02632, 295.415};
 double fitInfoCor_1p41[]={6.76549, 293.594, 1.75452, 303.695, 1.70365, 295.915, -33.1962, 293.617};
 double fitInfoCor_1p42[]={7.17846, 291.647, 1.69589, 307.029, 2.62399, 292.102, -64.4806, 287.527};
 double fitInfoCor_1p43[]={7.30847, 294.705, 1.06474, 305.909, 2.49376, 296.971, -95.6882, 286.831};
 double fitInfoCor_1p45[]={8.05342, 300.519, 1.16475, 311.915, 1.28182, 303.519, -154.935, 289.287};

 double fitInfoCor1_1p35[]=={6.19109, 307.316, 2.37869, 319.3, 1.99949, 303.388, 90.1678, 300.954};
 double fitInfoCor1_1p36[]={6.50096, 295.674, 1.22812, 308.256, 1.35798, 295.664, 55.7693, 293.507};
 double fitInfoCor1_1p37[]={6.71268, 296.456, 2.14134, 309.219, 1.22664, 297.038, 21.2338, 299.175};
 double fitInfoCor1_1p38[]={7.49636, 291.164, 1.04928, 305.633, 1.39565, 288.845, -10.675, 289.086};
 double fitInfoCor1_1p383[]={6.78718, 292.465, 1.78891, 303.332, 1.3372, 294.053, -20.3039, 295.072};
 double fitInfoCor1_1p384[]={6.35637, 291.505, 1.64396, 302.7, 1.79571, 293.359, -24.3407, 293.454};
 double fitInfoCor1_1p385[]={7.41085, 291.405, 1.06077, 304.886, 1.85429, 289.964, -26.5583, 287.355};
 double fitInfoCor1_1p39[]={6.90609, 292.052, 1.50888, 301.861, 1.29354, 290.81, -42.7543, 291.332};
 double fitInfoCor1_1p40[]={5.91825, 295.173, 1.76076, 303.524, 2.85856, 296.565, -75.3113, 291.612};
 double fitInfoCor1_1p41[]={7.08934, 296.533, 1.46158, 306.354, 2.05393, 299.245, -105.979, 290.33};
 double fitInfoCor1_1p42[]={7.49581, 297.298, 1.40377, 312.189, 2.95233, 298.288, -136.759, 284.872};
 double fitInfoCor1_1p43[]={7.62349, 302.828, 0.782514, 313.677, 2.83126, 305.53, -167.464, 284.685};
 double fitInfoCor1_1p45[]={8.35577, 313.123, 0.89035, 324.107, 1.63919, 316.538, -225.755, 288.005};

//R = (Rfit-Rmc)*Rmc.Unit(), p0=1, p1
 double fitInfo1Cor_1p35[]={5.64077, 321.524, 3.22485, 332.926, 0.904442, 317.019, 140.308,318.853};
 double fitInfo1Cor_1p36[]={5.94111, 307.432, 2.06887, 319.65, 0.265798, 306.63, 107.295,309.491};
 double fitInfo1Cor_1p37[]={6.16812, 305.299, 2.3687, 311.701, 0.119632, 305.081, 72.845,315.388};
 double fitInfo1Cor_1p38[]={6.95801, 297.273, 1.8761, 311.325, 0.294108, 294.332, 39.5382,310.801};
 double fitInfo1Cor_1p383[]={6.25337, 297.642, 2.63445, 308.14, 0.242115, 298.541, 30.2119,311.913};
 double fitInfo1Cor_1p384[]={5.81795, 296.399, 2.47519, 307.007, 0.713736, 297.523, 26.1897,310.717};
 double fitInfo1Cor_1p385[]={6.88338, 296.07, 1.85846, 309.301, 0.772469, 294.033, 23.8328,305.759};
 double fitInfo1Cor_1p39[]={6.3869, 295.319, 2.32596, 304.765, 0.199697, 293.438, 7.8707,308.309};
 double fitInfo1Cor_1p40[]={5.39453, 295.501, 2.57475, 303.568, 1.78653, 296.226, -25.1402,309.067};
 double fitInfo1Cor_1p41[]={6.55162, 294.214, 2.27693, 303.855, 0.965511, 296.371, -56.1317,306.774};
 double fitInfo1Cor_1p42[]={6.95274, 292.343, 2.22646, 307.016, 1.91301, 292.589, -87.4152,300.61};
 double fitInfo1Cor_1p43[]={7.09668, 295.422, 1.61536, 306.074, 1.75428, 297.431, -117.949,302.322};
 double fitInfo1Cor_1p45[]={7.82997, 301.138, 1.69691, 312.083, 0.528011, 304.146, -176.717,302.765};

 double PosMag[]={fitInfo_1p35[6],fitInfo_1p36[6],fitInfo_1p37[6],fitInfo_1p38[6],fitInfo_1p383[6],fitInfoCor_1p384[6],fitInfoCor_1p385[6],fitInfo_1p39[6],fitInfo_1p40[6],fitInfo_1p41[6],fitInfo_1p42[6],fitInfo_1p43[6],fitInfo_1p45[6]};
 double PosErr[]={fitInfo_1p35[7],fitInfo_1p36[7],fitInfo_1p37[7],fitInfo_1p38[7],fitInfo_1p383[7],fitInfoCor_1p384[7],fitInfoCor_1p385[7],fitInfo_1p39[7],fitInfo_1p40[7],fitInfo_1p41[7],fitInfo_1p42[7],fitInfo_1p43[7],fitInfo_1p45[7]};


 double PosMag1[]={fitInfo1_1p35[6],fitInfo1_1p36[6],fitInfo1_1p37[6],fitInfo1_1p38[6],fitInfo1_1p383[6],fitInfoCor_1p384[6],fitInfoCor_1p385[6],fitInfo1_1p39[6],fitInfo1_1p40[6],fitInfo1_1p41[6],fitInfo1_1p42[6],fitInfo1_1p43[6],fitInfo1_1p45[6]};
 double PosErr1[]={fitInfo1_1p35[7],fitInfo1_1p36[7],fitInfo1_1p37[7],fitInfo1_1p38[7],fitInfo1_1p383[7],fitInfoCor_1p384[7],fitInfoCor_1p385[7],fitInfo1_1p39[7],fitInfo1_1p40[7],fitInfo1_1p41[7],fitInfo1_1p42[7],fitInfo1_1p43[7],fitInfo1_1p45[7]};

 double PosMagCor[]={fitInfoCor_1p35[6],fitInfoCor_1p36[6],fitInfoCor_1p37[6],fitInfoCor_1p38[6],fitInfoCor_1p383[6],fitInfoCor_1p384[6],fitInfoCor_1p385[6],fitInfoCor_1p39[6],fitInfoCor_1p40[6],fitInfoCor_1p41[6],fitInfoCor_1p42[6],fitInfoCor_1p43[6],fitInfoCor_1p45[6]};
 double PosCorErr[]={fitInfoCor_1p35[7],fitInfoCor_1p36[7],fitInfoCor_1p37[7],fitInfoCor_1p38[7],fitInfoCor_1p383[7],fitInfoCor_1p384[7],fitInfoCor_1p385[7],fitInfoCor_1p39[7],fitInfoCor_1p40[7],fitInfoCor_1p41[7],fitInfoCor_1p42[7],fitInfoCor_1p43[7],fitInfoCor_1p45[7]};
 
 double PosMagCor1[]={fitInfoCor1_1p35[6],fitInfoCor1_1p36[6],fitInfoCor1_1p37[6],fitInfoCor1_1p38[6],fitInfoCor1_1p383[6],fitInfoCor1_1p384[6],fitInfoCor1_1p385[6],fitInfoCor1_1p39[6],fitInfoCor1_1p40[6],fitInfoCor1_1p41[6],fitInfoCor1_1p42[6],fitInfoCor1_1p43[6],fitInfoCor1_1p45[6]};
 double PosCor1Err[]={fitInfoCor1_1p35[7],fitInfoCor1_1p36[7],fitInfoCor1_1p37[7],fitInfoCor1_1p38[7],fitInfoCor1_1p383[7],fitInfoCor1_1p384[7],fitInfoCor1_1p385[7],fitInfoCor1_1p39[7],fitInfoCor1_1p40[7],fitInfoCor1_1p41[7],fitInfoCor1_1p42[7],fitInfoCor1_1p43[7],fitInfoCor1_1p45[7]};
 

 double PosMagCor2[]={fitInfo1Cor_1p35[6],fitInfo1Cor_1p36[6],fitInfo1Cor_1p37[6],fitInfo1Cor_1p38[6],fitInfo1Cor_1p383[6],fitInfo1Cor_1p384[6],fitInfo1Cor_1p385[6],fitInfo1Cor_1p39[6],fitInfo1Cor_1p40[6],fitInfo1Cor_1p41[6],fitInfo1Cor_1p42[6],fitInfo1Cor_1p43[6],fitInfo1Cor_1p45[6]};
 double PosCor2Err[]={fitInfo1Cor_1p35[7],fitInfo1Cor_1p36[7],fitInfo1Cor_1p37[7],fitInfo1Cor_1p38[7],fitInfo1Cor_1p383[7],fitInfo1Cor_1p384[7],fitInfo1Cor_1p385[7],fitInfo1Cor_1p39[7],fitInfo1Cor_1p40[7],fitInfo1Cor_1p41[7],fitInfo1Cor_1p42[7],fitInfo1Cor_1p43[7],fitInfo1Cor_1p45[7]};


 
 double vgErr[]={0,0,0,0,0,0,0,0,0,0};

// TGraph *gR = new TGraphErrors(13,vg,PosMag,vgErr,PosErr);
// gR->Draw("a*"); 
// TGraph *gRcor = new TGraphErrors(13,vg,PosMagCor,vgErr,PosCorErr);
// gRcor->SetLineColor(kRed);
// gRcor->SetMarkerColor(kRed);
// gRcor->SetMarkerStyle(22);
// gRcor->Draw("p");

 TGraph *gR = new TGraph(13,vg,PosMag);
 gR->GetXaxis()->SetTitle("water_RI");
 gR->GetYaxis()->SetTitle("|X_{fit}|-|X_{MC}| [mm]");
 gR->SetMarkerStyle(21);
 gR->SetMarkerSize(2);
 gR->Draw("ap");

 TGraph *gRcor = new TGraph(13,vg,PosMagCor2);
 gRcor->SetMarkerSize(2);
 gRcor->SetMarkerColor(kRed);
 gRcor->SetMarkerStyle(23);
 gRcor->Draw("p");
/*
 TGraph *gRcor1 = new TGraph(13,vg,PosMagCor);
 gRcor1->SetMarkerSize(2);
 gRcor1->SetMarkerColor(kBlue);
 gRcor1->SetMarkerStyle(22);
 gRcor1->Draw("p");

 TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
 legend->AddEntry(gR,"MPW fitted without drive correction","p");
 legend->AddEntry(gRcor,"With drive correction Xcor = 1*Xfit+p1*Ufit","p");
 legend->AddEntry(gRcor1,"With drive correction Xcor = p0*Xfit+p1*Ufit","p");
 legend->Draw("same");
*/
}
