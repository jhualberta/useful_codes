{
double range = 3;
const Int_t n = 30;
Double_t x[n+1],y[n+1];
Double_t rcut = 0.85*range;
Double_t dphi = TMath::TwoPi()/n;
for (Int_t i=0;i<n;i++) {
x[i] = rcut*TMath::Cos(i*dphi);
y[i] = rcut*TMath::Sin(i*dphi);
}
x[n] = x[0]; y[n] = y[0];
TCutG *mycut = new TCutG("mycut",n+1,x,y);
TF2 *f2 = new TF2("f2","1-x^2-y^2+(x^2+y^2)^2",-range,range,-range,range);
f2->Draw("surf1 [mycut]");
}
