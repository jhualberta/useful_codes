from ROOT import *
from array import array
import ctypes
x = [357.4794, 614.5129, 614.5129, 357.4794, 22.3105, -393.5785, -600.7243, -600.7243, -393.5785, 22.3105]
y = [529.9850, 176.2087, -176.2087, -529.985, -638.888, -503.7575, -218.6458, 218.6458, 503.7575, 638.8880]
z = [633.5017]

xx = array('f', x)
yy = array('f', y)

g1 = TGraph(10, xx, yy)
g1.GetXaxis().SetTitle("x [mm]")
g1.GetYaxis().SetTitle("y [mm]")

#l = TLatex();
#l.SetTextSize(0.03);
#l.SetTextFont(42);
#l.SetTextAlign(21);
#l.SetTextColor(kBlue);

g1.Draw("AP")

#def drawtext():
n = g1.GetN();
x1 = ctypes.c_double() 
x2 = ctypes.c_double()
height = ctypes.c_double(3.0)
labels = []
for i in range(n): 
   g1.GetPoint(i, x1, x2)
   #print x1, x2
   j = 30+i
   l = TText(x1, x2.value+height.value, str(j))
   labels.append(l)  
   l.Draw("same")

#for ll in labels:
#    ll.Draw("same")
#ex = TExec("ex","drawtext()")
#g1.GetListOfFunctions().Add(ex);
g1.SetMarkerStyle(8)
#g1.Draw("AP")
raw_input()
