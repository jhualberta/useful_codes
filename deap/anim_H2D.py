import sys
import glob
from rat import *
from ROOT import *
import numpy as np
from numpy import sqrt
import time
from time import sleep
import struct
#from run_numbers import runs
from datetime import datetime, time as datetime_time, timedelta

f = TFile("CombinedHistos_levelCuts_run28103to28849.root");##Combined_levelCutsNTP_run28103to28849.root");
h2d_0 = f.Get("H2_nSCBayes_rprompt60Bayes_0");	
h2d_1 = f.Get("H2_nSCBayes_rprompt60Bayes_1");	
h2d_2 = f.Get("H2_nSCBayes_rprompt60Bayes_2");	
h2d_3 = f.Get("H2_nSCBayes_rprompt60Bayes_3");	
h2d_4 = f.Get("H2_nSCBayes_rprompt60Bayes_4");	
h2d_5 = f.Get("H2_nSCBayes_rprompt60Bayes_5");	
h2d_6 = f.Get("H2_nSCBayes_rprompt60Bayes_6");	
h2d_7 = f.Get("H2_nSCBayes_rprompt60Bayes_7");	
h2d_8 = f.Get("H2_nSCBayes_rprompt60Bayes_7_neck");	
h2d_9 = f.Get("H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe");
h2d_10 = f.Get("H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe_mbR");	
h2d_11 = f.Get("H2_nSCBayes_rprompt60Bayes_7_str_qRatio");
h2d_12 = f.Get("H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ");
h2d_13 = f.Get("H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ_pulseG");

h2d_0.SetTitle("no cuts")
h2d_1.SetTitle("level 1: dtmTrigSrc&0x82 == 0")
h2d_2.SetTitle("level 2: calcut&0x31f8 == 0")
h2d_3.SetTitle("level 3: deltat > 20000")
h2d_4.SetTitle("level 4: numEarlyPulses <= 3")
h2d_5.SetTitle("level 5: subeventN==1")
h2d_6.SetTitle("level 6: 2250<eventTime < 2700")
h2d_7.SetTitle("level 7: qPE>60")
h2d_8.SetTitle("7-level && neckVetoN==0")
h2d_9.SetTitle("7-level && neckVetoN==0 && fmaxpe<0.4")
h2d_10.SetTitle("7-level && neckVetoN==0 && fmaxpe<0.4 && mbR<630")
h2d_11.SetTitle("7-level && str && qRatio")
h2d_12.SetTitle("7-level && str && qRatio && mbZ<550")
h2d_13.SetTitle("7-level && str && qRatio && mbZ<550 && pulseindexfirstgar>2")

h2d_0.GetXaxis().SetRangeUser(0,400);h2d_0.GetYaxis().SetRangeUser(0.5,1.0);
h2d_1.GetXaxis().SetRangeUser(0,400);h2d_1.GetYaxis().SetRangeUser(0.5,1.0);
h2d_2.GetXaxis().SetRangeUser(0,400);h2d_2.GetYaxis().SetRangeUser(0.5,1.0);
h2d_3.GetXaxis().SetRangeUser(0,400);h2d_3.GetYaxis().SetRangeUser(0.5,1.0);
h2d_4.GetXaxis().SetRangeUser(0,400);h2d_4.GetYaxis().SetRangeUser(0.5,1.0);
h2d_5.GetXaxis().SetRangeUser(0,400);h2d_5.GetYaxis().SetRangeUser(0.5,1.0);
h2d_6.GetXaxis().SetRangeUser(0,400);h2d_6.GetYaxis().SetRangeUser(0.5,1.0);
h2d_7.GetXaxis().SetRangeUser(0,400);h2d_7.GetYaxis().SetRangeUser(0.5,1.0);
h2d_8.GetXaxis().SetRangeUser(0,400);h2d_8.GetYaxis().SetRangeUser(0.5,1.0);
h2d_9.GetXaxis().SetRangeUser(0,400);h2d_9.GetYaxis().SetRangeUser(0.5,1.0);
h2d_10.GetXaxis().SetRangeUser(0,400);h2d_10.GetYaxis().SetRangeUser(0.5,1.0);
h2d_11.GetXaxis().SetRangeUser(0,400);h2d_11.GetYaxis().SetRangeUser(0.5,1.0);
h2d_12.GetXaxis().SetRangeUser(0,400);h2d_12.GetYaxis().SetRangeUser(0.5,1.0);
h2d_13.GetXaxis().SetRangeUser(0,400);h2d_13.GetYaxis().SetRangeUser(0.5,1.0);

h2d_0.GetYaxis().SetTitleOffset(1.5)
h2d_1.GetYaxis().SetTitleOffset(1.5)
h2d_2.GetYaxis().SetTitleOffset(1.5)
h2d_3.GetYaxis().SetTitleOffset(1.5)
h2d_4.GetYaxis().SetTitleOffset(1.5)
h2d_5.GetYaxis().SetTitleOffset(1.5)
h2d_6.GetYaxis().SetTitleOffset(1.5)
h2d_7.GetYaxis().SetTitleOffset(1.5)
h2d_8.GetYaxis().SetTitleOffset(1.5)
h2d_9.GetYaxis().SetTitleOffset(1.5)
h2d_10.GetYaxis().SetTitleOffset(1.5)
h2d_11.GetYaxis().SetTitleOffset(1.5)
h2d_12.GetYaxis().SetTitleOffset(1.5)
h2d_13.GetYaxis().SetTitleOffset(1.5)

ffroi = TFile("roi_802_days_11March2020_nsc_rp60.root");
mycut = ffroi.Get("roi_top55_gr");
mycut.SetName("mycut");
mycut.SetLineColor(kRed);

plotNum = 14;
c1 = TCanvas("c1","The HSUM example",1100,800);
c1.SetGrid();
#c1.Divide(2,1);
gBenchmark.Start("hsum")

h2d = []
h2d.append(h2d_0); h2d.append(h2d_1); h2d.append(h2d_2); h2d.append(h2d_3); h2d.append(h2d_4); h2d.append(h2d_5); h2d.append(h2d_6); h2d.append(h2d_7);
h2d.append(h2d_8); h2d.append(h2d_9); h2d.append(h2d_10); h2d.append(h2d_11); h2d.append(h2d_12); h2d.append(h2d_13);

slider = 0;
#gSystem.Unlink("levelCuts2D.gif");

gSystem.Unlink("smallRegion_levelCuts2D.gif");
###Fill histograms randomly
kUPDATE = 1;
gifcnt = 0;

### h2d_2.Draw()
for i in range(plotNum):
   if ((i%kUPDATE) == 0): 
     gPad.SetLogz();
     h2d[i].Draw("colz");
     mycut.Draw("same");
     c1.Update();
     ## if (slider) slider.SetRange(0,Float_t(i)/20.);
     c1.Modified();
     c1.Update();
     #if gROOT.IsBatch():
     #c1.Print("levelCuts2D.gif+200");## default gif refresh time = 10ms, now is 50*10 ms = 0.5 s
     c1.Print("smallRegion_levelCuts2D.gif+200");
     sleep(1)
## slider.SetRange(0,1);
#h2d[0].GetYaxis().SetLabelOffset(1.5)
h2d[0].Draw("sameaxis"); ## to redraw axis hidden by the fill area

## c1.Modified();
### make infinite animation by adding "++" to the file name
### if (gROOT.IsBatch()):
##c1.Print("levelCuts2D.gif++")
c1.Print("smallRegion_levelCuts2D.gif++")

##You can view the animated file hsumanim.gif with Netscape/IE or mozilla
gBenchmark.Show("hsum");