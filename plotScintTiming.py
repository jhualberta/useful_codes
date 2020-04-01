#!/usr/bin/env python
import os, sys, string, ROOT#, rat, ConfigParser,  optparse
import numpy as np
from ROOT import TF1, TH1F, TCanvas, gPad

#def func_timing( x ):
#  y1 = a1*(np.exp(-x/-t1)-np.exp(-x/rt))/(-t1-rt)
#  y2 = a2*(np.exp(-x/-t2)-np.exp(-x/rt))/(-t2-rt)
#  y3 = a3*(np.exp(-x/-t3)-np.exp(-x/rt))/(-t3-rt)
#  y4 = a4*(np.exp(-x/-t4)-np.exp(-x/rt))/(-t4-rt)
#  return y1+y2+y3+y4

# simplified timing function, with 3 components
def func_timing_default(tau, amp, x ):
  rt = 0.8
  t1, t2, t3, t4 = tau 
  a1, a2, a3, a4 = amp
  y1 = a1*(np.exp(-x/-t1)-np.exp(-x/rt))/(-t1-rt)
  y2 = a2*(np.exp(-x/-t2)-np.exp(-x/rt))/(-t2-rt)
  y3 = a3*(np.exp(-x/-t3)-np.exp(-x/rt))/(-t3-rt)
  y4 = a4*(np.exp(-x/-t4)-np.exp(-x/rt))/(-t4-rt)
  return y1+y2+y3+y4

def func_timing_3pars( tau, amp, x ):
  rt = 0.8
  a1, a2, a3 = amp
  t1, t2, t3 = tau
  y1 = a1*(np.exp(-x/-t1)-np.exp(-x/rt))/(-t1-rt)
  y2 = a2*(np.exp(-x/-t2)-np.exp(-x/rt))/(-t2-rt)
  y3 = a3*(np.exp(-x/-t3)-np.exp(-x/rt))/(-t3-rt)
  return y1+y2+y3

def func_timing_oxford( tau, amp, x ):
  a1, a2, a3, aOffset = amp
  rt, t1, t2, t3 = tau
  y1 = a1*(np.exp(-x/-t1)-np.exp(-x/rt))/(-t1-rt)
  y2 = a2*(np.exp(-x/-t2)-np.exp(-x/rt))/(-t2-rt)
  y3 = a3*(np.exp(-x/-t3)-np.exp(-x/rt))/(-t3-rt)
  y4 = aOffset*(np.exp(-x/rt))/rt # a small additional amplitude for rise time
  return y1+y2+y3+y4

binNum = 400
binstep = 1
hScint_timing = TH1F("hScint_timing","scintillator timing, 0.5g/L PPO, penn", 400, -100, 300)
hScint_timing_default = TH1F("hScint_timing_default","scintillator timing, default, penn", 400, -100, 300)

hScint_0p25_oxford = TH1F("hScint_0p25_oxford","scintillator timing, 0.25g/L PPO, oxford", 400, -100, 300)
hScint_0p5_oxford = TH1F("hScint_0p5_oxford","scintillator timing, 0.5g/L PPO, oxford", 400, -100, 300)
hScint_1p0_oxford = TH1F("hScint_1p0_oxford","scintillator timing, 1.0g/L PPO, oxford", 400, -100, 300)
hScint_2p0_oxford = TH1F("hScint_2p0_oxford","scintillator timing, 2.0g/L PPO, oxford", 400, -100, 300)
hScint_6p0_oxford = TH1F("hScint_6p0_oxford","scintillator timing, 6.0g/L PPO, oxford", 400, -100, 300)

hpenn_0p5 = TH1F("hpenn_0p5", "Penn scint timing 0.5 g/L", 400, -100,300)
hpenn_alpha_0p5 = TH1F("hpenn_alpha_0p5", "Penn #alpha scint timing 0.5 g/L", 400, -100,300)

hdefault = TH1F("hdefault", "e- scint timing 2.0 g/L", 400, -100,300)
hdefault_alpha = TH1F("hdefault_alpha", "#alpha scint timing 2.0 g/L", 400, -100,300)


hOxford = [hScint_0p25_oxford,hScint_0p5_oxford,hScint_1p0_oxford,hScint_2p0_oxford,hScint_6p0_oxford]
# 0.25, 0.5, 1, 2, 6 g/L
# rise time, t1, t2, t3
tau025 = [1.25, -8.1, -25.0, -68.2]
tau05 = [1.12, -7.2, -18.7, -49.1]
tau1 = [1.18, -5.5, -13.3, -40.9]
tau2 = [1.06, -4.2, -11.7, -48.9]
tau6 = [0.94, -2.5, -9.3, -46.0]

tauAlpha_default = [-4.79, -18.4, -92.0, -900.0]
ampAlpha_default = [0.427, 0.313, 0.157, 0.1027]
# Penn 0.5 g/L
penn_0p5_amp  = [0.553, 0.331, 0.116]
penn_0p5_tau =  [-7.19, -24.81, -269.87]
penn_alpha_0p5_amp = [0.574, 0.311, 0.115]
penn_alpha_0p5_tau = [-6.56, -23.82, -224.19]

penn_default_tau = [-4.88, -15.4, -66.0, -400.0]
penn_default_amp = [0.665, 0.218, 0.083, 0.0346]

penn_alpha_default_tau =  [-4.79, -18.4, -92.0, -900.0]
penn_alpha_default_amp = [0.427, 0.313, 0.157, 0.1027]

amp025 =  [0.292, 0.531, 0.139, 0.038]
amp05 =  [0.435, 0.404, 0.126, 0.035]
amp1 =  [0.456, 0.375, 0.133, 0.036]
amp2 =  [0.579, 0.278, 0.089, 0.054]
amp6 =  [0.637, 0.170, 0.086, 0.107]

dataTau = [tau025,tau05,tau1,tau2,tau6]
dataAmp = [amp025,amp05,amp1,amp2,amp6]


for j in range(binNum): # 1600 bins
  xx = -100 + j*binstep
  if xx>=0: # this timing spectrum is always positive
    yy = func_timing_3pars(penn_0p5_tau, penn_0p5_amp, xx)
    yy1 = func_timing_3pars(penn_alpha_0p5_tau, penn_alpha_0p5_amp, xx)

    yy2 = func_timing_default(penn_default_tau, penn_default_amp, xx) 
    yy3 = func_timing_default(penn_alpha_default_tau, penn_alpha_default_amp, xx) 

    hpenn_0p5.Fill(xx, yy)
    hpenn_alpha_0p5.Fill(xx, yy1)
    hdefault.Fill(xx, yy2)
    hdefault_alpha.Fill(xx, yy3)

  else:
    hpenn_0p5.Fill(xx,0)
    hpenn_alpha_0p5.Fill(xx,0)
    hdefault.Fill(xx, 0)
    hdefault_alpha.Fill(xx, 0)



#for j in range(binNum): # 1600 bins
#  xx = -100 + j*binstep
#  if xx>=0: # this timing spectrum is always positive
#    yy = func_timing_3pars(xx)
#    hScint_timing.Fill(xx, yy)
#  else:
#    hScint_timing.Fill(xx, 0)
#
#for j in range(binNum): # 1600 bins
#  xx = -100 + j*binstep
#  if xx>=0: # this timing spectrum is always positive
#    yy = func_timing_default(xx)
#    hScint_timing_default.Fill(xx, yy)
#  else:
#    hScint_timing_default.Fill(xx, 0)
#
## oxford measurements
#for item in range(0,5):
#  for j in range(binNum): # 1600 bins
#    xx = -100 + j*binstep
#    if xx>=0: # this timing spectrum is always positive
#      yy = func_timing_3pars(xx)
#      yy2 = func_timing_oxford(dataTau[item], dataAmp[item], xx)
#      hOxford[item].Fill(xx, yy2)
#    else:
#      hOxford[item].Fill(xx, 0)  

c1 = ROOT.TCanvas()
gPad.SetLogy()
#plot the shifted MPW with MPW

#hpenn_alpha_0p5.SetLineColor(2)
#hpenn_alpha_0p5.Draw()
#hpenn_0p5.Draw("sames")
hdefault_alpha.SetLineColor(2)
hdefault_alpha.Draw()
hdefault.Draw("sames")


#hOxford[0].SetLineColor(1)
#hOxford[1].SetLineColor(2)
#hOxford[2].SetLineColor(3)
#hOxford[3].SetLineColor(4)
#hOxford[4].SetLineColor(6)
#
#scale = hScint_timing.GetMaximum()
#scale2 = hOxford[1].GetMaximum()
##hScint_timing_oxford.Scale(scale/scale2)
##hScint_timing_oxford.Scale(1./hScint_timing_oxford.Integral())
##hScint_timing.Scale(1./hScint_timing.Integral())
#
#
#hOxford[4].Draw()
#for j in range(0,4):
#    hOxford[j].Draw("sames")

#hOxford[1].SetLineColor(2)
#hOxford[3].SetLineColor(4)
#hScint_timing.SetLineColor(2)
#hScint_timing_default.SetLineColor(4)
#hScint_timing.SetLineStyle(2)
#hScint_timing_default.SetLineStyle(2)
#
#hOxford[3].Draw()
#hScint_timing_default.Draw("sames")
#hOxford[1].Draw("sames")
#hScint_timing.Draw("sames")

raw_input("Press 'Enter' to exit")
