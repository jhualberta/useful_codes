#!/usr/bin/env python
import os, sys, string, ROOT, rat
import numpy as np
from ROOT import RAT
from ROOT import TF1, TH1F, TGraph, TCanvas, gPad
ndata = 7
speed = np.array([170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0])
# 0.5
#bias = np.array([-315.644,-197.462,-82.3972,35.3163,144.416,245.884,337.388])
#1.0
#bias = np.array([-298.667,-179.491,-60.2557,56.9042,167.542,266.731,357.572])
#2.0
bias =  np.array([-309.227,-185.951,-65.1289,53.1423,164.702,265.039, 354.629])
#6.0


linear = TGraph(ndata, speed, bias)
linear.SetMarkerStyle(21)
linear.Draw("AP*")
linear.Fit("pol1")
func =  linear.GetListOfFunctions().FindObject("pol1");
zero = func.GetX(0,speed[0],speed[ndata-1],1e-5)
print zero, 299.792/zero
raw_input("hold")


