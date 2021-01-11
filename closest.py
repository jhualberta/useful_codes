import numpy as np 
import ROOT
import sys, os, getopt
from ROOT import TCanvas, TTree, TFile, TH1F, AddressOf
from array import array
import random
import matplotlib.pyplot as plt

table = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 650, 600, 700, 525, 625, 600, 550, 600, 600, 675, 925, 1150, 725, 1000, 1100, 1000, 1000, 1025, 800, 1075, 1550, 1150, 1200, 1625, 1325, 1525, 1250, 1475, 1625, 1500, 1250, 1750, 1824, 1575, 2124, 1950, 1924, 1575, 1525, 1924, 1525, 2250, 2624, 2175, 2224, 2475, 2075, 2000, 2075, 2424, 2275, 2400, 2824, 2875, 2250, 2800, 2300, 2775, 2924, 2400, 3300, 3275, 3175, 2024, 3175, 2850, 3150, 3075, 3275, 3224, 3224, 3400, 3500, 3400, 3575, 3700, 3975, 3700, 3400, 3700]

fname = 'Merged_fEnergy_Rat6176_test_6and7MeVgamma.root'
binN = 240
ff = TFile(fname,"READ")
hEsum = ff.Get("hNCherenkov")
def closest(lst, K):
  lst = np.asarray(lst) 
  idx = (np.abs(lst - K)).argmin() 
  return [lst[idx],idx]


hSx = TH1F("hSx","",100,0,10)
for i in range(100):
  ycounts = hEsum.GetBinContent(i);
  nphoton = hEsum.GetXaxis().GetBinCenter(i);
  xx = closest(table, nphoton)[1]
  print xx, ycounts
  hSx.Fill(xx/10,ycounts)

hSx.Draw()
raw_input()
