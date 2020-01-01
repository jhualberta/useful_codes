#!/usr/bin/env python
import os, sys, string, ROOT, rat
import numpy as np
from ROOT import RAT
from ROOT import TF1, TH1F, TCanvas, gPad
dbMP = RAT.DB.Get()
dbMP.LoadAll(os.environ["GLG4DATA"], "FIT_MULTIPATH.ratdb")
mpScintName = ['labppo_0p25_oxford','labppo_0p5_oxford','labppo_1p0_oxford','labppo_2p0_oxford','labppo_6p0_oxford']
dbLink_MP = [dbMP.GetLink("FIT_MULTIPATH", kk) for kk in mpScintName]
binNum = 1600
#pdf_mpw_data = dbLink_MP.GetDArray("sPDF_multipathfit")
#fpdf_mpw = []
#hMPWpdf = TH1F("hMPWpdf","MPW pdf", binNum, -100, 300)
#for i in range(pdf_mpw_data.size()): # 1600 bins
#    fpdf_mpw.append(pdf_mpw_data[i])
#    hMPWpdf.SetBinContent(i+1, pdf_mpw_data[i])
# check the current existing scint pdf
pdf_scint_data = [dbLink_MP[i].GetDArray("sPDF_scint") for i in range(len(dbLink_MP))]
hscintOld = []
for i in range(len(dbLink_MP)):
   hscintOld.append(TH1F("hscint_"+mpScintName[i],"existed pdf in ratdb", binNum, -100, 300))
   for j in range(binNum):
     hscintOld[i].SetBinContent(j+1, pdf_scint_data[i][j])

pdf_MPW_shift = [dbLink_MP[i].GetDArray("sPDF_multipathfit_shift") for i in range(len(dbLink_MP))]
hMPWshift = []
for i in range(len(dbLink_MP)):
   hMPWshift.append(TH1F("hMPWshift_"+mpScintName[i],"existed pdf in ratdb", binNum, -100, 300))
   for j in range(binNum):
     hMPWshift[i].SetBinContent(j+1, pdf_MPW_shift[i][j])

c1 = ROOT.TCanvas()
gPad.SetLogy()
#plot the shifted MPW with MPW
#hMPWpdf.Draw()
#hMPWpdf_shift.Draw("sames")
#hscintOld[0].Draw()
#hscintOld[0].SetLineColor(1)
#for i in range(1,5):
#    hscintOld[i].SetLineColor(i+1)
#    hscintOld[i].Draw("sames")

hMPWshift[0].Draw()
hMPWshift[0].SetLineColor(1)
for i in range(1,5):
    hMPWshift[i].SetLineColor(i+1)
    hMPWshift[i].Draw("sames")

raw_input("Press 'Enter' to exit")

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetCanvasBorderMode(0)
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetPadColor(0)
ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
