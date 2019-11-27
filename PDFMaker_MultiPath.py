#!/usr/bin/env python
import os, sys, string, ROOT, rat, ConfigParser,  optparse
import numpy as np
from ROOT import RAT
from ROOT import TF1, TH1F, TCanvas, gPad
# Author Jie Hu - 24 Apr 2019 <jhu9@ualberta.ca>
# Create the pdfs for the MultiPath ScintWater(Partial) Fitter, Scint Fitter or Te-Scint Fitter
# It prints out the data of pdfs and make plots.
# usage: python PDFMaker_MultiPath.py
# enter "Enter" to quit

# Convolve the scint timing histogram with the MPW pdf histogram
# NOTE: Convolving the different regions of the two pdfs gives different shapes of the convolved histograms.
# The selection of convolving regions can be tuned here.

def pdfMaker(options):

  # Search OPTICS files in the defaults_MP.ini by options
  config = ConfigParser.ConfigParser()
  config.read("defaults_MP.ini")
  scintIndex = config.get("unloaded_scint", "default_material") #labppo by default
  testAlpha = False

  if options.opticsFile == None:
      options.opticsFile = config.get("unloaded_scint", "default_optics")
      # set by default
  elif options.opticsFile == "labppo_dda":
      options.opticsFile = config.get("labppo_dda", "default_optics")
      scintIndex = config.get("labppo_dda", "default_material")
      print "!!! Choose OPTICS file and scint index: ", options.opticsFile, scintIndex
  elif options.opticsFile == "te_dda":
      options.opticsFile = config.get("TeDiol_0p5_dda", "default_optics")
      scintIndex = config.get("TeDiol_0p5_dda", "default_material")
      print "!!! Choose OPTICS file and scint index: ", options.opticsFile, scintIndex
  elif options.opticsFile == "te_bismsb_dda":
      options.opticsFile = config.get("TeDiol_0p5_bismsb_dda", "default_optics")
      scintIndex = config.get("TeDiol_0p5_bismsb_dda", "default_material")
      print "!!! Choose OPTICS file and scint index: ", options.opticsFile, scintIndex
  if options.particle == None:
      print "use e- by default"
  elif options.particle == "alpha":
      print "!!! Choose alpha particle"
      testAlpha = True

  ## Extract the scintillator timing parameters
  dbOPTICS = RAT.DB.Get()
  # labppo by default,  add options for the other OPTICS
  dbOPTICS.LoadAll(os.environ["GLG4DATA"], options.opticsFile)
  dbLink_OPTICS = dbOPTICS.GetLink("OPTICS", scintIndex)
  # e- timing by default
  timing_params = dbLink_OPTICS.GetDArray("SCINTWAVEFORM_value1")
  amplitude = dbLink_OPTICS.GetDArray("SCINTWAVEFORM_value2")
  rise_time = dbLink_OPTICS.GetD("SCINT_RISE_TIME")
  # alpha timing
  if testAlpha:
      timing_params = dbLink_OPTICS.GetDArray("SCINTWAVEFORMalpha_value1")
      amplitude = dbLink_OPTICS.GetDArray("SCINTWAVEFORMalpha_value2")
      rise_time = dbLink_OPTICS.GetD("SCINT_RISE_TIME")

  timing_data = [timing_params, amplitude, rise_time]

  binNum = 1600 # use 1600 bins for pdf histograms
  ## Extract the MPW pdf as the PMT time response
  dbMP = RAT.DB.Get()
  dbMP.LoadAll(os.environ["GLG4DATA"],"FIT_MULTIPATH.ratdb")
  dbLink_MP = dbMP.GetLink("FIT_MULTIPATH","")
  pdf_mpw_data = dbLink_MP.GetDArray("sPDF_multipathfit")
  fpdf_mpw = []
  hMPWpdf = TH1F("hMPWpdf","MPW pdf", binNum, -100, 300)
  for i in range(pdf_mpw_data.size()): # 1600 bins
      fpdf_mpw.append(pdf_mpw_data[i])
      hMPWpdf.SetBinContent(i+1, pdf_mpw_data[i])
  
  # check the current existing scint pdf
  pdf_scint_data = dbLink_MP.GetDArray("sPDF_scint")
  hscintOld = TH1F("hscintOld","existed pdf in ratdb", binNum, -100, 300)
  for i in range(binNum):
      hscintOld.SetBinContent(i+1, pdf_scint_data[i])
 
  # Build-up the scint timing histogram
  t1, t2, t3, t4 = timing_data[0]
  a1, a2, a3, a4 = timing_data[1]
  rt = timing_data[2]
  print "timing constants: ", t1, t2, t3, t4
  print "amplitudes: ", a1, a2, a3, a4
 
  def func_timing( x ):
      y1 = a1*(np.exp(-x/-t1)-np.exp(-x/rt))/(-t1-rt)
      y2 = a2*(np.exp(-x/-t2)-np.exp(-x/rt))/(-t2-rt)
      y3 = a3*(np.exp(-x/-t3)-np.exp(-x/rt))/(-t3-rt)
      y4 = a4*(np.exp(-x/-t4)-np.exp(-x/rt))/(-t4-rt)
      return y1+y2+y3+y4
 
  hScint_timing = TH1F("hScint_timing","scintillator timing", binNum, -100, 300)
  binstep = float(400./binNum) # 0.25 ns step for a whole range of 400 ns
  for j in range(binNum): # 1600 bins
      xx = -100 + j*binstep
      if xx>=0: # this timing spectrum is always positive
          yy = func_timing(xx)
          hScint_timing.Fill(xx, yy)
      else:
          hScint_timing.Fill(xx, 0)

  ftiming = []
  for j in range(binNum): # 1600 bins
      ftiming.append( hScint_timing.GetBinContent(j+1) )

  print "Now convolving and coordinating the pdfs:"
  start_convl = -6.25 # instead of -6, tuned to get a best fit result.
  end_convl = 10 # late light are deweighted after 10 ns in the MPW pdf
  # Set the region of the MPW pdf for the convolution.
  loCutBin = hMPWpdf.FindBin(start_convl)
  hiCutBin = hMPWpdf.FindBin(end_convl)
  # Take the t>0 region for the scint timing spectrum.
  nonzeroBin = hScint_timing.FindBin(0.0)
 
  hConvl = TH1F("hConvl","convolved histogram", binNum, -100, 300)
 
  # the numpy.convolve(list1, list2) method can be used here and it is much faster. However, it will introduce larger bin numbers and make the bins of the histograms tricky.
  for j in range(binNum):
      convl = 0; xi = 0; sxi = 0; rxi = 0
      xx = -100 + j*binstep
      xi = -100
      for ii in range(binNum):
          ind = hMPWpdf.FindBin(xx - xi)
          sxi = hMPWpdf.GetBinContent(ind)
          if ind<loCutBin or ind>hiCutBin: # only convolve the selected region of the MPW pdf, ignoring the deweighting components in the MPW pdf
              sxi = 0
          if xi>=0: # only convolve the positive region of the scint timing spectrum
              rxi = func_timing(xi)
          else: 
              rxi = 0
          convl = convl+ rxi*sxi
          xi = xi+binstep
      hConvl.Fill(xx,convl)
  
  # De-weighting for early time [-100, -6] ns and late time [220, 300] ns
  # These regions are the same to the ET1D pdf
  earlyTime = start_convl
  lateTime = 220 # this value is fixed
  
  earlyBin = hConvl.FindBin(earlyTime)
  lateBin = hConvl.FindBin(lateTime)
  # the earlyBin and the earlyBin+1 have 0 bin contents, caused by the convolution, so starting at earlyBin+2
  earlyDeweight = hConvl.GetBinContent(earlyBin+2) # use the second bin content on the right to the last bin
  lateDeweight = hConvl.GetBinContent(lateBin)
  
  hScintPDF = TH1F("hScintPDF","convolved histogram, deweighted", binNum, -100, 300)
  for i in range(1, earlyBin+2):
      hScintPDF.SetBinContent(i, earlyDeweight)
 
  for i in range(earlyBin+2, lateBin):
      hScintPDF.SetBinContent(i, hConvl.GetBinContent(i))
 
  for i in range(lateBin, binNum+1):
      hScintPDF.SetBinContent(i, lateDeweight)
 
  # Plot and check histograms
  # Find the half-maximum bins of the two pdfs, shift the MPW pdf to the right to match with the half-maximum bins.
  # Calculate half-maximum
  startTime = -9.5 #NOTE: This is also a tuning parameter, defining from where to start to calculate the half maximum of the MPW pdf.
                   # This parameter affects how much the pdf needed to be shifted.
                   # the value -9.5 ns is the last step of the early noise level in the MPW pdf. The half maximum is calculate from here.
  mpwStartBin = hMPWpdf.FindBin(startTime)
  halfMax_mpw = 0.5*( hMPWpdf.GetMaximum() + hMPWpdf.GetBinContent(mpwStartBin) )
  # Looking for the bin close to get the half-maximum value, mpw
  halfval1 = []
  for i in range( mpwStartBin, hMPWpdf.GetMaximumBin() ):
      val = hScintPDF.GetBinContent(i)
      halfval1.append( abs(val-halfMax_mpw) )
  minIndex1 = halfval1.index(min(halfval1))
  halfMaxBin_mpw = minIndex1 + mpwStartBin
 
  # Looking for the bin close to get the half-maximum value, scint
  hScintPDF.Scale(hMPWpdf.GetMaximum()/hScintPDF.GetMaximum()) # scale the scint pdf to the same maximum to the MPW pdf
  startleveltemp = []
  startlevel = hMPWpdf.GetBinContent(mpwStartBin) #start at the level close to the MPW start bin. Then the value of the half-max is close to the MPW
  scintStartBin = hScintPDF.FindBin(earlyTime)
  for i in range( scintStartBin, hScintPDF.GetMaximumBin() ):
      val = hScintPDF.GetBinContent(i)
      startleveltemp.append(abs(val-startlevel))
  minIndexLevel = scintStartBin + startleveltemp.index(min(startleveltemp))
  scintStartBin1 = minIndexLevel
  halfMax_scint = 0.5*( hScintPDF.GetMaximum() + hScintPDF.GetBinContent(scintStartBin1) )
 
  halfval2 = []
  for i in range( scintStartBin, hScintPDF.GetMaximumBin() ):
      val = hScintPDF.GetBinContent(i)
      halfval2.append(abs(val-halfMax_scint))
  minIndex2 = halfval2.index( min(halfval2) )
  halfMaxBin_scint = minIndex2 + scintStartBin1
 
  shiftBin = abs( halfMaxBin_scint - halfMaxBin_mpw )
  print "check MPW startbin, scint startbin, halfMax_mpw, halfMax_scint, halfMaxBin_scint, halfMaxBin_mpw, shifted ", mpwStartBin, scintStartBin, halfMax_mpw, halfMax_scint, halfMaxBin_mpw, halfMaxBin_scint, shiftBin
  # move the MPW pdf to the right
  hMPWpdf_shift = TH1F("hMPWpdf_shift","MPW pdf", binNum, -100, 300)
  for i in range(shiftBin, binNum+1):
      hMPWpdf_shift.SetBinContent(i, hMPWpdf.GetBinContent(i-shiftBin))
 
  for i in range(0, shiftBin):
      hMPWpdf_shift.SetBinContent(i+1, hMPWpdf.GetBinContent(1)) # fill the left blanks with the noise level
 
  # Finally, both the pdfs are scaled to their integrals for the likelihood calculations.
  scalePDF = hMPWpdf_shift.Integral()
  hScintPDF.Scale(scalePDF/hScintPDF.Integral())
 
  sPDF_multipath_shift = []
  sPDF_scint = []
 
  for i in range(binNum):
      sPDF_multipath_shift.append(hMPWpdf_shift.GetBinContent(i+1))
      sPDF_scint.append(hScintPDF.GetBinContent(i+1))
 
  # Print out the data of the pdfs, which are put into FIT_MULTIPATH.ratdb
  print "Please put the following data into FIT_MULTIPATH.ratdb"
  print "/// MPW pdf shifted peak for partial fill"
  print "sPDF_multipathfit_shift: ", sPDF_multipath_shift, ","
  print "/// scint pdf, scint timing (Chicago) convolved with MPW pdf, -6 to 220 ns"
  print "sPDF_scint: ", sPDF_scint, ","
 
  # Print out the peaked time
  sPDF_multipath_shift_maxbin = sPDF_multipath_shift.index(max(sPDF_multipath_shift))
  sPDF_scint_maxbin = sPDF_scint.index(max(sPDF_scint))
  print "sPDF_multipath_shift is peaked at bin: ",  sPDF_multipath_shift_maxbin, str(-100 + sPDF_multipath_shift_maxbin*0.25) + " ns"
  print "sPDF_scint is peaked at bin: ", sPDF_scint_maxbin, str(-100 + sPDF_scint_maxbin*0.25) + " ns"
  # Plot and check histograms
  c1 = ROOT.TCanvas()
  c1.Divide(2, 2)
  c1.cd(1)
  gPad.SetLogy()
  #plot the shifted MPW with MPW
  hMPWpdf_shift.SetLineColor(2)
  hMPWpdf.Draw()
  hMPWpdf_shift.Draw("sames")
  c1.cd(2)
  #plot the scintillation timing spectrum
  gPad.SetLogy()
  hScint_timing.Draw()
  c1.cd(3)
  #compare the existing scint pdf with the created one
  gPad.SetLogy()
  hscintOld.Draw()
  hScintPDF.Draw("sames")
  c1.cd(4)
  #plot the shifted MPW with scint pdf
  gPad.SetLogy()
  hMPWpdf_shift.SetLineColor(4)
  hScintPDF.SetLineColor(2)
  hMPWpdf_shift.Draw()
  hScintPDF.Draw("sames")
 
  raw_input("Press 'Enter' to exit")
 
  ROOT.gROOT.SetStyle("Plain")
  ROOT.gStyle.SetCanvasBorderMode(0)
  ROOT.gStyle.SetPadBorderMode(0)
  ROOT.gStyle.SetPadColor(0)
  ROOT.gStyle.SetCanvasColor(0)
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetOptStat(0)

if __name__ == '__main__':
    parser = optparse.OptionParser(usage = "usage: %prog [options] target", version = "%prog 1.0")
    parser.add_option("-s", type = "string", dest = "opticsFile", help = "Load an OPTICS.ratdb, default='OPTICS_labppo.ratdb'", default = "OPTICS_labppo.ratdb")
    parser.add_option("-p", type = "string", dest = "particle", help = "Particle type to use (see generator documentation for available particles), default = 'e-'", default = "e-")
    (options, args) = parser.parse_args()
    pdfMaker(options)
