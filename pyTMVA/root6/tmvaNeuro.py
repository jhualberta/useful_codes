import ROOT
import sys, os, getopt
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array
path = os.getcwd()
file_list = os.listdir(path)

inputfile = "rawData.csv"
f0 = open(inputfile)
reader = csv.reader(f0)
csvdata = list(reader)
csvdata = csvdata[1:] #delete the table title
totalevents = len(csvdata)

# create a TNtuple
ntuple = ROOT.TNtuple("ntuple","neuro data","x1:x2:x3:x4:x5:x6:x7:x8:signal")#"brain:synek:signal")
gos = []
brain = []
synek = []
n20 =[]
indexIII = []
indexV = []
indexIII_V_I_III = []
four = []
gcs = []
lesion = []

count = 0
# generate 'signal' and 'background' distributions
for i in range(totalevents):
  goodData = True 
  for j in range(8):
    check1 = (csvdata[i][j]!='NaN')
    goodData = goodData*check1
    
  if goodData:
    valGos = float(csvdata[i][0])
    valBrain = float(csvdata[i][1])
    valSynek = float(csvdata[i][2])
    valN20 = float(csvdata[i][3])
    valIndexIII = float(csvdata[i][4])
    valIndexV = float(csvdata[i][5])
    valIndexII_V_I_III  = float(csvdata[i][6])
    valfour = float(csvdata[i][7])
    valGcs = float(csvdata[i][8])
    gos.append(valGos)
    brain.append(valBrain)
    synek.append(valSynek)
    ntuple.Fill(valBrain, valSynek, valN20, valIndexIII, valIndexV, valIndexII_V_I_III, valfour, valGcs, valGos)
    count += 1 

nTrain = int(count*0.7) # set 70% of the raw data for training
nTest = count - nTrain
print "total events:", count, "trained:",nTrain, "test:",nTest


# keeps objects otherwise removed by garbage collected in a list
gcSaver = []
 
# create a new TCanvas
gcSaver.append(ROOT.TCanvas())
 
# draw an empty 2D histogram for the axes
histo = ROOT.TH2F("histo","",1,-5,5,1,-5,5)
histo.Draw()
 
# draw the signal events in red
ntuple.SetMarkerColor(ROOT.kRed)
ntuple.Draw("x2:x1","signal > 0.5","same")
 
# draw the background events in blue
ntuple.SetMarkerColor(ROOT.kBlue)
ntuple.Draw("x2:x1","signal <= 0.5","same")

TMVA.Tools.Instance()

# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.
fout = ROOT.TFile("test.root","RECREATE")

factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=D;D;D;D;D;D;D;D",#I;D;P;G,D",
                                "AnalysisType=Classification"]
                                     ))
dataloader = TMVA.DataLoader("dataset")

dataloader.AddVariable("x1","F")
dataloader.AddVariable("x2","F")
dataloader.AddVariable("x3","F")
dataloader.AddVariable("x4","F")
dataloader.AddVariable("x5","F")
dataloader.AddVariable("x6","F")
dataloader.AddVariable("x7","F")
dataloader.AddVariable("x8","F")

dataloader.AddSignalTree(ntuple)
dataloader.AddBackgroundTree(ntuple)

# cuts defining the signal and background sample
sigCut = ROOT.TCut("signal > 0.5")
bgCut = ROOT.TCut("signal <= 0.5")

dataloader.PrepareTrainingAndTestTree(sigCut,   # signal events
                                   bgCut,    # background events
                                   ":".join([
                                        "nTrain_Signal=0",
                                        "nTrain_Background=0",
                                        "SplitMode=Random",
                                        "NormMode=NumEvents",
                                        "!V"
                                       ]))


method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", ":".join(["!H","!V", "NTrees=nTrain", "nEventsMin=nTest", "MaxDepth=3", "BoostType=AdaBoost", "AdaBoostBeta=0.5",
 "SeparationType=GiniIndex","nCuts=10", "PruneMethod=NoPruning", ]))
                                   
dataloader.TrainAllMethods()
dataloader.TestAllMethods()
dataloader.EvaluateAllMethods()

reader = ROOT.TMVA.Reader()
import array
varx1 = array.array('f',[0]) ; reader.AddVariable("x1",varx1)
varx2 = array.array('f',[0]) ; reader.AddVariable("x2",varx2)
varx3 = array.array('f',[0]) ; reader.AddVariable("x3",varx3)
varx4 = array.array('f',[0]) ; reader.AddVariable("x4",varx4)
varx5 = array.array('f',[0]) ; reader.AddVariable("x5",varx5)
varx6 = array.array('f',[0]) ; reader.AddVariable("x6",varx6)
varx7 = array.array('f',[0]) ; reader.AddVariable("x7",varx7)
varx8 = array.array('f',[0]) ; reader.AddVariable("x8",varx8)

reader.BookMVA("BDT","weights/TMVAClassification_BDT.weights.xml")
# create a new 2D histogram with fine binning
histo2 = ROOT.TH2F("histo2","",200,-5,5,200,-5,5)
 
# loop over the bins of a 2D histogram
for i in range(1,histo2.GetNbinsX() + 1):
    for j in range(1,histo2.GetNbinsY() + 1):
         
        # find the bin center coordinates
        varx1[0] = histo2.GetXaxis().GetBinCenter(i)
        varx2[0] = histo2.GetYaxis().GetBinCenter(j)
         
        # calculate the value of the classifier
        # function at the given coordinate
        bdtOutput = reader.EvaluateMVA("BDT")
         
        # set the bin content equal to the classifier output
        histo2.SetBinContent(i,j,bdtOutput)
 
gcSaver.append(ROOT.TCanvas())
histo2.Draw("colz")
 
# draw sigma contours around means
for mean, color in (
    ((1,1), ROOT.kRed), # signal
    ((-1,-1), ROOT.kBlue), # background
    ):
     
    # draw contours at 1 and 2 sigmas
    for numSigmas in (1,2):
        circle = ROOT.TEllipse(mean[0], mean[1], numSigmas)
        circle.SetFillStyle(0)
        circle.SetLineColor(color)
        circle.SetLineWidth(2)
        circle.Draw()
        gcSaver.append(circle)
 
ROOT.gPad.Modified()
# fill histograms for signal and background from the test sample tree
ROOT.TestTree.Draw("BDT>>hSig(22,-1.1,1.1)","classID == 0","goff")  # signal
ROOT.TestTree.Draw("BDT>>hBg(22,-1.1,1.1)","classID == 1", "goff")  # background
 
ROOT.hSig.SetLineColor(ROOT.kRed); ROOT.hSig.SetLineWidth(2)  # signal histogram
ROOT.hBg.SetLineColor(ROOT.kBlue); ROOT.hBg.SetLineWidth(2)   # background histogram
 
# use a THStack to show both histograms
hs = ROOT.THStack("hs","")
hs.Add(ROOT.hSig)
hs.Add(ROOT.hBg)
 
# show the histograms
gcSaver.append(ROOT.TCanvas())
hs.Draw()

rocgraph = TMVA.ROCCurve("dataset");
rocgraph.Draw()
