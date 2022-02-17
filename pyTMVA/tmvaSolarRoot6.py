import ROOT
import sys, os, getopt
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
from ROOT.TMVA import ROCCurve
import numpy as np
import csv
from array import array
path = os.getcwd()
file_list = os.listdir(path)

inputfile = "../MixedDataSet_MP_r200004to200658.root"
#MixedDataSet_MPW_SolarBi214_lowE_nhit15.root"#MixedMPWDataSet_solarBi214_lowE.root"
f0 = TFile(inputfile)

# create a TNtuple

methodString = "BDT" 
#methodString = "BDTF"
#methodString = "Likelihood"



ntuple = f0.Get("mixedBkg")
#ROOT.TNtuple("ntuple","solar+bkg","x1:x2:x3:x4:x5:x6:signal")
# beta14:thetaij:nhits:udotR:logL:scaleLogL:signal

beta14 = []
thetaij = []
nhits = []
energy = []
udotR = []
scaledLogL = []
signal = []

count = 0
# generate 'signal' and 'background' distributions
totalEvts = ntuple.GetEntries()
for i in range(totalEvts):
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
#ntuple.Draw("beta14:thetaij","signal > 0.0","same")
 
# draw the background events in blue
ntuple.SetMarkerColor(ROOT.kBlue)
#ntuple.Draw("beta14:thetaij","signal < 0.0","same")

TMVA.Tools.Instance()

# note that it seems to be mandatory to have an
# output file, just passing None to TMVA::Factory(..)
# does not work. Make sure you don't overwrite an
# existing file.
fout = ROOT.TFile("testSolar_results.root","RECREATE")

factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=D;D;D;D;D;D;",
                                #"Transformations=D;D;D;D;D;D;",#I;D;P;G,D",
                                "AnalysisType=Classification"]
                                     ))
# beta14:thetaij:nhits:udotR:logL:scaleLogL

dataloader = ROOT.TMVA.DataLoader("dataset")

#dataloader.AddVariable("itr","F")
#dataloader.AddVariable("beta14","F")
#dataloader.AddVariable("thetaij","F")
#dataloader.AddVariable("energy","F")
dataloader.AddVariable("klDiv","F")
dataloader.AddVariable("udotR","F")
dataloader.AddVariable("Utest","F")
dataloader.AddVariable("Gtest","F")
dataloader.AddVariable("zfactor","F")
dataloader.AddVariable("scaleLogL","F")

dataloader.SetBackgroundWeightExpression( "weight" )

#factory.AddVariable("itr","F")
#factory.AddVariable("beta14","F")
#factory.AddVariable("thetaij","F")
##factory.AddVariable("energyCor","F")
#factory.AddVariable("udotR","F")
##factory.AddVariable("logL","F")
##factory.AddVariable("scaleLogL","F")

dataloader.AddSignalTree(ntuple)
dataloader.AddBackgroundTree(ntuple)

# cuts defining the signal and background sample
sigCut = ROOT.TCut("signal > 0.0  ")
bgCut = ROOT.TCut("signal < 0.0")

dataloader.PrepareTrainingAndTestTree(sigCut,   # signal events
                                   bgCut,    # background events
                                   ":".join([
                                        "nTrain_Signal=0",
                                        "nTrain_Background=0",
                                        "SplitMode=Random",
                                        "NormMode=NumEvents",
                                        "!V"
                                       ]))
if methodString == "Likelihood":
   method = factory.BookMethod( dataloader, ROOT.TMVA.Types.kLikelihood, "Likelihood", 
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=100:NSmoothBkg[0]=100:NSmoothBkg[1]=100:NSmooth=1:NAvEvtPerBin=100" );

if methodString == "BDT":
   method = factory.BookMethod( dataloader, ROOT.TMVA.Types.kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )
#method = factory.BookMethod(dataloader,ROOT.TMVA.Types.kBDT, methodString, ":".join(["!H","!V","TransformOutput","PDFInterpol=Spline2","NTrees=nTrain", "nEventsMin=nTest", "MaxDepth=3", "BoostType=AdaBoost", "AdaBoostBeta=0.5","SeparationType=GiniIndex","nCuts=10", "PruneMethod=NoPruning", ]))


if methodString == "BDTF":
   method = factory.BookMethod( dataloader, ROOT.TMVA.Types.kBDT, "BDTF",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

reader = ROOT.TMVA.Reader()
import array

#dataloader.AddVariable("klDiv","F")
#dataloader.AddVariable("udotR","F")
#dataloader.AddVariable("Utest","F")
#dataloader.AddVariable("Gtest","F")
#dataloader.AddVariable("zfactor","F")
#dataloader.AddVariable("scaleLogL","F")

varx1 = array.array('f',[0]) ; reader.AddVariable("klDiv",varx1)
varx2 = array.array('f',[0]) ; reader.AddVariable("udotR",varx2)
varx3 = array.array('f',[0]) ; reader.AddVariable("Utest",varx3)
varx4 = array.array('f',[0]) ; reader.AddVariable("Gtest",varx4)
varx5 = array.array('f',[0]) ; reader.AddVariable("zfactor",varx5)
varx6 = array.array('f',[0]) ; reader.AddVariable("scaleLogL",varx6)

#reader.BookMVA(methodString,"weights/TMVAClassification_BDT.weights.xml")
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
        #bdtOutput = reader.EvaluateMVA(methodString)##
         
        # set the bin content equal to the classifier output
        #histo2.SetBinContent(i,j,bdtOutput)
 
gcSaver.append(ROOT.TCanvas())
#histo2.Draw("colz")
 
# draw sigma contours around means
#for mean, color in (
#    ((1,1), ROOT.kRed), # signal
#    ((-1,-1), ROOT.kBlue), # background
#    ):
#     
#    # draw contours at 1 and 2 sigmas
#    for numSigmas in (1,2):
#        circle = ROOT.TEllipse(mean[0], mean[1], numSigmas)
#        circle.SetFillStyle(0)
#        circle.SetLineColor(color)
#        circle.SetLineWidth(2)
#        circle.Draw()
#        gcSaver.append(circle)
 
ROOT.gPad.Modified()
# fill histograms for signal and background from the test sample tree
#ROOT.TestTree.Draw(methodString+">>hSig(22,-1.1,1.1)","classID == 0","gon")#"goff")  # signal
#ROOT.TestTree.Draw(methodString+">>hBg(22,-1.1,1.1)","classID == 1", "gon")#"goff")  # background
 
#ROOT.hSig.SetLineColor(ROOT.kRed); ROOT.hSig.SetLineWidth(2)  # signal histogram
#ROOT.hBg.SetLineColor(ROOT.kBlue); ROOT.hBg.SetLineWidth(2)   # background histogram
 
# use a THStack to show both histograms
#hs = ROOT.THStack("hs","")
#hs.Add(ROOT.hSig)
#hs.Add(ROOT.hBg)
 
# show the histograms
#gcSaver.append(ROOT.TCanvas())
#hs.Draw()
raw_input("enter")
#rocgraph = factory.GetROCCurve();
#auc = rocgraph.GetROCIntegral()
#print "AUC curve!! ", auc

#rocgraph.Draw()
