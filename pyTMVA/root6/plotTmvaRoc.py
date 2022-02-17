import ROOT
import sys, os, getopt
from ROOT import TMVA, TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array
path = os.getcwd()
file_list = os.listdir(path)


rocgraph = TMVA.GetROCCurve(dataset);
rocgraph.Draw()
