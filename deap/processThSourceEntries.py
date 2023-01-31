#/usr/bin/python
import sys,os
import shutil
import collections
import subprocess
from subprocess import *
from ROOT import *
from array import array
import numpy as np
from numpy import mean, std 

database = [(28103, 33562.2), (28105, 51166.7), (28110, 80554.3), (28119, 82929.4), (28123, 82706), (28132, 73977.5), (28137, 91422.5), (28147, 102921), (28160, 82727.9), (28164, 84241.9), (28173, 83011.8), (28188, 84856.2), (28192, 37543.3), (28193, 42662.2), (28203, 86068.1), (28208, 71545.9), (28218, 81955.7), (28222, 76137), (28225, 13144.8), (28234, 84796.8), (28239, 86912.7), (28248, 84832.1), (28252, 84456.4), (28261, 84567.6), (28266, 82621.2), (28275, 65775.2), (28279, 82917.7), (28288, 83391.6), (28294, 13550.4), (28295, 69157.3), (28310, 97976.1), (28314, 83849.4), (28323, 35333.2), (28325, 49740), (28330, 8242.77), (28331, 72731), (28343, 80154.1), (28348, 101548), (28357, 81105.3), (28362, 42269.8), (28363, 38520.8), (28372, 31747.3), (28374, 33233.2), (28375, 16737.9), (28379, 34726.6), (28381, 3552.17), (28383, 37236.6), (28393, 86581.7), (28398, 79217.8), (28409, 80981.9), (28413, 81907.9), (28422, 93353.4), (28427, 73278.8), (28436, 88788.7), (28440, 79075), (28449, 81333.6), (28454, 83239.7), (28464, 78085.2), (28468, 21203), (28472, 59867.8), (28481, 83806.3), (28486, 66983.8), (28487, 1466.59), (28512, 85700.5), (28516, 107079), (28525, 87512), (28530, 79275), (28545, 79737.2), (28549, 103215), (28565, 82179.1), (28570, 88141.2), (28579, 73621.4), (28583, 83467.9), (28593, 83246.2), (28595, 68020.6), (28608, 75458.5), (28624, 29596), (28632, 76061.5), (28637, 92863.4), (28655, 80712.5), (28659, 5483.64), (28660, 80960.5), (28670, 81830.4), (28676, 104177), (28685, 72406.1), (28689, 83168.9), (28702, 78084.5), (28708, 76079.3), (28730, 10423.1), (28734, 69218.3), (28743, 83800.5), (28748, 82275.2), (28758, 82719.1), (28763, 87651.6), (28786, 51708.9), (28795, 83852.8), (28800, 82276.1), (28809, 67379.8), (28813, 19081.8), (28815, 63646), (28824, 81088.6), (28829, 80787.7), (28838, 77165.1), (28844, 8456.76), (28847, 35769.3), (28849, 37971.4)]

print "there are ", len(database), "runs."

path=os.getcwd()
file_list = os.listdir(path)
fList = open('fileRuns.dat')
processfile = []
count = 0
for i in fList:
  processfile.append(i.rstrip())
  count = count+1

print "now processing ", count, " files."

ndata = len(database)

if count != ndata:
    print "data size not match."

runStart = 30627
runEnd = 32149

H2_qpe_fprompt_7 = TH2D( "H2_qpe_fprompt_7", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 );
H2_qpe_fprompt_cuts = TH2D( "H2_qpe_fprompt_cuts", "fmaxpe<0.75; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 );
H2_qpe_fprompt_cuts1 = TH2D( "H2_qpe_fprompt_cuts1", "fmaxpe<0.4; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 );
H2_qpe_fprompt_cutROI = TH2D( "H2_qpe_fprompt_cutROI", "MBR<630, fmaxpe<0.4; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 1000, 0.0, 1.0 );
H2_qpe_fprompt_cut2ROI = TH2D( "H2_qpe_fprompt_cut2ROI", "MBR<630, fmaxpe<0.4, cft2R<0.04, cfb3r<0.1, pulseindexfirstgar>2; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 1000, 0.0, 1.0 );


array_rates = array('d')
array_runs = array('d')

i = 0
count = 0
count481 = 0
count465 = 0
runLowRate = []
for fname in processfile:
    # Merged_31561.root
    number = fname[fname.index('run')+3:fname.index('run')+8]
    run = int(number)
    print run
    runcheck = database[i][0]
    if runcheck != run:
        print "data size bug"
        break
    ff = TFile(fname)
    tree = ff.Get("T")
    n_level7 = tree.GetEntries()
    liveTime = database[i][1]
    rate = float(n_level7)/liveTime
    array_rates.append(rate)
    array_runs.append(run)
    i += 1
    count += 1

print "strangeRate runs", runLowRate

print "Mean rates=", mean(array_rates), "+-", std(array_rates)

gRates = TGraph(count, array_runs, array_rates)
gRates.SetLineColor( 2 )
gRates.SetLineWidth( 3 )
gRates.SetMarkerColor( 2 )
gRates.SetMarkerStyle( 21 )
gRates.SetTitle( 'Total runs' )
gRates.GetXaxis().SetTitle( 'run number' )
gRates.GetYaxis().SetTitle( 'counts/sec' )
gRates.Draw( 'AP' )


raw_input()
