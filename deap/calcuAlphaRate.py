import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import rat 
from ROOT import *
from array import array

checkLevel = 2 

#import statistics
#362 remove bad: 31863, 31975, 31871, 31265
run = [30627,30636,30642,30646,30651,30655,30658,30663,30672,30676,30681,30686,30690,30695,30699,30700,30705,30706,30707,30708,30712,30717,30721,30726,30730,30731,30732,30735,30736,30741,30742,30743,30744,30746,30747,30751,30756,30760,30765,30769,30774,30778,30779,30784,30785,30789,30794,30799,30804,30808,30813,30814,30815,30819,30824,30826,30832,30837,30841,30846,30850,30855,30858,30859,30865,30870,30876,30881,30889,30894,30898,30903,30905,30909,30916,30922,30923,30926,30929,30936,30943,30948,30953,30961,30962,30963,30964,30965,30969,30974,30978,30983,30987,30992,31002,31003,31008,31012,31016,31020,31025,31029,31034,31038,31043,31053,31057,31063,31067,31072,31076,31081,31082,31083,31084,31085,31086,31094,31099,31103,31113,31114,31118,31123,31124,31125,31130,31132,31137,31141,31146,31150,31155,31159,31161,31162,31163,31166,31170,31171,31172,31182,31186,31192,31196,31201,31205,31206,31211,31216,31224,31228,31233,31241,31246,31250,31251,31256,31260,31265,31275,31279,31284,31285,31286,31287,31291,31296,31297,31301,31302,31304,31309,31313,31318,31322,31327,31331,31334,31339,31343,31351,31355,31360,31361,31365,31370,31374,31379,31382,31386,31391,31395,31403,31405,31409,31410,31411,31416,31420,31425,31429,31435,31437,31448,31452,31458,31464,31470,31476,31480,31485,31489,31495,31499,31505,31511,31527,31528,31533,31538,31542,31547,31552,31557,31561,31573,31577,31582,31586,31591,31595,31600,31601,31602,31606,31611,31612,31613,31614,31615,31619,31622,31627,31629,31633,31634,31636,31638,31643,31647,31652,31656,31662,31666,31671,31675,31676,31681,31684,31688,31693,31697,31702,31706,31711,31715,31722,31727,31728,31729,31731,31733,31735,31737,31739,31747,31752,31757,31761,31766,31770,31775,31779,31780,31786,31794,31797,31798,31806,31809,31813,31818,31822,31827,31831,31836,31840,31844,31854,31858,31863,31871,31874,31877,31882,31886,31891,31895,31902,31908,31916,31922,31923,31930,31936,31944,31945,31947,31949,31951,31953,31955,31957,31959,31964,31970,31975,31983,31987,31992,31998,32003,32008,32022,32026,32031,32035,32042,32046,32056,32061,32065,32070,32072,32076,32082,32086,32088,32093,32097,32101,32106,32110,32115,32119,32121,32131,32139,32144,32148,32149]

livetime = [62394.9,21943.2,60512.7,83095,73466.6,17714.3,71308.7,3222.85,61320.3,88991.8,13770.6,75230.4,73088.6,79612.2,77065.7,7243.32,270.931,2349.43,12.8163,85115,72370.8,82933.2,83841.8,81547,47860.3,5070.23,826.497,421.752,25834.8,44492.9,43.9193,375.716,1469.93,160.438,34032.4,80918,80770.1,83923.7,80991.2,82698.8,81060.1,1917.54,81483.2,10462.3,86314.2,67413.3,81395.5,93703.2,79808,68538.6,12928.7,60676.7,13759.9,74455,1257.71,79253.2,78631.8,80863.2,83606.7,80518.4,83136.3,16053.1,140.95,76768.8,78020.6,79269.3,84146.3,65718.6,66158.1,80542.9,84624.4,39209.5,46023.8,76824.5,84658.6,88179.1,603.69,5091.48,429.794,63734.3,74614.2,81372.9,87823.2,28011.4,41776.1,692.611,1.94703,17981.1,72118.9,76340.4,82887.5,80673.3,87531.3,76816.8,28606.1,2977.54,80772.1,83359.3,81241.1,84577.8,98278.8,64160.8,103650,61862.8,104204,78113.3,83888.3,87658.8,72147.5,81111,81036.7,42511.8,958.46,285.827,812.906,1109.15,11445.6,79323.4,87384.3,89181.2,41846.7,2771.61,79299.9,2943.55,152.458,79937.5,1481.39,74587.4,89393.7,78607.3,80092,83610.9,2995.24,6313.59,507.136,367.961,920.92,60959.9,80443,4994,350.646,77146.4,74367.7,78420.3,76285.1,79436.8,4500.37,84405.3,76995.5,83069.9,81395.2,79430.2,79124.8,83539.5,76698.6,81583.3,23124.6,59849.1,84634.3,2114.66,67712.1,83773.7,76968.3,1002.28,409.173,1423.84,78992.3,84169.3,601.46,644.341,1410.32,80525.4,80490.6,99462.5,65689.5,82664.2,79432.9,15921.6,73555.5,91010.1,79320.1,66580.8,83462,78369.2,471.164,113021,54163.9,90315.6,12996.8,52805.3,83768.5,84134.2,78587,57745.2,14236.4,16847.8,65733,6928.7,77870.6,82435.4,82884.4,86504.1,18081.6,64874.9,207.861,72645.4,79948.4,54685.5,26534.2,49061.4,85829.9,80827.9,84373.9,79515,81021.1,81670,79590.1,26750.7,36232.1,82923.8,80945.7,84652.2,79400.1,83805,96558.9,65201.9,64901,81517.9,80496.4,83101.7,97836.7,67382.6,76246.5,79.1604,4124.05,79427.8,6325.44,8384.6,61972.3,2348.86,3327.93,40107.3,44308.4,55020.7,25977.8,3807.79,392.928,38037.4,42678.9,73504.2,82833.5,107163,55914.9,80321.5,84165.5,80927.9,28971.5,54452.9,7504.87,65045.6,83508.2,80784.1,83865.4,90213.2,89099.3,71294.2,78193.4,25721.3,14059.6,84816.2,438.243,82467.9,85408.7,84410.4,84903.6,87808,282.282,73262.7,104412,58517.2,80014.6,83490.3,79524.2,3700.12,90475.8,69398.5,67437.8,2075.08,9097.7,497.467,64073.8,85912.3,78824.5,83901.8,84229.7,77986.6,79835.1,20804.2,24938.8,55094.3,82220.7,-0.998378,273.974,3952.62,70696.7,82191.5,82917.5,81082.5,102239,74058.9,71811.4,35627.8,696.547,80650,80510,83353.6,80115.6,84056.1,113585,69708.4,76701.7,77842,83571.5,94897.7,48948.4,64551,39549.7,-0.998378,67139.2,85555.8,78762.8,120602,41166.3,85709.8,60965.1,82341.3,81457.1,86798.2,78647.8,2529.34,68120,77803.5,82761,65507.3,26866.8,88313.9,65050.7,8902.8,71087.5,16292.7,62217.8,93708.8,78244.2,86164,69821.1,33461,15719.3,76597.7,77509,80069.2,11394.2,88916.9]

nRuns = len(run)
## live time after muon-veto correction

###print "Already removed bad rates, total t=",sum(livetime)/3600/24, "days", "total runs", len(livetime)

day = [float(k)/3600/24 for k in livetime]

print len(day)

rawdata = zip(run,day)
#362 -4 remove bad: 31265, 31863, 31871, 31975

deleteData = []
for data in rawdata:
   if data[0] == 31265: deleteData.append(data)
   if data[0] == 31863: deleteData.append(data)
   if data[0] == 31871: deleteData.append(data) 
   if data[0] == 31975: deleteData.append(data)

for data in deleteData:
   rawdata.remove(data)

#print len(rawdata)

n1 = []
n2 = []
fileList = []

count = 0
path=os.getcwd()
fList = open('list358runs.dat')
fileList = []
fileZombie = []
for i in fList:
   fileList.append(i.rstrip())
   # print i.rstrip()
   count = count+1
print "now checking", count, "root files."

bad_runs = []
### check livetime long enough!!!

countGood = 0
runGood = []
livetimeGood = []
i = 0
processedData = []
hnscbVsRprompt = TH2F("hnscbVsRprompt", "nSCBayes vs rprompt60Bayes", 2000, 0, 2000, 100, 0, 1)
hqpeVsFprompt = TH2F("hqpeVsFprompt", "qpe vs fprompt", 2000, 0, 2000, 100, 0, 1)

data_variables_single = []
data_variables_double = []

for data in rawdata:
     fname = "Merged_OptLowThresh_" + str(data[0]) + ".root"
     ff = TFile(fname)
     if ff.IsZombie():
          print "zombie file", fname 
          fileZombie.append(fname)
          continue

     if data[1]<0.01:
          ## print "check --- run with too short livetime", data[0], "run time", data[1], "day"
          continue
     t1 = ff.Get("T1")
     t2 = ff.Get("T2")

     qpe1 = array('f',[0]) #unsigned int
     nSCBayes1 = array('f',[0]) #unsigned double 
     fprompt1 = array('f',[0]) #unsigned double 
     fmaxpe1 = array('f',[0]) #unsigned double 
     rprompt1 = array('f',[0]) #unsigned double 
     eventTime1 = array('f',[0])
     subeventN1 = array('f',[0])
     nhit1 = array('i',[0])

     t1.SetBranchAddress("qpe",qpe1)
     t1.SetBranchAddress("fprompt",fprompt1)
     t1.SetBranchAddress("nSCBayes",nSCBayes1)
     t1.SetBranchAddress("rprompt60Bayes",rprompt1)
     t1.SetBranchAddress("eventTime", eventTime1)
     t1.SetBranchAddress("subeventN", subeventN1)
     t1.SetBranchAddress("nhit", nhit1)
  
     qpe2 = array('f',[0]) #unsigned int
     nSCBayes2 = array('f',[0]) #unsigned double 
     fprompt2 = array('f',[0]) #unsigned double 
     fmaxpe2 = array('f',[0]) #unsigned double 
     rprompt2 = array('f',[0]) #unsigned double 
     eventTime2 = array('f',[0])
     subeventN2 = array('f',[0])
     nhit2 = array('i',[0])

     t2.SetBranchAddress("qpe",qpe2)
     t2.SetBranchAddress("fprompt",fprompt2)
     t2.SetBranchAddress("nSCBayes",nSCBayes2)
     t2.SetBranchAddress("rprompt60Bayes",rprompt2)
     t2.SetBranchAddress("eventTime", eventTime2)
     t2.SetBranchAddress("subeventN", subeventN2)
     t2.SetBranchAddress("nhit", nhit2)

     hncluster1 = TH1F("hncluster1","",3,0,3)
     hncluster2 = TH1F("hncluster2","",3,0,3)

     ## level of cuts
     cutLevel4 = "numEarlyPulses <=3"          
     cutLevel5 = "&& subeventN <= 1"
     cutLevel6 = "&& eventTime > 2250 && eventTime < 2700"
     #cutLevel7 = "&& fmaxpe<0.4"
     cutLevel7 = "&& qpe>200"
     cutLevel8 = "&& fmaxpe<0.4"
     cutLevel9 = "&& neckVeto == 0"
     cutLevel10 = "&& sqrt(evtx*evtx+evty*evty+evtz*evtz)<830"
     cutLevel11 = "&& pulseindexfirstgar > 2"
     cutLevel12 = "&& cft2r < 0.04 && evtz < 550"
     cutLevel13 = "&& cfb3r<0.1"
     #cutLevel14 = ""## mb-tf2 Z
     #cutLevel15 = "" ## mb-tf2 R
     t1.Project("hncluster1", "ncluster")
     t2.Project("hncluster2", "ncluster")
     t1.Project("hnscbVsRprompt", "rprompt60Bayes:nSCBayes")
     t2.Project("hqpeVsFprompt", "fprompt:qpe")
     var = (nscb, rprompt, qpe, fprompt, nhit)

     if checkLevel == 4:
        t1.Project("hncluster1", "ncluster", cutLevel4)
        t2.Project("hncluster2", "ncluster", cutLevel4)
        if numEarlyPulses<=3:
            data_variables.append(var)

     if checkLevel == 5:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5)
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5)
        if numEarlyPulses<=3 && subeventN <= 1:
            data_variables.append(var)

     if checkLevel == 6:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6)
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6)

     if checkLevel == 7:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7)
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7)
     
     if checkLevel == 8: 
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8 )
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8 )
 
     if checkLevel == 9:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9 )
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9 )
 
     if checkLevel == 10:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10 )
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10 )

     if checkLevel == 11:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10+cutLevel11 )
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10+cutLevel11 )
     
     if checkLevel == 12:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10+cutLevel11+cutLevel12 )
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10+cutLevel11+cutLevel12 )
     
     if checkLevel == 13:
        t1.Project("hncluster1", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10+cutLevel11+cutLevel12+cutLevel13 )
        t2.Project("hncluster2", "ncluster", cutLevel4+cutLevel5+cutLevel6+cutLevel7+cutLevel8+cutLevel9+cutLevel10+cutLevel11+cutLevel12+cutLevel13 )

     n1val = hncluster1.Integral()
     n2val = hncluster2.Integral()
     ## n1val = hh.GetBinContent(hh.GetXaxis().FindBin(1.5))
     ## n2val = hh.GetBinContent(hh.GetXaxis().FindBin(2.5))
     n1.append(n1val)
     n2.append(n2val)
     runGood.append(data[0])
     livetimeGood.append(round(data[1],3)) ### !! unit:days
     processedData.append(data)
     countGood += 1

     ## print runname, n1val, n2val


print "Good runs for analysis:", countGood
#### odd runs
## 30681, 30686, 30743, 30751, 30785
print "after removing live time too short:", sum(livetimeGood), "days"

pars =[1.5, 0.0625, 1.0625, 0.5]
#haojie
##pars = [9/6.1, 1/6.1, 7/6.1, 1/3.05]
nAV = []
nTPB = []
runCut = []
badRuns = []
nSingle = array('d')
nDouble = array('d')
actualLiveTime = [] ### removing wrong calculated ones
bad_nSingle = array('d') ## still to save n1, n2 even it's bad for calculating nAV and nTPB
bad_nDouble = array('d')
for i in range(countGood):
    a = n1[i]
    b = n2[i]
    d = float(livetimeGood[i])
    #if d<0.2:
    #    print "run too short",d
    #    continue
    if not(a>0) or not(b>0):
        ## print "why nClusters == 0?? n1=", a, "n2=", b
        continue
    av = round( (pars[0]*b - pars[1]*a)/d, 3)
    tpb = round( (pars[2]*a - pars[3]*b)/d, 3)
    if av<0 or tpb<0:
        ## print "wrong calculations in AV rates, run", run[i], av, tpb, "days", d
        ## still saves nClusters even the derived nAV and nTPB are wrong
        bad_nSingle.append( round(a/d, 3) )
        bad_nDouble.append( round(b/d, 3) )
        badRuns.append(processedData[i][0])
        continue
    nAV.append(av)
    nTPB.append(tpb)
    nSingle.append(a/d)
    nDouble.append(b/d)

    bad_nSingle.append( round(a/d, 3) )
    bad_nDouble.append( round(b/d, 3) )

    runCut.append(processedData[i][0])
    actualLiveTime.append(processedData[i][1])

    if av>1000:
        print "n1, n2, nAV, ntpb",a/d,b/d,av,tpb, "time", d

print "total=",round(sum(day),2),"days"
print "after removing too short runs", sum(livetimeGood), "days"
print "after removing too short runs and wrong calculated runs", sum(actualLiveTime), "days"
print "average nAV rates over runs (event/day)",  round(np.mean(nAV) , 3), "+-", round(np.std(nAV) ,3)
print "convert to mHz", round(np.mean(nAV)/3600/24*1000 , 2), "+-", round(np.std(nAV)/3600/24*1000 ,2)
print "average nTPB rates over runs (event/day)", round(np.mean(nTPB), 3), "+-", round(np.std(nTPB),3)
print "convert to mHz", round(np.mean(nTPB)/3600/24*1000, 2), "+-", round(np.std(nTPB)/3600/24*1000,2)

print "Cluster average n1 rates over runs", round(np.mean(nSingle),3), "+-", round(np.std(nSingle),3)
print "Cluster convert to mHz", round(np.mean(nSingle)/3600/24*1000 , 2), "+-", round(np.std(nSingle)/3600/24*1000 ,2)
print "Cluster average n2 rates over runs", round(np.mean(nDouble),3), "+-", round(np.std(nDouble),3)
print "Cluster convert to mHz", round(np.mean(nDouble)/3600/24*1000, 2), "+-", round(np.std(nDouble)/3600/24*1000,2)

print "including wrong calculation!!!"
print "Cluster average n1 rates over runs", round(np.mean(bad_nSingle),3), "+-", round(np.std(bad_nSingle),3)
print "Cluster convert to mHz", round(np.mean(bad_nSingle)/3600/24*1000 , 2), "+-", round(np.std(bad_nSingle)/3600/24*1000 ,2)
print "Cluster average n2 rates over runs", round(np.mean(bad_nDouble),3), "+-", round(np.std(bad_nDouble),3)
print "Cluster convert to mHz", round(np.mean(bad_nDouble)/3600/24*1000, 2), "+-", round(np.std(bad_nDouble)/3600/24*1000,2)


#plt.scatter(runCut, nAV,alpha=0.8, marker='s')
#plt.scatter(runCut, nTPB, alpha=0.8, marker='o')
#plt.xlabel("run")
#plt.ylabel("events/day")
#plt.errorbar(runCut, nAV, yerr=np.std(nAV))
#plt.errorbar(runCut, nTPB, yerr=np.std(nTPB))
#plt.grid()
#plt.legend(["AV","TPB"], loc='upper right')#, loc='upper left')
#plt.show()



hnscbVsRprompt.Draw("colz")

raw_input("enter")
