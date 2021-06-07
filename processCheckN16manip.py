#/usr/bin/python
# sim_2p5MeVbeta_z0_splitZ-1000_1.root
import sys, os
import subprocess
import ROOT
from ROOT import TH1F, TFile, TTree
path=os.getcwd()
file_list = os.listdir(path)
fList = open('listMC.dat')
processfile = []
count = 0
processfile = []
for i in fList:
 data = []   
 line = i.rstrip() 
 data = [i.split(',')[0],i.split(',')[1],i.split(',')[2],i.split(',')[3]]
 processfile.append(data)
 count = count+1

for item in processfile:
    ff = TFile(item[0],"read")
    t1 = ff.Get("T")
    hx = TH1F("hx","",2000,-9000,9000)
    hy = TH1F("hy","",2000,-9000,9000)
    hz = TH1F("hz","",2000,-9000,9000)
    t1.Project("hx","posx","")
    t1.Project("hy","posy","")
    t1.Project("hz","posz","")
    meanx = hx.GetMean()    
    meany = hy.GetMean()    
    meanz = hz.GetMean()
    print item[0]
    print "Recon: ", meanx, ", ", meany, ", ", meanz-108
    print "manip: ", item[1], ", ", item[2], ", ", item[3]
    manipx = float(item[1])
    manipy = float(item[2])
    manipz = float(item[3])


    if (abs(meanx - manipx)>1000 or abs(meany - manipy)>1000 or abs(meanz - manipz)>1000 ):
      print "warning", meanx - manipx, meany - manipy, meanz - manipz
