import csv
import glob
import ROOT
from ROOT import *
import numpy as np
from numpy import sqrt, pi
path = r'tables/*.csv'
files = glob.glob(path)
print(files)
print("There are", len(files), "files. Don't forget to count the total area!!")

dataFile1 = []
dataFile2 = []
dataFile3 = []
warning_points = 100 # print out a warning if found more than this value

Harea = TH1F("Harea","",1000,0,1000);
Hradius = TH1F("Hradius","",100,0,50);
Hmass = TH1F("Hmass","",100,0,50);

Harea.GetXaxis().SetTitle("Radius of dust [#mu m^{2}]")
Harea.GetYaxis().SetTitle("Number of particles (/1 #mu m^{2})")

Hradius.GetXaxis().SetTitle("Radius of dust [#mu m]")
Hradius.GetYaxis().SetTitle("Number of particles (/0.5 #mu m)")

Hmass.GetXaxis().SetTitle("Radius of dust [#mu m]")
Hmass.GetYaxis().SetTitle("Total mass contribution (ng/0.5 #mu m)")


totalArea = 0.8*1.2; ## mm
## 2720 x 1824 pixels
## 0.22 um/px

pix2um = 0.22

countPic = 0
for fname in files:
    data1 = []
    data2 = []
    data3 = []
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            area = float(row['Area'])
            area = area*pix2um*pix2um
            radius = sqrt(area/pi)
            volume = 4./3*pi*pow(radius,3)
            mass = volume*0.001
            data1.append(area)
            data2.append(radius)
            data3.append(mass)
            Hmass.Fill(radius, mass)
            if mass>15:
                print("???? warning: find large mass particles>15 ng, m=", round(mass,2), "ng. please check this file: ","topEdge"+fname[fname.index('Edge')+4:fname.index('Edge')+6])
            if area>1000:
                print("???? warning: find large area particles>1000 um^2, S=", round(area,2), "mm^2. please check this file: ","topEdge"+fname[fname.index('Edge')+4:fname.index('Edge')+6])
        if len(data1)>warning_points:
            print("warning: find ", len(data1), ", more than 500 dusts, please check this file: ",fname)
    dataFile1.append(data1)
    dataFile2.append(data2)
    dataFile3.append(data3)
    countPic += 1
## print(dataFile1)
countMass = 0
for fdata in dataFile1:
    for val in fdata:
        Harea.Fill(val)

for fdata in dataFile2:
    for val in fdata:
        Hradius.Fill(val)

#for fdata in dataFile3:
#    for val in fdata:
#        Hmass.Fill(radius, val)

countDust = Hradius.Integral()
countMass = Hmass.Integral()
print("Analyzed " + str(countPic) + " pictures")
print("Counts " + str(countDust) + " dusts.")
print("mass " + str(round(countMass,4)) + "ng.")
### totalArea is in mm^2, so switch to 1./100 cm^2
print("Density "+str(round(countMass/(totalArea/100*countPic),4))+" ng/cm^2")

ff = TFile("save_processed.root","recreate")
ff.cd()
Harea.Write()
Hradius.Write()
Hmass.Write()

#print(dataFile)
