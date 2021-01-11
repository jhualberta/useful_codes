#/usr/bin/python
# sim_2p5MeVbeta_z0_splitZ-1000_1.root
import sys, os
import subprocess
path=os.getcwd()
file_list = os.listdir(path)
fList = open('fmpw.dat')
processfile = []
count = 0
for i in fList:
 processfile.append(i.rstrip())
 count = count+1

for fname in processfile:
 subprocess.call(["root","-q" ,"AnalysisSolarMPW.C+(\"" +fname+ "\")"])
