#/usr/bin/python
# sim_2p5MeVbeta_z0_splitZ-1000_1.root
import sys, os
import subprocess
path=os.getcwd()
file_list = os.listdir(path)
count = 0
for i in range(100,200):
 subprocess.call(["root","-q" ,"AnalyWaterTres_allRat6176_"+str(i)+".C+"])
