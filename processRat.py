#/usr/bin/python
# sim_2p5MeVbeta_z0_splitZ-1000_1.root
import sys, os
import subprocess
path=os.getcwd()
file_list = os.listdir(path)
fList = open('fList.dat')
processfile = []
count = 0
for i in fList:
 processfile.append(i.rstrip())
 count = count+1

for fname in processfile:
 #print "now process: "+ fname
 i1 = fname.index("Z")
 ## find the third "_"
 i2 = fname.index("_", fname.index("_",fname.index("_")+1)+1)
 ## find the fourth "_"
 i3 = fname.index("_",i2+1)
 print "z = " + fname[i1+1:i3] # print z position
# subprocess.call(["root","-lq" ,"\'getMean.C(\"" +fname+ "\")\' "])
 subprocess.call(["root","-lq" ,"getMean.C(\"" +fname+ "\")"])

