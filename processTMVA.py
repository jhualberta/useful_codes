#/usr/bin/python
# sim_2p5MeVbeta_z0_splitZ-1000_1.root
import sys, os
import subprocess
import shutil
path=os.getcwd()
file_list = os.listdir(path)
count = 0
NN = 5000 

for i in range(0,NN):
   filename = 'analyzeData.cc'
   ftemp = open(filename,'r+')
   flist=ftemp.readlines()
   flist[20] = '  TFile* histFile = new TFile(\"OutputMinus0p1_BDT'+str(i)+'_\"+fileName, \"RECREATE\");\n'
   flist[75]= '  TTree *T2 = (TTree*)inputFile->Get(\"T2_'+str(i)+'\");\n'
   flist[119]= '  TTree* sig = dynamic_cast<TTree*>(inputFile->Get(\"T2_'+str(i)+'\"));\n'
   ftemp=open(filename,'w+')
   ftemp.writelines(flist)
   ftemp.close()
   #subprocess.call(["mv","analyzeData_"+str(i)+".cc" ,"analyzeData.cc"])
   subprocess.call(["make"])
   subprocess.call(["./analyzeData"])
