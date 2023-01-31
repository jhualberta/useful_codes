import sys,os
import subprocess

path=os.getcwd()
file_list = os.listdir(path)
fgetfile = open('done.dat','r+')
fList = fgetfile.readlines()
print fList
count = 0
fileZombie = []
fileNoTree = []
for line in fList:
     fname = line.rstrip() ## remove '\n' at the end of filename!!!
     ### print "process", fname
     if count%100 == 0:
         print "has processed", count, "files."
     ff = TFile(fname)
     if ff.IsZombie():
          print "zombie file", fname
          fileZombie.append(fname)
     if not(ff.GetListOfKeys().Contains("T")):
          print "file has no TTree!!!", fname
          fileNoTree.append(fname)
     count += 1

print "these files are zombies!!", fileZombie
print "thest files have no tree!!", fileNoTree

