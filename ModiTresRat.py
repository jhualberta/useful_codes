#/usr/bin/python
import sys,os
import shutil

path=os.getcwd()
file_list = os.listdir(path)
fList = open('all.dat')
processfile = []
count = 0
for i in fList:
  processfile.append(i.rstrip())
# print i.rstrip()
  count = count+1

for i in range(count):
  shutil.copyfile("AnalyWaterTres_allRat6176.bb", "AnalyWaterTres_allRat6176_"+str(i)+".C")

for i in  range(count):
 print 'write', processfile[i]
 filename = 'AnalyWaterTres_allRat6176_'+str(i)+'.C'
 ftemp = open(filename,'r+')
 flist=ftemp.readlines()
 flist[25]='void AnalyWaterTres_allRat6176_'+str(i)+'()\n'
 flist[27]=' const char* filename = \"'+processfile[i]+'\";\n'

 ftemp=open(filename,'w+') 
 ftemp.writelines(flist)
 ftemp.close()
