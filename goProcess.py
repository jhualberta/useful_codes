#/usr/bin/python
import sys,os
import subprocess
path=os.getcwd()
file_list = os.listdir(path)
fList = open('filelist1.dat')
processfile = []
count = 0
for i in fList:
 processfile.append(i.rstrip())
 count = count+1

for j in range(1,count+1):
 print str(j) + ' ' + processfile[j-1]
 filename = 'analyPartialLeon.C'
 ftemp = open(filename,'r+')
 flist=ftemp.readlines()
 flist[35]= '   const char* filename = \"'+processfile[j-1]+'\";\n'
 ftemp=open(filename,'w+') 
 ftemp.writelines(flist)
 ftemp.close()
 subprocess.call(["root","-q" ,"analyPartialLeon.C+"])

fList.close()
