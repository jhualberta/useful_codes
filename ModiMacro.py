#/usr/bin/python
import sys,os
path=os.getcwd()
file_list = os.listdir(path)
fList = open('filelist.dat')
processfile = []
count = 0
for i in fList:
 processfile.append(i.rstrip())
# print i.rstrip()
 count = count+1

start = 61
for j in range(start,start+count):
 print str(j) + ' ' + processfile[j-start]
 filename = 'FitMulti'+str(j)+'.mac'
 ftemp = open(filename,'r+')
 flist=ftemp.readlines()
 flist[0]='/rat/inroot/load '+processfile[j-start]+'\n'
 flist[11]='/rat/procset file \"FitMultiWater_'+processfile[j-start]+'\"\n'
 flist[13]='/rat/inroot/read '+processfile[j-start]+'\n'
 ftemp=open(filename,'w+') 
 ftemp.writelines(flist)
 ftemp.close()
