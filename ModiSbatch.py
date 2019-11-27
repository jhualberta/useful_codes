#/usr/bin/python
import sys,os
path=os.getcwd()
file_list = os.listdir(path)
count = 60 

for j in range(1,count+1):
 filename = 'sbatch_'+str(j)+'.sh'
 ftemp = open(filename,'r+')
 flist=ftemp.readlines()
 flist[7]='rat FitMulti'+str(j)+'.mac\n'
 ftemp=open(filename,'w+') 
 ftemp.writelines(flist)
 ftemp.close()
