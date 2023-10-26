#Jie: search rootfiles for some runs in deap data
import sys,os
import subprocess

## put the run numbers want to check
runList = [28103, 28105, 28110]

all_ratVer = ['v5.9.6','v5.10.1','v5.10.3','v5.10.4','v5.10.5','v5.10.6','v5.10.8','v5.11.0','v5.11.2','v5.12.1','v5.12.2','v5.12.3','v5.12.6','v5.12.8','v5.13.1','v5.13.2','v5.14.0','v5.14.2','v5.14.2_noPLSC','v5.14.3_SCalFix','v5.14.3_SCalFix2','v5.15.0','v5.15.2','v5.15.3']

print "specify a rat version since 5.13?"
index0 = all_ratVer.index('v5.13.1') 
version_rat = all_ratVer[index0:]

fileType = 'ntp'
print "cal or ntp? checking: ", fileType 

fileList = []
get_ratVer = []
for run in runList:
    for ratVer in version_rat:
       try:        
           fileList.append(subprocess.check_output(["ls", "-1","/project/6004969/data/"+ratVer+"/"+fileType+"/run0"+str(run)]))
           print "!!! this run", run, "exists in ", ratVer
           get_ratVer.append(ratVer)
       except:
           print "this run", run, "not exist in ", ratVer
           continue  

fileList = [i.split('\n') for i in fileList]

exist_ratVer = list(set(get_ratVer))
print "exit rat version", exist_ratVer
print "access all the rootfiles you want"
for item in fileList:
    for nn in item:
         if nn == '': continue
         for ratVer in exist_ratVer: 
            print "/project/6004969/data/"+ratVer+"/"+fileType+"/run"+nn[nn.index('0'):nn.index('0')+6]+"/"+nn
