#/usr/bin/python
# print out E resolution, E scale from fitConvol output txt files
import sys, os
import subprocess
path=os.getcwd()
file_list = os.listdir(path)
fList = open('listOutput.dat')
processfile = []
count = 0
result = []
for i in fList:
 ftemp = open(i.rstrip(),'r+')
 flist=ftemp.readlines()
 for data in flist:
    result = [data.split(' ')[4],data.split(' ')[5],data.split(' ')[8],data.split(' ')[9]]
 processfile.append(result) 
 count = count+1
#
print "Esigma, EsigmaErr, Escale, EscaleErr"
for data in processfile:
    print data[2],data[3],data[0],data[1]
