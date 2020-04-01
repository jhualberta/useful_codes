#/usr/bin/python
import sys,os
path=os.getcwd()
file_list = os.listdir(path)
fList = open('fList_water2018Nov.txt')
for j in fList:
    if "root" in j:
        print j,
