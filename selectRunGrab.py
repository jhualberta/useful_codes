#/usr/bin/python
import sys,os
path=os.getcwd()
file_list = os.listdir(path)
fList = open('filelist.dat')
processfile = []
count = 0
fCheck = open('fileList_oct2018.txt')
name = []
for i in fCheck:
    f = i.split('.')
    name.append(f[0])

for i in fList:
  fname = i.split('.')
  for j in name:
    if fname[0] == j:
       print i,
       break

# print i.rstrip()
