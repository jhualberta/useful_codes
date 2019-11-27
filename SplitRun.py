#/usr/bin/python
import sys,os
import collections 
path=os.getcwd()
file_list = os.listdir(path)
fList = open('mpwList.dat')
processfile = []
count = 0
runList = []
rootfileList = []

for line in fList:
   processfile.append(line.rstrip())
   count = count+1
   indexBegin = line.index('r00')
   #print indexBegin,line[indexBegin],line[indexBegin+5],line[indexBegin+5:indexBegin+11]
   run = line[indexBegin+5:indexBegin+11] 
   runList.append(run)
   rootfileList.append(line)

counter = collections.Counter(runList)
#sort counter by run ID
sort_counter = sorted(counter.items())
mainrun = [sort_counter[i][0] for i in range(0,len(sort_counter))]
repeat = [sort_counter[i][1] for i in range(0,len(sort_counter))] 
print counter
countSplit = len(mainrun)

mainrun.sort()
for i in mainrun:
    print i

newfileList = []
countNewFile = 0 
iterflag = 0

for val in repeat:
  print val
  f = open("processList"+str(countNewFile+1)+".txt",'r+')
  for i in range(iterflag, iterflag+val):
     f.write(rootfileList[i])
  iterflag = iterflag + val 
  countNewFile = countNewFile+1
  f.close() 
if countSplit != countNewFile:
   print "check data, something wrong!"
print str(count)+" files in total"
print "splitted files "+str(countNewFile)
