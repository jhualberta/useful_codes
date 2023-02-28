#/usr/bin/python
import sys,os
import shutil
from datetime import datetime
from time import strptime
import pytz

path=os.getcwd()
file_list = os.listdir(path)
fList = open('u232_run_dates.dat')
processfile = []
count = 0
run = []
date = []
date_str = []
utc = []
kk = 0
for i in fList:
  processfile.append(i.rstrip())
  if kk%2 == 0:
      x1 = int(i.split('\t')[0])
      #print xx
      run.append(x1)
  elif kk%2 != 0:
      x2 = i.split('\t')[2]
      ss = x2.split(' ')[1:5]
      ### ss = ['Jul', '13', '10:53:05']
      mon = ss[0]
      mm = strptime(mon,'%b').tm_mon
      dd = ss[1]
      time = ss[2]
      yy = ss[3]
      recast = yy+'-'+str(mm)+'-'+dd+' '+time
      date_str.append(recast)
      naive = datetime.strptime(recast, "%Y-%m-%d %H:%M:%S")
      #timestamp = time.mktime(naive.timetuple()) + naive.microsecond/1e6
      #print utc
      date.append(naive)
  kk = kk+1

data = zip(run,date_str)

data = sorted(data, key = lambda x: x[0])
print data
print len(data)
print sorted(run)

realRun = [28103, 28105, 28110, 28119, 28123, 28132, 28137, 28147, 28160, 28164, 28173, 28188, 28192, 28193, 28203, 28208, 28218, 28222, 28225, 28234, 28239, 28248, 28252, 28261, 28266, 28275, 28279, 28288, 28294, 28295, 28310, 28314, 28323, 28325, 28330, 28331, 28343, 28348, 28357, 28362, 28363, 28372, 28374, 28375, 28379, 28381, 28383, 28393, 28398, 28409, 28413, 28422, 28427, 28436, 28440, 28449, 28454, 28464, 28468, 28472, 28481, 28486, 28487, 28512, 28516, 28525, 28530, 28545, 28549, 28565, 28570, 28579, 28583, 28593, 28595, 28608, 28624, 28632, 28637, 28655, 28659, 28660, 28670, 28676, 28685, 28689, 28702, 28708, 28730, 28734, 28743, 28748, 28758, 28763, 28786, 28795, 28800, 28809, 28813, 28815, 28824, 28829, 28838, 28844, 28847, 28849]
print set(run)-set(realRun)

for x in data:
    if x[0]!=28178 and x[0]!=28839:
       print x[0],x[1]










