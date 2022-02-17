import numpy as np
from numpy import log, sqrt
# zhihua zhou, watermelon book, P100
def ent(p):
  if p != 0:
    ent = p*np.log(p)/np.log(2)
  else:
    ent = 0
  return ent

densityData = [0.697, 0.774, 0.634, 0.608, 0.556, 0.403, 0.481, 0.437, 0.666, 0.243, 0.245, 0.343, 0.639, 0.657, 0.36, 0.593, 0.719]
output = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]
count = sum(output)
entTot = -(ent(count/17.) + ent(1.0 - count/17.))

rawpairs = zip(densityData, output)
rawpairs.sort(key=lambda x: x[0])
print "# of data: ", len(densityData)
N = len(densityData)
T = []
count = 0
for i in range(N-1):
    ta = (rawpairs[i][0]+rawpairs[i+1][0])/2
    tpair = (ta, rawpairs[i][1])
    T.append(tpair)
    if rawpairs[i][1]>0:
        count += 1

pt = []
Dv = []
for i in rawpairs:
  count0 = 0
  dv = 0
  for j in T:
      if i[0]>j[0]:
          count0 += i[1]
          dv += 1
  p = count0/17.
  ee = -(ent(p)+ent(1-p))
  pt.append(ee)
  Dv.append(dv)

sum0 = 0
gain = []
for i in range(17):
  sum0 += Dv[i]/17.*pt[i]
  gain.append(entTot-sum0) 
#print sum0
#print pt
print max(gain)
