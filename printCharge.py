#/usr/bin/python
# create q-ratio tabler for partialEnergy

import sys,os
path=os.getcwd()
file_list = os.listdir(path)
fList = open('in_x2500_y0_z.dat')
z = []
q = []
for line in fList:
    if line.__contains__('root'):
       ind1 = line.index('_')
       ind2 = line.index('_',ind1+1)
       ind3 = line.index('_',ind2+1)
       ind4 = line.index('_',ind3+1)
       ind5 = line.index('_',ind4+1)
       ind6 = line.index('_',ind5+1)
       z.append(int(line[line.index('z')+1:ind6]))
       #print line[line.index('z'):ind6],
    else:
       #print line.split(' ')[2]
       q.append(line.split(' ')[2])

result = zip(z,q)
result.sort(key=lambda x: x[0])
for i in result:
    print i[0],i[1],
