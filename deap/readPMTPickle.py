import pandas as pd
from numpy import *
from ROOT import *
## for rAV = 900, pmtID = 30
## 357.48 529.98 633.50

rAV = 900 #851
obj = pd.read_pickle(r'savePMTinfo.pickle')
# Get list of 'Course' values and find index of target value
#result = obj['95']
print "radius=",rAV

for item in obj:
    pmtid = item.keys()[0]
    #if pmtid<30 or pmtid>39: continue
    pmtPhi = item[pmtid][0]
    pmtCostheta = item[pmtid][1]
    vector = TVector3()
    vector.SetMagThetaPhi(rAV, arccos(pmtCostheta), pmtPhi )
    x = vector.X() 
    y = vector.Y() 
    z = vector.Z()
    #if x<300 or z<600: continue
    print pmtid, '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(x, y, z, pmtPhi, pmtCostheta)
