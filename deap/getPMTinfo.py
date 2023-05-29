import os
import couchdb
import sys
import glob
from rat import *
#from rat import PMTInfoUtil
from ROOT import *
import numpy as np
from numpy import *
import math
from datetime import datetime, time as datetime_time, timedelta
from array import array
global qc
PI = np.pi
## print "will read from files %s and write to file %s\n"%(sys.argv[1],sys.argv[2])
server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]
#Dir=sys.argv[1]
file0 = str(sys.argv[0])
fileName = os.path.basename(file0)


server = couchdb.Server("https://deimos.physics.carleton.ca:6984/")
db = server["deapdb"]

results = db.view('WebView/PMTPos', include_docs=True)
Phi = np.zeros(255)
CosTheta = np.zeros(255)
index = np.zeros(255,dtype = int)

pmtX = np.zeros(255)

## buil PMT look-up tables
for row in results:
    if row.doc["run_range"][0] == 0:
         #print(int(row.key), row.doc["locationPhi"])
         CosTheta[int(row.key)] = np.cos(row.doc["locationTheta"]*PI/180)
         Phi[int(row.key)] = row.doc["locationPhi"]*PI/180
         index[int(row.key)] = int(row.key)

pmtPos = []  ## Cartesian coordination
for x in zip(CosTheta, Phi):
    pmtPos.append(TVector3(np.sqrt(1.0-x[0]**2)*np.cos(x[1]),np.sqrt(1.0-x[0]**2)*np.sin(x[1]),x[0]))

for pmt in pmtPos:
    pmt.SetMag(900)
    if pmt.Z()>551:
        print pmt.X(), pmt.Y(), pmt.Z()






