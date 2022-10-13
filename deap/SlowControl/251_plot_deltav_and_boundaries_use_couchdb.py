#!/usr/bin/python
import os,sys
import couchdb
import ROOT
import time
from datetime import datetime

rn_min, rn_max = 15390, 16267
#rn_min, rn_max = 15390, 15440
deltav_txt_file = 'GArdata.txt'
ymax = 30.

## ___ useful funcs ___
def ODBtime_to_sec(s):
    return int(time.mktime(time.strptime(s,'%a %b %d %H:%M:%S %Y')))

## ___ init couchdb ___
serv_url = 'http://ddaq3.snolab.ca:5984'
dbna = 'deapdb'
serv = couchdb.Server(serv_url)
db = serv[dbna]

if __name__=='__main__':

    ## 2 steps:
    #  1. read in start/end times from horribly slow couchdb
    #  2. read in the deltav info (pressure vs tutc)

    ## ___ 1. read in the start/end times from couch ___
    print 'INFO: reading in start/end times from couchdb'
    print 'INFO:  rn_min=',rn_min
    print 'INFO:  rn_max=',rn_max

    tstarts = {}
    tends = {}

    print 'INFO: retrieving WebView/runinfoConfig...',
    results = db.view('WebView/runinfoConfig', reduce=False)
    print 'DONE'
    for row in results:
        rn = row.key
        if rn<rn_min or rn>rn_max: continue
        print 'INFO:  getting data for',rn

        entry = db[row.id]
        try:
            tstart = entry['dateTimeStart']
            tend = entry['dateTimeEnd']
        except KeyError:
            print '\nWARNING: cannot find dateTimeStart and dateTimeEnd for',rn,', will NOT include this run'
            print 'INFO: entry follows'
            print entry
            tstart = None
            tend = None

        if tstart==None:
            continue

        tstart_int = ODBtime_to_sec(tstart)
        tend_int = ODBtime_to_sec(tend)

        tstarts[rn] = tstart_int
        tends[rn] = tend_int


    goodrns = []
    for rn in range(rn_min,rn_max+1):
        if rn in tstarts.keys():
            goodrns.append(rn)
    goodrns.sort()

    ## ___ 2. read in the deltav info ___
    print 'INFO: reading in deltav info'
#    tutc_press = []

    gpres = ROOT.TGraph()

    f = file(deltav_txt_file, 'r')
    iline=0
    for l in f:
        iline+=1
        if not iline%1000:
            print 'INFO: processing line',iline
        if l.startswith('#'): continue
        ls = l.split()
        if len(l.split())!=2: continue
        if not l.startswith('2016'): continue
        press = float(ls[1])
        dvdate, dvtime = ls[0].split('T')
        dvtime = dvtime.rstrip('Z')
#        print dvdate,dvtime

        t = datetime.strptime(dvdate+' '+dvtime, '%Y-%m-%d %H:%M:%S')
        tepoch = datetime(1970,1,1)
        dt = (t-tepoch).total_seconds()# + fracsec + self.toffset

        for rn in goodrns:
#            print 'rn=',rn
            if dt>=tstarts[rn] and dt<=tends[rn] and tends[rn]-tstarts[rn]!=0.:
                rneff = rn+(dt-tstarts[rn])/(tends[rn]-tstarts[rn])
                gpres.SetPoint(gpres.GetN(), rneff, press)

    gpres.SetMarkerStyle(20)
    gpres.Draw('APL')
    au = raw_input('>')
