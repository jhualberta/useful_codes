#/usr/bin/python

# search for ".mac" in "run file"
# modify N16 macros, add the threshold, number of events

import sys,os
from decimal import Decimal
import shutil

path=os.getcwd() #path name string
newpath = 'newMacros'
#os.mkdir(newpath)

a=[]
#Need format like: -5861.1, -2524.1, -5000 
#Be careful to the whitespaces

f1=open('position.txt','r')

for line in f1.read().split('\n'):
  if line.strip(): #check line not empty
    #print line
    a=line.split(',')
    xValue = Decimal(a[0])
    yValue = Decimal(a[1])
    zValue = Decimal(a[2])
    xInt = int(xValue)
    yInt = int(yValue)
    zInt = int(zValue)
    TemplateFilename = 'template.mac'
    shutil.copy2(TemplateFilename, path+'/newname.mac') # copy file
    newMacroName = 'N16Source_x'+str(xInt)+'_y'+str(yInt)+'_z'+str(zInt)+'.mac'
    print 'Creating new macro: '+newMacroName
    os.rename("newname.mac", newMacroName)
    f2=open(newMacroName,'r+') # get Macro Template
    flist=f2.readlines()
    flist[0] = '# File: '+newMacroName+'\n'
    flist[4] = '#    X = '+str(xValue)+' mm'+'\n'
    flist[5] = '#    Y = '+str(yValue)+' mm'+'\n'
    flist[6] = '#    Z = '+str(zValue)+' mm'+'\n'
    flist[18] = '/rat/db/set GEO[N16Source] position'+' [ '+str(xValue)+', '+str(yValue)+', '+str(zValue)+' ]'+'\n'
    flist[32] = '/rat/procset file '+'\"SNOPMC_N16Source_x'+str(xInt)+'_y'+str(yInt)+'_z'+str(zInt)+'_1e4evts.root\"'	
    f2=open(newMacroName,'w+') # this line is required, otherwise all the file contents will be attached
    f2.writelines(flist)
    f2.close()
    #shutil.copy2(TemplateFilename, path+'/'+newpath) 
f1.close()

