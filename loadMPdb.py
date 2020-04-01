import os, sys, string, rat
from ROOT import RAT
db = RAT.DB.Get()
scintMaterialName = "labppo_0p5_oxford"
db.LoadAll(os.environ["GLG4DATA"], "FIT_MULTIPATH.ratdb")
link = db.GetLink("FIT_MULTIPATH", scintMaterialName)
try:
  a = link.GetD("water_RI_tuned")
  print a 
except:
    print "Failed to find the water_RI for MPW"

