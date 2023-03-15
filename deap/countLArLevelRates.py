''' This comes with RAT, so source rat.env first! '''
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

##21604, 75803 seconds = 0.877349 days
##18831, 103383 seconds = 1.19657 days
##22387, 73903.6 seconds = 0.855366 days
##22394, 74506.8 seconds = 0.862347 days

  TFile *fLar = new TFile("Combined_Sina_Lar_u232_run21604.root");
  TFile *f18831 = new TFile("Combined_Sina_Lar_u232_run18831.root");
  TFile *f22387 = new TFile("Combined_Sina_Lar_u232_run22387.root");
  TFile *f22394 = new TFile("Combined_Sina_Lar_u232_run22394.root");

  double time0 = 75803;
  double time1 = 103383;
  double time2 = 73903.6;
  double time3 = 74506.8; 
  





