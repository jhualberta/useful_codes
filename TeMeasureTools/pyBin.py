#!/usr/bin/env
# read binary
import ROOT
import sys, os, getopt
from ROOT import TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array
import struct

path = os.getcwd()
file_list = os.listdir(path)
inputfile = "ch0_LU151_UA-8_bin_attenu_Coin_16July.bin"
'''
format:
for each event:
board channel timeStamp Elong Eshort flag nSample waveform
2       8         2     2     2      4    252
'''
with open(inputfile, "rb") as file1:
  byte = file1.read(1)
  while byte != "":
    # Do stuff with byte.
    byte = file1.read(1)
    struct.unpack("i" * ((len(byte) -24) // 4), byte[20:-4])


#nsample = np.uint(0xfc) # take 252 (0xfc) samples of waveform
#np.uint(0xfc)














