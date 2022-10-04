from ROOT import *

fileList = ["Merged_30681.root","Merged_30686.root","Merged_30695.root","Merged_30705.root","Merged_30706.root","Merged_30707.root","Merged_30708.root","Merged_30717.root","Merged_30726.root","Merged_30741.root","Merged_30742.root","Merged_30743.root","Merged_30744.root","Merged_30746.root","Merged_30747.root","Merged_30751.root","Merged_30756.root","Merged_30760.root","Merged_30765.root","Merged_30769.root","Merged_30774.root","Merged_30784.root","Merged_30785.root","Merged_30804.root","Merged_30813.root","Merged_30815.root","Merged_30826.root","Merged_30837.root"]

n_double = []
n_single = []

print "total, single, double"
for files in fileList:
    ff = TFile(files)
    hist = ff.Get("nclusters")
    number = hist.GetEntries()
    single = hist.GetBinContent(2)
    double = number - single
    #print number, single, number-single
    print number single, double
    print "--"
    n_double.append(double)
    n_single.append(single)
    ff.Close()

print "double: ", n_double
print "single: ", n_single
