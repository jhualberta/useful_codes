#!/usr/bin/env python
'''
A python tool for evaluating the performance
of SNO+ fitter algorithms.

Author : Ed Leming <edward.leming@physics.ox.ac.uk>
'''
import ROOT 
# Stop ROOT stealing our useage messages
ROOT.PyConfig.IgnoreCommandLineOptions = True
import rat
import json
import time
import numpy as np

config = {}

class DimensionTool(object):
    '''
    Class to retreive some relavent parameters for the different possible dimensions
    '''
    def __init__(self, dimension='udotr_vs_rcube'):
        self.set_dimension(dimension)
        
    def set_dimension(self, dimension):
        if not dimension.lower() in ['udotr_vs_rcube', 'e_vs_rcube', 'ebias_vs_rcube',  "ebias_vs_efit", "rbias_vs_efit", "rbias_vs_rfit"]:
            raise Exception("AxisConverter: {0} is not a recognised dimension".format(dimension))
        self._dimension = dimension.lower()

    def get_axis_bounds_2d(self):
        if self._dimension in ['udotr_vs_rcube']:
            return config["rcube_min"], config["rcube_max"], config["udotr_min"], config["udotr_max"]
        elif self._dimension in ['e_vs_rcube']:
            return config["rcube_min"], config["rcube_max"], config["e_min"], config["e_max"]
        elif self._dimension in ['ebias_vs_rcube']:
            return config["rcube_min"], config["rcube_max"], config["ebias_min"], config["ebias_max"]
        elif self._dimension in ['ebias_vs_efit']:
            return config["e_min"], config["e_max"], config["ebias_min"], config["ebias_max"]
        elif self._dimension in ['rbias_vs_efit']:
            return config["e_min"], config["e_max"], -1*config["r_max"], config["r_max"]
        elif self._dimension in ['rbias_vs_rfit']:
            return config["r_min"], config["r_max"], -1*config["r_max"], config["r_max"]

    def get_histo_bounds_2d(self):
        if self._dimension in ['udotr_vs_rcube']:
            return config["rcube_min"], config["rcube_max"], config["udotr_min"], config["udotr_max"]
        elif self._dimension in ['e_vs_rcube']:
            return config["rcube_min"], config["rcube_max"], config["e_min"], config["e_max"]
        elif self._dimension in ['ebias_vs_rcube']:
            return config["rcube_min"], config["rcube_max"], config["ebias_min"], config["ebias_max"]
        elif self._dimension in ['ebias_vs_efit']:
            return config["e_min"], config["e_max"], config["ebias_min"], config["ebias_max"]
        elif self._dimension in ['rbias_vs_efit']:
            return config["e_min"], config["e_max"], -1*config["r_max"], config["r_max"]
        elif self._dimension in ['rbias_vs_rfit']:
            return config["r_min"], config["r_max"], -1*config["r_max"], config["r_max"]

    def get_nbins_2d(self):
        if self._dimension in ['udotr_vs_rcube']:
            return config["rcube_nbins"], config["udotr_nbins"]
        elif self._dimension in ['e_vs_rcube']:
            return config["rcube_nbins"], config["e_nbins"]
        elif self._dimension in ['ebias_vs_rcube']:
            return config["rcube_nbins"], config["ebias_nbins"]
        elif self._dimension in ['ebias_vs_efit']:
            return config["e_nbins"], config["ebias_nbins"]
        elif self._dimension in ['rbias_vs_efit']:
            return config["e_nbins"], config["r_nbins"]
        elif self._dimension in ['rbias_vs_rfit']:
            return config["r_nbins"], config["r_nbins"]

    def get_bin_width_2d(self):
        axisXMin, axisXMax, axisYMin, axisYMax = self.get_axis_bounds_2d()
        nbinsX, nbinsY = self.get_nbins_2d()
        return (axisXMax - axisXMin) / nbinsX, (axisYMax - axisYMin) / nbinsY

    def get_axis_title(self):
        if self._dimension in ['udotr_vs_rcube']:
            return ("R^{3}/R^{3}_{AV}", "u #cdot R")
        elif self._dimension.split('_')[0] in ['ebias']:
            return ("R^{3}/R^{3}_{AV}", "(E_{fit}-E_{MC})/E_{MC} [*100% MeV]")
        elif self._dimension.split('_')[0] in ['e']:
            return ("R^{3}/R^{3}_{AV}", "E_{fit} [% MeV]")

    def get_truth(self, mc):
        if self._dimension[0] in ['e']:
            return mc.GetMCParticle(0).GetKineticEnergy()
        elif self._dimension[0] in ['u']:
            return mc.GetMCParticle(0).GetPosition().Mag()
        elif self._dimension[0] in ['r']:
            return mc.GetMCParticle(0).GetPosition()
        else:
            return 0

    def get_2dplot(self, fitVertex, mc):
        if self._dimension in ['udotr_vs_rcube']:
            return ( pow( (fitVertex.GetPosition().Mag()/6005.0), 3 ), (fitVertex.GetPosition().Unit()*fitVertex.GetDirection().Unit()) )
        elif self._dimension in ['e_vs_rcube']:
            return ( pow( (fitVertex.GetPosition().Mag()/6005.0), 3 ), fitVertex.GetEnergy() )
        elif self._dimension in ['ebias_vs_rcube']:
            return ( pow( (fitVertex.GetPosition().Mag()/6005.0), 3 ), ((fitVertex.GetEnergy() - self.get_truth(mc)) / self.get_truth(mc))*100 )
        elif self._dimension in ['ebias_vs_efit']:
            return ( fitVertex.GetEnergy(), ((fitVertex.GetEnergy() - self.get_truth(mc)) / self.get_truth(mc))*100 ) 
        elif self._dimension in ['rbias_vs_efit']:
            return ( fitVertex.GetEnergy(), (fitVertex.GetPosition() - self.get_truth(mc)) * self.get_truth(mc).Unit() )
        elif self._dimension in ['rbias_vs_rfit']:
            return ( fitVertex.GetPosition().Mag(), (fitVertex.GetPosition() - self.get_truth(mc)) * self.get_truth(mc).Unit() )
        else:
            return (0,0)

def plot_2d(fname, dimension, outfile, fitter="waterFitter", nevents=1e5, verbose=False):
    '''
    Plot the 2D fitter results
    The performance will be evaluated for variable defined in plotting.cfg,
    all of which will be plotted as a function of 'dimension'.
    '''
    # Make a primary dimension tool to handle requests wrt the passed dimension
    primaryDimTool = DimensionTool(dimension)

    minBinX, maxBinX, minBinY, maxBinY = primaryDimTool.get_histo_bounds_2d()
    nbinX, nbinY = primaryDimTool.get_nbins_2d()
    nbinX = int(nbinX)
    nbinY = int(nbinY)

    ######################
    # Book some histograms
    # Loop over each requested dimension
    hist2d = {}
    secondaryDimTool = DimensionTool(dimension) # Make a secondary tool to handle conversions wrt internal dimensions

    for dim in config["plot2d"]:
        secondaryDimTool.set_dimension(dim)
        if dim in ['udotr_vs_rcube','e_vs_rcube', 'ebias_vs_rcube', 'ebias_vs_efit', 'rbias_vs_efit','rbias_vs_rfit']:
           hist2d[dim] = ROOT.TH2D("{0}".format(dim),
                                          "", nbinX, minBinX, maxBinX, nbinY, minBinY, maxBinY)
           if dim.split('_')[2] in ['rcube']:           
              hist2d[dim].GetXaxis().SetTitle("R^{3}/R_{AV}^{3}")
           elif dim.split('_')[2] in ['efit']:
              hist2d[dim].GetXaxis().SetTitle("e_{fit} [MeV]")
           elif dim.split('_')[2] in ['rfit']:
              hist2d[dim].GetXaxis().SetTitle("r_{fit} [mm]")
           if dim[0] in ['e']:
              if dim.split('_')[0] in ['ebias']:
                 hist2d[dim].GetYaxis().SetTitle(dim.split('_')[0]+"=(e_{fit}-e_{mc})/e_{mc} *100%")
              else:
                 hist2d[dim].GetYaxis().SetTitle(dim.split('_')[0]+" [MeV]")
           elif dim[0] in ['r']:
              if dim.split('_')[0] in ['rbias']:
                 hist2d[dim].GetYaxis().SetTitle(dim.split('_')[0]+"=(R_{fit}-R_{mc})*R_{mc}/|R_{mc}|")
           else:    
              hist2d[dim].GetYaxis().SetTitle(dim.split('_')[0])

    # Loop over the MC events
    simCounter, evCounter, fittedCounter = 0, 0, 0
    for ds,run in rat.dsreader(fname):
        if simCounter == 0:
            loopStart=time.time()
            print "Begining event loop..."
        simCounter = simCounter + 1
        for iev in range(0, ds.GetEVCount()):
            # Use retriggers?
            if config["use_retriggers"] == 0 and iev > 0:
                continue 

            # Increment counter
            evCounter = evCounter + 1

            # Get DS variables
            ev = ds.GetEV(iev)
            mc = ds.GetMC()

            # Get truth and add to nEvents histo
            truth = primaryDimTool.get_truth(mc)

            # Get fitter vertex
            try:
                fitResult = ev.GetFitResult(fitter)
                fitVertex = fitResult.GetVertex(0)
            except Exception as e:
                if verbose:
                    print "Simulated event {0:d}, GTID {1:d} : {2}".format(simCounter,
                                                                           ev.GetGTID(),
                                                                           e)
                continue
            #fitExists.Fill(truth)
            if (not fitVertex.ContainsPosition()) or (not fitVertex.ValidPosition()):
                continue
            if (not fitVertex.ContainsTime()) or (not fitVertex.ValidTime()):
                continue
            if (not fitVertex.ContainsEnergy()) or (not fitVertex.ValidEnergy()):
                continue
            if (not fitVertex.ContainsDirection()) or (not fitVertex.ValidDirection()):
                continue

            # Increment counter
            fittedCounter = fittedCounter + 1
            for dim in hist2d:
                secondaryDimTool.set_dimension(dim)
                try:
                  hist2d[dim].Fill( secondaryDimTool.get_2dplot(fitVertex, mc)[0], secondaryDimTool.get_2dplot(fitVertex, mc)[1] )
                except IndexError:
                    continue 

            # Print 
            if fittedCounter > 0 and fittedCounter % 5000 == 0:
                print "{0} events processed, [{1:.1f} s]".format(fittedCounter, time.time()-loopStart)
                loopStart=time.time()
            # If we're above threshold, break
            if fittedCounter >= nevents:
                break
        if fittedCounter >= nevents:
            break

    print "######################################"
    print "{0} simulated decays".format(simCounter)
    print "{0} triggered events".format(evCounter)
    print "{0} fitted events".format(fittedCounter)

    #######################################
    # Make summary plots and save to file
    outFile = ROOT.TFile(outfile, "RECREATE")

    # Fit all other build histograms and make plots
    can = ROOT.TCanvas("c1","c1")

    for dim in sorted(hist2d):
       axisXMin, axisXMax, axisYMin, axisYMax = secondaryDimTool.get_axis_bounds_2d()
       axisXBinWidth, axisYBinWidth = secondaryDimTool.get_bin_width_2d()
       nbins0, nbins1 = secondaryDimTool.get_nbins_2d()
       h2d = hist2d[dim]
       h2d.Write()
       # Format plots and write

    print "Results written to: {0}".format(outfile)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in. Will accept wildcard")
    parser.add_argument('dimension', type=str,
                        help="Variable to make 1D plots wrt [x, y, z, r, e]")
    parser.add_argument('-n', '--nevents', type=int, default=100000,
                        help="Maximum number of fitted events to evaluate from infile(s) [1e5]")
    parser.add_argument('-f', '--fitter', type=str, default="scintFitter",
                        help="Name of the fitter processor to be evaluated [scintFitter]")
    parser.add_argument('-o', '--outfile', type=str,
                        help="Path to outfile [./2d_plots_{dimension}.root]")
    parser.add_argument('-c', '--config', type=str, default="./plotting.cfg",
                        help="Path to plotting config file. If black wil use ./plotting.cfg")
    args = parser.parse_args()
    start = time.time()

    # ROOT stuff
    ROOT.gROOT.SetBatch(True)

    # Read in plotting config
    execfile(args.config, {}, config)

    outfile = args.outfile
    if not args.outfile:
        outfile = "./2d_plots_{0}.root".format(args.dimension)

    # Plot performance vs the requested dimension
    plot_2d(args.infile, args.dimension, outfile,
                                  fitter=args.fitter,
                                  nevents=args.nevents)
    print "######################################"
    print "Full script took {0:.1f} mins".format((time.time() - start)/60.)
