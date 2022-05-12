from ROOT import TFile, TTree, gRandom
from array import array
 
file = TFile("tree.root", 'recreate')
tree = TTree("tree_name", "tree title")
 
Nmax = 10
nParticles = array( 'i', [ 0 ] )
pt = array( 'd', Nmax*[ 0. ] )
tree.Branch( 'nParticles', nParticles, 'nParticles/I' )
tree.Branch( 'pt', pt, 'pt[nParticles]/D' )
 
for i in xrange(1000): # loop over events
    nParticles[0] = int(gRandom.Uniform()*10)
    for j in xrange(nParticles[0]): # loop over particles in this event
       pt[j] = gRandom.Gaus(20,2)
    tree.Fill()
 
file.Write("",TFile.kOverwrite)
file.Close()
