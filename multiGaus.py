import ROOT
import numpy as np
import time
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
import sys
#import fitmrq
from ROOT import TH2F, TH1F, TRandom


# create fake data
npoints = 40 # 49 data points
nfit = 9 # fit 3*3 of parameters
tolerance = 0.1

# for fitter converge
max_iter = 1000
NGOOD = 4 # if find 4 good fits, quit

xdata = []
for i in range(85,125):
   xdata.append(i)

ydata = [1.913521, 1.953769, 2.347435, 2.883654, 3.493567,4.047560, 4.337210, 4.364347, 4.563004, 5.054247,5.194183, 5.380521, 5.303213, 5.384578, 5.563983, 5.728500, 5.685752, 5.080029, 4.251809, 3.372246, 2.207432, 1.227541, 0.8597788,0.8220503,0.8046592, 0.7684097,0.7469761,0.8019787,0.8362375,0.8744895, 0.9143721,0.9462768,0.9285364,0.8954604,0.8410891, 0.7853871,0.7100883,0.6938808,0.7363682,0.7032954]#, 0.6029015, 0.5600163, 0.7477068, 1.188785, 1.938228, 2.602717, 3.472962, 4.465014, 5.177035]

# check data
hdata = TH1F("h","Example of several fits in subranges",npoints,85,134);
#for i in range(npoints):
#  hdata.SetBinContent(i+1,ydata[i])
#hdata.Draw()
#raw_input("enter to quit")

# Model function
K = 3 # 3 gaussian model
# y(x) = B0_0*exp(-((x-E_0)/G_0)*((x-E_0)/G_0))+ ...
#B_k

def fgauss(xx, pars): #3*3 = 9 pars
  dydpars = [ 0 for k in range(nfit)]
  par0 = pars[0:3]
  par1 = pars[3:6]
  par2 = pars[6:9]
  y = 0.0
  for i in range(nfit/3):
    arg = (xx-par1[i])/par2[i]
    ex = np.exp(-arg*arg)
    fac = par0[i]*ex*2*arg # derivative core
    y += par0[i]*ex # ymodel, calculate xdata[i] in model
    dydpars[3*i] = ex # B0, B1, B2
    dydpars[3*i+1] = fac/par2[i] # E0, E1, E2
    dydpars[3*i+2] = fac*arg/par2[i] # G0, G1, G2
  #print dydpars
  results = [y, dydpars]
  return results

# calculate chi2, beta 
# calculate alpha use (15.5.11)
def mrqcof(pars_trial):
  # initializing
  alpha = [ [0 for i in range(nfit)] for i in range(nfit) ]
  beta = [ 0 for i in range(nfit) ]
  chi2 = 0
  # loop the data
  for i in range(npoints):
    results = fgauss(xdata[i],pars_trial)
    ymodel = results[0]
    dydpars = results[1]
    dy = ydata[i] - ymodel
    sig2i = 1 # set sigma = 1 for convenience, otherwise should come from data
    j = 0
    for l in range(nfit):
       #if par_isfine[l]: # reserved for free and hold
       wt = dydpars[l]*sig2i
       k = 0
       for m in range(l+1):
         alpha[j][k] += alpha[j][k]+wt*dydpars[m]
         k = k+1
       beta[j] += dy*wt
       j = j+1
    chi2 += dy*dy*sig2i # for data loop
  for j in range(nfit):
    for k in range(j):
      alpha[k][j] = alpha[j][k]
  results = [chi2, alpha, beta]
  return results

# fitmrq mpars_start = [1., 1. , 1., 100., 100., 100., 20.,  20.,  20.]

alambda = 0.001
oneda = [ [1] for k in range(nfit) ] # nfit*1 vector
# trial parameters
alpha = [ [0 for i in range(nfit)] for i in range(nfit) ]
beta = [0 for kk in range(nfit)]
pars_start = [1., 1. , 1., 100., 100., 100., 20.,  20.,  20.]
#[5., 100., 2., 5., 100., 2.,5., 100., 2.] 

atry = pars_start
a = atry
results_start = mrqcof(pars_start)
chisqOld = results_start[0]
alpha = results_start[1]
beta = results_start[2]
cov = [[0 for i in range(nfit)] for kk in range(nfit)]
temp = [[0 for i in range(nfit)] for kk in range(nfit)]
da =  [0 for kk in range(nfit)]
nGoodfits = 0
for iter in range(max_iter):
    if nGoodfits == NGOOD:
       alambda = 0.0
       print "reach goodfits: ", nGoodfits
       break
    for j in range(nfit):
       for k in range(nfit):
          cov[j][k] = alpha[j][k]
       cov[j][j] = alpha[j][j]*(1.0+alambda)
       for k in range(nfit):
          temp[j][k] = cov[j][k]
       oneda[j][0] = beta[j] # delta a
    # solve linear equations: alpha'*da = beta
    # ie: cov*da = oneda, where temp = cov
    # do QR decomposition of temp to solve, eg: Ax = B, QRx = B, Q^TQ = I, Rx = Q^TB, y = Q^TB, Rx = y
    Acov = scipy.array(temp)
    Boneda = scipy.array(oneda)
    try:
      Qcov, Rcov = scipy.linalg.qr(Acov)
    except ValueError:
      print "can not solve da from cov matrix, break" 
      break
    solveX = scipy.linalg.solve(Rcov, np.dot(Qcov.T, Boneda))
    print Acov
    #print Qcov, Rcov 
    for j in range(nfit):
      oneda[j] = solveX[j]

    for j in range(nfit):
       for k in range(nfit):
          cov[j][k] = temp[j][k]
       da[j] = oneda[j][0]
    j = 0
    for l in range(nfit):
      atry[l] = a[l] + da[j]
      j = j + 1 
    #print atry, a, da
    results_new = mrqcof(atry)
    chisq = results_new[0]
    cov = results_new[1]
    da = results_new[2] # beta
    if abs(chisq - chisqOld)<max(tolerance, tolerance*chisq):
        nGoodfits += 1
    if chisq<chisqOld:
        alambda *= 0.1
        chisqOld = chisq
        for j in range(nfit):
            for k in range(nfit):
              alpha[j][k] = cov[j][k]
            beta[j] = da[j]    
        for l in range(nfit):
            a[l] = atry[l]
    else:
        alambda *= 10.0
        chisq = chisqOld
        if alambda>1e6:
            print "lambda range out"
            break
    print chisq, chisqOld
    print alambda
print "best fit results: ", a


#   alambda = 0.001
#   mfit = 0
