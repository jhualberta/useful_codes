import ROOT
import numpy as np
import time
import scipy
import scipy.linalg   # SciPy Linear Algebra Library
import sys
#import fitmrq
from ROOT import TH2F, TH1F, TRandom


# create fake data
npoints = 10 #data points
nfit = 9 # fit 3*3 of parameters
tolerance = 0.001

# for fitter converge
max_iter = 1000
NGOOD = 4 # if find 4 good fits, quit

xdata = []
#for i in range(85,125):
#   xdata.append(i)

for i in range(0,10):
   xdata.append(i)
#TF1 *g1    = new TF1("g1","gaus",85,95);
#TF1 *g2    = new TF1("g2","gaus",98,108);
#TF1 *g3    = new TF1("g3","gaus",110,121);
ydata = [2.970149501247504, 2.9757981498167134, 2.980861309137447, 2.9853359562474893, 2.98921941669298, 2.99250936719238, 2.995203837952819, 2.997301214635582, 2.9988002399680034, 2.9997000149995]

# check data
hdata = TH1F("hdata","Example of several fits in subranges", npoints, xdata[0], xdata[npoints-1]);
for i in range(npoints):
  hdata.SetBinContent(i+1,ydata[i])

# Model function
K = 3 # 3 gaussian model
# y(x) = B0_0*exp(-((x-E_0)/G_0)*((x-E_0)/G_0))+ ...

def fgauss(xx, pars): #3*3 = 9 pars
  dyda = [ 0 for k in range(nfit)]
  par0 = pars[0:3] #B0, B1, B2
  par1 = pars[3:6] #E0, E1, E2
  par2 = pars[6:9] #G0, G1, G2
  y = 0.0
  for i in range(nfit/3):
    arg = (xx-par1[i])/par2[i]
    ex = np.exp(-arg*arg)
    fac = par0[i]*ex*2*arg # derivative core
    y += par0[i]*ex # ymodel, calculate xdata[i] in model
    dyda[3*i] = ex # B0, B1, B2
    dyda[3*i+1] = fac/par2[i] # E0, E1, E2
    dyda[3*i+2] = fac*arg/par2[i] # G0, G1, G2
  # print y, dyda
  results = [y, dyda]
  return results

# calculate chi2, beta 
# calculate alpha use (15.5.11)
def mrqcof(pars_trial):
  # initializing
  alpha = [ [0 for i in range(nfit)] for i in range(nfit) ]
  beta = [ 0 for i in range(nfit) ]
  chisq = 0
  # loop the data
  for i in range(npoints):
    results = fgauss(xdata[i],pars_trial)
    ymodel = results[0]
    dyda = results[1]
    dy = ydata[i] - ymodel
    sig2i = 1.#/(100*100) # set sigma = 1 for convenience, otherwise should come from data
    j = 0
    for l in range(nfit):
       #if par_isfine[l]: # reserved for free and hold
       wt = dyda[l]*sig2i
       k = 0
       for m in range(l+1):
         alpha[j][k] += wt*dyda[m] # use (15.5.11)
         k = k+1
       beta[j] += dy*wt
       j = j+1
    chisq += dy*dy*sig2i # for data loop
  for j in range(1, nfit):
    for k in range(j):
      alpha[k][j] = alpha[j][k]
  results = [chisq, alpha, beta]
  return results

# fitmrq mpars_start = [1., 1. , 1., 100., 100., 100., 20.,  20.,  20.]

alambda = 0.001
oneda = [ [1] for k in range(nfit) ] # nfit*1 vector
# trial parameters
alpha = [ [0 for i in range(nfit)] for i in range(nfit) ]
beta = [0 for kk in range(nfit)]
pars_start = [1., 1., 1., 10., 8., 8., 99.,  99.,  99.]

atry = pars_start
a = atry
cov = [ [0 for i in range(nfit)] for kk in range(nfit) ]
temp = [ [0 for i in range(nfit)] for kk in range(nfit) ]
da =  [0 for i in range(nfit)]
results_start = mrqcof(a)
chisqOld = results_start[0]
alpha = results_start[1]
beta = results_start[2]
nGoodfits = 0
#print "alpha", alpha
#print "beta", beta
for iter in range(max_iter):
    if nGoodfits == NGOOD:
       alambda = 0.0
       print "reach goodfits: ", nGoodfits
       break
    for j in range(nfit):
       for k in range(nfit):
          cov[j][k] = alpha[j][k]
       cov[j][j] = alpha[j][j]*(1.0+alambda) # calculate alpha'
       for k in range(nfit):
          temp[j][k] = cov[j][k]
       oneda[j][0] = beta[j] # delta a
       
    #print('\n'.join([' '.join(['{:3}'.format(item) for item in row]) for row in cov])) # check cov
    #print oneda
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
    #print "A:", Acov
    #print "B:", Boneda
    #print "x:", solveX

    for j in range(nfit):
        oneda[j] = solveX[j]
    #print oneda
    #for j in range(nfit):
    #   #for k in range(nfit):
    #   #   cov[j][k] = temp[j][k]
        da[j] = oneda[j][0]
    j = 0
    for l in range(nfit):
      atry[l] = a[l] + da[j] #reserve for hold/free
      j = j + 1 
    # print "atry=", atry, "a=", a, "da=", da

    results_new = mrqcof(atry)
    chisq = results_new[0]
    #print "chiNew, chiOld", chisq, chisqOld
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
        if alambda>1e8:
            print "lambda range out"
            break
    print alambda
print "best fit results: ", a

hfit = TH1F("hfit","Example of several fits in subranges", npoints, xdata[0], xdata[npoints-1]);
i = 0
for xx in xdata:
  results_final = fgauss(xx, a) 
  y = results_final[0]
  hfit.SetBinContent(i+1,y)
  i = i+1

hfit.SetLineColor(2)
hfit.Draw()
hdata.Draw("sames")
raw_input("enter to quit")
