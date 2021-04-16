#  simpleFit.py
#  G. Cowan / RHUL Physics / October 2017
#  Simple program to illustrate least-squares fitting with curve_fit

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# set data values
x   = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
y   = np.array([1.3, 2.5, 2.6, 3.5, 3.1, 2.2, 1.7, 1.5])
sig = np.array([0.15, 0.2, 0.25, 0.3, 0.3, 0.2, 0.15, 0.15])


# define fit function (order of polynomial set using numPar below).
def func(x, *theta):
    m = len(theta)
    f = 0.0
    for i in range(m):
        f += theta[i]*pow(x,i)
    return f

# set default parameter values and do the fit
numPar = 3                     # set number of parameters here
p0 = np.array(numPar*[1.0])
thetaHat, cov = curve_fit(func, x, y, p0, sig, absolute_sigma=True)

# Retrieve minimized chi-squared, etc.
numPoints = len(x)
ndof = numPoints - numPar
chisq = sum(((y - func(x, *thetaHat))/sig)**2)
print "chisq = ", chisq, ",     ndof = ", ndof

# Print fit parameters and covariance matrix
print "\n", "Fitted parameters and standard deviations:"
sigThetaHat = np.sqrt(np.diag(cov))
for i in range(len(thetaHat)):
    print "thetaHat[", i, "] = ", thetaHat[i], "  +-  ", sigThetaHat[i]

print "\n", "i, j, cov[i,j], rho[i,j]:"
for i in range(len(thetaHat)):
    for j in range(len(thetaHat)):
        rho = cov[i][j] / (sigThetaHat[i]*sigThetaHat[j])
        print i, "  ", j, "  ", cov[i][j], "  ", rho

# Set up plot
matplotlib.rcParams.update({'font.size':18})     # set all font sizes
plt.clf()
fig, ax = plt.subplots(1,1)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.errorbar(x, y, yerr=sig, xerr=0, color='black', fmt='o', label='data')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$', labelpad=10)
xMin = 0
xMax = 10
yMin = 0
yMax = 5
plt.xlim(xMin, xMax)
plt.ylim(yMin, yMax)
xPlot = np.linspace(xMin, xMax, 100)        # enough points for a smooth curve
fit = func(xPlot, *thetaHat)
plt.plot(xPlot, fit, 'red', linewidth=2, label='fit result')

# Tweak legend
handles, labels = ax.get_legend_handles_labels()
handles = [handles[1], handles[0]]
labels = [labels[1], labels[0]]
handles = [handles[0][0], handles[1]]      # turn off error bar for data in legend
plt.legend(handles, labels, loc='upper right', fontsize=14, frameon=False)

# Make and store plot
plt.show()
plt.savefig("simpleFit.pdf", format='pdf')
