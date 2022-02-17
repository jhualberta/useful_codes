import numpy as np
from numpy import log, sqrt
# zhihua zhou, watermelon book, P75
def logp(p):
  r = p*log(p)/log(2)
  return r

a = [1,2,3,4,5,6,8,10,15]
nhard = 6
nsoft = 3
enthard = -( 1*(log(1)/log(2)))
entsoft = -(1./3*log(1./3)/log(2)+2./3*log(2./3)/log(2))
# 6 hard, all good
# 3 soft, 1 good
entTot = -(7./9*(log(7./9)/log(2))+2./9*(log(2./9)/log(2)))
#print entTot-(6./9*enthard + 3./9*entsoft)
# sound
entTot_sound = entTot
entLoud = -(logp(5./6)+logp(1./6))
entDumb = -(2./2*logp(1)) 
entCrispy = 0
print entTot_sound-(6./9*entLoud+2./9*entDumb+1./9*entCrispy)



