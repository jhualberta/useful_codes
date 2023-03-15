#/usr/bin/python
import operator
import sys,os
import collections
import subprocess
from subprocess import *
from ROOT import *
from rat import *
import numpy as np
from numpy import array
#### roi 55 data ############
npoints = 179
x = [0 for k in range(npoints)]
y = [0 for k in range(npoints)]
x[0]=111; y[0]=0.684439;
x[1]=112; y[1]=0.681856;
x[2]=113; y[2]=0.679527;
x[3]=114; y[3]=0.677336;
x[4]=115; y[4]=0.674967;
x[5]=116; y[5]=0.672673;
x[6]=117; y[6]=0.670425;
x[7]=118; y[7]=0.668411;
x[8]=119; y[8]=0.666301;
x[9]=120; y[9]=0.663939;
x[10]=121; y[10]=0.661738;
x[11]=122; y[11]=0.659478;
x[12]=123; y[12]=0.657319;
x[13]=124; y[13]=0.655317;
x[14]=125; y[14]=0.653208;
x[15]=126; y[15]=0.650971;
x[16]=127; y[16]=0.649363;
x[17]=128; y[17]=0.647897;
x[18]=129; y[18]=0.645601;
x[19]=130; y[19]=0.643595;
x[20]=131; y[20]=0.642117;
x[21]=132; y[21]=0.64023;
x[22]=133; y[22]=0.638307;
x[23]=134; y[23]=0.636658;
x[24]=135; y[24]=0.634565;
x[25]=136; y[25]=0.632896;
x[26]=137; y[26]=0.631477;
x[27]=138; y[27]=0.629822;
x[28]=139; y[28]=0.628284;
x[29]=140; y[29]=0.626446;
x[30]=141; y[30]=0.62475;
x[31]=142; y[31]=0.623116;
x[32]=143; y[32]=0.621598;
x[33]=144; y[33]=0.620073;
x[34]=145; y[34]=0.61844;
x[35]=146; y[35]=0.616892;
x[36]=147; y[36]=0.615514;
x[37]=148; y[37]=0.613992;
x[38]=149; y[38]=0.612619;
x[39]=150; y[39]=0.611277;
x[40]=151; y[40]=0.609653;
x[41]=152; y[41]=0.607869;
x[42]=153; y[42]=0.606611;
x[43]=154; y[43]=0.605346;
x[44]=155; y[44]=0.603965;
x[45]=156; y[45]=0.60277;
x[46]=157; y[46]=0.600999;
x[47]=158; y[47]=0.599367;
x[48]=159; y[48]=0.598301;
x[49]=160; y[49]=0.597359;
x[50]=161; y[50]=0.59612;
x[51]=162; y[51]=0.594866;
x[52]=163; y[52]=0.593814;
x[53]=164; y[53]=0.592567;
x[54]=165; y[54]=0.593081;
x[55]=166; y[55]=0.593607;
x[56]=167; y[56]=0.594146;
x[57]=168; y[57]=0.594698;
x[58]=169; y[58]=0.595264;
x[59]=170; y[59]=0.595844;
x[60]=171; y[60]=0.596439;
x[61]=172; y[61]=0.597048;
x[62]=173; y[62]=0.597673;
x[63]=174; y[63]=0.598314;
x[64]=175; y[64]=0.59897;
x[65]=176; y[65]=0.599644;
x[66]=177; y[66]=0.600205;
x[67]=178; y[67]=0.600636;
x[68]=179; y[68]=0.601075;
x[69]=180; y[69]=0.601523;
x[70]=181; y[70]=0.601978;
x[71]=182; y[71]=0.602441;
x[72]=183; y[72]=0.602913;
x[73]=184; y[73]=0.603394;
x[74]=185; y[74]=0.603883;
x[75]=186; y[75]=0.604381;
x[76]=187; y[76]=0.604888;
x[77]=188; y[77]=0.605404;
x[78]=189; y[78]=0.60593;
x[79]=190; y[79]=0.606464;
x[80]=191; y[80]=0.607008;
x[81]=192; y[81]=0.607562;
x[82]=193; y[82]=0.608125;
x[83]=194; y[83]=0.608698;
x[84]=195; y[84]=0.60928;
x[85]=196; y[85]=0.609873;
x[86]=197; y[86]=0.610281;
x[87]=198; y[87]=0.610642;
x[88]=199; y[88]=0.611006;
x[89]=199; y[89]=0.704787;
x[90]=198; y[90]=0.704616;
x[91]=197; y[91]=0.704443;
x[92]=196; y[92]=0.704267;
x[93]=195; y[93]=0.704089;
x[94]=194; y[94]=0.703907;
x[95]=193; y[95]=0.703724;
x[96]=192; y[96]=0.703538;
x[97]=191; y[97]=0.703349;
x[98]=190; y[98]=0.703159;
x[99]=189; y[99]=0.702966;
x[100]=188; y[100]=0.702771;
x[101]=187; y[101]=0.702574;
x[102]=186; y[102]=0.702375;
x[103]=185; y[103]=0.702174;
x[104]=184; y[104]=0.701972;
x[105]=183; y[105]=0.701768;
x[106]=182; y[106]=0.701562;
x[107]=181; y[107]=0.701355;
x[108]=180; y[108]=0.701147;
x[109]=179; y[109]=0.700938;
x[110]=178; y[110]=0.700728;
x[111]=177; y[111]=0.700517;
x[112]=176; y[112]=0.700305;
x[113]=175; y[113]=0.700093;
x[114]=174; y[114]=0.699875;
x[115]=173; y[115]=0.699654;
x[116]=172; y[116]=0.699432;
x[117]=171; y[117]=0.699211;
x[118]=170; y[118]=0.698991;
x[119]=169; y[119]=0.698772;
x[120]=168; y[120]=0.698553;
x[121]=167; y[121]=0.698336;
x[122]=166; y[122]=0.69812;
x[123]=165; y[123]=0.697905;
x[124]=164; y[124]=0.697692;
x[125]=163; y[125]=0.69748;
x[126]=162; y[126]=0.69727;
x[127]=161; y[127]=0.697062;
x[128]=160; y[128]=0.696856;
x[129]=159; y[129]=0.696653;
x[130]=158; y[130]=0.696451;
x[131]=157; y[131]=0.696252;
x[132]=156; y[132]=0.696055;
x[133]=155; y[133]=0.695861;
x[134]=154; y[134]=0.695669;
x[135]=153; y[135]=0.69548;
x[136]=152; y[136]=0.695294;
x[137]=151; y[137]=0.695111;
x[138]=150; y[138]=0.69493;
x[139]=149; y[139]=0.694752;
x[140]=148; y[140]=0.694576;
x[141]=147; y[141]=0.694404;
x[142]=146; y[142]=0.694234;
x[143]=145; y[143]=0.694066;
x[144]=144; y[144]=0.693901;
x[145]=143; y[145]=0.693738;
x[146]=142; y[146]=0.693577;
x[147]=141; y[147]=0.693417;
x[148]=140; y[148]=0.69326;
x[149]=139; y[149]=0.693104;
x[150]=138; y[150]=0.692948;
x[151]=137; y[151]=0.692794;
x[152]=136; y[152]=0.69264;
x[153]=135; y[153]=0.692486;
x[154]=134; y[154]=0.692331;
x[155]=133; y[155]=0.692176;
x[156]=132; y[156]=0.692019;
x[157]=131; y[157]=0.69186;
x[158]=130; y[158]=0.691699;
x[159]=129; y[159]=0.691534;
x[160]=128; y[160]=0.691366;
x[161]=127; y[161]=0.691194;
x[162]=126; y[162]=0.691017;
x[163]=125; y[163]=0.690834;
x[164]=124; y[164]=0.690645;
x[165]=123; y[165]=0.690449;
x[166]=122; y[166]=0.690245;
x[167]=121; y[167]=0.690033;
x[168]=120; y[168]=0.689804;
x[169]=119; y[169]=0.689565;
x[170]=118; y[170]=0.689315;
x[171]=117; y[171]=0.689054;
x[172]=116; y[172]=0.688781;
x[173]=115; y[173]=0.688496;
x[174]=114; y[174]=0.688198;
x[175]=113; y[175]=0.687887;
x[176]=112; y[176]=0.687562;
x[177]=111; y[177]=0.687223;
x[178]=111; y[178]=0.684439;

roi55 = TCutG("top55",npoints)
roi55.SetVarX("nscb");
roi55.SetVarY("rprompt");

for i in range(npoints):
    roi55.SetPoint(i,x[i], y[i])

### drawing Sideband for roi55
xline = []
yline = []
for i in range(npoints): ### save the points for top curve
    if y[i]>0.687:
        xline.append(float(x[i]))
        yline.append(float(y[i]))

nlinepoints = len(xline)
print "save points", nlinepoints ### 89 points!!!
nscbMin = min(xline) 
nscbMax = max(xline) 
idx1 = xline.index(nscbMin) ## left y vertical
idx2 = xline.index(nscbMax) ## right y vertical

yscbMin = yline[idx1] ## the y-value (rprompt) of the nscbMin
yscbMax = yline[idx2] ## the y-value (rprompt) of the nscbMax

pointsteps = 100 #nscbMax - nscbMin
print "ROI55 ##################"
print nscbMin, nscbMax, yscbMin, yscbMax
print xline
print yline

print "set steps", nscbMax - nscbMin 

## print leftx
## print lefty

allPoints = nlinepoints + pointsteps + pointsteps + pointsteps
cutSideBand55 = TCutG("sideband55", allPoints) ## for top55
cutSideBand55.SetVarX("nscb");
cutSideBand55.SetVarY("rprompt");

zipped = zip(xline, yline)
res = sorted( zipped, key = lambda xx: xx[0]) ## element in zipped is (x,y), sort by x or nscb

### sorted lists
xline = [xx[0] for xx in res] ## for zipped pair elements
yline = [xx[1] for xx in res]
#

### counter clock-wise to fill the CutG, bottom tilted line
## from x = 111 to 199
for i in range(nlinepoints):
    ### print "points #", i
    #if yline[i]<0.1:
    #    print "what??"
    #print xline[i], yline[i]
    cutSideBand55.SetPoint(i,xline[i], yline[i])

## right vertical line, x = 199, y from bottom to top
rightx = [round(nscbMax,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
righty = [round(kk,6) for kk in np.linspace(yscbMax, 1.0, pointsteps)]

for i in range(pointsteps):
    cutSideBand55.SetPoint(i+nlinepoints, rightx[i], righty[i])

## top horizontal line for sideband, from right to left
topx = [round(kk,0) for kk in np.linspace(nscbMax, nscbMin, pointsteps)]
topy = [1.0 for i in np.linspace(nscbMin, nscbMax, pointsteps)]
print "top", topx, topy
#print topy, topx

for i in range(pointsteps):
    cutSideBand55.SetPoint(i+nlinepoints+pointsteps, topx[i], topy[i])
    
## left vertical line, from top to bottom
leftx = [round(nscbMin,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
lefty = [round(kk,6) for kk in np.linspace(1.0, yscbMin, pointsteps)]

for i in range(pointsteps):
    cutSideBand55.SetPoint(i+nlinepoints+pointsteps + pointsteps, leftx[i], lefty[i])
### cutSideBand55.RemovePoint(0,0)


### roi_top30 data
npoints30 = 202+1 
x30 = [0 for k in range(npoints30)]
y30 = [0 for k in range(npoints30)]
x30[0]=99; y30[0]=0.717298
x30[1]=100; y30[1]=0.714353
x30[2]=101; y30[2]=0.711464
x30[3]=102; y30[3]=0.708478
x30[4]=103; y30[4]=0.705839
x30[5]=104; y30[5]=0.702971
x30[6]=105; y30[6]=0.700019
x30[7]=106; y30[7]=0.697562
x30[8]=107; y30[8]=0.694773
x30[9]=108; y30[9]=0.692056
x30[10]=109; y30[10]=0.689742
x30[11]=110; y30[11]=0.68699
x30[12]=111; y30[12]=0.684439
x30[13]=112; y30[13]=0.681856
x30[14]=113; y30[14]=0.679527
x30[15]=114; y30[15]=0.677336
x30[16]=115; y30[16]=0.674967
x30[17]=116; y30[17]=0.672673
x30[18]=117; y30[18]=0.670425
x30[19]=118; y30[19]=0.668411
x30[20]=119; y30[20]=0.666301
x30[21]=120; y30[21]=0.663939
x30[22]=121; y30[22]=0.661738
x30[23]=122; y30[23]=0.659478
x30[24]=123; y30[24]=0.657319
x30[25]=124; y30[25]=0.655317
x30[26]=125; y30[26]=0.653208
x30[27]=126; y30[27]=0.650971
x30[28]=127; y30[28]=0.649363
x30[29]=128; y30[29]=0.647897
x30[30]=129; y30[30]=0.645601
x30[31]=130; y30[31]=0.643595
x30[32]=131; y30[32]=0.642117
x30[33]=132; y30[33]=0.64023
x30[34]=133; y30[34]=0.638307
x30[35]=134; y30[35]=0.636658
x30[36]=135; y30[36]=0.634565
x30[37]=136; y30[37]=0.632896
x30[38]=137; y30[38]=0.631477
x30[39]=138; y30[39]=0.629822
x30[40]=139; y30[40]=0.628284
x30[41]=140; y30[41]=0.626446
x30[42]=141; y30[42]=0.62475
x30[43]=142; y30[43]=0.623116
x30[44]=143; y30[44]=0.621598
x30[45]=144; y30[45]=0.620073
x30[46]=145; y30[46]=0.61844
x30[47]=146; y30[47]=0.616892
x30[48]=147; y30[48]=0.615514
x30[49]=148; y30[49]=0.613992
x30[50]=149; y30[50]=0.612619
x30[51]=150; y30[51]=0.611277
x30[52]=151; y30[52]=0.609653
x30[53]=152; y30[53]=0.607869
x30[54]=153; y30[54]=0.606611
x30[55]=154; y30[55]=0.605346
x30[56]=155; y30[56]=0.603965
x30[57]=156; y30[57]=0.60277
x30[58]=157; y30[58]=0.600999
x30[59]=158; y30[59]=0.599367
x30[60]=159; y30[60]=0.598301
x30[61]=160; y30[61]=0.597359
x30[62]=161; y30[62]=0.59612
x30[63]=162; y30[63]=0.594866
x30[64]=163; y30[64]=0.593814
x30[65]=164; y30[65]=0.592567
x30[66]=165; y30[66]=0.593081
x30[67]=166; y30[67]=0.593607
x30[68]=167; y30[68]=0.594146
x30[69]=168; y30[69]=0.594698
x30[70]=169; y30[70]=0.595264
x30[71]=170; y30[71]=0.595844
x30[72]=171; y30[72]=0.596439
x30[73]=172; y30[73]=0.597048
x30[74]=173; y30[74]=0.597673
x30[75]=174; y30[75]=0.598314
x30[76]=175; y30[76]=0.59897
x30[77]=176; y30[77]=0.599644
x30[78]=177; y30[78]=0.600205
x30[79]=178; y30[79]=0.600636
x30[80]=179; y30[80]=0.601075
x30[81]=180; y30[81]=0.601523
x30[82]=181; y30[82]=0.601978
x30[83]=182; y30[83]=0.602441
x30[84]=183; y30[84]=0.602913
x30[85]=184; y30[85]=0.603394
x30[86]=185; y30[86]=0.603883
x30[87]=186; y30[87]=0.604381
x30[88]=187; y30[88]=0.604888
x30[89]=188; y30[89]=0.605404
x30[90]=189; y30[90]=0.60593
x30[91]=190; y30[91]=0.606464
x30[92]=191; y30[92]=0.607008
x30[93]=192; y30[93]=0.607562
x30[94]=193; y30[94]=0.608125
x30[95]=194; y30[95]=0.608698
x30[96]=195; y30[96]=0.60928
x30[97]=196; y30[97]=0.609873
x30[98]=197; y30[98]=0.610281
x30[99]=198; y30[99]=0.610642
x30[100]=199; y30[100]=0.611006
x30[101]=199; y30[101]=0.730395
x30[102]=198; y30[102]=0.730285
x30[103]=197; y30[103]=0.730173
x30[104]=196; y30[104]=0.730058
x30[105]=195; y30[105]=0.729948
x30[106]=194; y30[106]=0.729842
x30[107]=193; y30[107]=0.729734
x30[108]=192; y30[108]=0.729625
x30[109]=191; y30[109]=0.729513
x30[110]=190; y30[110]=0.7294
x30[111]=189; y30[111]=0.729284
x30[112]=188; y30[112]=0.729168
x30[113]=187; y30[113]=0.729049
x30[114]=186; y30[114]=0.728929
x30[115]=185; y30[115]=0.728808
x30[116]=184; y30[116]=0.728686
x30[117]=183; y30[117]=0.728562
x30[118]=182; y30[118]=0.728437
x30[119]=181; y30[119]=0.728312
x30[120]=180; y30[120]=0.728185
x30[121]=179; y30[121]=0.728058
x30[122]=178; y30[122]=0.72793
x30[123]=177; y30[123]=0.727802
x30[124]=176; y30[124]=0.727674
x30[125]=175; y30[125]=0.727545
x30[126]=174; y30[126]=0.727417
x30[127]=173; y30[127]=0.727288
x30[128]=172; y30[128]=0.72716
x30[129]=171; y30[129]=0.727032
x30[130]=170; y30[130]=0.726905
x30[131]=169; y30[131]=0.726778
x30[132]=168; y30[132]=0.726653
x30[133]=167; y30[133]=0.726528
x30[134]=166; y30[134]=0.726404
x30[135]=165; y30[135]=0.726282
x30[136]=164; y30[136]=0.726162
x30[137]=163; y30[137]=0.726043
x30[138]=162; y30[138]=0.725926
x30[139]=161; y30[139]=0.725811
x30[140]=160; y30[140]=0.725698
x30[141]=159; y30[141]=0.725587
x30[142]=158; y30[142]=0.725479
x30[143]=157; y30[143]=0.725374
x30[144]=156; y30[144]=0.725271
x30[145]=155; y30[145]=0.725171
x30[146]=154; y30[146]=0.725073
x30[147]=153; y30[147]=0.724979
x30[148]=152; y30[148]=0.724888
x30[149]=151; y30[149]=0.7248
x30[150]=150; y30[150]=0.724716
x30[151]=149; y30[151]=0.724634
x30[152]=148; y30[152]=0.724556
x30[153]=147; y30[153]=0.724481
x30[154]=146; y30[154]=0.72441
x30[155]=145; y30[155]=0.724341
x30[156]=144; y30[156]=0.724276
x30[157]=143; y30[157]=0.724214
x30[158]=142; y30[158]=0.724154
x30[159]=141; y30[159]=0.724098
x30[160]=140; y30[160]=0.724044
x30[161]=139; y30[161]=0.723992
x30[162]=138; y30[162]=0.723942
x30[163]=137; y30[163]=0.723895
x30[164]=136; y30[164]=0.723849
x30[165]=135; y30[165]=0.723804
x30[166]=134; y30[166]=0.72376
x30[167]=133; y30[167]=0.723716
x30[168]=132; y30[168]=0.723673
x30[169]=131; y30[169]=0.723629
x30[170]=130; y30[170]=0.723584
x30[171]=129; y30[171]=0.723538
x30[172]=128; y30[172]=0.72349
x30[173]=127; y30[173]=0.723439
x30[174]=126; y30[174]=0.723386
x30[175]=125; y30[175]=0.723329
x30[176]=124; y30[176]=0.723268
x30[177]=123; y30[177]=0.723202
x30[178]=122; y30[178]=0.72313
x30[179]=121; y30[179]=0.723052
x30[180]=120; y30[180]=0.722968
x30[181]=119; y30[181]=0.722876
x30[182]=118; y30[182]=0.722776
x30[183]=117; y30[183]=0.722668
x30[184]=116; y30[184]=0.722551
x30[185]=115; y30[185]=0.722424
x30[186]=114; y30[186]=0.722286
x30[187]=113; y30[187]=0.722138
x30[188]=112; y30[188]=0.721979
x30[189]=111; y30[189]=0.721808
x30[190]=110; y30[190]=0.721625
x30[191]=109; y30[191]=0.72143
x30[192]=108; y30[192]=0.721223
x30[193]=107; y30[193]=0.721004
x30[194]=106; y30[194]=0.720772
x30[195]=105; y30[195]=0.720528
x30[196]=104; y30[196]=0.720272
x30[197]=103; y30[197]=0.720004
x30[198]=102; y30[198]=0.719747
x30[199]=101; y30[199]=0.719481
x30[200]=100; y30[200]=0.719205
x30[201]=99; y30[201]=0.718921
x30[202]=99; y30[202]=0.717298

roi30 = TCutG("top30",npoints30)
roi30.SetVarX("nscb");
roi30.SetVarY("rprompt");
for i in range(npoints30):
    roi30.SetPoint(i,x30[i], y30[i])

### drawing Sideband for roi30
xline = []
yline = []
for i in range(npoints30): 
    if y30[i]>0.718:
        xline.append(float(x30[i]))
        yline.append(float(y30[i]))

nlinepoints = len(xline)
print "save points for SB top30", nlinepoints ###
nscbMin = min(xline) 
nscbMax = max(xline) 
idx1 = xline.index(nscbMin) ## left y vertical
idx2 = xline.index(nscbMax) ## right y vertical

y_nscbMin = yline[idx1] ## the y-value (rprompt) of the nscbMin
y_nscbMax = yline[idx2] ## the y-value (rprompt) of the nscbMax

pointsteps = 100 #nscbMax - nscbMin

print "ROI 30 ############### "
print "xmin, xmax, ymin, ymax, xline, yline"
print nscbMin, nscbMax, y_nscbMin, y_nscbMax
print xline
print yline

print "set steps", nscbMax - nscbMin 

## print leftx
## print lefty

allPoints = nlinepoints + pointsteps + pointsteps + pointsteps
cutSideBand30 = TCutG("sideband30", allPoints) ## for top30
cutSideBand30.SetVarX("nscb");
cutSideBand30.SetVarY("rprompt");

zipped = zip(xline, yline)
res = sorted(zipped, key = lambda x: x[0]) ### using nscb to sort the pair

print "!!!!!!!!", xline
### sorted lists
xline = [xx[0] for xx in res] ### zipped element!!
yline = [xx[1] for xx in res]

### !!counter clock-wise to fill the CutG, bottom tilted line
## from x = 111 to 199
for i in range(nlinepoints):
    print "top30 points#", i
    #if yline[i]<0.1:
    #    print "what??"
    print xline[i], yline[i]
    cutSideBand30.SetPoint(i,xline[i], yline[i])

## right vertical line, x is fixed, x= nscbMax, y from bottom to top
rightx = [round(nscbMax,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
righty = [round(kk,6) for kk in np.linspace(y_nscbMax, 1.0, pointsteps)]

for i in range(pointsteps):
    cutSideBand30.SetPoint(i+nlinepoints, rightx[i], righty[i])

## top horizontal line for sideband, y is fixed, y =1, x from right to left
topx = [round(kk,0) for kk in np.linspace(nscbMax, nscbMin, pointsteps)]
topy = [1.0 for i in np.linspace(nscbMin, nscbMax, pointsteps)]
print "top", "x=", topx, "y=", topy
#print topy, topx

for i in range(pointsteps):
    cutSideBand30.SetPoint(i+nlinepoints+pointsteps, topx[i], topy[i])
    
## left vertical line, x is fixed to x= nscbMin, y from top y=1 to bottom y = y_nscbMin
leftx = [round(nscbMin,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
lefty = [round(kk,6) for kk in np.linspace(1.0, y_nscbMin, pointsteps)]

for i in range(pointsteps):
    cutSideBand30.SetPoint(i+nlinepoints+pointsteps+pointsteps, leftx[i], lefty[i])
### cutSideBand30.RemovePoint(0,0)

##### roi05 data
npoints05 = 238+1
x05 = [0 for k in range(npoints05)]
y05 = [0 for k in range(npoints05)]
x05[0]=81; y05[0]=0.779134
x05[1]=82; y05[1]=0.774714
x05[2]=83; y05[2]=0.770641
x05[3]=84; y05[3]=0.767196
x05[4]=85; y05[4]=0.763081
x05[5]=86; y05[5]=0.75925
x05[6]=87; y05[6]=0.755826
x05[7]=88; y05[7]=0.752281
x05[8]=89; y05[8]=0.74921
x05[9]=90; y05[9]=0.745817
x05[10]=91; y05[10]=0.742119
x05[11]=92; y05[11]=0.738677
x05[12]=93; y05[12]=0.735369
x05[13]=94; y05[13]=0.732153
x05[14]=95; y05[14]=0.729144
x05[15]=96; y05[15]=0.726132
x05[16]=97; y05[16]=0.723314
x05[17]=98; y05[17]=0.720575
x05[18]=99; y05[18]=0.717298
x05[19]=100; y05[19]=0.714353
x05[20]=101; y05[20]=0.711464
x05[21]=102; y05[21]=0.708478
x05[22]=103; y05[22]=0.705839
x05[23]=104; y05[23]=0.702971
x05[24]=105; y05[24]=0.700019
x05[25]=106; y05[25]=0.697562
x05[26]=107; y05[26]=0.694773
x05[27]=108; y05[27]=0.692056
x05[28]=109; y05[28]=0.689742
x05[29]=110; y05[29]=0.68699
x05[30]=111; y05[30]=0.684439
x05[31]=112; y05[31]=0.681856
x05[32]=113; y05[32]=0.679527
x05[33]=114; y05[33]=0.677336
x05[34]=115; y05[34]=0.674967
x05[35]=116; y05[35]=0.672673
x05[36]=117; y05[36]=0.670425
x05[37]=118; y05[37]=0.668411
x05[38]=119; y05[38]=0.666301
x05[39]=120; y05[39]=0.663939
x05[40]=121; y05[40]=0.661738
x05[41]=122; y05[41]=0.659478
x05[42]=123; y05[42]=0.657319
x05[43]=124; y05[43]=0.655317
x05[44]=125; y05[44]=0.653208
x05[45]=126; y05[45]=0.650971
x05[46]=127; y05[46]=0.649363
x05[47]=128; y05[47]=0.647897
x05[48]=129; y05[48]=0.645601
x05[49]=130; y05[49]=0.643595
x05[50]=131; y05[50]=0.642117
x05[51]=132; y05[51]=0.64023
x05[52]=133; y05[52]=0.638307
x05[53]=134; y05[53]=0.636658
x05[54]=135; y05[54]=0.634565
x05[55]=136; y05[55]=0.632896
x05[56]=137; y05[56]=0.631477
x05[57]=138; y05[57]=0.629822
x05[58]=139; y05[58]=0.628284
x05[59]=140; y05[59]=0.626446
x05[60]=141; y05[60]=0.62475
x05[61]=142; y05[61]=0.623116
x05[62]=143; y05[62]=0.621598
x05[63]=144; y05[63]=0.620073
x05[64]=145; y05[64]=0.61844
x05[65]=146; y05[65]=0.616892
x05[66]=147; y05[66]=0.615514
x05[67]=148; y05[67]=0.613992
x05[68]=149; y05[68]=0.612619
x05[69]=150; y05[69]=0.611277
x05[70]=151; y05[70]=0.609653
x05[71]=152; y05[71]=0.607869
x05[72]=153; y05[72]=0.606611
x05[73]=154; y05[73]=0.605346
x05[74]=155; y05[74]=0.603965
x05[75]=156; y05[75]=0.60277
x05[76]=157; y05[76]=0.600999
x05[77]=158; y05[77]=0.599367
x05[78]=159; y05[78]=0.598301
x05[79]=160; y05[79]=0.597359
x05[80]=161; y05[80]=0.59612
x05[81]=162; y05[81]=0.594866
x05[82]=163; y05[82]=0.593814
x05[83]=164; y05[83]=0.592567
x05[84]=165; y05[84]=0.593081
x05[85]=166; y05[85]=0.593607
x05[86]=167; y05[86]=0.594146
x05[87]=168; y05[87]=0.594698
x05[88]=169; y05[88]=0.595264
x05[89]=170; y05[89]=0.595844
x05[90]=171; y05[90]=0.596439
x05[91]=172; y05[91]=0.597048
x05[92]=173; y05[92]=0.597673
x05[93]=174; y05[93]=0.598314
x05[94]=175; y05[94]=0.59897
x05[95]=176; y05[95]=0.599644
x05[96]=177; y05[96]=0.600205
x05[97]=178; y05[97]=0.600636
x05[98]=179; y05[98]=0.601075
x05[99]=180; y05[99]=0.601523
x05[100]=181; y05[100]=0.601978
x05[101]=182; y05[101]=0.602441
x05[102]=183; y05[102]=0.602913
x05[103]=184; y05[103]=0.603394
x05[104]=185; y05[104]=0.603883
x05[105]=186; y05[105]=0.604381
x05[106]=187; y05[106]=0.604888
x05[107]=188; y05[107]=0.605404
x05[108]=189; y05[108]=0.60593
x05[109]=190; y05[109]=0.606464
x05[110]=191; y05[110]=0.607008
x05[111]=192; y05[111]=0.607562
x05[112]=193; y05[112]=0.608125
x05[113]=194; y05[113]=0.608698
x05[114]=195; y05[114]=0.60928
x05[115]=196; y05[115]=0.609873
x05[116]=197; y05[116]=0.610281
x05[117]=198; y05[117]=0.610642
x05[118]=199; y05[118]=0.611006
x05[119]=199; y05[119]=0.772743
x05[120]=198; y05[120]=0.772746
x05[121]=197; y05[121]=0.772748
x05[122]=196; y05[122]=0.772749
x05[123]=195; y05[123]=0.772749
x05[124]=194; y05[124]=0.772747
x05[125]=193; y05[125]=0.772744
x05[126]=192; y05[126]=0.77274
x05[127]=191; y05[127]=0.772735
x05[128]=190; y05[128]=0.77273
x05[129]=189; y05[129]=0.772723
x05[130]=188; y05[130]=0.772716
x05[131]=187; y05[131]=0.772708
x05[132]=186; y05[132]=0.7727
x05[133]=185; y05[133]=0.772691
x05[134]=184; y05[134]=0.772682
x05[135]=183; y05[135]=0.772672
x05[136]=182; y05[136]=0.772663
x05[137]=181; y05[137]=0.772654
x05[138]=180; y05[138]=0.772645
x05[139]=179; y05[139]=0.772637
x05[140]=178; y05[140]=0.772629
x05[141]=177; y05[141]=0.772622
x05[142]=176; y05[142]=0.772616
x05[143]=175; y05[143]=0.77261
x05[144]=174; y05[144]=0.772606
x05[145]=173; y05[145]=0.772604
x05[146]=172; y05[146]=0.772602
x05[147]=171; y05[147]=0.772603
x05[148]=170; y05[148]=0.772605
x05[149]=169; y05[149]=0.772609
x05[150]=168; y05[150]=0.772616
x05[151]=167; y05[151]=0.772625
x05[152]=166; y05[152]=0.772636
x05[153]=165; y05[153]=0.77265
x05[154]=164; y05[154]=0.772667
x05[155]=163; y05[155]=0.772687
x05[156]=162; y05[156]=0.77271
x05[157]=161; y05[157]=0.772737
x05[158]=160; y05[158]=0.772767
x05[159]=159; y05[159]=0.772801
x05[160]=158; y05[160]=0.772839
x05[161]=157; y05[161]=0.772881
x05[162]=156; y05[162]=0.772927
x05[163]=155; y05[163]=0.772977
x05[164]=154; y05[164]=0.773031
x05[165]=153; y05[165]=0.77309
x05[166]=152; y05[166]=0.773153
x05[167]=151; y05[167]=0.773221
x05[168]=150; y05[168]=0.773294
x05[169]=149; y05[169]=0.773371
x05[170]=148; y05[170]=0.773453
x05[171]=147; y05[171]=0.773539
x05[172]=146; y05[172]=0.773631
x05[173]=145; y05[173]=0.773726
x05[174]=144; y05[174]=0.773826
x05[175]=143; y05[175]=0.773931
x05[176]=142; y05[176]=0.77404
x05[177]=141; y05[177]=0.774152
x05[178]=140; y05[178]=0.774269
x05[179]=139; y05[179]=0.77439
x05[180]=138; y05[180]=0.774514
x05[181]=137; y05[181]=0.774641
x05[182]=136; y05[182]=0.774772
x05[183]=135; y05[183]=0.774905
x05[184]=134; y05[184]=0.775041
x05[185]=133; y05[185]=0.775179
x05[186]=132; y05[186]=0.775319
x05[187]=131; y05[187]=0.775461
x05[188]=130; y05[188]=0.775604
x05[189]=129; y05[189]=0.775747
x05[190]=128; y05[190]=0.775892
x05[191]=127; y05[191]=0.776036
x05[192]=126; y05[192]=0.77618
x05[193]=125; y05[193]=0.776323
x05[194]=124; y05[194]=0.776465
x05[195]=123; y05[195]=0.776606
x05[196]=122; y05[196]=0.776745
x05[197]=121; y05[197]=0.776882
x05[198]=120; y05[198]=0.777017
x05[199]=119; y05[199]=0.777148
x05[200]=118; y05[200]=0.777276
x05[201]=117; y05[201]=0.777401
x05[202]=116; y05[202]=0.777522
x05[203]=115; y05[203]=0.777639
x05[204]=114; y05[204]=0.777752
x05[205]=113; y05[205]=0.77786
x05[206]=112; y05[206]=0.777964
x05[207]=111; y05[207]=0.778062
x05[208]=110; y05[208]=0.778156
x05[209]=109; y05[209]=0.778246
x05[210]=108; y05[210]=0.77833
x05[211]=107; y05[211]=0.77841
x05[212]=106; y05[212]=0.778486
x05[213]=105; y05[213]=0.778557
x05[214]=104; y05[214]=0.778624
x05[215]=103; y05[215]=0.778688
x05[216]=102; y05[216]=0.778749
x05[217]=101; y05[217]=0.778807
x05[218]=100; y05[218]=0.778864
x05[219]=99; y05[219]=0.778919
x05[220]=98; y05[220]=0.778974
x05[221]=97; y05[221]=0.779029
x05[222]=96; y05[222]=0.779085
x05[223]=95; y05[223]=0.779144
x05[224]=94; y05[224]=0.779205
x05[225]=93; y05[225]=0.779271
x05[226]=92; y05[226]=0.779341
x05[227]=91; y05[227]=0.779417
x05[228]=90; y05[228]=0.779499
x05[229]=89; y05[229]=0.779589
x05[230]=88; y05[230]=0.779687
x05[231]=87; y05[231]=0.779793
x05[232]=86; y05[232]=0.779909
x05[233]=85; y05[233]=0.780044
x05[234]=84; y05[234]=0.780221
x05[235]=83; y05[235]=0.780409
x05[236]=82; y05[236]=0.780609
x05[237]=81; y05[237]=0.780819
x05[238]=81; y05[238]=0.779134
roi05 = TCutG("top05",npoints05)
roi05.SetVarX("nscb");
roi05.SetVarY("rprompt");
for i in range(npoints05):
    roi05.SetPoint(i,x05[i], y05[i])
    
### drawing Sideband for roi30
xline = []
yline = []
###NOTE:!! This curve is complicated, need to be careful!
for i in range(npoints05): 
    if y05[i]>0.772:
        if x05[i]<88:
            if x05[i] < 84:
                if y05[i]>0.78:
                    xline.append(float(x05[i]))
                    yline.append(float(y05[i]))
            else:
                if y05[i]>0.7726:
                   xline.append(float(x05[i]))
                   yline.append(float(y05[i]))
        else:
            xline.append(float(x05[i]))
            yline.append(float(y05[i]))

nlinepoints = len(xline)
print "save points for SB top30", nlinepoints ###

xxx1 = array(xline)
yyy1 = array(yline)
g = TGraph(len(xline),xxx1,yyy1)
#g.Draw("AP*")
#raw_input()

print "ROI 05 ###############"
print "xmin, xmax, ymin, ymax, xline, yline"
nscbMin = min(xline) 
nscbMax = max(xline) 
idx1 = xline.index(nscbMin) ## left y vertical
idx2 = xline.index(nscbMax) ## right y vertical

yscbMin = yline[idx1] ## the y-value (rprompt) of the nscbMin
yscbMax = yline[idx2] ## the y-value (rprompt) of the nscbMax

pointsteps = 100 #nscbMax - nscbMin

print nscbMin, nscbMax, yscbMin, yscbMax
print "set steps", nscbMax - nscbMin 

## print leftx
## print lefty

allPoints = nlinepoints + pointsteps + pointsteps + pointsteps
cutSideBand05 = TCutG("sideband05", allPoints) ## for top05
cutSideBand05.SetVarX("nscb");
cutSideBand05.SetVarY("rprompt");

zipped = zip(xline, yline)
res = sorted(zipped, key = lambda xx: xx[0]) ### using x to sort

### sorted lists
xline = [xx[0] for xx in res]
yline = [xx[1] for xx in res]
#
xxx2 = array(xline)
yyy2 = array(yline)
g1 = TGraph(len(xline),xxx2,yyy2)
## g1.Draw("AP*")
#raw_input()

print "x=", len(xline), xline
print "y=", yline

### !!counter-clock-wise to fill the CutG, bottom tilted line
## from x = 81 to 199
for i in range(nlinepoints):
    print "points #", i
    #if yline[i]<0.1:
    #    print "what??"
    print xline[i], yline[i]
    cutSideBand05.SetPoint(i,xline[i], yline[i])

## right vertical line, x = 199, y from bottom to top
rightx = [round(nscbMax,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
righty = [round(kk,6) for kk in np.linspace(yscbMax, 1.0, pointsteps)]

for i in range(pointsteps):
    cutSideBand05.SetPoint(i+nlinepoints, rightx[i], righty[i])

## top horizontal line for sideband, from right to left
topx = [round(kk,0) for kk in np.linspace(nscbMax, nscbMin, pointsteps)]
topy = [1.0 for i in np.linspace(nscbMin, nscbMax, pointsteps)]
print "top", topx, topy
#print topy, topx

for i in range(pointsteps):
    cutSideBand05.SetPoint(i+nlinepoints+pointsteps, topx[i], topy[i])
    
## left vertical line, from top to bottom
leftx = [round(nscbMin,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
lefty = [round(kk,6) for kk in np.linspace(1.0, yscbMin, pointsteps)]

for i in range(pointsteps):
    cutSideBand05.SetPoint(i+nlinepoints+pointsteps + pointsteps, leftx[i], lefty[i])
### cutSideBand05.RemovePoint(0,0)    

fnew = TFile("saveSideBandROI.root","recreate")
fnew.cd()
cutSideBand55.Write()
roi55.Write()
cutSideBand30.Write()
roi30.Write()
cutSideBand05.Write()
roi05.Write()

## cutSideBand.Draw()
## raw_input()
