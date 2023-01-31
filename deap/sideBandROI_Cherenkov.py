#/usr/bin/python
import operator
import sys,os
import collections
import subprocess
from subprocess import *
from ROOT import *
from rat import *
import numpy as np
xline = []
yline = []
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

roi = TCutG("roi",npoints)
roi.SetVarX("nSCBayes");
roi.SetVarY("rprompt60Bayes");

for i in range(npoints):
    roi.SetPoint(i,x[i], y[i])

for i in range(npoints): 
    if y[i]>0.685:
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

print nscbMin, nscbMax, yscbMin, yscbMax
print xline
print yline

print "set steps", nscbMax - nscbMin 

## print leftx
## print lefty

allPoints = nlinepoints + pointsteps + pointsteps + pointsteps
cutSideBand = TCutG("cutSideBand", allPoints)
cutSideBand.SetVarX("nSCBayes");
cutSideBand.SetVarY("rprompt60Bayes");

zipped = zip(xline, yline)
res = sorted( zipped, key = lambda x: x[1])

### sorted lists
xline = [x[0] for x in res]
yline = [x[1] for x in res]
#

### counter clock-wise to fill the CutG, bottom tilted line
## from x = 111 to 199
for i in range(nlinepoints):
    print "points #", i
    #if yline[i]<0.1:
    #    print "what??"
    print xline[i], yline[i]
    cutSideBand.SetPoint(i,xline[i], yline[i])

## right vertical line, x = 199, y from bottom to top
rightx = [round(nscbMax,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
righty = [round(kk,6) for kk in np.linspace(yscbMax, 1.0, pointsteps)]

for i in range(pointsteps):
    cutSideBand.SetPoint(i+nlinepoints, rightx[i], righty[i])

## top horizontal line for side-band, from right to left
topx = [round(kk,0) for kk in np.linspace(nscbMax, nscbMin, pointsteps)]
topy = [1.0 for i in np.linspace(nscbMin, nscbMax, pointsteps)]
print "top", topx, topy
#print topy, topx

for i in range(pointsteps):
    cutSideBand.SetPoint(i+nlinepoints+pointsteps, topx[i], topy[i])


## left vertical line, from top to bottom
leftx = [round(nscbMin,0) for i in np.linspace(nscbMin, nscbMax, pointsteps)]
lefty = [round(kk,6) for kk in np.linspace(1.0, yscbMin, pointsteps)]

for i in range(pointsteps):
    cutSideBand.SetPoint(i+nlinepoints+pointsteps + pointsteps, leftx[i], lefty[i])

### cutSideBand.RemovePoint(0,0)

fnew = TFile("saveSideBandROI.root","recreate")
fnew.cd()
cutSideBand.Write()
roi.Write()

cutSideBand.Draw()
raw_input()
