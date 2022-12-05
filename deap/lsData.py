#/usr/bin/python
import sys,os
import subprocess as sp
path=os.getcwd()
file_list = os.listdir(path)
#/usr/bin/python
import sys,os
import subprocess
path=os.getcwd()
file_list = os.listdir(path)

runList = [30681, 30686, 30695, 30705, 30706, 30707, 30708, 30717, 30726, 30741, 30742, 30743, 30744, 30746, 30747, 30751, 30756, 30760, 30765, 30769, 30774, 30784, 30785, 30804, 30813, 30815, 30826, 30837]
path = "/project/6004969/data/v5.14.0/cal/"
subprocess.call(["cd",path])
basename = os.path.basename(path)
print basename
allFile = []
for run in runList:
      p = sp.Popen(["ls "+path+"/*"+str(run)+"*"], shell=True, stdout=sp.PIPE, universal_newlines=True)
      out = p.communicate()
      for fname in out:
          if fname != None and fname != "\n":
              #print fname
              allFile.append(path+fname)

for fname in allFile:
    print fname
