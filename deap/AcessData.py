import sys,os
import subprocess

#runList = [30681, 30686, 30695, 30705, 30706, 30707, 30708, 30717, 30726, 30741, 30742, 30743, 30744, 30746, 30747, 30751, 30756, 30760, 30765, 30769, 30774, 30784, 30785, 30804, 30813, 30815, 30826, 30837]

## runList = [27696, 27697, 27699, 27704, 27708, 27713, 27717, 27722, 27726, 27732, 27747, 27752, 27756, 27765, 27770, 27779, 27783, 27793, 27798, 27807, 27811, 27822, 27827, 27842, 27844, 27851, 27863, 27867, 27876, 27881, 27890, 27894, 27899, 27903, 27908, 27909, 27910, 27911, 27927, 27936, 27940, 27950, 27951, 27952, 27953, 27954, 27955, 27965, 27969, 27975, 27979, 27980, 27981, 27982, 27987, 27988, 27989, 27994, 27999, 28008, 28013, 28023, 28026, 28030, 28039, 28044, 28046, 28055, 28059, 28066, 28067, 28071, 28076, 28085, 28089, 28090, 28092, 28103, 28105, 28110, 28119, 28123, 28132, 28137, 28147, 28160, 28164, 28173, 28178, 28188, 28192, 28193, 28203, 28208, 28218, 28222, 28225, 28234, 28239, 28248, 28252, 28261, 28266, 28275, 28279, 28288, 28294, 28295, 28310, 28314, 28323, 28325, 28330, 28331, 28343, 28348, 28357, 28362, 28363, 28372, 28374, 28375, 28379, 28381, 28383, 28393, 28398, 28409, 28413, 28422, 28427, 28436, 28440, 28449, 28454, 28464, 28468, 28472, 28481, 28486, 28487, 28488, 28512, 28516, 28525, 28530, 28545, 28549, 28565, 28570, 28579, 28583, 28593, 28595, 28608, 28624, 28632, 28637, 28649, 28653, 28655, 28659, 28660, 28670, 28676, 28685, 28689, 28702, 28708, 28730, 28734, 28743, 28748, 28758, 28763, 28786, 28795, 28800, 28809, 28813, 28815, 28824, 28829]

runList = [28103, 28105, 28110, 28119, 28123, 28132, 28137, 28147, 28160, 28164, 28173, 28188, 28192, 28193, 28203, 28208, 28218, 28222, 28225, 28234, 28239, 28248, 28252, 28261, 28266, 28275, 28279, 28288, 28294, 28295, 28310, 28314, 28323, 28325, 28330, 28331, 28343, 28348, 28357, 28362, 28363, 28372, 28374, 28375, 28379, 28381, 28383, 28393, 28398, 28409, 28413, 28422, 28427, 28436, 28440, 28449, 28454, 28464, 28468, 28472, 28481, 28486, 28487, 28512, 28516, 28525, 28530, 28545, 28549, 28565, 28570, 28579, 28583, 28593, 28595, 28608, 28624, 28632, 28637, 28655, 28659, 28660, 28670, 28676, 28685, 28689, 28702, 28708, 28730, 28734, 28743, 28748, 28758, 28763, 28786, 28795, 28800, 28809, 28813, 28815, 28824, 28829]

fileList = []
for run in runList:
    #subprocess.call(["cd","/project/6004969/data/v5.14.0/cal/run0"+str(run)])
    #path=os.getcwd()
    #file_list = os.listdir(path)
    try:
        fileList.append(subprocess.check_output(["ls", "-1","/project/6004969/data/v5.14.0/ntp/run0"+str(run)]))
    except:
        print "this run", run, "not exist"
        continue  

fileList = [i.split('\n') for i in fileList]

for item in fileList:
    for nn in item:
         if nn == '': continue
         print "/project/6004969/data/v5.14.0/ntp/run"+nn[nn.index('0'):nn.index('0')+6]+"/"+nn
