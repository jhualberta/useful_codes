import json, glob, os
read_files = glob.glob("*.ratdb")
output_list = []
position = []
pos_err = []

# check whether a run is invalid (empty json files)
validfiles = []
for f in read_files:
  check = (os.stat(f).st_size == 0)
  if check == True:
      print "invalid runs ", f
  else:
      validfiles.append(f)

for i in range(0, len(validfiles) ):
  with open(validfiles[i], "r") as infile:
    d = json.load(infile)
    pos = d['position']
    position.append(pos)
    print "manip position", pos[0], pos[1], pos[2]
    infile.close()

