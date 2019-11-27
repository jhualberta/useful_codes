import json, glob, os

read_files = glob.glob("*.ratdb")
output_list = []
rawlivetime = []
livetime = []
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
    rawlivetime.append(d['raw_livetime'])    
    livetime.append(d['livetime'])
    infile.close()

print "raw live time", sum(rawlivetime), "seconds", sum(rawlivetime)/60/60/24, "days"
print "live time", sum(livetime), "seconds", sum(livetime)/60/60/24, "days"
