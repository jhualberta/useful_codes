import json, glob, os, sys
# ratdb -s "postgres://snoplus:password@pgsql.snopl.us:5400/ratdb" dump_table LIVETIME 200004-207476
# need to change *.ratdb to *.json
path=os.getcwd()
file_list = os.listdir(path)

read_files = glob.glob("*.json")
fList = open('LIVETIME_200004-207476.json')
output_list = []
rawlivetime = []
livetime = []
# check whether a run is invalid (empty json files)
validfiles = []

for i in fList:
    d = json.loads(i)
    rawlivetime.append(d['raw_livetime'])    
    livetime.append(d['livetime'])

print "raw live time", sum(rawlivetime), "seconds", sum(rawlivetime)/60/60/24, "days"
print "live time", sum(livetime), "seconds", sum(livetime)/60/60/24, "days"
