rm listScan.dat
ls MPWsolarTest*>>listScan.dat
sed 's/_/ /g' listScan.dat>>listScan.sort | mv listScan.sort listScan.dat
sort -k 4 listScan.dat>>temp.dat | mv temp.dat listScan.dat
sed 's/ /_/g' listScan.dat>>temp.dat | mv temp.dat listScan.dat
