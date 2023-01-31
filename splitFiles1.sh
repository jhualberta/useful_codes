split -dl 50 -a 3 --additional-suffix=.sh allProcessRuns.sh processA1_

k=0
for f in processA1_0*; 
  do mv $f processA1_$k.sh;
	  k=$((k+1))
	  echo "Renamed processA1_"$k
done
