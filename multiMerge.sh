rm Merge*.dat

split -dl 900 -a 2 --additional-suffix=.dat allMerge.dat Merge_ 
k=0
for f in Merge_0*;
  do mv $f Merge_$k.dat;
          k=$((k+1))
          echo "Renamed Merge_"$k
done

echo "total merge list is:"$k

i=0
#for i in $(seq 0 $k)
for ff in Merge_*;
do
   echo "processing file Merge_"$i
   mv $ff MergeList.dat
   root -q newLoadData.C
   mv Merged_Resol.root merged_$i.root
   rm MergeList.dat
   i=$((i+1))
done
