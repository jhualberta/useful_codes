#!/bin/bash
N=1
cat $1 | while read line
do
    echo "$line" > file$N
    N=`expr $N + 1`
done
