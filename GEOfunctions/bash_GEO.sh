#!/bin/bash

# Bash script to preprocess GEO matrix files to be read into R
# First input is the matrix file

FN=$1
FNoutSer=${FN/.txt/_series.txt}
FNoutSam=${FN/.txt/_samples.txt}
FNoutData=${FN/.txt/_data.txt}
echo $FNoutSer
cat $FN | grep !Series > $FNoutSer
echo $FNoutSam
cat $FN | grep !Sample > $FNoutSam
echo $FNoutData
cat $FN | egrep "ID_REF|G[0-9]{6}" > $FNoutData
