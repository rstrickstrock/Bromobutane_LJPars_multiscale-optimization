#! /bin/bash
result_cwd=$1

cat $result_cwd"/this_prediction.csv" | cut -f2 -d ' '
