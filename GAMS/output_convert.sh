#!/bin/bash

# Example:
# You need to compare the upper bound of MIU before any second optimisation has been done
#   ./output_convert.sh d2013RRR_Stern-output.txt MIU upper first
# Then copy/pasta into jupyter.

file=$1
variable=$2

case "$3" in
   "upper")
     cmd=4
     ;;
   "lower")
     cmd=2
     ;;
   "level")
     cmd=3
     ;;
   "marginal")
     cmd=5
     ;;
esac

case "$4" in
    "first")
        data=$(rg -m 1 -A 63 $variable $1 | tail -n +5 | tr -s ' ' | cut -d ' ' -f $cmd | tr '\n' ', ' | sed 's/EPS/eps\(\)/g; s/\.\,/NaN,/g')
        ;;
    "last")
        data=$(tac $1 | rg -m 1 -B 63 $variable | tac | tail -n +5 | tr -s ' ' | cut -d ' ' -f $cmd | tr '\n' ', ' | sed 's/EPS/eps\(\)/g; s/\.\,/NaN,/g')
        ;;
esac

echo '['"${data:0:-1}"']'

