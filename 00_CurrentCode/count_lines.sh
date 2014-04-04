#!/bin/bash

let i=`cat *.cpp *.h | wc -l`
let k=`cat MakeParamsList.cpp | wc -l`
let l=`ls -l *.cpp *.h | wc -l`-1
let j=$i-$k-$l*20-3
echo " "
echo "GeneStream project has $j lines of code, comments and doc strings."
echo " "
