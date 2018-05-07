#!/bin/bash

echo
echo "Running test calculation on 25 TCRs..."
cmd="../../bin/pgen -i tcrs.txt -d ../../db/"
echo
echo $cmd "> test_output.txt"
$cmd > test_output.txt

echo
echo "Comparing output to expected output..."
cmd="diff expected_output.txt test_output.txt"
echo
echo $cmd
$cmd

echo
echo "Did that last command generate any big differences in pvalues?"
echo "(No output means diff didn't find any differences)"
