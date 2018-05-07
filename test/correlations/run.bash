#!/bin/bash

echo
echo "# Running test calculation"
cmd="../../bin/correlations -d ../../db -m matrix_small.txt -f cmv_features.txt -p 1e-3 -q 1e-2"
echo
echo $cmd "> test_output.txt"
$cmd > test_output.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output.txt test_output.txt"
echo
echo $cmd
$cmd




