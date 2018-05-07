#!/bin/bash

echo
echo "Compute tcrdist distance matrix:"

cmd="../../bin/tcrdists -i tcrs.txt -d ../../db"

echo
echo $cmd "> test_output.txt"
$cmd > test_output.txt

echo
echo "Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output.txt test_output.txt"
echo
echo $cmd
$cmd

