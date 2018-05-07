#!/bin/bash

echo
echo "This just computes the neighbors:"

cmd="../../bin/neighbors -s 2 -i example1_matrix.txt -p 1e-4 -b example1_subject_bias.txt -f example1_features.txt -d ../../db"

echo
echo $cmd "> test_output1.txt"
$cmd > test_output1.txt 

echo
echo "Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output1.txt test_output1.txt"
echo
echo $cmd
$cmd

echo
echo "This command uses the neighbor relationships for DBSCAN clustering"
cmd="../../bin/neighbors -s 2 -i example1_matrix.txt -p 1e-4 -b example1_subject_bias.txt -f example1_features.txt -c -m 4 -z 4 -d ../../db"

echo
echo $cmd "> test_output2.txt"
$cmd > test_output2.txt

echo
echo "Now look for differences in the cluster output lines versus expected. Because the"
echo "output logfile contains the Z_CO scores and these are based on random shuffling,"
echo "the expected and observed logfiles will differ (but the Z_CO scores should be pretty"
echo "similar)."
echo
echo
echo "Comparing output to expected output (this command should produce no output):"
echo
echo "diff <(grep ^cluster expected_output2.txt) <(grep ^cluster test_output2.txt)"
diff <(grep ^cluster expected_output2.txt) <(grep ^cluster test_output2.txt)

