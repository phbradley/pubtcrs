#!/bin/bash

echo
echo "# Compute tcrdist distance matrix:"

cmd="../../bin/tcrdists --tcrs_file1 Vbeta_family_tcrs.txt --database ../../db"

echo
echo $cmd "> test_output1.txt"
$cmd > test_output1.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output1.txt test_output1.txt"
echo
echo $cmd
$cmd

echo
echo "# Compute alpha chain tcrdist distance matrix, with tcrs defined at the allele level:"

cmd="../../bin/tcrdists --tcrs_file1 allele_tcrs_alpha.txt --database ../../db"

echo
echo $cmd "> test_output2.txt"
$cmd > test_output2.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output2.txt test_output2.txt"
echo
echo $cmd
$cmd

echo
echo "# Compute beta chain tcrdist distance matrix, with tcrs defined at the allele level, comparing tcrs in"
echo "# one file to tcrs in another. Will compute nearest-neighbor (NN) distances for tcrs in file1 with respect"
echo "# to the tcrs in file2"

cmd="../../bin/tcrdists --tcrs_file1 allele_tcrs_beta.txt --tcrs_file2 allele_tcrs_beta_ref.txt --database ../../db"

echo
echo $cmd "> test_output3.txt"
$cmd > test_output3.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output3.txt test_output3.txt"
echo
echo $cmd
$cmd

