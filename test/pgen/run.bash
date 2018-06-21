#!/bin/bash

echo
echo "# Running test calculation on 25 TCRs defined at the V-family level and CDR3 protein sequence..."
cmd="../../bin/pgen -i tcrs.txt -d ../../db/"
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
echo "# Running test calculation on 412 TCRs defined at the V/J allele-level and CDR3 protein sequence..."
cmd="../../bin/pgen -i protseq_tcrs.tsv -d ../../db/ -o test_output2.txt"
echo
echo $cmd
$cmd

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output2.txt test_output2.txt"
echo
echo $cmd
$cmd
echo

echo "# Running test calculation on 412 TCRs defined at the V/J allele-level and CDR3 nucleotide sequence..."
cmd="../../bin/pgen -i nucseq_tcrs.tsv -d ../../db/ -o test_output3.txt"
echo
echo $cmd
$cmd

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_output3.txt test_output3.txt"
echo
echo $cmd
$cmd
