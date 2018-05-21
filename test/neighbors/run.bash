#!/bin/bash

echo
echo "## This just computes the neighbors:"

cmd="../../bin/neighbors --similarity_mode 2 --matrix1 example1_matrix.txt --nbr_pval_threshold 1e-4 --subject_bias example1_subject_bias.txt --feature_mask example1_features.txt --database ../../db"

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
echo "## This command uses the neighbor relationships for DBSCAN clustering"
cmd="../../bin/neighbors --similarity_mode 2 --matrix1 example1_matrix.txt --nbr_pval_threshold 1e-4 --subject_bias example1_subject_bias.txt --feature_mask example1_features.txt --cluster --min_core_nbrcount 4 --min_cluster_size 4 --database ../../db"

echo
echo $cmd "> test_output2.txt"
$cmd > test_output2.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
echo
echo "diff <(grep ^cluster: expected_output2.txt) <(grep ^cluster: test_output2.txt)"
diff <(grep ^cluster: expected_output2.txt) <(grep ^cluster: test_output2.txt)

echo
echo "## This command uses the neighbor relationships for DBSCAN clustering and reads DBSCAN "
echo "## parameters from a database file."
cmd="../../bin/neighbors --similarity_mode 2 --matrix1 example1_matrix.txt --subject_bias example1_subject_bias.txt --feature_mask example1_features.txt --cluster --database ../../db"

echo
echo $cmd "> test_output3.txt"
$cmd > test_output3.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
echo
echo "diff <(grep ^cluster: expected_output3.txt) <(grep ^cluster: test_output3.txt)"
diff <(grep ^cluster: expected_output3.txt) <(grep ^cluster: test_output3.txt)

echo
echo "## This command performs DBSCAN clustering using the TCRdist measure, using default neighbor threshold"
echo "## and precomputed DBSCAN paramters read from a database file."
cmd="../../bin/neighbors --similarity_mode 3 --tcrs example4_tcrs.txt --cluster --database ../../db"

echo
echo $cmd "> test_output4.txt"
$cmd > test_output4.txt

echo
echo "# Comparing output to expected output (this command should produce no output):"
echo
echo "diff expected_output4.txt test_output4.txt"
diff expected_output4.txt test_output4.txt

