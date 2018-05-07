#!/bin/bash

echo
echo "## ==============================================="
echo "## test pgen"
echo "cd pgen/"
cd pgen/
./run.bash

echo
echo "## ==============================================="
echo "## test tcrdists"
echo "cd ../tcrdists/"
cd ../tcrdists/
./run.bash

echo
echo "## ==============================================="
echo "## test correlations"
echo "cd ../correlations/"
cd ../correlations/
./run.bash

echo
echo "## ==============================================="
echo "## test neighbors"
echo "cd ../neighbors/"
cd ../neighbors/
./run.bash
echo
cd ../


