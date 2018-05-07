#!/bin/bash

echo "==============================================="
echo "test pgen"
cd pgen/
./run.bash

echo "==============================================="
echo "test tcrdists"
cd ../tcrdists/
./run.bash

echo "==============================================="
echo "test correlations"
cd ../correlations/
./run.bash

echo "==============================================="
echo "test neighbors"
cd ../neighbors/
./run.bash
echo "==============================================="
cd ../


