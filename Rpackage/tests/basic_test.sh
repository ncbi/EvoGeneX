#!/bin/bash

# Create results directory if it doesn't exists
if [ ! -d "tests/results/" ]; then
    mkdir "tests/results/"
fi

# Run the different versions
Rscript tests/run_local.R
Rscript tests/run_release.R

# Print results
echo "Diff between brown slow:"
diff tests/results/local_brown_slow.csv tests/results/release_brown_slow.csv
echo "Diff between brown fast:"
diff tests/results/local_brown_fast.csv tests/results/release_brown_fast.csv
echo "Diff between single regime evog slow:"
diff tests/results/local_single_evog_slow.csv tests/results/release_single_evog_slow.csv
echo "Diff between single regime evog fast:"
diff tests/results/local_single_evog_fast.csv tests/results/release_single_evog_fast.csv
echo "Diff between two regime evog slow:"
diff tests/results/local_two_evog_slow.csv tests/results/release_two_evog_slow.csv
echo "Diff between two regime evog fast:"
diff tests/results/local_two_evog_fast.csv tests/results/release_two_evog_fast.csv
echo "Done"