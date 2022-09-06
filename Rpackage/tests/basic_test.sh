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
diff tests/results/local_brown_slow.csv tests/curr_brown_slow.csv
echo "Diff between brown fast:"
diff tests/results/local_brown_fast.csv tests/curr_brown_fast.csv
echo "Diff between evog slow:"
diff tests/results/local_evog_slow.csv tests/curr_evog_slow.csv
echo "Diff between evog fast:"
diff tests/results/local_evog_fast.csv tests/curr_evog_fast.csv
echo "Done"