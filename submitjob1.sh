#!/bin/bash
sbatch -p RM-shared -t 3:00:00 \
--mail-type=ALL --mail-user=martip8@rpi.edu \
-A phy200008p -o output/job1.output -e job1.err -J "775p0" \
--ntasks=61 ./job1.job 
