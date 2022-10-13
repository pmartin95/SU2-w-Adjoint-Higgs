#!/bin/bash
sbatch -p RM-shared -t 12:00:00 \
--mail-type=ALL --mail-user=martip8@rpi.edu \
-A tra160032p -o output/job1.output -e job1.err -J "775p0" \
--ntasks=1 ./job1.job 
