#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate coloc_reporter
cd /home/ec2-user/software/COLOC-reporter
Rscript LAVA_EAS.R
