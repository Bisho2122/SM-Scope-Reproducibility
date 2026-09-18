#!/bin/bash

# Declare arrays

# declare -a metric_feat=("correlation" "euclidean")
# declare -a weight_factors=("1511" "5511" "1515")
# declare -a dim_out=(10 30 50)
# declare -a k_feat=(15 50)

declare -a metric_feat=("euclidean")
declare -a weight_factors=("5511")
declare -a dim_out=(50)
declare -a k_feat=(15)

job_name="SG"

# Nested for loops
for i in "${dim_out[@]}"; do
  for j in "${k_feat[@]}"; do
    for k in "${metric_feat[@]}"; do
      for l in "${weight_factors[@]}"; do
            err_filename="err/${job_name}_${i}_${j}_${k}_${l}.err"
            out_filename="out/${job_name}_${i}_${j}_${k}_${l}.out"
            sbatch -e $err_filename -o $out_filename SpatialGlue_job.sh $i $j $k $l
            sleep 1
      done
    done
  done
done