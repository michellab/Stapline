#!/bin/bash

eq_jobid=$(sbatch --parsable min_eq.sh)

echo "Submitted min_eq.sh with Job ID $eq_jobid"

sbatch --dependency=afterok:$eq_jobid prod.sh
echo "Submitted prod.sh (dependent on min_eq.sh)"

sbatch --dependency=afterok:$eq_jobid prod2.sh
echo "Submitted prod2.sh (dependent on min_eq.sh)"

sbatch --dependency=afterok:$eq_jobid prod3.sh
echo "Submitted prod3.sh (dependent on min_eq.sh)"
