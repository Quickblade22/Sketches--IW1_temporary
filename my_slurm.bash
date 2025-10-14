#!/bin/bash
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=16          # 16 CPU cores per task
#SBATCH --mem=64g                   # 64GB memory allocation
#SBATCH --partition=rleap_cpu  # Partition (queue) to use
#SBATCH --output=/work/rleap1/aaditya_mehta/%container_run.txt  # Output log file



# Change to the directory containing your scripts
cd /work/rleap1/aaditya_mehta/Sketches--IW1_temporary  # UPDATE THIS PATH

# Run the Python script that handles container execution
python3 running_container.py && echo "Job finished at: $(date)"