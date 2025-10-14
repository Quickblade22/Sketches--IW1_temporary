#!/bin/bash
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=16          # 16 CPU cores per task
#SBATCH --mem=64g                   # 64GB memory allocation
#SBATCH --partition=rleap_cpu       # Partition (queue) to use
#SBATCH --output=output.txt  # Output log file
#SBATCH --error=error.err



# Run the Python script inside the container using exec with proper directory navigation
echo "Running Python script in container using exec..."
apptainer exec \
    --bind .:/work/rleap1/aaditya_mehta/Sketches--IW1_temporary \
    mycontainer.sif \
    sh -c "cd src && make clean && python3 running.py /work/rleap1/aaditya_mehta/Sketches--IW1_temporary/Adventure.bin 10000 5000000 --compile-first --output-file=output.txt"

# Check if the execution was successful
if [ $? -eq 0 ]; then
    echo "Container execution completed successfully!"
else
    echo "Container execution failed."
    exit 1
fi

echo "Script completed successfully!"
