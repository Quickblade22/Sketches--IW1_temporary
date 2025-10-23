#!/bin/bash
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --cpus-per-task=32          # 16 CPU cores per task
#SBATCH --mem=128g                   # 64GB memory allocation
#SBATCH --partition=rleap_cpu       # Partition (queue) to use
#SBATCH --output=output3.txt  # Output log file
#SBATCH --error=error3.err



# Run the Python script inside the container using exec with proper directory navigation
echo "Running Python script in container using exec..."
apptainer exec \
    --bind .:/work/rleap1/aaditya_mehta/Sketches--IW1_temporary \
    mycontainer2.sif \
    sh -c "cd src && make clean && echo cpu 32 138 50000 700000 >> output2.txt && python3 running.py /work/rleap1/aaditya_mehta/Sketches--IW1_temporary/Adventure.bin 50000 700000 --additional-args=' --features 4' --compile-first --output-file=output3.txt"

# Check if the execution was successful
if [ $? -eq 0 ]; then
    echo "Container execution completed successfully!"
else
    echo "Container execution failed."
    exit 1
fi

echo "Script completed successfully!"
