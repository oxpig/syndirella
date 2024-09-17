#!/bin/bash

# Check if the user provided JOB_LIST, JOB_DIR, and LOG_DIR as arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 job_list_file job_directory log_directory"
    exit 1
fi

# Use the first argument as the JOB_LIST file
JOB_LIST="$1"

# Use the second argument as the directory containing the job scripts
JOB_DIR="$2"

# Use the third argument as the directory for log files
LOG_DIR="$3"

# Ensure the log directory exists
if [ ! -d "$LOG_DIR" ]; then
    echo "Log directory $LOG_DIR does not exist."
    exit 1
fi

# Extract the base name of the job list file (without the directory path)
JOB_LIST_BASENAME=$(basename "$JOB_LIST")

# Replace the file extension with .log and use the log directory
LOG_FILE="${LOG_DIR}/${JOB_LIST_BASENAME%.*}.log"

# Clear the log file if it already exists
> "$LOG_FILE"

# Loop through each job script in the list
while IFS= read -r job_script; do
    # Full path to the job script
    full_path="$JOB_DIR/$job_script"

    if [[ -f "$full_path" ]]; then
        # Submit the job and capture the job ID
        job_id=$(sbatch "$full_path" | awk '{print $4}')

        # Log the job ID with the corresponding script name
        echo "$job_script: $job_id" >> "$LOG_FILE"

        # Wait for 30 seconds before submitting the next job
        sleep 30
    else
        echo "Job script $full_path not found!" >> "$LOG_FILE"
    fi
done < "$JOB_LIST"
