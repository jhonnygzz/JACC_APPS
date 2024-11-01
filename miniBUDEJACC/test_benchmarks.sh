#!/bin/bash
# HOW TO RUN:
#chmod +x test_benchmarks.sh
#./test_benchmarks.sh

# Array of values to test
WGSIZES=(128)
PPWIS=(4 8 16 32 64 128)

# Backend to use (cuda, amdgpu, or threads)
BACKEND="cuda"

# Create results directory if it doesn't exist
mkdir -p results

# Create a CSV file for results
RESULTS_FILE="results/benchmark_results_${BACKEND}_$(date +%Y%m%d_%H%M%S).csv"
echo "wgsize,ppwi,kernel_time_ms,avg_time_ms,interactions_per_sec,gflops,gfinsts" > "$RESULTS_FILE"

# Run benchmarks for each combination
for wgsize in "${WGSIZES[@]}"; do
    for ppwi in "${PPWIS[@]}"; do
        echo "Running benchmark with wgsize=$wgsize ppwi=$ppwi"
        
        # Run the benchmark and capture output
        output=$(julia --project=. src/JACCBUDE.jl -w "$wgsize" -p "$ppwi")
        
        # Extract metrics using grep and sed to get just the number
        kernel_time=$(echo "$output" | grep "Kernel time:" | sed 's/.*: *\([0-9.]*\).*/\1/')
        avg_time=$(echo "$output" | grep "Average time:" | sed 's/.*: *\([0-9.]*\).*/\1/')
        interactions=$(echo "$output" | grep "Interactions/s:" | sed 's/.*: *\([0-9.]*\).*/\1/')
        gflops=$(echo "$output" | grep "GFLOP/s:" | sed 's/.*: *\([0-9.]*\).*/\1/')
        gfinsts=$(echo "$output" | grep "GFInst/s:" | sed 's/.*: *\([0-9.]*\).*/\1/')
        
        # Check if we got all values
        if [ -z "$kernel_time" ] || [ -z "$avg_time" ] || [ -z "$interactions" ] || [ -z "$gflops" ] || [ -z "$gfinsts" ]; then
            echo "Warning: Failed to extract some metrics for wgsize=$wgsize ppwi=$ppwi"
            echo "Output was:"
            echo "$output"
            echo "Skipping this combination"
            continue
        fi
        
        # Save results to CSV
        echo "$wgsize,$ppwi,$kernel_time,$avg_time,$interactions,$gflops,$gfinsts" >> "$RESULTS_FILE"
        
        # Print progress
        echo "Completed: wgsize=$wgsize ppwi=$ppwi"
        echo "Last result: $kernel_time,$avg_time,$interactions,$gflops,$gfinsts"
        echo "----------------------------------------"
    done
done

echo "Benchmarking complete! Results saved to $RESULTS_FILE"