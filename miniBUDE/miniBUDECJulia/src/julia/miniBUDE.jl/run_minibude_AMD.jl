#!/usr/bin/env julia

function run_minibude(ppwi::Int, wgsize::Int)
    start_time = time()
    cmd = `julia --project=AMDGPU src/AMDGPU.jl -p $ppwi -w $wgsize --deck ../../../data/bm2`
    
    # Create a unique log file for each run
    log_file = "minibude_p$(ppwi)_w$(wgsize)AMD.log"
    
    println("Running miniBUDE with $ppwi ppwi and $wgsize wgsize")
    println("Logging output to $log_file")
    
    # Run the command and redirect output to log file
    open(log_file, "w") do io
        process = run(pipeline(cmd, stdout=io, stderr=io), wait=false)
        pid = getpid(process)
        
        while process_running(process)
            if (time() - start_time) > 150
                try
                    run(`kill -9 $pid`)
                catch e
                    println("Error killing process $pid: $e")
                end
                
                println("Timeout reached for p=$ppwi, w=$wgsize, deleting log file...")
                close(io)
                rm(log_file)
                return
            end
            sleep(0.1)
        end
    end
    
    # Small delay after completion
    sleep(5)
end

function main()
    # Define parameter ranges
    ppwi = [1, 2, 4, 8, 16, 32, 64, 128]
    wgsize = [1, 2, 4, 8, 16, 32, 64, 128]
    
    # Create results directory if it doesn't exist
    mkpath("minibude_results")
    
    # Change to the directory containing miniBUDE.jl
    cd(dirname(@__FILE__))
    
    # Track start time
    start_time = time()
    
    # Run all combinations
    for p in ppwi
        for w in wgsize
            try
                run_minibude(p, w)
                println("Completed run with p=$p, w=$w")
            catch e
                println("Error running with p=$p, w=$w: $e")
            end
        end
    end
    
    # Calculate and print total runtime
    total_time = time() - start_time
    println("\nAll runs completed in $(round(total_time/60, digits=2)) minutes")
    
end

# Run the main function
main()