#!/bin/bash

experiments_directory="experiments"
mkdir -p "$experiments_directory"

declare -a settings=(
    '-r 1 -c 5 -m 5 -p 5 -l 1'
    '-r 2 -c 5 -m 5 -p 5 -l 1'
    '-r 4 -c 5 -m 5 -p 5 -l 1'
    '-r 5 -c 5 -m 5 -p 5 -l 1'
    '-r 1 -c 10 -m 5 -p 5 -l 1'
    '-r 1 -c 20 -m 5 -p 5 -l 1'
    '-r 1 -c 40 -m 5 -p 5 -l 1'

)

evaluations=10000

declare -a datasetfiles=("banknotes.txt" "diabetes.txt")
declare -a methods=(1 2)

for method in "${methods[@]}"; do
    for datasetfile in "${datasetfiles[@]}"; do
        for setting in "${settings[@]}"; do
            # strip the file extension from the dataset name
            dataset_name=$(echo "$datasetfile" | sed 's/\.txt$//')

            settings_directory="$experiments_directory/t${method}d${dataset_name}${setting//[- ]}"
            mkdir -p "$settings_directory"

            for i in {1..30}; do
                experiment_directory="$settings_directory/experiment_$i"
                mkdir -p "$experiment_directory"

                population_size=$(echo "$setting" | sed -n 's/.*-p \([0-9]*\).*/\1/p')
                generations=$((evaluations / population_size))

                args=("-g" "$generations" "$setting")
                args=("-f" "$experiment_directory/best_chromosome.chr" "${args[@]}")
                args=("-d" "$datasetfile" "${args[@]}")
                args=("-t" "$method" "${args[@]}")

                echo "Running experiment $i with the following arguments:"
                echo "-d $datasetfile -t $method $setting"

    #            run_command="./cmake-build-debug/bin"
                run_command="./cmake-build-debug/BIN_projekt.exe"
                for arg in "${args[@]}"; do
                    run_command="$run_command $arg"
                done

    #            echo "$run_command"
                ($run_command > "$experiment_directory/output.txt") &

                # limit the number of parallel jobs
                if ((i % 10 == 0)); then
                    wait
                fi
            done
        done
    done
done