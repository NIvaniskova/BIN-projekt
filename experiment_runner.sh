#!/bin/bash

experiments_directory="experiments"
mkdir -p "$experiments_directory"

declare -a settings=(
    '-r 1 -c 15 -m 5 -p 5 -l 15'
    #'-r 1 -c 15 -m 5 -p 5 -l 3'
    '-r 3 -c 5 -m 5 -p 5 -l 15'
    #'-r 3 -c 5 -m 5 -p 5 -l 3'
    '-r 1 -c 15 -m 5 -p 10 -l 15'
    #'-r 1 -c 15 -m 5 -p 10 -l 3'
    '-r 3 -c 5 -m 5 -p 10 -l 15'
    #'-r 3 -c 5 -m 5 -p 10 -l 3'
    # todo
    # přidej další nastavení
)

evaluations=10000
j = 0
for setting in "${settings[@]}"; do
    settings_directory="$experiments_directory/$j"
    mkdir -p "$settings_directory"

    for i in {1..30}; do
        experiment_directory="$settings_directory/experiment_$i"
        mkdir -p "$experiment_directory"

        population_size=$(echo "$setting" | sed -n 's/.*-p \([0-9]*\).*/\1/p')
        generations=$((evaluations / population_size))

        args=("-g" "$generations" $setting)

        # todo
        # DONE v cgp.h přidej přepínač toho, kam se ukládá nejlepší chromozom
        # DONE a v main.cpp oddělej tu podmínku, že se ukládá jenom někdy, ulož ho vždycky
        # a nastav tu cetu na ten argument, ja mu tady dávám přepínač f,
        # tak pak dej jenom case f... v tom cgp.h
        args=("-f" "$experiment_directory/best_chromosome.chr" $args)

        echo "Running experiment $i with arguments:"
        echo "${args[@]}"
        echo "./cmake-build-debug/BIN_projekt.exe" "${args[@]}" # > "$experiment_directory/output.txt"
        # (./cmake-build-debug/BIN_projekt.exe "${args[@]}" > "$experiment_directory/output.txt") &

        # limit the number of parallel jobs
        if ((i % 10 == 0)); then
            wait
        fi
    done
    j++
done
