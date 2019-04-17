#!/bin/bash

for job in *.job; do
    [ -e "$job" ] || continue
    echo "sbatch '${job}'"
    sbatch "${job}"
done
