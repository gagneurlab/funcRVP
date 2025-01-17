#!/bin/bash

# Define a embedding list and version number
emb_list=("omics_pops")
version="v111cov"

while read trait
do
    for emb in "${emb_list[@]}"
    do
        for genotype in deepRVAT 
        do
            sbatch submit_hyperopt_job.sh $trait $emb $genotype $version
            # sbatch submit_hyperopt_mean_job.sh $trait $emb $genotype v5cov
        done
    done
done < phenotype_list_41.txt