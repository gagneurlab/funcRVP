#!/bin/bash

version="v1NEWsplitHO"
# ukb_ver="inv_train"
test_size='0.25'

while read trait
do
    for emb in omics_pops
    do
        for genotype in deepRVAT
        do
            sbatch submit_funcrvp_hyperopt.sh $trait $emb $genotype $version $test_size
            # sbatch submit_funcrvp_hyperopt.sh $trait $emb $genotype $version $test_size
            # sbatch submit_funcrvp_hyperopt.sh $trait $emb $genotype $version $ukb_ver
            # sbatch submit_funcrvp_hyperopt.sh $trait $emb $genotype v5cov
        done
    done
done < phenotype_list_41.txt
