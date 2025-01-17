#!/bin/bash

test_size='0.25'
genotype='deeprvat'
emb_state='--shuffled_emb'
# version='v1NEWsplit'

while read trait
do
    sbatch submit_funcrvp_training.sh $trait pops_exp $genotype $test_size v1NEWsplit
    # sbatch -p urgent submit_funcrvp_training.sh $trait omics_pops_exp $genotype $test_size v1NEWsplit
    # sbatch -p urgent submit_funcrvp_training.sh $trait string $genotype $test_size v1NEWsplit
    # sbatch -p urgent submit_funcrvp_training.sh $trait string_exp $genotype $test_size v1NEWsplit
    # sbatch -p urgent submit_funcrvp_training.sh $trait omics_pops $genotype $test_size v1NEWsplit_randEmb $emb_state '0'
    # sbatch submit_funcrvp_training.sh $trait enformer_small deeprvat $test_size v1NEWsplit
done < phenotype_list_41.txt
