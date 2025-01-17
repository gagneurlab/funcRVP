#!/bin/bash

while read trait
do
    sbatch -p urgent submit_arch_job.sh $trait pops deepRVAT v112arch
    sbatch -p urgent submit_arch_job.sh $trait omics deepRVAT v113arch
    sbatch -p urgent submit_arch_job.sh $trait omics_pops_exp deepRVAT v114arch
done < phenotype_list_41.txt





# Define a embedding list and version number
# emb_list=("pops")
# version="v112arch"


# while read trait
# do
#     for emb in "${emb_list[@]}"
#     do
#         for genotype in deepRVAT 
#         do
#             sbatch submit_arch_job.sh $trait $emb $genotype $version
#             # sbatch submit_arch_job_shuffledpheno.sh $trait $emb $genotype $version 1234
#             # sbatch submit_arch_job_shuffledpheno.sh $trait $emb $genotype $version 111
#             # sbatch submit_arch_job_shuffledpheno.sh $trait $emb $genotype $version 225
#         done
#     done
# done < phenotype_list_41.txt