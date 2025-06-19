Okay, that's a standard and very useful technique! You want to define the "genetic space" using the diverse 1000 Genomes Project (1KG) samples and then see where your UK Biobank (UKB) samples fall within that predefined space.

This involves:

1.  Performing PCA **only** on the 1KG samples using the common SNPs.
2.  Saving the results (specifically the SNP loadings or the PCA model information).
3.  Using that saved information to calculate the PC scores for the UKB samples based on their genotypes at the same common SNPs.

Here's how you can do this using **PLINK 2** with PGEN files:

**Assumptions:**

*   You have 1KG data as a single PGEN fileset: `<1kg_prefix>.pgen`, `<1kg_prefix>.pvar[.zst]`, `<1kg_prefix>.psam`.
*   You have UKB data as a single PGEN fileset: `<ukb_prefix>.pgen`, `<ukb_prefix>.pvar[.zst]`, `<ukb_prefix>.psam`.
*   Both are on the **same reference genome build**.
*   You have your `common_snplist.txt`.
*   You want PGEN outputs where applicable.

**Workflow & Commands:**

**Phase 1: Prepare Reference (1KG) and Target (UKB) Datasets**

*   **Step 1a: Extract Common SNPs from 1000 Genomes (Reference)**

    ```bash
    plink2 \
        --pfile <1kg_prefix> \
        --extract common_snplist.txt \
        --make-pgen \
        --out 1kg_common_snps
    ```

*   **Step 1b: LD Prune the 1KG Reference Dataset**
    *   PCA should be run on a set of relatively independent SNPs.

    ```bash
    # Generate list of SNPs to keep after pruning 1KG data
    plink2 \
        --pfile 1kg_common_snps \
        --indep-pairwise 1000kb 50 0.1 \
        # Adjust window(1000kb), step(50), r^2(0.1) as needed
        --out 1kg_pruning

    # Create the pruned 1KG PGEN fileset
    plink2 \
        --pfile 1kg_common_snps \
        --extract 1kg_pruning.prune.in \
        --make-pgen \
        --out 1kg_common_snps_pruned
    ```
    *   **Important:** Note the list of SNPs that survived pruning (`1kg_pruning.prune.in`). You'll need this exact list for the UKB data.

*   **Step 1c: Extract *Pruned* Common SNPs from UK Biobank (Target)**
    *   You need the UKB genotypes **only** for the SNPs used in the final 1KG PCA.
    *   **Crucially, ensure alleles match the 1KG reference.** Mismatched or flipped alleles will invalidate the projection. We use `--ref-allele force` referencing the 1KG `.pvar` file.

    ```bash
    plink2 \
        --pfile <ukb_prefix> \
        --extract 1kg_pruning.prune.in \
        # Use the list of SNPs KEPT after pruning 1KG
        --ref-allele force 1kg_common_snps_pruned.pvar.zst 3 4 \
        # Force UKB alleles (REF=col 4) to match 1KG pruned SNPs (ID=col 3)
        # Check .pvar format if not standard CHROM POS ID REF ALT...
        --make-pgen \
        --out ukb_target_snps
    ```
    *   **Action:** Carefully check the log file (`ukb_target_snps.log`). Note any variants that `--ref-allele force` could not reconcile or had to remove. If many variants fail, investigate the reason (e.g., strand issues, different alleles recorded). The projection accuracy depends on having the correct alleles for the vast majority of SNPs.

**Phase 2: Calculate Reference PCA and Project Target Data**

*   **Step 2a: Calculate PCA on 1KG Reference using `approx`**
    *   The `--pca approx` modifier calculates PCs and saves necessary info for later projection using the `--out` prefix.

    ```bash
    plink2 \
        --pfile 1kg_common_snps_pruned \
        --pca approx 20 \
        # Calculate top 20 PCs for 1KG, adjust number as needed
        --out 1kg_ref_pca
    ```
    *   This generates `1kg_ref_pca.eigenvec` (PCs for 1KG samples), `1kg_ref_pca.eigenval`, and internal information associated with the `1kg_ref_pca` prefix used for projection.

*   **Step 2b: Project UKB Samples onto 1KG PCs**
    *   Use `--pca approx` again, but this time provide the target dataset (`ukb_target_snps`) and tell it to use the reference PCA information calculated in the previous step.

    ```bash
    plink2 \
        --pfile ukb_target_snps \
        --pca approx 20 read-ref-info=1kg_ref_pca \
        # Read the PCA model info from the previous run
        --out ukb_projected_pca
    ```

**Output:**

*   `1kg_ref_pca.eigenvec`: Contains the principal components for the **1000 Genomes samples only**.
*   `1kg_ref_pca.eigenval`: Eigenvalues from the 1000 Genomes PCA.
*   `ukb_projected_pca.eigenvec`: This is the key output file. It contains the **projected principal component scores for the UK Biobank samples**, placing them into the PC space defined by 1000 Genomes. The columns will be FID, IID, PC1_proj, PC2_proj, ..., PC20_proj.

**Important Considerations:**

1.  **Allele Harmonization:** The success of the projection hinges *critically* on Step 1c (`--ref-allele force`). Mismatched alleles between the reference PCA calculation and the target data will lead to incorrect projections. Scrutinize the log from that step.
2.  **SNP Set:** Ensure you use the *exact* same set of SNPs (`1kg_pruning.prune.in`) for both the final 1KG PCA (Step 2a) and the UKB data used for projection (Step 1c -> Step 2b).
3.  **Relatedness in UKB:** Projection itself doesn't require removing relateds from UKB, but if you later use these PCs for association testing within UKB, you'll likely need to account for relatedness there. Consider filtering relateds when creating `ukb_target_snps` if desired.
4.  **Interpretation:** Plotting the projected PCs (`ukb_projected_pca.eigenvec`) alongside the reference PCs (`1kg_ref_pca.eigenvec`), colored by population labels (1KG populations and UKB ancestry groups/subset labels), will show how UKB individuals relate to global populations defined by 1KG.
5.  **Alternative Tools:** While PLINK 2's `approx` method is convenient, tools like `GCTA` (`--pca-project`) are also specifically designed for this task and might offer different options or robustness checks.

This workflow allows you to analyze your UKB samples within the context of the established 1000 Genomes genetic variation structure.