# Thinking:

Let me think through the next steps after the RNA regional filtering completes.

Current state:

- We have 130 million RNA reads from 2,711 cells with barcodes attached
- We're filtering these to HLA and IGH regions using BWA alignment to extended references
- We'll get a filtered FASTQ with reads that map to our regions + unmapped reads

After regional filtering completes, we need to:

1. Align filtered reads to BOTH linear reference and pangenome
2. Call variants from both alignments
3. Compare variant sets
4. Link variants to cell types
5. Visualize results

Let me think through each step in detail...

For alignment:

- Linear: BWA-MEM to the original HLA/IGH reference sequences
- Graph: GraphAligner to the pangenome graphs we created

For variant calling:

- Linear: Could use bcftools, FreeBayes, or DeepVariant
- Graph: vg call on the GAM output

For single-cell analysis:

- We need to create a count matrix from our alignments
- Cluster cells to identify cell types
- Link variants to specific cells/clusters

For visualization:

- UMAP of cells colored by variant presence
- Comparison of variant counts
- Coverage plots showing where pangenome finds extra variants

Key challenges:

- Preserving cell barcode information through all steps
- Ensuring fair comparison between linear and graph alignment
- Having enough reads per cell in our small regions for meaningful analysis

---

# Response:

## Next Steps Roadmap - Post Regional Filtering

Let me think through this systematically:

### Phase 1: Quality Check Filtered Data

Once `align_to_regions.sh` completes, we need to verify:

1. **How many reads survived filtering?**
    
    - Total reads that mapped to HLA region
    - Total reads that mapped to IGH region
    - Unmapped reads we're keeping (potential alt haplotype reads)
    - Distribution across cells - do we still have good coverage?
2. **Key concern**: With only ~5Mb regions out of 3Gb genome, we might have very few reads
    
    - Expected: Maybe 0.1-0.5% of reads will map to our regions
    - From 130M reads → possibly 130K-650K filtered reads
    - Across 2,711 cells → 50-250 reads per cell in regions
    - This might be too sparse for robust single-cell analysis

### Phase 2: Parallel Alignment Strategy

**Linear Reference Alignment:**

- Use BWA-MEM on filtered reads → original `hla_hg38.fa` and `igh_hg38.fa`
- Output: BAM files with preserved CB tags
- This gives us baseline variant calls

**Pangenome Graph Alignment:**

- Use GraphAligner on same filtered reads → `HLA_region_full.vg` and `IGH_region_full.vg`
- Need to preserve cell barcodes through alignment
- Output: GAM format (Graph Alignment Map)
- Convert GAM → BAM for variant calling

### Phase 3: Variant Calling

**Critical consideration**: With sparse data per cell, we need to:

- Call variants at the **population level** (all cells together) first
- Then genotype individual cells
- Cannot reliably call variants per individual cell with ~100 reads

**Linear variants:**

- Use bcftools mpileup + call (simple, fast)
- Or FreeBayes for better indel calling
- Generate VCF with all variants found

**Graph variants:**

- Use `vg call` on GAM file
- Extract novel variants not in linear reference
- Generate comparable VCF

### Phase 4: Single-Cell Analysis Approach

Given sparse data, we need a modified strategy:

1. **Pseudo-bulk by cell type first**
    
    - Use existing cell type annotations from the h5 file
    - Or do quick clustering based on full transcriptome (if available)
    - Pool reads by cell type for better coverage
2. **Variant-cell matrix creation**
    
    - For each variant: count supporting reads per cell
    - Create sparse matrix: variants × cells
    - Most entries will be 0 or 1 read
3. **Statistical approach**
    
    - Can't do per-cell genotyping reliably
    - Instead: "Which cell types are enriched for variant X?"
    - Use aggregation and statistical tests

### Phase 5: Integration with ATAC Data

Using the fragments file:

1. Identify which variants fall in open chromatin regions
2. Check if novel pangenome variants are in accessible regions
3. Link to cell-type specific accessibility

### Phase 6: Key Comparisons and Visualizations

**Primary comparison:**

- Total variants: Linear vs Pangenome
- Novel variants found only in pangenome
- Which are in coding regions vs regulatory regions?

**Cell type associations:**

- Which variants show cell-type specific patterns?
- Even with sparse data, strong effects might be visible

**Critical plots:**

1. Venn diagram of variants (linear vs pangenome)
2. Coverage comparison plot for a example region
3. Cell type enrichment heatmap for novel variants
4. ATAC accessibility overlay on variant positions

### Key Decision Points:

1. **If filtered data is too sparse** (< 50 reads/cell):
    
    - Might need to focus on pseudo-bulk analysis only
    - Or demonstrate proof-of-concept on most covered region
2. **If no novel variants found**:
    
    - Check if graph alignment is working properly
    - Might need to adjust variant calling parameters
    - Could focus on improved mapping quality instead
3. **Computational considerations**:
    
    - GraphAligner might be slow on 100k+ reads
    - May need to subsample or run in batches
    - GAM files can be large

The overall goal remains: **Show that pangenome alignment finds variants missed by linear reference, particularly in the highly polymorphic HLA region.**

Copy