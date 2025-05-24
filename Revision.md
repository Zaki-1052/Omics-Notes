# Thinking:

Let me carefully think through this project step by step.

Starting point:

- We have 10x multiome data (RNA + ATAC from the same cells)
- We have a pangenome graph for two regions (HLA and IGH)
- We have linear reference sequences for the same regions

The core hypothesis:

- The pangenome contains sequences missing from the linear reference
- When we align reads to the pangenome vs linear reference, we'll find MORE variants
- These additional variants can be linked to specific cell types

The workflow should be:

1. Take ALL the reads from the 10x data
2. Split into two parallel pipelines:
    - Pipeline A: Align to LINEAR reference → call variants
    - Pipeline B: Align to PANGENOME graph → call variants
3. Compare the variants found in each pipeline
4. Use the cell barcode information to link variants to specific cell types
5. Show that pangenome finds more variants and these are biologically meaningful

Current issue with the script:

- We're filtering reads BEFORE aligning to the pangenome
- This defeats the purpose! We want to see what the pangenome can capture that the linear reference might miss
- We should align ALL reads to both references

The filtering approach only makes sense if:

- We have computational limitations
- We want to focus on specific regions due to time constraints

But even then, we should:

1. Filter to get reads that MIGHT map to our regions (using a permissive threshold)
2. Align these to BOTH linear and pangenome references
3. Compare variant calls

The biological significance:

- ATAC data is crucial - it shows which regions are accessible in each cell type
- Variants in accessible regions are more likely to be functional
- By combining RNA + ATAC + variants, we can show:
    - Which variants are in open chromatin (ATAC)
    - Which affect gene expression (RNA)
    - Which cell types are affected

---

# Response:

Looking at the roadmap and project goals, I need to reconsider the entire workflow. Let me think through this carefully:

## What We're Starting With:

- **10x Multiome data**: Paired RNA + ATAC reads from the same individual cells (~3,000 cells)
- **Two reference types**:
    1. Linear reference (traditional hg38) for HLA and IGH regions
    2. Pangenome graph containing multiple haplotypes for the same regions
- **Cell barcodes**: Each read is tagged with which cell it came from

## The Core Goal:

**Demonstrate that the pangenome captures more genetic variation than the linear reference, and link these variants to specific cell types.**

## The Correct Workflow:

### Current Script Problem:

The script is doing **pre-filtering** - only keeping reads that map to our regions on the linear reference, then planning to realign these to the pangenome. **This is backwards!**

### Why This Is Wrong:

1. **We're biasing toward the linear reference** - if a read doesn't map well to linear reference but would map to an alternate haplotype in the pangenome, we'd miss it
2. **We lose the comparison** - we need to align the SAME reads to both references to fairly compare variant discovery

### What We Should Do:

```
Original Reads (ALL cells, ALL genomic regions)
           |
           ├─→ Align to Linear Reference (BWA-MEM)
           |          ↓
           |    Call Variants (DeepVariant)
           |          ↓
           |    Variants Set A
           |
           └─→ Align to Pangenome Graph (GraphAligner)  
                      ↓
                Call Variants (vg call)
                      ↓
                Variants Set B
                      
Compare: Set B should contain Set A + additional variants
```

Copy Code

### The Biologically Correct Approach:

1. **Option A: Full Genome Alignment (Ideal but computationally intensive)**
    
    - Align ALL reads to full linear genome
    - Align ALL reads to full pangenome
    - Focus analysis on HLA/IGH regions
2. **Option B: Smart Pre-filtering (Practical compromise)**
    
    - Use a **permissive** filter to get reads that MIGHT map to our regions
    - Include unmapped reads and reads mapping nearby
    - Then align this subset to both references

### For Option B, the filtering should be:

```
# Get reads mapping to extended regions (± 1Mb buffer)
samtools view -b all_reads.bam chr6:27500000-34500000 chr14:104600000-108300000 -o region_candidates.bam

# ALSO include unmapped reads (they might map to pangenome alternate sequences)
samtools view -b -f 4 all_reads.bam -o unmapped.bam

# Merge and convert to FASTQ for realignment
samtools merge candidates.bam region_candidates.bam unmapped.bam
```

Copy Code

## Why ATAC Data Matters:

1. **Chromatin accessibility** = which DNA regions are "open" in each cell type
2. **Functional variants** are more likely in open chromatin
3. We can show:
    - Variant X is in open chromatin in T cells but not B cells
    - Therefore, Variant X likely affects T cell function
    - The pangenome found Variant X but linear reference missed it

## The Final Biological Story:

"We found [N] additional variants using the pangenome that were missed by the linear reference. Of these, [M] variants were in cell-type-specific open chromatin regions, suggesting functional relevance. For example, variant rs##### was found in open chromatin specifically in CD8+ T cells and correlated with increased expression of gene GENE1."

Should we revise our approach to properly compare linear vs pangenome alignment on the same set of reads?

Copy