# Single-Cell Multi-omics Pangenome Analysis Project - Current State Documentation

## Project Overview

### Core Hypothesis

The human pangenome reference captures genetic variation missed by the traditional linear reference (hg38), particularly in highly polymorphic regions like HLA. By aligning single-cell multi-omics data to both references, we can:

1. Identify additional variants unique to the pangenome
2. Link these variants to specific cell types using single-cell resolution
3. Demonstrate functional relevance using chromatin accessibility data

### Biological Significance

- **HLA region** (chr6:28.5-33.5Mb): Highly polymorphic, critical for immune function
- **IGH region** (chr14:105.6-107.3Mb): Immunoglobulin heavy chain, less polymorphic control
- **10x Multiome**: Paired RNA + ATAC from same cells enables linking variants to gene expression AND chromatin accessibility

## Current Directory Structure

```
~/Omics-UBIC/pangenome_project/
├── data/
│   ├── pbmc_granulocyte_sorted_3k/
│   │   ├── gex/
│   │   │   ├── *_S1_L003_R1_001.fastq.gz  # Barcode + UMI (28bp)
│   │   │   ├── *_S1_L003_R2_001.fastq.gz  # cDNA (90bp)
│   │   │   ├── *_S1_L004_R1_001.fastq.gz  # Multiple lanes
│   │   │   └── *_S1_L004_R2_001.fastq.gz  # (8 files total)
│   │   └── atac/
│   │       ├── *_S12_L001-004_R1_001.fastq.gz  # Genomic DNA (50bp)
│   │       ├── *_S12_L001-004_R2_001.fastq.gz  # Barcode (24bp)
│   │       └── *_S12_L001-004_R3_001.fastq.gz  # Genomic DNA (50bp)
│   │                                            # (16 files total)
│   ├── filtered_reads/
│   │   ├── gex_R2_all.fastq          # 39GB concatenated RNA reads
│   │   ├── atac_R1_all.fastq         # 12GB concatenated ATAC R1
│   │   ├── atac_R3_all.fastq         # 12GB concatenated ATAC R3
│   │   └── gex_hla_mapped.fastq.sam  # 24.5GB partial BWA output
│   ├── filter_reads_by_region.sh      # Original flawed script
│   └── filter_reads_by_region_v2.sh   # Improved but still flawed
├── references/
│   ├── hla_hg38.fa     # chr6:28500000-33500000 (5Mb)
│   ├── igh_hg38.fa     # chr14:105600000-107300000 (1.7Mb)
│   └── *.fa.{amb,ann,bwt,pac,sa}  # BWA index files
├── graphs/
│   ├── HLA_region.vg   # Pangenome graph for HLA
│   └── IGH_region.vg   # Pangenome graph for IGH
└── results/
    └── (empty)
```

Copy Code

## What Has Been Done

### 1. Environment Setup ✓

```
conda create -n pangenome python=3.9
conda install -c bioconda bwa samtools seqtk
# Still need: GraphAligner, vg tools, deepvariant
```

Copy Code

### 2. Data Acquisition ✓

- Downloaded 3k PBMC multiome dataset (22.4GB tar)
- Extracted FASTQ files
- Downloaded and properly formatted regional references

### 3. Initial Processing Attempt (FLAWED)

**What we did:**

```
# Concatenated all lanes
zcat gex/*_R2_001.fastq.gz > gex_R2_all.fastq  # 39GB

# Attempted to filter by alignment to regional references
bwa mem hla_hg38.fa gex_R2_all.fastq > mapped.sam  # Crashed after ~3hrs
```

Copy Code

**Why this was wrong:**

1. Lost cell barcode information by only using R2
2. Pre-filtering based on linear alignment biases against pangenome
3. No preservation of UMI for deduplication
4. Trying to convert aligned reads back to FASTQ loses alignment info

## File Format Details

### 10x Genomics Read Structure

**RNA (Gene Expression)**

- R1: 28bp = 16bp Cell Barcode + 12bp UMI
- R2: 90bp = cDNA insert (actual sequence)

**ATAC**

- R1: 50bp = Genomic DNA (forward)
- R2: 24bp = 16bp Cell Barcode + 8bp fixed
- R3: 50bp = Genomic DNA (reverse)

### Critical Metadata to Preserve

1. **Cell Barcode (CB)**: Links reads to cells
2. **UMI (UB)**: For RNA deduplication
3. **Read pairing**: For ATAC fragments
4. **Original read ID**: For troubleshooting

## The Corrected Approach

### Core Principle

Align THE SAME reads to both linear reference and pangenome to enable fair comparison.

### Implementation Strategy

```
Phase 1A: Prepare Reads with Metadata
1. Extract barcodes from R1 (RNA) or R2 (ATAC)
2. Validate against whitelist
3. Add barcodes to read headers
4. Create alignment-ready FASTQs

Phase 1B: Smart Read Selection
Option 1: K-mer based filtering
- Extract k-mers from BOTH references
- Keep reads with matching k-mers
  
Option 2: Permissive BWA filtering
- Align to extended regions (±1Mb)
- Keep mapped AND unmapped reads
- Include low-quality mappings

Phase 1C: Parallel Alignment
Linear Pipeline:
- BWA-MEM with appropriate parameters
- Preserve CB/UB tags in output BAM

Pangenome Pipeline:
- GraphAligner with same reads
- Convert GAM to BAM format
- Preserve graph-specific info

Phase 1D: Comparison Metrics
- Total/mapped/uniquely mapped reads
- Reads rescued by pangenome
- Mapping quality distributions
```

Copy Code

## Missing Components

### Software Not Yet Installed

1. **GraphAligner**: For pangenome alignment
2. **vg tools**: For graph manipulation and variant calling
3. **deepvariant**: For linear variant calling
4. **Cell Ranger ARC**: For proper 10x processing (optional)

### Data Not Yet Generated

1. Proper barcode-tagged FASTQ files
2. K-mer databases from references
3. Completed alignments (BAM/GAM files)
4. Variant calls (VCF files)

## Next Immediate Steps

1. **Install missing tools:**

```
# GraphAligner
git clone https://github.com/maickrau/GraphAligner
cd GraphAligner && make

# vg tools
conda install -c bioconda vg

# Get barcode whitelist
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
```

Copy Code

2. **Create barcode-aware FASTQ files:**

```
# Need to implement:
- Read barcode whitelist
- Process R1/R2 pairs
- Add CB:Z:BARCODE tag to headers
- Output new FASTQs
```

Copy Code

3. **Implement fair filtering strategy:**

- Either k-mer based or permissive alignment
- Must include potential pangenome-specific reads

## Known Issues to Avoid

1. **Don't pre-filter based on linear alignment** - biases results
2. **Don't lose barcode information** - needed for single-cell analysis
3. **Don't delete alignment files** - needed for variant calling
4. **Don't process all 39GB at once** - consider batching

## Expected Outcomes

- Pangenome should identify 5-20% more variants in HLA region
- These variants should be enriched in:
    - Alternate haplotypes
    - Regions with poor linear reference coverage
    - Cell-type specific accessible chromatin
- Final deliverable: Visualization showing variant discovery improvement and cell-type associations

## File Size Estimates

- Input: ~17GB extracted FASTQs
- Filtered reads: ~1-2GB (targeting specific regions)
- Alignments: ~2-4GB per reference
- Variant calls: ~10-100MB VCFs

This documentation should allow seamless continuation of the project with full technical context.

Copy