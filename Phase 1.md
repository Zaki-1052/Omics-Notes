# Thinking:

The user is asking me to think very carefully about the implementation approach for Phase 1. Let me reconsider what we have and what we need to do.

Starting materials:

1. 10x multiome FASTQ files (RNA + ATAC)
    - RNA: R1 (barcode+UMI), R2 (cDNA)
    - ATAC: R1 (genomic), R2 (barcode), R3 (genomic)
2. Linear reference (hla_hg38.fa, igh_hg38.fa) - just the regions
3. Pangenome graph (HLA_region.vg, IGH_region.vg)

The key insight: We need to compare variant discovery between linear and graph alignment using THE SAME reads.

Current problem: We're pre-filtering based on linear alignment, which biases against pangenome-specific sequences.

Better approach:

1. Since we only have regional references (not full genome), we need a smart filtering strategy
2. We should cast a wide net to capture reads that might map to alternate haplotypes
3. Include unmapped reads in case they map to pangenome-specific sequences

The biological workflow should preserve:

- Cell barcode information throughout
- Paired-end relationships for ATAC
- UMI information for deduplication

Let me think through the detailed implementation...

---

# Response:

## Detailed Implementation Plan for Phase 1

After thinking deeply about this, here's the correct approach:

### Starting Materials Assessment:

```
Data we have:
├── 10x Multiome FASTQs (~3k cells)
│   ├── RNA: R1 (barcode+UMI), R2 (cDNA sequence)
│   └── ATAC: R1 (genomic), R2 (barcode), R3 (genomic)
├── Linear references (region-specific)
│   ├── hla_hg38.fa (chr6:28.5-33.5Mb)
│   └── igh_hg38.fa (chr14:105.6-107.3Mb)
└── Pangenome graphs
    ├── HLA_region.vg
    └── IGH_region.vg
```

Copy Code

### Key Realization:

Since we only have **regional references** (not whole genome), we can't do a traditional "align to whole genome first" approach. We need a different strategy.

### Proposed Implementation:

#### Step 1: Intelligent Read Selection

```
Problem: We need reads that MIGHT map to our regions, without biasing toward linear reference

Solution: Multi-stage filtering
1. Quick k-mer based filtering
   - Extract k-mers from BOTH linear reference AND pangenome graph
   - Keep reads sharing k-mers with either reference
   - This captures reads that might map to alternate haplotypes

2. Preserve cell information
   - Keep barcode-UMI mappings
   - Maintain read pairing for ATAC

Pseudocode:
create_kmer_database(hla_linear, hla_pangenome, k=21) -> hla_kmers
create_kmer_database(igh_linear, igh_pangenome, k=21) -> igh_kmers

for read in all_reads:
    if shares_kmers(read, hla_kmers, min_matches=2):
        add_to_hla_candidate_reads(read, preserve_barcode=True)
    if shares_kmers(read, igh_kmers, min_matches=2):
        add_to_igh_candidate_reads(read, preserve_barcode=True)
```

Copy Code

#### Step 2: Parallel Alignment Strategy

```
For each region (HLA, IGH):
├── Linear Alignment Pipeline
│   ├── RNA reads:
│   │   └── bwa mem -M (mark supplementary) 
│   │       → SAM with CB/UB tags preserved
│   ├── ATAC reads:
│   │   └── bwa mem -M (paired-end mode)
│   │       → SAM with CB tags preserved
│   └── Output: region_linear.bam (sorted, indexed)
│
└── Pangenome Alignment Pipeline
    ├── RNA reads:
    │   └── GraphAligner with single-end mode
    │       → GAM format with barcode metadata
    ├── ATAC reads:
    │   └── GraphAligner with paired-end mode
    │       → GAM format with barcode metadata
    └── Output: region_pangenome.gam

Key: Same input reads → two different alignments
```

Copy Code

#### Step 3: Alignment Post-Processing

```
Linear pipeline:
- Sort by position
- Add read groups
- Mark duplicates (using UMI for RNA)
- Index BAM files

Pangenome pipeline:
- Convert GAM to BAM (projected onto linear coords where possible)
- Sort and index
- Preserve graph-specific mappings in tags
```

Copy Code

#### Step 4: Critical Quality Metrics

```
For fair comparison, track:
1. Total input reads
2. Reads mapped (linear vs graph)
3. Reads mapped uniquely vs multi-mapped
4. Reads mapped to alt haplotypes (pangenome only)
5. Unmapped reads rescued by pangenome

Create comparison matrix:
                Linear  Pangenome  Improvement
Total reads      X        X          -
Mapped          Y1       Y2         +Z%
High quality    Q1       Q2         +W%
Alt haplotypes   0       A          +A
```

Copy Code

### Biological Significance Tracking:

#### Preserve Throughout Pipeline:

1. **Cell Barcode (CB)** - Links every read to its cell of origin
2. **UMI (UB)** - For RNA deduplication
3. **Read Position** - For variant calling
4. **Mapping Quality** - For filtering

#### Why This Matters:

```
Later phases need:
- CB tag → assign variants to cells → link to cell types
- Position → call variants → compare variant sets
- ATAC accessibility → functional interpretation
```

Copy Code

### Modified Workflow Diagram:

```
10x Multiome FASTQs
        |
   [K-mer Filter]  ← Uses kmers from BOTH references
        |
  Candidate Reads
   /          \
Linear Ref   Pangenome
  (BWA)     (GraphAligner)
   |            |
  BAM          GAM
   |            |
[Preserve CB]  [Preserve CB]
   |            |
Variants_L   Variants_P
   \          /
    Compare:
    - Total variants
    - Unique to pangenome
    - Cell type distribution
```

Copy Code

### Key Differences from Original Approach:

1. **No pre-alignment filtering** - We use k-mer matching instead
2. **Preserve all metadata** - Cell barcodes stay with reads throughout
3. **Fair comparison** - Same input reads to both pipelines
4. **Capture alt haplotypes** - Include reads that might only map to pangenome

### Expected Outcomes:

- Pangenome should capture 5-20% more variants in HLA region (highly polymorphic)
- Fewer additional variants in IGH (less polymorphic)
- Additional variants should cluster in specific cell types
- ATAC data will show if variants are in regulatory regions

This approach ensures we fairly demonstrate the pangenome's advantages while maintaining biological interpretability through single-cell resolution.

Copy