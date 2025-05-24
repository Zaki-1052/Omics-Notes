## Revised Project Roadmap: Practical Implementation

### Phase 0: Critical Infrastructure Fix (Must Complete First!)

#### 0.1 Get Working Pangenome Graphs

**Option A: Download Pre-built Graphs (Fastest)**

```
Sources to check:
- Human Pangenome Reference Consortium (HPRC)
- UCSC Genome Browser pangenome tracks
- Existing VG graph repositories

What we need:
- HLA region graph (chr6:28.5-33.5Mb)
- IGH region graph (chr14:105.6-107.3Mb)
- Or whole chromosome graphs we can subset
```

Copy Code

**Option B: Build Minimal Graphs (If pre-built unavailable)**

```
Required:
- Download 3-5 alternate HLA haplotypes from IPD-IMGT/HLA
- Get reference sequence for our regions
- Use vg construct to build graph
- Minimum viable graph: reference + 2-3 alternate haplotypes
```

Copy Code

#### 0.2 Verify Graph Integrity

```
Essential checks:
- vg stats -z HLA_region.vg (should show thousands of nodes)
- vg paths -L -x HLA_region.vg (should list multiple haplotypes)
- vg view HLA_region.vg | head (should see actual graph structure)
```

Copy Code

### Phase 1: Simplified 10x Data Processing

#### 1.1 Get Proper Tools

```
Priority installs:
1. UMI-tools (for barcode extraction)
   conda install -c bioconda umi_tools
   
2. VG toolkit (for graph work)
   conda install -c bioconda vg
   
3. GraphAligner (for alignment)
   git clone and compile

4. Get proper barcode whitelist
   Already have arc_barcode_whitelist.txt.gz ✓
```

Copy Code

#### 1.2 Extract and Tag Barcodes (Keep It Simple)

```
For RNA (GEX):
- Input: R1 (barcode+UMI) + R2 (cDNA)
- Extract first 16bp from R1 as cell barcode
- Validate against whitelist
- Add as tag to R2 header
- Output: barcode-tagged R2 FASTQ

For ATAC:
- Input: R1 (genomic) + R2 (barcode) + R3 (genomic)
- Extract first 16bp from R2 as cell barcode
- Add to both R1 and R3 headers
- Keep paired relationship
- Output: barcode-tagged R1/R3 FASTQs
```

Copy Code

### Phase 2: Smart Regional Filtering

#### 2.1 Create Extended Regional References

```
Why: Cast a wider net to catch reads that might map to alternate haplotypes

HLA extended: chr6:27000000-35000000 (8Mb total)
IGH extended: chr14:104000000-109000000 (5Mb total)

Extract from full hg38 reference
```

Copy Code

#### 2.2 Permissive Pre-filtering

```
Strategy: One-pass alignment to extended regions

1. Align barcode-tagged reads to extended regions with BWA
   - Use very permissive settings (-k 15, -w 100)
   - Keep secondary alignments
   
2. Extract:
   - All mapped reads (any quality)
   - Their mate pairs
   - Unmapped reads with mapped mates
   
3. Convert back to FASTQ preserving barcodes
```

Copy Code

### Phase 3: Parallel Alignment Comparison

#### 3.1 Prepare Inputs

```
Same filtered, barcode-tagged FASTQ files for both pipelines:
- GEX: ~500MB-1GB (from 39GB original)
- ATAC: ~200-500MB paired files
```

Copy Code

#### 3.2 Linear Reference Pipeline

```
1. BWA-MEM alignment to original regional references
   - Standard parameters
   - Output: SAM with barcode tags preserved
   
2. Convert to sorted BAM
   - samtools sort
   - samtools index
   
3. Basic variant calling
   - bcftools mpileup + call (faster than deepvariant)
   - Simple SNPs/indels
   - Output: linear_variants.vcf
```

Copy Code

#### 3.3 Pangenome Pipeline

```
1. GraphAligner to pangenome graphs
   - Appropriate parameters for short reads
   - Output: GAM format
   
2. Project to reference coordinates
   - vg surject to get BAM
   - Keep graph-specific mappings in tags
   
3. Graph-aware variant calling
   - vg call on original GAM
   - Output: graph_variants.vcf
```

Copy Code

### Phase 4: Single-Cell Analysis (Simplified)

#### 4.1 Quick Cell Clustering

```
Use existing tools on subset data:
1. Create basic count matrix from aligned reads
2. Simple clustering with scanpy:
   - Basic QC (cells with >100 genes)
   - Normalize and log transform
   - PCA → Leiden clustering
   - UMAP for visualization
   
No need for sophisticated integration initially
```

Copy Code

#### 4.2 Variant-Cell Association

```
Simple approach:
1. For each variant, count reads per cell
2. Create variant-by-cell matrix
3. Overlay on UMAP
4. Calculate enrichment per cluster
```

Copy Code

### Phase 5: Demonstration Outputs

#### 5.1 Core Comparisons

```
Key metrics table:
- Total variants found: Linear vs Graph
- Unique to graph
- Variants per cell type
- Coverage improvements
```

Copy Code

#### 5.2 Visualization

```
Essential plots:
1. Venn diagram of variants
2. UMAP colored by specific variant presence
3. Bar chart of variants per cell type
4. Coverage plot for one example region
```

Copy Code

## Critical Path & Time Estimates

**Day 1:**

- Morning: Fix pangenome graphs (Phase 0) ← BLOCKING ISSUE
- Afternoon: Install tools, process barcodes

**Day 2:**

- Morning: Regional filtering
- Afternoon: Run alignments

**Day 3:**

- Morning: Variant calling
- Afternoon: Basic single-cell analysis

**Day 4:**

- Create visualizations
- Write up results

## Key Simplifications from Previous Plan:

1. **No complex k-mer filtering** - Just use permissive BWA alignment
2. **No TensorQTL/advanced QTL** - Simple association is enough for proof of concept
3. **Use bcftools instead of deepvariant** - Faster, sufficient for demonstration
4. **Basic scanpy clustering** - Don't need perfect cell type annotation
5. **Focus on clear demonstration** - Not publication-ready analysis

## Next Immediate Actions:

1. **Check if pangenome graphs can be downloaded** - Search HPRC/UCSC
2. **If not, identify source haplotypes** - Need FASTA files for HLA alternate alleles
3. **Install vg toolkit** - Required for any graph work
4. **Test graph construction on small example** - Ensure we can build graphs

The most critical realization: **Without working pangenome graphs, nothing else matters.** Everything else is just standard bioinformatics that we know works.

Copy