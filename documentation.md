# Single-Cell Multi-omics Pangenome Analysis Project - Complete Documentation

## Project Overview

**Goal**: Demonstrate that aligning single-cell multi-omics data to the human pangenome reference captures more genetic variants than traditional linear reference alignment, and link these variants to specific cell types.

**Datasets**:

- 10x Genomics Multiome (paired RNA + ATAC-seq from same cells)
- Target regions: HLA (chr6:28.5-33.5Mb) and IGH (chr14:105.6-107.3Mb)

## Initial State and Critical Issues Discovered

### Starting Files

```
~/Omics-UBIC/pangenome_project/
├── data/
│   ├── pbmc_granulocyte_sorted_3k/          # Raw FASTQs (Aug 2020)
│   │   ├── gex/  (RNA: L003, L004)
│   │   └── atac/ (ATAC: L001-L004)
│   ├── pbmc_granulocyte_sorted_10k*.h5      # Cell Ranger output (Apr 2021)
│   └── pbmc_granulocyte_sorted_10k*.tsv.gz  # ATAC fragments (Apr 2021)
├── references/
│   ├── hla_hg38.fa  (5.1MB)
│   └── igh_hg38.fa  (1.5MB)
└── graphs/
    ├── HLA_region.vg (0 bytes - EMPTY!)
    └── IGH_region.vg (0 bytes - EMPTY!)
```

Copy Code

### Critical Issue #1: Empty Pangenome Graphs

The pangenome graph files were 0 bytes - completely empty. Without these, the entire project premise was impossible.

**Resolution**:

```
# Downloaded HPRC pangenome graphs
cd ~/Omics-UBIC/pangenome_project/pangenome_source
aws s3 cp --no-sign-request s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.xg .
aws s3 cp --no-sign-request s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gbwt .

# Extracted regional graphs
cd ../graphs
vg chunk -x ../pangenome_source/hprc-v1.0-mc-grch38.xg \
         -g ../pangenome_source/hprc-v1.0-mc-grch38.gbwt \
         -p GRCh38.chr6:28500000-33500000 -c 10000 \
         > HLA_region_full.vg

vg chunk -x ../pangenome_source/hprc-v1.0-mc-grch38.xg \
         -g ../pangenome_source/hprc-v1.0-mc-grch38.gbwt \
         -p GRCh38.chr14:105600000-107043000 -c 10000 \
         > IGH_region_full.vg

# Created indexes
vg index -x HLA_region_full.xg HLA_region_full.vg
vg index -x IGH_region_full.xg IGH_region_full.vg
```

Copy Code

**Result**:

- HLA graph: 287,996 nodes, 397,041 edges (23.6MB)
- IGH graph: 206,451 nodes, 290,182 edges (14MB)

### Critical Issue #2: Dataset Mismatch

Discovered we had:

- 3k raw FASTQ files (no Cell Ranger output)
- 10k Cell Ranger output (no raw FASTQs)
- These were DIFFERENT experiments from different dates

**Initial failed approach**:

- Concatenated all reads: `zcat gex/*_R2_001.fastq.gz > gex_R2_all.fastq` (39GB!)
- Used 10k barcodes with 3k data - resulted in almost no matches

## Correct Data Processing Pipeline

### Step 1: Obtain Matching 3k Cell Ranger Output

```
cd ~/Omics-UBIC/pangenome_project/data

# Downloaded 3k Cell Ranger output
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz
```

Copy Code

### Step 2: Extract Valid Cell Barcodes

```
# extract_3k_barcodes.py
import h5py
with h5py.File('pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5', 'r') as f:
    barcodes = f['matrix/barcodes'][:]
    print(f'Number of 3k cells: {len(barcodes)}')  # Result: 2,711 cells
    with open('cell_barcodes_3k.txt', 'wb') as out:
        for bc in barcodes:
            out.write(bc + b'\n')
```

Copy Code

### Step 3: Create Barcode Extraction Scripts

**RNA Barcode Extraction** (`extract_barcodes_rna.py`):

```
#!/usr/bin/env python3
import gzip
import sys

# Read cell barcodes
valid_barcodes = set()
with open('cell_barcodes_3k.txt', 'rb') as f:
    for line in f:
        bc = line.strip().decode()
        if bc.endswith('-1'):
            bc = bc[:-2]
        valid_barcodes.add(bc)

print(f"Loaded {len(valid_barcodes)} valid barcodes", file=sys.stderr)

# Process RNA reads
r1_file = sys.argv[1]  # Barcode + UMI
r2_file = sys.argv[2]  # cDNA sequence

with gzip.open(r1_file, 'rt') as r1, gzip.open(r2_file, 'rt') as r2:
    while True:
        # Read 4 lines from each file
        r1_header = r1.readline().strip()
        r1_seq = r1.readline().strip()
        r1_plus = r1.readline().strip()
        r1_qual = r1.readline().strip()
        
        r2_header = r2.readline().strip()
        r2_seq = r2.readline().strip()
        r2_plus = r2.readline().strip()
        r2_qual = r2.readline().strip()
        
        if not r1_header:
            break
            
        # Extract barcode (first 16bp of R1) and UMI (next 12bp)
        barcode = r1_seq[:16]
        umi = r1_seq[16:28]
        
        # Check if valid cell
        if barcode in valid_barcodes:
            # Add barcode and UMI to R2 header
            new_header = f"{r2_header} CB:Z:{barcode} UB:Z:{umi}"
            print(new_header)
            print(r2_seq)
            print(r2_plus)
            print(r2_qual)
```

Copy Code

**ATAC Barcode Extraction** (`extract_barcodes_atac.py`):

```
#!/usr/bin/env python3
import gzip
import sys

# Read cell barcodes
valid_barcodes = set()
with open('cell_barcodes_3k.txt', 'rb') as f:
    for line in f:
        bc = line.strip().decode()
        if bc.endswith('-1'):
            bc = bc[:-2]
        valid_barcodes.add(bc)

print(f"Loaded {len(valid_barcodes)} valid barcodes", file=sys.stderr)

# Process ATAC reads
r1_file = sys.argv[1]  # Genomic DNA forward
r2_file = sys.argv[2]  # Barcode
r3_file = sys.argv[3]  # Genomic DNA reverse

with gzip.open(r1_file, 'rt') as r1, gzip.open(r2_file, 'rt') as r2, gzip.open(r3_file, 'rt') as r3:
    while True:
        # Read 4 lines from each file
        r1_header = r1.readline().strip()
        r1_seq = r1.readline().strip()
        r1_plus = r1.readline().strip()
        r1_qual = r1.readline().strip()
        
        r2_header = r2.readline().strip()
        r2_seq = r2.readline().strip()
        r2_plus = r2.readline().strip()
        r2_qual = r2.readline().strip()
        
        r3_header = r3.readline().strip()
        r3_seq = r3.readline().strip()
        r3_plus = r3.readline().strip()
        r3_qual = r3.readline().strip()
        
        if not r1_header:
            break
            
        # Extract barcode (first 16bp of R2)
        barcode = r2_seq[:16]
        
        # Check if valid cell
        if barcode in valid_barcodes:
            # Add barcode to headers
            read_id = r1_header.split()[0]
            # Output R1
            print(f"{read_id}/1 CB:Z:{barcode}")
            print(r1_seq)
            print(r1_plus)
            print(r1_qual)
            # Output R3 as R2
            print(f"{read_id}/2 CB:Z:{barcode}")
            print(r3_seq)
            print(r3_plus)
            print(r3_qual)
```

Copy Code

### Step 4: Process All Reads with Barcodes

**Note**: Encountered disk quota issue. Solution: Remove old concatenated files:

```
rm -f filtered_reads/gex_R2_all.fastq      # 39GB
rm -f filtered_reads/atac_R1_all.fastq     # 12GB
rm -f filtered_reads/atac_R3_all.fastq     # 12GB
rm -f filtered_reads/gex_hla_mapped.bam    # 4.9GB
rm -f filtered_reads/gex_hla_mapped.fastq.sam  # 23GB
```

Copy Code

**Processing script** (`reprocess_3k_barcodes.sh`):

```
#!/bin/bash
set -e

# Process RNA
for lane in L003 L004; do
    python3 extract_barcodes_rna.py \
        pbmc_granulocyte_sorted_3k/gex/pbmc_granulocyte_sorted_3k_S1_${lane}_R1_001.fastq.gz \
        pbmc_granulocyte_sorted_3k/gex/pbmc_granulocyte_sorted_3k_S1_${lane}_R2_001.fastq.gz \
        > rna_barcoded_${lane}.fastq
done
cat rna_barcoded_L*.fastq > rna_barcoded_all.fastq
rm rna_barcoded_L*.fastq

# Process ATAC
for lane in L001 L002 L003 L004; do
    python3 extract_barcodes_atac.py \
        pbmc_granulocyte_sorted_3k/atac/pbmc_granulocyte_sorted_3k_S12_${lane}_R1_001.fastq.gz \
        pbmc_granulocyte_sorted_3k/atac/pbmc_granulocyte_sorted_3k_S12_${lane}_R2_001.fastq.gz \
        pbmc_granulocyte_sorted_3k/atac/pbmc_granulocyte_sorted_3k_S12_${lane}_R3_001.fastq.gz \
        > atac_barcoded_${lane}.fastq
done
cat atac_barcoded_L*.fastq > atac_barcoded_all.fastq
rm atac_barcoded_L*.fastq
```

Copy Code

### Step 5: Regional Filtering with BWA

Create extended references and align to filter reads:

```
# Create extended references (±1-2Mb buffer)
cd ../references
samtools faidx /home/odyssey-comp-07/Omics-UBIC/hg38.fa chr6:27000000-35000000 > hla_extended.fa
samtools faidx /home/odyssey-comp-07/Omics-UBIC/hg38.fa chr14:104000000-109000000 > igh_extended.fa
bwa index hla_extended.fa
bwa index igh_extended.fa

cd ../data
# Align with permissive settings
bwa mem -t 4 -k 15 -w 100 -M ../references/hla_extended.fa rna_barcoded_all.fastq | samtools view -b > rna_hla_extended.bam
bwa mem -t 4 -k 15 -w 100 -M ../references/igh_extended.fa rna_barcoded_all.fastq | samtools view -b > rna_igh_extended.bam

# Extract mapped + unmapped reads
samtools view -b -F 4 rna_hla_extended.bam > rna_hla_mapped.bam
samtools view -b -F 4 rna_igh_extended.bam > rna_igh_mapped.bam
samtools view -b -f 4 rna_hla_extended.bam > rna_unmapped.bam

# Merge and convert back to FASTQ
samtools merge -f rna_regions_all.bam rna_hla_mapped.bam rna_igh_mapped.bam rna_unmapped.bam
samtools fastq rna_regions_all.bam > rna_filtered.fastq
```

Copy Code

## Current Project State

### Successfully Completed:

1. ✅ Downloaded and extracted HPRC pangenome graphs for HLA/IGH regions
2. ✅ Obtained matching 3k Cell Ranger output files
3. ✅ Created barcode extraction scripts for 10x data
4. ✅ Identified and resolved dataset mismatch issue

### Ready for Next Steps:

1. Regional filtering of barcoded reads (BWA alignment to extended regions)
2. Parallel alignment to linear reference vs pangenome graph
3. Variant calling and comparison
4. Single-cell analysis and variant-cell type association

### Key Files Created:

- `HLA_region_full.vg/xg` - Pangenome graph with 287k nodes
- `IGH_region_full.vg/xg` - Pangenome graph with 206k nodes
- `cell_barcodes_3k.txt` - 2,711 valid cell barcodes
- `extract_barcodes_rna.py` - RNA barcode extraction script
- `extract_barcodes_atac.py` - ATAC barcode extraction script

### Tools Installed:

- vg toolkit (for pangenome analysis)
- BWA, samtools (for alignment)
- GraphAligner (compiled but not yet used)

This documentation provides the exact steps to reproduce the analysis up to the current point, including all issues encountered and their resolutions.