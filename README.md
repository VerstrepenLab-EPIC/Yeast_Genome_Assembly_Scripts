# Genome Assembly Post-processing & Annotation Toolkit

This repository provides a streamlined workflow for post-processing, quality control, and functional annotation of genome assemblies.

It is primarily designed for **fungal genome assemblies**, but can be adapted for other eukaryotic genomes.

---

# Dependencies

The following software must be installed and available in your `$PATH`:

- minimap2  
- purge_dups  
- RagTag  
- Chromeister  
- BUSCO  
- meryl  
- Merqury  
- bwa  
- samtools  
- DIAMOND  

---

## Overview

The toolkit includes scripts for:

1. **Assembly refinement and scaffolding**
   - Self-alignment  
   - Duplication purging  
   - Haplotype separation  
   - Reference-guided scaffolding  

2. **Assembly quality control**
   - BUSCO completeness assessment  
   - Merqury k-mer–based evaluation  
   - Illumina read mapping statistics  

3. **Functional annotation**
   - DIAMOND-based annotation against UniProt  
   - UniProt-like FASTA header formatting  

4. **Reference database download**
   - Automated retrieval of UniProt fungal datasets  

---

# 1️⃣ run_purge_ragtag.sh

## Purpose

Performs:

- Assembly self-alignment  
- Duplication purging (`purge_dups`)  
- Haplotype separation  
- Reference-guided scaffolding (RagTag)  
- Scaffold renaming  
- Chromeister scoring (before and after scaffolding)  

## Usage

```bash
bash run_purge_ragtag.sh \
    assembly.fasta \
    reference.fasta \
    output_directory \
    SAMPLE_PREFIX \
    32
```

### Arguments

| Argument   | Description |
|------------|-------------|
| `assembly` | Input assembly FASTA (can be gzipped) |
| `reference` | Reference genome FASTA (can be gzipped) |
| `outdir` | Output directory |
| `prefix` | Sample name / prefix |
| `threads` | Number of threads (default: 32) |

### Main Outputs

```
final_<SAMPLE>/
├── <SAMPLE>.purged.final.fasta
├── <SAMPLE>.purged.name_map.tsv
├── <SAMPLE>.hap.final.fasta        (if haplotigs exist)
└── <SAMPLE>.hap.name_map.tsv
```

---

# 2️⃣ qc_busco_merqury_flagstat.py

## Purpose

Performs assembly quality control:

- BUSCO (auto-lineage: eukaryota)  
- Merqury (k-mer completeness and QV)  
- Illumina read mapping (`bwa mem`)  
- `samtools flagstat`  

Automatically skips completed steps unless `--force` is specified.

## Usage

```bash
python qc_busco_merqury_flagstat.py \
    --fasta final_SAMPLE/SAMPLE.purged.final.fasta \
    --fq1 reads_R1.fastq.gz \
    --fq2 reads_R2.fastq.gz \
    --prefix SAMPLE \
    --cores 32
```

### Optional Parameters

| Parameter | Description |
|-----------|-------------|
| `--outdir` | Custom output directory (default: `QC_<prefix>`) |
| `--k` | k-mer size for meryl (default: 21) |
| `--force` | Re-run steps even if outputs already exist |

---

# 3️⃣ download_uniprot-fungi.sh

## Purpose

Downloads UniProtKB FASTA sequences via the UniProt REST API for fungal taxa.

Default taxonomy:

```
taxonomy_id:4751  (Fungi)
```

The query can be modified directly inside the script if needed.

## Usage

```bash
bash download_uniprot-fungi.sh
```

### Output Directory

```
refdb/
├── uniprotkb_saccharomycetales.fasta.gz
└── uniprotkb_fungi.fasta
```

---

# 4️⃣ annotate_uniprotlike.py

## Purpose

Annotates predicted protein FASTA files using DIAMOND hits against UniProt.

Produces:

- UniProt-like FASTA headers  
- Query → best hit mapping table  

## Required DIAMOND Output Format

The DIAMOND output must contain:

```
qseqid pident length qlen qcovhsp evalue bitscore stitle
```

### Example DIAMOND Command

```bash
diamond blastp \
    -q predicted_proteins.faa \
    -d uniprot_db.dmnd \
    -o hits.tsv \
    --outfmt 6 qseqid pident length qlen qcovhsp evalue bitscore stitle \
    --max-target-seqs 1 \
    --evalue 1e-5
```

## Usage

```bash
python annotate_uniprotlike.py \
    --fasta predicted_proteins.faa \
    --hits hits.tsv \
    --out_fasta annotated.faa \
    --out_map annotation_map.tsv \
    --organism "Saccharomyces cerevisiae" \
    --taxid 4932
```

---

# Maintainer

**Verstrepen Lab — KU Leuven**  
https://verstrepenlab.sites.vib.be/en  

Contact: michael.abrouk@kuleuven.be  
