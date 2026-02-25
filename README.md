# Genome Assembly Post-processing & Annotation Toolkit

This repository contains a set of scripts for:

1.  Assembly purging and reference-guided scaffolding\
2.  Assembly quality control (BUSCO, Merqury, read mapping)\
3.  Functional annotation using UniProt + DIAMOND\
4.  Downloading UniProt fungal reference databases

The workflow is designed primarily for fungal genome assemblies but can
be adapted for other eukaryotes.

------------------------------------------------------------------------

# Workflow Overview (Recommended Order)

The scripts are intended to be used in the following order:

1)  run_purge_ragtag.sh\
2)  qc_busco_merqury_flagstat.py\
3)  download_uniprot-fungi.sh\
4)  annotate_uniprotlike.py

------------------------------------------------------------------------

# 1️⃣ run_purge_ragtag.sh

## Purpose

Performs:

-   Self-alignment
-   Duplication purging (purge_dups)
-   Haplotype separation
-   Reference-guided scaffolding (RagTag)
-   Scaffold renaming
-   Chromeister scoring (before and after scaffolding)

## Usage

``` bash
bash run_purge_ragtag.sh \
    assembly.fasta \
    reference.fasta \
    output_directory \
    SAMPLE_PREFIX \
    32
```

### Arguments

  Argument    Description
  ----------- -----------------------------------------
  assembly    Input assembly FASTA (can be gzipped)
  reference   Reference genome FASTA (can be gzipped)
  outdir      Output directory
  prefix      Sample name / prefix
  threads     Number of threads (default: 32)

### Main Outputs

    final_<SAMPLE>/
     ├── <SAMPLE>.purged.final.fasta
     ├── <SAMPLE>.purged.name_map.tsv
     ├── <SAMPLE>.hap.final.fasta (if haplotigs exist)
     └── <SAMPLE>.hap.name_map.tsv

------------------------------------------------------------------------

# 2️⃣ qc_busco_merqury_flagstat.py

## Purpose

Runs assembly quality control:

-   BUSCO (auto-lineage eukaryota)
-   Merqury (k-mer completeness and QV)
-   Illumina read mapping (bwa mem)
-   samtools flagstat

Includes automatic skipping if results already exist.

## Usage

``` bash
python qc_busco_merqury_flagstat.py \
    --fasta final_SAMPLE/SAMPLE.purged.final.fasta \
    --fq1 reads_R1.fastq.gz \
    --fq2 reads_R2.fastq.gz \
    --prefix SAMPLE \
    --cores 32
```

Optional parameters:

    --outdir   Custom output directory (default: QC_<prefix>)
    --k        k-mer size for meryl (default: 21)
    --force    Re-run even if outputs exist

------------------------------------------------------------------------

# 3️⃣ download_uniprot-fungi.sh

## Purpose

Downloads UniProtKB FASTA sequences via REST API for fungal taxa.

Default taxonomy: taxonomy_id:4751 (Fungi)

You may modify the query inside the script if needed.

## Usage

``` bash
bash download_uniprot-fungi.sh
```

Outputs are written to:

    refdb/
     ├── uniprotkb_saccharomycetales.fasta.gz
     └── uniprotkb_fungi.fasta

------------------------------------------------------------------------

# 4️⃣ annotate_uniprotlike.py

## Purpose

Annotates predicted protein FASTA files using DIAMOND hits against
UniProt.

Produces:

-   UniProt-like FASTA headers
-   Mapping table of query → best hit

## Expected DIAMOND Format

The DIAMOND output must contain:

qseqid pident length qlen qcovhsp evalue bitscore stitle

Example DIAMOND command:

``` bash
diamond blastp \
    -q predicted_proteins.faa \
    -d uniprot_db.dmnd \
    -o hits.tsv \
    --outfmt 6 qseqid pident length qlen qcovhsp evalue bitscore stitle \
    --max-target-seqs 1 \
    --evalue 1e-5
```

## Usage

``` bash
python annotate_uniprotlike.py \
    --fasta predicted_proteins.faa \
    --hits hits.tsv \
    --out_fasta annotated.faa \
    --out_map annotation_map.tsv \
    --organism "Saccharomyces cerevisiae" \
    --taxid 4932
```

------------------------------------------------------------------------

# Dependencies

-   minimap2
-   purge_dups
-   RagTag
-   Chromeister
-   BUSCO
-   meryl
-   Merqury
-   bwa
-   samtools
-   DIAMOND

------------------------------------------------------------------------

# Suggested Repository Structure

    .
    ├── run_purge_ragtag.sh
    ├── qc_busco_merqury_flagstat.py
    ├── annotate_uniprotlike.py
    ├── download_uniprot-fungi.sh
    ├── README.md
    └── LICENSE

------------------------------------------------------------------------

# Maintainer

Verstrepen Lab --- KU Leuven\
https://verstrepenlab.sites.vib.be/en\
Contact: michael.abrouk@kuleuven.be
