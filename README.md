# Genome assembly processing and annotation scripts

## Install the environment

Install Conda (Miniforge or Miniconda) first. From this directory, create and activate the shared environment using the included `INSTALL.yaml`:

```bash
conda env create -f INSTALL.yaml
conda activate yeast_assembly
```

If these files are under `Scripts/`, run `cd Scripts` first, or use `conda env create -f Scripts/INSTALL.yaml` from the repository root. The Zenodo draft also contains an identical environment file at its root; create the environment once.

The environment targets **Linux x86_64**, including Linux under WSL2. It supplies Python 3.11, the analysis tools, R and GNU shell utilities. Main tool versions are pinned in the YAML; transitive packages are resolved by Conda, so this is an environment specification rather than a complete lockfile. The recorded Linux dependency solve passed; a complete installation and genome-scale run in this shared environment have not been validated.

Reference genomes, sequencing reads, predicted proteins, UniProt sequences and BUSCO lineage datasets are separate inputs. BUSCO's auto-lineage run needs access to its downloadable datasets. The annotation converter uses Python's standard library; DIAMOND is included for the preceding protein search.

Check the command interfaces after activation:

```bash
bash run_purge_ragtag.sh --help
python qc_busco_merqury_flagstat.py --help
python annotate_uniprotlike.py --help
bash run_tidk_telomeres.sh --help
```

## Included files and workflow

| File | Purpose |
| --- | --- |
| `INSTALL.yaml` | Shared Conda environment specification |
| `run_purge_ragtag.sh` | Duplicate-sequence purging, reference-guided scaffolding and Chromeister plots |
| `qc_busco_merqury_flagstat.py` | Genome BUSCO, Merqury and paired-end Illumina mapping QC |
| `annotate_uniprotlike.py` | Convert predicted protein headers using precomputed DIAMOND hits |
| `run_tidk_telomeres.sh` | Self-contained candidate telomeric-repeat screen for one genome |
| `download_uniprot-fungi.sh` | Optional UniProt fungal FASTA download helper; requires curl separately |

Start with de novo assembly contigs, run purging/scaffolding if appropriate, then assess the resulting scaffolds with QC and telomere screening. Protein prediction is a separate upstream step: provide a predicted protein FASTA and run DIAMOND before the annotation converter. Protein prediction is not included. The optional `download_uniprot-fungi.sh` helper downloads UniProt fungal protein sequences before database construction.

These are generic reusable scripts. They do not reproduce every stage of the project's historical analysis pipelines. Use fresh output directories, retain input checksums and software/database versions, and inspect results and logs before selecting an assembly.

## Assembly purging and RagTag scaffolding

```bash
bash run_purge_ragtag.sh \
    assembly.fasta reference.fasta results SAMPLE 32
```

Arguments are assembly FASTA, reference FASTA, output directory, sample prefix and optional thread count (default 32). Both FASTAs may be gzip-compressed. The output directory must be new or empty. Prefixes must start with a letter or digit and contain only letters, digits, underscores, dots or hyphens.

The script self-aligns split contigs with minimap2, runs `purge_dups`, extracts sequences with `get_seqs`, and scaffolds the retained and available hap sequences with RagTag. It does not calculate long-read depth for purging. The inherited cutoff values are 5/35/90; evaluate them for the input dataset. Override them through `PURGE_LOW`, `PURGE_MID` and `PURGE_HIGH`, with integer values satisfying `0 <= LOW < MID < HIGH`.

Purged scaffolds are sorted by length when the reference contains more than 20 sequences; otherwise they follow reference order, with unmatched sequences retained at the end. Final identifiers become `SAMPLE_scaff1`, `SAMPLE_scaff2`, etc.; name maps record the previous identifiers. Keep purged and hap FASTAs separate because their renamed identifiers can overlap.

Main outputs under `results/`:

```text
SAMPLE_pipeline.log
purged.fa                       # contigs before RagTag
hap.fa                          # extracted sequences or temporary placeholder
purged_SAMPLE_ragtag/            # RagTag FASTA, AGP and other outputs
final_SAMPLE/
    SAMPLE.purged.final.fasta
    SAMPLE.purged.name_map.tsv
    SAMPLE.hap.final.fasta       # only when hap scaffolding was produced
    SAMPLE.hap.name_map.tsv      # only when hap scaffolding was produced
```

Intermediate alignments, BED files and Chromeister matrices/score outputs are also retained. Chromeister uses the reference as query (x-axis), the assembly as database (y-axis), and dimension 2000. The initial comparison uses the input assembly; subsequent comparisons use final purged and available hap scaffolds.

Tools must be on PATH: `minimap2`, `split_fa`, `purge_dups`, `get_seqs`, `samtools`, `ragtag.py`, `CHROMEISTER` and `Rscript`. `compute_score.R` is located on PATH or beside the resolved Chromeister executable. Set `CHROMEISTER_SCORE_SCRIPT` to its full path if discovery fails. `PYTHON_BIN` overrides the Python interpreter.

If the cleaned duplicate BED is empty, the script retains the original assembly. If `get_seqs` fails or produces no usable purged FASTA, it logs a warning and also falls back to the original assembly; the final success message alone does not establish successful purging. A temporary `hap_empty` sequence containing one `N` is used when no usable hap output exists and is excluded from final hap scaffolding. Its presence is not evidence of removed biological sequence.

This script does not perform FCS-GX screening, mitochondrial detection or junction-support validation. RagTag joins are reference-guided hypotheses; Chromeister agreement alone does not establish their correctness.

## Assembly QC

```bash
python qc_busco_merqury_flagstat.py \
    --fasta results/final_SAMPLE/SAMPLE.purged.final.fasta \
    --fq1 reads_R1.fastq.gz --fq2 reads_R2.fastq.gz \
    --prefix SAMPLE --outdir QC_SAMPLE --cores 32 --k 21
```

Use an uncompressed assembly FASTA and paired Illumina FASTQs (plain or gzip-compressed). Defaults are 32 cores, k-mer length 21 and output directory `QC_<prefix>`. `--force` reruns stages and replaces existing results. BWA writes index files beside the assembly FASTA, so provide a writable working copy.

The script runs BUSCO in genome mode with `--auto-lineage-euk`, counts Illumina k-mers with meryl, runs Merqury, and maps reads using `bwa mem` piped to `samtools sort`, followed by `samtools flagstat`.

Outputs include:

- `busco/SAMPLE-BUSCO/`: BUSCO summaries and run outputs.
- `merqury/SAMPLE.meryl`: read k-mer database.
- `merqury/SAMPLE.merqury_out.*`: Merqury QV, completeness and spectra outputs.
- `mapping/SAMPLE.sorted.bam` and `mapping/SAMPLE.flagstat.txt`.
- `logs/`: command logs.

Existing results are skipped using file-presence checks. Those checks do not verify input hashes or complete output integrity; use a new directory/prefix when inputs or parameters change, and inspect partial runs before resuming. The `--cores` setting is passed to BUSCO and mapping tools; meryl's command has no explicit thread limit.

Auto-lineage can select different marker panels and does not reproduce a fixed-lineage comparison automatically. Record the selected lineage/version. `flagstat` reports several quantities: overall mapped reads and properly paired reads are distinct metrics. Merqury results depend on the independent reads, coverage and k-mer settings; low coverage can limit their interpretation.

## Optional UniProt download

The additional `download_uniprot-fungi.sh` helper uses `curl`, gzip, grep and sed. `curl` is not explicitly listed in the shared environment specification. Install it if unavailable:

```bash
conda install -n yeast_assembly -c conda-forge curl
bash download_uniprot-fungi.sh
```

Run this helper from a fresh working directory. It writes `refdb/uniprotkb_saccharomycetales.fasta.gz` and `refdb/uniprotkb_fungi.fasta`. Despite the compressed filename, its query is `taxonomy_id:4751` (Fungi), not a Saccharomycetales-only query. It follows paginated API responses and then decompresses the result.

The current helper truncates an existing compressed output and overwrites the decompressed file. It has no shell fail-fast setting, HTTP-status failure handling, retries or atomic output replacement. A failed/partial download can therefore leave output files behind. Inspect curl messages, validate the gzip archive with `gzip -t`, and verify that the decompressed FASTA is nonempty before database construction. It does not record a database release, response headers or checksums; retain those separately when reproducibility matters. Its syntax was checked, but a live download was not run in this review.

## Protein annotation from DIAMOND hits

Provide a plain protein FASTA with unique nonempty identifiers and a UniProt protein FASTA obtained separately. Keep the database release, query and download date in your provenance. Build and search a local database:

```bash
diamond makedb --in uniprot.fasta --db uniprot_db

diamond blastp \
    --query predicted_proteins.faa --db uniprot_db \
    --out hits.tsv \
    --outfmt 6 qseqid pident length qlen qcovhsp evalue bitscore stitle \
    --max-target-seqs 1 --evalue 1e-5 --threads 32

python annotate_uniprotlike.py \
    --fasta predicted_proteins.faa --hits hits.tsv \
    --out_fasta annotated.faa --out_map annotation_map.tsv \
    --organism "Saccharomyces cerevisiae" --taxid 4932 \
    --keep_uniprot_accession_tag
```

The converter does not run DIAMOND or download databases. Its hits file must be tab-separated, without a header, with exactly the eight columns shown above and UniProt-style subject titles. Rows with fewer than eight fields are skipped; malformed numeric fields raise an error.

Default filters are identity >=30% (`--min_pident 30`), query coverage >=70% (`--min_qcov 0.70`, expressed as a fraction), and E-value <=1e-5 (`--max_evalue 1e-5`). Among retained hits, the highest bit score wins, with lower E-value breaking ties. Searching for only one target can leave a query without an assignment if that hit fails the converter's filters; retain more hits upstream if alternatives are required.

Outputs are a protein FASTA with UniProt-like headers and a query-to-hit TSV. Queries without an accepted, parseable hit are labelled hypothetical proteins. Leading underscores are removed from query identifiers; collisions after normalization receive `__dupN` suffixes. The TSV records raw and normalized identifiers. Sequences remain unchanged. Output directories must already exist, and existing output files are overwritten; use paths distinct from inputs.

Copied `sp`/`tr`, gene-name, `PE` and `SV` fields describe the database hit. They do not establish UniProt deposition, experimental evidence or validated function for the query protein. Treat these assignments as sequence-similarity annotations and retain the mapping table.

## Candidate telomeric-repeat screening

```bash
bash run_tidk_telomeres.sh assembly.fasta -o SAMPLE_telomeres
```

Run one genome per invocation with a new or empty output directory. Plain or gzip-compressed FASTA is accepted. The Python implementation is embedded in the shell script; no sibling Python file is needed. TIDK 0.2.7 and Python are supplied by the environment. `PYTHON_BIN` selects another interpreter; `RAYON_NUM_THREADS` can limit TIDK threads.

Defaults:

| Option | Meaning | Default |
| --- | --- | --- |
| `-l` | Retain sequences strictly longer than this length | 200000 bp |
| `-m`, `-x` | Discovery repeat-length range | 5–30 bp |
| `-t` | TIDK discovery threshold | 20 |
| `-d` | Discovery fraction at each sequence end | 0.05 |
| `-n` | Unique primitive/canonical candidates to validate; 0 means all | 10 |
| `-e` | Terminal validation window | 10000 bp |
| `-c` | Minimum exact tandem copies | 5 |
| `-w` | TIDK search output window | 10000 bp |

Discovery is independent for each genome. A candidate passes when it has at least one detected end and a greater terminal than internal exact-repeat density. Each qualifying array must lie wholly inside a terminal window. Discovery threshold 20 and the five-copy validation threshold serve different purposes. Overlapping phase/strand arrays merge only when their union remains an uninterrupted exact repeat.

Main outputs are `assembly_summary.tsv`, `telomere_counts_by_sequence.tsv`, `best_telomere_candidate.txt`, `sequence_inclusion.tsv`, `candidate_motifs.tsv`, `telomere_validation_summary.tsv`, `telomere_validation_by_sequence.tsv`, `repeat_arrays.tsv`, `run_config.json` and `complete.json`. TIDK discovery/search outputs and logs are retained too. Array coordinates are 1-based inclusive, with lengths, complete copies and distances to the actual sequence ends.

`NO_TELOMERIC_REPEAT_IDENTIFIED` and `NO_SEQUENCES_ABOVE_THRESHOLD` are valid screen outcomes. Counts describe candidate terminal repeat arrays, not confirmed physical telomeres, chromosomes or assembly completeness. A negative screen does not establish biological absence.

## Review and validation

The four scripts pass syntax/CLI checks and small portability smoke tests covering argument validation, spaced paths, annotation hits/no hits, and single-genome telomere detection with the strict length boundary. These checks do not constitute a full purging/QC run or validation of every scientific output. The environment's recorded Linux Conda dependency solve passed; end-to-end installation remains untested.

For GitHub, keep this README, `INSTALL.yaml` and all four core scripts together; include the optional download helper if you intend to provide it. A repository license is not included here; add the license selected by the authors before distributing under explicit reuse terms.

## Maintainer

Verstrepen Lab — KU Leuven  
https://verstrepenlab.sites.vib.be/en  
Contact: michael.abrouk@kuleuven.be
