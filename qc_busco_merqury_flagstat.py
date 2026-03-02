#!/usr/bin/env python3
"""
qc_busco_merqury_flagstat.py

Inputs:
  --fasta   path to final assembly (purged final)
  --fq1     Illumina R1 fastq(.gz)
  --fq2     Illumina R2 fastq(.gz)
  --prefix  output prefix
  --outdir  output directory (default: ./QC_<prefix>)
  --cores   cores/threads for BUSCO + mapping + samtools (default: 32)
  --k       meryl k (default: 21)
  --force   re-run even if outputs exist

Workflow:
  1) BUSCO:
     busco -i FASTA -c CORES -m geno -f --auto-lineage-euk --out_path OUTDIR/busco -o <prefix>-BUSCO
     Skips if short_summary*.txt exists.

  2) Merqury:
     meryl k=K count fq1 fq2 output OUTDIR/merqury/<prefix>.meryl
     merqury.sh <prefix>.meryl FASTA <prefix>.merqury_out   (NOTE: outprefix is RELATIVE basename)
     Runs with cwd=OUTDIR/merqury and ensures OUTDIR/merqury/logs exists.
     Skips if <outprefix>.qv and <outprefix>.spectra-cn.hist exist in merqury dir.

  3) Read mapping + QC:
     - Build bwa index for FASTA if missing: bwa index FASTA
     - Map paired-end reads: bwa mem -t CORES FASTA fq1 fq2
     - Sort: samtools sort -@ CORES -o OUTDIR/mapping/<prefix>.sorted.bam
     - Flagstat: samtools flagstat -@ CORES OUTDIR/mapping/<prefix>.sorted.bam > OUTDIR/mapping/<prefix>.flagstat.txt
     Skips if sorted BAM and flagstat file exist (and BAM is non-empty).

Important:
  Some merqury.sh installs prepend "logs/" to the provided outprefix. Passing an absolute
  outprefix breaks logging (creates logs//abs/path/...). This script always passes a basename.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from datetime import datetime


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def die(msg: str, code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def check_file(p: Path, label: str) -> None:
    if not p.exists():
        die(f"{label} not found: {p}")
    if p.is_dir():
        die(f"{label} is a directory, expected a file: {p}")


def which_or_die(prog: str) -> str:
    path = shutil.which(prog)
    if not path:
        die(f"Required program not found in PATH: {prog}")
    return path


def run_cmd(cmd: list[str], log_path: Path, cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{ts()}] Running:\n  " + " ".join(cmd))
    print(f"[{ts()}] Log: {log_path}")
    with open(log_path, "w") as log:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if proc.returncode != 0:
        die(f"Command failed (exit={proc.returncode}). See log: {log_path}")


def run_pipeline_bwa_mem_sort(
    bwa: str,
    samtools: str,
    fasta: Path,
    fq1: Path,
    fq2: Path,
    bam_out: Path,
    threads: int,
    log_path: Path,
    cwd: Path | None = None,
) -> None:
    """
    Run:
      bwa mem -t threads fasta fq1 fq2 | samtools sort -@ threads -o bam_out
    Captures stderr from both commands into log_path.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    bam_out.parent.mkdir(parents=True, exist_ok=True)

    print(f"[{ts()}] Running pipeline:\n  bwa mem | samtools sort")
    print(f"[{ts()}] BAM: {bam_out}")
    print(f"[{ts()}] Log: {log_path}")

    with open(log_path, "w") as log:
        # bwa mem writes alignments to stdout; messages to stderr
        bwa_cmd = [
            bwa, "mem",
            "-t", str(threads),
            str(fasta),
            str(fq1),
            str(fq2),
        ]
        sort_cmd = [
            samtools, "sort",
            "-@", str(threads),
            "-o", str(bam_out),
            "-",  # read SAM/BAM from stdin
        ]

        bwa_p = subprocess.Popen(
            bwa_cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=log,
            text=False,  # binary stream
        )
        assert bwa_p.stdout is not None

        sort_p = subprocess.Popen(
            sort_cmd,
            cwd=str(cwd) if cwd else None,
            stdin=bwa_p.stdout,
            stdout=subprocess.DEVNULL,
            stderr=log,
            text=False,
        )

        # allow bwa to receive SIGPIPE if samtools exits early
        bwa_p.stdout.close()

        sort_rc = sort_p.wait()
        bwa_rc = bwa_p.wait()

    if bwa_rc != 0 or sort_rc != 0:
        die(f"Pipeline failed (bwa exit={bwa_rc}, samtools sort exit={sort_rc}). See log: {log_path}")


def run_flagstat(
    samtools: str,
    bam: Path,
    out_txt: Path,
    threads: int,
    log_path: Path,
    cwd: Path | None = None,
) -> None:
    """
    Run:
      samtools flagstat -@ threads bam > out_txt
    Writes stdout to out_txt, stderr to log_path.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    out_txt.parent.mkdir(parents=True, exist_ok=True)

    cmd = [samtools, "flagstat", "-@", str(threads), str(bam)]
    print(f"[{ts()}] Running:\n  " + " ".join(cmd))
    print(f"[{ts()}] Output: {out_txt}")
    print(f"[{ts()}] Log: {log_path}")

    with open(log_path, "w") as log, open(out_txt, "w") as out:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=out,
            stderr=log,
            text=True,
        )
    if proc.returncode != 0:
        die(f"samtools flagstat failed (exit={proc.returncode}). See log: {log_path}")


# -------------------------
# Success checks
# -------------------------
def busco_success(busco_out_path: Path, run_name: str) -> bool:
    run_dir = busco_out_path / run_name
    if not run_dir.exists():
        return False
    return any(run_dir.rglob("short_summary*.txt"))


def meryl_success(meryl_db: Path) -> bool:
    if not meryl_db.exists() or not meryl_db.is_dir():
        return False
    try:
        next(meryl_db.iterdir())
        return True
    except StopIteration:
        return False


def merqury_success(merqury_dir: Path, out_prefix: str) -> bool:
    qv = merqury_dir / f"{out_prefix}.qv"
    spectra = merqury_dir / f"{out_prefix}.spectra-cn.hist"
    return qv.exists() and spectra.exists()


def bwa_index_success(fasta: Path) -> bool:
    """
    BWA index creates files alongside the fasta:
      .amb .ann .bwt .pac .sa
    """
    exts = [".ամբ", ".ann", ".bwt", ".pac", ".sa"]  # guard; .amb is ASCII, but keep simple below
    # Use correct expected extensions:
    exts = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    return all((Path(str(fasta) + ext)).exists() for ext in exts)


def bam_nonempty(p: Path) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size > 0


def mapping_success(bam: Path, flagstat_txt: Path) -> bool:
    return bam_nonempty(bam) and flagstat_txt.exists() and flagstat_txt.is_file() and flagstat_txt.stat().st_size > 0


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run BUSCO + Merqury + bwa mem mapping (sorted BAM + flagstat), with skip-if-done."
    )
    ap.add_argument("--fasta", required=True, help="Final (purged) assembly FASTA")
    ap.add_argument("--fq1", required=True, help="Illumina FASTQ R1 (can be .gz)")
    ap.add_argument("--fq2", required=True, help="Illumina FASTQ R2 (can be .gz)")
    ap.add_argument("--prefix", required=True, help="Output prefix (used in file/folder names)")
    ap.add_argument("--outdir", default=None, help="Output directory (default: ./QC_<prefix>)")
    ap.add_argument("--cores", type=int, default=32, help="Cores/threads (default: 32)")
    ap.add_argument("--k", type=int, default=21, help="k-mer size for meryl (default: 21)")
    ap.add_argument("--force", action="store_true", help="Force re-run even if outputs exist")
    args = ap.parse_args()

    fasta = Path(args.fasta).resolve()
    fq1 = Path(args.fq1).resolve()
    fq2 = Path(args.fq2).resolve()

    check_file(fasta, "FASTA")
    check_file(fq1, "FASTQ R1")
    check_file(fq2, "FASTQ R2")

    prefix = args.prefix.strip()
    if not prefix:
        die("Prefix cannot be empty")

    outdir = Path(args.outdir).resolve() if args.outdir else Path(f"./QC_{prefix}").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    logs_dir = outdir / "logs"
    busco_out_path = outdir / "busco"
    merqury_dir = outdir / "merqury"
    mapping_dir = outdir / "mapping"
    merqury_dir.mkdir(parents=True, exist_ok=True)
    mapping_dir.mkdir(parents=True, exist_ok=True)

    # Tools
    busco = which_or_die("busco")
    meryl = which_or_die("meryl")
    merqury_sh = which_or_die("merqury.sh")
    bwa = which_or_die("bwa")
    samtools = which_or_die("samtools")

    print("============================================================")
    print(f"[{ts()}] QC workflow: BUSCO + Merqury + bwa mem (sorted BAM + flagstat)")
    print(f"Prefix:   {prefix}")
    print(f"FASTA:    {fasta}")
    print(f"FASTQ1:   {fq1}")
    print(f"FASTQ2:   {fq2}")
    print(f"Outdir:   {outdir}")
    print(f"Cores:    {args.cores}")
    print(f"k-mer:    {args.k}")
    print(f"Force:    {args.force}")
    print("============================================================")

    # -------------------------
    # 1) BUSCO
    # -------------------------
    busco_run_name = f"{prefix}-BUSCO"
    busco_log = logs_dir / f"{prefix}.busco.log"
    busco_out_path.mkdir(parents=True, exist_ok=True)

    if not args.force and busco_success(busco_out_path, busco_run_name):
        print(f"[{ts()}] BUSCO appears complete -> skipping (found short_summary in {busco_out_path / busco_run_name})")
    else:
        busco_cmd = [
            busco,
            "-i", str(fasta),
            "-c", str(args.cores),
            "-m", "geno",
            "-f",
            "--auto-lineage-euk",
            "--out_path", str(busco_out_path),
            "-o", busco_run_name,
        ]
        run_cmd(busco_cmd, busco_log, cwd=outdir)

    # -------------------------
    # 2) Merqury (meryl + merqury.sh)
    # -------------------------
    meryl_db = merqury_dir / f"{prefix}.meryl"
    meryl_log = logs_dir / f"{prefix}.meryl.log"

    if not args.force and meryl_success(meryl_db):
        print(f"[{ts()}] meryl DB appears present -> skipping ({meryl_db})")
    else:
        if meryl_db.exists():
            print(f"[{ts()}] Removing existing meryl DB dir: {meryl_db}")
            shutil.rmtree(meryl_db)
        meryl_cmd = [
            meryl,
            f"k={args.k}",
            "count",
            str(fq1),
            str(fq2),
            "output",
            str(meryl_db),
        ]
        run_cmd(meryl_cmd, meryl_log, cwd=merqury_dir)

    (merqury_dir / "logs").mkdir(parents=True, exist_ok=True)

    out_prefix = f"{prefix}.merqury_out"  # MUST be relative/basename
    merqury_log = logs_dir / f"{prefix}.merqury.log"

    if not args.force and merqury_success(merqury_dir, out_prefix):
        print(f"[{ts()}] Merqury appears complete -> skipping (found {out_prefix}.qv and {out_prefix}.spectra-cn.hist)")
    else:
        merqury_cmd = [
            merqury_sh,
            str(meryl_db),   # absolute ok
            str(fasta),      # absolute ok
            out_prefix       # MUST be relative/basename
        ]
        run_cmd(merqury_cmd, merqury_log, cwd=merqury_dir)

    # -------------------------
    # 3) bwa mem mapping + samtools sort + flagstat
    # -------------------------
    bwa_index_log = logs_dir / f"{prefix}.bwa_index.log"
    bwa_map_log = logs_dir / f"{prefix}.bwa_mem_sort.log"
    flagstat_log = logs_dir / f"{prefix}.flagstat.log"

    bam_out = mapping_dir / f"{prefix}.sorted.bam"
    flagstat_txt = mapping_dir / f"{prefix}.flagstat.txt"

    if not args.force and mapping_success(bam_out, flagstat_txt):
        print(f"[{ts()}] Mapping/flagstat appears complete -> skipping ({bam_out}, {flagstat_txt})")
    else:
        # Index FASTA if needed
        if not args.force and bwa_index_success(fasta):
            print(f"[{ts()}] bwa index appears present -> skipping")
        else:
            bwa_index_cmd = [bwa, "index", str(fasta)]
            run_cmd(bwa_index_cmd, bwa_index_log, cwd=outdir)

        # Run bwa mem | samtools sort
        if args.force and bam_out.exists():
            print(f"[{ts()}] Removing existing BAM: {bam_out}")
            bam_out.unlink()

        run_pipeline_bwa_mem_sort(
            bwa=bwa,
            samtools=samtools,
            fasta=fasta,
            fq1=fq1,
            fq2=fq2,
            bam_out=bam_out,
            threads=args.cores,
            log_path=bwa_map_log,
            cwd=outdir,
        )

        # Flagstat
        if args.force and flagstat_txt.exists():
            print(f"[{ts()}] Removing existing flagstat: {flagstat_txt}")
            flagstat_txt.unlink()

        run_flagstat(
            samtools=samtools,
            bam=bam_out,
            out_txt=flagstat_txt,
            threads=args.cores,
            log_path=flagstat_log,
            cwd=outdir,
        )

    print("============================================================")
    print(f"[{ts()}] DONE")
    print(f"BUSCO output path:   {busco_out_path}  (run name: {busco_run_name})")
    print(f"Merqury directory:   {merqury_dir}")
    print(f"Merqury out prefix:  {out_prefix}  (outputs in merqury dir)")
    print(f"Mapping directory:   {mapping_dir}")
    print(f"Sorted BAM:          {bam_out}")
    print(f"Flagstat:            {flagstat_txt}")
    print(f"Logs:                {logs_dir}")
    print("============================================================")


if __name__ == "__main__":
    main()
