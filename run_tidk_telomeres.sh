#!/usr/bin/env bash
# Portable single-genome version of the selected source, with validated screening fixes.
# Requires Python >=3.10 and TIDK on PATH; no Conda installation is assumed.
set -euo pipefail
PYTHON_BIN="${PYTHON_BIN:-python3}"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || { echo "ERROR: Python 3 is required" >&2; exit 1; }
TEL_SCREEN_SCRIPT=$(mktemp "${TMPDIR:-/tmp}/tidk_screen.XXXXXXXX.py")
trap 'rm -f "$TEL_SCREEN_SCRIPT"' EXIT
cat > "$TEL_SCREEN_SCRIPT" <<'TIDK_SCREEN_PYTHON'
#!/usr/bin/env python3
"""Exploratory exact-repeat screen; retained sequences are not assumed chromosomes."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read_fasta(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    name, chunks, seen = None, [], set()
    with opener(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    if not chunks:
                        raise ValueError(f'Empty FASTA sequence: {name}')
                    yield name, ''.join(chunks).upper()
                name = line[1:].strip()
                if not name or name.split()[0] in seen:
                    raise ValueError('Empty or duplicate FASTA identifier')
                seen.add(name.split()[0])
                chunks = []
            else:
                if name is None or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', line):
                    raise ValueError('Invalid FASTA sequence or missing header')
                chunks.append(line)
    if name is not None:
        if not chunks:
            raise ValueError(f'Empty FASTA sequence: {name}')
        yield name, ''.join(chunks).upper()


def write_tsv(path, rows, fields):
    with open(path, 'w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def primitive(motif):
    for k in range(1, len(motif) + 1):
        if len(motif) % k == 0 and motif[:k] * (len(motif) // k) == motif:
            return motif[:k]


def variants(motif):
    motif = primitive(motif)
    rc = motif.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return sorted({s[i:] + s[:i] for s in (motif, rc) for i in range(len(s))})


def parse_candidates(path, top_n):
    candidates = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != 2 or header[0] != 'canonical_repeat_unit' or not header[1].startswith('count_repeat_runs_gt_'):
            raise ValueError(f'Unsupported TIDK explore header: {header}')
        for row in reader:
            raw = row[header[0]].upper()
            if not re.fullmatch('[ACGT]+', raw):
                raise ValueError(f'Invalid motif: {raw}')
            motif = min(variants(raw))
            score = int(row[header[1]])
            candidates[motif] = max(score, candidates.get(motif, 0))
    ranked = sorted(candidates.items(), key=lambda x: (-x[1], len(x[0]), x[0]))
    if top_n:
        ranked = ranked[:top_n]
    return [dict(candidate_rank=i, motif=m, repeat_length=len(m), explore_score=s)
            for i, (m, s) in enumerate(ranked, 1)]


def arrays(sequence, motif, copies):
    """Merge overlapping phase/strand representations only when their union is exact.

    Coordinates are 0-based half-open internally; union retains partial terminal
    units, so complete copies are floor(array length / primitive motif length).
    Adjacent nonoverlapping arrays are not joined.
    """
    spans = sorted((m.start(), m.end()) for v in variants(motif)
                   for m in re.finditer('(?:' + v + '){' + str(copies) + ',}', sequence))
    merged = []
    for start, end in spans:
        if merged and start < merged[-1][1]:
            union_end = max(end, merged[-1][1])
            union = sequence[merged[-1][0]:union_end]
            unit = union[:len(motif)]
            if union == (unit * ((len(union) + len(unit) - 1) // len(unit)))[:len(union)]:
                merged[-1][1] = union_end
                continue
        merged.append([start, end])
    return merged


def analyse(sequence, motif, copies, end_size):
    length = len(sequence)
    spans = arrays(sequence, motif, copies)
    # Require the full qualifying array to lie in the terminal window.
    left = [s for s in spans if s[1] <= end_size]
    right = [s for s in spans if s[0] >= length - end_size]
    internal = [s for s in spans if s not in left and s not in right]
    choose = lambda spans, right=False: min(spans, key=lambda s: (length-s[1] if right else s[0], -(s[1]-s[0]))) if spans else None
    def covered_bp(intervals):
        # Distinct exact arrays can overlap by a few bases at an interruption.
        # Count those bases once when estimating regional density.
        end = -1
        total = 0
        for start, stop in sorted(intervals):
            total += max(0, stop - max(start, end))
            end = max(end, stop)
        return total
    return spans, choose(left), choose(right, True), covered_bp(left+right), covered_bp(internal)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('fasta', type=Path)
    p.add_argument('-o', type=Path, default=Path('tidk_telomeres'), dest='outdir')
    p.add_argument('-l', '--min-sequence-length', type=int, default=200000, help='Retain sequences STRICTLY longer than this (default: 200000)')
    for flag, dest, default, helptext in [('-m','minimum',5,'Minimum discovery repeat length'),('-x','maximum',30,'Maximum discovery repeat length'),('-n','top_n',10,'Unique candidates to validate; 0 means all'),('-e','end_size',10000,'Terminal window in bp'),('-w','window',10000,'TIDK search window in bp'),('-c','copies',5,'Minimum exact tandem copies'),('-t','threshold',20,'TIDK discovery threshold (strictly greater than)')]:
        p.add_argument(flag, dest=dest, type=int, default=default, help=helptext + f' (default: {default})')
    p.add_argument('-d', dest='distance', type=float, default=0.05, help='Fraction searched at EACH sequence end (default: 0.05)')
    a = p.parse_args(argv)
    if any(getattr(a,k) <= 0 for k in ['minimum','maximum','end_size','window','copies']) or a.top_n < 0 or a.threshold < 0 or a.min_sequence_length < 0 or a.minimum > a.maximum or not 0 < a.distance <= 0.5:
        p.error('Invalid numeric parameter or range')
    tidk = shutil.which('tidk')
    if not tidk:
        p.error('tidk unavailable: activate the telomere Conda environment')
    seqs = list(read_fasta(a.fasta))
    if not seqs:
        p.error('Empty FASTA')
    retained = [(name,seq) for name,seq in seqs if len(seq) > a.min_sequence_length]
    if any(len(seq) <= 2*a.end_size for _,seq in retained):
        p.error('Terminal windows overlap for retained sequences; reduce -e or increase -l')
    out = a.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        p.error(f'Output directory is not empty: {out}; choose a new directory')
    version = subprocess.check_output([tidk,'--version'], text=True).strip()
    config = {k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    config.update(input_path=str(a.fasta.resolve()), input_sha256=sha256(a.fasta), tidk=tidk, tidk_version=version, python=sys.version, script_sha256=sha256(__file__), started_utc=datetime.now(timezone.utc).isoformat(), rayon_threads=os.environ.get('RAYON_NUM_THREADS','tool default'), commands=[])
    def save_config():
        (out/'run_config.json').write_text(json.dumps(config,indent=2)+'\n')
    save_config()
    inclusion = [dict(sequence=name.split()[0],sequence_length=len(seq),included=int(len(seq)>a.min_sequence_length),reason='LENGTH_GT_THRESHOLD' if len(seq)>a.min_sequence_length else 'LENGTH_LE_THRESHOLD',header=name) for name,seq in seqs]
    write_tsv(out/'sequence_inclusion.tsv',inclusion,['sequence','sequence_length','included','reason','header'])
    filtered = out/'screened_sequences.fasta'
    with filtered.open('w') as f:
        for name,seq in retained:
            f.write('>'+name+'\n')
            for i in range(0,len(seq),80): f.write(seq[i:i+80]+'\n')
    config['screened_fasta_sha256']=sha256(filtered)
    def run(cmd, stdout, stderr):
        config['commands'].append(cmd); save_config()
        with open(stdout,'w') as so, open(stderr,'w') as se:
            subprocess.run(cmd,stdout=so,stderr=se,check=True)
    candidates=[]
    if retained:
        run([tidk,'explore','--minimum',str(a.minimum),'--maximum',str(a.maximum),'--threshold',str(a.threshold),'--distance',str(a.distance),str(filtered)],out/'tidk_explore.tsv',out/'tidk_explore.log')
        candidates=parse_candidates(out/'tidk_explore.tsv',a.top_n)
    write_tsv(out/'candidate_motifs.tsv',candidates,['candidate_rank','motif','repeat_length','explore_score'])
    summaries, details, all_arrays = [], {}, []
    total_bp=sum(len(seq) for _,seq in retained)
    for c in candidates:
        motif=c['motif']; rows=[]; terminal_bp=internal_bp=0
        for name,seq in retained:
            spans,left,right,tbp,ibp=analyse(seq,motif,a.copies,a.end_size)
            terminal_bp+=tbp; internal_bp+=ibp
            row=dict(sequence=name.split()[0],sequence_length=len(seq),motif=motif,left_detected=int(left is not None),right_detected=int(right is not None),detected_ends=int(left is not None)+int(right is not None))
            for side,span in [('left',left),('right',right)]:
                row.update({f'{side}_start_1based':span[0]+1 if span else '',f'{side}_end_1based':span[1] if span else '',f'{side}_array_bp':span[1]-span[0] if span else '',f'{side}_complete_copies':(span[1]-span[0])//len(motif) if span else '',f'{side}_distance_bp':(span[0] if side=='left' else len(seq)-span[1]) if span else ''})
            rows.append(row)
            for start,end in spans:
                region='LEFT' if end<=a.end_size else 'RIGHT' if start>=len(seq)-a.end_size else 'INTERNAL_OR_WINDOW_BOUNDARY'
                all_arrays.append(dict(motif=motif,sequence=name.split()[0],start_1based=start+1,end_1based=end,array_bp=end-start,complete_copies=(end-start)//len(motif),distance_left_bp=start,distance_right_bp=len(seq)-end,region=region))
        terminal_density=terminal_bp/(2*a.end_size*len(retained))
        internal_density=internal_bp/(total_bp-2*a.end_size*len(retained))
        ends=sum(r['detected_ends'] for r in rows)
        # Explicit screening rule, not a biological validation or confidence score.
        supported=ends>0 and terminal_density>internal_density
        summaries.append(dict(**c,supported=int(supported),terminal_ends_detected=ends,both_ends=sum(r['detected_ends']==2 for r in rows),one_end=sum(r['detected_ends']==1 for r in rows),zero_ends=sum(r['detected_ends']==0 for r in rows),terminal_repeat_bp=terminal_bp,internal_repeat_bp=internal_bp,terminal_density=terminal_density,internal_density=internal_density))
        details[motif]=rows
    summaries.sort(key=lambda r:(-r['supported'],-r['terminal_ends_detected'],-r['both_ends'],-(r['terminal_density']-r['internal_density']),-r['terminal_repeat_bp'],r['repeat_length'],r['motif']))
    summary_fields=['candidate_rank','motif','repeat_length','explore_score','supported','terminal_ends_detected','both_ends','one_end','zero_ends','terminal_repeat_bp','internal_repeat_bp','terminal_density','internal_density']
    write_tsv(out/'telomere_validation_summary.tsv',summaries,summary_fields)
    detail_fields=['sequence','sequence_length','motif','left_detected','right_detected','detected_ends']+[f'{side}_{field}' for side in ['left','right'] for field in ['start_1based','end_1based','array_bp','complete_copies','distance_bp']]
    write_tsv(out/'telomere_validation_by_sequence.tsv',[r for c in summaries for r in details[c['motif']]],detail_fields)
    write_tsv(out/'repeat_arrays.tsv',all_arrays,['motif','sequence','start_1based','end_1based','array_bp','complete_copies','distance_left_bp','distance_right_bp','region'])
    best=next((s for s in summaries if s['supported']),None)
    status='CANDIDATE_IDENTIFIED' if best else 'NO_TELOMERIC_REPEAT_IDENTIFIED' if retained else 'NO_SEQUENCES_ABOVE_THRESHOLD'
    if best:
        final=[dict(r,screen_status=status) for r in details[best['motif']]]
    else:
        final=[dict(sequence=name.split()[0],sequence_length=len(seq),motif='',left_detected=0,right_detected=0,detected_ends=0,screen_status=status) for name,seq in retained]
    write_tsv(out/'telomere_counts_by_sequence.tsv',final,detail_fields+['screen_status'])
    assembly=dict(screen_status=status,best_motif=best['motif'] if best else '',input_sequences=len(seqs),retained_sequences=len(retained),excluded_sequences=len(seqs)-len(retained),candidate_motifs_tested=len(candidates),detected_ends=sum(r['detected_ends'] for r in final),possible_sequence_ends=2*len(retained),both_ends=sum(r['detected_ends']==2 for r in final),one_end=sum(r['detected_ends']==1 for r in final),zero_ends=sum(r['detected_ends']==0 for r in final))
    write_tsv(out/'assembly_summary.tsv',[assembly],list(assembly))
    (out/'best_telomere_candidate.txt').write_text('\n'.join(f'{k}: {v}' for k,v in assembly.items())+'\nCounts refer to candidate exact repeat arrays within the terminal window, not confirmed telomeres or chromosomes. No detection does not establish biological absence.\n')
    search=out/'search'; search.mkdir()
    for c in candidates:
        prefix=f"candidate_{c['candidate_rank']:02d}_{c['motif']}"
        run([tidk,'search','--string',c['motif'],'--window',str(a.window),'--output',prefix,'--dir',str(search),'--extension','tsv',str(filtered)],search/(prefix+'.stdout.log'),search/(prefix+'.stderr.log'))
    config['finished_utc']=datetime.now(timezone.utc).isoformat(); save_config()
    assert sha256(a.fasta)==config['input_sha256'], 'Input changed during analysis'
    (out/'complete.json').write_text(json.dumps(assembly,indent=2)+'\n')
    print(json.dumps(assembly),flush=True)


if __name__=='__main__':
    main()
TIDK_SCREEN_PYTHON
"$PYTHON_BIN" "$TEL_SCREEN_SCRIPT" "$@"
