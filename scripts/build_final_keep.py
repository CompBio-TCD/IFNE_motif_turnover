#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def parse_args():
    p=argparse.ArgumentParser(description="Build final species whitelist from QC and biological filters.")
    p.add_argument('--qc-kept', required=True, help='QC *.kept.txt; one species ID per line.')
    p.add_argument('--keep-list', action='append', default=[], help='Whitelist to intersect; may be repeated.')
    p.add_argument('--exclude-list', action='append', default=[], help='Species IDs to remove; may be repeated.')
    p.add_argument('--sanity-tsv', default=None, help='Optional maf sanity TSV; status=FAIL taxa are removed.')
    p.add_argument('--out', required=True, help='Output final *.kept.txt.')
    return p.parse_args()

def load_ids(path):
    path=Path(path).expanduser().resolve()
    if not path.exists(): raise FileNotFoundError(path)
    ids=[]
    with path.open() as fh:
        for raw in fh:
            line=raw.strip()
            if not line or line.startswith('#'): continue
            ids.append(line.split()[0])
    return ids

def main():
    a=parse_args()
    qc_order=load_ids(a.qc_kept)
    current=set(qc_order)
    print(f'QC-passing species: {len(current)}')
    for path in a.keep_list:
        keep=set(load_ids(path)); before=len(current); current &= keep
        print(f'After keep-list {Path(path).name}: {len(current)} (removed {before-len(current)})')
    if a.sanity_tsv:
        df=pd.read_csv(Path(a.sanity_tsv).expanduser().resolve(), sep='\t')
        if not {'species','status'}.issubset(df.columns):
            raise ValueError('sanity TSV must contain species and status columns')
        failed=set(df.loc[df['status'].astype(str).str.upper().eq('FAIL'),'species'].astype(str))
        before=len(current); current -= failed
        print(f'After removing sanity FAIL taxa: {len(current)} (removed {before-len(current)})')
    for path in a.exclude_list:
        exc=set(load_ids(path)); before=len(current); current -= exc
        print(f'After exclusion list {Path(path).name}: {len(current)} (removed {before-len(current)})')
    final=[sp for sp in qc_order if sp in current]
    out=Path(a.out).expanduser().resolve(); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(''.join(f'{sp}\n' for sp in final))
    print(f'Final retained species: {len(final)}')
    print(f'Wrote: {out}')

if __name__=='__main__': main()
