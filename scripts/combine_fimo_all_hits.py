#!/usr/bin/env python3
import argparse,re
from pathlib import Path
import pandas as pd
import yaml

def parse_args():
    p=argparse.ArgumentParser(description="Combine raw FIMO TSVs; use FIMO q-values; preserve all significant positions.")
    p.add_argument('--fimo-dir',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--config',default='../config/config.yaml')
    p.add_argument('--p-thresh',type=float,default=None)
    p.add_argument('--q-thresh',type=float,default=None)
    return p.parse_args()

def parse_filename(path):
    m=re.match(r'(?P<species>[^_]+)__(?P<tf>.+)__(?P<source>jaspar|homococo)\.tsv$',Path(path).name)
    return None if not m else (m.group('species'),m.group('tf'),m.group('source'))

def main():
    a=parse_args(); cfg={}
    cp=Path(a.config).expanduser().resolve()
    if cp.exists():
        with cp.open() as fh: cfg=yaml.safe_load(fh) or {}
    pth=a.p_thresh if a.p_thresh is not None else float(cfg.get('fimo',{}).get('p_value',1e-4))
    qth=a.q_thresh if a.q_thresh is not None else float(cfg.get('fimo',{}).get('q_value',0.05))
    rows=[]
    for path in sorted(Path(a.fimo_dir).expanduser().resolve().glob('*__*__*.tsv')):
        parsed=parse_filename(path)
        if parsed is None: continue
        species,tf,source=parsed
        try: df=pd.read_csv(path,sep='\t')
        except pd.errors.EmptyDataError: continue
        if df.empty: continue
        if 'p-value' not in df.columns or 'q-value' not in df.columns:
            raise ValueError(f'{path}: expected FIMO p-value and q-value columns')
        df['p-value']=pd.to_numeric(df['p-value'],errors='coerce')
        df['q-value']=pd.to_numeric(df['q-value'],errors='coerce')
        df=df.dropna(subset=['p-value','q-value'])
        df=df[(df['p-value']<=pth)&(df['q-value']<=qth)].copy()
        if df.empty: continue
        df.insert(0,'species',species); df.insert(1,'TF',tf); df.insert(2,'source',source)
        rows.append(df)
    if rows:
        out=pd.concat(rows,ignore_index=True)
        keys=[c for c in ['species','TF','source','sequence_name','start','stop','strand','p-value','q-value'] if c in out.columns]
        out=out.drop_duplicates(subset=keys)
    else: out=pd.DataFrame()
    op=Path(a.out).expanduser().resolve(); op.parent.mkdir(parents=True,exist_ok=True); out.to_csv(op,sep='\t',index=False)
    print(f'p threshold: {pth}'); print(f'q threshold: {qth}'); print(f'Significant hits written: {len(out)}'); print(f'Wrote: {op}')

if __name__=='__main__': main()
