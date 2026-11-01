"""Resolve shared frozen inputs without duplicating their bytes in the archive."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def candidate_path(name):
    rows=json.loads((ROOT/'cases.json').read_text(encoding='utf-8'))
    matches=[r for r in rows if r['candidate']==name]
    if len(matches)!=1:raise ValueError('Unknown or ambiguous candidate: '+name)
    row=matches[0];p=(ROOT/row['input_path']).resolve()
    if not p.is_relative_to(ROOT.parent.resolve()):raise ValueError('Input path escapes the package')
    if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Candidate hash mismatch: '+name)
    return p
def candidate_paths():
    return [candidate_path(r['candidate']) for r in json.loads((ROOT/'cases.json').read_text(encoding='utf-8')) if r['quadratic']]
