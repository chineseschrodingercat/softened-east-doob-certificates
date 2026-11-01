"""Check every shipped byte before executing a frozen-certificate replay."""
import hashlib,json
from pathlib import Path
def verify_manifest(root):
    if not __debug__:raise RuntimeError('Do not use Python -O or PYTHONOPTIMIZE for verification')
    root=Path(root).resolve();manifest=json.loads((root/'MANIFEST.json').read_text(encoding='utf-8'))
    for name,item in manifest.items():
        p=(root/name).resolve()
        if not p.is_relative_to(root):raise ValueError('Manifest path escapes package')
        data=p.read_bytes()
        if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:raise ValueError('File mismatch: '+name)
    return len(manifest)

if __name__=='__main__':
    print(json.dumps({'status':'PASS','files':verify_manifest(Path(__file__).resolve().parents[1])}))
