"""Package integrity and fresh-output utilities, Python standard library only."""
from pathlib import Path
import hashlib,json
MODULE=Path(__file__).resolve().parents[1]
PACKAGE=MODULE.parent

def check_manifest():
    if not __debug__:raise RuntimeError("Do not disable certificate assertions")
    p=PACKAGE/"MANIFEST.json"
    manifest=json.loads(p.read_text(encoding="utf-8-sig"))
    for name,entry in manifest.items():
        q=(PACKAGE/name).resolve()
        assert q.is_relative_to(PACKAGE.resolve()),name
        raw=q.read_bytes()
        assert len(raw)==entry["bytes"] and hashlib.sha256(raw).hexdigest()==entry["sha256"],name
    return {"files_checked":len(manifest),"manifest_sha256":hashlib.sha256(p.read_bytes()).hexdigest()}

def new_output(name):
    dest=(MODULE/name).resolve()
    if not dest.is_relative_to(MODULE.resolve()) or dest==MODULE.resolve():
        raise ValueError("Output must be a new subdirectory of monthus/")
    if dest.exists():raise FileExistsError("Refusing to overwrite existing output: "+str(dest))
    dest.mkdir(parents=True)
    return dest
