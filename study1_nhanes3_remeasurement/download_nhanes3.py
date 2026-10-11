import hashlib, json, sys, urllib.request
from pathlib import Path
DATA = Path(__file__).resolve().parent / 'studies/nhanes3_testretest/data'

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()
bad = []
for e in json.loads((DATA / 'download_manifest.json').read_text(encoding='utf-8-sig')):
    dest = DATA / e['file']
    if dest.exists() and sha(dest) == e['sha256']:
        print('ok', e['file'])
        continue
    print('downloading', e['file'], f"({e['bytes'] / 1000000.0:,.0f} MB)", flush=True)
    req = urllib.request.Request(e['url'], headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, 'wb') as out:
        for c in iter(lambda: r.read(1 << 20), b''):
            out.write(c)
    ok = sha(dest) == e['sha256']
    print('verified' if ok else 'HASH MISMATCH', e['file'])
    if not ok:
        bad.append(e['file'])
print(f'{17 - len(bad)} of 17 files verified')
sys.exit(1 if bad else 0)
