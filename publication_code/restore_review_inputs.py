"""Restore excluded inputs from a lawfully held complete review directory.
Usage: python publication_code/restore_review_inputs.py /path/to/Glutamicum_Preprint_Draft
No network access; each file must match its archived SHA-256 before copying.
"""
from pathlib import Path
import csv, hashlib, shutil, sys
root=Path(__file__).resolve().parents[1]
source=Path(sys.argv[1]).resolve()
rows=list(csv.DictReader((root/'excluded_inputs.csv').open()))
for row in rows:
    p=source/row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
        raise SystemExit('Missing or nonmatching archived input: '+row['path'])
for row in rows:
    dest=root/row['path'];dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source/row['path'],dest)
print('Restored',len(rows),'hash-verified inputs; these remain excluded from version control.')
