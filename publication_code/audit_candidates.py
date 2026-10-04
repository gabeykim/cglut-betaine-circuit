"""Independent descriptive audit of supplied P1 files. Python standard library only.
Run: python audit_candidates.py path/to/P1_directory path/to/output_directory
This does not reproduce the original candidate generator or biological predictions.
"""
import csv, json, hashlib, sys, platform, re
from pathlib import Path
from collections import Counter, defaultdict

def read(p):
    with p.open(newline='') as f: return list(csv.DictReader(f))

def main(src, out):
    src, out = Path(src), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows=read(src/'proposal_97_NOT_FOR_ORDER.csv')
    base=read(src/'BASELINE_order_manifest_final_v4.csv')
    bm={r['construct_id']:r for r in base}
    parts=defaultdict(list)
    for r in read(src/'proposal_components.csv'): parts[r['proposal_id']].append(r)
    fa={}; current=None
    for line in (src/'proposal_97_NOT_FOR_ORDER.fasta').read_text().splitlines():
        if line.startswith('>'):
            current=line[1:].split(' | ')[0]; assert current not in fa
            fa[current]=''
        elif line.strip(): fa[current]+=line.strip()
    assert len(rows)==97 and len(base)==95
    assert len({r['proposal_id'] for r in rows})==97
    assert len({r['full_sequence'] for r in rows})==97
    assert set(fa)=={r['proposal_id'] for r in rows}
    metrics=[]; change_counts=Counter(); stage_counts=Counter()
    for r in rows:
        pid=r['proposal_id']; seq=r['full_sequence']; parent=bm[r['baseline_id']]
        assert set(seq)<=set('ACGT') and fa[pid]==seq
        assert len(seq)==int(r['length_bp'])
        gc=100*(seq.count('G')+seq.count('C'))/len(seq)
        assert abs(gc-float(r['GC_percent']))<=0.0051
        pos=0
        for c in parts[pid]:
            a,b=int(c['start']),int(c['end'])
            assert a==pos and b-a==len(c['sequence']) and seq[a:b]==c['sequence']
            pos=b
        assert pos==len(seq)
        change=r['change']; change_counts[change]+=1
        orig=parent['full_sequence']
        if change=='unchanged': assert seq==orig
        elif change=='conditional_direct_ATG_TC': assert seq==orig+'TC'
        elif change=='conditional_SD1':
            c=next(c for c in parts[pid] if c['component']=='proposed:SD1')
            a,b=int(c['start']),int(c['end'])
            assert c['sequence']=='AAAGGAGGACAAC' and seq[:a]+seq[b:]==orig
        elif change=='new_matched_scrambled_control':
            partner=next(x for x in rows if x['proposal_id']==r['baseline_id']+'__P1')
            c=next(c for c in parts[pid] if c['component']=='operator:matched_scrambled_opA')
            a,b=int(c['start']),int(c['end']); ps=partner['full_sequence']
            assert seq[:a]==ps[:a] and seq[b:]==ps[b:] and len(seq)==len(ps)
            assert Counter(seq[a:b])==Counter(ps[a:b]) and b-a==16
        else: raise AssertionError(change)
        panel=parent['panel']
        stage='B2 EnvZ' if panel=='kinase_junction_alt_backbone' else ('B1 PhoQ' if panel.startswith('kinase_') else 'A output')
        stage_counts[stage]+=1
        metrics.append({'proposal_id':pid,'baseline_id':r['baseline_id'],'stage':stage,'baseline_panel':panel,'change':change,'length_bp':len(seq),'gc_percent':round(gc,4),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'identity_components_and_edit_checks':'PASS'})
    assert stage_counts=={'A output':73,'B1 PhoQ':18,'B2 EnvZ':6}
    legacy=json.loads((src/'verification.json').read_text())
    for rel,h in legacy['source_sha256'].items(): assert hashlib.sha256((src/('BASELINE_'+Path(rel).name)).read_bytes()).hexdigest()==h
    summary={'audit_scope':'Independent CSV/FASTA/component/edit consistency audit; not biological validation or original pipeline reproduction','python_version':platform.python_version(),'baseline_records':len(base),'candidate_records':len(rows),'unique_ids':97,'unique_sequences':97,'stage_counts':dict(stage_counts),'change_counts':dict(change_counts),'length_min_bp':min(x['length_bp'] for x in metrics),'length_max_bp':max(x['length_bp'] for x in metrics),'gc_min_percent':min(x['gc_percent'] for x in metrics),'gc_max_percent':max(x['gc_percent'] for x in metrics),'all_records_passed':True,'baseline_checksums_match_archived_verification':True,'original_scanner_rerun':False,'folding_rerun':False,'experimental_data_present':False}
    with (out/'candidate_audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
    (out/'audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'input_sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(src.iterdir()) if p.is_file()},indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main(*sys.argv[1:3])
