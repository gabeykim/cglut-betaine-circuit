"""Build a reversible proposal; never writes to the source project.
Usage: python3 build_revision_proposal.py /path/to/Glutamicum /path/to/output
Requires the original project's scripts and cached literature.
"""
import csv, json, sys, hashlib, shutil, re
from pathlib import Path
from collections import Counter, defaultdict
from dataclasses import asdict

repo, out = map(lambda s: Path(s).resolve(), sys.argv[1:3])
assert repo != out and repo not in out.parents
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(repo / 'scripts'))
import seqchecks as S
from calibration_constructs import henke_primers

def read(p):
    with p.open(newline='') as f: return list(csv.DictReader(f))
def write(name, rows):
    with (out/name).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

sources=[repo/'designs'/n for n in ['order_manifest_final_v4.csv','order_manifest_final_v4_components.csv']]
hashes={str(p.relative_to(repo)):sha(p) for p in sources}
base=read(sources[0]); comp=read(sources[1]); assert len(base)==95
byid={r['construct_id']:r for r in base}; assert len(byid)==95
maps=defaultdict(list)
for c in comp:
    cid=c['construct_id']; a,b=int(c['start']),int(c['end'])
    maps[cid].append((c['component'],byid[cid]['full_sequence'][a:b],c['provenance']))
for cid, parts in maps.items():
    assert ''.join(p[1] for p in parts)==byid[cid]['full_sequence']
sd='AAAGGAGGACAAC'
assert sd==json.loads((repo/'data/cache/translational_sources.json').read_text())['approach_B_source']['sd_motif']
scramble_id=next(x for x in byid if x.startswith('V2-OTH-02_'))
scramble=next(p[1] for p in maps[scramble_id] if p[0]=='operator:scrambled_opA')
primers,_=henke_primers()
reporter_prefix='ATGAGTAAAGGAGAAGAACTTTTCA'
assert primers['HN49'].endswith(reporter_prefix)

proposed=[]; components=[]; edits=[]; details=[]; checks=[]
def emit(cid,parent,parts,scope,change,is_new=False):
    seq=''.join(p[1] for p in parts); orig=byid[parent]['full_sequence']
    assert set(seq)<=set('ACGT')
    proposed.append(dict(proposal_id=cid,baseline_id=parent,scope=scope,status='DRAFT_NOT_FOR_ORDER',change=change,
                         length_bp=len(seq),GC_percent=round(S.gc(seq),2),full_sequence=seq))
    pos=0
    for label,s,prov in parts:
        components.append(dict(proposal_id=cid,component=label,start=pos,end=pos+len(s),sequence=s,provenance=prov))
        pos+=len(s)
    assert pos==len(seq)
    # A coordinate map identifies which new positions still correspond to original bases.
    from difflib import SequenceMatcher
    matcher=SequenceMatcher(None,orig,seq,autojunk=False)
    mapping={}
    for tag,a,b,c,d in matcher.get_opcodes():
        if tag=='equal': mapping.update({j:a+j-c for j in range(c,d)})
        else: edits.append(dict(proposal_id=cid,baseline_id=parent,operation=tag,baseline_start=a,baseline_end=b,
                                proposal_start=c,proposal_end=d,old_sequence=orig[a:b],new_sequence=seq[c:d]))
    sites=S.site_scan(seq); oldsites={(x.enzyme,x.start,x.end,x.strand) for x in S.site_scan(orig)}
    gained=[]
    for x in sites:
        positions=[mapping.get(j) for j in range(x.start,x.end)]
        continuous=all(p is not None for p in positions) and positions==list(range(positions[0],positions[0]+len(positions)))
        if not continuous or (x.enzyme,positions[0],positions[-1]+1,x.strand) not in oldsites: gained.append(asdict(x))
    known=[]; k=0
    for label,s,_ in parts:
        if label.startswith('forecistron:'): known.append(k)
        k+=len(s)
    output=not byid[parent]['panel'].startswith('kinase_')
    boundary=seq+reporter_prefix if output else None
    boundary_sites=[asdict(x) for x in S.site_scan(boundary) if x.start<len(seq)<x.end] if boundary else []
    basic=S.basic_checks(seq)
    gg=[x for x in sites if x.enzyme in S.GOLDEN_GATE]
    assert not gg, (cid,gg)
    coupled=any(p[0]=='TAATG_frag' for p in parts)
    if coupled:
        assert boundary[len(seq)-2:len(seq)+3]=='TAATG'
        assert (len(seq)-2-known[0])%3==0
    if 'conditional_SD1' in change:
        a=known[0]; assert seq[a-13:a]==sd
        assert seq[a:a+3]=='ATG'
        assert a-(a-13+sd.index('GGAGG')+5)==5
        assert seq[a:]==orig[a-13:]
    checks.append(dict(proposal_id=cid,length_bp=len(seq),GC_percent=basic['GC_percent'],
                       golden_gate_hits=len(gg),all_35_enzyme_hits=len(sites),gained_enzyme_hits=len(gained),
                       homopolymer_runs_gt6=len(basic['homopolymers_gt6']),
                       sigA_like_hits_both_strands=len(S.sigA_scan(seq,both_strands=True)),
                       reporter_boundary_test='primer_prefix_only' if output else 'not_applicable',
                       reporter_boundary_enzyme_hits=len(boundary_sites),
                       coupling_frame='preserved' if coupled else 'not_applicable',
                       folding='NOT_RECOMPUTED',full_plasmid='NOT_AVAILABLE'))
    details.append(dict(proposal_id=cid,sites=[asdict(x) for x in sites],gained_sites=gained,
                        sigA=[asdict(x) for x in S.sigA_scan(seq,both_strands=True)],
                        sd_atg_legacy_scan=[asdict(x) for x in S.sd_atg_scan(seq,known_starts=known)],
                        host_rm=[asdict(x) for x in S.host_rm_scan(seq)],basic=basic,
                        reporter_prefix_boundary_sites=boundary_sites))

plain=[]; psoda=[]
for r in base:
    cid=r['construct_id']; parts=list(maps[cid]); labels=[p[0] for p in parts]
    scope='main'; change='unchanged'
    if 'core:Psyn' in labels and not any(x.startswith('forecistron:') for x in labels) and not cid.startswith('CAL-'):
        assert parts[-1][0]=='core:Psyn'
        parts.append(('proposed:direct_ATG_spacer_TC','TC','Designed option: use only with direct reporter ATG; no extra vector scar.'))
        plain.append(cid); change='conditional_direct_ATG_TC'
        core_start=sum(len(p[1]) for p in parts[:-2]); core=parts[-2][1]
        hx_last=core.index('TATAAT')+5
        assert len(''.join(p[1] for p in parts))-(core_start+hx_last)-1==7
    if cid.startswith(('V2-CX-18_','V2-CX-19_')):
        idx=next(i for i,p in enumerate(parts) if p[0].startswith('forecistron:'))
        parts.insert(idx,('proposed:SD1',sd,'Liu 2023 SD2 motif reused here as a designed SD1; not experimentally validated.'))
        psoda.append(cid); scope='conditional_PsodA'; change='conditional_SD1'
    if cid.startswith(('V2-OTH-13_','V2-OTH-14_','V2-OTH-15_')): scope='optional_Ptuf'
    if r['panel']=='kinase_junction_alt_backbone': scope='deferred_EnvZ'
    emit(cid+'__P1',cid,parts,scope,change)
assert len(plain)==39 and len(psoda)==2
for number,prefix in [(1,'V2-SUB-04_'),(2,'V2-CX-06_')]:
    parent=next(x for x in byid if x.startswith(prefix))
    existing=[x for x in components if x['proposal_id']==parent+'__P1']
    parts=[(x['component'],x['sequence'],x['provenance']) for x in existing]
    positions=[i for i,p in enumerate(parts) if p[0]=='operator:opA_phoP']; assert len(positions)==1
    i=positions[0]; assert len(parts[i][1])==len(scramble) and sorted(parts[i][1])==sorted(scramble)
    parts[i]=('operator:matched_scrambled_opA',scramble,'Existing V2-OTH-02 scramble transferred into the exact parent context; binding loss unvalidated.')
    emit(f'NEW-CTRL-{number:02d}__P1',parent,parts,'main','new_matched_scrambled_control',True)
    proposed_parent=next(x['full_sequence'] for x in proposed if x['proposal_id']==parent+'__P1')
    a=sum(len(x[1]) for x in parts[:i]); b=a+len(scramble)
    assert proposed[-1]['full_sequence'][:a]==proposed_parent[:a]
    assert proposed[-1]['full_sequence'][b:]==proposed_parent[b:]
assert len(proposed)==97 and len({r['full_sequence'] for r in proposed})==97
assert Counter(r['scope'] for r in proposed)=={'main':86,'conditional_PsodA':2,'optional_Ptuf':3,'deferred_EnvZ':6}
write('proposal_97_NOT_FOR_ORDER.csv',proposed)
write('proposal_components.csv',components)
write('changes_from_baseline.csv',edits)
write('recomputed_checks.csv',checks)
(out/'scan_details.json').write_text(json.dumps(details,indent=2)+'\n')
(out/'proposal_97_NOT_FOR_ORDER.fasta').write_text(''.join('>'+r['proposal_id']+' | '+r['scope']+' | DRAFT_NOT_FOR_ORDER\n'+r['full_sequence']+'\n' for r in proposed))
for p in sources: shutil.copy2(p,out/('BASELINE_'+p.name))
assert hashes=={str(p.relative_to(repo)):sha(p) for p in sources}
facts=dict(baseline_count=95,proposed_count=97,baseline_records_unchanged=54,
           conditional_TC_edits=39,conditional_SD1_edits=2,new_controls=2,
           scope_counts=dict(Counter(r['scope'] for r in proposed)),
           source_sha256=hashes,reporter_prefix=reporter_prefix,
           gained_site_findings=[{'id':x['proposal_id'],'sites':x['gained_sites']} for x in details if x['gained_sites']],
           boundary_site_findings=[{'id':x['proposal_id'],'sites':x['reporter_prefix_boundary_sites']} for x in details if x['reporter_prefix_boundary_sites']],
           positive_control_fixtures=S.positive_control_check(force=True))
(out/'verification.json').write_text(json.dumps(facts,indent=2)+'\n')
shutil.copy2(Path(__file__),out/'build_revision_proposal.py')
print(json.dumps({k:v for k,v in facts.items() if k not in ['source_sha256','positive_control_fixtures']},indent=2))
