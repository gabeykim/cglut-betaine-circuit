from pathlib import Path
import csv, json, shutil, zipfile, hashlib, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; WORK=ROOT
OUT=ROOT
REL=ROOT
FIG=REL/'figures'; FIG.mkdir(exist_ok=True)
read=lambda p:list(csv.DictReader(p.open()))
rec=read(ROOT/'project/results/identity_rec.csv')
natives=list(dict.fromkeys(r['native_label'] for r in rec))
donors=['Sal_PhoP','Ec_OmpR']
values=np.array([[float(next(r['pct_identity_aligned_columns'] for r in rec if r['native_label']==n and r['donor_label']==d)) for d in donors] for n in natives])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(8.6,4.2),gridspec_kw={'width_ratios':[1.6,1]})
im=ax[0].imshow(values,cmap='Blues',vmin=20,vmax=55,aspect='auto')
ax[0].set_xticks([0,1],['Salmonella PhoP','E. coli OmpR'])
ax[0].set_yticks(range(len(natives)),[n.replace('Cg_','') for n in natives])
for i in range(len(natives)):
    for j in range(2):ax[0].text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if values[i,j]>44 else 'black')
ax[0].set_title('A  REC aligned-column identity (%)',loc='left',fontweight='bold')
ax[1].bar(['Original\nflagged set','Retained on\nrescoring'],[304,82],color=['#465C6B','#A0B5C1'],width=.6)
for i,n in enumerate([304,82]):ax[1].text(i,n+7,str(n),ha='center')
ax[1].set_ylim(0,345);ax[1].set_ylabel('Sites in the original subset');ax[1].set_title('B  Background sensitivity',loc='left',fontweight='bold')
ax[1].spines[['right','top']].set_visible(False)
fig.tight_layout(pad=1.5)
for ext in ['png','pdf','svg']:fig.savefig(FIG/f'figure_1.{ext}',dpi=300)
plt.close(fig)
with (FIG/'figure_1_identity_values.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['native_regulator']+donors);w.writerows([[n]+list(v) for n,v in zip(natives,values)])
(FIG/'figure_1_site_counts.csv').write_text('comparison,count\noriginal_flagged_subset,304\nretained_after_intergenic_rescoring,82\n')
fig,ax=plt.subplots(1,2,figsize=(8.6,3.2),gridspec_kw={'width_ratios':[1,1.5]})
ax[0].bar(['Output','PhoQ','EnvZ'],[73,18,6],color=['#465C6B','#9CAEA9','#C2B6A4'])
for i,n in enumerate([73,18,6]):ax[0].text(i,n+2,str(n),ha='center')
ax[0].set_ylim(0,87);ax[0].set_ylabel('Fragments');ax[0].set_title('A  Retained collection',loc='left',fontweight='bold')
labels=['Unchanged','Conditional TC extension','Conditional SD1 insertion','New matched controls']
ax[1].barh(labels[::-1],[54,39,2,2][::-1],color='#829AA7')
for i,n in enumerate([54,39,2,2][::-1]):ax[1].text(n+.7,i,str(n),va='center')
ax[1].set_xlim(0,61);ax[1].set_xlabel('Fragments');ax[1].set_title('B  Relationship to the 95-fragment baseline',loc='left',fontweight='bold')
for a in ax:a.spines[['right','top']].set_visible(False)
fig.tight_layout(pad=1.5)
for ext in ['png','pdf','svg']:fig.savefig(FIG/f'figure_2.{ext}',dpi=300)
plt.close(fig)

