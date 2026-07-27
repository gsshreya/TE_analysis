#!/usr/bin/env python3
"""Block-wise, haplotype-level sharing in +/-500 kb around the CABYR Alu.

For the six carriers, compares sharing in each block between:
  (1) their ALU-bearing haplotype  (the founder haplotype)
  (2) their other, non-ALU haplotype (built-in control)
and gives a haplotype-level (12-haplotype) identity view over the window.

Shows: sharing is specific to the ALU haplotype and decays outside the core
(non-core blocks), while the non-ALU haplotypes stay at background sharing.
"""
import gzip, os, itertools, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE=os.path.dirname(os.path.abspath(__file__))
VCF=os.environ.get('CABYR_VCF','/path/to/hpc/cabyr/te_region_cabyr.vcf.gz')
FIG=os.path.join(HERE,'..','figures','new_fig'); os.makedirs(FIG,exist_ok=True)
RES=os.path.join(HERE,'..','results'); os.makedirs(RES,exist_ok=True)
INS=24159972; FLANK=500_000; BLOCK=100_000
LO,HI=INS-FLANK,INS+FLANK
CAR=['DE27046720QK','DJ41072338LK','FD54524847UE','VI92225701RH','WE47272522PX','ZG03084144JQ']

def info_cm(info):
    for kv in info.split(';'):
        if kv.startswith('CM='): return float(kv[3:])
    return np.nan

# ---------- parse VCF (restricted to the +/-500 kb window) ----------
pos=[]; rows=[]; alu_gt=None; samples=None
with gzip.open(VCF,'rt') as fh:
    for line in fh:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples=line.rstrip('\n').split('\t')[9:]
            cix=[samples.index(c) for c in CAR]; continue
        f=line.rstrip('\n').split('\t'); p=int(f[1]); gts=f[9:]
        if p==INS and 'INS:ME:ALU' in f[4]:
            alu_gt=[gts[i].split(':')[0] for i in cix]; continue
        if p<LO or p>HI: continue
        if len(f[3])!=1 or len(f[4])!=1: continue
        row=[]; ok=True
        for g in gts:
            a=g.split(':')[0]
            if '|' not in a: ok=False;break
            x,y=a.split('|')
            if x not in '01' or y not in '01': ok=False;break
            row+=[int(x),int(y)]
        if not ok: continue
        pos.append(p); rows.append(row)
# ALU GT sits at INS; if window parse skipped it, re-read just that record
if alu_gt is None:
    with gzip.open(VCF,'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if f[1]==str(INS) and 'INS:ME:ALU' in f[4]:
                alu_gt=[f[9:][i].split(':')[0] for i in cix]; break
pos=np.array(pos)
H=np.array(rows,dtype=np.int8).reshape(len(pos),len(samples),2)
n=len(pos)
A=[0 if gt.replace('/','|').split('|')[0]=='1' else 1 for gt in alu_gt]   # ALU hap index per carrier
aluH=np.stack([H[:,cix[k],A[k]]     for k in range(6)],axis=1)             # n x 6 (founder hap)
altH=np.stack([H[:,cix[k],1-A[k]]   for k in range(6)],axis=1)            # n x 6 (control hap)
print(f"{n} SNPs in chr18:{LO:,}-{HI:,} (+/-{FLANK//1000} kb); block={BLOCK//1000} kb")

# ---------- per-block sharing ----------
def mean_pair_identity(mat, mask):
    sub=mat[mask]
    if sub.shape[0]==0: return np.nan
    vals=[(sub[:,i]==sub[:,j]).mean() for i,j in itertools.combinations(range(6),2)]
    return float(np.mean(vals))
def all6_frac(mat, mask):
    sub=mat[mask]
    if sub.shape[0]==0: return np.nan
    agree=np.array([ (sub[k]==sub[k][0]).all() for k in range(sub.shape[0]) ])
    return float(agree.mean())

edges=list(range(LO,HI+1,BLOCK))
tbl=[]
print(f"\n{'block(kb from INS)':22}{'nSNP':>6}{'ALU_pairID':>11}{'alt_pairID':>11}{'ALU-alt':>9}{'ALU_all6':>10}")
for b0,b1 in zip(edges[:-1],edges[1:]):
    mask=(pos>=b0)&(pos<b1)
    a=mean_pair_identity(aluH,mask); c=mean_pair_identity(altH,mask); a6=all6_frac(aluH,mask)
    rel0=(b0-INS)//1000; rel1=(b1-INS)//1000
    label=f"[{rel0:+d},{rel1:+d})"
    core=' *core*' if (b0<INS<b1 or (b0>=INS-200000 and b1<=INS+400000)) else ''
    d = (a-c) if (not np.isnan(a) and not np.isnan(c)) else np.nan
    print(f"{label:22}{int(mask.sum()):6d}{a:11.3f}{c:11.3f}{d:9.3f}{a6:10.3f}{core}")
    tbl.append((label,b0,b1,int(mask.sum()),a,c,d,a6))

# window-wide summary
mask_all=(pos>=LO)&(pos<=HI)
print(f"\nwhole +/-500 kb: ALU pairID={mean_pair_identity(aluH,mask_all):.3f}  "
      f"alt pairID={mean_pair_identity(altH,mask_all):.3f}  ALU all-6 frac={all6_frac(aluH,mask_all):.3f}")

with open(os.path.join(RES,'hap_block_sharing_500kb.tsv'),'w') as f:
    f.write("block\tstart\tend\tnSNP\tALU_pairID\talt_pairID\tALU_minus_alt\tALU_all6frac\n")
    for r in tbl: f.write("\t".join(str(x) for x in r)+"\n")

# ---------- FIGURE 1: block sharing decay (ALU vs non-ALU) ----------
mids=[(t[1]+t[2])/2/1000 - INS/1000 for t in tbl]
fig,ax=plt.subplots(figsize=(8,4.2))
ax.plot(mids,[t[4] for t in tbl],'o-',color='#c0392b',label='ALU haplotype (founder)')
ax.plot(mids,[t[5] for t in tbl],'s--',color='#2980b9',label='non-ALU haplotype (control)')
ax.axvline(0,color='k',lw=0.8,ls=':'); ax.set_xlabel('distance from insertion (kb)')
ax.set_ylabel('mean pairwise identity (6 carriers, 15 pairs)')
ax.set_title('CABYR: block-wise haplotype sharing, +/-500 kb'); ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'fig_block_sharing_500kb.png'),dpi=150); plt.close(fig)

# ---------- FIGURE 2: haplotype-level identity (12 carrier haplotypes) ----------
allhap=np.concatenate([aluH,altH],axis=1)             # n x 12
labels=[f"{c[:4]}:ALU" for c in CAR]+[f"{c[:4]}:alt" for c in CAR]
sub=allhap[mask_all]
Idn=np.zeros((12,12))
for i in range(12):
    for j in range(12): Idn[i,j]=(sub[:,i]==sub[:,j]).mean()
fig,ax=plt.subplots(figsize=(7,6))
im=ax.imshow(Idn,cmap='magma',vmin=0.5,vmax=1.0)
ax.set_xticks(range(12)); ax.set_xticklabels(labels,rotation=90,fontsize=7)
ax.set_yticks(range(12)); ax.set_yticklabels(labels,fontsize=7)
ax.set_title('Haplotype-level pairwise identity (+/-500 kb)\n6 ALU vs 6 non-ALU haplotypes')
fig.colorbar(im,label='fraction of SNPs identical'); fig.tight_layout()
fig.savefig(os.path.join(FIG,'fig_hap_identity_500kb.png'),dpi=150); plt.close(fig)

# ---------- FIGURE 3: clustering of the 12 carrier haplotypes ----------
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
D=1.0-Idn; np.fill_diagonal(D,0.0); D=(D+D.T)/2
Z=linkage(squareform(D,checks=False),method='average')
fig,(axd,axh)=plt.subplots(1,2,figsize=(12,5.2),gridspec_kw={'width_ratios':[1,1.25]})
dn=dendrogram(Z,labels=labels,orientation='left',ax=axd,color_threshold=0.15)
axd.set_title('Haplotype clustering (+/-500 kb)\n1 - identity, average linkage'); axd.set_xlabel('distance (1 - identity)')
order=dn['leaves'][::-1]
Ir=Idn[np.ix_(order,order)]
im=axh.imshow(Ir,cmap='magma',vmin=0.5,vmax=1.0)
axh.set_xticks(range(12)); axh.set_xticklabels([labels[o] for o in order],rotation=90,fontsize=7)
axh.set_yticks(range(12)); axh.set_yticklabels([labels[o] for o in order],fontsize=7)
axh.set_title('identity matrix, clustered order')
fig.colorbar(im,ax=axh,label='fraction of SNPs identical'); fig.tight_layout()
fig.savefig(os.path.join(FIG,'fig_hap_cluster_500kb.png'),dpi=150); plt.close(fig)
print("wrote figures/new_fig/fig_hap_cluster_500kb.png")

# block of 6 ALU-haps vs the rest, as numbers
alu_block=Idn[:6,:6][np.triu_indices(6,1)].mean()
alt_block=Idn[6:,6:][np.triu_indices(6,1)].mean()
cross=Idn[:6,6:].mean()
print(f"\nhaplotype-level mean identity (+/-500 kb): ALU-ALU={alu_block:.3f}  "
      f"alt-alt={alt_block:.3f}  ALU-alt(cross)={cross:.3f}")
print("wrote results/hap_block_sharing_500kb.tsv, figures/new_fig/fig_block_sharing_500kb.png, fig_hap_identity_500kb.png")
