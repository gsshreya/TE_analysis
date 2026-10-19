#!/usr/bin/env python3
"""Founder-haplotype analysis of the CABYR AluYa insertion (chr18:24,159,972).

Samples
Population   : DR_NGH_1_01 (the only community carrying the insertion).
Carriers     : read from the phased ALU genotype. Main analysis = unrelated carriers;
sensitivity analysis = all carriers (adds the carrier who is a child
of another carrier).
Null pool    : unrelated DR_NGH_1_01 individuals homozygous for the reference allele.

Window         : +/-500 kb around the insertion; phased biallelic SNPs only.
Identity       : fraction of SNPs at which two haplotypes carry the same allele.

Analyses
1. Within carriers (Fig S10B): carriers' ALU-bearing vs non-ALU haplotypes,
pairwise identity over +/-500 kb, average-linkage clustering.
2. Population clustering (Fig S10A): carriers' ALU haplotypes clustered with every
    null-pool haplotype in the shared core (the contiguous region around the
    nsertion where all ALU haplotypes carry the same allele). Tests:
    clade test  - size of the smallest clade containing all ALU haplotypes vs the
    same for k random null haplotypes (one per individual), R draws;
    repeated over the full +/-500 kb window (window not chosen from data);
3. Block-wise null (Fig S11): mean pairwise identity of the ALU haplotypes in ten
    100 kb blocks vs R draws of k null haplotypes (one per individual). Empirical
    p = (1 + #null >= obs)/(R + 1); Bonferroni over 10 blocks (alpha = 0.005).
    The carriers' non-ALU haplotypes are placed in the same null as a control.
"""
import gzip, os, sys, csv, itertools, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist, squareform

HERE=os.path.dirname(os.path.abspath(__file__))
VCF =os.environ.get('CABYR_POP_VCF','/path/to/hpc/cabyr/te_window_ALLsamples.vcf.gz')
META=os.environ.get('GI_META','/path/to/hpc/metadata_Jan_2025_noeths.csv')
OUT =os.environ.get('OUT_DIR',os.path.join(HERE,'outputs')); os.makedirs(OUT,exist_ok=True)
sys.path.insert(0,os.environ.get('FIG_STYLE_DIR',os.path.join(HERE,'..','Figures')))
import fig_style as S                                    # repo-wide figure style

POP='DR_NGH_1_01'; INS=24159972; FLANK=500_000; BLOCK=100_000
R=int(os.environ.get('NULL_R',10000)); SEED=11
LO,HI=INS-FLANK,INS+FLANK
C_ALU,C_ALT,C_NULL='#D55E00','#0072B2','0.80'            # ALU / non-ALU / null-pool haplotypes

class Tee:                                               # print to screen and to the log
    def __init__(self,path): self.f=open(path,'w'); self.o=sys.stdout
    def write(self,s): self.o.write(s); self.f.write(s)
    def flush(self): self.o.flush(); self.f.flush()
sys.stdout=Tee(os.path.join(OUT,'cabyr_founder_analysis.log'))
rng=np.random.default_rng(SEED)
print(f"CABYR founder analysis | VCF={VCF}\nmetadata={META} | R={R} | seed={SEED}\n")

# 1. samples
meta={r['IID']:r for r in csv.DictReader(open(META))}
in_pop=lambda s: s in meta and meta[s]['Code']==POP
is_child=lambda s: meta[s]['Child'].strip()=='Yes'

# 2. load window 
pos=[]; rows=[]; alu_gt=None
with gzip.open(VCF,'rt') as fh:
    for line in fh:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples=line.rstrip('\n').split('\t')[9:]
            keep=[i for i,s in enumerate(samples) if in_pop(s)]; ks=[samples[i] for i in keep]; continue
        f=line.rstrip('\n').split('\t'); p=int(f[1])
        if p==INS and 'INS:ME:ALU' in f[4]:
            alu_gt={samples[i]:f[9+i].split(':')[0].replace('/','|') for i in keep}; continue
        if p<LO or p>HI or len(f[3])!=1 or len(f[4])!=1: continue
        g=[f[9+i].split(':')[0] for i in keep]
        if any(len(a)!=3 or a[1]!='|' or a[0] not in '01' or a[2] not in '01' for a in g): continue
        pos.append(p); rows.append([int(a[0]) for a in g]+[int(a[2]) for a in g])
assert alu_gt is not None, 'ALU record not found at the insertion position'
pos=np.array(pos); n=len(pos); N=len(ks)
Hm=np.array(rows,dtype=np.int8)                          # SNPs x 2N (hap0 of all samples, then hap1)
hap=lambda s,h: Hm[:,ks.index(s)+h*N]

carriers=[s for s in ks if '1' in alu_gt[s]]
assert all(alu_gt[c] in ('0|1','1|0') for c in carriers), 'expected heterozygous carriers only'
A={c:(0 if alu_gt[c]=='1|0' else 1) for c in carriers}  # index of the ALU-bearing haplotype
car_u=[c for c in carriers if not is_child(c)]
pool=[s for s in ks if not is_child(s) and alu_gt[s]=='0|0']
P=np.stack([hap(s,h) for s in pool for h in (0,1)],1)   # null-pool haplotypes, 2 per individual
print(f"{POP}: {N} samples in VCF; {len(carriers)} carriers ({len(car_u)} unrelated; "
      f"excluded as offspring: {sorted(set(carriers)-set(car_u))}); null pool {len(pool)} "
      f"unrelated non-carriers = {P.shape[1]} haplotypes")
print(f"window chr18:{LO:,}-{HI:,}: {n} phased biallelic SNPs\n")

def alu_mat(cars):  return np.stack([hap(c,A[c])   for c in cars],1)
def alt_mat(cars):  return np.stack([hap(c,1-A[c]) for c in cars],1)
def null_cols(k):   # k distinct pool individuals, one random haplotype each, R draws -> column indices into P
    ind=np.array([rng.choice(len(pool),k,replace=False) for _ in range(R)])
    return 2*ind+rng.integers(0,2,size=ind.shape)

#  3. block-wise null (Fig S11)
edges=list(range(LO,HI+1,BLOCK))
masks=[(pos>=b0)&(pos<b1) for b0,b1 in zip(edges[:-1],edges[1:])]
blab=[f"{(b0-INS)//1000:+d} to {(b1-INS)//1000:+d}" for b0,b1 in zip(edges[:-1],edges[1:])]
W=np.stack([m/m.sum() for m in masks])                   # block averaging weights
def block_identity(counts,k):                            # mean pairwise identity per block from allele counts
    pid=(counts*(counts-1)+(k-counts)*(k-counts-1))/(k*(k-1))
    return W@pid
def block_test(cars):
    k=len(cars); cols=null_cols(k)
    cnt=np.zeros((n,R),dtype=np.int16)
    for j in range(k): cnt+=P[:,cols[:,j]]
    obs=block_identity(alu_mat(cars).sum(1).astype(float),k)
    alt=block_identity(alt_mat(cars).sum(1).astype(float),k)
    nul=block_identity(cnt.astype(float),k)
    p=(1+(nul>=obs[:,None]).sum(1))/(R+1); pct=(nul<alt[:,None]).mean(1)*100
    q=np.percentile(nul,[50,95,99],axis=1)
    return dict(k=k,obs=obs,alt=alt,p=p,alt_pct=pct,q=q)
blocks={'unrelated':block_test(car_u),'all':block_test(carriers)}

#  4. clustering + clade test (Fig S10A)
def core_window(M,tol=1):
    """contiguous SNP run around the insertion where all ALU haplotypes agree
    (a run ends after more than `tol` consecutive discordant SNPs)."""
    ci=int(np.argmin(np.abs(pos-INS))); agree=(M==M[:,[0]]).all(1)
    def run(step):
        v=0; last=k=ci
        while 0<=k<n:
            if agree[k]: last=k; v=0
            else:
                v+=1
                if v>tol: break
            k+=step
        return last
    return pos[run(-1)],pos[run(1)]
def clade_test(cars,mask):
    k=len(cars); X=np.concatenate([alu_mat(cars),P],1)[mask].T     # haplotypes x SNPs
    Z=linkage(pdist(X,'hamming'),'average'); _,nodes=to_tree(Z,rd=True)
    leafsets=[frozenset(nd.pre_order()) for nd in nodes]
    size=lambda ids: min(len(s) for s in leafsets if set(ids)<=s)
    obs=size(range(k)); nul=np.array([size(k+c) for c in null_cols(k)])
    return dict(k=k,X=X,Z=Z,obs=obs,null_med=float(np.median(nul)),p=(1+(nul<=obs).sum())/(R+1))
c0,c1=core_window(alu_mat(car_u)); core=(pos>=c0)&(pos<=c1)
full=np.ones(n,bool)
clade={(cs,wn):clade_test(cars,m) for cs,cars in (('unrelated',car_u),('all',carriers))
       for wn,m in (('core',core),('full',full))}
# exclusivity: identity of each null haplotype to the ALU consensus in the core
Acore=alu_mat(car_u)[core]; cons=(Acore.mean(1)>=0.5).astype(np.int8)
alu_to_cons=(Acore==cons[:,None]).mean(0); null_to_cons=(P[core]==cons[:,None]).mean(0)
pw_null=1-squareform(pdist(P[core].T,'hamming')); np.fill_diagonal(pw_null,0)
n_tight=int(((pw_null>=0.99).sum(1)>=len(car_u)-1).sum())

#  5. within-carrier clustering (Fig S10B)
WC=np.concatenate([alu_mat(car_u),alt_mat(car_u)],1).T   # 2k x SNPs, +/-500 kb
Iwc=1-squareform(pdist(WC,'hamming'))
k=len(car_u); iu=np.triu_indices(k,1)
wc_alu=Iwc[:k,:k][iu].mean(); wc_alt=Iwc[k:,k:][iu].mean(); wc_x=Iwc[:k,k:].mean()

print(f"--- Within carriers (+/-500 kb, {k} unrelated carriers, {k*(k-1)//2} pairs) ---")
print(f"mean pairwise identity: ALU-ALU {wc_alu:.3f} | non-ALU-non-ALU {wc_alt:.3f} | ALU-non-ALU {wc_x:.3f}\n")
print(f"--- Shared core: chr18:{c0:,}-{c1:,} ({(c0-INS)/1000:+.0f} to {(c1-INS)/1000:+.0f} kb; {int(core.sum())} SNPs) ---")
print(f"ALU haplotypes vs their consensus: min identity {alu_to_cons.min():.3f}")
print(f"null-pool haplotypes vs ALU consensus: max identity {null_to_cons.max():.3f}; "
      f">=0.99: {(null_to_cons>=0.99).sum()} of {P.shape[1]}")
print(f"null-pool haplotypes with >={len(car_u)-1} near-identical (>=0.99) partners in the core: {n_tight} of {P.shape[1]}\n")
print("--- Clade test: smallest clade containing all ALU haplotypes (leaves) ---")
for (cs,wn),r in clade.items():
    print(f"{cs:9} carriers (k={r['k']}), {wn:4} window: observed {r['obs']} (non-carriers inside: {r['obs']-r['k']}); "
          f"null median {r['null_med']:.0f}; p = {r['p']:.2g}")
for cs,b in blocks.items():
    print(f"\n--- Block-wise null, {cs} carriers (k={b['k']}); Bonferroni alpha = 0.005 ---")
    print(f"{'block (kb)':14}{'nSNP':>6}{'ALU':>8}{'non-ALU':>9}{'null med':>10}{'null 99%':>10}{'p':>9}{'non-ALU %ile':>14}")
    for i,l in enumerate(blab):
        print(f"{l:14}{int(masks[i].sum()):6d}{b['obs'][i]:8.3f}{b['alt'][i]:9.3f}{b['q'][0,i]:10.3f}"
              f"{b['q'][2,i]:10.3f}{b['p'][i]:9.4f}{b['alt_pct'][i]:13.0f}%")

with open(os.path.join(OUT,'table_block_null.tsv'),'w') as f:
    f.write("carrier_set\tk\tblock_kb\tnSNP\tALU_identity\tnonALU_identity\tnull_median\tnull_p95\tnull_p99\tp\tnonALU_percentile\n")
    for cs,b in blocks.items():
        for i,l in enumerate(blab):
            f.write(f"{cs}\t{b['k']}\t{l}\t{int(masks[i].sum())}\t{b['obs'][i]:.4f}\t{b['alt'][i]:.4f}\t"
                    f"{b['q'][0,i]:.4f}\t{b['q'][1,i]:.4f}\t{b['q'][2,i]:.4f}\t{b['p'][i]:.4g}\t{b['alt_pct'][i]:.1f}\n")
with open(os.path.join(OUT,'table_clustering.tsv'),'w') as f:
    f.write("carrier_set\tk\twindow\tnSNP\tobs_clade_size\tnoncarriers_in_clade\tnull_median_clade_size\tp\n")
    for (cs,wn),r in clade.items():
        f.write(f"{cs}\t{r['k']}\t{wn}\t{r['X'].shape[1]}\t{r['obs']}\t{r['obs']-r['k']}\t{r['null_med']:g}\t{r['p']:.4g}\n")
    f.write(f"\n# core window chr18:{c0}-{c1}; ALU min identity to consensus {alu_to_cons.min():.4f}; "
            f"max null identity to consensus {null_to_cons.max():.4f}; null haplotypes with >={len(car_u)-1} "
            f"near-identical partners {n_tight}/{P.shape[1]}\n")

# 7. figures 
# Figure S10
with matplotlib.rc_context(rc=matplotlib.rcParamsDefault):
    matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42})
    r=clade[('unrelated','core')]; kk=r['k']
    fig=plt.figure(figsize=(13,9.3))
    gs=fig.add_gridspec(2,1,height_ratios=[1,1.1],hspace=0.2)
    # ---- A: carriers' ALU haplotypes among all non-carrier haplotypes ----
    axA=fig.add_subplot(gs[0])
    dnA=dendrogram(r['Z'],ax=axA,no_labels=True,color_threshold=0.03,above_threshold_color='#bbbbbb')
    for xi,leaf in enumerate(dnA['leaves']):
        if leaf<kk: axA.plot(5+10*xi,0,marker='^',color='#d6371f',ms=9,clip_on=False,zorder=5)
    axA.set_ylabel('Haplotype distance',fontsize=12); axA.set_xticks([])
    axA.legend(handles=[plt.Line2D([],[],marker='^',ls='',color='#d6371f',ms=9,label='Carrier Alu Haplotypes')],
               loc='upper right',frameon=False,fontsize=13,handletextpad=0.2)
    fig.text(0.04,0.9,'A',fontsize=26,va='top')
    # ---- B: carriers' ALU vs non-ALU haplotypes ----
    gsB=gs[1].subgridspec(1,2,width_ratios=[1,1],wspace=0.6)
    axd=fig.add_subplot(gsB[0]); axh=fig.add_subplot(gsB[1])
    labd=[f"{c[:4]}:ALU" for c in car_u]+[f"{c[:4]}:alt" for c in car_u]
    labh=[f"{c[:4]}:Alu" for c in car_u]+[f"{c[:4]}:alt" for c in car_u]
    D=1-Iwc; np.fill_diagonal(D,0)
    Zb=linkage(squareform(D,checks=False),'average')
    dnB=dendrogram(Zb,labels=labd,orientation='left',ax=axd,leaf_font_size=14)
    axd.set_xlabel('Distance (1-identity)',fontsize=14); axd.tick_params(axis='x',labelsize=12)
    order=dnB['leaves'][::-1]                              # top-to-bottom, matching the dendrogram
    axh.imshow(Iwc[np.ix_(order,order)],cmap='magma',vmin=0.5,vmax=1.0,aspect='auto')
    im=axh.images[0]
    axh.set_xticks(range(len(order))); axh.set_xticklabels([labh[o] for o in order],rotation=90,fontsize=9)
    axh.set_yticks(range(len(order))); axh.set_yticklabels([labh[o] for o in order],fontsize=9)
    cb=fig.colorbar(im,ax=axh,fraction=0.05,pad=0.04); cb.set_label('Fraction of SNPs identical',fontsize=12)
    fig.text(0.04,0.47,'B',fontsize=26,va='top')
    for ext in ('pdf','png'):
        fig.savefig(os.path.join(OUT,f'Figure_S10.{ext}'),dpi=150,bbox_inches='tight')
    plt.close(fig)

S.set_style()
# Figure S11
b=blocks['unrelated']; mids=[(b0+b1)/2/1000-INS/1000 for b0,b1 in zip(edges[:-1],edges[1:])]
fig,ax=plt.subplots(figsize=(S.COL15,S.COL15*0.72))
ax.fill_between(mids,b['q'][0],b['q'][2],color='0.85',linewidth=0,label=f'{POP} non-carrier\nnull, median-99th')
ax.fill_between(mids,b['q'][0],b['q'][1],color='0.70',linewidth=0,label='null, median-95th')
ax.plot(mids,b['q'][0],color='0.35')
ax.plot(mids,b['obs'],'o-',color=C_ALU,ms=3,label='ALU haplotype')
ax.plot(mids,b['alt'],'s--',color=C_ALT,ms=3,label='non-ALU haplotype')
for x,pv in zip(mids,b['p']):
    if pv<0.005: ax.text(x,1.012,'*',ha='center',va='bottom')
ax.axvline(0,color='k',lw=0.6,ls=':'); ax.set_ylim(top=1.06)
ax.set_xlabel('Distance from insertion (kb)'); ax.set_ylabel('Mean pairwise identity')
ax.yaxis.grid(True,linestyle='--',alpha=0.4); ax.set_axisbelow(True)
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.2),ncol=2,frameon=False)
S.save(fig,os.path.join(OUT,'Figure_S11')); plt.close(fig)
