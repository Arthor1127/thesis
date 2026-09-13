# usage: python viz_sweep_trends.py <dir-with-shard_*.csv>   (default /tmp/sw)
import glob, csv, datetime
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows=[]
for f in glob.glob(__import__('sys').argv[1] + '/shard_*.csv') if len(__import__('sys').argv)>1 else glob.glob('/tmp/sw/shard_*.csv'):
    for r in csv.DictReader(open(f)):
        if (r.get('Phi0_per_mA') or '').strip():
            try:
                d={k:float(r[k]) for k in ['route_clearance_um','segment_length_um','squid_gap_um','arm_semi_length_um','arm_width_um','Phi0_per_mA','n_tet','n_ksp','seconds']}
                d['conv']=r.get('converged')=='True'; rows.append(d)
            except: pass
n=len(rows)
phi=np.array([r['Phi0_per_mA'] for r in rows]); sec=np.array([r['seconds'] for r in rows])
tet=np.array([r['n_tet'] for r in rows]); ksp=np.array([r['n_ksp'] for r in rows])
AX=['route_clearance_um','squid_gap_um','arm_semi_length_um','arm_width_um','segment_length_um']
SH={'route_clearance_um':'route_clearance','squid_gap_um':'squid_gap','arm_semi_length_um':'arm_semi_length','arm_width_um':'arm_width','segment_length_um':'segment_length'}
SB={'route_clearance_um':'wire-array gap [um]','squid_gap_um':'loop height [um]','arm_semi_length_um':'half loop width [um]','arm_width_um':'lead thickness [um]','segment_length_um':'inter-loop lead [um]'}
FULL={'route_clearance_um':[30,33,36,40],'squid_gap_um':[14.,19.2,24.4,29.6,34.8,40.],'arm_semi_length_um':[8,12,16,20,24],'arm_width_um':[3,4.5,6,8],'segment_length_um':[0,3,6,9,12]}
X={a:np.array([r[a] for r in rows]) for a in AX}
at=lambda a,v: phi[np.abs(X[a]-v)<1e-6]
INK,INK2,GRID,ACC,FIT='#1a1a2e','#6b6b83','#dcdce4','#3a6ea5','#c0392b'
plt.rcParams.update({'font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.linewidth':0.8})

lx=lambda a: np.log(X[a]+3.0) if a=='segment_length_um' else np.log(X[a])
M=np.column_stack([np.ones(n)]+[lx(a) for a in AX])
coef,*_=np.linalg.lstsq(M,np.log(phi),rcond=None)
logpred=M@coef; r2log=1-np.var(np.log(phi)-logpred)/np.var(np.log(phi)); pred=np.exp(logpred)
exps=dict(zip(AX,coef[1:]))
print("=== GLOBAL SCALING FIT  Phi0/mA ~ prefactor * prod x_a^p_a ===")
print(f"  prefactor exp(c0) = {np.exp(coef[0]):.3e}")
for a in AX: print(f"  {SH[a]:16} p = {exps[a]:+.2f}")
print(f"  log-space R^2 = {r2log:.3f}   median|rel err| = {np.median(np.abs(pred/phi-1))*100:.1f}%")

print("\n=== PER-AXIS 1-D FITS (per-value means, covered values) ===")
afit={}
for a in AX:
    vs=[v for v in FULL[a] if len(at(a,v))>=25]
    if len(vs)<2: continue
    ms=np.array([at(a,v).mean() for v in vs]); vs=np.array(vs,float)
    lin=np.polyfit(vs,ms,1); r2l=1-np.var(ms-np.polyval(lin,vs))/np.var(ms)
    pw=np.polyfit(np.log(vs+(3 if a=='segment_length_um' else 0)),np.log(ms),1)
    afit[a]=(lin,pw,r2l)
    print(f"  {SH[a]:16} linear {lin[0]:+.3f}/um {lin[1]:+.1f} (R2={r2l:.3f})   power x^{pw[0]:+.2f}")

print("\n=== COST FITS ===")
A2=np.column_stack([np.ones(n),tet*ksp/1e12]); (t0,t1),*_=np.linalg.lstsq(A2,sec,rcond=None)
r2t=1-np.var(sec-A2@[t0,t1])/np.var(sec)
print(f"  seconds = {t0:.0f} + {t1:.1f}*(n_tet*n_ksp/1e12)   R2={r2t:.3f}")
Bg=np.column_stack([np.ones(n)]+[X[a] for a in ['squid_gap_um','arm_semi_length_um','arm_width_um','segment_length_um']])
tc,*_=np.linalg.lstsq(Bg,tet,rcond=None); r2m=1-np.var(tet-Bg@tc)/np.var(tet)
print(f"  n_tet = {tc[0]/1e6:.2f}M + {tc[1]/1e3:.0f}k*gap + {tc[2]/1e3:.0f}k*semi + {tc[3]/1e3:.0f}k*aw + {tc[4]/1e3:.0f}k*seg   R2={r2m:.3f}")
plog=np.polyfit(np.log(tet),np.log(sec),1); print(f"  seconds ~ n_tet^{plog[0]:.2f}")

fig=plt.figure(figsize=(14,13),facecolor='white')
gs=fig.add_gridspec(4,15,height_ratios=[1.0,0.92,0.95,1.05],hspace=0.75,wspace=2.6,left=0.115,right=0.965,top=0.92,bottom=0.05)
ymax=max(np.percentile(at(a,v),92) for a in AX for v in FULL[a] if len(at(a,v)))*1.05

for i,a in enumerate(AX):
    ax=fig.add_subplot(gs[0,i*3:i*3+3]); vals=np.array(FULL[a],float)
    mean=np.array([at(a,v).mean() if len(at(a,v)) else np.nan for v in vals])
    lo=np.array([np.percentile(at(a,v),10) if len(at(a,v)) else np.nan for v in vals])
    hi=np.array([np.percentile(at(a,v),90) if len(at(a,v)) else np.nan for v in vals])
    cnt=[len(at(a,v)) for v in vals]
    ax.fill_between(vals,lo,hi,color=ACC,alpha=0.15,lw=0)
    ax.plot(vals,mean,'o',color=INK,ms=6,zorder=4)
    if a in afit:
        lin,pw,r2l=afit[a]; xx=np.linspace(vals.min(),vals.max(),60)
        yy=np.exp(np.polyval(pw,np.log(xx+(3 if a=='segment_length_um' else 0))))
        ax.plot(xx,yy,'-',color=FIT,lw=2,zorder=3)
        ax.text(0.04,0.94,f'$\\propto x^{{{pw[0]:+.2f}}}$',transform=ax.transAxes,fontsize=9,color=FIT,va='top',weight='bold')
    for v,m,c in zip(vals,mean,cnt):
        if c<25 and np.isfinite(m): ax.plot(v,m,'o',mfc='none',mec=FIT,ms=13,mew=1.8,zorder=5)
    ax.set_title(SH[a],fontsize=10.5,weight='bold',pad=11)
    ax.text(0.5,1.02,SB[a],transform=ax.transAxes,ha='center',fontsize=7.5,color=INK2)
    ax.set_ylim(0,ymax); sp=(vals.max()-vals.min())*0.08; ax.set_xlim(vals.min()-sp,vals.max()+sp)
    ax.grid(True,color=GRID,lw=0.7); ax.set_axisbelow(True); ax.spines[['top','right']].set_visible(False)
    (ax.set_ylabel('mean Phi0/mA\n(band 10-90%)',fontsize=8.5) if i==0 else ax.set_yticklabels([]))

fig.text(0.5,0.965,f'SQUID-array flux threading - trends & fits - partial sweep {n}/2400',ha='center',fontsize=13,weight='bold')
fig.text(0.5,0.94,f'checked {datetime.date.today()} ~20:00 -03  .  red curve = 1-D power fit  .  red ring = <25 samples  .  each panel averages the other 4 axes',ha='center',fontsize=8.2,color=INK2)

axE=fig.add_subplot(gs[1,0:7]); order=sorted(AX,key=lambda a:exps[a])
axE.barh(range(len(order)),[exps[a] for a in order],color=['#c0392b' if exps[a]<0 else '#2e8b57' for a in order],height=0.6)
for k,a in enumerate(order):
    axE.text(exps[a]+(0.04 if exps[a]>=0 else -0.04),k,f'{exps[a]:+.2f}',va='center',ha='left' if exps[a]>=0 else 'right',fontsize=9)
axE.set_yticks(range(len(order))); axE.set_yticklabels([SH[a] for a in order]); axE.axvline(0,color=INK2,lw=0.8)
axE.set_xlabel('fitted exponent p   in   Phi0/mA ~ prod (axis)^p',fontsize=8.5)
axE.set_title(f'Global multiplicative scaling fit  (log-space R2 = {r2log:.2f})',fontsize=10,weight='bold',pad=8)
axE.grid(True,axis='x',color=GRID,lw=0.7); axE.set_axisbelow(True); axE.spines[['top','right']].set_visible(False)

axP=fig.add_subplot(gs[1,8:12]); axP.scatter(phi,pred,s=8,color=ACC,alpha=0.5,lw=0)
lim=[0,phi.max()*1.05]; axP.plot(lim,lim,'-',color=INK2,lw=1); axP.set_xlim(lim); axP.set_ylim(lim)
axP.set_xlabel('measured Phi0/mA',fontsize=8.5); axP.set_ylabel('scaling-fit Phi0/mA',fontsize=8.5)
axP.set_title(f'fit vs data  (median err {np.median(np.abs(pred/phi-1))*100:.0f}%)',fontsize=9.5,weight='bold',pad=8)
axP.grid(True,color=GRID,lw=0.7); axP.set_axisbelow(True); axP.spines[['top','right']].set_visible(False)

axT=fig.add_subplot(gs[2,0:7]); axT.axis('off')
best=max(rows,key=lambda r:r['Phi0_per_mA']); nconv=sum(1 for r in rows if not r['conv'])
axT.text(0,1.0,
 f'SCALING   Phi0/mA  ~  {np.exp(coef[0]):.2e}  *\n'
 f'   route_clearance^{exps["route_clearance_um"]:+.2f} * squid_gap^{exps["squid_gap_um"]:+.2f} * arm_semi^{exps["arm_semi_length_um"]:+.2f}\n'
 f'   * arm_width^{exps["arm_width_um"]:+.2f} * (seg+3)^{exps["segment_length_um"]:+.2f}\n'
 f'   log-R2 {r2log:.2f}   median err {np.median(np.abs(pred/phi-1))*100:.0f}%\n\n'
 f'COST      seconds = {t0:.0f} + {t1:.1f}*(n_tet*n_ksp/1e12)   R2 {r2t:.2f}\n'
 f'          n_tet = {tc[0]/1e6:.2f}M + {tc[2]/1e3:.0f}k*arm_semi + {tc[3]/1e3:.0f}k*arm_width\n'
 f'                  + {tc[1]/1e3:.0f}k*squid_gap + {tc[4]/1e3:.0f}k*seg   R2 {r2m:.2f}\n'
 f'          seconds ~ n_tet^{plog[0]:.1f}\n\n'
 f'STATE     {n}/2400 ({100*n/2400:.0f}%)  .  Phi0/mA {phi.min():.1f}->{phi.max():.1f}  .  not-conv {nconv} ({100*nconv/n:.1f}%)\n'
 f'BEST      clr {best["route_clearance_um"]:.0f} seg {best["segment_length_um"]:.0f} gap {best["squid_gap_um"]:.0f} '
 f'semi {best["arm_semi_length_um"]:.0f} aw {best["arm_width_um"]:.1f} -> {best["Phi0_per_mA"]:.1f} (grid corner)\n'
 f'NOTE      segment_length still ~seg-0 only -> its fit is preliminary',
 fontsize=8.2,va='top',family='monospace',linespacing=1.5)

axC=fig.add_subplot(gs[2,8:13]); w=tet*ksp/1e12
axC.scatter(w,sec,s=8,color=ACC,alpha=0.45,lw=0); ww=np.linspace(w.min(),w.max(),50)
axC.plot(ww,t0+t1*ww,'-',color=FIT,lw=2)
axC.set_xlabel('n_tet * n_ksp / 1e12',fontsize=8.5); axC.set_ylabel('solve seconds',fontsize=8.5)
axC.set_title(f'cost ~ mesh x iterations  (R2 {r2t:.2f})',fontsize=9.5,weight='bold',pad=8)
axC.grid(True,color=GRID,lw=0.7); axC.set_axisbelow(True); axC.spines[['top','right']].set_visible(False)

pairs=[('squid_gap_um','arm_semi_length_um'),('squid_gap_um','arm_width_um'),('arm_semi_length_um','arm_width_um')]
vmax=float(np.nanpercentile(phi,99)); ims=[]
for j,(xa,ya) in enumerate(pairs):
    axh=fig.add_subplot(gs[3,j*5:j*5+5]); xv,yv=FULL[xa],FULL[ya]
    Z=np.full((len(yv),len(xv)),np.nan)
    for iy,vy in enumerate(yv):
        for ix,vx in enumerate(xv):
            m=phi[(np.abs(X[xa]-vx)<1e-6)&(np.abs(X[ya]-vy)<1e-6)]
            if len(m): Z[iy,ix]=m.max()
    im=axh.imshow(Z,origin='lower',aspect='auto',cmap='viridis',vmin=0,vmax=vmax); ims.append(im)
    axh.set_xticks(range(len(xv))); axh.set_xticklabels([f'{v:g}' for v in xv],fontsize=7)
    axh.set_yticks(range(len(yv))); axh.set_yticklabels([f'{v:g}' for v in yv],fontsize=7)
    axh.set_xlabel(SH[xa],fontsize=8); axh.set_ylabel(SH[ya],fontsize=8)
    for iy in range(len(yv)):
        for ix in range(len(xv)):
            if np.isfinite(Z[iy,ix]):
                axh.text(ix,iy,f'{Z[iy,ix]:.0f}',ha='center',va='center',fontsize=6.8,color='white' if Z[iy,ix]<vmax*0.5 else '#101015')
    if j==1: axh.set_title('Pairwise landscape - best Phi0/mA per cell (max over the other 3 axes)',fontsize=10,weight='bold',pad=12)
cb=fig.colorbar(ims[-1],ax=fig.axes[-3:],fraction=0.02,pad=0.03); cb.set_label('max Phi0/mA',fontsize=8)

out=(__import__('sys').argv[1] if len(__import__('sys').argv)>1 else '.') + '/sweep_trends.png'
fig.savefig(out,dpi=140,facecolor='white'); print('\nwrote',out)
