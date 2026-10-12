"""Literal full-source nine-column critical normal on the third RR axis.

Chart p13=1, Rbar3=1; no p or D inverse. The exact ell=kappa e
coordinates and all first kappa/t terms are retained. This file alone
does not assert rank nine or exclude a lift.
"""
from sage.all import *
from pathlib import Path
import json,time,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));K=qf['ring'].base_ring();a,r,H,Z=K.gens();q=1+r;m=1+a*r
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));O=qf['ring']
P=PolynomialRing(K,names=('cs_x','cs_z','cs_w','cs_e0','cs_e1','cs_e2','cs_delta'))
x,z,w,e0,e1,e2,de=P.gens();es=vector(P,(e0,e1,e2));vals=(x,P.one(),z,w,P.zero(),P.zero(),P.zero(),P.zero(),P.one())
N=PolynomialRing(K,names=P.variable_names()+('cs_t','cs_kappa'));nv=N.gens();nx,nz,nw,ne0,ne1,ne2,nd,nt,nk=nv
ps=(nx,N.one(),nz,nw,(nt*nz-nk*ne2)/beta,nx*nw*nt/beta+nk*(nx*ne1-ne0)/alpha)
rs=(alpha/beta*nw*nt+nk*ne1,nt,N.one())
ck('exact chart ell0',rs[0]*ps[0]-alpha*rs[2]*ps[5]==nk*ne0)
ck('exact chart ell1',rs[0]*ps[1]-alpha/beta*rs[1]*ps[3]==nk*ne1)
ck('exact chart ell2',rs[1]*ps[2]-beta*rs[2]*ps[4]==nk*ne2)
def cp(f):return O({tuple(map(int,e))+(0,)*(9-len(e)):K(c) for e,c in f.dict().items()})
def layer(f,ij):return {tuple(map(int,e)):K(c) for e,c in cp(f).dict().items() if (sum(e[:6]),sum(e[6:]))==ij}
def jet(dd,cv=lambda c:K(c)):
    # Order: value, dp31, dp32, dRR1, dRR2. Skip every zero-support
    # monomial BEFORE mapping the original coefficient parameters.
    out=[P.zero() for unused in range(5)]
    for ex,c in dd.items():
        e=tuple(map(int,ex));e=e+(0,)*(9-len(e))
        for l,j in enumerate((None,4,5,6,7)):
            ee=list(e);fac=1
            if j is not None:
                fac=ee[j]
                if not fac:continue
                ee[j]-=1
            if any(n and not val for n,val in zip(ee,vals)):continue
            out[l]+=fac*cv(c)*prod(val**n for n,val in zip(ee,vals))
    return out
def jp(f):return jet(cp(f).dict())
def addjets(*parts):return [sum((item[j] for item in parts),P.zero()) for j in range(5)]
def timesjet(f,j):return [f*g for g in j]
def normal(jet0,extra):
    value,dp4,dp5,dr0,dr1=jet0
    jt=z/beta*dp4+x*w/beta*dp5+alpha*w/beta*dr0+dr1
    jk=-e2/beta*dp4+(x*e1-e0)/alpha*dp5+e1*dr0+extra
    return [value.derivative(v) for v in P.gens()]+[jt,jk]
off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];off=vector(O,[cp(f) for f in off]);D=matrix(O,qf['D']);C=matrix(O,qf['C']);DC=D*C;Dd=D*off
rows=[];boundary=[];labels=[];extras=[]
for j,R0 in enumerate(qf['full_R_sources']):
    Rf=cp(R0);jj=addjets(jp(Dd[j]),timesjet(de,jet(layer(Rf,(1,2)))),*[timesjet(es[k],jp(DC[j,k])) for k in range(3)])
    ext=de*jet(layer(Rf,(2,1)))[0]
    boundary.append(jj[0]);rows.append(normal(jj,ext));labels.append(('R',j));extras.append(ext)
Q=cp(qf['full_Q0']);B=list(map(cp,qf['Q0_L_cofactors']))
jj=addjets(jet(layer(Q,(2,1))),*[timesjet(es[k],jp(B[k])) for k in range(3)])
ext=jet(layer(Q,(3,0)))[0]+de*jet(layer(Q,(0,2)))[0]
boundary.append(jj[0]);rows.append(normal(jj,ext));labels.append(('Q0',0));extras.append(ext)
mod=load(str(p/'independent_delta_H_global_L_module.sobj'));old=mod['original_field'];Hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));sp=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'));pool=full['polynomial_pool'];ph=old.ring().hom([a,a*q/m,H/r,q*Z/r],K);cache={}
def cv(c):
    c=old(c)
    if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
    return cache[c]
def oldface(dag,ij):
    dd={}
    for sign,lref,rref in dag:
        for ex,c in pool[lref].items():
            for ey,b in pool[rref].items():
                e=tuple(int(v)+int(u) for v,u in zip(ex,ey))
                if (sum(e[:6]),sum(e[6:]))==ij:dd[e]=dd.get(e,old.zero())+sign*c*b
    return {e:c for e,c in dd.items() if c}
for j,it in enumerate(Hs['items']):
    jj=addjets(jet(it['H41_original'],cv),timesjet(de,jet(sp['H13_original'][j],cv)),*[timesjet(es[k],jet(it['L_cofactors_original'][k],cv)) for k in range(3)])
    h50=oldface(full['records'][j]['full_signed_product_DAG'],(5,0));h22=oldface(full['records'][j]['full_signed_product_DAG'],(2,2))
    ext=jet(h50,cv)[0]+de*jet(h22,cv)[0]
    boundary.append(jj[0]);rows.append(normal(jj,ext));labels.append(('H',j));extras.append(ext)
    print('H normal',j,'terms',[len(f.dict()) for f in rows[-1]],'seconds',round(time.monotonic()-st,2),flush=True)
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));el=load(str(p/'independent_delta_even_L_module.sobj'))
for j,it in enumerate(ev['items']):
    jj=addjets(jp(it['pure_p4']),*[timesjet(es[k],jp(it['A_p_quadratic'][k])) for k in range(3)])
    jj[0]+=K(it['square_unit'])*es[int(it['square_component'])]**2
    ext=de*(sum((jp(el['items'][j]['B_RR_linear'][k])[0]*es[k] for k in range(3)),P.zero())+jet(layer(it['full_source'],(2,1)))[0])
    boundary.append(jj[0]);rows.append(normal(jj,ext));labels.append(('even',int(it['Q_index'])));extras.append(ext)
oldbd=load(str(p/'independent_inactive_single_RR_boundary.sobj'));OP=oldbd['ring'];evaxis=OP.hom([x,P.one(),z,w,e0,e1,e2,de],P)
for j,f in enumerate(oldbd['items'][2]['rows']+oldbd['items'][2]['even']):ck('all original boundary row '+str(j),boundary[j]==evaxis(f))
graph=load(str(p/'independent_inactive_single_RR_even_graph.sobj'))['items'][2]['graph'];images=[evaxis(graph.get(v,v)) for v in OP.gens()]
gp=P.hom([x,z,w,images[4],images[5],images[6],images[7]],P)
mat=matrix(P,rows);mapped=matrix(P,28,9,[gp(f) for f in mat.list()]);bmap=[gp(f) for f in boundary]
for j in range(25):ck('all 25 genuine graph equations '+str(j),bmap[j]==0)
payload={'source_files':['independent_delta_subcritical_Q0_full_source.sobj','independent_ratio_global_R_derivative_degreefall.sobj','independent_delta_subcritical_even_full_source.sobj','independent_delta_even_L_module.sobj'],'ring':P,'chart_ring':N,'chart_p':ps,'chart_Rbar':rs,'boundary':boundary,'normal_matrix':mat,'column_order':['x','z','w','e0','e1','e2','delta','t','kappa'],'labels':labels,'kappa_native_extra':extras,'graph_images':[gp(v) for v in P.gens()],'graph_normal_matrix':mapped,'graph_boundary':bmap,'candidate_equations':(w*w-a*r,x*z-10*m,3*(x*x+z*z)-x*z),'checks':checks,'scope':__doc__}
stem='independent_critical_single_RR3_normal';tmp=p/(stem+'.tmp.sobj');save(payload,str(tmp),compress=False);os.replace(tmp,p/(stem+'.sobj'))
(p/(stem+'.json')).write_text(json.dumps({'status':'SOURCES_BUILT','checks':len(checks),'rows':28,'columns':9,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('DONE',len(checks),time.monotonic()-st,flush=True)
