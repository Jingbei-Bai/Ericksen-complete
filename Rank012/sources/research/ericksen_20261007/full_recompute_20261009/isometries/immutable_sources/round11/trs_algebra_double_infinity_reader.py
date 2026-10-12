"""Exact full-source audit excluding B=infinity,L=infinity for p!=0."""
exec(open(__file__.replace('trs_algebra_double_infinity_reader.py','trs_algebra_frobenius_infinity_leading.py')).read().split('items=[]')[0])
lp=K.ring().gen(6);lr=R.gen(5);checks=[];records=[];corners={}
def check(label,cond):
 if not cond:raise AssertionError(label)
 checks.append(label)
def lcL(x):
 return nf(x.numerator().coefficient({lr:x.numerator().degree(lr)})/x.denominator().coefficient({lr:x.denominator().degree(lr)}))
for odd,pol in enumerate(d['frobenius'][1]):
 for mon,x in sorted(pol.items()):
  if not x:continue
  nl=int(x.numerator().degree(lp));dl=int(x.denominator().degree(lp))
  nn=poly_lead(x.numerator());dd=poly_lead(x.denominator())
  zl=-nn[0]+dd[0]-odd;rl=-nl+dl+odd
  dllead=int(dd[1].numerator().degree(lr)-dd[1].denominator().degree(lr))
  nllead=int(nn[1].numerator().degree(lr)-nn[1].denominator().degree(lr))
  tag=str((odd,mon))
  check('denominator_joint_corner '+tag,dllead==dl)
  check('numerator_z_bound '+tag,zl>=2)
  check('numerator_r_bound '+tag,rl>=-2)
  dvtx=lcL(dd[1])
  check('protected_denominator '+tag,all(nf(ff) in [v,j,p*p-v,v-p*p] for ff,ee in dvtx.factor()))
  record={'odd':odd,'pressure':mon,'z_bound':int(zl),'r_bound':int(rl),'denominator_vertex':str(dvtx.factor())}
  if zl==2 and rl==-2 and nllead==nl:
   zz,c,parity=rational_lead(x);corner=lcL(c)
   corners[(odd,mon)]=(corner,parity)
   record['contributes_corner']=True;record['corner_without_radical']=str(corner)
  records.append(record)
I1formal=(v**3+v*I2-j-P)/v**2
check('exactly_two_corner_terms',set(corners)=={(0,(0,0)),(1,(1,0))})
check('constant_corner',corners[(0,(0,0))]==(nf(-6*I1formal*p*p/(j*v)),0))
check('KC_corner',corners[(1,(1,0))]==(nf(7*p*p/(j*j)),1))
payload={'status':'PASS_ALL_PITCH_NONZERO_L_INFINITY','checks':len(checks),'seconds':time.time()-START,'records':records,
 'corner':'-6 I1 p^2/(j v) +/- 7 i KC p^2/(j sqrt(v(v-p^2)))',
 'norm_squared':'(p^4/j^2)[36 I1^2/v^2+49 KC^2/(v(v-p^2))] > 0',
 'scope':'For every fixed positive simple spectrum, arbitrary real KC,KB,p!=0 and v>p^2, the full C1 source has no complex algebraic branch with B->infinity and L->infinity at any relative rate. All four radical signs are covered.'}
(OUT/'trs_algebra_double_infinity_reader.json').write_text(json.dumps(payload,indent=2))
save({'checks':checks,'records':records,'corners':corners,'source_file':'trs_algebra_regular_frobenius.sobj'},str(OUT/'trs_algebra_double_infinity_certificate.sobj'))
signal.alarm(0);print(json.dumps({'status':payload['status'],'checks':len(checks),'seconds':time.time()-START}),flush=True)

