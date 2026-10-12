"""Independent full archived Frobenius replay of the second Newton face."""
import sys
sys.argv=[__file__,'7','1']
exec(open(__file__.replace('trs_algebra_newton_full_source_reader.py','trs_algebra_newton_blowup.py')).read().split('cols={}')[0])
source=load(str(OUT/'trs_algebra_regular_frobenius.sobj'))
old=load(str(OUT/'trs_algebra_newton_blowup_1.sobj'))
cc=[]
for row,pair in enumerate(source['frobenius']):
 out=SR(0)
 for odd,pol in enumerate(pair):
  for mon,coef in sorted(pol.items()):
   term=evrat(coef)*KC**mon[0]*KB**mon[1]
   if odd:term*=U
   out+=term
   print('source',row,odd,mon,'seconds',round(time.time()-START,2),flush=True)
 out=cl(out)
 checks=[]
 for degree in range(6):
  assert nk(out[degree]-K(old['residuals']['C'+str(row)][degree]))==0
  checks.append(degree)
 cc.append(out);print('ROW_REPLAY_PASS',row,checks,flush=True)
CB=cl(cc[0]+l*cc[1]/(v+l*l))
assert nk(CB[5]-K(old['residuals']['CB'][5]))==0
payload={'status':'PASS_FULL_FROBENIUS_SOURCE_REPLAY','checks':13,'seconds':time.time()-START,'precision':int(CB.precision_absolute())}
(OUT/'trs_algebra_newton_full_source_reader.json').write_text(json.dumps(payload,indent=2));print(payload,flush=True)
signal.alarm(0)

