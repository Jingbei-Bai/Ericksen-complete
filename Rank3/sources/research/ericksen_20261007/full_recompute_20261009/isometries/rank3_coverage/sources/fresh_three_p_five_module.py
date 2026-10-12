"""Extract the entire five-generator source matrix from fresh full H2 normals."""
from sage.all import *
from pathlib import Path
p=Path(__file__).absolute().parent
d=load(str(p/'independent_subcritical_single_RR_three_p_second_normal.sobj'))
assert d['matrix'].dimensions()==(28,9),'Complete second normal requires all 28 source rows and 9 coordinates'
E=d['ring'];M=d['matrix'].matrix_from_rows_and_columns(range(25),(2,3,4,6,8))
save({'ring':E,'matrix':M,'source_columns':(2,3,4,6,8),'scope':__doc__},str(p/'independent_subcritical_single_RR_three_p_five_module.sobj'),compress=False)
print('FULL H2 FIVE MODULE REBUILT',M.dimensions(),flush=True)
