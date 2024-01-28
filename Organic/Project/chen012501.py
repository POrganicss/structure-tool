import __init__
from Tool.File import File
from Functions.Mechanism_Calculation import Mechanism_Calculation as ms
from Tool.Code import Code



# gjf = Xyz.togjf(Test,OPT)
# print(gjf)
# Mechanism_Calculation.gjf_optimize(gjf,'cstest01',2)
OPT1 = Code.set(opt='',nproc=12, additional='em=gd3bj',method='B3LYP/6-31g(d)', name='Pd01', mixset=[['Pd','Br'], ['SDD']])

ms.xyzs_optimize([],'Vis',OPT1)
