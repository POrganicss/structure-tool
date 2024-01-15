import os
import sys

path = os.getcwd()
path = path.replace("\\", "/")+'/Organic'
sys.path.append(path+'/Format')
sys.path.append(path+'/Tool')
sys.path.append(path+'/Applications')
sys.path.append(path+'/Functions')
sys.path.append(path+'/resources')
from tqdm import tqdm
from Data import *
from Verify import *
from Code import *
from Gaussian import *
from Log import *
from Xyz import *
from Gjf import *
from Smile import *
from File import *
from Datatransmission import *
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, wait

class Mechanism_Calculation:

    def gjf_optimize(index, name, gjf, proname, para, solvent=None):
        def run_gaussian(gjf, proname, name, additional='', opt_settings=''):
            para['name'] = name
            xyz = Gaussian.run(gjf, proname, para['name'])
            para['opt'] = para['opt'].replace('calcfc', 'calcall') + opt_settings
            if 'freq' not in para['additional']:
                para['additional'] += ' freq'
            if additional not in para['additional']:
                para['additional'] += ' ' + additional
            return xyz

        if index == 2:
            para['name'] = name + '_01'
            xyz = run_gaussian(gjf, proname, para['name'], 'int=superfine', ',tight,maxstep=1,notrust')
            para['name'] = name + '_02'
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para['name'], 1)
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]

        elif index == 3:
            para['name'] = name + '_01'
            xyz = run_gaussian(gjf, proname, para['name'], 'int=superfine')
            para['name'] = name + '_02'
            gjf = Xyz.togjf(xyz, Code.topara(para))
            xyz = run_gaussian(gjf, proname, para['name'], 'freq', ',tight,maxstep=1,notrust')
            para['name'] = name + '_03'
            para['solvent'] = solvent
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para['name'], 1)
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]


    def gjf_optimize(gjf: str, proname: str, index=1, solvent=''):
        Verify.Is([gjf, proname, index, solvent], [str, str, int, str])
        para = Gjf.getparameters(gjf)
        name = para['name']
        
        if index == 1:
            if 'freq' not in para['additional']:
                para['additional'] = para['additional']+' freq'
            xyz = Gjf.toxyz(gjf)
            
            gjf = Xyz.togjf(xyz, Code.topara(para))
            
            result = Gaussian.run(gjf, proname, para['name'], 1)
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]
        elif index == 2:
            para['name'] = name+'_01'
            xyz = Gaussian.run(gjf, proname, para['name'])
            para['opt'] = para['opt'].replace(
                'calcfc', 'calcall')+',tight,maxstep=1,notrust'
            if 'freq' not in para['additional']:
                para['additional'] = para['additional']+' freq'
            if 'int=superfine' not in para['additional']:
                para['additional'] = para['additional']+' int=superfine'
            para['name'] = name+'_02'
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para['name'], 1)
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]
        elif index == 3:
            para['name'] = name+'_01'
            xyz = Gaussian.run(gjf, proname, para['name'])
            para['opt'] = para['opt'].replace(
                'calcfc', 'calcall')+',tight,maxstep=1,notrust'
            if 'int=superfine' not in para['additional']:
                para['additional'] = para['additional']+' int=superfine'
            para['additional'] = para['additional']+'int=superfine'
            para['name'] = name+'_02'
            gjf = Xyz.togjf(xyz, Code.topara(para))
            xyz = Gaussian.run(gjf, proname, para['name'])
            para['opt'] = ''
            if 'freq' not in para['additional']:
                para['additional'] = para['additional']+' freq'
            para['name'] = name+'_03'
            para['solvent'] = solvent
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para['name'], 1)
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]
    
    def gjf_optimizes(item):
            Mechanism_Calculation.gjf_optimize(gjf=item[0], proname=item[1], index=item[2], solvent=item[3])
    
    def gjfs_optimize(gjfs, proname, indexs, solvents):
        pronames = [proname]*len(gjfs)
        arguments_list = list(zip(gjfs, pronames, indexs, solvents))
        max_workers = 128
        total_tasks = len(arguments_list)
        max_workers = min(max_workers, total_tasks) if max_workers else total_tasks
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交初始任务集（最多为max_workers）
            futures = {executor.submit(Mechanism_Calculation.gjf_optimizes, element): (
                element) for element in arguments_list[:max_workers]}
            running_tasks = max_workers

            # 在任务完成时对其进行处理
            results = []
            with tqdm(total=total_tasks, desc="Processing tasks") as pbar:
                for completed_future in concurrent.futures.as_completed(futures):
                    try:
                        # 检索已完成任务的结果
                        result = completed_future.result()
                        results.append(result)

                        pbar.update(1)  # 更新进度条

                        # 如果还有剩余的参数，请提交新任务
                        if running_tasks < total_tasks:
                            new_args = arguments_list[running_tasks]
                            new_future = executor.submit(
                                Mechanism_Calculation.gjf_optimizes, *new_args)
                            futures[new_future] = new_args
                            running_tasks += 1
                    except Exception as e:
                        # 处理任务执行期间引发的异常
                        print(f"Exception: {e}")

        return results
    
    def xyz_optimize(xyz: str, proname: str, OPT: dict, index=1, solvent=''):
        Verify.Is([xyz, proname, OPT, index, solvent],
                  [str, str, dict, int, str])

        gjf = Xyz.togjf(xyz, OPT)
        return Mechanism_Calculation.gjf_optimize(gjf, proname, index, solvent)

    def xyzs_optimize(xyzs: list, proname: str, OPT: dict, charges=[],names=[],indexs=[],solvents=[]):
        Verify.Is([xyzs, proname, OPT, charges, names, indexs],
                  [list, str, dict, list, list, list])
        
        if indexs == [] or len(indexs) < len(xyzs):
            indexs = [1]*len(xyzs)
            
        if solvents == [] or len(solvents) < len(xyzs):
            if OPT['solvent'] != '':
                solvents = [OPT['solvent']]*len(xyzs)
            else:
                solvents = ['']*len(xyzs)
            
        gjfs = Xyz.togjfs(xyzs, OPT, charges,names)
        return Mechanism_Calculation.gjfs_optimize(gjfs, proname, indexs,solvents)
        #return Mechanism_Calculation.multrungjf(items)

    def smile_optimize(smile: str, proname: str, OPT: dict, index=1, solvent=''):
        Verify.Is([smile, proname, OPT, index, solvent],
                  [str, str, dict, int, str])
        try:
            xyz = Smile.toxyz(smile)
        except:
            print('错误分子：'+smile)
            # messagebox.showinfo('错误提示', '分子式错误')
        else:
            return Mechanism_Calculation.xyz_optimize(xyz, proname, OPT, index, solvent)

    def smiles_optimize(smiles: list, proname: str, OPT: dict, charges=[], names=[], indexs=[], solvents=[]):
        Verify.Is([smiles, proname, OPT, charges,names,indexs],
                  [list, str, dict, list,list,list])
        xyzs = Smile.toxyzs(smiles)
        return Mechanism_Calculation.xyzs_optimize(xyzs, proname, OPT, charges, names, indexs, solvents)

    def getitems(gjfs: list, proname: str, indexs: list, solvent=''):
        items = []
        for gjf, index in zip(gjfs, indexs):
            items.append([gjf, proname, index, solvent])
        return items

    def xyz_irc(xyz, proname, OPT):
        print()


# gjf = Xyz.togjf(Test,OPT)
# print(gjf)
# Mechanism_Calculation.gjf_optimize(gjf,'cstest01',2)
OPT1 = Code.set(opt='',nproc=12, additional='em=gd3bj',method='B3LYP/6-31g(d)', name='Pd01', mixset=[['Pd','Br'], ['SDD']])

Mechanism_Calculation.xyzs_optimize([],'Vis',OPT1)
