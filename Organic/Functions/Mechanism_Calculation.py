import __init__

from Tool.Verify import *
from Tool.Code import *
from Tool.File import *
from Tool.Datatransmission import *

from Format.Log import *
from Format.Xyz import *
from Format.Gjf import *
from Format.Smile import *

from Applications.Gaussian import *

from tqdm import tqdm
from Data import *
import concurrent.futures

class Mechanism_Calculation:
    
    def gjf_optimize(gjf: str, proname: str, index=1, solvent=''):
        Verify.Is([gjf, proname, index, solvent], [str, str, int, str])
        para = Gjf.getparameters(gjf)
        name = Code.getname(para['name']) #给文件加上日期后缀
        
        if index == 1:
            #直接进行优化
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
            
            #进行一第一次优化
            para['name'] = name+'_01'
            if 'freq' in para['additional']:
                para['additional'] = para['additional'].replace('freq', ' ').replace('  ',' ')
            xyz = Gaussian.run(gjf, proname, para['name'])
            
            
            #第二次优化的参数设置
            para['opt'] = para['opt'].replace(
                'calcfc', 'calcall')+',tight,maxstep=1,notrust'
            if 'freq' not in para['additional']:
                para['additional'] = para['additional']+' freq'
            if 'int=superfine' not in para['additional']:
                para['additional'] = para['additional']+' int=superfine'
            para['name'] = name+'_02'
            
            #进行第二次几何优化
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para['name'], 1)
            
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]
        elif index == 3:
            #进行一第一次优化
            para['name'] = name+'_01'
            if 'freq' in para['additional']:
                para['additional'] = para['additional'].replace('freq', ' ').replace('  ',' ')
            xyz = Gaussian.run(gjf, proname, para['name'])
            
            #第二次优化的参数设置
            para['opt'] = para['opt'].replace(
                'calcfc', 'calcall')+',tight,maxstep=1,notrust'
            if 'int=superfine' not in para['additional']:
                para['additional'] = para['additional']+' int=superfine'
            para['additional'] = para['additional']+'int=superfine'
            para['name'] = name+'_02'
            
            #进行第二次几何优化
            gjf = Xyz.togjf(xyz, Code.topara(para))
            xyz = Gaussian.run(gjf, proname, para['name'])
            
            #第三次几何优化的参数设置
            para['opt'] = ''
            if 'freq' not in para['additional']:
                para['additional'] = para['additional']+' freq'
            para['name'] = name+'_03'
            para['solvent'] = solvent
            
            #进行第三次几何优化
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para['name'], 1)
            
            if result[0]:
                return result[1], result[2]
            else:
                print('计算失败')
                return result[1]
    
    def gjfs_optimize(gjfs, proname, indexs, solvents):
        def gjf_optimizes(item):# 定义一个内部函数，用于处理单个任务
            # 调用 Mechanism_Calculation.gjf_optimize 函数处理单个任务
            Mechanism_Calculation.gjf_optimize(gjf=item[0], proname=item[1], 
                                               index=item[2], solvent=item[3])
        
        # 生成多线程任务参数
        # 创建一个 pronames 列表，长度与 gjfs 相同，每个元素都是 proname
        pronames = [proname]*len(gjfs)
        # 将四个列表打包成一个元组列表，每个元组包含一个任务的所有参数
        arguments_list = list(zip(gjfs, pronames, indexs, solvents))
        # 设置最大工作线程数
        max_workers = 128
        # 获取总任务数
        total_tasks = len(arguments_list)
        # 确保最大工作线程数不超过总任务数
        max_workers = min(max_workers, total_tasks) if max_workers else total_tasks
        
        
        # 创建一个 ThreadPoolExecutor 实例
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交初始任务集（最多为max_workers）
            futures = {executor.submit(gjf_optimizes, element): (
                element) for element in arguments_list[:max_workers]}
            running_tasks = max_workers

            # 在任务完成时对其进行处理
            results = []
            # 创建一个进度条
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
        except Exception:
            print('错误分子：'+smile)
            # messagebox.showinfo('错误提示', '分子式错误')
        else:
            return Mechanism_Calculation.xyz_optimize(xyz, proname, OPT, index, solvent)

    def smiles_optimize(smiles: list, proname: str, OPT: dict, charges=[], names=[], indexs=[], solvents=[]):
        Verify.Is([smiles, proname, OPT, charges,names,indexs],
                  [list, str, dict, list,list,list])
        xyzs = Smile.toxyzs(smiles)
        return Mechanism_Calculation.xyzs_optimize(xyzs, proname, OPT, charges, names, indexs, solvents)

