import datetime
import re
import math


#提供gaussian运算参数的设置
class Code:

    # 作为总的参数生成命令，根据函数输入生成命令代码
    @staticmethod  
    def set(nproc=16, mem='55GB', name='proj', charge=0, spin='1', opt='calcfc', method='B3LYP/6-31g', additional='', solvent='',mixset=''):
        
        code = {'opt': opt, 'method': method,
                'solvent': solvent, 'additional': additional,
                'mixset':mixset}

        para = {'nproc': str(nproc), 'mem': mem, 'name': name,
                'charge': str(charge), 'spin': str(spin),
                'code': Code.getcode(code),
                'mixset':mixset,
                'method':method,
                'solvent':solvent}
        return para
  
    # 根据参数字典生成命令代码 
    def getcode(code:dict):
        Code = '# '
        # 生成opt部分
        if code.get('opt','') != '': #判断逻辑为opt不为空则加到code中，否则不进行opt计算
            Code=Code+'opt=('+code['opt']+") "
        
        # 生成method部分
        if code['mixset']!='':
            Code=Code+code['method'].split('/')[0]+'/'+'genecp '
        else:
            Code=Code+code['method']+ ' '
            
        # 生成溶剂的的部分
        if code['solvent']!='':
            Code = Code+'scrf=(smd,solvent=' +code['solvent']+') '
        Code = Code+code.get('additional','')+' '
        return Code


    # 将修改后的代码重新生成参数字典
    def topara(para:dict):
        fpara = {'nproc': str(para['nproc']), 'mem': para['mem'], 'name': para['name'],
                 'charge': str(para['charge']), 'spin': str(para['spin']), 'mixset': para['mixset'], 'method': para['method']
                 
                 }
        fpara['code']=Code.getcode(para)
        return fpara
    
    # 将gjf命令解析拆分

    def resvoling(para, opt='-1', method='-1', additional='-1', solvent='-1'):

        units = para['code'].replace('# ', '').replace('freq', '').split(' ')
        code = {}
        code['additional'] = ''
        code['solvent'] = ''
        for unit in units:
            if 'opt' in unit:
                code['opt'] = unit.replace('opt=(', '').replace(')', '')
            elif 'scrf' in unit:
                code['solvent'] = unit.replace(
                    'scrf=(smd,solvent=', '').replace(')', '')
            else:
                code['additional'] = code['additional']+unit

        if opt != '-1':
            code['opt'] = opt
        if method != '-1':
            code['method'] = method
        if additional != '-1':
            code['additional'] = additional
        if solvent != '-1':
            code['solvent'] = solvent

        
        return code

        # codes = {'job': 0, 'Basiset': 0, 'DFT': 0}
        # jobs = ['sp', 'opt', 'freq', 'IRC', 'scan']
        # Basisets = ['3-21g', '6-31g', '6-311g', 'sdd', 'genecp',
        #             'defTZVP', 'def2TZVP', 'def2TZVPP', 'def2QZVP']
        # DFTs = ['B3lyp', 'pbe1pbe', 'wb97xd', 'm062x']
        # for job in jobs:
        #     if job in code:

   
    def get(para: dict):
        fpara = {'nproc': str(para['nproc']), 'mem': para['mem'], 'name': para['name'],
                 'charge': str(para['charge']), 'spin': str(para['spin']), 'method': para['method'],
                 'mixset':para['mixset']
                 }
        fpara.update(Code.resvoling(para))
        return fpara
    
      
    #将代码生成参数字典
    def tocode(order:str):
        units = order.replace('# ', '').replace('', '').split(' ')
        
        code = {}
        code['additional'] = ''
        code['solvent'] = ''
        
        for unit in units:
            if 'opt' in unit:
                code['opt'] = unit.replace('opt=(', '').replace(')', '')
            elif '/' in unit:
                code['method'] = unit
            elif 'scrf=' in unit:
                code['solvent'] = unit.replace(
                    'scrf=(smd,solvent=', '').replace(')', '')
            else:
                code['additional'] = code['additional']+' '+unit
        return code
    

    # 获取gjf中的前后命令

    def getcode1(gjf):
        elements = ['Si', 'Cl', 'Br', 'C', 'N', 'O', 'F', 'P', 'S', 'I']
        number = re.match(suffixsnumber, smile)
        code = []
        glines = gjf.readline()
        for gline in glines:
            if '#' in gline:
                code.append(gline+'\n')
            elif "" in gline:
                print()
        print()

    # 根据xyz赋予命令行生成gjf文件
    def getopt(charge=0, opt='calcfc', method='B3LYP/6-31g', additional='', solvent='', nproc='16', mem='55GB', name='proj', spin='1',mixset=''):
        OPT = []
        OPT1 = {'nproc': nproc, 'mem': mem, 'name': name,
                'charge': str(charge), 'spin': spin,  # 价态和自旋多重度
                # 计算指令
                'code': '# opt=('+opt+') '+method+' '+additional,
                'mixset':mixset,
                'method': method
                }

        if solvent != '':
            OPT1['code'] = OPT1['code']+' scrf=(smd,solvent='+solvent+')'
            
        if mixset !='':
            OPT1['code']=''

        OPT2 = {'nproc': '16', 'mem': '55GB', 'name': 'proj',
                'charge': str(charge), 'spin': '1',  # 价态和自旋多重度
                'code': '# opt=calcfc freq B3LYP/6-31g em=gd3bj'  # 计算指令
                }

        OPT3 = {'nproc': '16', 'mem': '55GB', 'name': 'proj',
                'charge': str(charge), 'spin': '1',  # 价态和自旋多重度
                # 计算指令
                'code': '# opt=calcfc freq B3LYP/6-311g(d,p) em=gd3bj,'
                }

        OPT4 = {'nproc': '16', 'mem': '55GB', 'name': 'proj',
                'charge': str(charge), 'spin': '1',  # 价态和自旋多重度
                'code': '# opt=calcfc freq wB97XD/def2TZVP,'  # 计算指令
                }

        OPT.append(OPT2)
        OPT.append(OPT3)
        OPT.append(OPT4)

        return OPT1, OPT

    def getts(charge=0, opt='calcfc,ts,noeigen', method='B3LYP/6-31g', additional='', solvent='', nproc='16', mem='55GB', name='proj', spin='1'):
        TS = []
        TS1 = {'nproc': nproc, 'mem': mem, 'name': name,
               'charge': str(charge), 'spin': spin,  # 价态和自旋多重度
               # 计算指令
               'code': '# opt=('+opt+') freq '+method+' '+additional
               }
        if solvent != '':
            TS1['code'] = TS1['code']+' scrf=(smd,solvent='+solvent+')'

        TS2 = {'nproc': '16', 'mem': '55GB', 'name': 'proj',
               'charge': str(charge), 'spin': '1',  # 价态和自旋多重度
               # 计算指令
               'code': '# opt=(calcfc,ts,noeigen) freq wB97XD/def2TZVP'
               }
        TS3 = {'nproc': '16', 'mem': '55GB', 'name': 'proj',
               'charge': str(charge), 'spin': '1',  # 价态和自旋多重度
               # 计算指令
               'code': '# opt=(calcfc,ts,noeigen,gdiis) freq wB97XD/def2TZVP'
               }
        TS4 = {'nproc': '16', 'mem': '55GB', 'name': 'proj',
               'charge': str(charge), 'spin': '1',  # 价态和自旋多重度
               # 计算指令
               'code': '# opt=(calcall,ts,noeigen,gdiis) int=superfine wB97XD/def2TZVP'
               }
        TS.append(TS2)
        TS.append(TS3)
        TS.append(TS4)
        return TS1, TS

    # 根据要求生成sets
    def getsets():
        print()

    # 根据要求生成%cpu内容 nmber:需要的总核数
    def cpu(number, start=0):
        ends = number+start-1
        return str(start)+'-'+str(ends), ends+1

    # 根据要求平均分配cpu核数 nmber:需要的组数
    def cpusets(number, start=0, sum=128):
        sets = Code.nprocsets(number, sum)
        cpusets = []
        for set in sets:
            r, start = Code.cpu(set, start)
            cpusets.append(r)
        return cpusets

    # 根据要求平均分配cpu核数 nmber:需要的组数
    def nprocsets(number, sum=128):
        sets = []
        if number >= 128:
            for i in range(number):
                sets.append('1')
            return sets

        return Code.sets(number, sum)[0]

    # 根据要求平均分配cpu核数 nmber:需要的组数
    def memsets(number, sum=450):
        sets = []
        if number >= 128:
            for i in range(number):
                sets.append('3GB')
            return sets

        memsets = Code.sets(number, sum)[0]
        sets = []
        for mem in memsets:
            sets.append(str(mem)+'GB')
        return sets

    # 根据要求平均分配资源数，nmber:需要的组数，sum:总资源数
    def sets(number, sum):
        sets = []

        i = math.floor(sum/number)
        k = sum/number-i
        s = 0  # 已经使用的内存
        for j in range(number):
            med = math.floor(i*(j+1)+k*j-s+1)
            sets.append(med)
            s = s+med
        if s > sum:
            sets[number-1] = sets[number-1]-1
            s = s-1
        return sets, s

    # 获取文件名
    def getname(name: str):
        date = str(datetime.date.today()).split('-')
        name_postfix = date[0][2:]+date[1]+date[2]
        return name+'_'+name_postfix


    def allocate_resources(total_cores, total_memory, num_tasks):
        # 初始化每个任务的资源分配
        allocations = [{'cores': 0, 'memory': 0} for _ in range(num_tasks)]

        # Step 1: 均匀分配核心
        cores_per_task = total_cores // num_tasks
        for allocation in allocations:
            allocation['cores'] += cores_per_task
        
        # Step 2: 随机分配剩余的核心
        remaining_cores = total_cores % num_tasks
        while remaining_cores > 0:
            allocation = random.choice(allocations)
            allocation['cores'] += 1
            remaining_cores -= 1

        # Step 3: 分配内存
        memory_per_core = total_memory / total_cores
        for allocation in allocations:
            allocation['memory'] = allocation['cores'] * memory_per_core

        # 生成Gaussian所需的资源字符串
        resource_strings = []
        for allocation in allocations:
            resource_string = f"%NProcShared={allocation['cores']}\n%Mem={int(allocation['memory'])}GB"
            resource_strings.append(resource_string)
        
        return resource_strings