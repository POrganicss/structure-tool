import __init__
import datetime
import re
import math


#提供gaussian运算参数的设置
class Code:

    # 根据函数输出生成命令代码
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

    def get(para: dict):
        fpara = {'nproc': str(para['nproc']), 'mem': para['mem'], 'name': para['name'],
                 'charge': str(para['charge']), 'spin': str(para['spin']), 'method': para['method']
                 }
        fpara['code'] = Code.getcode(para)

        return fpara
    
    def getpara(nproc='16', mem='55GB', name='proj', charge=0, spin='1', opt='calcfc', method='B3LYP/6-31g', additional='', solvent=''):
        
        code = {'opt': opt, 'method': method,
                'solvent': solvent, 'additional': additional}

        para = {'nproc': str(nproc), 'mem': mem, 'name': name,
                'charge': str(charge), 'spin': str(spin),
                'code': Code.getcode(code)}
        return para
    

    # 根据参数字典生成命令代码
    def getcode(code:dict):
        Code = '# '
        if code['opt']!='':
            Code=Code+'opt=('+code['opt']+") "
        if code['mixset']!='':
            Code=Code+code['method'].split('/')[0]+'/'+'genecp '
        else:
            Code=Code+code['method']
        if code['solvent']!='':
            Code = Code+' scrf=(smd,solvent=' +code['solvent']+') '
        Code = Code+code['additional']  
        return Code

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
    
    # 将代码生成参数字典
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
            elif '/' in unit:
                code['method'] = unit
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

        return Code.getopt(charge=para['charge'],
                           opt=code['opt'],
                           method=code['method'],
                           additional=code['additional'],
                           solvent=code['solvent'],
                           nproc=para['nproc'],
                           mem=para['mem'],
                           name=para['name'],
                           spin=para['spin']
                           )[0]

        # codes = {'job': 0, 'Basiset': 0, 'DFT': 0}
        # jobs = ['sp', 'opt', 'freq', 'IRC', 'scan']
        # Basisets = ['3-21g', '6-31g', '6-311g', 'sdd', 'genecp',
        #             'defTZVP', 'def2TZVP', 'def2TZVPP', 'def2QZVP']
        # DFTs = ['B3lyp', 'pbe1pbe', 'wb97xd', 'm062x']
        # for job in jobs:
        #     if job in code:

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
    def getname(name):
        date = str(datetime.date.today()).split('-')
        name_postfix = date[0][2:]+date[1]+date[2]
        return name[0]+name_postfix+str(name[1])
