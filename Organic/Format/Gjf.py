import __init__
from tkinter import messagebox
from File import *
from Tool.Code import Code
from Tool.Verify import Verify

class Gjf:
    #所有的参数：name,mem,nproc,code,mixset
    
    # 从gjf中获取xyz坐标文件
    @staticmethod
    def toxyz(gjf):
        lines = gjf.splitlines()  # 将输入数据拆分为行
        xyzs = []
        xyz = []  # 初始化存储容器，用来存储xyz结果的list容器
        for line in lines:
            string = line.split()
            if len(string) == 4 and not Verify.isnum(string[0]) and Verify.isnums(string[1:4]):
                xyz.append(line)
            elif len(xyz) != 0:
                xyzs.append('\n'.join(xyz))
                xyz = []
        if len(xyzs) <= 1:
            return xyzs[0]
        else:
            return xyzs

    # 获取gjf中的所有参数
    def getcode(gjf):
        
        print()
    
    def getparameters(gjf):
        chunk=gjf.replace('\n\n\n','\n\n').split('\n\n')
        para = {}
        if len(chunk[0].splitlines())>3:
            para.update(Gjf.getfunctional(chunk[0]))
        else:
            messagebox.showinfo('gjf文件错误,命令行少于3行')
            
        if len(chunk[2].splitlines()) > 1:
            para.update(Gjf.getchargeandspin(chunk[2]))
        else:
            messagebox.showinfo('gjf文件错误,未发现电荷和自旋多重度')
        if len(chunk)>4:
            paramenters, sets=Gjf.getmixset(chunk[3])
            para.update(paramenters)
            para['code'] = para['code'].replace('genecp',sets)
        else:
            para['mixset']=''
        para.update(Code.tocode(para['code']))
        del para['code']
        return para
    
    def getfunctional(chunk):
        readline=chunk.splitlines()
        paramenters = {}
        for i, line in enumerate(readline):
            line = line.strip()
            if '%chk' in line:
                paramenters['name'] = line.split("=")[1].split(".")[0]
            elif '%mem' in line:
                paramenters['mem'] = line.split("=")[1]
            elif '%nproc' in line:
                paramenters['nproc'] = line.split("=")[1]
            elif '#' in line:
                paramenters['code'] = line
        return paramenters
         
    def getchargeandspin(chunk):
        readline=chunk.splitlines()
        paramenters = {}
        for line in readline:
            if len(line.split()) == 2 and Gjf.isnums(line.split()):
                paramenters['charge'] = line.split()[0]
                paramenters['spin'] = line.split()[1]
        return paramenters
    
    def getmixset(chunk):
        readline=chunk.splitlines()
        sets = readline[1]
        
        paramenters = {}
        
        elements=[]
        MIXSET=[]
        for i in range(1, int(len(readline)/3)):
            element = readline[i*3].split()
            element.remove('0')
            elements.append(element)
            MIXSET.append(readline[i*3+1])
        elements.append(MIXSET)
        paramenters['mixset'] = elements
            
        return paramenters, sets
    
    # 获取gjf中的其他参数
    def getname(gjf, name):
        readline = gjf.splitlines()
        for line in readline:
            if name in line:
                return line