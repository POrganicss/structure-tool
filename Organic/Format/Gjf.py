import __init__
from tkinter import messagebox
from Tool.Code import Code
from Tool.Verify import Verify
from Tool.File import File
class Gjf:
    #所有的参数：name,mem,nproc/cpu,code,xyz,mixset
    #其中基组的存储形式sets=[[['C'，'H','O'],'6-31g'],[['Pd'],'SDD']]
    #读取形式：elements=[['C'，'H','O'],['Pd']],mixset=['6-31g','SDD']
    
    # 从gjf中获取xyz坐标文件
    @staticmethod
    def toxyz(gjf):
        """
        将 Gaussian 输入文件格式（GJF）转换为 XYZ 格式。

        参数:
            gjf (str): Gaussian 输入文件的内容。

        返回:
            str: XYZ 格式的坐标。或者list类型:多个分子的xyz坐标
        """
        
        # 分割 GJF 内容为行
        lines = gjf.splitlines()
        xyzs = []  # 存储所有的 XYZ 块
        xyz = []   # 存储当前的 XYZ 块

        # 遍历每一行
        for line in lines:
            string = line.split()

            # 检查行是否符合 XYZ 块的格式要求
            if len(string) == 4 and not Verify.isnum(string[0]) and Verify.isnums(string[1:4]):
                xyz.append(line)
            elif len(xyz) != 0:
                # 如果当前行不符合要求但已经有坐标存在，则将当前 XYZ 块添加到 'xyzs' 列表中
                xyzs.append('\n'.join(xyz))
                xyz = []

        # 如果只有一个 XYZ 块，则返回它；否则返回所有 XYZ 块
        if len(xyzs) <= 1:
            return xyzs[0]
        else:
            return xyzs

    # 从gjf中获取除坐标外的所有参数
    def getcode(gjf):
        
        def get_functional(chunk):
            readline=chunk.splitlines()
            paramenters = {}
            for line in readline:
                line = line.strip()
                if '%chk' in line:
                    paramenters['name'] = line.split("=")[1].split(".")[0]
                elif '%mem' in line:
                    paramenters['mem'] = line.split("=")[1]
                elif '%nproc' in line:
                    paramenters['nproc'] = line.split("=")[1]
                elif '%cpu' in line:
                    paramenters['cpu'] = line.split("=")[1]
                elif '#' in line:
                    paramenters['code'] = line
                elif len(line.split()) == 2 and Verify.isnums(line.split()):
                    paramenters['charge'] = line.split()[0]
                    paramenters['spin'] = line.split()[1]
            return paramenters
     
        def get_mixset(chunk):######待进一步修改
            lines=chunk.strip().splitlines()
            
            elements=[]#存储类似elements=[['C'，'H','O'],['Pd']]
            MIXSET=[]#存储类似mixset=['6-31g','SDD']
            
            elements.append(lines[0].split())
            elements[0].remove('0')
            MIXSET.append(lines[1])
            
            for i in range(1, int(len(lines)/3)):
                element = lines[i*3].split()
                element.remove('0')
                elements.append(element)
                MIXSET.append(lines[i*3+1])
            
            return {'mixset':[[elements[0],MIXSET[0]],[elements[1],MIXSET[1]]]}
        
        # 对输入文件进行预处理
        chunk=gjf.replace('\n\n\n','\n\n').split('\n\n')
        
        # 获取坐标区域之前的内容
        Header=chunk[0]+'\n'+chunk[2].splitlines()[0]+'\n'
        
        # 获取坐标区域之后的内容
        Footer=chunk[3]
        
        # 定义一个字典，用于存储所有的参数
        para = {}
        
        # 将坐标区域之前的内容添加到字典中
        para.update(get_functional(Header))
        
        # 将坐标区域之后的内容添加到字典中
        para.update(get_mixset(Footer))
        return para

gjf=File.getdata('D:\\BOrganic\\Desktop\\test.gjf')
result=Gjf.toxyz(gjf)
print(result[2])
