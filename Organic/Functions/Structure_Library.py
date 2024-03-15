import os
import sys
import time
path = os.getcwd()
path = path.replace("\\", "/")+'/Organic'
sys.path.append(path+'/Format')
sys.path.append(path+'/Tool')
sys.path.append(path+'/Applications')
sys.path.append(path+'/Functions')
from Code import *
from Gaussian import Gaussian
from Log import *
from Xyz import *
from Gjf import *
from Smile import *

class Structure_Library:
    # projectname：项目名称
    # smile：骨架结构
    # positions：需要取代的位置
    # fragments：与位置对应的分子碎片，需要逐个单独定义
    # brackets：该分子碎片是否为取代基，为取代基则为1，不为取代基则为0，默认为取代基
    def getlibrary(projectname: str, smile: str, positions: list, fragments: list, brackets=[-1]):

        smiles = Smile.getsmiless(
            smile, positions, fragments, brackets)  # 建立smile小分子库

        paramenters = {'nproc': '128', 'mem': '400GB', 'charge': '0', 'spin': '1',
                       'code': '# opt pm6',
                       'name': 'proj'}

        # 获取分子的立体信息，并将其转换为xyz文件
        xyzs = Gaussian.runsmiles(smiles, projectname, paramenters)

        sdfs = Xyz.tosdfs(xyzs)  # 将xyz文件转换为sdf文件

        File.save(sdfs, path+'/all.sdf')  # 保存sdf文件

    def getmultlibrary(projectname: str, smile: str, positions: list, fragments: list, brackets: list):
        print('------------gjf生成中------------')
        smiles = Smile.getsmiless(
            smile, positions, fragments, brackets)  # 建立smile小分子库
        paramenters = {'charge': '0', 'spin': '1',
                       'code': '# opt pm6', 'name': 'proj','mixset':''}

        # 获取分子的立体信息，并将其转换为xyz文件
        Gaussian.runmultsmiles(smiles, projectname, paramenters)

        xyzs = []
        for i in range(len(smiles)):
            path2 = path+'/Temp/'+Code.getname(['proj', str(i+1)])+'.xyz'
            xyz = File.read(path2)
            xyzs.append(xyz)
            os.remove(path2)

        sdfs = Xyz.tosdfs(xyzs)  # 将xyz文件转换为sdf文件
        File.save(sdfs, path+'/all.sdf')  # 保存sdf文件
        
    def getmult(projectname: str, smile: str, positions: list, fragments: list, brackets: list):
        print('------------gjf生成中------------')
        smiles = Smile.getsmiless(
            smile, positions, fragments, brackets)  # 建立smile小分子库


        paramenters = {'charge': '0', 'spin': '1',
                       'code': '# opt b3lyp', 'name': 'proj'}
        
        items=Gaussian.getitemss(smiles,projectname,paramenters) #打包内容
        
        print(items)
        
        
        print('------------gjf生成完成，开始计算------------')
        result=Gaussian.multrunsmiless(items)
        print(result)

