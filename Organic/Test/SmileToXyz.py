import __init__

from rdkit.Chem import AllChem
from rdkit import Chem
from openbabel import openbabel
from Tool.File import *
class Toxyz:
    
    def rdkitoxyztwo(smile):
    
        # 将smile转换为分子对象
        mol = Chem.MolFromSmiles(smile)
        
        if mol is not None:
            # 为分子计算3D坐标
            AllChem.Compute2DCoords(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)  # 使用随机种子以确保可复现性
            
            # 获取原子坐标
            conf = mol.GetConformer()
            coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            
            # 构建xyz格式字符串
            xyz_string = f"{mol.GetNumAtoms()}\n\n"
            for i, coord in enumerate(coords):
                xyz_string += f"{mol.GetAtomWithIdx(i).GetSymbol()} {coord.x} {coord.y} {coord.z}\n"
            
            return xyz_string
        else:
            return None
        
    # 将smiles转换成xyz文件
    def rdkitoxyz(smile):

        moldh = Chem.MolFromSmiles(smile)

        # 补上所有氢原子
        mol = Chem.AddHs(moldh)

        # 生成初始3D构象
        AllChem.EmbedMolecule(mol, randomSeed=1)

        # 利用分子力场方法进行能量最小化优化
        # AllChem.MMFFOptimizeMolecule(mol)

        # 获取分子电荷数
        # Charge = Chem.rdmolops.GetFormalCharge(mol)

        # 输出优化后的3D结构xyz文件
        xyz = Chem.MolToXYZBlock(mol)
        xyzs = xyz.splitlines()
        return "\n".join(xyzs[2:])
    
    
    def openbabeltoxyz(smile):
        File.save(smile, path+'/temp02.smi')

        conv = openbabel.OBConversion()  # 使用openbabel模块中的OBConversion函数，用于文件格式转换的

        # 输入需要转换的文件的名字，以及定义转换后文件的文件名
        conv.OpenInAndOutFiles(path+'/temp02.smi', path+'/temp02.xyz')
        conv.SetInAndOutFormats("smi", "xyz")  # 定义转换文件前后的格式
        conv.Convert()  # 执行转换操作
        conv.CloseOutFile()  # 转换完成后关闭转换后的文件，完成转换

        xyz = File.getdata(path+'/temp02.xyz')
        os.remove(path+'/temp02.xyz')
        return xyz



Toxyz.rdkitoxyztwo()