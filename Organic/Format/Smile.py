import time
import __init__
from Tool.File import *
from openbabel import openbabel
import re
from rdkit.Chem import AllChem
from rdkit import Chem
from openbabel import OBMol, OBConversion
from openbabel import openbabel, pybel
from Tool.Executor import Executor

class Smile:

    def smile_divide(smile):
        smile_units = []
        prefixes = ["=", "/", "\\", "@",'#',"[","]"]
        elements = ["Si", "Cl", "Br", "C", "N", "O", "F", "P", "S", "I"]
        suffixsnumber = r"\d+"
        chirality = "@@|@|\\|/"  # 添加手性表示
        
        def match_and_append(pattern, s):
            match = re.match(pattern, s)
            if match:
                unit.append(match.group(0))
                return s.replace(match.group(0), "", 1), True
            return s, False
        
        while smile:
            unit = []
            
            print('-----------1------------')
            # 匹配键型
            for prefix in prefixes:
                if smile.startswith(prefix):
                    unit.append(prefix)
                    smile = smile[len(prefix):]
                    break


            print('-----------2------------')
            # 匹配元素
            for element in elements:
                if smile.startswith(element):
                    unit.append(element)
                    smile = smile[len(element):]
                    break
            
            print('-----------3------------')
            # 匹配手性
            for chir in chirality.split('|'):
                if smile.startswith(chir):
                    unit.append(chir)
                    smile = smile[len(chir):]
                    break
            print('-----------4------------')
            # 匹配键型
            for prefix in prefixes:
                if smile.startswith(prefix):
                    smile, matched = match_and_append(suffixsnumber, smile[1:])
                    if matched:
                        unit.append(prefix)
                        break
            print('-----------5------------')
            # 匹配后缀数字
            smile, matched = match_and_append(suffixsnumber, smile)
            
            if not matched:
                # 匹配后缀括号部分
                while smile.startswith(("(", "[")):
                    bracket = []
                    bracket.append(smile[0])
                    smile = smile[1:]
                    i = 1
                    while i > 0:
                        if smile.startswith(("(", "[")):
                            bracket.append(smile[0])
                            smile = smile[1:]
                            i = i + 1
                        elif smile.startswith((")", "]")):
                            bracket.append(smile[0])
                            smile = smile[1:]
                            i = i - 1
                        else:
                            bracket.append(smile[0])
                            smile = smile[1:]
                    unit.append("".join(bracket))
            
            unit="".join(unit)
            print(unit)
            smile_units.append(unit)
        return smile_units

    def getsmiles(skeleton, positions, Fragments, brackets=[-1]):
        
        def addgroups(result, position, Fragments, bracket):
            def reorder_smiles(smiles, start_num=1):
                # 查找所有的环编号（一位数）
                numbers = re.findall(r'(\d)', smiles)
                
                # 去重并保持顺序
                unique_numbers_ordered = sorted(set(numbers), key=numbers.index)
                
                # 创建从旧编号到新编号的映射，从start_num开始
                number_mapping = {num: str(i + start_num) for i, num in enumerate(unique_numbers_ordered)}
                
                # 用新编号替换旧编号
                new_smiles = ''.join(number_mapping.get(char, char) for char in smiles)
                
                # 计算最大编号
                max_num = start_num + len(unique_numbers_ordered) - 1
                return new_smiles, max_num
            
            Results = []
            result,maxnum=reorder_smiles(result)
            frg1, frg2 =result[:position], result[position:]
            
            frg1, frg2 = frg1.strip(), frg2.strip()
            
            new_Fragments=[]
            for Fragment in Fragments:
                new_Fragments.append(reorder_smiles(Fragment,maxnum)[0])
            
            # for Fragment in new_Fragments:
            #     if Fragment.strip()=='':
            #         Results.append("".join(frg1 + frg2))
            #     else:
            #         Fragment = f"({Fragment})" if bracket == 1 else Fragment
            #         Results.append("".join(frg1 + Fragment + frg2))
            
            
            Results = [frg1 + (f"({Fragment})" if bracket == 1 else Fragment) 
                    + frg2 if Fragment.strip() else frg1 + frg2
                        for Fragment in new_Fragments
                        ]        
                    
            return Results
        if len(brackets) == 1 and brackets[0] == -1:
            brackets = [1]*len(positions)
            
        new_Fragments={}
        new_brackets={}
        for i,position in enumerate(positions):
            new_Fragments[position]=Fragments[i]
            new_brackets[position]=brackets[i]
        
        Result = [skeleton]
        positions.sort(reverse=True)
        for position in positions:
            num=len(Result)
            Result=Executor.ThreadExecutor(addgroups, Result, [position]*num, [new_Fragments[position]]*num, [new_brackets[position]]*num)
            Result=[item for sublist in Result for item in sublist]
            print(len(Result))
        return Result
    
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
        try:
            # Create a molecule from SMILES string
            mol = pybel.readstring("smi", smile)

            # Add hydrogens
            mol.OBMol.AddHydrogens()

            # Generate 3D coordinates
            mol.make3D()

            # Convert to XYZ format
            xyz_data = mol.write("xyz")
            
            return xyz_data

        except Exception as e:
            print(str(smile)+'无法转换成xyz坐标')
            # Handle exceptions and provide useful error information
            return f"Error during conversion: {str(e)}"

    def rdkitosdf(smile):
        mol = Chem.MolFromSmiles(smile)

        if mol is not None:
            # Generate 3D coordinates
            AllChem.Compute2DCoords(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)

            # Convert molecule to SDF format string
            sdf_string = Chem.MolToMolBlock(mol)

            return sdf_string
        else:
            return "Invalid SMILES input."
    
    def openbabeltosdf(smile):
        try:
            # Create a molecule from SMILES string
            mol = pybel.readstring("smi", smile)

            # Add hydrogens
            mol.OBMol.AddHydrogens()

            # Generate 3D coordinates
            mol.make3D()

            # Convert to SDF format
            sdf_data = mol.write("sdf")
            
            return sdf_data

        except Exception as e:
            print(str(smile)+'无法转换成SDF格式')
            # Handle exceptions and provide useful error information
            return f"Error during conversion: {str(e)}"
    
    def getsdf(smiles):
        return Executor.ThreadExecutor(Smile.openbabeltosdf,smiles)

