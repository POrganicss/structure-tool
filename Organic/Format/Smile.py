import __init__
from Tool.File import *
from openbabel import openbabel
import re
from rdkit.Chem import AllChem
from rdkit import Chem
from openbabel import OBMol, OBConversion
from openbabel import openbabel, pybel

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
     
    def getsmile(Result, position, Fragments, bracket):
        # 按照最大连接数修改碎片数值
        def modifyfrag(frag, num):
            if num == -1:
                return frag
            for i in range(9, 0, -1):
                if str(i) in frag:
                    frag = frag.replace(str(i), str(i+num-1), 2)
            return frag
        
        Results = []
        for result in Result:
            smile_units = Smile.smile_divide(result)
            frg1, frg2 ="".join(smile_units[:position-1]), "".join(smile_units[position-1:])
            digits = [int(char) for char in result if char.isdigit()]
            maxnum = max(digits) + 1 if digits else -1
            
            for Fragment in Fragments:
                Fragment = modifyfrag(Fragment, maxnum)
                
                Fragment = f"({Fragment})" if bracket == 1 else Fragment
                Results.append("".join(frg1 + Fragment + frg2))
        return Results

    # 计算效率公式 单核单线程：
    def getsmiless(skeleton, positions, Fragments, brackets=[-1]):
        # 将某一部分与位置建立联系
        def getconnect(positions, fragments):
            for i in range(len(positions)-1):
                for j in range(len(positions)-i-1):
                    if positions[j] < positions[j+1]:
                        positions[j], positions[j+1] = positions[j+1], positions[j]
                        fragments[j], fragments[j+1] = fragments[j+1], fragments[j]
            Fragments = {}
            for i in range(len(positions)):
                Fragments[positions[i]] = fragments[i]
            return Fragments

        if len(brackets) == 1 and brackets[0] == -1:
            brackets = [1 for _ in range(len(positions))]
        Fragments = getconnect(positions, Fragments)
        brackets = getconnect(positions, brackets)
        Result = [skeleton]
        for position in positions:
            Result = Smile.getsmile(Result, position,
                                    Fragments[position], brackets[position]
                                    )
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

    def toxyz(smile):
        try:
            return Smile.openbabeltoxyz(smile)
        except Exception :
            return Smile.rdkitoxyz(smile)
    
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
        # 创建一个 Open Babel 分子对象
        ob_mol = openbabel.OBMol()
        
        # 将 SMILES 转换为 Open Babel 分子对象
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("smi", "sdf")
        conv.ReadString(ob_mol, smile)

        # 将 Open Babel 分子对象转换为 SDF 格式字符串
        sdf_string = conv.WriteString(ob_mol)

        return sdf_string
    
    
    # 将一些smiles转换成xyz文件
    def toxyzs(smiles):
        xyzs = []
        for smile in smiles:
            xyzs.append(Smile.toxyz(smile))
        return xyzs

    def getfragments(valence, total=0, fragment=[]):
        Elements = {
            1: "H",
            2: "He",
            3: "Li",
            4: "Be",
            5: "B",
            6: "C",
            7: "N",
            8: "O",
            9: "F",
            10: "Ne",
            11: "Na",
            12: "Mg",
            13: "Al",
            14: "Si",
            15: "P",
            16: "S",
            17: "Cl",
            18: "Ar",
            19: "K",
            20: "Ca",
            21: "Sc",
            22: "Ti",
            23: "V",
            24: "Cr",
            25: "Mn",
            26: "Fe",
            27: "Co",
            28: "Ni",
            29: "Cu",
            30: "Zn",
            33: "As",
            35: "Br",
            44: "Ru",
            45: "Rh",
            46: "Pd",
            47: "Ag",
            49: "In",
            53: "I",
            55: "Cs",
            74: "W",
            75: "Re",
            77: "Ir",
            78: "Pt",
            79: "Au",
            82: "Pb",
        }
        elements = {1: ["", "F", "Cl", "Br", "I"],
                    2: ["O", "S"], 3: ["N"], 4: ["C"]}

        if total == 1:
            for i in range(1, 4):
                fragment = Smile.getconnect(i, "0", fragment)
            return fragment
        if valence == 1:
            for key in elements:
                for i in range(len(elements[key])):
                    fragment.append(elements[key][i])
            # fragments=['','N','C','O','S','F','Cl','Br','I']
            """elif valence==2:
            for key in elements:
                Keys=list(elements.keys)
                
                for Key in Keys:
                    for i in range(len(elements[key])):
                        for j in range(len(elements[Key])):
                            if key!=1:
                                fragment.append(elements[key][i]+elements[Key][j])"""
        elif valence == 2:
            fragment = [
                "CC",
                "C=C",
                "NN",
                "OO",
                "CN",
                "C#N",
                "NC",
                "N#C",
                "OC",
                "CO",
                "C=O",
                "ON",
                "NO",
                "N=O",
            ]
        elif valence == 3:
            fragment = [
                "CCC",
                "C(C)C",
                "C=CC",
                "CC=C",
                "C(C)=C",
                "C#CC",
                "CC#C",
                "CCO",
                "OCC",
                "COC",
                "C(O)C",
                "CC=O",
                "OC=C",
                "OC#C",
                "CCN",
                "CNN",
                "NCC",
                "CNC",
                "C(N)C",
                "N(C)C",
                "CC#N",
                "C=NN",
                "NC=C",
                "N(=O)=O",
                "C1CC1",
            ]
        elif valence == 4:
            fragment = [
                "CCCC",
                "C(C)CC",
                "CC(C)C",
                "C(C)(C)C",
                "C=CCC",
                "C=CC=C",
                "CC=CC",
                "CCC=C",
                "C#CCC",
                "CC#CC",
                "CCC#C",
                "C(=C)CC",
                "C(C)=CC",
                "C(C)C=C",
                "C(C)C#C",
                "C=C(C)C",
                "CC(=C)C",
            ]
        return fragment

    
