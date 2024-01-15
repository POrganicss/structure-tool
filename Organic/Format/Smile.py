from File import *
from openbabel import openbabel
import re
from rdkit.Chem import AllChem
from rdkit import Chem
import os
import sys
path = os.getcwd()
path = path.replace("\\", "/")+'/Organic'
sys.path.append(path+'/Format/')
sys.path.append(path+'/Tool/')
sys.path.append(path+'/Applications/')
sys.path.append(path+'/Functions/')


class Smile:

    # smile左侧切出一个小单元结构
    def smile_unit(smile):
        unit = []
        prefixs = ["=", "/", "\\", "@"]
        elements = ["Si", "Cl", "Br", "C", "N", "O", "F", "P", "S", "I"]
        suffixsnumber = r"\d+"  # 数字

        # 匹配键型
        for prefix in prefixs:
            if smile.find(prefix) == 0:
                unit.append(prefix)
                smile = smile.replace(prefix, "", 1)
                break

        # 匹配元素
        for element in elements:
            if smile.find(element) == 0:
                unit.append(element)
                smile = smile.replace(element, "", 1)
                break

        # 匹配键型
        for prefix in prefixs:
            if smile.find(prefix) == 0:
                number = re.match(suffixsnumber, smile[1])
                if number is not None:
                    unit.append(prefix)
                    smile = smile.replace(prefix, "", 1)
                    break

        # 匹配后缀数字
        number = re.match(suffixsnumber, smile)
        if number is not None:
            unit.append(number.group(0))
            smile = smile.replace(number.group(0), "", 1)

        # 匹配后缀括号部分
        while smile.find("(") == 0:
            bracket = []
            bracket.append(smile[0])
            smile = smile.replace(smile[0], "", 1)
            i = 1
            while i > 0:
                if smile.find("(") == 0:
                    bracket.append(smile[0])
                    smile = smile.replace(smile[0], "", 1)
                    i = i + 1
                elif smile.find(")") == 0:
                    bracket.append(smile[0])
                    smile = smile.replace(smile[0], "", 1)
                    i = i - 1
                else:
                    bracket.append(smile[0])
                    smile = smile.replace(smile[0], "", 1)
            unit.append("".join(bracket))
        return "".join(unit), smile

    # 将smile分解成小单元
    def smile_divide(smile):
        smile_units = []
        while smile != "":
            smile_unit, smile = Smile.smile_unit(smile)
            smile_units.append(smile_unit)
            if smile == "":
                break
        return smile_units

    # 选择性对分子碎片加括号
    def addbracket(Fragment, bracket):
        if Fragment != "" and bracket == 1:
            Fragment = "(" + Fragment + ")"
        elif Fragment != "" and bracket == 0:
            return Fragment
        return Fragment

    # 将smile的小单元按照确定位置进行重新组合
    def keletsonsplit(smile_units, position):
        return "".join(smile_units[:position-1]), "".join(smile_units[position-1:])

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

    # 获得最大连接数
    def getmaxcon(smile):
        for i in range(9, 0, -1):
            if str(i) in smile:
                return i + 1
        return -1

    # 按照最大连接数修改碎片数值
    def modifyfrag(frag, num):
        if num == -1:
            return frag
        for i in range(9, 0, -1):
            if str(i) in frag:
                frag = frag.replace(str(i), str(i+num-1), 2)
        return frag

    def getsmile(Result, position, Fragments, bracket):
        Results = []
        for result in Result:
            smile_units = Smile.smile_divide(result)
            frg1, frg2 = Smile.keletsonsplit(smile_units, position)
            maxnum = Smile.getmaxcon(result)
            for Fragment in Fragments:
                Fragment = Smile.modifyfrag(Fragment, maxnum)
                Fragment = Smile.addbracket(Fragment, bracket)
                Results.append("".join(frg1 + Fragment + frg2))
        return Results

    # 计算效率公式 单核单线程：time=exp(0.8951*log3(x)-9.3175)
    def getsmiless(skeleton, positions, Fragments, brackets=[-1]):
        if len(brackets) == 1 and brackets[0] == -1:
            brackets = Smile.getbrackets(len(positions))
        Fragments = Smile.getconnect(positions, Fragments)
        brackets = Smile.getconnect(positions, brackets)
        Result = [skeleton]
        for position in positions:
            Result = Smile.getsmile(Result, position,
                                    Fragments[position], brackets[position]
                                    )
        return Result

    # 将smiles转换成xyz文件
    def toxyz(smile):
        print(smile)
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
        File.tofile(smile, path+'/temp02.smi')

        conv = openbabel.OBConversion()  # 使用openbabel模块中的OBConversion函数，用于文件格式转换的

        # 输入需要转换的文件的名字，以及定义转换后文件的文件名
        conv.OpenInAndOutFiles(path+'/temp02.smi', path+'/temp02.xyz')
        conv.SetInAndOutFormats("smi", "xyz")  # 定义转换文件前后的格式
        conv.Convert()  # 执行转换操作
        conv.CloseOutFile()  # 转换完成后关闭转换后的文件，完成转换

        xyz = File.getdata(path+'/temp02.xyz')
        os.remove(path+'/temp02.xyz')
        return xyz

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

    def SmilesToSdf(smiles, names, ind=0):
        mols = []

        for i in range(len(smiles)):
            moldh = Chem.MolFromSmiles(smiles[i])

            try:
                mol = Chem.AddHs(moldh)
            except Exception as e:
                print(smiles[i] + "\n")
                print(e)

            # 生成初始3D构象
            AllChem.EmbedMolecule(mol, randomSeed=1)

            # 利用分子力场方法进行能量最小化优化

            AllChem.MMFFOptimizeMolecule(mol)
            mols.append(mol)

        """ writer = Chem.SDWriter('CS'+str(ind)+'.sdf')
        writer.SetProps(['LOGP', 'MW'])
    
        for i, mol in enumerate(mols):
            mw = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            mol.SetProp('MW', '%.2f' %(mw))
            mol.SetProp('LOGP', '%.2f' %(logp))
            mol.SetProp('_Name', 'cs'+names[i])
            writer.write(mol)
        writer.close() """

        return mols

    def getbrackets(number):
        li = []
        for i in range(number):
            li.append(1)
        return li
