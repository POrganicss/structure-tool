import re
import __init__
from Format.Xyz import Xyz
from Format.Smile import Smile

def smile_divide(smile):
    
    def match_bond(smile,unit):
        prefixes = ["=", "/", "\\", "@",'#',"[","]"]
        for prefix in prefixes:
            if smile.startswith(prefix):
                unit.append(prefix)
                smile = smile[len(prefix):]
                return smile,unit
        return smile,unit
    
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
        print("smile:"+smile)
        # 匹配后缀数字
        smile, matched = match_and_append(suffixsnumber, smile)
        
        #if not matched:
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
        print("unit:"+unit)
        smile_units.append(unit)
        print(smile_units)
    return smile_units
smile1='=C2[C@]3(C#N)O[C@H](CO[P@@](OC4=CC=CC=C4)(N[C@H](C(OCC(CC)CC)=O)C)=O)[C@@H](O)[C@H]3O'
smile='C1=CC=CC=C1'

#xyz1=Smile.rdkitoxyz(smile)
#xyz2=Smile.openbabeltoxyz(smile)
#print(Xyz.toxyz(xyz2))
smile_divide(smile)


