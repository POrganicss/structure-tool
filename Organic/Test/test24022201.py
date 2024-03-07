# 将代码中的注释转换为中文
import re

def optimized_smile_divide(smile):
    # 定义 SMILES 解析的模式和结构元素
    prefixes = ["=", "/", "\\", "@", '#', "[", "]"]
    elements = ["Si", "Cl", "Br", "C", "N", "O", "F", "P", "S", "I"]
    suffix_number_pattern = re.compile(r"\d+")
    chirality = ["@@", "@", "\\", "/"]  # 手性表示

    # 辅助函数，用于模式匹配和字符串操作
    def match_and_append(pattern, s, unit):
        match = re.match(pattern, s)
        if match:
            unit.append(match.group(0))
            return s.replace(match.group(0), "", 1), True
        return s, False
    
    # 主解析逻辑
    
    
    smile_units = []
    while smile:
        unit = []
        # 匹配化学前缀（键、手性）
        for prefix in prefixes:
            if smile.startswith(prefix):
                unit.append(prefix)
                smile = smile[len(prefix):]

        # 匹配元素（包括双字母元素）
        for element in elements:
            if smile.startswith(element):
                unit.append(element)
                smile = smile[len(element):]

        # 匹配手性符号
        for chir in chirality:
            if smile.startswith(chir):
                unit.append(chir)
                smile = smile[len(chir):]

        # 匹配带有后缀数字的元素
        smile, matched = match_and_append(suffix_number_pattern, smile, unit)
        if not matched:
            # 处理由括号指示的嵌套结构
            if smile.startswith(("(", "[")):
                nested_structure = [smile[0]]
                smile = smile[1:]
                depth = 1
                while depth > 0:
                    char = smile[0]
                    smile = smile[1:]
                    nested_structure.append(char)
                    if char in "([":
                        depth += 1
                    elif char in ")]":
                        depth -= 1
                unit.append("".join(nested_structure))

        # 将单元组件合并成一个字符串，并添加到结果列表中
        if unit:  # 避免添加空单元
            smile_units.append("".join(unit))

    return smile_units

# 以下用于测试优化后的函数（用户决定是否运行）
# test_smile = "CC1=NC=CC2=C1NC3=CC=C(C(C)=O)C=C32"
# print(optimized_smile_divide(test_smile))
