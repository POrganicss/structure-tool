import __init__
import re
from Tool.File import File
from Tool.Executor import Executor
import fitz
class Pdf:

    def getnum(line):
        # 使用正则表达式提取数字，包括小数和负数
        numbers = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
        # 将字符串数字转换为实际数字类型
        return numbers[-1]  
    
    def get_pages(path):
        def extract_text_from_page(page_num):
            return doc.load_page(page_num).get_text()
        with fitz.open(path) as doc:
            page_text=Executor.ThreadExecutor(extract_text_from_page,range(doc.page_count))
            #doc.close()
        return page_text
    
    def get_score(path):
        print('1')
        try:
            pages = Pdf.get_pages(path)
        except Exception as e:
            print(f"访问PDF页面时出错: {e}")
            return []
        print('2')
        scores = []
        for page in pages:
            lines = page.splitlines()
            score = {}

            # 改进的数据提取逻辑，增加了错误处理
            for i, line in enumerate(lines):
                try:
                    if "Molecule Name" in line:
                        score["Molecule Name"] = lines[i + 1].strip()
                    elif "Molecular Weight" in line:
                        score["Molecular Weight"] = lines[i + 1].strip()
                    elif "XLogP" in line:
                        score["XLogP"] = lines[i + 1].strip()
                    elif "PSA" in line:
                        score["PSA"] = lines[i + 1].strip()
                    elif "Heavy Atoms" in line:
                        score["Heavy Atoms"] = lines[i + 1].strip()
                    elif "Acceptor Count" in line:
                        score["Acceptor Count"] = lines[i + 1].strip()
                    elif "Donor Count" in line:
                        score["Donor Count"] = lines[i + 1].strip()
                    elif "Chelator Count" in line:
                        score["Chelator Count"] = lines[i + 1].strip()
                except IndexError:
                    # 处理标签后数据缺失的情况
                    continue

            # 改进的分数读取逻辑，增加了错误处理
            for line in lines:
                try:
                    if "Total Score" in line:
                        score["Total Score"] = Pdf.getnum(line)
                    elif "Shape" in line:
                        score["Shape"] = Pdf.getnum(line)
                    elif "Hydrogen Bond" in line:
                        score["Hydrogen Bond"] = Pdf.getnum(line)
                    elif "Protein Desolvation" in line:
                        score["Protein Desolvation"] = Pdf.getnum(line)
                    elif "Ligand Desolvation" in line:
                        score["Ligand Desolvation"] = Pdf.getnum(line)
                except ValueError:
                    # 处理getnum无法将分数转换为float的情况
                    continue

            # 验证分数并处理缺失数据
            try:
                # 检查score字典中是否存在必要的键，否则跳过这个分数
                required_keys = ['Shape', 'Hydrogen Bond', 'Protein Desolvation', 'Ligand Desolvation', 'Total Score']
                if not all(key in score for key in required_keys):
                    continue  # 如果缺少关键数据，则跳过此分数

                # 计算总分并进行验证
                _sum = sum(float(score[key]) for key in required_keys[:-1])  # 从总分计算中排除'Total Score'
                if abs(_sum - float(score["Total Score"])) <= 0.01:
                    scores.append(score)
                else:
                    # 如果总分不匹配，则打印错误详情
                    print(f"分子错误: {score.get('Molecule Name', '未知')}, 检测到分数不匹配。")
                    for key in required_keys:
                        print(f"{key} 数据: {score.get(key, 'N/A')}")
            except Exception as e:
                # 处理分数验证中的意外错误
                print(f"分数验证过程中出现意外错误: {e}")

        return scores


   