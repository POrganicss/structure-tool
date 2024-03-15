import os
import __init__
from Format.Pdb import Pdb
from Format.Sdf import Sdf
from Tool.Datatransmission import LocalCommand
from Tool.File import File
from Format.Pdf import Pdf
import os
import shlex
class Openeye:
    #分子对接
    def run(path, Ligands_name, Protein_name, nproc):
        # 生成openeye命令
        command = []
        # -----------------------新建项目文件夹-----------------------#
        File.create_files(path, "fixpka_oeomega", "hybrid", "report")
        File.copy(path, os.path.join(path, "fixpka_oeomega"), f"{Ligands_name}.sdf")
        File.copy(path, os.path.join(path, "hybrid"), f"{Protein_name}.pdb")
        
        # -----------------分别对小分子和蛋白质进行处理----------------#
        Sdf.get_oeomega(os.path.join(path, "fixpka_oeomega"), Ligands_name, nproc=nproc)
        File.cp(path, "fixpka_oeomega", "hybrid", f"{Ligands_name}_fixpka_oeomega.oeb.gz")
        Pdb.get_receptor(os.path.join(path, "hybrid"), Protein_name)
        
        # -----------------进行小分子对接----------------#
        hybrid_cmd = (
            f"hybrid -mpi_np {nproc} -receptor {shlex.quote(os.path.join(path, 'hybrid', f'{Protein_name}_clean_receptor.oeb.gz'))} "
            f"-dbase {shlex.quote(os.path.join(path, 'hybrid', f'{Ligands_name}_fixpka_oeomega.oeb.gz'))} "
            f"-docked_molecule_file {shlex.quote(os.path.join(path, 'hybrid', f'{Ligands_name}_{Protein_name}hybrid_docked_molecule.sdf'))} "
            f"-undocked_molecule_file {shlex.quote(os.path.join(path, 'hybrid', f'{Ligands_name}_{Protein_name}hybrid_undocked_molecule.sdf'))} "
            f"-score_file {shlex.quote(os.path.join(path, 'hybrid', 'hybrid_score.txt'))} "
            f"-report_file {shlex.quote(os.path.join(path, 'hybrid', 'hybrid_report.txt'))} "
            f"-settings_file {shlex.quote(os.path.join(path, 'hybrid', 'hybrid_settings.param'))} "
            f"-status_file {shlex.quote(os.path.join(path, 'hybrid', 'hybrid_status.txt'))} && "
        )
        docking_report_cmd = (
            f"docking_report -docked_poses {shlex.quote(os.path.join(path, 'hybrid', f'{Ligands_name}_{Protein_name}hybrid_docked_molecule.sdf'))} "
            f"-receptor {shlex.quote(os.path.join(path, 'hybrid', f'{Protein_name}_clean_receptor.oeb.gz'))} "
            f"-report {shlex.quote(os.path.join(path, 'hybrid', 'report.pdf'))}\n"
        )

        # 将命令添加到命令列表
        command.append(hybrid_cmd)
        command.append(docking_report_cmd)

        # 执行命令
        LocalCommand.execute_command(''.join(command))
        print("对接结果生成完成！")

    def get_score(path):
        print('1')
        try:
            # 假设Pdf类在其他地方定义，并具有必要的方法
            pages = Pdf.get_pages(path)
            print(pages)
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

    def trans_score(scores:list):
        T_scores = {
        "Molecule Name": [],
        "Molecular Weight": [],
        "XLogP": [],    
        "PSA": [],
        "Heavy Atoms": [],
        "Acceptor Count": [],
        "Donor Count": [],
        "Chelator Count": [],
        "Total Score": [],
        "Shape": [],
        "Hydrogen Bond": [],
        "Protein Desolvation": [],
        "Ligand Desolvation": []
        }
        for score in scores:
            T_scores["Molecule Name"].append(score["Molecule Name"])
            T_scores["Molecular Weight"].append(score["Molecular Weight"])
            T_scores["XLogP"].append(score["XLogP"])
            T_scores["PSA"].append(score["PSA"])
            T_scores["Heavy Atoms"].append(score["Heavy Atoms"])
            T_scores["Acceptor Count"].append(score["Acceptor Count"])
            T_scores["Donor Count"].append(score["Donor Count"])
            T_scores["Chelator Count"].append(score["Chelator Count"])
            T_scores["Total Score"].append(score["Total Score"])
            T_scores["Shape"].append(score["Shape"])
            T_scores["Hydrogen Bond"].append(score["Hydrogen Bond"])
            T_scores["Protein Desolvation"].append(score["Protein Desolvation"])
            T_scores["Ligand Desolvation"].append(score["Ligand Desolvation"])
        
        return T_scores

    def pyrun():
        from openeye import oechem, oedocking

        # 读取受体分子
        protein = oechem.OEGraphMol()
        oechem.OEReadMolecule(oechem.oemolistream("protein.pdb"), protein)

        # 读取配体分子
        ligand = oechem.OEGraphMol()
        oechem.OEReadMolecule(oechem.oemolistream("ligand.sdf"), ligand)

        # 创建一个OEDocking对象
        dock = oedocking.OEDock()

        # 初始化受体
        dock.Initialize(protein)

        # 对配体进行对接，并获取对接结果
        score = dock.DockMultiConformerMolecule(ligand)

        # 保存对接结果
        oechem.OEWriteMolecule(oechem.oemolostream("docked_ligand.sdf"), ligand)


         
path='D:\\BOrganic\\Documents\\Tencent Files\\1028209087\\FileRecv\\docking_report.pdf'
print(Openeye.get_score(path))