import __init__
import os
import shlex
from Tool.Executor import Executor
class Openeye:
    #分子对接
    def run(path, Ligands_name, Protein_name, nproc):
        from Tool.Datatransmission import LocalCommand
        from Tool.File import File
        from Format.Sdf import Sdf
        from Format.Pdb import Pdb
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

    def get_score(path): # type: ignore
        from Format.Pdf import Pdf
        try:
            if os.path.isfile(path):
                pages = Pdf.get_pages(path)
        except Exception as e:
            print(f"访问PDF页面时出错: {e}")
            return []
        T_scores = {key: [] for key in ["Molecule Name", "Molecular Weight", "XLogP", "PSA", "Heavy Atoms",
                                        "Acceptor Count", "Donor Count", "Chelator Count", "Total Score",
                                        "Shape", "Hydrogen Bond", "Protein Desolvation", "Ligand Desolvation"]}
        def getpage(page):
            lines = page.splitlines()
            score = {}
            # 提取和存储分子属性
            for i, line in enumerate(lines):
                for key in T_scores.keys():
                    if key in line:
                        tf,score_value=Pdf.getnum(line)
                        if not tf:
                            score[key] = lines[i + 1].strip()
                            break  # 跳出内循环
                        else:
                            score[key] = score_value
                            break  # 跳出内循环
            return score
        scores=Executor.ThreadExecutor(getpage,pages)
        
        # 校验并存储分数数据
        required_keys = ['Shape', 'Hydrogen Bond', 'Protein Desolvation', 'Ligand Desolvation', 'Total Score']
        for score in scores:
            for key in T_scores.keys():
                T_scores[key].append(score[key])
            try:
                if all(key in score for key in required_keys):
                    # 计算总分并进行验证
                    _sum = sum(float(score[key]) for key in required_keys[:-1])  # 排除'Total Score'
                    if abs(_sum - float(score["Total Score"])) > 0.2:
                        # 总分不匹配，打印错误详情
                        print("总的误差：" + str(abs(_sum - float(score["Total Score"]))))
                        print(f"分子错误: {score.get('Molecule Name', '未知')}, 检测到分数不匹配。")
                        for key in required_keys:
                            print(f"{key} 数据: {score.get(key, 'N/A')}")
            except ValueError as e:
                # 处理分数转换中的错误
                print(f"分数转换过程中出现错误: {e}")

        return scores,T_scores
       
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
