import os
import __init__
from Format.Pdb import Pdb
from Format.Sdf import Sdf
from Tool.Datatransmission import LocalCommand
from Tool.File import File

class Openeye:
    #分子对接
    def run(path, Ligands_name, Protein_name, nproc):

        # 生成openeye命令
        command = []
        # -----------------------新建项目文件夹-----------------------#
        File.create_files(path, "fixpka_oeomega", "hybrid", "report")
        File.copy(path, os.path.join(path, "fixpka_oeomega"), Ligands_name+".sdf")
        File.copy(path,os.path.join(path, "hybrid"), Protein_name+".pdb")
        
        
        # -----------------进行电荷平衡并生成分子3D构象----------------#
        Sdf.get_oeomega(os.path.join(path, "fixpka_oeomega"),Ligands_name,nproc=nproc)
        
        File.cp(path,"fixpka_oeomega","hybrid",Ligands_name+"_fixpka_oeomega.oeb.gz")
        
        # ------------------生成对接口袋并进行分子对接-----------------#
        Pdb.get_receptor(os.path.join(path, "hybrid"),Protein_name)
        
        command.append(
            "hybrid "
            + " -mpi_np "
            + str(nproc)
            + " -receptor "
            + os.path.join(path, "hybrid", Protein_name+"_clean_receptor.oeb.gz")
            + " -dbase "
            + os.path.join(path, "hybrid", Ligands_name+"_fixpka_oeomega.oeb.gz")
            + " -docked_molecule_file "
            + os.path.join(path, "hybrid", Ligands_name+"_"+Protein_name+"hybrid_docked_molecule.sdf")
            + " -undocked_molecule_file "
            + os.path.join(path, "hybrid", Ligands_name+"_"+Protein_name+"hybrid_undocked_molecule.sdf")
            + " -score_file "
            + os.path.join(path, "hybrid", "hybrid_score.txt")
            + " -report_file "
            + os.path.join(path, "hybrid", "hybrid_report.txt")
            + " -settings_file "
            + os.path.join(path, "hybrid", "hybrid_settings.param")
            + " -status_file "
            + os.path.join(path, "hybrid", "hybrid_status.txt")
            + "\n"
            
        )
        LocalCommand.execute_command(''.join(command))
        print("分子对接完成！")
        command = []
        
        # ------------------------生成对接结果------------------------#
        File.copy(os.path.join(path, "hybrid"), os.path.join(path, "report"), Ligands_name+"_"+Protein_name+"hybrid_docked_molecule.sdf", Protein_name+"_clean_receptor.oeb.gz")
        
        command.append(
            "docking_report"
            + " -docked_poses "
            + os.path.join(path, "hybrid", Ligands_name+"_"+Protein_name+"hybrid_docked_molecule.sdf")
            + " -receptor "
            + os.path.join(path, "report", Protein_name+"_clean_receptor.oeb.gz")
            + " -report "
            + os.path.join(path, "report", "report.pdf")
            + "\n"
        )
        LocalCommand.execute_command(''.join(command))
        print("对接结果生成完成！")
        command = []
        # -----------------------获取对接结果------------------------#

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


