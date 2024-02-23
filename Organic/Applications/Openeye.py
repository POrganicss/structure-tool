import os
import __init__
from Format.Pdb import Pdb
from Tool.Datatransmission import LocalCommand
from Tool.File import File


class Openeye:
    
    def run(proname, path, Ligands, Protein, nproc):
        path = os.path.join(path, proname)
        
        if not os.path.exists(path):
            os.makedirs(path)

        File.tofile(path + "\\Ligands.sdf", Ligands)  # 将sdf文件写入文件
        File.tofile(path + "\\Protein.pdb", Protein)  # 将pdb文件写入文件

        # 生成openeye命令
        command = []
        # -----------------------新建项目文件夹-----------------------#
        File.create_files(path,"fixpka", "oeomega", "hybrid", "report")
        File.copy(path, os.path.join(path, "fixpka"), "Ligands.sdf")
        File.copy(path,os.path.join(path, "hybrid"), "Protein.pdb")
        
        
        # -----------------进行电荷平衡并生成分子3D构象----------------#
        command.append(
            "fixpka "
            + " -in "
            + os.path.join(path, "fixpka", "Ligands.sdf")
            + " -out "
            + os.path.join(path, "fixpka", "Ligands_fixpka.oeb.gz")
            + " && "
        )
        
        command.append(
            "copy "
            + os.path.join(path, "fixpka", "Ligands_fixpka.oeb.gz")
            + " "
            + os.path.join(path, "oeomega")
            + " && "
        )
        
        command.append(
            "oeomega "
            + " pose "
            + " -mpi_np "
            + str(nproc)
            + " -in "
            + os.path.join(path, "oeomega", "Ligands_fixpka.oeb.gz")
            + " -out "
            + os.path.join(path, "oeomega", "Ligands_fixpka_oeomega.oeb.gz")
            + " -log "
            + os.path.join(path, "oeomega", "oeomega_pose.log")
            + " -progress "
            + "log"
            + "\n"
        )
        LocalCommand.execute_command(''.join(command))
        print("3D构象生成完成！")
        command = []
        
        # ------------------生成对接口袋并进行分子对接-----------------#
        command.append(
            "pdb2receptor "
            + " -pdb "
            + os.path.join(path, "hybrid", "Protein.pdb")
            + " -ligand_residue "
            + Pdb.get_ligand_residue(Protein)
            + " -receptor "
            + os.path.join(path, "hybrid", "receptor.oeb.gz")
            + " && "
        )
        
        command.append(
            "copy "
            + os.path.join(path, "oeomega", "Ligands_fixpka_oeomega.oeb.gz")
            + " "
            + os.path.join(path, "hybrid")
            + " && "
        )
        command.append(
            "hybrid "
            + " -mpi_np "
            + str(nproc)
            + " -receptor "
            + os.path.join(path, "hybrid", "receptor.oeb.gz")
            + " -dbase "
            + os.path.join(path, "hybrid", "Ligands_fixpka_oeomega.oeb.gz")
            + " -docked_molecule_file "
            + os.path.join(path, "hybrid", "hybrid_docked_molecule.sdf")
            + " -undocked_molecule_file "
            + os.path.join(path, "hybrid", "hybrid_undocked_molecul.sdf")
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
        File.copy(os.path.join(path, "hybrid"), os.path.join(path, "report"), "hybrid_docked_molecule.sdf", "receptor.oeb.gz")
        
        command.append(
            "docking_report"
            + " -docked_poses "
            + os.path.join(path, "hybrid", "hybrid_docked_molecule.sdf")
            + " -receptor "
            + os.path.join(path, "report", "receptor.oeb.gz")
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


def openeye_test():
    # Test24022301为项目名称
    path="D:\\BOrganic\\Desktop\\cs01" # 项目路径
    Ligands=File.getdata(path+"/CS-7-YA.sdf") #读取配体Ligands.sdf文件
    Protein=File.getdata(path+"/4agd_DeW.pdb") #读取蛋白质Protein.pdb文件
    print(type(Ligands))
    # 6, nproc数为并行计算的核心数

    Openeye.run("Test24022301", path, Ligands, Protein, nproc=6)
    
openeye_test()
