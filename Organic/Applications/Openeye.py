import os
import __init__
from Tool.Datatransmission import LocalCommand
from Tool.File import File


class Openeye:
    def run(proname, path, sdf, pdb, cpu):
        path = os.path.join(path, proname)
        os.makedirs(path)

        File.tofile(path + "/Temp.sdf", sdf)  # 将sdf文件写入文件
        File.tofile(path + "/Temp.pdb", pdb)  # 将pdb文件写入文件

        # 生成openeye命令
        command = []
        # -----------------------新建项目文件夹-----------------------#
        command.append("mkdir " + os.path.join(path, "fixpka") + "\n")
        command.append("mkdir " + os.path.join(path, "oeomega") + "\n")
        command.append("mkdir " + os.path.join(path, "hybrid") + "\n")
        command.append("mkdir " + os.path.join(path, "report") + "\n")

        # -----------------进行电荷平衡并生成分子3D构象----------------#
        command.append(
            "cp "
            + os.path.join(path, "Temp.sdf")
            + " "
            + os.path.join(path, "fixpka")
            + "\n"
        )
        command.append(
            "fixpka "
            + " -in "
            + os.path.join(path, "fixpka", "Temp.sdf")
            + " -out "
            + os.path.join(path, "fixpka", "Temp_fixpka.oeb.gz")
            + "\n"
        )

        command.append(
            "cp "
            + os.path.join(path, "fixpka", "Temp_fixpka.oeb.gz")
            + " "
            + os.path.join(path, "oeomega")
            + "\n"
        )
        command.append(
            "oeomega rocs "
            + " -in "
            + os.path.join(path, "oeomega", "Temp_fixpka.oeb.gz")
            + " -out "
            + os.path.join(path, "oeomega", "Temp_fixpka_oeomega.oeb.gz")
            + "\n"
        )
        LocalCommand.execute_command(command)
        print("3D构象生成完成！")
        command = []
        
        # ------------------生成对接口袋并进行分子对接-----------------#
        command.append(
            "cp "
            + os.path.join(path, "Temp.pdb")
            + " "
            + os.path.join(path, "hybrid")
            + "\n"
        )
        command.append(
            "make_receptor "
            + " -p "
            + os.path.join(path, "hybrid", "Temp.pdb")
            + " -o "
            + os.path.join(path, "hybrid", "receptor.oeb")
            + "\n"
        )

        command.append(
            "cp "
            + os.path.join(path, "oeomega", "Temp_fixpka_oeomega.oeb.gz")
            + " "
            + os.path.join(path, "hybrid")
            + "\n"
        )
        command.append(
            "hybrid "
            + " -receptor "
            + os.path.join(path, "hybrid", "receptor.oeb")
            + " -dbase "
            + os.path.join(path, "hybrid", "Temp_fixpka_oeomega.oeb.gz")
            + "\n"
        )
        
        LocalCommand.execute_command(command)
        print("对接口袋生成完成！")
        command = []
        
        # ------------------------生成对接结果------------------------#
        command.append(
            "cp "
            + os.path.join(path, "hybrid", "hybrid.oeb.gz")
            + " "
            + os.path.join(path, "report")
            + "\n"
        )
        command.append(
            "cp "
            + os.path.join(path, "hybrid", "receptor.oeb")
            + " "
            + os.path.join(path, "report")
            + "\n"
        )
        command.append(
            "docking_report"
            + " -mpi_np "
            + str(cpu)
            + " -docked_poses "
            + os.path.join(path, "hybrid", "hybrid.oeb.gz")
            + " -receptor "
            + os.path.join(path, "report", "receptor.oeb")
            + " -report "
            + os.path.join(path, "report", "report.pdf")
            + "\n"
        )
        LocalCommand.execute_command(command)
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
