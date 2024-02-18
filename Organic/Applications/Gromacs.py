import os
import __init__
from Tool.Datatransmission import LocalCommand as LC
from Tool.File import File
class Gromacs:
    def run(proname,protein,ligand):
        path = os.path.join(path, proname)
        os.makedirs(path)

        File.tofile(path + "/Protein.pdb", protein)  # 将sdf文件写入文件
        File.tofile(path + "/ligand.pdb", ligand)  # 将pdb文件写入文件

        # 生成gromacs命令
        command = []
        # -----------------------新建项目文件夹-----------------------#
        command.append("mkdir " + os.path.join(path, "1.1Protein") + "\n")
        command.append("mkdir " + os.path.join(path, "1.2Ligand") + "\n")
        command.append("mkdir " + os.path.join(path, "1.3Complex") + "\n")
        
        command.append("mkdir " + os.path.join(path, "2.1Solvate") + "\n")
        command.append("mkdir " + os.path.join(path, "2.2Ions") + "\n")
        
        command.append("mkdir " + os.path.join(path, "3.1Minimization") + "\n")
        command.append("mkdir " + os.path.join(path, "3.2NVTequilibration") + "\n")
        command.append("mkdir " + os.path.join(path, "3.3NPTequilibration") + "\n")
        
        command.append("mkdir " + os.path.join(path, "4.0MD") + "\n")
        command.append("mkdir " + os.path.join(path, "4.1Analysis") + "\n")
        
        # 生成gromacs命令
        command.append(
            "cp "
            + os.path.join(path, "Protein.pdb")
            + " "
            + os.path.join(path, "1.1protein")
            + "\n"
        )
        
        command.append(
            "gmx pdb2gmx -f "
            + os.path.join(path, "1.1protein", "Protein.pdb")
            + " -ter -ff charmm36-jul2022.ff -water tip3p"
            + "\n"
        )
        
        outinfo=LC.execute_command(command)
        if 'a' in outinfo:
            outinfo=LC.execute_command('8')
            if 'b' in outinfo:
                outinfo=LC.execute_command('1')
                if 'c' in outinfo:
                    outinfo=LC.execute_command('0')
                    if 'd' in outinfo:
                        outinfo=LC.execute_command('0')
        
        
        command.append(
            "cp "
            + os.path.join(path, "ligand.pdb")
            + " "
            + os.path.join(path, "1.2Ligand")
            + "\n"
        )
        
        command.append(
            "perl "
            + os.path.join(path, "sort_mol2_bonds.pl")
            + " 7w.mol2 "
            + os.path.join(path, "1.2Ligand",'7w_fix.mol2')
            + "\n"
        )