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
        File.create_files(path, "1.1Protein", "1.2Ligand", "1.3Complex", 
                          "2.1Solvate", "2.2Ions", "3.1Minimization",
                          "3.2NVTequilibration","3.3NPTequilibration",
                          "4.0MD","4.1Analysis")
        
        # -----------------------分别处理蛋白质与配体------------------------#
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
        
        # -----------------------溶剂化与离子化------------------------#
        
        command.append("gmx editconf"
                       +" -f "
                       + os.path.join(path,"2.1Solvate","complex.gro")
                       +" -o "
                       +os.path.join(path,"2.1Solvate","complex_newbox.gro")
                       +" -c -d 1.0 -bt dodecahedron ")
        command.append("gmx solvate"
                + " -cp "
                + os.path.join(path,"2.1Solvate","complex_newbox.gro")
                + " -cs "
                + "spc216.gro"
                + " -p "
                + os.path.join(path,"2.1Solvate","topol.top")
                + " -o "
                + os.path.join(path,"2.1Solvate","complex_solv.gro"))
        command.append("python cgenff_charmm2gmx.py ")
        LC.execute_command(command)
        File.cp(path,"2.1Solvate","2.2Ions",
                "solv.gro","topol.top","ions.itp","posre.itp")
        command.append("gmx grompp"
                       + " -f "
                       + os.path.join(path,"2.2Ions","ions.mdp")
                       + " -c "
                       + os.path.join(path,"2.2Ions","solv.gro")
                       + " -p "
                       + os.path.join(path,"2.2Ions","topol.top")
                       + " -o "
                       + os.path.join(path,"2.2Ions","ions.tpr")
                       + " -maxwarn 20")

        command.append("gmx genion"
                       + " -s "
                       + os.path.join(path,"2.2Ions","ions.tpr")
                       + " -o "
                       + os.path.join(path,"2.2Ions","solv_ions.gro")
                       + " -p "
                       + os.path.join(path,"2.2Ions","topol.top")
                       + " -pname NA -nname CL -neutral")
        outinfo=LC.execute_command(command)
        if 'a' in outinfo:
            outinfo=LC.execute_command('15')
        # -----------------------能量最小化与PVT,NVT平衡------------------------#
        File.cp(path,"2.2Ions","3.1minimization",
                "solv_ions.gro","topol.top","ions.itp","posre.itp")
        
        command.append("gmx grompp"
                       + " -f "
                       + os.path.join(path,"2.2Ions","em.mdp")
                       + " -c "
                       + os.path.join(path,"2.2Ions","solv_ions.gro")
                       + " -p "
                       + os.path.join(path,"2.2Ions","topol.top")
                       + " -o "
                       + os.path.join(path,"2.2Ions","em.tpr")
                       + " -maxwarn 20")
        command.append("gmx mdrun"
                       + " -v -deffnm em -ntmpi 1")
        
        File.cp(path,"3.1minimization","3.2nvtequilibration",
                "em.gro","topol.top","ions.itp","posre.itp")
        command.append("gmx make_ndx"
                        +"-f "
                        +os.path.join(path,"3.2nvtequilibration","em.gro")
                        +"-o "
                        +os.path.join(path,"3.2nvtequilibration","index_em.ndx"))
        outinfo=LC.execute_command(command)
        if 'a' in outinfo:
            outinfo=LC.execute_command('0 & ! a H* ')
            outinfo=LC.execute_command('q')
            
        command.append("gmx genrestr"
                       +" -f "
                       +os.path.join(path,"3.2nvtequilibration","em.gro")
                       +" -n "
                       +os.path.join(path,"3.2nvtequilibration","index_em.ndx")
                       +" -o "
                       +os.path.join(path,"3.2nvtequilibration","posre_em.itp")
                       +" -fc 1000 1000 1000")
        outinfo=LC.execute_command(command)
        if 'a' in outinfo:
            outinfo=LC.execute_command('3')
            
        command.append("gmx make_ndx"
                       +" -f "
                       +os.path.join(path,"3.2nvtequilibration","em.gro")
                       +" -o "
                       +os.path.join(path,"3.2nvtequilibration","index.ndx"))
        if 'a' in outinfo:
            outinfo=LC.execute_command('1 | 13')
            outinfo=LC.execute_command('q')
        command.append("gmx grompp"
                       +" -f "
                       +os.path.join(path,"3.2nvtequilibration","nvt.mdp")
                       +" -c "
                       +os.path.join(path,"3.2nvtequilibration","em.gro")
                       +" -r "
                       +os.path.join(path,"3.2nvtequilibration","em.gro")
                       +" -p "
                       +os.path.join(path,"3.2nvtequilibration","topol.top")
                       +" -n "
                       +os.path.join(path,"3.2nvtequilibration","index.ndx")
                       +" -o "
                       +os.path.join(path,"3.2nvtequilibration","nvt.tpr")
                       +" -maxwarn 20")
        command.append("gmx mdrun"
                       +" -v -deffnm nvt -ntmpi 1")
        
        File.cp(path,"3.2nvtequilibration","3.3nptequilibration",
                "em.gro","topol.top","ions.itp","posre.itp")
        command.append("gmx grompp"
                       +" -f "
                       +os.path.join(path,"3.3nptequilibration","npt.mdp")
                       +" -c "
                       +os.path.join(path,"3.3nptequilibration","nvt.gro")
                       +" -r "
                       +os.path.join(path,"3.3nptequilibration","nvt.gro")
                       +" -p "
                       +os.path.join(path,"3.3nptequilibration","topol.top")
                       +" -n "
                       +os.path.join(path,"3.3nptequilibration","index.ndx")
                       +" -o "
                       +os.path.join(path,"3.3nptequilibration","npt.tpr")
                       +" -maxwarn 21")
        command.append("gmx mdrun"
                       +" -v -deffnm npt -ntmpi 1")
        # -----------------------动力学模拟------------------------#
        File.cp(path,"3.2nvtequilibration","3.3nptequilibration",
                "em.gro","topol.top","ions.itp","posre.itp")
        command.append("gmx grompp"
                       +" -f "
                       +os.path.join(path,"4.0MD","md.mdp")
                       +" -c "
                       +os.path.join(path,"4.0MD","npt.gro")
                       +" -p "
                       +os.path.join(path,"4.0MD","topol.top")
                       +" -n "
                       +os.path.join(path,"4.0MD","index.ndx")
                       +" -o "
                       +os.path.join(path,"4.0MD","md_0_10.tpr")
                       +" -maxwarn 17")
        command.append("gmx mdrun"
                       +" -v -deffnm md_0_10 -ntmpi 1")
        