import os
import __init__
from Format.Pdb import Pdb
from Format.Sdf import Sdf
from Tool.Datatransmission import LocalCommand
from Tool.File import File
from Format.Pdf import Pdf

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
        
        #hybrid -mpi_np 4 -receptor 3jyu.oeb -dbase Ligands_fixpka_oeomega.oeb.gz -docked_molecule_file hybrid_docked_molecule.sdf
        
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
        
        # docking_report -docked_poses hybrid_docked_molecule.sdf -receptor 3jyu.oeb
        
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

    def get_score(path):
        pages=Pdf.get_pages(path)
        scores=[]
        for page in pages:
            lines=page.splitlines()
            score={}
            
            #--------------------读取基本参数--------------------#
            for i,line in enumerate(lines):
                if "Molecule Name" in line:#读取分子名称
                    score["Molecule Name"]=lines[i+1]
                elif "Molecular Weight" in line:#读取分子重量
                    score["Molecular Weight"]=lines[i+1]
                elif "XLogP" in line:
                    score["XLogP"]=lines[i+1]  
                elif "PSA" in line:
                    score["PSA"]=lines[i+1]    
                elif "Heavy Atoms" in line:
                    score["Heavy Atoms"]=lines[i+1]
                elif "Acceptor Count" in line:
                    score["Acceptor Count"]=lines[i+1] 
                elif "Donor Count" in line:
                    score["Donor Count"]=lines[i+1]  
                elif "Chelator Count" in line:
                    score["Chelator Count"]=lines[i+1] 
                    
                #--------------------读取评分--------------------#
            for line in lines:     
                if "Total Score" in line:
                    score["Total Score"]=Pdf.getnum(line)
                elif "Shape" in line:
                    score["Shape"]=Pdf.getnum(line)
                elif "Hydrogen Bond" in line:
                    score["Hydrogen Bond"]=Pdf.getnum(line)
                elif "Protein Desolvation" in line:
                    score["Protein Desolvation"]=Pdf.getnum(line)
                elif "Ligand Desolvation" in line:
                    score["Ligand Desolvation"]=Pdf.getnum(line)
            _sum=float(score["Shape"])+float(score["Hydrogen Bond"])+float(score["Protein Desolvation"])+float(score["Ligand Desolvation"])
            
            if _sum-float(score["Total Score"])<=0.01:
                scores.append(score)
            else:
                print("分子名称："+score["Molecule Name"]+" , 分数统计有误")
                print("Shape数据:" +score["Shape"])
                print("Hydrogen Bond数据:" +score["Hydrogen Bond"])
                print("Protein Desolvation数据:" +score["Protein Desolvation"])
                print("Ligand Desolvation数据:" +score["Ligand Desolvation"])
                print("Total Score数据:" +score["Total Score"])
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



