import __init__
import os
import glob
from Tool.Datatransmission import LocalCommand as LC
from Tool.File import File
class Pdb:
    def get_receptor(path,Protein_name):
        if not os.path.exists(os.path.join(path, "Receptor")):
            os.makedirs(os.path.join(path, "Receptor"))
        Protein=File.getdata(os.path.join(path,Protein_name+'.pdb'))
        print(Pdb.get_ligand_residue(Protein))
        if Pdb.get_ligand_residue(Protein) is not None:
            print(Protein_name)
            LC.execute_command(
                "pdb2receptor "
                + " -pdb "
                + os.path.join(path,  Protein_name+'.pdb')
                + " -ligand_residue "
                + Pdb.get_ligand_residue(Protein)
                + " -receptor "
                + os.path.join(path, "Receptor",  Protein_name+"_receptor.oeb.gz")
                +'\n'
            )
    
    def get_receptors(path):
        Proteins = File.getfile_name(path,'pdb')
        for Protein in Proteins:
            Pdb.get_receptor(path,Protein)
            
    def get_ligand_residue(pdb):

        exclude_list=['HOH',"SO4","PO4","CL"]
        lines = pdb.splitlines()
        residues = {}  # 创建一个字典来存储每个残基名称及其对应的原子数
        for line in lines:
            if line.startswith("HETATM"):
                columns = line.split()
                residue_name = columns[3]
                if residue_name in exclude_list:
                    continue  # 如果残基名称在排除列表中，则跳过
                if residue_name not in residues:
                    residues[residue_name] = 0
                residues[residue_name] += 1
        # 如果residues为空，则表示没有找到HETATM头，返回None
        if not residues:
            return None
        # 找到原子数最多的残基名称
        print(residues)
        max_residue_name = max(residues, key=residues.get)
        return max_residue_name
    
    def protein_deal(Protein):
        #-----------清除额外数据-----------#
        residues = {}  # 创建一个字典来存储每个残基名称及其对应的原子数
        exclude_list=['HOH',"SO4","PO4","CL"]
        lines = Protein.splitlines()
        new_lines=[]
        for line in lines:
            if line.startswith("HETATM"): 
                columns = line.split()
                residue_name = columns[3]
                if residue_name not in exclude_list:
                    if residue_name not in residues:
                        residues[residue_name] = 0
                    residues[residue_name] += 1
                    new_lines.append(line)
            else:
                new_lines.append(line)
        max_residue_name = max(residues, key=residues.get)
        Protein = '\n'.join(new_lines)
        print(max_residue_name)
        #-----------读取重复单元-----------#
        Nosidues = {}  # 创建一个字典来存储每个残基名称及其对应的原子数
        lines = Protein.splitlines()
        new_lines=[]
        for line in lines:
            if line.startswith("HETATM"): 
                columns = line.split()
                residue_name = columns[3]
                No_name=columns[4]
                if residue_name == max_residue_name:
                    if No_name not in Nosidues:
                        Nosidues[No_name] = 0
                    Nosidues[No_name] += 1
        
        print(Nosidues)
        if Nosidues:
            max_value = max(Nosidues.values())
            max_Nosidues_name = next(key for key, value in Nosidues.items() if value == max_value)
        else:
            max_Nosidues_name = None
        print(max_Nosidues_name)
        #-----------删除重复单元-----------#
      
        for line in lines:
            if line.startswith("HETATM"): 
                columns = line.split()
                residue_name = columns[4]
                if residue_name != max_Nosidues_name:  # 修改此处的判断条件
                    new_lines.append(line)
            else:
                new_lines.append(line)
        # 将修改后的行重新组合成字符串
        new_protein = '\n'.join(new_lines)
        return new_protein
        
    def clean_protein(Protein):
        # 清除额外数据
        exclude_list = ['HOH', "SO4", "PO4", "CL"]
        residue_counts = {}  # 存储每个残基名称及其对应的原子数
        duplicate_counts = {}  # 存储每个残基名称及其对应的重复单元数
        
        lines = Protein.splitlines()
        cleaned_lines = []
        
        for line in lines:
            if line.startswith("HETATM"): 
                columns = line.split()
                residue_name = columns[3]
                if residue_name not in exclude_list:
                    # 更新残基计数
                    if residue_name not in residue_counts:
                        residue_counts[residue_name] = 0
                    residue_counts[residue_name] += 1
                    
                    # 更新重复单元计数
                    duplicate_name = columns[4]
                    if residue_name not in duplicate_counts:
                        duplicate_counts[residue_name] = {}
                    if duplicate_name not in duplicate_counts[residue_name]:
                        duplicate_counts[residue_name][duplicate_name] = 0
                    duplicate_counts[residue_name][duplicate_name] += 1
                    
                    cleaned_lines.append(line)
            else:
                cleaned_lines.append(line)
        
        # 获取最大残基名称和最大重复单元名称
        max_residue_name = max(residue_counts, key=residue_counts.get)
        max_duplicate_name = None
        if duplicate_counts:
            max_duplicate_value = max(max(sub_dict.values()) for sub_dict in duplicate_counts.values())
            max_duplicate_name = next((key for sub_dict in duplicate_counts.values() for key, value in sub_dict.items() if value == max_duplicate_value), None)
        
        # 删除重复单元
        cleaned_lines = [line for line in cleaned_lines if not (line.startswith("HETATM") and line.split()[4] == max_duplicate_name)]
        
        # 将修改后的行重新组合成字符串
        cleaned_protein = '\n'.join(cleaned_lines)
        
        return cleaned_protein,max_residue_name
    
path="D:\\BOrganic\\Desktop\\cs01\Protein\\W\\3g0f.pdb"
Protein=File.getdata(path)
new_protein=Pdb.clean_protein(Protein)

file_path="D:\\BOrganic\\Desktop\\cs01\Protein\\W\\3g0f_RR.pdb"
File.tofile(file_path,new_protein)