import __init__
import os
from Tool.Datatransmission import LocalCommand as LC
from Tool.File import File
class Pdb:
    
    
    def get_receptor(path,Protein_name):
        #处理逻辑：1.去除蛋白质中的无机小分子
        #2.检查是否又两个相同的配体小分子
        #3.如果有，则进一步检查是否又重复的多肽链然后删除重复的多肽链
        
        Protein=File.read(os.path.join(path,Protein_name+'.pdb'))
        new_protein,max_Nosidues_name=Pdb.clean(Protein)
        
        
        #chain=Pdb.simplify_chain(os.path.join(path,Protein_name+'.pdb'))
        Protein_name=Protein_name+"_clean"
        new_protein,max_Nosidues_name=Pdb.clean(Protein)
        File.save(os.path.join(path,Protein_name+'.pdb'),new_protein)
        print(max_Nosidues_name)
        if Pdb.get_ligand_residue(Protein) is not None:
            LC.execute_command(
                "pdb2receptor "
                + " -pdb "
                + os.path.join(path,  Protein_name+'.pdb')
                + " -ligand_residue "
                + max_Nosidues_name
                + " -receptor "
                + os.path.join(path,  Protein_name+"_receptor.oeb.gz")
                +'\n'
            )
            
        #return File.read(os.path.join(path,Protein_name+"_receptor.oeb.gz"))
        
    def get_receptors(path,Protein_names:list):
        #Protein_names = File.getfile_name(path,'pdb')#获取该路径下所有pdb文件名字
        for Protein_name in Protein_names:#遍历所有pdb文件名
            Pdb.get_receptor(path,Protein_name)
            
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
        max_residue_name = max(residues, key=residues.get)
        return max_residue_name
    
    def simplify_chain(Protein_path):
        import random
        import os
        from Bio import PDB

        # 创建PDB解析器
        parser = PDB.PDBParser(QUIET=True)

        # 解析PDB文本
        structure = parser.get_structure("pdb_structure", Protein_path)

        # 遍历每一个链，记录是否包含小分子
        chains_with_ligands = []
        chains_without_ligands = []
        for model in structure:
            for chain in model:
                has_ligand = any(atom.element == "H" for atom in chain.get_atoms() if atom.id == "HETATM")
                if has_ligand:
                    chains_with_ligands.append(chain)
                else:
                    chains_without_ligands.append(chain)

        # 如果PDB中只有一条链，则不进行修改，直接返回原始PDB文本
        if len(chains_with_ligands) + len(chains_without_ligands) == 1:
            return File.read(Protein_path)

        # 优先选择包含小分子的链，如果没有，则随机选择一个链
        if chains_with_ligands:
            selected_chain = random.choice(chains_with_ligands)
        else:
            selected_chain = random.choice(chains_without_ligands)

        # 创建PDB输出对象
        io = PDB.PDBIO()

        # 设置输出文件的结构为选定的链
        io.set_structure(selected_chain)

        # 创建一个内存中的文件对象，用于保存PDB文本
        temp_out_file = io.save("temp_out.pdb")

        # 读取临时输出的PDB文件内容
        with open("temp_out.pdb", "r") as temp_out_file:
            new_protein = temp_out_file.read()
        
        # 移除临时文件
        os.remove("temp_out.pdb")
        
        return new_protein
    
    #蛋白质处理程序，1.去除
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

    def clean(Protein):
        residues={}
        exclude_list=['HOH',"SO4","PO4","CL"]
        lines = Protein.splitlines()
        new_lines=[]
        for line in lines:
            if line.startswith("HETATM"): 
                columns = line.split()
                residue_name=columns[3]
                if columns[3] not in exclude_list:
                    if residue_name not in residues:
                        residues[residue_name] = 0
                    residues[residue_name] += 1
                    new_lines.append(line)
            else:
                new_lines.append(line)
        if residues:
            max_Nosidues_name = next(iter(residues.keys()))  # 或者使用 next(iter(Nosidues.items()))[0]
        else:
            max_Nosidues_name = None
        return '\n'.join(new_lines),max_Nosidues_name


#gmx pdb2gmx -f 3HTB_clean.pdb -o 3HTB_processed.gro -ter

    def get_top(path,Protein_name):
        command=[]
        protein_path=os.path.join(path, Protein_name+'.pdb')
        out_path=os.path.join(path, Protein_name+'_processed.gro')
        
        command.append("gmx pdb2gmx"
                       +" -f "
                       + protein_path
                       + " -o "
                       + out_path+"_processed.gro"
                       + " -ter"
                       )
        
        #command1="cd D:\\BOrganic\\Desktop\\Protein_Ligand\\getop"
        command2="gmx pdb2gmx -f 3HTB_clean.pdb -o 3HTB_clean_processed.gro_processed.gro -ter"
        #print(command)
        #child=LC.submit(command1)
        #child.interact()
        #print(child.before)
        
        #child.sendline(command2)
        #print(child.before)
        
        os.chdir("D:\\BOrganic\\Desktop\\Protein_Ligand\\getop")

        child=LC.submit(command2)
        while True:  # 继续循环直到进程结束或你决定中断
            child.expect('\n', timeout=1)  # 假设每次期待新行的出现，设置超时为1秒
            print(child.before.decode('utf-8'))  # 打印自上次匹配以来的所有输出
        child.expect('\n')
        print(child.before)
        child.terminate()
path="D:\\BOrganic\\Desktop\\Protein_Ligand\\getop"
Pdb.get_top(path,"3HTB_clean")