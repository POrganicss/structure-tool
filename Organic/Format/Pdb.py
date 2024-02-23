import __init__

class Pdb:
    def get_ligand_residue(pdb):
        lines = pdb.splitlines()
        residues = {}  # 创建一个字典来存储每个残基名称及其对应的原子数
        for line in lines:
            if line.startswith("HETATM"):
                columns = line.split()
                residue_name = columns[3]
                if residue_name not in residues:
                    residues[residue_name] = 0
                residues[residue_name] += 1
        # 找到原子数最多的残基名称
        max_residue_name = max(residues, key=residues.get)
        return max_residue_name