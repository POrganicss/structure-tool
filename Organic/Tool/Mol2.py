from rdkit import Chem

def generate_atom_combinations(smiles, target_atom):
    """
    生成给定分子和原子的所有可能组合。

    参数：
    - smiles (str): 分子的 SMILES 表示。
    - target_atom (str): 目标原子的符号。

    返回：
    - list: 包含所有可能组合的 RDKit Mol 对象的列表。
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        all_combinations = []

        # 1. 与目标原子以不同键连接，考虑键的性质和立体异构体
        for bond_type in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            for bond_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                new_mol = add_atom_with_bond(mol, target_atom, bond_type, bond_stereo)
                all_combinations.append(new_mol)

        # 2. 插入到分子中不同键的位置，考虑键的性质和立体异构体
        for bond_type in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            for bond_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                new_mol = insert_atom_into_bond(mol, target_atom, bond_type, bond_stereo)
                all_combinations.append(new_mol)

        # 3. 与分子中不同原子形成键连接，考虑键的性质和立体异构体
        for bond_type in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            for bond_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                new_mol_list = connect_atom_to_other_atoms(mol, target_atom, bond_type, bond_stereo)
                all_combinations.extend(new_mol_list)

        # 4. 考虑芳香键的情况
        for bond_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
            new_mol = add_aromatic_atom(mol, target_atom, bond_stereo)
            all_combinations.append(new_mol)

        # 5. 插入到环中的不同位置
        for bond_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
            new_mol = insert_atom_into_ring(mol, target_atom, bond_stereo)
            all_combinations.append(new_mol)

        # 6. 与环外的原子形成键连接
        for bond_stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
            new_mol_list = connect_atom_to_external_atoms(mol, target_atom, bond_stereo)
            all_combinations.extend(new_mol_list)

        return all_combinations

    else:
        print("无法生成分子。")

# 其他辅助函数的定义...

def add_atom_with_bondss(mol, target_atom, bond_type, bond_stereo=None):
    """
    将目标原子添加到分子中，并与目标原子以指定键连接。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_type (Chem.BondType): 连接键的类型。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - RDKit Mol 对象: 包含新原子的分子。
    """
    new_mol = Chem.RWMol(mol)
    target_atom_index = new_mol.GetNumAtoms()

    # 添加新原子
    new_atom = Chem.Atom(target_atom)
    new_atom_idx = new_mol.AddAtom(new_atom)

    # 添加新键
    new_mol.AddBond(target_atom_index, new_atom_idx, bond_type)

    if bond_stereo:
        new_mol.GetBondBetweenAtoms(target_atom_index, new_atom_idx).SetStereo(bond_stereo)

    # 更新分子，重建其结构
    new_mol = new_mol.GetMol()
    Chem.SanitizeMol(new_mol)

    return new_mol



def add_atom_with_bond(mol, target_atom_index, target_atom, bond_type, bond_stereo=None):
    """
    将目标原子添加到分子中，并与目标原子以指定键连接。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom_index (int): 与新原子形成键的原子的索引。
    - target_atom (str): 目标原子的符号。
    - bond_type (Chem.BondType): 连接键的类型。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - RDKit Mol 对象: 包含新原子的分子。
    """
    new_mol = Chem.RWMol(mol)

    # 添加新原子
    new_atom = Chem.Atom(target_atom)
    new_atom_idx = new_mol.AddAtom(new_atom)

    # 确保不会尝试添加自连接
    if target_atom_index == new_atom_idx:
        raise ValueError("Cannot add a bond between an atom and itself.")

    # 添加新键
    new_mol.AddBond(target_atom_index, new_atom_idx, bond_type)

    if bond_stereo:
        new_mol.GetBondBetweenAtoms(target_atom_index, new_atom_idx).SetStereo(bond_stereo)

    # 更新分子，重建其结构
    new_mol = new_mol.GetMol()
    Chem.SanitizeMol(new_mol)

    return new_mol


def add_atom_with_bonds(mol, target_atom, bond_type, bond_stereo):
    """
    将目标原子添加到分子中，并与目标原子以指定键连接。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_type (Chem.BondType): 连接键的类型。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - RDKit Mol 对象: 包含新原子的分子。
    """
    new_mol = Chem.AddHs(mol)
    target_atom_index = new_mol.GetSubstructMatch(Chem.MolFromSmiles(target_atom))[0]

    new_atom = Chem.Atom(target_atom)
    new_mol.AddAtom(new_atom)

    new_bond = Chem.BondTypeToBondType(bond_type)
    new_mol.AddBond(target_atom_index, mol.GetNumAtoms(), new_bond)

    new_mol.GetBondBetweenAtoms(target_atom_index, mol.GetNumAtoms()).SetStereo(bond_stereo)

    return new_mol

def insert_atom_into_bond(mol, target_atom, bond_type, bond_stereo):
    """
    将目标原子插入到分子中不同键的位置。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_type (Chem.BondType): 连接键的类型。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - RDKit Mol 对象: 包含插入原子的分子。
    """
    new_mol = Chem.AddHs(mol)

    for i in range(new_mol.GetNumBonds()):
        new_atom = Chem.Atom(target_atom)
        new_mol.AddAtom(new_atom)

        new_bond = Chem.BondTypeToBondType(bond_type)
        new_mol.AddBond(new_mol.GetNumAtoms() - 1, i, new_bond)

        new_mol.GetBondBetweenAtoms(new_mol.GetNumAtoms() - 1, i).SetStereo(bond_stereo)

    return new_mol

def connect_atom_to_other_atoms(mol, target_atom, bond_type, bond_stereo):
    """
    将目标原子与分子中不同原子以指定键连接。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_type (Chem.BondType): 连接键的类型。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - list: 包含所有连接组合的 RDKit Mol 对象的列表。
    """
    new_mol_list = []
    target_atom_index = mol.GetSubstructMatch(Chem.MolFromSmiles(target_atom))[0]

    for i in range(mol.GetNumAtoms()):
        if i != target_atom_index:
            new_mol = Chem.AddHs(mol)
            new_bond = Chem.BondTypeToBondType(bond_type)
            new_mol.AddBond(target_atom_index, i, new_bond)
            new_mol.GetBondBetweenAtoms(target_atom_index, i).SetStereo(bond_stereo)
            new_mol_list.append(new_mol)

    return new_mol_list

def add_aromatic_atom(mol, target_atom, bond_stereo):
    """
    将目标原子添加到分子中，并与目标原子以芳香键连接。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - RDKit Mol 对象: 包含新原子的分子。
    """
    new_mol = Chem.AddHs(mol)
    target_atom_index = new_mol.GetSubstructMatch(Chem.MolFromSmiles(target_atom))[0]

    new_atom = Chem.Atom(target_atom)
    new_mol.AddAtom(new_atom)

    new_bond = Chem.BondType.AROMATIC
    new_mol.AddBond(target_atom_index, mol.GetNumAtoms(), new_bond)

    new_mol.GetBondBetweenAtoms(target_atom_index, mol.GetNumAtoms()).SetStereo(bond_stereo)

    return new_mol

def insert_atom_into_ring(mol, target_atom, bond_stereo):
    """
    将目标原子插入到分子中不同环的位置。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - RDKit Mol 对象: 包含插入原子的分子。
    """
    new_mol = Chem.AddHs(mol)

    for ring in mol.GetRingInfo().BondRings():
        for i in range(len(ring)):
            new_atom = Chem.Atom(target_atom)
            new_mol.AddAtom(new_atom)

            new_bond = Chem.BondType.SINGLE  # 在环中插入原子使用单键
            new_mol.AddBond(new_mol.GetNumAtoms() - 1, ring[i], new_bond)

            new_mol.GetBondBetweenAtoms(new_mol.GetNumAtoms() - 1, ring[i]).SetStereo(bond_stereo)

    return new_mol

def connect_atom_to_external_atoms(mol, target_atom, bond_stereo):
    """
    将目标原子与分子中不同原子以指定键连接。

    参数：
    - mol (RDKit Mol 对象): 输入分子。
    - target_atom (str): 目标原子的符号。
    - bond_stereo (Chem.BondStereo): 连接键的立体异构体。

    返回：
    - list: 包含所有连接组合的 RDKit Mol 对象的列表。
    """
    new_mol_list = []
    target_atom_index = mol.GetSubstructMatch(Chem.MolFromSmiles(target_atom))[0]

    for i in range(mol.GetNumAtoms()):
        if i not in mol.GetRingInfo().AtomRings():  # 不连接到环内的原子
            new_mol = Chem.AddHs(mol)
            new_bond = Chem.BondType.SINGLE
            new_mol.AddBond(target_atom_index, i, new_bond)
            new_mol.GetBondBetweenAtoms(target_atom_index, i).SetStereo(bond_stereo)
            new_mol_list.append(new_mol)

    return new_mol_list





# 示例用法：
molecular_formula_example = "CCO"
target_atom_example = "O"
result_molecules = generate_atom_combinations(molecular_formula_example, target_atom_example)

# 打印每个生成的分子的 SMILES 字符串
for i, mol in enumerate(result_molecules):
    smiles = Chem.MolToSmiles(mol)
    print(f"Molecule {i+1}: {smiles}")
