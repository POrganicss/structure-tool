
from matplotlib import pyplot as plt
import mplcursors
import numpy as np
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem
from File import File
#提供绘图工具
class Draw:
    
    def draw_PESC(data, ENERGY_LEVEL_LENGTH=8.1, LINE_LENGTH=12.0):
        
        xylines = [[None] * (3 * len(data[i])) for i in range(len(data))]

        for i in range(len(data)):
            for j in range(len(data[i])):
                xylines[i][3 * j] = ((ENERGY_LEVEL_LENGTH +
                                    LINE_LENGTH) * float(j), data[i][j])
                xylines[i][3 * j + 1] = ((ENERGY_LEVEL_LENGTH + LINE_LENGTH)
                                        * float(j) + ENERGY_LEVEL_LENGTH / 2.0, data[i][j])
                xylines[i][3 * j + 2] = (ENERGY_LEVEL_LENGTH +
                                        (ENERGY_LEVEL_LENGTH + LINE_LENGTH) * float(j), data[i][j])

        mol = ob.OBMol()
        
        for _ in range(len(xylines)*len(xylines[0])):
            atom = mol.NewAtom()
            atom.SetAtomicNum(6)  # Atomic number for carbon
            mol.AddAtom(atom)
        atom_index = 0
        for points in xylines:
            for point in points:
                atom = mol.GetAtom(atom_index + 1)
                atom.SetVector(point[0], point[1], 0)
                atom_index += 1
           
        for i in range(1, mol.NumAtoms()):
            mol.AddBond(i, i + 1, 1)  # 1 represents a single bond
        
        sum_i=0
        for i in range(1, len(xylines)+1):
            index = (len(xylines[i-1]) - 1)+sum_i
            sum_i = index
            mol.DeleteBond(mol.GetBond(index))

        return mol
    
    def draw_line(*args, labels={'xlabel': 'x', 'ylabel': 'y','zlabel': 'z', 'label': 'Draws', 'title': 'title','marker':'x'}, limit=''):
        for i,arg in enumerate(args):
            print("第"+str(i+1)+'个数据的个数: '+str(len(arg)))
        
        if len(args)==1:
            x = list(range(1, len(args[0]) + 1))
            plt.plot(x, args[0], marker=labels['marker'],
                     label=labels['label'])
            # 添加基准线
            if limit != '' and isinstance(limit, list) and round(len(limit)/2) == len(limit)/2:
                for i in range(int(len(limit)/2)):
                    if min(args[0]) - float(limit[i*2]) < limit[i*2]*limit[i*2+1]:
                        plt.axhline(
                            y=float(limit[i*2]), color='r', linestyle='--', label='y='+str(limit[i*2]))
            # 在x轴上以5为间隔设置刻度
            plt.xticks(np.array(x))
            # 设置标签和标题
            plt.xlabel(labels['xlabel'])
            plt.ylabel(labels['ylabel'])
            plt.title(labels['title'])
            # 在底部添加数据点的总数
            plt.annotate(f'Data Sizes: {len(args[0])}', xy=(
                0.5, -0.1), ha='center', va='center', xycoords='axes fraction')
            # 显示图例
            plt.legend()
            mplcursors.cursor(hover=True)
        elif len(args)==2:
            #plt.plot(args[0], args[1], marker='x', label=labels['label'])
            plt.scatter(args[0], args[1], marker='x', label=labels['label'])
            # 在x轴上以5为间隔设置刻度
            plt.xticks(np.array(args[0]))
            # 设置标签和标题
            plt.xlabel(labels['xlabel'])
            plt.ylabel(labels['ylabel'])
            plt.title(labels['title'])
            # 在底部添加数据点的总数
            plt.annotate(f'Data Sizes: {len(args[0])}', xy=(
                0.5, -0.1), ha='center', va='center', xycoords='axes fraction')
            # 显示图例
            plt.legend()
            mplcursors.cursor(hover=True)
        elif len(args)==3:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            #ax.plot(args[0], args[1], args[2])
            #ax.plot_surface(args[0], args[1], args[2], cmap='viridis')
            ax.scatter(args[0], args[1], args[2])
            
            ax.set_xlabel(labels['xlabel'])
            ax.set_ylabel(labels['ylabel'])
            ax.set_zlabel(labels['zlabel'])
            ax.set_title(labels['title'])
        plt.show()

    def draw_mol(smile):
        moldH = Chem.MolFromSmiles(smile)
        mol=Chem.AddHs(moldH)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        img = Chem.Draw.MolToImage(mol, size=(
            1500, 1500), kekulize=True, highlightBonds=[0, 1])
        img.show()

    def draw_dline():
        print()
        
a=Draw.draw_PESC([[0,22.02,-5.09,-8.35,-18.34,-37.14,-4.34,-15.82]])

File.save('test.cdxml',a)
