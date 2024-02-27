import __init__
from Applications.Openeye import Openeye


path="D:\\BOrganic\\Desktop\\cs01\\Test24022301" # 项目路径
Ligands_name="CS-7-YA" #读取配体Ligands.sdf文件,后缀名sdf要小写
Protein_name="4agd_DeW" #读取蛋白质Protein.pdb文件,后缀名pdb要小写
nproc=6 #nproc数为并行计算的核心数
#配体文件和蛋白质文件都要在项目路径下面

Openeye.run(path, Ligands_name, Protein_name, nproc)
