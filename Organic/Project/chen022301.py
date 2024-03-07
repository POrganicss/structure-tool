import __init__
from Applications.Openeye import Openeye
from Format.Sdf import Sdf
from Tool.File import File

#-------计算参数的设置-------#
path="D:\\BOrganic\\Desktop\\Protein_Ligand\\Openeye" # 项目路径
Ligands_name="Ligands" #读取配体Ligands.sdf文件,后缀名sdf要小写
Protein_name="3jyu" #读取蛋白质Protein.pdb文件,后缀名pdb要小写
nproc=6 #nproc数为并行计算的核心数
#配体文件和蛋白质文件都要在项目路径下面

#-------开始计算-------#
Openeye.run(path, Ligands_name, Protein_name, nproc)
