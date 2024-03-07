import os
import __init__
from Format.Pdb import Pdb
from Format.Sdf import Sdf
from Tool.Datatransmission import LocalCommand
from Tool.File import File

class Filter_protein_targets:
    
    def Filter(path):
        Protein_names = File.getfile_name(path,'pdb')#获取该路径下所有pdb文件名字
        Pdb.get_receptors(os.path.join(path,"Protein"),Protein_names)#获取所有蛋白质的受体
        

