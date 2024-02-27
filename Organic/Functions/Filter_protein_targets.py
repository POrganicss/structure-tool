import os
import __init__
from Format.Pdb import Pdb
from Format.Sdf import Sdf
from Tool.Datatransmission import LocalCommand
from Tool.File import File

class Filter_protein_targets:
    
    def Filter(path,):
        Pdb.get_receptors(os.path.join(path,Protein_name+'.pdb'))

