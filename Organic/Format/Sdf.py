import __init__
import os
from Tool.Datatransmission import LocalCommand as LC
from Tool.File import File


class Sdf:
    
    ''' Available Modes:
    classic - The original customizable omega2 interface
    macrocycle - Conformer generation for macrocycles
    rocs - Optimal conformer generation for ROCS
    pose - Optimal conformer generation for molecular alignment
    and pose prediction by docking
    dense - Optimal conformer generation for FREEFORM '''
    
    def get_oeomega(path, Ligand_name, Mode='pose',nproc=1):
        command = []

        command.append(
            "fixpka "
            + " -in "
            + os.path.join(path, Ligand_name + ".sdf")
            + " -out "
            + os.path.join(path, Ligand_name + "_fixpka.oeb.gz")
            + " && "
        )
        command.append(
            "oeomega "
            + Mode
            + " -mpi_np "
            + str(nproc)
            + " -in "
            + os.path.join(path, Ligand_name + "_fixpka.oeb.gz")
            + " -out "
            + os.path.join(path, Ligand_name + "_fixpka_oeomega.oeb.gz")
            + " -log "
            + os.path.join(path, "oeomega_pose.log")
            + " -progress "
            + "log"
            + "\n"
        )
        LC.execute_command("".join(command))
        
        #return File.getdata(os.path.join(path, Ligand_name + "_fixpka_oeomega.oeb.gz"))
