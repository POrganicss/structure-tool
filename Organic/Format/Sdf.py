import shlex
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
        
        #return File.read(os.path.join(path, Ligand_name + "_fixpka_oeomega.oeb.gz"))

    def getoeomega(path, Ligand_name, Mode='pose',nproc=1):
        # 使用shlex.quote安全地引用路径和文件名
        safe_path = shlex.quote(path)
        safe_ligand_name = shlex.quote(Ligand_name)
        input_file = os.path.join(safe_path, safe_ligand_name + ".sdf")
        output_file = os.path.join(safe_path, safe_ligand_name + "_fixpka.oeb.gz")
        omega_output_file = os.path.join(safe_path, safe_ligand_name + "_fixpka_oeomega.oeb.gz")
        log_file = os.path.join(safe_path, f"oeomega_{Mode}.log")

        # 构建命令，注意Mode参数应该如何正确传递取决于oeomega的实际参数规范
        commands = [
            f"fixpka -in {input_file} -out {output_file} && ",
            f"oeomega -mode {Mode} -mpi_np {nproc} -in {output_file} 
            -out {omega_output_file} -log {log_file} -progress log\n"
        ]
        # 执行命令
        LC.execute_command("".join(commands))
        
    def renames(ligands,names:list):
        modules=ligands.split("$$$$")
        new_modules=[]
        for i,module in enumerate(modules):
            lines=module.splitlines()
            if len(lines)>1:
                if "-MTS-"in lines[1]:
                    lines[0]=names[i]
                elif "-MTS-"in lines[2]:
                    lines[1]=names[i]
                new_modules.append("\n".join(lines))
        return "$$$$".join(new_modules)+"$$$$"+'\n'
    
     