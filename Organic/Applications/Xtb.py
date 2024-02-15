import __init__
from Tool.Datatransmission import LocalCommand


class Xtb:
    def run(self,xyz,inp,command):
        # 生成xtb命令
        xyz_name='xyz.xyz'
        inp_name='scan.inp'
        xtb_command=[]
        xtb_command.append('xtb '+xyz_name+' ')
        xtb_command.append('--input '+inp_name+' ')
        xtb_command.append('--opt ')
        xtb_command.append('--gfn2 ')
        xtb_command.append('--parallel '+str(4))
        
        # 执行xtb命令
        LocalCommand.execute_command(''.join(xtb_command))
        # 移动xtb文件
        LocalCommand.execute_command('mv ' + self.xtb_name + ' ' + self.path)
        return self.path + '/' + self.xtb_name