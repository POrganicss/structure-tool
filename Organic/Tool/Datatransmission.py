import paramiko
import subprocess
from gromacs_py import *

#通过SSH进行数据和文件的传输
class Datatransmission:
    def Commandtransmission(content):
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(hostname='192.168.3.98', port=22, username='root', password='shuoxing614',
                    look_for_keys=True)
        _, stdout, stderr = ssh.exec_command(content)
        out=stdout.read().decode()
        err=stderr.read().decode()
        ssh.close()

    def Command(content):
                                                # 使用Popen创建进程，并与进程进行复杂的交互
        proc = subprocess.Popen(content,        # cmd特定的查询空间的命令
            stdin=None,                         # 标准输入 键盘
            stdout=subprocess.PIPE,             # -1 标准输出（演示器、终端) 保存到管道中以便进行操作
            stderr=subprocess.PIPE,             # 标准错误，保存到管道
            shell=True).wait()
        #outinfo, errinfo = proc.communicate()   # 获取输出和错误信息
        #return outinfo.decode('gbk'), errinfo.decode('gbk')
        return proc

    def Filetransmission(remove_path, local_path):
        SSH = paramiko.Transport(('192.168.3.98', 22))
        SSH.connect(username='root', password='shuoxing614')
        sftp = paramiko.SFTPClient.from_transport(SSH)
        
        sftp.get(remove_path, local_path)
        sftp.close()
        SSH.close()
