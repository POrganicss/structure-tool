import os
import paramiko
import subprocess
from gromacs_py import *


# 通过SSH进行数据和文件的传输
class Datatransmission:
    hostname = "192.168.3.98"
    port = 22
    username = "root"
    password = "shuoxing614"

    def Commandtransmission(content):
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(
            Datatransmission.hostname,
            Datatransmission.port,
            Datatransmission.username,
            Datatransmission.password,
            look_for_keys=True,
        )
        _, stdout, stderr = ssh.exec_command(content)
        out = stdout.read().decode()
        err = stderr.read().decode()
        ssh.close()


    def Filetransmission(remove_path, local_path):
        SSH = paramiko.Transport((Datatransmission.hostname, Datatransmission.port))
        SSH.connect(Datatransmission.username, Datatransmission.password)
        sftp = paramiko.SFTPClient.from_transport(SSH)
        
        sftp.get(os.path.normpath(remove_path), os.path.normpath(local_path))
        sftp.close()
        SSH.close()

    def exist(path):

        # 创建 SSH 客户端
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # 连接到远程服务器
        ssh.connect(
            Datatransmission.hostname,
            Datatransmission.port,
            Datatransmission.username,
            Datatransmission.password,
        )

        # 创建 SFTP 客户端
        sftp = ssh.open_sftp()

        try:
            # 尝试获取远程文件或文件夹的属性
            sftp.stat(os.path.normpath(path))
            return True
        except IOError as _:
            # 如果 stat 调用失败，那么文件或文件夹可能不存在
            return False
        finally:
            # 关闭连接
            sftp.close()
            ssh.close()


class LocalCommand:
    @staticmethod
    def execute_command(command):
        # 使用Popen创建进程，并与进程进行复杂的交互
        proc = subprocess.Popen(
            command,  # cmd特定的查询空间的命令
            stdin=subprocess.PIPE,  # 标准输入 键盘
            stdout=subprocess.PIPE,  # -1 标准输出（演示器、终端) 保存到管道中以便进行操作
            stderr=subprocess.PIPE,  # 标准错误，保存到管道
            shell=True,
        )
        outinfo, errinfo = proc.communicate()  # 获取输出和错误信息
        #print(outinfo.decode("gbk"))
        #print(errinfo.decode("gbk"))
        return outinfo.decode("gbk")


    def submit(command):
        import wexpect
        """
        执行系统命令并允许与进程进行交互。
        
        :param initial_command: 要执行的初始命令字符串。
        :return: 一个wexpect的子程序对象，可用于后续与进程的交互。
        """
        if not command:  # 检查命令是否为空
            raise ValueError("命令不能为空")
        
        # 启动子进程
        child = wexpect.spawn(command, timeout=None)  # 设置timeout=None使其不会超时
        return child  # 返回子进程对象以供进一步交互
    
    
    @staticmethod
    def exists(path):
        return os.path.exists(path)
