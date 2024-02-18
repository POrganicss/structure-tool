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

    def Command(content):
        # 使用Popen创建进程，并与进程进行复杂的交互
        proc = subprocess.Popen(
            content,  # cmd特定的查询空间的命令
            stdin=None,  # 标准输入 键盘
            stdout=subprocess.PIPE,  # -1 标准输出（演示器、终端) 保存到管道中以便进行操作
            stderr=subprocess.PIPE,  # 标准错误，保存到管道
            shell=True,
        ).wait()
        # outinfo, errinfo = proc.communicate()   # 获取输出和错误信息
        # return outinfo.decode('gbk'), errinfo.decode('gbk')
        return proc

    def Filetransmission(remove_path, local_path):
        SSH = paramiko.Transport(("192.168.3.98", 22))
        SSH.connect(username="root", password="shuoxing614")
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
        return outinfo.decode("gbk")

    @staticmethod
    def exists(path):
        return os.path.exists(os.path.normpath(path))
