import os
import paramiko
import subprocess
from gromacs_py import *
import wexpect

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

    def execute_command(command, input_queue, output_queue):
        import wexpect
        from queue import Queue, Empty
        import threading
        """
        执行命令并实时处理输入和输出。
        
        :param command: 要执行的命令字符串。
        :param input_queue: 一个队列，用于从外部发送输入到命令。
        :param output_queue: 一个队列，用于将命令的输出发送到外部。
        """
        child = wexpect.spawn(command, timeout=None)

        def send_input():
            while True:
                try:
                    # 等待外部输入，这里的超时时间可以根据实际情况进行调整
                    line = input_queue.get(timeout=0.1)
                    child.sendline(line)
                except Empty:
                    # 如果没有外部输入，继续等待
                    continue

        # 创建一个线程用于监听外部输入并发送到子进程
        input_thread = threading.Thread(target=send_input)
        input_thread.daemon = True  # 设为守护线程，确保主进程退出时子线程也会被结束
        input_thread.start()

        # 主循环，用于读取命令的输出
        while True:
            try:
                line = child.readline(timeout=0.1)  # 这里的超时时间可以根据实际情况进行调整
                if line:
                    # 将输出放到输出队列
                    output_queue.put(line.decode('utf-8'))
                if child.eof():
                    # 如果子进程结束，退出循环
                    break
            except wexpect.TIMEOUT:
                # 超时意味着暂时没有输出，继续等待
                continue

        # 等待输入线程结束
        input_thread.join()

    def submit():
        import subprocess
        import threading
        from queue import Queue, Empty
        import time

        # 初始化队列
        input_queue = Queue()
        output_queue = Queue()

        # 启动子进程
        proc = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        # 定义读取输出的线程函数
        def read_output():
            while True:
                # 读取一行输出
                output = proc.stdout.readline()
                if output:
                    output_queue.put(output)
                else:
                    # 如果没有输出，可能子进程已结束
                    break

        # 定义发送输入的线程函数
        def send_input():
            while True:
                try:
                    # 从外部获取输入
                    line = input_queue.get(timeout=0.1)
                    # 向子进程发送输入
                    proc.stdin.write(line + '\n')
                    proc.stdin.flush()
                except Empty:
                    # 如果没有输入，继续循环等待
                    continue

        # 创建和启动线程
        output_thread = threading.Thread(target=read_output, daemon=True)
        output_thread.start()

        input_thread = threading.Thread(target=send_input, daemon=True)
        input_thread.start()

        # 在这里，您可以模拟外部输入，例如：
        # input_queue.put('echo Hello World')

        # 等待线程结束（在实际情况下，您可能需要更复杂的逻辑来决定何时停止）
        input_thread.join()
        output_thread.join()

        # 在脚本结束时关闭子进程
        proc.terminate()

    @staticmethod
    def exists(path):
        return os.path.exists(path)
