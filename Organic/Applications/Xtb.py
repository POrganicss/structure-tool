import os
import tempfile
import __init__
from Tool.Datatransmission import LocalCommand
from Tool.File import File


class Xtb:
    def run(xyz, inp, command):
        # 创建临时文件
        with tempfile.NamedTemporaryFile(
            delete=True, suffix=".xyz"
        ) as xyz_file, tempfile.NamedTemporaryFile(
            delete=True, suffix=".inp"
        ) as inp_file:
            # 写入xyz和inp内容
            xyz_file.write(xyz.encode())
            inp_file.write(inp.encode())

            # 确保内容被写入磁盘
            xyz_file.flush()
            inp_file.flush()

            # 生成xtb命令
            xtb_command = []
            xtb_command.append("xtb " + xyz_file.name + " ")
            xtb_command.append("--input " + inp_file.name + " ")
            xtb_command.append("--parallel " + str(4) + " ")
            xtb_command.append(command)

            # 执行xtb命令
            LocalCommand.execute_command("".join(xtb_command))

        # 获取临时文件的文件名（不含路径和后缀）
        result_name = os.path.splitext(os.path.basename(xyz_file.name))[0] + ".log"
        # 读取xtb结果文件
        result = File.read(result_name + ".log发")

        # 删除结果文件
        os.remove(result_name)

        return result
