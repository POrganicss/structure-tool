import os
import __init__

from Format.Gjf import *
from Smile import *
from Xyz import *
from Log import *
from Datatransmission import *
from Code import *
from File import *
from multiprocessing.pool import ThreadPool
import time
from alive_progress import alive_bar
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor, wait


class Gaussian:  # 提交gaussian任务的包，其中不涉及具体的任务类型
    path = "/root/chen/Compute"
    sxyzs = []

    # 自动提交任务并监控任务的施行，并生成坐标数据
    def run(gjf, proname, name, index=0):
        # ----------运算指令的编写----------
        code = []
        code.append("cd " + Gaussian.path + "\n")
        if Datatransmission.exist(Gaussian.path + "/" + proname) == False:
            code.append("mkdir " + proname + "\n")
        code.append("cd " + proname + "\n")

        while Datatransmission.exist(Gaussian.path + "/" + proname + "/" + name):
            match = re.match(r"(\D+)(\d+)$", name)
            if match:
                prefix = match.group(1)
                number = match.group(2)
            else:
                print("文件夹后面没有数字")
            name = prefix + str(int(number) + 1)

        code.append("mkdir " + name + "\n")
        code.append("cd " + name + "\n")

        code.append("touch " + name + ".gjf" + "\n")

        code.append('echo -e "' + gjf + '" > ' + name + ".gjf " + "\n")
        code.append("nohup time g16 " + name + ".gjf" + "&" + "\n")

        Datatransmission.Commandtransmission("".join(code))

        # ----------运算情况的判断----------
        if index == 0:  # 默认情况，只获取xyz文件
            return Gaussian.getxyz(proname, name)

        elif (
            index == 1
        ):  # 监控计算的结果，但成功时。输出True和能量和坐标，当计算失败时，输出False和坐标
            if Gaussian.getsuccess(proname, name):
                return (
                    True,
                    Gaussian.getxyz(proname, name),
                    Gaussian.getenergy(proname, name),
                )
            else:
                return False, Gaussian.getxyz(proname, name)

        elif index == 3:  # 监控计算时的情况
            while True:
                # 获取计算状态
                status = get_status(hostname, username, password, status_command)

                # 分析计算状态
                correct = is_correct(status)

                # 如果计算不正确，停止计算
                if not correct:
                    stop_job(hostname, username, password, stop_command)
                    break

                # 等待一段时间后再次检查状态
                time.sleep(60)

    def getxyz(proname, name):
        # ----------调用服务器端的Log_SSH.py运行----------
        code = []
        chkpath = Gaussian.path + "/" + proname + "/" + name + "/" + name + ".chk"
        code.append("formchk " + chkpath + "\n")

        code.append("python3 /root/Log_xyz.py ")
        code.append(
            Gaussian.path + "/" + proname + "/" + name + "/" + name + ".log" + "\n"
        )
        Datatransmission.Commandtransmission("".join(code))

        # ----------将服务器运算完的数据导入到本地----------
        path1 = Gaussian.path + "/" + proname + "/" + name + "/" + name + ".xyz"
        path2 = os.getcwd() + "/Temp/" + name + ".xyz"
        Datatransmission.Filetransmission(path1, path2)
        xyz = File.read(path2)
        # os.remove(path2)
        return xyz

    def getenergy(proname, name):
        # ----------调用服务器端的Log_energy.py运行----------
        code = []
        code.append("python3 /root/Log_engergy.py ")
        code.append(Gaussian.path + "/" + proname + "/" + name + "/" + name + ".log")
        Datatransmission.Commandtransmission("".join(code))

        # ----------将服务器运算完的数据导入到本地----------
        path1 = Gaussian.path + "/" + proname + "/" + name + "/" + name + ".eng"
        path2 = path + "/Temp/" + name + ".eng"
        Datatransmission.Filetransmission(path1, path2)
        energy = File.read(path2)

        # 输出Gibbs free energy
        G = energy.splitlines()[0].split(":")[1]
        print(name + " G:" + G)
        # os.remove(path2)
        return energy, energy.splitlines()[6].split(":")[1]

    def getsuccess(proname, name):
        # ----------调用服务器端的Log_success.py运行----------
        code = []
        code.append("python3 /root/Log_success.py ")
        code.append(
            Gaussian.path + "/" + proname + "/" + name + "/" + name + ".log" + "\n"
        )
        Datatransmission.Commandtransmission("".join(code))

        # ----------将服务器运算完的数据导入到本地----------
        path1 = Gaussian.path + "/" + proname + "/" + name + "/" + name + ".res"
        path2 = path + "/Temp/" + name + ".result"

        path3 = "C:\\Users\\10282\\gitee\\structure-tool\\Temp\\BozPhos.result"
        Datatransmission.Filetransmission(path1, path2)
        xyz = File.read(path2)
        os.remove(path2)

        # ----------将服务器运算完的数据导入到本地----------
        scode = []
        scode.append(
            "rm -f " + Gaussian.path + proname + "/" + name + "/" + name + ".res"
        )
        Datatransmission.Commandtransmission("".join(code))
        if "Error termination" in xyz:
            return False
        elif "Normal termination" in xyz:
            return True
        else:
            return False  # 意外错误

    # 将smile转化为gjf文件并提交到远程服务器，然后监控任务的施行，结束时导出需要的数据
    def runsmile(smile, proname, paramenters):
        try:
            xyz = Smile.toxyz(smile)
        except:
            print("错误分子：" + smile)
            # messagebox.showinfo('错误提示', '分子式错误')
        else:
            gjf = Xyz.togjf(xyz, paramenters)
            xyz = Gaussian.run(gjf, proname, paramenters["name"])
            return xyz

    def runsmiles(smiles, proname, paramenters):
        name = paramenters["name"]
        xyzs = []
        with alive_bar(len(smiles)) as bar:
            for i, smile in enumerate(smiles, start=1):  # 遍历任务
                paramenters["name"] = Code.getname([name, str(i)])
                xyz = Gaussian.runsmile(smile, proname, paramenters)
                xyzs.append(xyz)
                time.sleep(1)  # 假设这代码部分需要1s
                bar()  # 显示进度

        return xyzs

    def getitems(gjfs: list, proname: str, indexs: list):
        items = []
        for gjf, index in zip(gjfs, indexs):
            items.append([gjf, proname, index])
        return items

    def getitemss(gjfs: list, proname: str, paramenters: dict):
        items = []
        for gjf in gjfs:
            items.append([gjf, proname, paramenters])
        return items

    # 封装多线程调用
    def runmultsmile(item):
        xyz = Gaussian.run(item[0], item[1], item[2])
        return xyz

    # 封装多线程参数
    """ def getitems(smiles:list,proname:str,paramenters:dict):
        jobnumber=len(smiles)
        memsets=Code.memsets(jobnumber)
        nprocsets=Code.nprocsets(jobnumber)
        items=[]  
        for i,smile in enumerate(smiles,start=1):  # 遍历任务
            paramenters['name']=Code.getname(['proj',str(i)])
            paramenters['nproc']=nprocsets[i-1]
            paramenters['mem']=memsets[i-1]
            xyz=Smile.toxyz(smile)
            gjf=Xyz.togjf(xyz,paramenters)
            items.append([gjf,proname,paramenters['name']])
            
        return items
     """

    def runmulgjf(items):
        return Gaussian.run(items[0], items[1], items[2], 1)

    # 将多个smile转化为gjf文件并提交到远程服务器，然后监控任务的施行，结束时导出需要的数据
    def runmultsmiles(smiles, proname, paramenters):
        # paramenters={'cpu':'0-29','mem':'50GB','charge':'0','spin':'1',
        #'code':'# opt pm6',
        #'name':'proj'}
        jobnumber = len(smiles)
        memsets = Code.memsets(jobnumber)
        nprocsets = Code.nprocsets(jobnumber)
        items = []
        for i, smile in enumerate(smiles, start=1):  # 遍历任务
            paramenters["name"] = Code.getname(["proj", str(i)])
            paramenters["nproc"] = nprocsets[i - 1]
            paramenters["mem"] = memsets[i - 1]
            xyz = Smile.toxyz(smile)
            gjf = Xyz.togjf(xyz, paramenters)
            items.append([gjf, proname, paramenters["name"]])
            print(smile)

        print("------------gjf生成完成，开始计算------------")

        if jobnumber > 128:
            pool = ThreadPool(128)
            pool.map(Gaussian.runmultsmile, items)
            pool.close()  # 关闭线程池，不再接受线程
            pool.join()  # 等待线程池中线程全部执行完
        elif jobnumber < 128:
            pool = ThreadPool(jobnumber)
            pool.map(Gaussian.runmultsmile, items)
            pool.close()  # 关闭线程池，不再接受线程
            pool.join()  # 等待线程池中线程全部执行完

        # with alive_bar(jobnumber) as bar:
        #     for i,smile in enumerate(smiles,start=1):  # 遍历任务
        #         paramenters['name']=Gaussian.getname(paramenters['name'],str(i))
        #         paramenters['nproc']=nprocsets[i-1]
        #         paramenters['mem']=memsets[i-1]

        #         if jobnumber>128:
        #             pool = ThreadPool(128)
        #             items=[smiles,proname,paramenters]

        #             pool.map(Gaussian.runsmiles,items)
        #             pool.close()  # 关闭线程池，不再接受线程
        #             pool.join()  # 等待线程池中线程全部执行完

        #         if jobnumber<128:
        #             pool = ThreadPool(jobnumber)
        #             items=[smiles,proname,paramenters]
        #             pool.map(Gaussian.runsmiles,items)
        #             pool.close()  # 关闭线程池，不再接受线程
        #             pool.join()  # 等待线程池中线程全部执行完

        #         xyz=Gaussian.runsmile(smile,proname,paramenters)
        #         xyzs.append(xyz)

        #         time.sleep(1) # 假设这代码部分需要1s
        #         bar()  # 显示进度

        # # xyzs=[]
        # # for i,smile in enumerate(smiles,start=1):
        # #     xyz=Gaussian.runsmile(smile,[name,str(i)],proname)
        # #     xyzs.append(xyz)

    def runmultsmiles2(smiles, proname, paramenters):
        with futures.ThreadPoolExecutor(128) as executor:
            res = executor.map(Gaussian.runsmiles, smiles)

        return res

    # 利用多线程技术加速运算速度
    def multrunsmiles(smiles, proname, name="proj"):
        jobnumber = len(smiles)
        if jobnumber > 128:
            pool = ThreadPool(128)
            pool.map(Gaussian.runsmiles, smiles)

            pool.close()  # 关闭线程池，不再接受线程
            pool.join()  # 等待线程池中线程全部执行完

        if jobnumber < 128:
            pool = ThreadPool(jobnumber)
            paramenters = {
                "cpu": "0-127",
                "mem": "8GB",
                "charge": "0",
                "spin": "1",
                "code": "# opt freq B3LYP/6-31g em=gd3bj",
            }
            items = [smiles, proname, paramenters]
            pool.map(Gaussian.runsmiles, items)
            pool.close()  # 关闭线程池，不再接受线程
            pool.join()  # 等待线程池中线程全部执行完

    def multrunsmiless(items):

        thread_num = min(127, len(items[0]))

        executor = ThreadPoolExecutor(max_workers=128)  # 实例化线程池，thread_num个线程
        fs = []  # future列表
        xyzs = []

        for item in items:
            fs.append(executor.submit(Gaussian.runmultsmile, item))  # 提交任务

        wait(fs)  # 等待计算结束
        executor.shutdown()  # 销毁线程池

        return xyzs
