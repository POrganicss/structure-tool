import __init__

from Tool.Verify import *
from Tool.Code import *
from Tool.File import *
from Tool.Datatransmission import *

from Format.Log import *
from Format.Xyz import *
from Format.Gjf import *
from Format.Smile import *
from Tool.Executor import Executor
from Applications.Gaussian import *

from tqdm import tqdm
from Data import *
import concurrent.futures


class Mechanism_Calculation:

    def gjf_optimize(gjf: str, proname: str, index=1, solvent=""):
        Verify.Is([gjf, proname, index, solvent], [str, str, int, str])
        para = Gjf.getparameters(gjf)
        name = Code.getname(para["name"])  # 给文件加上日期后缀

        if index == 1:
            # 直接进行优化
            if "freq" not in para["additional"]:
                para["additional"] = para["additional"] + " freq"
            xyz = Gjf.toxyz(gjf)
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para["name"], 1)

            if result[0]:
                return result[1], result[2]
            else:
                print("计算失败")
                return result[1]

        elif index == 2:

            # 进行一第一次优化
            para["name"] = name + "_01"
            if "freq" in para["additional"]:
                para["additional"] = (
                    para["additional"].replace("freq", " ").replace("  ", " ")
                )
            xyz = Gaussian.run(gjf, proname, para["name"])

            # 第二次优化的参数设置
            para["opt"] = (
                para["opt"].replace("calcfc", "calcall") + ",tight,maxstep=1,notrust"
            )
            if "freq" not in para["additional"]:
                para["additional"] = para["additional"] + " freq"
            if "int=superfine" not in para["additional"]:
                para["additional"] = para["additional"] + " int=superfine"
            para["name"] = name + "_02"

            # 进行第二次几何优化
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para["name"], 1)

            if result[0]:
                return result[1], result[2]
            else:
                print("计算失败")
                return result[1]
        elif index == 3:
            # 进行一第一次优化
            para["name"] = name + "_01"
            if "freq" in para["additional"]:
                para["additional"] = (
                    para["additional"].replace("freq", " ").replace("  ", " ")
                )
            xyz = Gaussian.run(gjf, proname, para["name"])

            # 第二次优化的参数设置
            para["opt"] = (
                para["opt"].replace("calcfc", "calcall") + ",tight,maxstep=1,notrust"
            )
            if "int=superfine" not in para["additional"]:
                para["additional"] = para["additional"] + " int=superfine"
            para["additional"] = para["additional"] + "int=superfine"
            para["name"] = name + "_02"

            # 进行第二次几何优化
            gjf = Xyz.togjf(xyz, Code.topara(para))
            xyz = Gaussian.run(gjf, proname, para["name"])

            # 第三次几何优化的参数设置
            para["opt"] = ""
            if "freq" not in para["additional"]:
                para["additional"] = para["additional"] + " freq"
            para["name"] = name + "_03"
            para["solvent"] = solvent

            # 进行第三次几何优化
            gjf = Xyz.togjf(xyz, Code.topara(para))
            result = Gaussian.run(gjf, proname, para["name"], 1)

            if result[0]:
                return result[1], result[2]
            else:
                print("计算失败")
                return result[1]

    def gjfs_optimize(gjfs, proname, indexs, solvents):
        results = Executor.ThreadExecutor(
            Mechanism_Calculation.gjf_optimize,
            gjfs,
            [proname] * len(gjfs),
            indexs,
            solvents,
        )
        return results

    def xyz_optimize(xyz: str, proname: str, OPT: dict, index=1, solvent=""):
        Verify.Is([xyz, proname, OPT, index, solvent], [str, str, dict, int, str])

        gjf = Xyz.togjf(xyz, OPT)
        return Mechanism_Calculation.gjf_optimize(gjf, proname, index, solvent)

    def xyzs_optimize(
        xyzs: list,
        proname: str,
        OPT: dict,
        charges=[],
        names=[],
        indexs=[],
        solvents=[],
    ):
        Verify.Is(
            [xyzs, proname, OPT, charges, names, indexs],
            [list, str, dict, list, list, list],
        )

        if indexs == [] or len(indexs) < len(xyzs):
            indexs = [1] * len(xyzs)

        if solvents == [] or len(solvents) < len(xyzs):
            if OPT["solvent"] != "":
                solvents = [OPT["solvent"]] * len(xyzs)
            else:
                solvents = [""] * len(xyzs)

        gjfs = Xyz.togjfs(xyzs, OPT, charges, names)
        return Mechanism_Calculation.gjfs_optimize(gjfs, proname, indexs, solvents)
        # return Mechanism_Calculation.multrungjf(items)

    def smile_optimize(smile: str, proname: str, OPT: dict, index=1, solvent=""):
        Verify.Is([smile, proname, OPT, index, solvent], [str, str, dict, int, str])
        try:
            xyz = Smile.openbabeltoxyz(smile)
        except Exception:
            print("错误分子：" + smile)
            # messagebox.showinfo('错误提示', '分子式错误')
        else:
            return Mechanism_Calculation.xyz_optimize(xyz, proname, OPT, index, solvent)

    def smiles_optimize(
        smiles: list,
        proname: str,
        OPT: dict,
        charges=[],
        names=[],
        indexs=[],
        solvents=[],
    ):
        Verify.Is(
            [smiles, proname, OPT, charges, names, indexs],
            [list, str, dict, list, list, list],
        )
        xyzs = Smile.openbabeltoxyzs(smiles)
        return Mechanism_Calculation.xyzs_optimize(
            xyzs, proname, OPT, charges, names, indexs, solvents
        )
