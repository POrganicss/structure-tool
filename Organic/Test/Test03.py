from __future__ import print_function  # 导入未来模块的print_function，确保与Python 2和3的兼容性
import sys  # 导入sys模块，用于系统特定的参数和函数
from openeye import oechem  # 导入openeye包中的oechem模块
from openeye import oedocking  # 导入openeye包中的oedocking模块

def main(argv=[__name__]):  # 定义主函数，接受argv作为可选参数
    itf = oechem.OEInterface(InterfaceData)  # 使用InterfaceData创建OEInterface对象
    oedocking.OEDockMethodConfigure(itf, "-method")  # 配置对接方法为：
    oedocking.OESearchResolutionConfigure(itf, "-resolution")  # 配置搜索分辨率

    if not oechem.OEParseCommandLine(itf, argv):  # 解析命令行参数
        return 1

    imstr = oechem.oemolistream(itf.GetString("-in"))  # 打开输入分子文件流
    omstr = oechem.oemolostream(itf.GetString("-out"))  # 打开输出分子文件流
    receptor = oechem.OEGraphMol()  # 创建一个空的图分子用于受体

    if not oedocking.OEReadReceptorFile(receptor, itf.GetString("-receptor")):  # 读取受体文件
        oechem.OEThrow.Fatal("Unable to read receptor")

    dockMethod = oedocking.OEDockMethodGetValue(itf, "-method")  # 获取对接方法值
    dockResolution = oedocking.OESearchResolutionGetValue(itf, "-resolution")  # 获取搜索分辨率值
    dock = oedocking.OEDock(dockMethod, dockResolution)  # 创建一个对接对象
    dock.Initialize(receptor)  # 使用受体初始化对接过程

    for mcmol in imstr.GetOEMols():  # 遍历输入流中的每个多构象分子
        print("docking", mcmol.GetTitle())  # 打印正在对接的分子标题
        dockedMol = oechem.OEGraphMol()  # 创建一个空的图分子用于对接后的分子
        dock.DockMultiConformerMolecule(dockedMol, mcmol)  # 对多构象分子进行对接
        sdtag = oedocking.OEDockMethodGetName(dockMethod)  # 获取对接方法的名称
        oedocking.OESetSDScore(dockedMol, dock, sdtag)  # 将对接分数设置为SD标签
        dock.AnnotatePose(dockedMol)  # 注释对接构象
        oechem.OEWriteMolecule(omstr, dockedMol)  # 将对接后的分子写入输出流

    return 0

InterfaceData = """
!PARAMETER -receptor
!ALIAS -rec
!TYPE string
!REQUIRED true
!LEGAL_VALUE *.oeb
!LEGAL_VALUE *.oeb.gz
!BRIEF A receptor file the molecules pass to the -in flag will be docked to
!END
!PARAMETER -in
!TYPE string
!REQUIRED true
!BRIEF Multiconformer file of molecules to be docked.
!END
!PARAMETER -out
!TYPE string
!REQUIRED true
!BRIEF Docked molecules will be written to this file
!END
"""

if __name__ == "__main__":  # 检查脚本是否作为主程序运行
    sys.exit(main(sys.argv))  # 以主函数的返回值退出脚本
''' 


1.导入模块和定义函数：
    从openeye模块中导入oechem和oedocking模块。
    定义了一个名为main的函数，该函数接受一个参数argv，默认为[__name__]。

2.配置接口参数和命令行解析：
    创建了一个OEInterface对象itf，并使用预定义的InterfaceData进行初始化。
    使用oedocking.OEDockMethodConfigure和oedocking.OESearchResolutionConfigure配置了对接方法和搜索分辨率。
    调用oechem.OEParseCommandLine函数解析命令行参数，并判断是否成功。如果失败则返回1，表示程序异常退出。

3.打开输入输出流和读取受体文件：
    使用oechem.oemolistream打开输入分子文件流imstr，使用oechem.oemolostream打开输出分子文件流omstr。
    创建一个空的图分子对象receptor，并尝试读取受体文件，如果读取失败则抛出异常并退出程序。

4.初始化对接过程：
    获取对接方法和搜索分辨率的值，并使用这些值创建一个对接对象dock。
    使用受体对象receptor初始化对接过程。
    
5.对每个多构象分子进行对接：
    遍历输入流中的每个多构象分子mcmol。
    对当前分子进行对接，将对接后的构象存储在dockedMol中。
    获取对接方法的名称并设置为SD标签。
    对对接后的构象进行注释，并将其写入输出流。

6.返回结果：
    主函数返回0，表示程序正常退出。
整体而言，该代码通过命令行参数配置对接过程，并针对每个输入的多构象分子执行对接操作，最后将对接后的构象写入输出文件中。 '''