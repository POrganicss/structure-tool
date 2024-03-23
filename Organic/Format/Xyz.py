import __init__
from Tool.Compute import *
from Tool.Verify import Verify
from Tool.File import *
from Tool.Code import *
from Format.Gjf import *
from openbabel import openbabel
from scipy.spatial.distance import cdist
import numpy as np

class Xyz:

    # 根据xyz坐标获取gjf文件mixset=[]
    def togjf(xyzs, para):
        g16_input = []
        g16_input.append('%chk='+para['name']+'.chk\n')
        g16_input.append('%nproc='+str(para['nproc'])+"\n")
        g16_input.append("%mem=" + para['mem']+"\n")
        g16_input.append(para['code']+"\n")

        # 匹配单坐标
        if isinstance(xyzs, str):
            g16_input.append("\n")
            g16_input.append(para['name']+'\n')
            g16_input.append("\n")
            g16_input.append(
                str(para['charge']) + " " + para['spin']+'\n')
            g16_input.append(Xyz.initialize(xyzs))
            g16_input.append("\n")

        # 匹配多坐标
        elif isinstance(xyzs, list):
            for xyz in xyzs:
                # 反应物坐标
                g16_input.append("\n")
                g16_input.append('Reaction or Product or TS'+'\n')
                g16_input.append("\n")
                g16_input.append(
                    para['charge'] + " " + para['spin']+'\n')
                g16_input.append(Xyz.initialize(xyz))
                g16_input.append("\n")

        if para['mixset'] != '':
            def justelement(xyz, eles):
                base_elements = Xyz.getelements(xyz)
                exclusive_sets = []

                for ele in eles:
                    exclusive_set = set(ele)
                    for e in ele:
                        if e in base_elements:
                            base_elements.remove(e)
                        else:
                            exclusive_set.remove(e)
                    exclusive_sets.append(' '.join(list(exclusive_set)))

                return ' '.join(list(base_elements)), exclusive_sets
            
            basets, exsets = justelement(
                xyzs, para['mixset'][:-1])
            g16_input.append(basets+' 0'+"\n")
            g16_input.append(para['method'].split('/')[1]+"\n")
            g16_input.append('****'+"\n")
            g16_input.append('****'+"\n")

            for i, st in enumerate(para['mixset'][-1]):
                g16_input.append(exsets[i]+' 0'+"\n")
                g16_input.append(st+"\n")
                g16_input.append('****'+"\n")

            g16_input.append("\n")

            for i, st in enumerate(para['mixset'][-1]):
                g16_input.append(exsets[i]+' 0'+"\n")
                g16_input.append(st+"\n")
                g16_input.append('****'+"\n")
            g16_input.remove('****'+"\n")
        g16_input.pop()
        g16_input.append("\n")
        g16_input.append("\n")
        return ''.join(g16_input)

    # 将一些xyz坐标分别转换为gjf文件
    def togjfs(xyzs, paramenters, charges, names):
        if charges is None or len(charges) < len(xyzs):
            charges = [0]*len(xyzs)
        bools = True
        if names is None or len(names) < len(xyzs):
            bools = False
        name = paramenters['name']
        gjfs = []
        for i, xyz in enumerate(xyzs):
            if bools:
                paramenters['name'] = names[i]
            else:
                paramenters['name'] = Code.getname([name, str(i+1)])
            paramenters['charge'] = charges[i]
            gjfs.append(Xyz.togjf(xyz, paramenters))

        return gjfs

    # 获取xyz数值字典
    def getcoords(xyz):
        contents = xyz.strip().splitlines()
        atoms_info = []
        atom_coords = []
        for content in contents:
            atoms = content.split()
            atoms_info.append(atoms[0])
            atom_coords.append(atoms[1:])
        return atoms_info, atom_coords

    # 将三维坐标转换为Gaussian中Z矩阵的格式

    def getmatrix(xyz, rvar=True, avar=True, dvar=True):
        matrix = []
        atomnames, xyzs = Xyz.getcoords(xyz)
        xyzarr = np.zeros([len(xyzs), 3])
        
        for i, xyz in enumerate(xyzs):
            xyzarr[i][0] = xyz[0]
            xyzarr[i][1] = xyz[1]
            xyzarr[i][2] = xyz[2]
            
        distmat = cdist(xyzarr, xyzarr)

        npart, ncoord = xyzarr.shape
        rlist = []  # list of bond lengths
        alist = []  # list of bond angles (degrees)
        dlist = []  # list of dihedral angles (degrees)
        if npart > 0:
            # Write the first atom
            matrix.append(atomnames[0])

            if npart > 1:
                # and the second, with distance from first
                n = atomnames[1]
                rlist.append(distmat[0][1])
                matrix.append('{:<3s} {:>4d}  {:11s}'.format(n, 1, 'R1'))

                if npart > 2:
                    n = atomnames[2]
                    rlist.append(distmat[0][2])
                    alist.append(Compute.getangle(
                        xyzarr[2], xyzarr[0], xyzarr[1]))
                    matrix.append(
                        '{:<3s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(n, 1, 'R2', 2, 'A1'))
                    if npart > 3:
                        for i in range(3, npart):
                            n = atomnames[i]
                            rlist.append(distmat[i-3][i])
                            alist.append(Compute.getangle(
                                xyzarr[i], xyzarr[i-3], xyzarr[i-2]))
                            dlist.append(Compute.getdihedral(
                                xyzarr[i], xyzarr[i-3], xyzarr[i-2], xyzarr[i-1]))
                            matrix.append('{:3s} {:>4d}  {:11s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(
                                n, i-2, 'R{:<4d}'.format(i), i-1, 'A{:<4d}'.format(i-1), i, 'D{:<4d}'.format(i-2)))
        matrix.append('')
        for i in range(npart-1):
            matrix.append('R{:<4d}   {:>11.5f}'.format(i+1, rlist[i]))
        for i in range(npart-2):
            matrix.append('A{:<4d}   {:>11.5f}'.format(i+1, alist[i]))
        for i in range(npart-3):
            matrix.append('D{:<4d}   {:>11.5f}'.format(i+1, dlist[i]))
        return '\n'.join(matrix)

    def xyztosdf(oldpath, newpath):
        conv = openbabel.OBConversion()  # 使用openbabel模块中的OBConversion函数，用于文件格式转换的
        conv.OpenInAndOutFiles(oldpath, newpath)  # 输入需要转换的文件的名字，以及定义转换后文件的文件名
        conv.SetInAndOutFormats("xyz", "sdf")  # 定义转换文件前后的格式
        conv.Convert()  # 执行转换操作
        conv.CloseOutFile()  # 转换完成后关闭转换后的文件，完成转换

    # 将xyz转换为sdf文件
    def tosdf(xyz):
        xyz = Xyz.toxyz(xyz)
        File.save(xyz, path+'/temp02.xyz')
        # Xyz.xyztosdf(path+'/temp02.xyz',path+'/temp02.sdf')

        conv = openbabel.OBConversion()  # 使用openbabel模块中的OBConversion函数，用于文件格式转换的
        # 输入需要转换的文件的名字，以及定义转换后文件的文件名
        conv.OpenInAndOutFiles(path+'/temp02.xyz', path+'/temp02.sdf')
        conv.SetInAndOutFormats("xyz", "sdf")  # 定义转换文件前后的格式
        conv.Convert()  # 执行转换操作
        conv.CloseOutFile()  # 转换完成后关闭转换后的文件，完成转换

        sdf = File.read(path+'/temp02.sdf')[:-1]
        # print(sdf)
        os.remove(path+'/temp02.sdf')
        return sdf

    def tosdfs(xyzs):
        sdfs = []
        for xyz in xyzs:
            sdfs.append(Xyz.tosdf(xyz))
        os.remove(path+'/temp02.xyz')
        return sdfs

    # 将xyz转换为<xyz文件>
    def toxyz(xyz):
        line = xyz.splitlines()
        return ''.join(str(len(line))+'\n'+'toxyz'+'\n'+xyz)

    # 将<xyz文件>转换成xyz
    def dexyz(xyz):
        line = xyz.splitlines()
        content = []
        for i in range(2, len(line)):
            if not line[i]:
                break
            content.append(line[i]+'\n')
        return ''.join(content)

    # 初始化xyz格式
    def initialize(xyz):
        # 给字符串前后加定量空格
        def addblank(st, num, dx=1):
            if dx == 1:
                for _ in range(num):
                    st = st+' '
            if dx == -1:
                for _ in range(num):
                    st = ' '+st
            return st

        contents = xyz.strip().splitlines()
        res = []
        for content in contents:
            me = content.split()
            if float(me[1]) < 0:
                res.append(''+addblank(me[0], 18-len(me[0]))+me[1])
            elif float(me[1]) >= 0:
                res.append(''+addblank(me[0], 19-len(me[0]))+me[1])
                
            if float(me[2]) < 0:
                res.append(addblank(me[2], 3, -1))
            elif float(me[2]) >= 0:
                res.append(addblank(me[2], 4, -1))

            if float(me[3]) < 0:
                res.append(addblank(me[3], 3, -1))
            elif float(me[3]) >= 0:
                res.append(addblank(me[3], 4, -1))
            res.append('\n')
        return ''.join(res)

   
    def getpoints(xyz, *atoms):
        xyzs = Xyz.todata(xyz)[1]
        if not atoms:
            # 如果 atoms 为空，则返回所有坐标
            return xyzs
        # 检查每个原子索引，并返回相应的坐标
        return [xyzs[atom - 1] for atom in atoms]

    def getelements(xyz):  # 获取xyz中的所有元素
        elements = set()
        for line in xyz.splitlines():
            string = line.split()
            if len(string) == 4 and not Verify.isnum(string[0]) and Verify.isnums(string[1:]):
                elements.add(string[0])
        return elements


    # 获取参数
    def getparameter(xyzs, *args):
        unique_xyzs = []

        for xyz in xyzs:
            lines = xyz.splitlines()
            coords = [list(map(float, lines[coord - 1].split()[1:4]))
                    for coord in args]
            unique_xyzs.append(coords)

        result = []

        for coords in unique_xyzs:
            if len(args) == 2:
                distance = -Compute.getdistance(coords[0], coords[1])
                result.append(distance)
            elif len(args) == 3:
                angle = Compute.getangle(coords[0], coords[1], coords[2])
                result.append(angle)
            elif len(args) == 4:
                dihedral = Compute.getdihedral(
                    coords[0], coords[1], coords[2], coords[3])
                result.append(dihedral)

        return result
            

        