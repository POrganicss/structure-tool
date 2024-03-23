import math
import numpy as np
import sympy as sp
from scipy.spatial.transform import Rotation
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.integrate import cumtrapz

import numpy as np
import importlib
#提供数学的计算工具
class Compute:
    def get_lib():
        """动态选择数值计算库。"""
        try:
            cp = importlib.import_module('cupy')  # 尝试导入CuPy库
            return cp if cp.is_available() else np  # 如果CuPy可用，则返回CuPy，否则返回NumPy
        except ImportError:
            return np  # 如果导入CuPy失败，则退回到NumPy
    
    # 获取任意两个原子间的距离
    def getdistance(coords1, coords2):
        """计算两个原子之间的欧几里得距离。"""
        xp = Compute.get_lib()  # 确定使用NumPy或CuPy
        squared_diff = xp.sum((coord2 - coord1)**2 for coord1, coord2 in zip(coords1, coords2))  # 计算坐标差的平方和
        distance = xp.sqrt(squared_diff)  # 计算距离的平方根
        return distance.get() if xp != np else distance  # 如果使用CuPy，则需要使用get()方法获取值

    # 获取任意三点构成角的角度
    def getangle(coords1, coords2, coords3):
        """计算由三个原子形成的角度。"""
        xp = Compute.get_lib()  # 确定使用NumPy或CuPy
        vector_a = xp.array(coords1) - xp.array(coords2)  # 计算第一个向量
        vector_b = xp.array(coords3) - xp.array(coords2)  # 计算第二个向量
        dot_prod = xp.dot(vector_a, vector_b)  # 计算两个向量的点积
        cross_prod_magnitude = xp.linalg.norm(xp.cross(vector_a, vector_b))  # 计算两个向量的叉积的模长
        cos_angle = dot_prod / (xp.linalg.norm(vector_a) * xp.linalg.norm(vector_b))  # 计算余弦值
        sin_angle = cross_prod_magnitude / (xp.linalg.norm(vector_a) * xp.linalg.norm(vector_b))  # 计算正弦值
        angle_rad = xp.arctan2(sin_angle, cos_angle)  # 通过arctan2得到角度的弧度值
        return xp.degrees(angle_rad).get() if xp != np else xp.degrees(angle_rad)  # 将弧度转换为度并返回


    # 获取任意四点构成的二面角度
    def getdihedral(coords1, coords2, coords3, coords4):
        """计算由四个原子形成的二面角。"""
        xp = Compute.get_lib()  # 确定使用NumPy或CuPy
        
        vector_ji = xp.array(coords2) - xp.array(coords1)  # 计算第一个向量
        vector_kj = xp.array(coords3) - xp.array(coords2)  # 计算第二个向量
        vector_lk = xp.array(coords4) - xp.array(coords3)  # 计算第三个向量
        
        normal_vector1 = xp.cross(vector_ji, vector_kj)  # 计算第一个平面的法向量
        normal_vector1 /= xp.linalg.norm(normal_vector1)  # 归一化第一个法向量
        normal_vector2 = xp.cross(vector_lk, vector_kj)  # 计算第二个平面的法向量
        normal_vector2 /= xp.linalg.norm(normal_vector2)  # 归一化第二个法向量
        
        cross_product_m1 = xp.cross(normal_vector1, vector_kj)  # 计算辅助向量
        cross_product_m1 /= xp.linalg.norm(vector_kj)  # 归一化辅助向量
        
        x_dot = xp.dot(normal_vector1, normal_vector2)  # 计算法向量间的点积
        y_dot = xp.dot(cross_product_m1, normal_vector2)  # 计算辅助向量与法向量的点积
        
        dihedral_rad = xp.arctan2(y_dot, x_dot)  # 计算二面角的弧度值
        
        dihedral_deg = -180.0 - xp.degrees(dihedral_rad)  # 将弧度转换为度，并调整范围
        dihedral_deg = (dihedral_deg + 180.0) % 360.0 - 180.0  # 确保二面角的值在-180到180度之间
        return dihedral_deg.get() if xp != np else dihedral_deg  # 返回二面角的度数值
    
    
    def getmatrix(xyz, rvar=True, avar=True, dvar=True):
        from Format.Xyz import Xyz  # 导入Xyz类
        from scipy.spatial.distance import cdist  # 从SciPy库导入cdist函数，用于计算距离矩阵
        
        xp = Compute.get_lib()  # 获取合适的数值计算库（NumPy或CuPy）
        atomnames, xyzs = Xyz.getcoords(xyz)  # 从Xyz对象获取原子名和坐标
        xyzarr = xp.zeros([len(xyzs), 3])  # 初始化一个适当大小的数组，用于存储坐标

        # 转换数据格式：确保每个原子的坐标都被正确填充到xyzarr数组中
        for i, xyz in enumerate(xyzs):
            xyzarr[i, :] = xp.asarray(xyz)  # 转换每个坐标为适用的数值数组并存储

        # 使用cdist计算所有原子之间的距离矩阵，注意需要将数据转换为NumPy数组，因为cdist不接受CuPy数组
        distmat = cdist(xyzarr.get() if xp != np else xyzarr, xyzarr.get() if xp != np else xyzarr)

        # 初始化列表，用于存储分子结构参数
        rlist, alist, dlist = [], [], []  # 分别用于存储键长、键角和二面角
        matrix = []  # 用于存储输出结果的矩阵

        npart, ncoord = xyzarr.shape  # 获取原子总数和坐标维度
        

        if npart > 0:  # 如果有原子存在
            matrix.append(atomnames[0])  # 添加第一个原子名称
            if npart > 1:  # 如果至少有两个原子
                rlist.append(distmat[0][1])  # 添加第一个原子到第二个原子的距离
                matrix.append(f'{atomnames[1]:<3s} {1:>4d}  {"R1":>11s}')  # 格式化并添加第二个原子信息

                if npart > 2:  # 如果至少有三个原子
                    rlist.append(distmat[0][2])  # 添加第一个原子到第三个原子的距离
                    alist.append(Compute.getangle(xyzarr[2], xyzarr[0], xyzarr[1]))  # 计算并添加角度信息
                    matrix.append(f'{atomnames[2]:<3s} {1:>4d}  {"R2":>11s} {2:>4d}  {"A1":>11s}')  # 格式化并添加第三个原子信息

                    if npart > 3:  # 如果至少有四个原子
                        for i in range(3, npart):  # 遍历剩余的原子
                            rlist.append(distmat[i-3][i])  # 添加距离信息
                            alist.append(Compute.getangle(xyzarr[i], xyzarr[i-3], xyzarr[i-2]))  # 添加角度信息
                            dlist.append(Compute.getdihedral(xyzarr[i], xyzarr[i-3], xyzarr[i-2], xyzarr[i-1]))  # 添加二面角信息
                            matrix.append(f'{atomnames[i]:3s} {i-2:>4d}  {"R{i}":>11s} {i-1:>4d}  {"A{i-1}":>11s} {i:>4d}  {"D{i-2}":>11s}')  # 格式化并添加原子信息

        matrix.append('')  # 添加空行作为分隔符
        for i in range(npart-1):
            matrix.append('R{:<4d}   {:>11.5f}'.format(i+1, rlist[i]))
        for i in range(npart-2):
            matrix.append('A{:<4d}   {:>11.5f}'.format(i+1, alist[i]))
        for i in range(npart-3):
            matrix.append('D{:<4d}   {:>11.5f}'.format(i+1, dlist[i]))


        return '\n'.join(matrix)  # 将结果列表转换为字符串并返回
   
        
    def align_complex_structure(simple_coords, simple_connectivity_matrix, complex_coords, complex_connectivity_matrix, core_indices, adjusted_core_coords):
                
        def calculate_rotation_matrix(source, target):
            # 计算旋转矩阵
            rotation_matrix, _ = Rotation.align_vectors(source, target)
            return rotation_matrix.as_matrix()
                
        # 提取核心区域的坐标
        simple_core_coords = simple_coords[core_indices]
        complex_core_coords = complex_coords[core_indices]

        # 计算质心
        simple_center = np.mean(simple_core_coords, axis=0)
        complex_center = np.mean(complex_core_coords, axis=0)

        # 平移，使得质心重合
        translation_vector = complex_center - simple_center
        complex_coords = complex_coords.astype(np.float64) - translation_vector

        # 计算旋转矩阵
        rotation_matrix = calculate_rotation_matrix(simple_core_coords, complex_core_coords)

        # 应用平移和旋转到整个复杂结构
        aligned_complex_coords = np.dot(complex_coords - complex_center, rotation_matrix.T) + simple_center

        # 将核心区域调整为新的坐标
        aligned_complex_coords[core_indices] = adjusted_core_coords

        return aligned_complex_coords
    
    # 计算数据的一阶导数 
    def getderivative(y, dx=1):
        """
        计算给定数据点的数值导数，假设x值是等间距的。

        :param y: 输入数据点的y坐标，numpy数组形式
        :param dx: x值之间的间隔，默认为1
        :return: 导数值的numpy数组
        """
        x = np.arange(len(y)) * dx
        f = interp1d(x, y, kind='cubic', fill_value="extrapolate")
        dydx = np.array([derivative(f, xi, dx=1e-6) for xi in x])
        
        return dydx

    def getantiderivative(y, dx=1):
        """
        计算给定数据点的数值原函数，假设x值是等间距的。

        :param y: 输入数据点的y坐标，numpy数组形式
        :param dx: x值之间的间隔，默认为1
        :return: 原函数值的numpy数组
        """
        x = np.arange(len(y)) * dx
        Y_int = cumtrapz(y, x, initial=0)
        
        return Y_int
        
    def integrate_derivative(derivative_data, initial_value=0):
        original_data = [initial_value]

        n = len(derivative_data)
        delta_x = 1  # 假设离散间隔为1，可以根据实际情况调整

        for i in range(n - 2):  # 修改这里的范围
            x_i = i * delta_x
            x_i1 = (i + 1) * delta_x

            integral_value = (
                delta_x / 6) * (derivative_data[i] + 4 * derivative_data[i + 1] + derivative_data[i + 2]) + original_data[-1]
            original_data.append(integral_value)

        return original_data

    def getoriginal(derivative_data, initial_value=0):
        """
        根据导数数据计算原函数数据

        参数:
        - derivative_data: 包含导数数据的列表
        - initial_value: 初始值，默认为0

        返回:
        一个包含原函数数据的列表
        """
        x = sp.symbols('x')
        integral_expr = initial_value

        original_data = [initial_value]
        for derivative in derivative_data:
            integral_expr += float(derivative)
            original_data.append(sp.integrate(integral_expr, x))

        return original_data
    
    def getderivative(original_data, n=1):
        """
        计算原函数的 n 阶导数数据

        参数:
        - original_data: 包含原函数数据的列表
        - n: 阶数

        返回:
        一个包含 n 阶导数数据的列表
        """
        x = sp.symbols('x')
        original_expr = sp.interpolate(original_data, x)
        nth_derivative_expr = sp.diff(original_expr, x, n)

        # 将符号表达式转换为可调用的函数
        nth_derivative_func = sp.lambdify(x, nth_derivative_expr, 'numpy')

        # 生成 x 值，以便计算导数
        x_values = np.linspace(0, len(original_data) - 1, len(original_data))

        # 计算 n 阶导数数据
        nth_derivative_data = nth_derivative_func(x_values)

        return nth_derivative_data

    def getstatistics(data):
        """
        计算给定数据的均值、方差和标准差

        参数:
        - data: 包含数据的列表

        返回:
        一个包含均值、方差和标准差的元组
        """
        mean_value = np.mean(data)
        variance_value = np.var(data)
        std_dev_value = np.std(data)

        return mean_value, variance_value, std_dev_value

    def renormalization(*args):
        index = len(args[0])
        
        for i in range(index-1):
            for j in range(index-i-1):
                if args[0][j] < args[0][j+1]:
                    args[0][j], args[0][j+1] = args[0][j+1], args[0][j]
                    args[1][j], args[1][j+1] = args[1][j+1], args[1][j]
                    if len(args)==3:
                        args[2][j], args[2][j+1] = args[2][j+1], args[2][j]
        for i in len(args[0]):
            while args[0][i]-args[0][i+1]<0.00001:
                if args[-1][i]<args[0][i+1]:
                    del args[0][i+1]
                    del args[1][i+1]
                    if len(args)==3:
                        del args[2][i+1]
        print(len(args[0]))
        return args



