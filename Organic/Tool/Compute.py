import math
import mplcursors
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import multiprocessing as mp


class Compute:

    # 获取任意两个原子间的距离
    def getdistance(coords1, coords2):
        # 使用zip将两个坐标列表的对应元素组合在一起，然后计算平方差，最后求和
        squared_diff = sum((p2 - p1)**2 for p1,p2 in zip(coords1, coords2))
        
        # 使用math.sqrt计算欧几里得距离
        distance = math.sqrt(squared_diff)
        
        return distance

    # 获取任意三点构成角的角度
    def getangle(coords1, coords2, coords3):
        # 计算两个向量
        vector1 = np.array(coords1) - np.array(coords2)
        vector2 = np.array(coords3) - np.array(coords2)
        
        # 计算点积和叉积
        dot_product = np.dot(vector1, vector2)
        cross_product = np.linalg.norm(np.cross(vector1, vector2))
        
        # 计算角度
        cos_theta = dot_product / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
        sin_theta = cross_product / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
        
        # 使用反三角函数计算角度值
        theta_rad = np.arctan2(sin_theta, cos_theta)
        theta_deg = np.degrees(theta_rad)
        
        return theta_deg 

    # 获取任意四点构成的二面角度
    def getdihedral(coords1, coords2, coords3, coords4):
        # 计算连接向量
        vector_ji = np.array(coords2) - np.array(coords1)
        vector_kj = np.array(coords3) - np.array(coords2)
        vector_lk = np.array(coords4) - np.array(coords3)

        # 计算法向量并归一化
        normal_vector_1 = np.cross(vector_ji, vector_kj)
        normal_vector_1 /= np.linalg.norm(normal_vector_1)

        normal_vector_2 = np.cross(vector_lk, vector_kj)
        normal_vector_2 /= np.linalg.norm(normal_vector_2)

        # 计算辅助向量并归一化
        cross_product_m1 = np.cross(normal_vector_1, vector_kj)
        cross_product_m1 /= np.linalg.norm(vector_kj)

        # 计算点乘
        dot_product_x = np.dot(normal_vector_1, normal_vector_2)
        dot_product_y = np.dot(cross_product_m1, normal_vector_2)

        # 计算夹角
        dihedral_angle_rad = np.arctan2(dot_product_y, dot_product_x)
        dihedral_angle_deg = -180.0 - np.degrees(dihedral_angle_rad)
        
        # 调整夹角范围
        dihedral_angle_deg = (dihedral_angle_deg + 180.0) % 360.0 - 180.0
        return dihedral_angle_deg
    
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
