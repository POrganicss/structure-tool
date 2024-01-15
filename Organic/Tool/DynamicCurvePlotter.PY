import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class DynamicCurvePlotter:
    def __init__(self):
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], marker='x', linestyle='-')
        self.data_points = []

        # 设置图表标题和轴标签
        plt.title("Dynamic Curve Plot")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")

        # 使用FuncAnimation创建动画，设置save_count为一个合适的值
        self.ani = FuncAnimation(
            self.fig, self.update, init_func=self.init, blit=True, save_count=10)

    def init(self):
        self.line.set_data([], [])
        return self.line,

    def update(self, frame):
        for _ in range(100):
            # 调用其他函数计算新的数据点
            new_data_point = self.calculate_new_data_point()

            # 将新数据点添加到列表中
            self.data_points.append(new_data_point)

            # 更新线条对象
            x, y = zip(*self.data_points)
            self.line.set_data(x, y)
        return self.line,

    def calculate_new_data_point(self):
        # 在这里调用其他函数计算新的数据点
        # 这里简单地以时间作为 x 轴，随机数作为 y 轴
        import random
        x = len(self.data_points) + 1
        y = random.randint(1, 10)
        return x, y

    def show_plot(self):
        plt.show()


# 示例用法
plotter = DynamicCurvePlotter()
plotter.show_plot()
