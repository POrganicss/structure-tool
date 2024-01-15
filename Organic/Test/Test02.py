import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random

def plot():
    def update_plot(frame, x_data, y_data, line, ax):
        # 模拟新数据
        new_x = frame
        new_y = random.randint(0, 1) * frame

        # 更新数据
        x_data.append(new_x)
        y_data.append(new_y)

        # 设置数据
        line.set_xdata(x_data)
        line.set_ydata(y_data)

        # 动态调整坐标轴范围
        ax.set_xlim(0, max(x_data) + 5)  # 增加一些空间，使得图形更美观
        ax.set_ylim(0, max(y_data) + 2)

        # 显式设置坐标轴标签
        ax.set_xlabel('X轴标签')
        ax.set_ylabel('Y轴标签')

        # 显式调用绘图更新
        ax.figure.canvas.draw()

        return line,


    # 创建初始空白图表
    fig, ax = plt.subplots()
    x_data, y_data = [], []
    line, = ax.plot(x_data, y_data)  # 创建一个空曲线

    # 设置坐标轴标签
    ax.set_xlabel('X轴标签')
    ax.set_ylabel('Y轴标签')

    # 创建动画
    animation = FuncAnimation(fig, update_plot, fargs=(
        x_data, y_data, line, ax), interval=5, blit=True)

    plt.show()  # 显示图表

plot()