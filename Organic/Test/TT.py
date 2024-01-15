from math import cos, pi, sin
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QOpenGLWidget
from PyQt5.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLUT import *
from openbabel import pybel

class MolViewer(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.mol = None

    def load_molecule(self, input_file):
        mol = next(pybel.readfile('xyz', input_file))
        self.mol = mol

    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)

        glShadeModel(GL_FLAT)
        glClearColor(1.0, 1.0, 1.0, 1.0)

    def paintGL(self):
        if self.mol:
            self.mol.OBMol.SetDimension(3)
            self.mol.OBMol.AddHydrogens()
            self.mol.OBMol.ConnectTheDots()
            self.mol.OBMol.PerceiveBondOrders()
            self.mol.OBMol.EmbedMultipleConfs()
            self.mol.OBMol.CalcSASA()
            gl_data = self.mol.write("sdf").splitlines()
            gl_data = '\n'.join(gl_data[3:-1])
            gl_format = self.mol.write("sdf").splitlines()
            gl_format = gl_format[1]
            self.mol.OBMol.DeleteHydrogens()
            gl_data += " " + gl_format
            self.render_molecule(gl_data)

    def render_molecule(self, gl_data):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # 设置视图
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, self.width() / self.height(), 1, 100)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        gluLookAt(0, 0, -15, 0, 0, 0, 0, 1, 0)

        # 渲染分子结构
        glCallList(glGenLists(1))
        glNewList(1, GL_COMPILE)
        self.draw_benzene_ring()
        glEndList()
        glCallList(1)

        glFlush()

    def draw_benzene_ring(self):
        radius = 1.5
        sides = 6

        glBegin(GL_POLYGON)
        glColor3f(1, 0, 0)  # 设置颜色为红色
        for i in range(sides):
            angle = 2.0 * pi * i / sides
            x = radius * cos(angle)
            y = radius * sin(angle)
            glVertex3f(x, y, 0.0)
        glEnd()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setGeometry(100, 100, 800, 600)

        layout = QVBoxLayout()
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        self.mol_viewer = MolViewer(self)
        layout.addWidget(self.mol_viewer)

        input_xyz_file = 'benzene.xyz'  # 替换为您的XYZ文件名
        self.mol_viewer.load_molecule(input_xyz_file)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
