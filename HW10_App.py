import sys
import os #HW10 addition
import numpy as np
from pathlib import Path
from HW10_OOP_stem import trendlineClass as trendLn

import PyQt5.QtWidgets as qtw

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import  NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class ApplicationWindow(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.FilePath=os.getcwd() #place to store directory for data file dialog box.  os.getcwd() get current working directory

        #region manually build the main window and widgets rather than using QtDesigner
        self.layout =qtw.QVBoxLayout() #main layout for the widgets
        self.setLayout(self.layout)
        self.gbInput=qtw.QGroupBox("Input") #group box for input controls
        self.gbOutput=qtw.QGroupBox() #groupbox for output display
        self.gbInput.setSizePolicy(qtw.QSizePolicy.Fixed,qtw.QSizePolicy.Fixed)

        self.layout.addWidget(self.gbInput) #add input group box to main layout
        self.layout.addWidget(self.gbOutput) #add output group box to main layout

        self.gbInputLayout=qtw.QGridLayout() #grid layout for input group box
        self.gbOutputLayout=qtw.QGridLayout() #grid layout for output group box

        #set up controls for input
        self.lbl_Filename=qtw.QLabel("Filename:")
        self.tb_Filename=qtw.QLineEdit("Data File 1.txt")
        self.rb_LinearFit=qtw.QRadioButton("Linear Fit",checked=True)
        self.rb_QuadraticFit=qtw.QRadioButton("Quadratic Fit")
        self.rb_CubicFit=qtw.QRadioButton("Cubic Fit")
        self.rb_ExponentialFit=qtw.QRadioButton("Exponential Fit") #added this for HW10
        self.rb_AllThree=qtw.QRadioButton("All") #modified this for HW10
        self.btn_LoadAndCalculate=qtw.QPushButton("Load and Calculate")

        #add controls to input group box
        self.gbInput.setLayout(self.gbInputLayout)
        self.gbInputLayout.addWidget(self.lbl_Filename,1,1)
        self.gbInputLayout.addWidget(self.tb_Filename,1,2,1,3)
        self.gbInputLayout.addWidget(self.rb_LinearFit,2,1)
        self.gbInputLayout.addWidget(self.rb_QuadraticFit,2,2)
        self.gbInputLayout.addWidget(self.rb_CubicFit,2,3)
        self.gbInputLayout.addWidget(self.rb_ExponentialFit,2,4) #added this for HW10
        self.gbInputLayout.addWidget(self.rb_AllThree,2,5)
        self.gbInputLayout.addWidget(self.btn_LoadAndCalculate,3,2,1,2)

        #set up controls for output
        self.lbl_Eq=qtw.QLabel("Equation")
        self.tb_Eq=qtw.QLineEdit("")

        #set up the canvas for the matplotlib display
        self.canvas = FigureCanvasQTAgg(Figure(figsize=(10, 20),tight_layout=True, frameon=True))
        self.toolbar=NavigationToolbar(self.canvas,self) #HW10 addition

        #add controls to output group box
        self.gbOutput.setLayout(self.gbOutputLayout)
        self.gbOutputLayout.addWidget(self.lbl_Eq,1,1)
        self.gbOutputLayout.addWidget(self.tb_Eq,1,2)
        self.gbOutputLayout.addWidget(self.canvas,2,1,1,2)
        self.gbOutputLayout.addWidget(self.toolbar,3,1,1,2) #HW10 addition

        self.ax = self.canvas.figure.add_subplot()
        #endregion

        # create a trendline object
        x = np.array([0.05, 0.11, 0.15, 0.31, 0.46, 0.52, 0.70, 0.74, 0.82, 0.98, 1.17])
        y = np.array([0.956, 1.09, 1.332, 0.717, 0.771, 0.539, 0.378, 0.370, 0.306, 0.242, 0.104])
        self.TL=trendLn(x, y)

        #setup signals and slots
        self.btn_LoadAndCalculate.clicked.connect(self.OpenDataFile)
        self.rb_LinearFit.clicked.connect(self.DoPlot)
        self.rb_QuadraticFit.clicked.connect(self.DoPlot)
        self.rb_CubicFit.clicked.connect(self.DoPlot)
        self.rb_ExponentialFit.clicked.connect(self.DoPlot) #HW10 addition
        self.rb_AllThree.clicked.connect(self.DoPlot)

        self.DoPlot()  #note that linear fit is the default upon running program

        self.resize(1000,1000)

    def DoPlot(self):
        axes=self.ax
        axes.clear()
        #added TexEq to return for displaying a nicely formatted equation in the legend
        if self.rb_LinearFit.isChecked():
            linx,liny,linRSq,lineq, linTexEq, linlbl=self.TL.PlotLeastSquares(1, axes=self.ax, showpoints=True, npoints=500)
            self.tb_Eq.setText(lineq)
        elif self.rb_QuadraticFit.isChecked():
            quadx,quady,quadRSq,quadeq, quadTexEq, quadlbl=self.TL.PlotLeastSquares(2, axes=self.ax, showpoints=True, npoints=500)
            self.tb_Eq.setText(quadeq)
        elif self.rb_CubicFit.isChecked():
            cubx,cuby,cubRSq,cubeq, cubTexEq, cubelbl=self.TL.PlotLeastSquares(3, axes=self.ax, showpoints=True, npoints=500)
            self.tb_Eq.setText(cubeq)
        elif self.rb_ExponentialFit.isChecked(): #HW10 addition
            expx, expy,expRSq, expeq, expTexEq, explbl=self.TL.PlotLeastSquares(3, axes=self.ax, showpoints=True, npoints=500, exp=True)
            self.tb_Eq.setText(expeq)
        elif self.rb_AllThree.isChecked():
            linx,liny,linRSq,lineq, linTexEq, linlbl=self.TL.GetPlotInfo(1, npoints=500)
            quadx,quady,quadRSq,quadeq, quadTexEq, quadlbl=self.TL.GetPlotInfo(2, npoints=500)
            cubx,cuby,cubRSq,cubeq,cubeTexEq, cubelbl=self.TL.GetPlotInfo(3, npoints=500)
            expx,expy,expRSq,expeq,expTexEq, explbl=self.TL.GetPlotInfo(3, npoints=500, exp=True)

            axes.clear()
            axes.plot(linx, liny, linestyle='dashed', color='black', linewidth='2', label="" + linlbl +" {"+linTexEq+"}"+ r'($R^2={:0.3f}$)'.format(linRSq))
            axes.plot(quadx, quady, linestyle='dotted', color='black', linewidth='2', label="" + quadlbl+" {"+quadTexEq+"}" + r'($R^2={:0.3f}$)'.format(quadRSq))
            axes.plot(cubx, cuby, linestyle='dashdot', color='black', linewidth='2', label="" + cubelbl +" {"+cubeTexEq+"}"+ r'($R^2={:0.3f}$)'.format(cubRSq))
            axes.plot(expx, expy, linestyle='solid', color='black', linewidth='2', label="" + explbl+" {"+expTexEq+"}" + r'($R^2={:0.3f}$)'.format(expRSq))
            axes.plot(self.TL.x, self.TL.y, linestyle='none', marker='o', markerfacecolor='orange', markeredgecolor='black', markersize=10)
            axes.grid(axis='both')
            axes.tick_params(axis='both', direction='in', grid_linewidth=1, grid_linestyle='dashed', grid_alpha=0.5)
            axes.set(xlabel='X values', ylabel='Y values')
            axes.legend( fontsize=18) #HW10 addition of fontsize=18
            self.tb_Eq.clear()

        self.canvas.draw()

    def OpenDataFile(self): #illustrates use of QFileDialog to open a data file
        fname=qtw.QFileDialog.getOpenFileName(self, 'Open file',self.FilePath, 'text files (*.txt)')
        self.FilePath=str(Path(fname[0]).parents[0])+'\\'
        self.tb_Filename.setText(fname[0])
        self.TL.x , self.TL.y=np.loadtxt(fname[0], skiprows=1, unpack=True)
        self.DoPlot()

if __name__ == "__main__":
    qapp = qtw.QApplication(sys.argv)
    app = ApplicationWindow()
    app.show()
    qapp.exec_()