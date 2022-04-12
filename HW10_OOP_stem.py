import numpy as np
import matplotlib.pyplot as pyplot
import scipy as sp
import scipy.optimize as opt

from math import *
import math

class trendlineClass():
    def __init__(self, xdata, ydata):
        """
        The class should be constructed with some x and y data points to be fit
        :param xdata:
        :param ydata:
        """
        self.x=xdata
        self.y=ydata

    def RSquared(self, a, exp=False):
        '''
        To calculate the coefficient of determination (R**2) value for a set of x,y data and fitting function
        :param a:  the coefficients for the fitting function
        :param exp:  is it an exponential fit?
        :return: the R**2 value
        '''
        AvgY=np.mean(self.y) #calculates the average value of y
        yHat=np.array(self.yFit(a,exp))  #calculates the y values at self.x points for the fitting function
        SSRes=np.sum((yHat-self.y)**2)
        SSTot=np.sum((self.y-AvgY)**2)
        RSq = 1-SSRes/SSTot
        return RSq

    def yFit(self, a, exp=False):
        return [self.Exp(x,a) if exp else self.Poly(x,a) for x in self.x]

    def Exp(self, xval, a):
        """
        Calculate the value of an exponential fit to the data (i.e., a[0]+a[1]*exp(a[2]*xval))
        :param xval: the x value where I want to evaluate the function
        :param a: the coefficients for the fit
        :return: the value of the function
        """
        #$JES MISSING CODE HERE$
        pass

    def SSE(self, fn):
        """
        Calculates the sum of squared errors (SSE) for a fitting function to the data stored in self.x and self.y
        :param fn: the function  (e.g., a+b*exp(c*x))
        :return:
        """
        iterable = (fn(x) for x in self.x)
        yhat=np.fromiter(iterable,float)  # the value of the fitting function at each value of self.x
        sse=np.sum((self.y-yhat)**2)
        return sse

    def Poly(self,xval, a):
        '''
        calculates the value for a polynomial given a value for x and the coefficients of the polynomial.
        :param a:  the list of coefficients from numpy polyfit
        :return:  value of the polynomial at xval
        '''
        y=np.poly1d(a)
        return y(xval)

    def LeastSquares(self, power):
        '''
        Calculates the coefficients for a polynomial of degree power to best
        fit a data set (x,y) using the least squares approach.
        :param x: the independent variable of the data set
        :param y: the value of the function evaluated at each value of x
        :param power: the degree of the polynomial
        :return: the array of coefficients from numpy.polyfit
        '''
        a=np.polyfit(self.x, self.y, power)
        return a

    def OptimizeFit(self):
        """
        I need to create an objective function for scipy.optimize.minimize to minimize.  In this case of fitting the
        data with an exponential equation, the objective function should be the sum of squared errors (SSE) for the
        exponential fit vs. the self.x,self.y data.
        :return: the coefficients for the exponential fit that minimizes the objective function
        """
        #$JES MISSING CODE HERE$
        pass

    def GetEq(self,aa, exp=False):
        """
        This function creates an equation for either a polynomial fit or an exponential fit.
        :param aa: an array of coefficients.  I reverse the array if it is for a polynomial.
        :param exp: a boolean to indicate if I'm working with an exponential fit or not
        :return: a string equation, a TeX equation, and a label
        """
        a=np.flip(aa) if not exp else aa
        strEq=""
        strTex="y = "
        lbl="const"
        if len(a)==2: lbl="linear"
        if len(a)==3 and exp == False: lbl="quadratic"
        if len(a)==3 and exp == True: lbl="exponential"
        if len(a)==4: lbl="cubic"

        if exp==False: #HW10 addition: do this if not exp.  Also build TeX string for equation.

            for i in range(len(a)):
                if i==0:
                    strEq+="{:0.4f}".format(a[i])
                    strTex+="%s"%("{:0.3f}".format(a[i]))
                else:
                    if a[i]>=0.0:
                        strEq+="+"
                        strTex+="+"
                    if i>1:
                        strEq+="{:0.4f}*x^{:0d}".format(a[i],i)
                        strTex+="$%s{\cdot}x^%s$"%("{:0.4f}".format(a[i]),"{:0d}".format(i))
                    else:
                        strEq+="{:0.4f}*x".format(a[i])
                        strTex += "$%s{\cdot}x$" % ("{:0.4f}".format(a[i]))
        else: #HW10 addition:  do this if exp = True
            # a+b*exp(c*x)
            strEq="{:0.4f}".format(a[0])
            strEq += "+" if a[1]>=0 else ""
            strEq += "{:0.4f}*exp({:0.4f}*x)".format(a[1],a[2])
            strTex="{:0.4f}".format(a[0])+"+" if a[1]>0 else ""
            strTex+=r"%s${\cdot}e^{%s{\cdot}x}$"%("{:0.4f}".format(a[1]),"{:0.4f}".format(a[2]))
        return strEq, strTex, lbl

    def GetPlotInfo(self, power=1, npoints=500, exp=False):
        """
        This function generates x,y data for plotting the fitting curve for poly or exponential fit.
        :param power: the power of a polynomial
        :param npoints: number of x,y points to generate
        :param exp: boolean True if exponential fit, False if polynomial fit.
        :return: xvalues array, yvalues array, R squared, string equation, TeX equation
        """
        Xmin = min(self.x)
        Xmax = max(self.x)
        Ymin = min(self.y)
        Ymax = max(self.y)
        dX = 1.0 * (Xmax - Xmin) / npoints

        if exp: #HW10 addition
            a=self.OptimizeFit()
        else:
            a = self.LeastSquares(power)
        eq, TexEq, lbl = self.GetEq(a, exp=exp)

        xvals = []
        yvals = []
        for i in range(npoints):
            xvals.append(Xmin + i * dX)
            yvals.append(self.Exp(xvals[i],a) if exp else self.Poly(xvals[i], a))
        RSq = self.RSquared(a, exp)
        return xvals,yvals,RSq,eq,TexEq, lbl

    def PlotLeastSquares(self, power, axes=None, showpoints=True, npoints=500, exp=False):
        """
        This creates the plots for the polynomial or exponential fit
        :param power: power of the polynomial fit
        :param axes: matplotlib axes object
        :param showpoints: show the self.x, self.y data points (True)
        :param npoints: how many points to plot for the fit
        :param exp: is it an exponential fit? (True) default is poly (False)
        :return: xvals, yvals, RSq, eq, TexEq, lbl
        """
        QTPlotting=True
        if axes == None:
            axes=pyplot.subplot()
            QTPlotting=False

        xvals,yvals,RSq,eq,TexEq, lbl=self.GetPlotInfo(power=power,npoints=npoints, exp=exp)

        axes.plot(xvals,yvals,linestyle='dashed',color='black',linewidth='2',label=""+lbl+"{"+TexEq+"}"+r' ($R^2={:0.3f}$)'.format(RSq))
        if showpoints: axes.plot(self.x,self.y,linestyle='none',  marker='o', markerfacecolor='white', markeredgecolor='black', markersize=10)
        axes.grid(axis='both')
        axes.tick_params(labelsize=18,axis='both', direction='in', grid_linewidth=1, grid_linestyle='dashed', grid_alpha=0.5)
        #axes.set(xlabel='X values',ylabel='Y values')
        axes.set_ylabel('Y values', fontsize=18) #HW10 addition
        axes.set_xlabel('X values', fontsize=18) #HW10 addition
        axes.legend(fontsize=18,loc='upper right') #HW10 addition

        title=""
        if QTPlotting==False: pyplot.show()
        return xvals, yvals, RSq, eq, TexEq, lbl

def main():
    x = np.array([0.05, 0.11, 0.15, 0.31, 0.46, 0.52, 0.70, 0.74, 0.82, 0.98, 1.17])
    y = np.array([0.956, 1.09, 1.332, 0.717, 0.771, 0.539, 0.378, 0.370, 0.306, 0.242, 0.104])

    LS=trendlineClass(x, y)

    linx, liny, RSqLin, lineq, linTexEq, linlbl=LS.PlotLeastSquares(1,showpoints=True, npoints=500)

    cubx,cuby, RSqCub, cubeq, cubTexEq, cublbl=LS.PlotLeastSquares(3,showpoints=True, npoints=500)

    expx, expy, RSqExp, expeq, expTexEq, explbl = LS.PlotLeastSquares(power=3, exp=True)

    pyplot.plot(linx,liny, linewidth=2, linestyle='dashed', color='black', label=r'Linear fit ($R^2={:0.3f}$)'.format(RSqLin))
    pyplot.plot(cubx,cuby, linewidth=2, linestyle='dotted', color='black', label='Cubic fit ($R^2={:0.3f}$)'.format(RSqCub))
    pyplot.plot(expx,expy, linewidth=2, linestyle='solid', color='black', label='Exp fit ($R^2={:0.3f}$)'.format(RSqExp))
    pyplot.plot(x, y, linestyle='none', marker='o', markersize=10, markerfacecolor='white', markeredgecolor='black', label='Data')
    pyplot.xlabel('X values')
    pyplot.ylabel('Y values')
    pyplot.legend()

    pyplot.grid(axis='both')
    pyplot.tick_params(axis='both',direction='in', grid_linewidth=1, grid_linestyle='dashed',grid_alpha=0.5)
    pyplot.show()
    x = np.array([[10, 20, 30], [40, 50, 60]])
    y = np.array([[100], [200]])
    print(np.append(x, y, axis=1))
    #Use proper titles, labels and legends

if __name__=="__main__":
    main()