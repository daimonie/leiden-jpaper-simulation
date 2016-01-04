#Regular calculation modules.
import numpy as np 
import scipy as sp 
#Allows a debug-output stream.
import sys as sys 
#Physical constants list.
from scipy.constants import *
#Time differences.
import time as time  
#Command line arguments.
import argparse as argparse  

#griddata to format data
from matplotlib.mlab import griddata 

#3d surf plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

#Commandline arguments instruction.


parser	= argparse.ArgumentParser(prog="Surface Plot",
  description = "Surface plot of data file")  
parser.add_argument('-f', '--filename', help='Data file.', action='store', type = str)   
parser.add_argument('-m', '--mode', help='What do we want to plot.', action='store', type = str)   
parser.add_argument('-s', '--save', help='Save or show on screen?.', action='store', type = str, default = "plot")   
args	= parser.parse_args() 


filename    = args.filename
mode        = args.mode
save        = args.save

print "Plotting from file [%s], mode [%s] " % (filename, mode)

file_handler = open( filename, "r" );

data = np.genfromtxt(file_handler, dtype=None, usecols=range(0,10)); #excluding the symtype col

xlabel = "x"
ylabel = "y" 
title = "Data file: %s" % filename
if (mode == "specific_heat"):
    xlabel = "beta"
    ylabel = "J1/J3" 
    title = "Specific heat [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,3]
    zdata = data[:,2]
elif (mode == "energy"):
    xlabel = "beta"
    ylabel = "J1/J3" 
    title = "Energy [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,3]
    zdata = data[:,1]
elif (mode == "order_one"):
    xlabel = "beta"
    ylabel = "J1/J3" 
    title = "Order one [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,3]
    zdata = data[:,6]
elif (mode == "order_two"):
    xlabel = "beta"
    ylabel = "J1/J3" 
    title = "Order two [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,3]
    zdata = data[:,7]
elif (mode == "chi_one"):
    xlabel = "beta"
    ylabel = "J1/J3" 
    title = "Chi one [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,3]
    zdata = data[:,8]
elif (mode == "chi_two"):
    xlabel = "beta"
    ylabel = "J1/J3" 
    title = "Chi two [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,3]
    zdata = data[:,9]

else:
    raise Exception("Incorrect mode.");


lin_x = np.linspace(min(xdata), max(xdata))
lin_y = np.linspace(min(ydata), max(ydata))

x, y = np.meshgrid(lin_x, lin_y)
z = griddata(xdata, ydata, zdata, lin_x, lin_y, interp='linear')



fig = plt.figure()
ax = fig.gca(projection='3d') 

norm_y = y/y.max();

surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.summer, linewidth=1, vmin=z.min(), vmax=z.max()) 

ax.view_init(60, 160)  

fig.colorbar(surf)
plt.xlabel( xlabel );
plt.ylabel( ylabel ); 
plt.title( title );
if(save == "plot"):
    plt.show()
else:
    fig.savefig(filename)