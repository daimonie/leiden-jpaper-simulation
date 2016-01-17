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
parser.add_argument('-f', '--filename', help='Data file.', action='store', type = str, default = 'exception')   
parser.add_argument('-m', '--mode', help='What do we want to plot.', action='store', type = str, default='xy')    
parser.add_argument('-c', '--clip', help='Clip.', action='store', type = float, default = 40.0)   

parser.add_argument('--xnumber', help='Column representing x-axis?', action='store', type = int, default = 0)   
parser.add_argument('--ynumber', help='Column representing y-axis?', action='store', type = int, default = 3)   
parser.add_argument('--znumber', help='Column representing z-axis?', action='store', type = int, default = 2)   
parser.add_argument('--wnumber', help='Column representing z-axis?', action='store', type = int, default = 4)   

parser.add_argument('-x', '--xlabel', help='Label for Horizontal axis?.', action='store', type = str, default = r"$\beta$")   
parser.add_argument('-y', '--ylabel', help='Label for Vertical axis?.', action='store', type = str, default = r"$J_1$" )   
parser.add_argument('-z', '--zlabel', help='Title?.', action='store', type = str, default = r"$J_1$ run" )   

parser.add_argument('--phi', help='first angle for view_init', action='store', type=int, default=90);
parser.add_argument('--theta', help='first angle for view_init', action='store', type=int, default=180);
parser.add_argument('--cols', help='Columns in data file, required for the order params', action='store', type=int, default=10);
args	= parser.parse_args() 


filename    = args.filename
mode        = args.mode 
clip_size   = args.clip

xnumber     = args.xnumber
ynumber     = args.ynumber
znumber     = args.znumber
wnumber     = args.wnumber

xlabel      = args.xlabel
ylabel      = args.ylabel
zlabel      = args.zlabel

phi         = args.phi
theta       = args.theta
cols        = args.cols

if filename == "exception":
    raise Exception("You have to pass a data file.")

print "Plotting from file [%s], mode [%s], labels (%s, %s)" % (filename, mode, xlabel, ylabel)

file_handler = open( filename, "r" );

data = np.genfromtxt(file_handler, dtype=None, usecols=range(0,cols)); #excluding the symtype col

filename = filename.replace("_", " - ")
    

title = "Data file: %s" % filename

if (mode == "collapse"): 
    title = "%s [%s]" % (zlabel, filename)
    xdata = data[:,xnumber]
    ydata = data[:,ynumber]
    zdata = (data[:,znumber] + data[:,wnumber])/2
elif (mode == "xy"): 
    title = "%s [%s]" % (zlabel, filename)
    xdata = data[:,xnumber]
    ydata = data[:,ynumber]
    zdata = data[:,znumber]
else:
    raise Exception("Incorrect mode.");
 
lin_x = np.linspace(min(xdata), max(xdata))
lin_y = np.linspace(min(ydata), max(ydata))

x, y = np.meshgrid(lin_x, lin_y)
z = griddata(xdata, ydata, zdata, lin_x, lin_y, interp='linear')

z = np.clip(z, 0, clip_size)


np.seterr('ignore')
fig = plt.figure(figsize=(20, 10))
ax = fig.gca(projection='3d') 

[xx,yy] = np.meshgrid(lin_x, lin_y);
zz = xx * 0;

#python -W ignore::FutureWarning surf_plot.py -f maris061_data_d2d_large.txt -m specific_heat -s plot -c 15

 
print "Setting (min,max) = (%.3f, %.3f) for colour scheme." % (z.min(), z.max()*0.80)
ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap=cm.afmhot, linewidth=1, vmin=z.min(), vmax=z.max()*0.80)  
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.afmhot, linewidth=1, vmin=z.min(), vmax=z.max()*0.80)  


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.view_init(phi, theta)   

fig.colorbar(surf)
plt.xlabel( xlabel ,fontsize=30);
plt.ylabel( ylabel ,fontsize=30); 
plt.title( title ,fontsize=20);
 
plt.show() 