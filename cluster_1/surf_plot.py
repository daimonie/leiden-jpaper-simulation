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
parser.add_argument('-c', '--clip', help='Clip.', action='store', type = float, default = 40.0)   
parser.add_argument('-j', '--jnumber', help='What J to plot?', action='store', type = int, default = 3)   
parser.add_argument('-x', '--xlabel', help='Label for Horizontal axis?.', action='store', type = str, default = r"$\beta$")   
parser.add_argument('-y', '--ylabel', help='Label for Vertical axis?.', action='store', type = str, default = r"$J_1$" )   
parser.add_argument('--phi', help='first angle for view_init', action='store', type=int, default=90);
parser.add_argument('--theta', help='first angle for view_init', action='store', type=int, default=180);
parser.add_argument('--cols', help='Columns in data file, required for the order params', action='store', type=int, default=10);
args	= parser.parse_args() 


filename    = args.filename
mode        = args.mode
save        = args.save
clip_size   = args.clip
jnumber     = args.jnumber
xlabel      = args.xlabel
ylabel      = args.ylabel
phi         = args.phi
theta       = args.theta
cols        = args.cols
 

print "Plotting from file [%s], mode [%s], labels (%s, %s)" % (filename, mode, xlabel, ylabel)

file_handler = open( filename, "r" );

data = np.genfromtxt(file_handler, dtype=None, usecols=range(0,cols)); #excluding the symtype col

filename = filename.replace("_", " - ")
    

title = "Data file: %s" % filename

if (mode == "specific_heat"): 
    title = "Specific heat [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,2]
elif (mode == "energy"): 
    title = "Energy [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,1]
elif (mode == "order_one"): 
    title = "Order one [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,6] 
elif (mode == "chi_one"): 
    title = "Chi one [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,7] 
elif (mode == "order_two"): 
    title = "Order two [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,8] 
elif (mode == "chi_two"): 
    title = "Chi two [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,9] 
elif (mode == "order_three"): 
    title = "Order Three [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,10] 
elif (mode == "chi_three"): 
    title = "Chi Three [%s]" % filename
    xdata = data[:,0]
    ydata = data[:,jnumber]
    zdata = data[:,11] 

else:
    raise Exception("Incorrect mode.");

if (save == "temperature"):
    xlabel = r"$T$"
    xdata = np.reciprocal(xdata);
    
    
lin_x = np.linspace(min(xdata), max(xdata))
lin_y = np.linspace(min(ydata), max(ydata))

x, y = np.meshgrid(lin_x, lin_y)
z = griddata(xdata, ydata, zdata, lin_x, lin_y, interp='linear')

z = np.clip(z, 0, clip_size)


np.seterr('ignore')
fig = plt.figure(figsize=(20, 10))
ax = fig.gca(projection='3d') 
 
print "Setting (min,max) = (%.3f, %.3f) for colour scheme." % (z.min(), z.max()*0.80)
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.afmhot, linewidth=1, vmin=z.min(), vmax=z.max()*0.80)  


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.view_init(phi, theta)   

fig.colorbar(surf)
plt.xlabel( xlabel ,fontsize=30);
plt.ylabel( ylabel ,fontsize=30); 
plt.title( title ,fontsize=20);

if(save == "plot" or save == "temperature"):
    plt.show()
else: 
    raise Exception("saving doesn't support latex. Just do it manually, you slacker.");