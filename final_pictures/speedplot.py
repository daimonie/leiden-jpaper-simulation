
import numpy as np
import matplotlib.pyplot as plt  
from matplotlib import cm

def color (number):
    
    return cm.summer( number** 2 *1.0/ 49.)


x = np.arange(1, 11)*100
y1 = np.array([19083262,31901685,46221187,60233434,74500158,88506607,101903320,115750292,129435472,143059204]) / 1000000.
y2 = np.array([12214158,23549816,34530753,45572280,56736176,68704006,79301541,89948718,100284716,110765216]) / 1000000.
y3 = y2/2.
y4 = y2/3.
y5 = y2/4.
y6 = y2/48.
y7 = y2/64.

opac = 1.0

plt.fill_between(x, y1, 0*y1, facecolor=color(1),    alpha=opac, label='original')
plt.fill_between(x, y2, 0*y1, facecolor=color(2),      alpha=opac, label='1 core')
plt.fill_between(x, y3, 0*y1, facecolor=color(3),     alpha=opac, label='2 cores')
plt.fill_between(x, y4, 0*y1, facecolor=color(4),    alpha=opac, label='3 cores')
plt.fill_between(x, y5, 0*y1, facecolor=color(5),   alpha=opac, label='4 cores')  
plt.fill_between(x, y6, 0*y1, facecolor=color(6),   alpha=opac, label='48 cores') 
plt.fill_between(x, y7, 0*y1, facecolor=color(7),   alpha=opac, label='64 cores') 


plt.scatter(x, y1, label='original',    s=50, c=color(1))
plt.scatter(x, y2, label='1 core',      s=50, c=color(2))
plt.scatter(x, y3, label='2 cores',     s=50, c=color(3))
plt.scatter(x, y4, label='3 cores',     s=50, c=color(4))
plt.scatter(x, y5, label='4 cores',     s=50, c=color(5))  
plt.scatter(x, y6, label='48 cores',    s=50, c=color(6)) 
plt.scatter(x, y7, label='64 cores',    s=50, c=color(7)) 
plt.grid(True)
plt.legend()

plt.xlabel("Samples")
plt.ylabel("Time [s]")
plt.show()
