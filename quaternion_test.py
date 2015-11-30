import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

num_bins = 20
N = 100000 #number of quaternions to generate
#planning.cs.uiuc.edu/node198.html
randoms = np.random.rand(3, N) 

quaternions = np.zeros((4,N))

quaternions[0,:] = np.sqrt(1-randoms[0,:]) * np.sin(2*np.pi*randoms[1,:])
quaternions[1,:] = np.sqrt(1-randoms[0,:]) * np.cos(2*np.pi*randoms[1,:])

quaternions[2,:] = np.sqrt( randoms[1,:]) * np.sin( 2 * np.pi * randoms[2,:])
quaternions[3,:] = np.sqrt( randoms[1,:]) * np.cos( 2 * np.pi * randoms[2,:])
  
for i in range(0,4):
        plt.subplot(2,2,i);
        plt.hist(quaternions[i,:], num_bins, facecolor='green', alpha=0.5) 
        plt.xlabel('random number')
        plt.ylabel('times rolled')
        plt.title('Histogram of $u^%d$' % i)



# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()