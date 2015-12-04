import numpy as np;
 
#U_x[511]: np.matrix([[-1.000,-1.000, 1.000], [1.000, -1.000, 1.000], [1.000, -1.000, 1.000]])
#R[511]: np.matrix([[0.800,0.353, -0.485], [0.571, -0.197, 0.797], [0.185, -0.915, -0.359]])
#R[510]: np.matrix([[-0.400,0.539, -0.741], [0.539, -0.516, -0.666], [-0.741, -0.666, -0.084]])
#comparison:
         #diag   0.395   0.000   0.000
        #Rfoo     0.395   1.505   0.073


0
uxi = np.matrix([[-1.000,0.000, 0.000], [0.000, 1.000, 0.000], [0.000, 0.000, -1.000]])
ri = np.matrix([[-0.000,-0.951, 0.309], [-0.588, 0.250, 0.769], [-0.809, -0.182, -0.559]])
rj = np.matrix([[0.400,-0.539, -0.741], [0.539, -0.516, 0.666], [-0.741, -0.666, 0.084]]) 
new = np.matrix([[0.741,-0.285, 0.607], [-0.671, -0.325, 0.667], [0.007, -0.902, -0.432]])
old = np.matrix([[0.283,0.697, 0.659], [0.940, -0.067, -0.334], [0.189, -0.714, 0.674]])



print old
print new
print "---"

final = uxi * ri * rj;
print final

print "\n"
#Now we need to do this manually
manual = np.zeros((3,3));
for i in range(0,3):
	for j in range (0,3):
		for l in range(0,3):
			for k in range(0,3):
				manual[i,j] += uxi[i, l] * ri[l,k] * rj[k,j]
				ii = i*3 + j
				#print "Added %.3f to %d" % (uxi[i, l] * ri[l,k] * rj[k,j], ii)
				
print manual
print "Formula seems correct.";