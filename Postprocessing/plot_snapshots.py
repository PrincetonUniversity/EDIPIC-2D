import math
import struct
import numpy as np
import sys
import glob
import re
import os
from subprocess import call
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D




#np.printoptions(threshold=100)
def read_next_float(f):
	return struct.unpack('f', f.read(4))[0]

def read_float_arr(n, f):
	return [read_next_float(f) for i in range(0, n)]


#Prompting the user for the variable name
data_name = str(input("Enter the snapshot paramter name (e.g. Ne/EX,TXi):\n") or sys.exit())

#Prompting the user for the number of steps between plots
plot_per = int(input("Enter the number of steps between plots (default is 1):\n") or "1")

#Prompting the user for zmin
uzmin = input("Enter zmin for the surface plot (default is the minimum of each data set):\n") or False
uzmax = input("Enter zmax for the surface plot (default is the maximum of each data set):\n") or False

if uzmin != False:
	uzmin = float(uzmin)
	
if uzmax != False:
	uzmax = float(uzmax)	


print_dir = data_name + "_plots"


#sys.exit()
#Getting the number of files for this plasma parameter
file_list = glob.glob("*" + data_name + "*.bin")
nfiles = len(file_list)

#Getting the file numbers
tlist = [int(re.findall('\d+',file_list[i])[0]) for i in range(nfiles)]

#Ordering the list
file_list = np.array(file_list)[np.argsort(tlist)]
tlist = np.sort(tlist)



#Setting up the directory for plots, over-write if existing
if (os.path.isdir(print_dir)):
        #os.system('rm -r -I ' + print_image_dir)
	os.system('rm -r ' + print_dir)

call(['mkdir',print_dir])


#Plot settings
label_size = 30
tick_size = 18
cbar_size = 14



#Plotting the data
for i in range(0,nfiles,plot_per):
	
	#Getting the file name and printing it
	fname = str(file_list[i])
	print(fname)

	#Opening the file
	f = open(fname,'rb')

	#Getting the number of cells in the x-direction
	nx = int(read_next_float(f))

	#Re-setting the pointer to the head of the file
	f.seek(0)

	#Getting the number of nodes in the y-direction from the file size and nx
	ny = int(os.path.getsize(fname)/4/(nx+1)-1)

	#Reading the first line of the matrix
	DAT = read_float_arr(nx+1,f)

	#Stacking all other lines of the matrix
	for j in range(ny):
		
		DAT = np.vstack([DAT,read_float_arr(nx+1,f)])

	#Getting the X-grid
	X = DAT[0,1:]

	#Getting the Y-grid
	Y = DAT[1:,0]

	#Refining the data
	DAT = DAT[1:,1:]


	#Setting zmin and zmax if not set by the user	
	if uzmin == False:
		zmin = np.amin(np.amin(DAT))
	else:
		zmin = uzmin
		
	if uzmax == False:
		zmax = np.amax(np.amax(DAT))
	else:
		zmax = uzmax
	
	#Creating the vector of levels	
	zlev = np.linspace(zmin,zmax,40)

	#continue
	#Initializing the figure
	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	#plotting
	cs = plt.contourf(X, Y, DAT, zlev)

	plt.xlabel(r"$x$", fontsize=label_size)
	plt.ylabel(r"$y$", fontsize=label_size)

	plt.xlim([X[0],X[-1]])
	plt.ylim([Y[0],Y[-1]])
	
	plt.tick_params(labelsize=tick_size)

	cb = plt.colorbar(cs)
	cb.ax.tick_params(labelsize=cbar_size)

	#Shifting axis slightly
	pos1 = ax.get_position()
	pos2 = [pos1.x0+0.02,pos1.y0+0.02,pos1.width,pos1.height]
	ax.set_position(pos2)

	pos1 = cb.ax.get_position()
	pos2 = [pos1.x0+0.02,pos1.y0+0.02,pos1.width,pos1.height]
	cb.ax.set_position(pos2)


	plt.savefig(print_dir + "/%d.png" % tlist[i])

	#plt.show()

	plt.close(fig)






