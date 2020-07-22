import numpy as np
import loaders
import matplotlib.pyplot as plt
import os
import argparse
import h5py
import glob
from array import *


#This script needs to:
# 1.) Loop over all run directories 
# 2.) Loop over all files in each directory
# 3.) Grab the data
# 4.) Calculate the (cumulative) precession 
# 5.) Filter out runs that crashed
# 6.) Plot e vs. q (MM08 and KD06) and precession vs. q


phase_total_vs_q = []
mm08_total_e_vs_q = []
q = []
directory_list = glob.glob('/Users/anshsehgal/Mara3/results/q*')
directory_list = sorted(directory_list)

for directory in directory_list:
	file_list = os.listdir(directory)
	file_list = sorted(file_list)
	mass = float(directory[33:37])
	q.append(mass/10.)
	counter = 0
	single_mm08_e_global = []
	single_mm08_phase_global =[]
	for file in file_list:
		filename = directory + '/' + file_list[counter]
		if len(file_list) < 50:
			#print("CRASH")
			break
		else:
			sigma_vs_r        = h5py.File(filename, 'r')['sigma_vs_r']
			mm08_vs_r         = h5py.File(filename, 'r')['mm08_vs_r']
			mm08_e_vs_r       = h5py.File(filename, 'r')['mm08_e_vs_r']
			mm08_phase_vs_r   = h5py.File(filename, 'r')['mm08_phase_vs_r']
			mm08_e_global     = (h5py.File(filename, 'r')['mm08_e_global'][...]).tolist()
			mm08_phase_global = (h5py.File(filename, 'r')['mm08_phase_global'][...]).tolist()
			kd06_e_vs_r       = h5py.File(filename, 'r')['kd06_e_vs_r']
			kd06_e_global     = h5py.File(filename, 'r')['kd06_e_global']
			#get_eccentricity_measures(filename)
			single_mm08_e_global.append(mm08_e_global)
			single_mm08_phase_global.append(mm08_phase_global)
			counter += 1
	avg_mm08_e_global = np.average(single_mm08_e_global)
	avg_mm08_phase_global = np.average(single_mm08_phase_global)
	mm08_total_e_vs_q.append(avg_mm08_e_global)
	phase_total_vs_q.append(avg_mm08_phase_global)

i = np.where(~np.isnan(mm08_total_e_vs_q))
print i 

#print mm08_total_e_vs_q

#plt.figure()
#plt.plot(q_plot,e_plot)
#plt.show()



