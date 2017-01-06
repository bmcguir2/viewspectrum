#!/usr/bin/env python

#############################################################
#						Revision History					#
#############################################################

# reads in an spcat catalog and displays a spectrum

# 1.0 - Project start
# 1.1 - Include temperature variation
# 1.2 - Include option to write spectrum
# 1.3 - Include option to simulate Gaussians
# 1.4 - Include ability to underplot actual spectra
# 1.5 - streamlining gaussian tomfoolery
# 1.6 - enabling frequency limits
# 1.7 - adding ability to apply a VLSR offset to the simulated spectrum
# 1.8 - adding ability to write a simulation out from interactive session
# 2.0 - dynamically-updating plots
# 3.0 - store and plot multiple species, switches to requiring ipython
# 3.1 - restore from save file
# 3.2 - add filter by error limit

#############################################################
#							Preamble						#
#############################################################

import os, sys, argparse, math

#Python version check

if sys.version_info.major != 3:

	print("This code is written in Python 3.  It will not execute in Python 2.7.  Exiting (sorry).")
	
	quit()

import numpy as np
from numpy import exp as exp
import time as tm
import random
import warnings
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import itertools
import matplotlib.lines as mlines
from datetime import datetime, date, time
#warnings.filterwarnings('error')

version = 3.2

h = 6.626*10**(-34) #Planck's constant in J*s
k = 1.381*10**(-23) #Boltzmann's constant in J/K
kcm = 0.69503476 #Boltzmann's constant in cm-1/K
ckm = 2.998*10**5 #speed of light in km/s

#############################################################
#							Defaults 						#
#############################################################

first_run = True

auto_update = False

error_limit = float('inf')

thermal = float('inf') #initial default cutoff for optically-thick lines (i.e. don't touch them unless thermal is modified.)

T = 300 #temperature for simulations.  Default is 300 K.

catalog_file = None #catalog file to load in.  Needs to be a string.

vlsr = 0.0 #vlsr offset applied to simulation.  Default is 0 km/s.

S = 1 #linear scaling factor applied to simulation.  Default is 1.

ll = float('-inf') #lower limit for the simulation range.  Default is none.

ul = float('inf') #upper limit for the simulation range.  Default is none.

spec = None	#file of a laboratory or observational spectrum to load in for comparison.  Default is none.

dV = 5.0 #linewidth of the simulation.  Default is 5.0 km/s.
	
CT = 300.0 #temperature the catalog is simulated at.  Default is 300 K.

gauss = True #toggle for simulating Gaussians or a stick spectrum.  Default is True.

#############################################################
#							Functions						#
#############################################################
	
	
#read_cat reads the catalog file in

def read_cat(catalog_file):

	'''
	Reads in a catalog file line by line
	'''

	my_array = []

	try:
		with open(catalog_file) as input:
	
			for line in input:
		
				my_array.append(line)	
	except TypeError:
		print('Specify a catalog file with catalog_file = \'x\'')
		return			
			
	return my_array	
	
#trim_raw_array takes the raw_array and trims it at the specified ll and ul

def trim_raw_array(x):

	'''
	takes the raw array and trims it at the specified ll and ul
	'''

	tmp_array = np.arange(len(x),dtype=np.float)	#create a temporary array to splice out the frequency
	
	for line in range(len(x)):
	
		tmp_array[line] = float(str(x[line][:13]).strip())
		
	
	try:
		i = np.where(tmp_array > ll)[0][0]	#get the index of the first value above the lower limit
	except IndexError:
		i = 0								#if the catalog begins after the lower limit
	try:
		i2 = np.where(tmp_array > ul)[0][0]	#get the index of the first value above the upper limit
	except IndexError:
		i2 = len(tmp_array)					#if the catalog ends before the upper limit is reached
		
	trimmed_array = []

	for y in range(i,i2):
	
		trimmed_array.append(x[y])
		
	return trimmed_array				
	
#fix_qn fixes quantum number issues

def fix_qn(qnarray,line,old_qn):

	'''
	fixes quantum number issues arising from the use of alphabet characters to represent numbers in spcat
	'''

	new_qn = 000
			
	if 'A' in old_qn:
		
		new_qn = 100 + int(old_qn[1])
		
	if 'B' in old_qn:
		
		new_qn = 110 + int(old_qn[1])	
		
	if 'C' in old_qn:
		
		new_qn = 120 + int(old_qn[1])		

	if 'D' in old_qn:
		
		new_qn = 130 + int(old_qn[1])
		
	if 'E' in old_qn:
		
		new_qn = 140 + int(old_qn[1])
		
	if 'F' in old_qn:
		
		new_qn = 150 + int(old_qn[1])
		
	if 'G' in old_qn:
		
		new_qn = 160 + int(old_qn[1])
		
	if 'H' in old_qn:
		
		new_qn = 170 + int(old_qn[1])				
		
	if 'I' in old_qn:
		
		new_qn = 180 + int(old_qn[1])	
		
	if 'J' in old_qn:
		
		new_qn = 190 + int(old_qn[1])
		
	if 'K' in old_qn:
		
		new_qn = 200 + int(old_qn[1])
		
	if 'L' in old_qn:
		
		new_qn = 210 + int(old_qn[1])
		
	if 'M' in old_qn:
		
		new_qn = 220 + int(old_qn[1])	
		
	if 'N' in old_qn:
		
		new_qn = 230 + int(old_qn[1])	
		
	if 'O' in old_qn:
		
		new_qn = 240 + int(old_qn[1])
		
	if 'P' in old_qn:
		
		new_qn = 250 + int(old_qn[1])
		
	if 'Q' in old_qn:
		
		new_qn = 260 + int(old_qn[1])	
		
	if 'R' in old_qn:
		
		new_qn = 270 + int(old_qn[1])
		
	if 'S' in old_qn:
		
		new_qn = 280 + int(old_qn[1])
		
	if 'T' in old_qn:
		
		new_qn = 290 + int(old_qn[1])	
		
	if 'U' in old_qn:
		
		new_qn = 300 + int(old_qn[1])	
		
	if 'V' in old_qn:
		
		new_qn = 310 + int(old_qn[1])
		
	if 'W' in old_qn:
		
		new_qn = 320 + int(old_qn[1])	
		
	if 'X' in old_qn:
		
		new_qn = 330 + int(old_qn[1])	
		
	if 'Y' in old_qn:
		
		new_qn = 340 + int(old_qn[1])	
		
	if 'Z' in old_qn:
		
		new_qn = 350 + int(old_qn[1])																																											
				
	qnarray[line] = int(new_qn)			
	
# splices the catalog file appropriately, then populates a numpy array with the data

def splice_array(x):

	'''
	splices the catalog file appropriately, then populates a numpy array with the data
	'''

	frequency = np.arange(len(x),dtype=np.float)
	error = np.arange(len(x),dtype=np.float)
	logint = np.arange(len(x),dtype=np.float)
	dof = np.arange(len(x),dtype=np.int)
	elower = np.arange(len(x),dtype=np.float)
	gup = np.arange(len(x),dtype=np.int)
	tag = np.arange(len(x),dtype=np.int)
	qnformat = np.arange(len(x),dtype=np.int)
	qn1 = np.arange(len(x),dtype=object)
	qn2 = np.empty(len(x),dtype=object)
	qn3 = np.empty(len(x),dtype=object)
	qn4 = np.empty(len(x),dtype=object)
	qn5 = np.empty(len(x),dtype=object)
	qn6 = np.empty(len(x),dtype=object)
	qn7 = np.empty(len(x),dtype=object)
	qn8 = np.empty(len(x),dtype=object)
	qn9 = np.empty(len(x),dtype=object)
	qn10 = np.empty(len(x),dtype=object)
	qn11 = np.empty(len(x),dtype=object)
	qn12 = np.empty(len(x),dtype=object)

	for line in range(len(x)):
	
		frequency[line] = float(str(x[line][:13]).strip())
		error[line] = float(str(x[line][13:21]).strip())
		logint[line] = float(str(x[line][21:29]).strip())
		dof[line] = int(str(x[line][29:31]).strip())
		elower[line] = float(str(x[line][31:41]).strip())
		try:
			gup[line] = int(str(x[line][41:44]).strip()) if str(x[line][41:44]).strip() else ''
		except ValueError:
			fix_qn(gup,line,str(x[line][41:44]))
		tag[line] = int(str(x[line][44:51]).strip())
		qnformat[line] = int(str(x[line][51:55]).strip())
		try:
			qn1[line] = int(x[line][55:57]) if str(x[line][54:57]).strip() else ''
		except ValueError:
			fix_qn(qn1,line,str(x[line][55:57]))	
		try:		
			qn2[line] = int(x[line][57:59]) if str(x[line][57:59]).strip() else ''
		except ValueError:
			fix_qn(qn2,line,str(x[line][57:59]))
		try:				
			qn3[line] = int(x[line][59:61]) if str(x[line][59:61]).strip() else ''
		except ValueError:
			fix_qn(qn3,line,str(x[line][59:61]))
		try:
			qn4[line] = int(x[line][61:63]) if str(x[line][61:63]).strip() else ''
		except ValueError:
			fix_qn(qn4,line,str(x[line][61:63]))
		try:
			qn5[line] = int(x[line][63:65]) if str(x[line][63:65]).strip() else ''
		except ValueError:
			fix_qn(qn5,line,str(x[line][63:65]))
		try:
			qn6[line] = int(x[line][65:67]) if str(x[line][65:67]).strip() else ''
		except ValueError:
			fix_qn(qn6,line,str(x[line][65:67]))
		try:
			qn7[line] = int(x[line][67:69]) if str(x[line][67:69]).strip() else ''
		except ValueError:
			fix_qn(qn7,line,str(x[line][67:69]))
		try:
			qn8[line] = int(x[line][69:71]) if str(x[line][69:71]).strip() else ''
		except ValueError:
			fix_qn(qn8,line,str(x[line][69:71]))
		try:
			qn9[line] = int(x[line][71:73]) if str(x[line][71:73]).strip() else ''
		except ValueError:
			fix_qn(qn9,line,str(x[line][71:73]))
		try:
			qn10[line] = int(x[line][73:75]) if str(x[line][73:75]).strip() else ''
		except ValueError:
			fix_qn(qn10,line,str(x[line][73:75]))
		try:
			qn11[line] = int(x[line][75:77]) if str(x[line][75:77]).strip() else ''
		except ValueError:
			fix_qn(qn11,line,str(x[line][75:77]))
		try:
			qn12[line] = int(x[line][77:]) if str(x[line][77:]).strip() else ''
		except ValueError:
			fix_qn(qn12,line,str(x[line][77:]))		
			

	return frequency,error,logint,dof,elower,gup,tag,qnformat,qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8,qn9,qn10,qn11,qn12
	
#det_qns determines how many qns represent each state

def det_qns(qn7,qn8,qn9,qn10,qn11,qn12):

	'''
	determines how many qns represent each state
	'''
	
	qns = 1
	
	try:
		int(qn8[0])
		qns = 2
	except ValueError:
		pass
	except IndexError:
		print('There are no lines in the chosen range.  The program is about to crash, sorry.  I will fix this at some point.')
		print('Catalog file is: {}.' .format(catalog_file))
		
	try:
		int(qn9[0])
		qns = 3
	except ValueError:
		pass	
		
	try:
		int(qn10[0])
		qns = 4
	except ValueError:
		pass	

	try:
		int(qn11[0])
		qns = 5
	except ValueError:
		pass
		
	try:
		int(qn12[0])
		qns = 6
	except ValueError:
		pass		

	return qns

#calc_q will dynamically calculate a partition function whenever needed at a given T.  The catalog file used must have enough lines in it to fully capture the partition function, or the result will not be accurate for Q.
	
def calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,catalog_file):

	'''
	Dynamically calculates a partition function whenever needed at a given T.  The catalog file used must have enough lines in it to fully capture the partition function, or the result will not be accurate for Q.  This is perfectly fine for the *relative* intensities of lines for a given molecule used by this program.  However, absolute intensities between molecules are not remotely accurate.  
	'''

	Q = np.float64(0.0) #Initialize a float for the partition function
	
	if catalog_file=='acetone.cat':
	
		Q = 2.91296*10**(-7)*T**6 - 0.00021050085*T**5 + 0.05471337*T**4 - 5.5477*T**3 + 245.28*T**2 - 2728.3*T + 16431 #Hard code for Acetone
		
	elif catalog_file=='sh.cat':
	
		Q = 0.000000012549467*T**4 - 0.000008528126823*T**3 + 0.002288160909445*T**2 + 0.069272946237033*T + 15.357239728157400
 #Hard code for SH.  Completely unreliable below 2.735 K or above 300 K.

	elif catalog_file=='h2s.cat':
	
		Q = -0.000004859941547*T**3 + 0.005498622332982*T**2 + 0.507648423477309*T - 1.764494755639740 #Hard code for H2S.  Completely unreliable below 2.735 K or above 300 K.
	
	elif catalog_file=='hcn.cat':
	
		Q = -1.64946939*10**-9*T**4 + 4.62476813*10**-6*T**3 - 1.15188755*10**-3*T**2 + 1.48629408*T + .386550361
		
	elif catalog_file=='hc9n.cat':
	
		Q = 71.730808*T + 0.040659
		
	elif catalog_file=='hc7n.cat':
	
		Q = 36.949992*T + 0.135605
		
	elif catalog_file.lower()=='methanol.cat' or catalog_file.lower()=='ch3oh.cat':
	
		Q = 4.83410*10**-11*T**6 - 4.04024*10**-8*T**5 + 1.27624*10**-5*T**4 - 1.83807*10**-3*T**3 + 2.05911*10**-1*T**2 + 4.39632*10**-1*T -1.25670
		
	elif catalog_file.lower()=='13methanol.cat' or catalog_file.lower()=='13ch3oh.cat':
	
		Q = 0.000050130*T**3 + 0.076540934*T**2 + 4.317920731*T - 31.876881967
		
	elif catalog_file.lower()=='c2n.cat' or catalog_file.lower()=='ccn.cat':
		
		Q = 1.173755*10**(-11)*T**6 - 1.324086*10**(-8)*T**5 + 5.99936*10**(-6)*T**4 - 1.40473*10**(-3)*T**3 + 0.1837397*T**2 + 7.135161*T + 22.55770
	
	else:
	
		nstates = elower.size #Determine the number of total states in the raw cat file
	
		combined_array = np.empty(shape=(nstates,qns+1)) #Set up an array that has [J, ka, kc, Elower]

		if (qns == 1):
	
			for i in range(nstates): #Fill that array with [J, ka, kc, Elower]
		
				combined_array[i][0] = qn7[i]
				combined_array[i][1] = elower[i]

		if (qns == 2):
	
			for i in range(nstates): #Fill that array with [J, ka, kc, Elower]
		
				combined_array[i][0] = qn7[i]
				combined_array[i][1] = qn8[i]
				combined_array[i][2] = elower[i]
	
		if (qns == 3):
	
			for i in range(nstates): #Fill that array with [J, ka, kc, Elower]
		
				combined_array[i][0] = qn7[i]
				combined_array[i][1] = qn8[i]
				combined_array[i][2] = qn9[i]
				combined_array[i][3] = elower[i]
			
		if (qns == 4):
	
			for i in range(nstates): #Fill that array with [J, ka, kc, QN10, Elower]
		
				combined_array[i][0] = qn7[i]
				combined_array[i][1] = qn8[i]
				combined_array[i][2] = qn9[i]
				combined_array[i][3] = qn10[i]
				combined_array[i][4] = elower[i]			

		if (qns == 5):
	
			for i in range(nstates): #Fill that array with [J, ka, kc, QN10, QN11, Elower]
		
				combined_array[i][0] = qn7[i]
				combined_array[i][1] = qn8[i]
				combined_array[i][2] = qn9[i]
				combined_array[i][3] = qn10[i]
				combined_array[i][4] = qn11[i]
				combined_array[i][5] = elower[i]	
			
		if (qns == 6):
	
			for i in range(nstates): #Fill that array with [J, ka, kc, QN10, QN11, QN12, Elower]
		
				try:
					combined_array[i][0] = qn7[i]
				except ValueError:
					print('I choked at index {}.' .format(i))
					quit()
				combined_array[i][1] = qn8[i]
				combined_array[i][2] = qn9[i]
				combined_array[i][3] = qn10[i]
				combined_array[i][4] = qn11[i]
				combined_array[i][5] = qn12[i]
				combined_array[i][6] = elower[i]									
		
		temp = list(set(map(tuple,combined_array))) #Don't know HOW this works, but it does: sorts through the array and removes all duplicate entries, so that only a single entry remains for each set of quantum numbers.
	
		ustates = len(temp) #Number of unique lower states
	
		e_index = qns
	
		for i in range(ustates):
	
			J = temp[i][0] #Extract a J value from the list
			E = temp[i][qns] #Extract its corresponding energy
		
			Q += (2*J+1)*exp(np.float64(-E/(kcm*T))) #Add it to Q
			
		#result = [Q,ustates] #can enable this function to check the number of states used in the calculation, but note that this will break calls to Q further down that don't go after element 0.
	
	return Q

#scale_temp scales the simulated intensities to the new temperature

def scale_temp(int_sim,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT,catalog_file):

	'''
	Converts linear intensities at one temperature to another.
	'''

	scaled_int = np.empty_like(int_sim)
	
	Q_T = calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,catalog_file)
	Q_CT = calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,CT,catalog_file)
	
	scaled_int = int_sim * (Q_CT/Q_T) * exp(-(((1/T)-(1/CT))*elower)/0.695)
	
# 	for i in range(len(scaled_int)):
# 	
# 		if catalog[1][i] > error_limit:
# 		
# 			scaled_int[i] = 0.0

	return scaled_int

#convert_int converts the intensity to not log units

def convert_int(logint):

	'''
	Converts catalog logarithmic intensity units to linear ones
	'''

	intensity = np.empty_like(logint)	
	
	intensity = 10**(logint)
	
	return intensity
	
#simulates Gaussian profiles after intensities are simulated.

def sim_gaussian(int_sim,freq,linewidth):

	'''
	Simulates Gaussian profiles for lines, after the intensities have been calculated.  Tries to be smart with how large a range it simulates over, for computational resources.  Includes a thermal cutoff for optically-thick lines.
	'''
	
	freq_gauss = []

	for x in range(int_sim.shape[0]):
	
		l_f = linewidth*freq[x]/ckm #get the FWHM in MHz
	
		min_f = freq[x] - 10*l_f #get the frequency 10 FWHM lower
		max_f = freq[x] + 10*l_f #get the frequency 10 FWHM higher
	
		res_pnts = l_f / 15 #determine a resolution element (15 points across the line)
	
		freq_line = np.arange(min_f,max_f,res_pnts)
	
		freq_gauss.extend(freq_line)
	
	int_gauss = [0] * len(freq_gauss)
	
	freq_gauss.sort()
	
	start_time = tm.time()
	
	alerted = False

	for x in range(int_sim.shape[0]):
	
		telapsed = tm.time() - start_time
		
		if telapsed > 5 and alerted == False:
		
			tstep = telapsed/x
			
			ttotal = (tstep * int_sim.shape[0])/60
		
			print('You have asked for a computationally-expensive simulation.  Either wait for it to finish, or narrow up your frequency range by setting ll or ul.')
			print('A (probably very poor) estimation of the time remaining is: {:.2f} minutes.' .format(ttotal))
			
			alerted = True 
	
		l_f = linewidth*freq[x]/ckm #get the FWHM in MHz
	
		c = l_f/2.35482

		int_gauss += int_sim[x]*exp(-((freq_gauss - freq[x])**2/(2*c**2)))
		
	int_gauss[int_gauss > thermal] = thermal	
	
	return(freq_gauss,int_gauss)

#write_spectrum writes out the current freq and intensity to output_file

def write_spectrum(x,output_file):

	'''
	Will write out the simulation frequency and intensity for the currently active simulation ('current'), the summed spectrum from sum_stored() ('sum'), or any stored simulation 'x' to output_file (which must be given as a string).  
	'''
	
	if x == 'current':
	
		freq_tmp = freq_sim
		int_tmp = int_sim
		
	elif x == 'sum':
	
		freq_tmp = freq_sum
		int_tmp = int_sum
		
	else:
	
		freq_tmp = sim[x].freq_sim
		int_tmp = sim[x].int_sim

	if gauss == True:
	
		for h in range(len(freq_tmp)):

			if (h == 0): #write out the results to a out_file
			
				with open(output_file, 'w') as output: 
					
					output.write('{} {}\n' .format(freq_tmp[0],int_tmp[0]))
						
			else:
				
				with open(output_file, 'a') as output:
					
					output.write('{} {}\n' .format(freq_tmp[h],int_tmp[h]))
						
			h += 1		
	
	else:
	
		for h in range(freq_tmp.shape[0]):

			if (h == 0): #write out the results to a out_file
				
				with open(output_file, 'w') as output: 
					
					output.write('{} {}\n' .format(freq_tmp[0],int_tmp[0]))
						
			else:
				
				with open(output_file, 'a') as output:
					
					output.write('{} {}\n' .format(freq_tmp[h],int_tmp[h]))
						
			h += 1				

#run_sim runs the simulation.  It's a meta routine, so that we can update later

def run_sim(freq,intensity,T,dV,S):

	'''
	Runs a full simulation accounting for the currently-active T, dV, S, and vlsr values, as well as any thermal cutoff for optically-thick lines
	'''

	int_temp = scale_temp(intensity,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT,catalog_file)

	int_temp *= S
	
	freq_tmp = np.copy(freq)
	
	int_temp[int_temp > thermal] = thermal	
	
	if gauss == True:

		freq_sim,int_sim = sim_gaussian(int_temp,freq_tmp,dV)
		
	else:
	
		#int_temp[int_temp > thermal] = thermal
		freq_sim = freq_tmp
		int_sim = int_temp
		
	return freq_sim,int_sim
	
#mod_T changes the temperature, re-simulates, and re-plots	
	
def modT(x):

	'''
	Modifies the temperature value of the current simulation to the value given
	'''

	try:
		float(x)
	except:
		print('x needs to be a number.')
		return
		
	global T,freq_sim,int_sim
		
	T = float(x)
		
	freq_tmp = np.copy(frequency)
	
	freq_tmp += (-vlsr)*freq_tmp/ckm		
	
	freq_sim,int_sim = run_sim(freq_tmp,intensity,T,dV,S)
	
	clear_line('current')
	
	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)
		
	ax.legend()
	fig.canvas.draw()
	
	save_results('last.results')
	
#modS changes the scaling, re-simulates, and re-plots

def modS(x):

	'''
	Modifies the scaling value of the current simulation to the value given
	'''

	try:
		float(x)
	except:
		print('x needs to be a number.')
		return
		
	global S,freq_sim,int_sim
		
	S = x
		
	freq_tmp = np.copy(frequency)
	
	freq_tmp += (-vlsr)*freq_tmp/ckm		
	
	freq_sim,int_sim = run_sim(freq_tmp,intensity,T,dV,S)
	
	clear_line('current')
	
	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)
		
	ax.legend()
	fig.canvas.draw()
	
	save_results('last.results')
	
#moddV changes the velocity width, re-simulates, and re-plots

def moddV(x):

	'''
	Modifies the dV value of the current simulation to the value given
	'''

	try:
		float(x)
	except:
		print('dV needs to be a number.')
		return
		
	global dV,freq_sim,int_sim
	
	dV = x
	
	freq_tmp = np.copy(frequency)
	
	freq_tmp += (-vlsr)*freq_tmp/ckm		
		
	freq_sim,int_sim = run_sim(freq_tmp,intensity,T,dV,S)
	
	clear_line('current')
	
	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)
		
	ax.legend()
	fig.canvas.draw()
	
	save_results('last.results')
	
#modVLSR changes the LSR velocity, re-simulates, and re-plots

def modVLSR(x):

	'''
	Modifies the vlsr value of the current simulation to the value given
	'''

	try:
		float(x)
	except:
		print('vlsr needs to be a number.')
		return	
		
	global vlsr,freq_sim,int_sim,frequency
	
	vlsr = x
		
	freq_tmp = np.copy(frequency)
	
	freq_tmp += (-vlsr)*freq_tmp/ckm		
	
	freq_sim,int_sim = run_sim(freq_tmp,intensity,T,dV,S)
	
	clear_line('current')
	
	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)
	
	ax.legend()
	fig.canvas.draw()
	
	save_results('last.results')		

#modV is an alias for modVLSR

def modV(vlsr):

	'''
	Modifies the vlsr value of the current simulation to the value given
	'''

	modVLSR(vlsr)	
	
#make_plot makes the plot!

def make_plot():

	'''
	Generates a plot of the currently-active molecular simulation, as well as any laboratory data or observations which are loaded.  This will *not* restore any overplots from a previously-closed plot.
	'''

	global fig,ax,line1,line2,freq_sim,intensity,int_sim

	plt.ion()	

	fig = plt.figure()
	ax = fig.add_subplot(111)

	minorLocator = AutoMinorLocator(5)
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Intensity (Probably Arbitrary)')

	plt.locator_params(nbins=4) #Use only 4 actual numbers on the x-axis
	ax.xaxis.set_minor_locator(minorLocator) #Let the program calculate some minor ticks from that

	ax.get_xaxis().get_major_formatter().set_scientific(False) #Don't let the x-axis go into scientific notation
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	
	try:
		freq_obs[0]
		lines['obs'] = ax.plot(freq_obs,int_obs,color = 'black',label='obs',zorder=0)
	except:
		pass

	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)	

	ax.legend()
	fig.canvas.draw()
	
	save_results('last.results')
	
#obs_off turns off the observations
	
def obs_off():

	'''
	turns off the laboratory data or observations on the plot
	'''

	try:
		clear_line('obs')
		save_results('last.results')
		return
	except:
		print('The observations are already off.  You can turn them on with obs_on()')
		return
	
	try:
		lines['obs']
		save_results('last.results')
	except:
		print('There are no observations loaded into the program to turn off.  Load in obs with read_obs()')
		return	
			
#obs_on turns on the observations
	
def obs_on():

	'''
	turns on the laboratory data or observations on the plot
	'''	
	
	try:
		lines['obs']
		return
	except:
		pass

	try:
		lines['obs'] = 	ax.plot(freq_obs,int_obs,color = 'black',label='obs',zorder=0)
		ax.legend()
		fig.canvas.draw()
		save_results('last.results')
	except:
		print('There are no observations loaded into the program to turn on.  Load in obs with read_obs()')
		return
			
#read_obs reads in observations or laboratory spectra and populates freq_obs and int_obs

def read_obs(x):

	'''
	reads in observations or laboratory spectra and populates freq_obs and int_obs
	'''

	global spec

	spec = x

	obs = read_cat(x)
	
	global freq_obs,int_obs

	freq_obs = []
	int_obs = []

	for x in range(len(obs)):

		freq_obs.append(float(obs[x].split()[0]))
		int_obs.append(float(obs[x].split()[1].strip('\n')))
		
	try:		
		lines['obs'] = 	ax.plot(freq_obs,int_obs,color = 'black',label='obs',zorder=0)
	except:
		return
		
#close closes the currently open plot

def close():

	'''
	Closes the currently open plot window.
	'''
	
	plt.close()	
	
	save_results('last.results')
	
#store saves the current simulation parameters for recall later.  *Not* saved as a Gaussian. 'x' must be entered as a string with quotes.

def store(x):

	'''
	saves the current simulation parameters for recall later.  *Not* saved as a Gaussian. 'x' must be entered as a string with quotes.
	'''
	
	sim[x] = Molecule(x,catalog_file,elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,CT,vlsr,frequency,freq_sim,intensity,int_sim)
	
	if auto_update == True:
	
		sum_stored()
		overplot_sum()
		
	if any(line for line in ax.lines if line.get_label()==x):
		overplot(x)	
	
	save_results('last.results') 
	
#recall wipes the current simulation and re-loads a previous simulation that was stored with store(). 'x' must be entered as a string with quotes. This will close the currently-open plot.

def recall(x):

	'''
	wipes the current simulation and re-loads a previous simulation that was stored with store(). 'x' must be entered as a string with quotes. This will close the currently-open plot.
	'''

	save_results('last.results')

	global elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,vlsr,frequency,freq_sim,intensity,int_sim,current,catalog_file

	current = sim[x].name
	elower = sim[x].elower
	eupper = sim[x].eupper
	qns = sim[x].qns
	logint = sim[x].logint
	qn7 = sim[x].qn7
	qn8 = sim[x].qn8
	qn9 = sim[x].qn9
	qn10 = sim[x].qn10
	qn11 = sim[x].qn11
	qn12 = sim[x].qn12
	S = sim[x].S
	dV = sim[x].dV
	T = sim[x].T
	vlsr = sim[x].vlsr
	frequency = sim[x].frequency
	freq_sim = sim[x].freq_sim
	intensity = sim[x].intensity
	int_sim = sim[x].int_sim
	catalog_file = sim[x].catalog_file

	try:
		clear_line('current')
	except:
		pass
		
	tmp_freq = np.copy(frequency)
	
	tmp_freq += (-vlsr)*tmp_freq/ckm
	
	freq_sim,int_sim=run_sim(tmp_freq,intensity,T,dV,S)	
		
	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)	
		
	try:
		plt.get_fignums()[0]	
		ax.legend()
		fig.canvas.draw()
	except:	
		make_plot()	
		
	save_results('last.results')	
	

#overplot overplots a previously-stored simulation on the current plot in a color other than red, but does not touch the simulation active in the main program. 'x' must be entered as a string with quotes.

def overplot(x):

	'''
	overplots a previously-stored simulation on the current plot in a color other than red, but does not touch the simulation active in the main program. 'x' must be entered as a string with quotes.
	'''

	global elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,vlsr,frequency,freq_sim,intensity,int_sim,current,fig,ax

	#store the currently-active simulation in a temporary holding cell

	freq_temp = freq_sim 
	int_temp = int_sim 
	
	#pull the old simulation out of storage
	
	freq_sim = sim[x].freq_sim
	int_sim = sim[x].int_sim

	if any(line for line in ax.lines if line.get_label()==sim[x].name):
		line_color = [line.get_color() for line in ax.lines if line.get_label()==sim[x].name][0]
		line_style = [line.get_linestyle() for line in ax.lines if line.get_label()==sim[x].name][0]
		clear_line(sim[x].name)
	else:
		line_color = next(colors)
		line_style = '-'		

#Not working at the moment
		
# 	if any(line for line in ax.lines if line.get_color()==line_color):
# 		match = False
# 		i = 1
# 		while match == False:
# 			#first we just try changing the colors.
# 			line_color = next(colors)
# 			if any(line for line in ax.lines if line.get_color()==line_color) == False:
# 				match = True
# 			line = [line for line in ax.lines if line.get_color()==line_color][0]
# 			line_style_tmp = line.get_linestyle()
# 			line_style = next(styles)
# 			#if the line style hasn't been used before for that color, we have a good match and we exit
# 			if line_style != line_style_tmp:
# 				match = True
# 			#otherwise we go to the next line style and try again.  We loop long enough to try all the combinations, and then if we can't find any unique combinations, we go back and just re-use the color w/ the basic line.  Too bad.  Add more colors.
# 			else:
# 				line_style = next(styles)
# 				if i > 4:
# 					line_style = '-'
# 					match = True			
	
	lines[sim[x].name] = ax.plot(freq_sim,int_sim,color = line_color, linestyle=line_style, label = sim[x].name)
	
	ax.legend()
	fig.canvas.draw()
	
	freq_sim = freq_temp
	int_sim = int_temp
	
	save_results('last.results')
		
#load_mol loads a new molecule into the system.  Make sure to store the old molecule simulation first, if you want to get it back.  The current graph will be updated with the new molecule.  Catalog file must be given as a string, without the *.cat as usual.  Simulation will begin with the same T, dV, S, vlsr as previous, so change those first if you want.

def load_mol(x):

	'''
	loads a new molecule into the system.  Make sure to store the old molecule simulation first, if you want to get it back.  The current graph will be updated with the new molecule.  Catalog file must be given as a string, without the *.cat as usual.  Simulation will begin with the same T, dV, S, vlsr as previous, so change those first if you want.
	'''

	global frequency,logint,qn7,qn8,qn9,qn10,qn11,qn12,elower,eupper,intensity,qns,catalog,catalog_file,fig,current,fig,ax,freq_sim,int_sim,first_run
	
	
	if first_run == False:
		try:
			clear_line('current')
		except:
			pass
	
	current = x
	
	catalog_file = x
	
	raw_array = read_cat(catalog_file)

	trimmed_array = trim_raw_array(raw_array)
	
	try:
		trimmed_array[0]
	except:
		print('There were no lines in the frequency range for {}.' .format(x))
		return False

	catalog = splice_array(trimmed_array)

	frequency = np.copy(catalog[0])
	logint = np.copy(catalog[2])
	qn7 = np.asarray(catalog[14])
	qn8 = np.asarray(catalog[15])
	qn9 = np.asarray(catalog[16])
	qn10 = np.asarray(catalog[17])
	qn11 = np.asarray(catalog[18])
	qn12 = np.asarray(catalog[19])
	elower = np.asarray(catalog[4])

	eupper = np.empty_like(elower)

	eupper = elower + frequency/29979.2458

	qns = det_qns(qn7,qn8,qn9,qn10,qn11,qn12) #figure out how many qns we have for the molecule

	intensity = convert_int(logint)
	
	tmp_freq = np.copy(frequency)
	
	tmp_freq += (-vlsr)*tmp_freq/ckm
	
	freq_sim,int_sim=run_sim(tmp_freq,intensity,T,dV,S)
	
	#check if a plot is already open.  If it is, nothing happens.  If there are no plots open, plt.get_fignums()[0] is empty and throws an IndexError, which the exception catches and just makes a fresh plot.  If it's the first time the program has been run since it was loaded, it ignores this check and just makes a new plot.  
	
	if first_run == True:
		make_plot()
		first_run = False
		return True
	else:
		try:
			plt.get_fignums()[0]
		except:	
			make_plot()
			return True
		
	#if there is a plot open, we just update the current simulation
	
	if gauss == False:

		lines['current'] = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red',label='current',zorder=500) #Plot sticks from TA down to 0 at each point in freq.

	else:

		lines['current'] = ax.plot(freq_sim,int_sim,color = 'red',label='current',zorder=500)	

	ax.legend()
	fig.canvas.draw()	
	
	save_results('last.results')
	
	return True
		
#clear_line removes a line labeled 'x' from the current plot window.  x must be a string

def clear_line(x):

	'''
	removes a line labeled 'x' from the current plot window.  x must be a string
	'''	
	
	try:
		line = lines.pop(x)
	except:
		return
	
	try:
		line.remove()
		ax.legend()
		fig.canvas.draw()
	except:
		line[0].remove()
		ax.legend()
		fig.canvas.draw()		


#save_results prints out a file with all of the parameters for each molecule *which has been stored.*  It will not print the active molecule unless it has been stored. 'x' must be a string, and the output will go to there.

def save_results(x):

	'''
	save_results prints out a file that contains at the top the output of status(), followed by all of the parameters for each molecule which has been stored.  This file can then be used with restore() to restore the parameter space and keep going.
	'''

	with open(x, 'w') as output:
	
		output.write('viewspectrum.py version {}\n' .format(version))		
		output.write('saved: {}\n\n' .format(datetime.now().strftime("%m-%d-%y %H:%M:%S")))
		
		output.write('#### Active Simulation ####\n\n')
	
		output.write('molecule:\t{}\n' .format(current))
		output.write('obs:\t{}\n' .format(spec))
		output.write('T:\t{} K\n' .format(T))
		output.write('S:\t{:.2f}\n' .format(S))
		output.write('dV:\t{:.2f} km/s\n' .format(dV))
		output.write('VLSR:\t{:.2f} km/s\n' .format(vlsr))
		output.write('ll:\t{} MHz\n' .format(ll))
		output.write('ul:\t{} MHz\n' .format(ul))
		output.write('CT:\t{} K\n' .format(CT))
		output.write('gauss:\t{}\n' .format(gauss))
		output.write('catalog_file:\t{}\n' .format(catalog_file))
		output.write('thermal:\t{} K\n\n' .format(thermal))
	
		output.write('#### Stored Simulations ####\n\n')
		
		output.write('Molecule\tT(K)\tS\tdV\tvlsr\tCT\tcatalog_file\n')
	
		for molecule in sim:
		
			output.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\n' .format(sim[molecule].name,sim[molecule].T,sim[molecule].S,sim[molecule].dV,sim[molecule].vlsr,sim[molecule].CT,sim[molecule].catalog_file))
			
		output.write('\n#### Active Graph Status ####\n\n')
		
		output.write('Label\tColor\tStyle\n')
		
		for line in ax.lines:
		
			output.write('{}\t{}\t{}\n' .format(line.get_label(),line.get_color(),line.get_linestyle()))	
	
#status prints the current status of the program and the various key variables.

def status():

	'''
	prints the current status of the program and the various key variables
	'''

	print('Current Molecule:\t {}' .format(current))
	print('Current Catalog:\t {}' .format(catalog_file))
	print('Current lab or observation: \t {}' .format(spec))
	print('T: \t {} K' .format(T))
	print('S: \t {}' .format(S))
	print('dV: \t {} km/s' .format(dV))
	print('VLSR: \t {} km/s' .format(vlsr))
	print('ll: \t {} MHz' .format(ll))
	print('ul: \t {} MHz' .format(ul))
	print('CT: \t {} K' .format(CT))
	print('gauss: \t {}' .format(gauss))
		
#sum_stored creates a combined spectrum of all stored molecule simulations and stores it in freq_sum and int_sum.  This might take a while.  It's done from scratch because the frequency axes in freq_sim stored for each molecule will not necessarily be the same, so co-adding is difficult without re-gridding, which well, maybe later.	

def sum_stored():

	'''
	Creates a combined spectrum of all stored molecule simulations and stores it in freq_sum and int_sum.  This might take a while.
	'''
	global freq_sum, int_sum
	
	freq_gauss = []

	for x in sim:
	
		tmp_freq = np.copy(sim[x].frequency)
		
		tmp_freq += (-sim[x].vlsr)*tmp_freq/ckm	

		for y in range(len(sim[x].intensity)):
	
			l_f = sim[x].dV*tmp_freq[y]/ckm #get the FWHM in MHz
	
			min_f = tmp_freq[y] - 10*l_f #get the frequency 10 FWHM lower
			max_f = tmp_freq[y] + 10*l_f #get the frequency 10 FWHM higher
	
			res_pnts = l_f / 15 #determine a resolution element (15 points across the line)
	
			freq_line = np.arange(min_f,max_f,res_pnts)
	
			freq_gauss.extend(freq_line)
	
	int_gauss = [0] * len(freq_gauss)
	
	freq_gauss.sort()
	
	del tmp_freq
		
	for x in sim:
	
		tmp_freq = np.copy(sim[x].frequency)
		
		tmp_freq += (-sim[x].vlsr)*tmp_freq/ckm	
		
		tmp_int = np.copy(sim[x].intensity)
		
		tmp_int = scale_temp(tmp_int,sim[x].qns,sim[x].elower,sim[x].qn7,sim[x].qn8,sim[x].qn9,sim[x].qn10,sim[x].qn11,sim[x].qn12,sim[x].T,sim[x].CT,sim[x].catalog_file)
		
		tmp_int *= sim[x].S
		
		tmp_int[tmp_int > thermal] = thermal
	
		for y in range(len(tmp_int)):
		
			l_f = sim[x].dV*tmp_freq[y]/ckm #get the FWHM in MHz
			
			c = l_f/2.35482

			int_gauss += tmp_int[y]*exp(-((freq_gauss - tmp_freq[y])**2/(2*c**2)))
			
	int_gauss[int_gauss > thermal] = thermal
	
	freq_sum = freq_gauss
	int_sum = int_gauss		

#overplot_sum overplots the summed spectrum of all stored molecules as created by sum_stored() on the current plot, in green.

def overplot_sum():

	'''
	Overplots the summed spectrum of all stored molecules as created by sum_stored() on the current plot, in green.
	'''	
	
	line_color = '#00FF00'
	
	if any(line for line in ax.lines if line.get_label()=='sum'):
		clear_line('sum')
	
	lines['sum'] = ax.plot(freq_sum,int_sum,color = line_color, label = 'sum', gid='sum', linestyle = '-',zorder=25)
	
	ax.legend()
	fig.canvas.draw()
	
#restore restores the state of the program from a save file, loading all stored spectra into memory, loads the previously active simulation into current, and restores the last active graph. x is a string with the filename of the restore file. The catalog files must be present, and named to match those in the save file.

def restore(x):

	'''
	restores the state of the program from a save file, loading all stored spectra into memory, loads the previously active simulation into current, and restores the last active graph. x is a string with the filename of the restore file. The catalog files must be present, and named to match those in the save file.
	'''

	global frequency,logint,qn7,qn8,qn9,qn10,qn11,qn12,elower,eupper,intensity,qns,catalog,catalog_file,fig,current,fig,ax,freq_sim,int_sim,T,dV,S,vlsr,ll,ul,CT,gauss,first_run,thermal,sim
	
	#empty out the previously-stored simulations
	
	sim = {}

	#read in the save file as an array line by line
	
	restore_array = []
	
	try:
		with open(x) as input:
			for line in input:
				restore_array.append(line)
	except TypeError:
		print('Input must be a string in \'\'.')
		return
	except FileNotFoundError:
		print('This is not the file you are looking for.')
		return
	
	#check if the file that was read it is actually a savefile for viewspectrum.py
	
	if restore_array[0].split()[0] != 'viewspectrum.py':
		print('The file is not a viewspectrum.py save file, has been altered, or was created with an older version of the program.  In any case, I can\'t read it, sorry.')
		return	
		
	#let's grab the date and time, so we can print some nice information to the terminal
	
	restore_date = restore_array[1].split()[1]
	restore_time = restore_array[1].split()[2]	
		
	#separate out the sections into their own arrays

	active_array = []	
	stored_array = []
	graph_array = []
	
	#figure out where each section starts
	
	active_index = restore_array.index('#### Active Simulation ####\n')
	stored_index = restore_array.index('#### Stored Simulations ####\n')
	graph_index = restore_array.index('#### Active Graph Status ####\n')
	
	for i in range(active_index+2,stored_index-1):
	
		active_array.append(restore_array[i])
		
	for i in range(stored_index+3,graph_index-1):
	
		stored_array.append(restore_array[i])
		
	for i in range(graph_index+3,len(restore_array)):
	
		graph_array.append(restore_array[i])
		
	#just to be safe, let's set the upper limits, lower limits, gaussian toggles, and thermal values now.
	
	
	if active_array[9].split('\t')[1].strip('\n') == 'True':
		gauss = True
	else:
		gauss = False
	ll = float(active_array[6].split('\t')[1].strip(' MHz\n'))
	ul = float(active_array[7].split('\t')[1].strip(' MHz\n'))
	thermal = float(active_array[11].split('\t')[1].strip(' K\n'))
	
	#OK, now time to do the hard part.  As one always should, let's start with the middle part of the whole file, and load and then store all of the simulations.
	
	for i in range(len(stored_array)):
	
		name = stored_array[i].split('\t')[0]
		T = float(stored_array[i].split('\t')[1])
		S = float(stored_array[i].split('\t')[2])
		dV = float(stored_array[i].split('\t')[3])
		vlsr = float(stored_array[i].split('\t')[4])
		CT = float(stored_array[i].split('\t')[5])
		catalog_file = str(stored_array[i].split('\t')[6]).strip('\n').strip()
		
# 		try:
		first_run = True
		try:
			foo = load_mol(catalog_file)
			if foo == False:
				raise ValueError()
		except ValueError:
			continue
# 		except FileNotFoundError:
# 			print('I was unable to locate the catalog file {} for molecule entry {}.  Sorry.' .format(catalog_file,name))
# 			return
			
		store(name)
		
		close()
	
	#Now we move on to loading in the currently-active molecule
	
	catalog_file = active_array[10].split('\t')[1].strip('\n')
	try:
		obs = active_array[1].split('\t')[1].strip('\n')
		read_obs(obs)
	except:
		pass
	T = float(active_array[2].split('\t')[1].strip(' K\n'))
	S = float(active_array[3].split('\t')[1].strip('\n'))
	dV = float(active_array[4].split('\t')[1].strip(' km/s\n'))
	vlsr = float(active_array[5].split('\t')[1].strip(' km/s\n'))
	CT = float(active_array[8].split('\t')[1].strip(' K\n'))
	current = active_array[0].split('\t')[1].strip('\n')
	
	try:
		first_run = True
		load_mol(catalog_file)
	except FileNotFoundError:
		print('I was unable to locate the catalog file {} for molecule entry {}.  Sorry.' .format(catalog_file,name))
		return
		
	#And finally, overplot anything that was on the plot previously
	
	for i in range(len(graph_array)):
	
		name = graph_array[i].split('\t')[0]
		line_color = graph_array[i].split('\t')[1]
		line_style = graph_array[i].split('\t')[2].strip('\n')
		
		if name == 'obs':
		
			continue
			
		elif name == 'current':
		
			continue
			
		elif name == 'sum':
		
			sum_stored()
			overplot_sum()
			continue			
			
		else:
			
			lines[name] = ax.plot(sim[name].freq_sim,sim[name].int_sim,color = line_color, linestyle=line_style, label = name)	
			
		ax.legend()
		fig.canvas.draw()
		
	#If we made it here, we were successful, so let's print out what we did
	
	print('Successfully restored from file {} which was saved on {} at {}.' .format(x,restore_date,restore_time))	
	
#fix_legend allows you to change the legend to meet its size needs.

def fix_legend(x,lsize):

	'''
	Modifies the legend to have x columns.  lsize must be a string. 'xx-small' 'x-small' and 'small' all work.
	'''

	plt.legend(ncol=x,prop={'size':lsize})
	
	fig.canvas.draw()	
	
#############################################################
#							Classes for Storing Results		#
#############################################################	

class Molecule(object):

	def __init__(self,name,catalog_file,elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,CT,vlsr,frequency,freq_sim,intensity,int_sim):
	
		self.name = name
		self.catalog_file = catalog_file
		self.elower = elower
		self.eupper = eupper
		self.qns = qns
		self.logint = logint
		self.qn7 = qn7
		self.qn8 = qn8
		self.qn9 = qn9
		self.qn10 = qn10
		self.qn11 = qn11
		self.qn12 = qn12
		self.S = S
		self.dV = dV
		self.T = T
		self.CT = CT
		self.vlsr = vlsr
		self.frequency = frequency
		self.freq_sim = freq_sim
		self.intensity = intensity
		self.int_sim = int_sim

	
#############################################################
#							Run Program						#
#############################################################

sim = {} #dictionary to hold stored simulations

lines = {} #dictionary to hold matplotlib lines

freq_obs = [] #to hold laboratory or observational spectra
int_obs = []

freq_sim = [] #to hold simulated spectra
int_sim = [] 

freq_sum = [] #to hold combined spectra
int_sum = []

current = catalog_file

colors = itertools.cycle(['#ff8172','#514829','#a73824','#7b3626','#a8ac87','#8c8e64','#974710','#d38e20','#ce9a3a','#ae7018','#ac5b14','#64350f','#b18f59','#404040','#791304','#1f2161','#171848','#3082fe','#2c5b5e','#390083','#5c65f8','#6346fa','#3c3176','#1cf6ba','#c9bcf0','#90edfc','#3fb8ee','#b89b33','#e7d17b'])
styles = itertools.cycle(['-','--','-.',':'])
