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
import time
import random
import warnings
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import itertools
#warnings.filterwarnings('error')

version = 3.0

h = 6.626*10**(-34) #Planck's constant in J*s
k = 1.381*10**(-23) #Boltzmann's constant in J/K
kcm = 0.69503476 #Boltzmann's constant in cm-1/K
ckm = 2.998*10**5 #speed of light in km/s

#############################################################
#							Defaults 						#
#############################################################

T = 300 #temperature for simulations.  Default is 300 K.

catalog_file = None #catalog file to load in.  Needs to be a string.

vlsr = 0.0 #vlsr offset applied to simulation.  Default is 0 km/s.

S = 1 #linear scaling factor applied to simulation.  Default is 1.

ll = float('-inf') #lower limit for the simulation range.  Default is none.

ul = float('inf') #upper limit for the simulation range.  Default is none.

spec = None	#file of a laboratory or observational spectrum to load in for comparison.  Default is none.

dV = 5.0 #linewidth of the simulation.  Default is 5.0 km/s.

R = 0.5 #resolution of the simulation if using Gaussians.  Default is 0.5 km/s.
	
CT = 300.0 #temperature the catalog is simulated at.  Default is 300 K.

gauss = True #toggle for simulating Gaussians or a stick spectrum.  Default is True.

#############################################################
#							Functions						#
#############################################################
	
	
#read_cat reads the catalog file in

def read_cat(catalog_file):

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

def trim_raw_array(raw_array):

	tmp_array = np.arange(len(raw_array),dtype=np.float)	#create a temporary array to splice out the frequency
	
	for line in range(len(raw_array)):
	
		tmp_array[line] = float(str(raw_array[line][:13]).strip())
		
	i = np.where(tmp_array > ll)[0][0]		#get the index of the first value above the lower limit
	try:
		i2 = np.where(tmp_array > ul)[0][0]	#get the index of the first value above the upper limit
	except IndexError:
		i2 = len(tmp_array)					#if the catalog ends before the upper limit is reached
		
	else:
	
		print('You somehow invoked trim_raw_array, but ll and ul are gone.  Stop that.')
		return
		
	trimmed_array = []

	for x in range(i,i2):
	
		trimmed_array.append(raw_array[x])
		
	return trimmed_array				
	
#fix_qn fixes quantum number issues

def fix_qn(qnarray,line,old_qn):

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

def splice_array(raw_array):

	frequency = np.arange(len(raw_array),dtype=np.float)
	error = np.arange(len(raw_array),dtype=np.float)
	logint = np.arange(len(raw_array),dtype=np.float)
	dof = np.arange(len(raw_array),dtype=np.int)
	elower = np.arange(len(raw_array),dtype=np.float)
	gup = np.arange(len(raw_array),dtype=np.int)
	tag = np.arange(len(raw_array),dtype=np.int)
	qnformat = np.arange(len(raw_array),dtype=np.int)
	qn1 = np.arange(len(raw_array),dtype=object)
	qn2 = np.empty(len(raw_array),dtype=object)
	qn3 = np.empty(len(raw_array),dtype=object)
	qn4 = np.empty(len(raw_array),dtype=object)
	qn5 = np.empty(len(raw_array),dtype=object)
	qn6 = np.empty(len(raw_array),dtype=object)
	qn7 = np.empty(len(raw_array),dtype=object)
	qn8 = np.empty(len(raw_array),dtype=object)
	qn9 = np.empty(len(raw_array),dtype=object)
	qn10 = np.empty(len(raw_array),dtype=object)
	qn11 = np.empty(len(raw_array),dtype=object)
	qn12 = np.empty(len(raw_array),dtype=object)

	for line in range(len(raw_array)):
	
		frequency[line] = float(str(raw_array[line][:13]).strip())
		error[line] = float(str(raw_array[line][13:21]).strip())
		logint[line] = float(str(raw_array[line][21:29]).strip())
		dof[line] = int(str(raw_array[line][29:31]).strip())
		elower[line] = float(str(raw_array[line][31:41]).strip())
		gup[line] = int(str(raw_array[line][41:44]).strip())
		tag[line] = int(str(raw_array[line][44:51]).strip())
		qnformat[line] = int(str(raw_array[line][51:55]).strip())
		try:
			qn1[line] = int(raw_array[line][55:57]) if str(raw_array[line][54:57]).strip() else ''
		except ValueError:
			fix_qn(qn1,line,str(raw_array[line][55:57]))	
		try:		
			qn2[line] = int(raw_array[line][57:59]) if str(raw_array[line][57:59]).strip() else ''
		except ValueError:
			fix_qn(qn2,line,str(raw_array[line][57:59]))
		try:				
			qn3[line] = int(raw_array[line][59:61]) if str(raw_array[line][59:61]).strip() else ''
		except ValueError:
			fix_qn(qn3,line,str(raw_array[line][59:61]))
		try:
			qn4[line] = int(raw_array[line][61:63]) if str(raw_array[line][61:63]).strip() else ''
		except ValueError:
			fix_qn(qn4,line,str(raw_array[line][61:63]))
		try:
			qn5[line] = int(raw_array[line][63:65]) if str(raw_array[line][63:65]).strip() else ''
		except ValueError:
			fix_qn(qn5,line,str(raw_array[line][63:65]))
		try:
			qn6[line] = int(raw_array[line][65:67]) if str(raw_array[line][65:67]).strip() else ''
		except ValueError:
			fix_qn(qn6,line,str(raw_array[line][65:67]))
		try:
			qn7[line] = int(raw_array[line][67:69]) if str(raw_array[line][67:69]).strip() else ''
		except ValueError:
			fix_qn(qn7,line,str(raw_array[line][67:69]))
		try:
			qn8[line] = int(raw_array[line][69:71]) if str(raw_array[line][69:71]).strip() else ''
		except ValueError:
			fix_qn(qn8,line,str(raw_array[line][69:71]))
		try:
			qn9[line] = int(raw_array[line][71:73]) if str(raw_array[line][71:73]).strip() else ''
		except ValueError:
			fix_qn(qn9,line,str(raw_array[line][71:73]))
		try:
			qn10[line] = int(raw_array[line][73:75]) if str(raw_array[line][73:75]).strip() else ''
		except ValueError:
			fix_qn(qn10,line,str(raw_array[line][73:75]))
		try:
			qn11[line] = int(raw_array[line][75:77]) if str(raw_array[line][75:77]).strip() else ''
		except ValueError:
			fix_qn(qn11,line,str(raw_array[line][75:77]))
		try:
			qn12[line] = int(raw_array[line][77:]) if str(raw_array[line][77:]).strip() else ''
		except ValueError:
			fix_qn(qn12,line,str(raw_array[line][77:]))		
			

	return frequency,error,logint,dof,elower,gup,tag,qnformat,qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8,qn9,qn10,qn11,qn12
	
#det_qns determines how many qns represent each state

def det_qns(qn7,qn8,qn9,qn10,qn11,qn12):

	qns = 1
	
	try:
		int(qn8[0])
		qns = 2
	except ValueError:
		pass
		
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
	
def calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T):

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

def scale_temp(int_sim,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT):

	scaled_int = np.empty_like(int_sim)
	
	Q_T = calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T)
	Q_CT = calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,CT)
	
	scaled_int = int_sim * (Q_CT/Q_T) * exp(-(((1/T)-(1/CT))*elower)/0.695)

	return scaled_int

#convert_int converts the intensity to not log units

def convert_int(logint):

	intensity = np.empty_like(logint)	
	
	intensity = 10**(logint)
	
	return intensity
	
#simulates Gaussian profiles after intensities are simulated.

def sim_gaussian(int_sim,frequency,linewidth):
	
	freq_gauss = []

	for x in range(int_sim.shape[0]):
	
		l_f = linewidth*frequency[x]/3E5 #get the FWHM in MHz
	
		min_f = frequency[x] - 10*l_f #get the frequency 10 FWHM lower
		max_f = frequency[x] + 10*l_f #get the frequency 10 FWHM higher
	
		res_pnts = l_f / 15 #determine a resolution element (15 points across the line)
	
		freq_line = np.arange(min_f,max_f,res_pnts)
		int_line = np.arange(min_f,max_f,res_pnts)
	
		freq_gauss.extend(freq_line)
	
	int_gauss = [0] * len(freq_gauss)
	
	freq_gauss.sort()

	for x in range(int_sim.shape[0]):
	
		l_f = linewidth*frequency[x]/3E5 #get the FWHM in MHz
	
		c = l_f/2.35482

		int_gauss += int_sim[x]*exp(-((freq_gauss - frequency[x])**2/(2*c**2)))
		
	int_gauss[int_gauss > thermal] = thermal	
	
	return(freq_gauss,int_gauss)

#write_spectrum writes out the current freq and intensity to output_file

def write_spectrum(output_file):

	if gauss == True:
	
		for h in range(len(freq_sim)):

			if (h == 0): #write out the results to a out_file
			
				with open(output_file, 'w') as output: 
					
					output.write('{} {}\n' .format(freq_sim[0],int_sim[0]))
						
			else:
				
				with open(output_file, 'a') as output:
					
					output.write('{} {}\n' .format(freq_sim[h],int_sim[h]))
						
			h += 1		
	
	else:
	
		for h in range(freq_sim.shape[0]):

			if (h == 0): #write out the results to a out_file
				
				with open(output_file, 'w') as output: 
					
					output.write('{} {}\n' .format(freq_sim[0],int_sim[0]))
						
			else:
				
				with open(output_file, 'a') as output:
					
					output.write('{} {}\n' .format(freq_sim[h],int_sim[h]))
						
			h += 1				

#run_sim runs the simulation.  It's a meta routine, so that we can update later

def run_sim(frequency,intensity,T,dV,S):

	int_temp = scale_temp(intensity,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT)

	int_temp *= S
	
	if gauss == True:

		freq_sim,int_sim = sim_gaussian(int_temp,frequency,dV)
		
	else:
	
		int_temp[int_temp > thermal] = thermal
		freq_sim = frequency
		int_sim = int_temp
		
	return freq_sim,int_sim
	
#mod_T changes the temperature, re-simulates, and re-plots	
	
def modT(x):

	try:
		float(x)
	except:
		print('x needs to be a number.')
		return
		
	global T,freq_sim,int_sim
		
	T = float(x)
		
	freq_sim,int_sim = run_sim(frequency,intensity,T,dV,S)
	
	line1.set_ydata(int_sim)
	line1.set_xdata(freq_sim)
		
	fig.canvas.draw()
	
#modS changes the scaling, re-simulates, and re-plots

def modS(x):

	try:
		float(x)
	except:
		print('x needs to be a number.')
		return
		
	global S,freq_sim,int_sim
		
	S = x
		
	freq_sim,int_sim = run_sim(frequency,intensity,T,dV,S)
	
	line1.set_ydata(int_sim)
	line1.set_xdata(freq_sim)
		
	fig.canvas.draw()
	
#moddV changes the velocity width, re-simulates, and re-plots

def moddV(x):

	try:
		float(x)
	except:
		print('dV needs to be a number.')
		return
		
	global dV,freq_sim,int_sim
	
	dV = x
		
	freq_sim,int_sim = run_sim(frequency,intensity,T,dV,S)
	
	line1.set_ydata(int_sim)
	line1.set_xdata(freq_sim)
		
	fig.canvas.draw()
	
#modVLSR changes the LSR velocity, re-simulates, and re-plots

def modVLSR(x):

	try:
		float(x)
	except:
		print('vlsr needs to be a number.')
		return	
		
	global vlsr,freq_sim,int_sim,frequency
	
	vlsr = x
		
	frequency = np.copy(catalog[0])
	
	frequency += (-vlsr)*frequency/ckm		
	
	freq_sim,int_sim = run_sim(frequency,intensity,T,dV,S)
	
	line1.set_ydata(int_sim)
	line1.set_xdata(freq_sim)
	
	fig.canvas.draw()		

#modV is an alias for modVLSR

def modV(vlsr):

	modVLSR(vlsr)	
	
#make_plot makes the plot!

def make_plot():

	global fig,ax,line1,line2,freq_sim,intensity,int_sim
	
	freq_sim,int_sim=run_sim(frequency,intensity,T,dV,S)

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
	
	if spec:

		line2, = ax.plot(freq_obs,int_obs,color = 'black')

	if gauss == False:

		line1, = ax.vlines(freq_sim,0,int_sim,linestyle = '-',color = 'red') #Plot sticks from TA down to 0 at each point in freq.

	else:

		line1, = ax.plot(freq_sim,int_sim,color = 'red')	

	fig.canvas.draw()
	
#obs_off turns off the observations
	
def obs_off():

	if spec == None:
	
		print('There are no observations to turn off!')
		return
		
	line2.set_xdata([])
	line2.set_ydata([])
			
	fig.canvas.draw()	
	
#obs_on turns on the observations
	
def obs_on():

	if spec == None:
	
		print('There are no observations to turn on!')
		return	
		
	line2.set_xdata(freq_obs)
	line2.set_ydata(int_obs)
			
	fig.canvas.draw()	
	
# setup runs the initial setup for a molecule
# 
# def setup():
# 
# 	raw_array = read_cat(catalog_file)
# 
# 	raw_array = trim_raw_array(raw_array)
# 
# 	catalog = splice_array(raw_array)
# 
# 	frequency = np.copy(catalog[0])
# 	logint = np.copy(catalog[2])
# 	qn7 = np.asarray(catalog[14])
# 	qn8 = np.asarray(catalog[15])
# 	qn9 = np.asarray(catalog[16])
# 	qn10 = np.asarray(catalog[17])
# 	qn11 = np.asarray(catalog[18])
# 	qn12 = np.asarray(catalog[19])
# 	elower = np.asarray(catalog[4])
# 
# 	eupper = np.empty_like(elower)
# 
# 	eupper = elower + frequency/29979.2458
# 
# 	qns = det_qns(qn7,qn8,qn9,qn10,qn11,qn12) #figure out how many qns we have for the molecule
# 
# 	intensity = convert_int(logint)
# 	
# 	if spec:
# 	
# 		freq_obs,int_obs = read_obs()
# 
# 	return frequency,logint,qn7,qn8,qn9,qn10,qn11,qn12,elower,eupper,intensity,qns,catalog
	
#read_obs reads in observations or laboratory spectra and populates freq_obs and int_obs

def read_obs():

	obs = read_cat(spec)
	
	global freq_obs,int_obs

	freq_obs = []
	int_obs = []

	for x in range(len(obs)):

		freq_obs.append(float(obs[x].split()[0]))
		int_obs.append(float(obs[x].split()[1].strip('\n')))	
	
#store saves the current simulation parameters for recall later.  *Not* saved as a Gaussian. 'x' must be entered as a string with quotes.

def store(x):

	setup()
	
	intensity = convert_int(logint)
	
	freq_sim,int_sim = run_sim(frequency,intensity,T,dV,S)
	
	sim[x] = Molecule(x,catalog_file,elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,vlsr,frequency,freq_sim,intensity,int_sim) 
	
#recall wipes the current simulation and re-loads a previous simulation that was stored with store(). 'x' must be entered as a string with quotes. This will close the currently-open plot.

def recall(x):

	global elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,vlsr,frequency,freq_sim,intensity,int_sim,current
	
	plt.close()

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
		
	make_plot()	

#overplot overplots a previously-stored simulation on the current plot in a color other than red, but does not touch the simulation active in the main program. 'x' must be entered as a string with quotes.

def overplot(x):

	global elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,vlsr,frequency,freq_sim,intensity,int_sim,current

	#store the currently-active simulation in a temporary holding cell

	freq_temp = freq_sim 
	int_temp = int_sim 
	
	#pull the old simulation out of storage
	
	freq_sim = sim[x].freq_sim
	int_sim = sim[x].int_sim
	
	line_color = next(colors)
	
	ax.plot(freq_sim,int_sim,color = line_color)
	
	fig.canvas.draw()
	
	freq_sim = freq_temp
	int_sim = int_temp
		
#load_mol loads a new molecule into the system.  Make sure to store the old molecule simulation first, if you want to get it back.  The current graph will be updated with the new molecule.  Catalog file must be given as a string, without the *.cat as usual.  Simulation will begin with the same T, dV, S, vlsr as previous, so change those first if you want.

def load_mol(x):

	global frequency,logint,qn7,qn8,qn9,qn10,qn11,qn12,elower,eupper,intensity,qns,catalog,catalog_file
	
	catalog_file = x
	
	frequency,logint,qn7,qn8,qn9,qn10,qn11,qn12,elower,eupper,intensity,qns,catalog = setup()
	
	frequency += (-vlsr)*frequency/ckm
	
	freq_sim,int_sim=run_sim(frequency,intensity,T,dV,S)
	
	line1.set_ydata(int_sim)
	line1.set_xdata(freq_sim)
	
	fig.canvas.draw()

#save_results prints out a file with all of the parameters for each molecule *which has been stored.*  It will not print the active molecule unless it has been stored. 'x' must be a string, and the output will go to there.

def save_results(x):

	'''
	'''

	with open(x, 'w') as output:
	
		output.write('ll {}\n' .format(ll))
		output.write('ul {}\n' .format(ul))
		output.write('obs {}\n' .format(spec))
	
		output.write('Molecule \t T(K) \t S \t dV \t vlsr \t catalog_file\n')
	
		for molecule in sim:
		
			output.write('{} \t {} \t {} \t {} \t {} \t {}\n' .format(sim[molecule].name,sim[molecule].T,sim[molecule].S,sim[molecule].dV,sim[molecule].vlsr,sim[molecule].catalog_file))
	
#status prints the current status of the program and the various key variables.

def status():

	'''
	prints the current status of the program and the various key variables
	'''

	print('Current Molecule or Catalog:\t {}' .format(current))
	print('T: \t {} K' .format(T))
	print('S: \t {}' .format(S))
	print('dV: \t {} km/s' .format(dV))
	print('VLSR: \t {} km/s' .format(vlsr))
	print('R: \t {} km/s' .format(R))
	print('ll: \t {} MHz' .format(ll))
	print('ul: \t {} MHz' .format(ul))
	print('CT: \t {} K' .format(CT))
	print('gauss: \t {}' .format(gauss))
		
	
#############################################################
#							Classes for Storing Results		#
#############################################################	

class Molecule(object):

	def __init__(self,name,catalog_file,elower,eupper,qns,logint,qn7,qn8,qn9,qn10,qn11,qn12,S,dV,T,vlsr,frequency,freq_sim,intensity,int_sim):
	
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

thermal = float('inf') #initial default cutoff for optically-thick lines (i.e. don't touch them unless thermal is modified.)

current = catalog_file

colors = itertools.cycle(['#9400D3','#4B0082','#0000FF','#00FF00','#FFFF00','#FF7F00'])

# frequency,logint,qn7,qn8,qn9,qn10,qn11,qn12,elower,eupper,intensity,qns,catalog = setup()
# 
# make_plot()
