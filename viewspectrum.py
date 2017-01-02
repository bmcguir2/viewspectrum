#!/usr/bin/python

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

#############################################################
#							Preamble						#
#############################################################

import os, sys, argparse, math
import numpy as np
from numpy import exp as exp
import time
import random
import warnings
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
#warnings.filterwarnings('error')

version = 1.8

h = 6.626*10**(-34) #Planck's constant in J*s
k = 1.381*10**(-23) #Boltzmann's constant in J/K
kcm = 0.69503476 #Boltzmann's constant in cm-1/K
ckm = 2.998*10**5 #speed of light in km/s

#############################################################
#							Read Input						#
#############################################################

parser = argparse.ArgumentParser(description='Displays a spectrum from a catalog file at the temperature at which the catalog file was simulated.  Can optionally scale with temperature.  Relative intensities will be accurate - absolute intensities will depend on the completeness of the catalog, which is used to calculate Q.')
parser.add_argument('basename', help='Catalog file basename of the form: basename.cat  (i.e., do not include .cat)')
parser.add_argument('-ll', type=float, help='lower frequency limit (MHz)')
parser.add_argument('-ul', type=float, help='upper frequency limit (MHz)')
parser.add_argument('-T', type=float, help='Rotational Temperature (K)')
parser.add_argument('-CT', type=float, help='Temperature catalog file was simulated at. Default is 300 K.')
parser.add_argument('-W', '--write', action='store_true', help='write the resulting spectrum to an ascii file (basename.spec)')
parser.add_argument('-NP', '--noplot', action='store_true', help='skip plotting')
parser.add_argument('-G', '--gauss', action='store_true', help='Simulate Gaussians.  Default linewidth is 5 km/s, with a resolution of 0.5 km/s.  Modify this with the -dV and -R options.  Warning: this may take a very, very long time.')
parser.add_argument('-dV', type=float, help='velocity width for simulated Gaussians in km/s.  Default is 5.')
parser.add_argument('-R', type=float, help='Resolution for simulated Gaussians in km/s.  Default is 0.5.  Currently does not work.')
parser.add_argument('-lsr', type=float, help='Apply a VLSR correction in km/s to the simulated spectrum.')
parser.add_argument('-spec', type=str, help='Underplot a lab or astronomical spectrum in a given filename.')
parser.add_argument('-specS', type=float, help='An optional intensity scaling factor to apply to the simulation.')
args = parser.parse_args()

catalog_file = '{}.cat' .format(args.basename.strip('.')) #Catalog file
output_file = '{}.spec' .format(args.basename.strip('.')) #Output file
	
T = args.T #Rotational temperature
basename = args.basename.strip('.')

vlsr = 0.0

if args.lsr:

	vlsr = args.lsr

if args.ll:

	ll = args.ll #lower limit to plot

if args.ul:

	ul = args.ul #upper limit to plot

if args.spec:

	spec = args.spec
	
if args.specS:

	specS = args.specS

dV = 5.0

if args.dV:

	dV = args.dV
	
R = 0.5

if args.R:

	R = args.R

if args.CT:

	CT = args.CT #Temperature catalog file was simulated at, if not 300 K.
	
else:

	CT = 300.0

#############################################################
#							Functions						#
#############################################################

#read_cat reads the catalog file in

def read_cat(catalog_file):

	my_array = []

	with open(catalog_file) as input:
	
		for line in input:
		
			my_array.append(line)
			
	return my_array	
	
#trim_raw_array takes the raw_array and trims it at the specified ll and ul

def trim_raw_array(raw_array):

	tmp_array = np.arange(len(raw_array),dtype=np.float)	#create a temporary array to splice out the frequency
	
	for line in range(len(raw_array)):
	
		tmp_array[line] = float(str(raw_array[line][:13]).strip())
		
	if args.ll and args.ul:
	
		i = np.where(tmp_array > ll)[0][0]		#get the index of the first value above the lower limit
		i2 = np.where(tmp_array > ul)[0][0]		#get the index of the first value above the upper limit
			
	elif args.ll:
	
		i = np.where(tmp_array > ll)[0][0]		#get the index of the first value above the lower limit
		i2 = len(raw_array)						#go to the upper limit of the catalog
		
	elif args.ul:
	
		i = 0									#start at the beginning of the catalog
		i2 = np.where(tmp_array > ul)[0][0]		#get the index of the first value above the upper limit
		
	else:
	
		print('You somehow invoked trim_raw_array without providing any limits.  Stop that.')
		
		quit()
		
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
	
	if basename=='acetone':
	
		Q = 2.91296*10**(-7)*T**6 - 0.00021050085*T**5 + 0.05471337*T**4 - 5.5477*T**3 + 245.28*T**2 - 2728.3*T + 16431 #Hard code for Acetone
		
	if basename=='sh':
	
		Q = 0.000000012549467*T**4 - 0.000008528126823*T**3 + 0.002288160909445*T**2 + 0.069272946237033*T + 15.357239728157400
 #Hard code for SH.  Completely unreliable below 2.735 K or above 300 K.

	if basename=='h2s':
	
		Q = -0.000004859941547*T**3 + 0.005498622332982*T**2 + 0.507648423477309*T - 1.764494755639740 #Hard code for H2S.  Completely unreliable below 2.735 K or above 300 K.
	if basename=='methanol' or basename=='CH3OH' or basename=='ch3oh' or basename=='Methanol':
	
		Q = 4.83410*10**-11*T**6 - 4.04024*10**-8*T**5 + 1.27624*10**-5*T**4 - 1.83807*10**-3*T**3 + 2.05911*10**-1*T**2 + 4.39632*10**-1*T -1.25670
	#Q = (-2*10**(-11)*T**6 + 1.4816*10**(-8)*T**5 - 3.649*10**(-6)*T**4 + 0.0002657*T**3 + 0.0518856*T**2 + 2.6652559*T + 3.0420215) #hard code for unresolved HNCS
	
	#Q = 8.098178*T**1.486328 #hard code for ethylene oxide
	
	return Q

#scale_temp scales the intensities to the new temperature

def scale_temp(int,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT):

	scaled_int = np.empty_like(int)
	
	Q_T = calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T)
	Q_CT = calc_q(qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,CT)
	
	scaled_int = int * (Q_CT/Q_T) * exp(-(((1/T)-(1/CT))*elower)/0.695)

	return scaled_int

#convert_int converts the intensity to not log units

def convert_int(logint):

	int = np.empty_like(logint)	
	
	int = 10**(logint)
	
	return int
	
#simulates Gaussian profiles after intensities are simulated.

def sim_gaussian(TA,frequency,linewidth):
	
	freq_gauss = []

	for x in range(TA.shape[0]):
	
		l_f = linewidth*frequency[x]/3E5 #get the FWHM in MHz
	
		min_f = frequency[x] - 10*l_f #get the frequency 10 FWHM lower
		max_f = frequency[x] + 10*l_f #get the frequency 10 FWHM higher
	
		res_pnts = l_f / 15 #determine a resolution element (15 points across the line)
	
		freq_line = np.arange(min_f,max_f,res_pnts)
		int_line = np.arange(min_f,max_f,res_pnts)
	
		freq_gauss.extend(freq_line)
	
	int_gauss = [0] * len(freq_gauss)
	
	freq_gauss.sort()

	for x in range(TA.shape[0]):
	
		l_f = linewidth*frequency[x]/3E5 #get the FWHM in MHz
	
		c = l_f/2.35482

		int_gauss += TA[x]*exp(-((freq_gauss - frequency[x])**2/(2*c**2)))
	
	return(freq_gauss,int_gauss)

#plot_spectrum will show a stick spectrum of whatever frequency and intensity values given to it.
	
def plot_spectrum(freq,TA):

	minorLocator = AutoMinorLocator(5)
	fig, ax = plt.subplots()
	plt.xlabel('Frequency (MHz)')
		
	plt.ylabel('Intensity (Probably Arbitrary)')
		
	plt.locator_params(nbins=4) #Use only 4 actual numbers on the x-axis
	ax.xaxis.set_minor_locator(minorLocator) #Let the program calculate some minor ticks from that
		
	ax.get_xaxis().get_major_formatter().set_scientific(False) #Don't let the x-axis go into scientific notation
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	
	fig.patch.set_facecolor('white')
	
	if args.gauss == False:
			
		plt.vlines(freq,0,TA,linestyle = '-',color = 'red') #Plot sticks from TA down to 0 at each point in freq.
		
	else:
	
		plt.plot(freq,TA,color = 'red')

	plt.show(block=False) #Plot it!
	
	
#plot_spectrum_obs will show a spectrum of whatever frequency and intensity values given to it overplotted on a lab or observational spectrum
	
def plot_spectrum_obs(freq,TA,freq_obs,int_obs):

	minorLocator = AutoMinorLocator(5)
	fig, ax = plt.subplots()
	plt.xlabel('Frequency (MHz)')
		
	plt.ylabel('Intensity (Probably Arbitrary)')
		
	plt.locator_params(nbins=4) #Use only 4 actual numbers on the x-axis
	ax.xaxis.set_minor_locator(minorLocator) #Let the program calculate some minor ticks from that
		
	ax.get_xaxis().get_major_formatter().set_scientific(False) #Don't let the x-axis go into scientific notation
	ax.get_xaxis().get_major_formatter().set_useOffset(False)
	
	fig.patch.set_facecolor('white')
	
	plt.plot(freq_obs,int_obs,color = 'black')
	
	if args.gauss == False:
			
		plt.vlines(freq,0,TA,linestyle = '-',color = 'red') #Plot sticks from TA down to 0 at each point in freq.
		
	else:
	
		plt.plot(freq,TA,color = 'red',gid = id)
	
	plt.show(block=False) #Plot it!
		

#write_spectrum writes the spectrum to the output file.
	
def write_spectrum(frequency,intensity,output_file):

	if args.gauss:
	
		for h in range(len(frequency)):

			if (h == 0): #write out the results to a out_file
			
				with open(output_file, 'w') as output: 
					
					output.write('{} {}\n' .format(frequency[0],intensity[0]))
						
			else:
				
				with open(output_file, 'a') as output:
					
					output.write('{} {}\n' .format(frequency[h],intensity[h]))
						
			h += 1		
	
	else:
	
		for h in range(frequency.shape[0]):

			if (h == 0): #write out the results to a out_file
				
				with open(output_file, 'w') as output: 
					
					output.write('{} {}\n' .format(frequency[0],intensity[0]))
						
			else:
				
				with open(output_file, 'a') as output:
					
					output.write('{} {}\n' .format(frequency[h],intensity[h]))
						
			h += 1				

#run_sim runs the simulation.  It's a meta routine, so that we can update later

def run_sim(frequency,intensity,T,dV,specS):

	freq = frequency

	if args.T:

		int = scale_temp(intensity,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT)
	
	if args.specS:

		int *= specS
	
	if args.gauss:

		freq,int = sim_gaussian(int,frequency,dV)
		
	return freq,int

#make_plots is a meta function for making the plots appear

def make_plots(freq,int):

	if args.noplot:

		pass
	
	elif args.spec:

		plot_spectrum_obs(freq,int,freq_obs,int_obs)	

	else:
	
		plot_spectrum(freq,int)

#update_T updates the temperature, redoes the simulation, and spits everything back out again.

def update_T(new_temp):
	
	global T,freq,int
	
	T = new_temp
	
	freq,int = run_sim(frequency,intensity,T,dV,specS)
	
	make_plots(freq,int)
	
#update_dV updates the linewidth, redoes the simulation, and spits everything back out again.

def update_dV(new_dV):
	
	global dV,freq,int
	
	dV = new_dV
	
	freq,int = run_sim(frequency,intensity,T,dV,specS)
	
	make_plots(freq,int)
	
#update_specS updates the scaling, redoes the simulation, and spits everything back out again.

def update_specS(new_specS):
	
	global specS,freq,int 
	
	specS = new_specS
	
	freq,int = run_sim(frequency,intensity,T,dV,specS)
	
	make_plots(freq,int)
	
#update_vlsr updates the vlsr, redoes the simulation, and spits everything back out again.

def update_vlsr(new_vlsr):
	
	global vlsr,freq,int,frequency
	
	vlsr = new_vlsr
	
	freq_temp = np.copy(frequency)
	
	freq_temp += (-vlsr)*freq_temp/ckm
	
	freq,int = run_sim(freq_temp,intensity,T,dV,specS)
	
	make_plots(freq,int)			
	
	

#############################################################
#							Run Program						#
#############################################################

#start_time = time.time()

raw_array = read_cat(catalog_file)

if args.ll or args.ul:

	raw_array = trim_raw_array(raw_array)

catalog = splice_array(raw_array)

frequency = np.empty_like(catalog[0])

for x in range(len(frequency)):

	frequency[x] = catalog[0][x]

logint = catalog[2]
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

if args.lsr:

	frequency += (-vlsr)*frequency/ckm

if args.spec:

	obs = read_cat(spec)
	
	freq_obs = []
	int_obs = []
	
	for x in range(len(obs)):
	
		freq_obs.append(obs[x].split()[0])
		int_obs.append(obs[x].split()[1].strip('\n'))	

freq,int=run_sim(frequency,intensity,T,dV,specS)
		
# freq = frequency		
# 
# if args.T:
# 
# 	int = scale_temp(int,qns,elower,qn7,qn8,qn9,qn10,qn11,qn12,T,CT)
# 	
# if args.specS:
# 
# 	int *= specS
# 	
# if args.gauss:
# 
# 	freq,int = sim_gaussian(int,frequency,dV)

if args.write:

	write_spectrum(freq,int,output_file)

make_plots(freq,int)
	
print('For a list of available commands, type \'help.\'')

quit_con = False

while quit_con == False:

	choice = input('> ')
	
	choice = choice.lower()
	
	quitlist = ['q','quit','e','exit']
	helplist = ['h', 'help']
	
	if any(choice == x for x in helplist):
	
		print('The program is now in interactive mode.  The following options are available:\n')
		
		print('S X:\t Update the scaling factor (specS) to X.  Re-draws the graph.\n')
		
		print('T X:\t Update the temperature (T) to X.  Re-draws the graph.\n')
		
		print('V X:\t Update the VLSR (vlsr) to X.  Re-draws the graph.\n')
		
		print('dV X:\t Update the linewidth to X. Re-draws the graph.\n')
		
		print('G:\t Re-draws the graph if accidentally closed.\n')
		
		print('W X:\t Exports the current simulation to an ascii file named X.\n')
		
		print('q:\t Quits the program.\n')
		
		print('h:\t Displays this message.\n')
	
	elif any(choice == x for x in quitlist):
	
		quit_con = True
		
	elif 't' in choice[0]:
	
		plt.close()
	
		try:
			new_temp = float(choice.split()[1])
		except IndexError:
			try:
				new_temp = float(choice[1:])
			except:
				print('Unrecognized command.  Type \'help\' for a list of commands.')
		
		update_T(new_temp)
		
	elif 'dv' in choice:
	
		plt.close()
	
		try:
			new_dV = float(choice.split()[1])
		except IndexError:
			try:
				new_dV = float(choice[1:])
			except:
				print('Unrecognized command.  Type \'help\' for a list of commands.')
		
		update_dV(new_dV)
		
	elif 's' in choice[0]:
	
		plt.close()
		
		try:
			new_specS = float(choice.split()[1])
		except IndexError:
			try:
				new_specS = float(choice[1:])
			except:
				print('Unrecognized command.  Type \'help\' for a list of commands.')
		
		update_specS(new_specS)
		
	elif 'v' in choice[0]:
	
		plt.close()
	
		try:
			new_vlsr = float(choice.split()[1])
		except IndexError:
			try:
				new_vlsr = float(choice[1:])
			except:
				print('Unrecognized command.  Type \'help\' for a list of commands.')
		
		update_vlsr(new_vlsr)		
		
	elif 'g' in choice[0]:
	
		plt.close()
		
		make_plots(freq,int)	
		
	elif 'w' in choice[0]:
	
		try:
			output_file = str(choice.split()[1])
			write_spectrum(freq,int,output_file)
		except:
			print('Sorry, I couldn\'t quite parse that.  Please remember a space between the \'w\' and the desired output file name.')
			
	else:
	
		print('Unrecognized command.  Type "help" for a list of commands.')			
			
#ttotal = end_time - start_time

#print('\nTook {:.1f} s to run program.' .format(ttotal))