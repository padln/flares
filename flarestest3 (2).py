###########################################################################################
###########################################################################################
###                                                                                     ###
###                             ~~~Cool Code for Finding Flares~~~                      ###
###                                       ~~(CCFF)~~                                    ###
###                                                                                     ###
###########################################################################################
###########################################################################################

# here we import the light curve data as well as import some functions and other misc.
import numpy as np
import numpy as numpy
from scipy.signal import lombscargle
import matplotlib.pyplot as plt
import scipy.signal as signal
from astropy.timeseries import LombScargle
import scipy
from scipy.signal import find_peaks
from numpy import sin
from scipy import optimize
from multiprocessing import Pool
import statistics
import sys
import math
from scipy import interpolate

pi = 3.14159265358979323846

plt.ion()
fig, (ax1, ax2) = plt.subplots(2)



xu1, yu1 = np.loadtxt("KIC_6791060_all_jump_removed.dat", unpack = True, usecols = [0,2])


frequency, power = LombScargle(xu1, yu1).autopower()

periods = 1/frequency

#        here we graph the data and the LombScarge (this will probbably be deleted later)
#plt.plot(periods, power)
#plt.plot(xu1, yu1)
#plt.show()

# nyquist_frequency is frequency/2
#######################################################################################



#######################################################################################
#  here we calculate the sampling frequency and nyquist frequency and period
sampling_frequency =  len(xu1)/(xu1[len(xu1)-1]-xu1[0])
nyquist_frequency = sampling_frequency/2
nyquist_period = 1/nyquist_frequency

#print ("nyquist_period:", nyquist_period)
#print ("nyquist_frequency:", nyquist_frequency)
#######################################################################################



#######################################################################################
#program for finding peaks and making sure they are large enough to be significant
all_peaks = []
periods_counter = 0
for periods_values in range(len(periods)-1):
	periods_counter += 1
	if ((((power[periods_counter] - power[periods_counter - 1]) > 0.2*power[periods_counter-1] ) \
	or ((power[periods_counter] - power[periods_counter - 1]) > 0.2*power[periods_counter-1] )) \
	and ((power[periods_counter] > power[periods_counter - 1]) \
	and (power[periods_counter] > power[periods_counter + 1])) \
	and (power[periods_counter] > 0.01)):
		all_peaks.append(periods[periods_counter])
#######################################################################################



#######################################################################################
#program to ensure peaks are far enough apart to be distinct periods
peaks_counter = 0
distinct_peaks = []
for peaks in range(len(all_peaks)-2):
	peaks_counter += 1
	if  (all_peaks[peaks_counter] - all_peaks[peaks_counter+1]) > int(0.005):
		distinct_peaks.append(all_peaks[peaks_counter])

peaks_counter2 = 0
distinct_peaks2 = []
for peaks in range(len(distinct_peaks)-2):
	peaks_counter2 += 1
	if  (distinct_peaks[peaks_counter2] - distinct_peaks[peaks_counter2+1]) > int(0.005):
		distinct_peaks2.append(distinct_peaks[peaks_counter2])
#######################################################################################



#######################################################################################
# defining printing functions to check answers for error (printlist)
def printlist(lst):
	for a in range(len(lst)-1):
		print (lst[a])

def printlist2(lst , lst2):
	for a in range(len(lst)-1):
		print (lst[a], lst2[a])

def printlistnumbered(lst):
	n = 1
	for a in range(len(lst)-1):
		print (n, ")" ," " , lst[a], sep = '')
		n +=1

#print ("all peak:")
#print (printlist(all_peaks))
#print (" ")
#print ("distinct peaks:")
#print (printlist(distinct_peaks2))
#######################################################################################



#######################################################################################
#function that takes the average value of a list ######################################
def average(lst):
	total = 0
	for i in lst:
		total += i
	average = total/ len(lst)
	return average
#######################################################################################



#######################################################################################
#program that stores that parameters of the function in a dynamic parameters dictionary
parameters_guess = {}           #used to store the parameters of the sine function


parameter_index = 0

for peak in distinct_peaks2:
	if (peak > nyquist_period):
		parameter_index = parameter_index + 1
		#using form An*sin(2*pi*an*x+bn)
		parameters_guess['A' + str(parameter_index)]='1'
		parameters_guess['a' + str(parameter_index)]= peak
		parameters_guess['b' + str(parameter_index)]='1'


print ("parameters guess: ", parameters_guess)
#######################################################################################



#######################################################################################
# function to calculate the standard deviation of a list of deviations
def standard_deviation(lst):
	total = 0
	for p in lst:
		total += p**2
	if len(lst) > 0:
		sd = (total / len(lst))**(0.5)
		#print ("total:", total)
		#print ("n=",  len(lst))

		return sd
#######################################################################################



#######################################################################################
#fits sin equations to the data in sets of *length* data points and deviation and sd###
length = 40     #can be adjusted for different fittings,length of running_sd calculator
x_value = 2650
initial = x_value #this is where the program starts in the data set
running_sd = []

for x_values in range(len(xu1)-64800):
	ax2.clear()
	ax1.clear()
	def test_func(x, A, b, A2, b2, Z):

		#y = parameters['A1']*sin(2*pi*parameters['a1']*x+paramters['b1']) + paramters['A2']*sin(2*pi*paramters['a2']*x+paramters['b2']) + Z
		return A*sin(2*pi*(1/parameters_guess['a1'])*x+b) + A2*sin(2*pi*(1/parameters_guess['a2'])*x+b2) + Z

	params, params_covariance = optimize.curve_fit(test_func, xu1[x_value:x_value+length], yu1[x_value:x_value+length],
		                                                       p0=[2, 2, 2, 2, 2], maxfev = 100000)
	ax1.scatter(xu1, yu1, label='Data', color = 'green')
	ax1.plot(xu1[x_value: x_value+length], test_func(xu1[x_value:x_value+length], params[0], params[1], params[2], params[3], params[4]), label='Fitted function', color = "blue")
	deviation = []
	for point in range(x_value, x_value+length):
		deviation.append( yu1[point] - ( params[0]*sin(2*pi*(1/parameters_guess['a1'])*point+params[1]) + params[2]*sin(2*pi*(1/parameters_guess['a2'])*point+params[3]) + params[4] ) )
	running_sd.append( standard_deviation(deviation) )
	x_value+=(1)
	ax1.set_xlim([xu1[x_value]-0.5, xu1[x_value+length]+0.5])
	#print (x_value-2743, len(running_sd))
	ax2.scatter(xu1[0:(x_value - initial)] , running_sd, 2, color = 'black', linewidth = 2)
	ax2.axhline(y=statistics.median(running_sd))
	running_sd_deviation = []
	for i in running_sd:
		running_sd_deviation.append(abs(i-average(running_sd)))
	ax2.axhline(y=statistics.median(running_sd) + 2.5* standard_deviation(running_sd_deviation), color = 'red')
	ax2.axhline(y=statistics.median(running_sd) - 2.5* standard_deviation(running_sd_deviation), color = 'red')
	#print (3* standard_deviation(running_sd))
	fig.canvas.draw()
	fig.canvas.flush_events()

#######################################################################################



#######################################################################################
#running standard deviation dependant on *length* and flares###########################

flares_not_3 = []
#^list of points that pass the standard deviation test but
#may not necesarrily be long enough to be a flare
flares_starttime = []
flares_starttime_index = []
flares_endtime = []
flares_endtime_index = []
flare_lengths = []
flares_dataPoints_x = [] #this is a list containing a list of all the data points in a flare
flares_dataPoints_y = [] #this is a list containing a list of all the data points in a flare
i=0
for k in range(len(running_sd)-length):
	if running_sd[i] > average(running_sd) + 2.5*standard_deviation(running_sd_deviation):
		#^this ensures that the flares are actually noteworthy (pass the 2.5sigma threshold)
		flares_not_3.append(xu1[i])
		for p in range(20,0,-1): #here the code has found a point that meets the 2.5sigma threshold
		#then the code looks 20 points ahead to see if that point meets the threshold too
		#if that point also meets the threshold it records it as the flares_endtime
		#if +20 does not meet the threshold it will look at +19, then +18, so on
			if running_sd[i+length+p-5] > average(running_sd) + 2.5*standard_deviation(running_sd_deviation):
				flares_starttime.append(xu1[i+initial+length-3])
				flares_starttime_index.append(i+initial+length-3)
				flares_endtime.append(xu1[i+p+initial+length])
				flares_endtime_index.append(i+p+initial+length)
				flare_lengths.append(p)
				flares_dataPoints_x.append([])
				flares_dataPoints_y.append([])
				for d in range(flares_starttime_index[len(flares_starttime_index)-1], flares_endtime_index[len(flares_endtime_index)-1]):
				#here we say "take every data point between a flares_starttime and flares_endtime_index"
				#and add them as a list to a list containing all the point that are in a flare
					flares_dataPoints_x[len(flares_dataPoints_x)-1].append(xu1[d])
					flares_dataPoints_y[len(flares_dataPoints_y)-1].append(yu1[d])
				i = i + length + p - 5
				break

	i = i+1
	#if running_sd[i] > average(running_sd) - 3*standard_deviation(running_sd):
	#	flares_not_3.append(xu1[i])
print ('flares_not_3', flares_not_3)
print ('flares_starttime', flares_starttime)
print ('flares_starttime_index', flares_starttime_index)
print ('flares_endtime', flares_endtime)
print ('flares_endtime_index', flares_endtime_index)
print ('flare_lengths',flare_lengths)
#######################################################################################



#######################################################################################
# here we will remove the flares from the data
xu3 = xu1
yu3 = yu1
for n in range(len(flares_starttime)):
	for i in range(flares_starttime_index[n], flares_endtime_index[n]):
		xu3 = np.delete(xu3, i)
		yu3 = np.delete(yu3, i)

print("xu1 Dimention = ",len(xu1))
print("xu3 Dimention = ",len(xu3))

#######################################################################################


'''THIS IS THE OLD METHOD: (I JUST LEFT IT HERE INCASE IT'S NEEDED LATER)
(IT DOESNT WORK, I DONT KNOW WHY)
#######################################################################################
# here we will use the new flareless data and interpolate it
f3 = interpolate.interp1d(xu3, yu3) #saves interpolated data as a function f3

xu3_new = np.arange(xu3[0], xu3[len(xu3)-2], 0.0001) #creates a whole bunch of x-values very close together

#gets y-values for those x-values to be used in reimann sum
yu3_new = f3(xu3_new)   # use interpolation function returned by `interp1d`

#ax1.plot(xu3, yu3, 'o', xu3_new, y_u3new, '-')

#######################################################################################



#######################################################################################
# here we interpolate the original lightcurve data with the flares
f1 = interpolate.interp1d(xu1, yu1)

xu1_new = np.arange(xu1[0], xu1[len(xu1)-2], 0.0001)

yu1_new = f1(xu1_new)


#######################################################################################



#######################################################################################
# here we calculate the area between the flareless curve and the flarecurve

energy = [] #this is a list that contains the area (proportional to energy) of each flare
for q in range(len(flares_starttime)):

	sector = flares_endtime[q] - flares_starttime[q]
	limit = 1000 #this is the number of dA's we are deviding the flare area into
	area = 0
	for p in range(limit): #the p is the counter varibale for the dA
		#print (f2(flares_starttime[q] + (p * sector) / limit))
		#print (f1(flares_starttime[q]+ (p * sector) / limit))
		#print ( f2(flares_starttime[q]+ (sector)/limit) - f1(flares_starttime[q]+ (sector)/ limit) )
		dA = (sector/limit)(f1[flares_starttime[q]+p]-f3[flares_starttime[q]+p])
		area += dA
		#print("dA =", dA)
	energy.append(area)
print ("energy", energy)
#######################################################################################
'''


#######################################################################################
#this gets us a list of the data points around a flare so we can fit a curve to it
#then using the fitted curve we can predict what the flareless behavior would be
points_around_flares_without_flares_x = []
points_around_flares_without_flares_y = []
for k in range(len(flares_starttime_index)):
	empty_list = []
	for i in range(flares_starttime_index[k]-(2*length),flares_starttime_index[k]-2):
		empty_list.append(xu1[i])
	for i in range(flares_endtime_index[k],flares_endtime_index[k]+(2*length)):
		empty_list.append(xu1[i])
	points_around_flares_without_flares_x.append(empty_list)
	
	empty_list = []
	for i in range(flares_starttime_index[k]-(2*length),flares_starttime_index[k]-2):
		empty_list.append(yu1[i])
	for i in range(flares_endtime_index[k],flares_endtime_index[k]+(2*length)):
		empty_list.append(yu1[i])	
	points_around_flares_without_flares_y.append(empty_list)
	
#######################################################################################



#######################################################################################
#defines a function that fits a sin graph to a curve, nothing fancy
def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = numpy.array(tt)
    yy = numpy.array(yy)
    ff = numpy.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(numpy.fft.fft(yy))
    guess_freq = abs(ff[numpy.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * numpy.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*numpy.pi)
    fitfunc = lambda t: A * numpy.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": numpy.max(pcov), "rawres": (guess,popt,pcov)}
   
#params[0]*sin(2*pi*(1/parameters_guess['a1'])*point+params[1]) + params[2]*sin(2*pi*(1/parameters_guess['a2'])*point+params[3]) + params[4]

#######################################################################################

#######################################################################################  
energy = [] 
flares_detrended_x = []
flares_detrended_y = []
flare_rises_dataset_x = []
flare_decays_dataset_x = []
flare_rises_dataset_y = []
flare_decays_dataset_y = []
flare_rises_dataset_y_detrended = []
flare_decays_dataset_y_detrended = []
for k in range(len(flares_starttime_index)):
	active_flare_dataset_x = flares_dataPoints_x[k]
	active_flare_dataset_y = flares_dataPoints_y[k]
	active_flareless_dataset_x = points_around_flares_without_flares_x[k]
	active_flareless_dataset_y = points_around_flares_without_flares_y[k]
	#############################
	# this fits a curve to the points around the flareless neighborhood
	def test_func(x, A, b, A2, b2, Z):

		#y = parameters['A1']*sin(2*pi*parameters['a1']*x+paramters['b1']) + paramters['A2']*sin(2*pi*paramters['a2']*x+paramters['b2']) + Z
		return A*sin(2*pi*(1/parameters_guess['a1'])*x+b) + A2*sin(2*pi*(1/parameters_guess['a2'])*x+b2) + Z

	params, params_covariance = optimize.curve_fit(test_func, active_flareless_dataset_x[0:len(active_flareless_dataset_x)], active_flareless_dataset_y[0:len(active_flareless_dataset_y)],
		                                                       p0=[2, 2, 2, 2, 2], maxfev = 100000)
	
	#############################
	#here we sum up dA to get the total area of each flare
	area = 0
	for t in range(len(active_flare_dataset_y)-1):
		dA = (   (active_flare_dataset_x[t+1]-active_flare_dataset_x[t])*( (active_flare_dataset_y[t]-test_func((t+1),params[0],params[1],params[2],params[3],params[4])) + (0.5)*(active_flareless_dataset_y[t+1]-active_flareless_dataset_y[t]) + (0.5)*(test_func((t+1),params[0],params[1],params[2],params[3],params[4])-test_func(t,params[0],params[1],params[2],params[3],params[4])) ) )
		area += dA
	energy.append(area)
	#############################
	oneflare_detrended_x = []
	oneflare_detrended_y = []
	for t in range(flares_starttime_index[k]-(2*length),flares_endtime_index[k]+(2*length)):
		detrend_value = test_func(xu1[t],params[0],params[1],params[2],params[3],params[4])
		oneflare_detrended_y.append( yu1[t] -  detrend_value )
		oneflare_detrended_x.append( xu1[t] )
	flares_detrended_x.append(oneflare_detrended_x)
	flares_detrended_y.append(oneflare_detrended_y)
	#############################
	#next we will split the flare into post peak and pre peak data 
	one_flare_rises_dataset_x = []
	one_flare_decays_dataset_x = []
	one_flare_rises_dataset_y = []
	one_flare_decays_dataset_y = []
	index = active_flare_dataset_y.index(max(active_flare_dataset_y))
	for t in range(0,index+1):
		one_flare_rises_dataset_y.append(active_flare_dataset_y[t])
		one_flare_rises_dataset_x.append(active_flare_dataset_x[t])
	for t in range(index, len(active_flare_dataset_y)):
		one_flare_decays_dataset_y.append(active_flare_dataset_y[t])
		one_flare_decays_dataset_x.append(active_flare_dataset_x[t])
	flare_rises_dataset_x.append(one_flare_rises_dataset_x)
	flare_decays_dataset_x.append(one_flare_decays_dataset_x)
	flare_rises_dataset_y.append(one_flare_rises_dataset_y)
	flare_decays_dataset_y.append(one_flare_decays_dataset_y)
	#############################
	one_flare_rises_dataset_y_detrended = []
	one_flare_decays_dataset_y_detrended = []
	for t in range(len(one_flare_rises_dataset_y)):
		detrend_value = test_func(xu1[t],params[0],params[1],params[2],params[3],params[4])
		one_flare_rises_dataset_y_detrended.append(one_flare_rises_dataset_y[t] - detrend_value)
	for t in range(len(one_flare_decays_dataset_y)):
		detrend_value = test_func(xu1[t],params[0],params[1],params[2],params[3],params[4])
		one_flare_decays_dataset_y_detrended.append(one_flare_decays_dataset_y[t] - detrend_value)
	flare_rises_dataset_y_detrended.append(one_flare_rises_dataset_y_detrended)	
	flare_decays_dataset_y_detrended.append(one_flare_decays_dataset_y_detrended)
#######################################################################################
print(flare_rises_dataset_y)
print(flare_decays_dataset_y)

print(flares_dataPoints_x[0])
print(flares_dataPoints_y[0])

#then we will fit an exponential line to that data to parameterize t_r and t_d




#original_stdout = sys.stdout # Save a reference to the original standard output

#with open('filename.txt', 'w') as f:
#    sys.stdout = f # Change the standard output to the file we created.
#    print(printlist2(running_sd, xu1[2700:len(xu1)-62400]))
#    sys.stdout = original_stdout # Reset the standard output to its original value

#plt.savefig("flares.png")

#print ("Running SD:")
#printlistnumbered(running_sd)

#ax2.scatter(xu1[0:945] , running_sd, 2)
