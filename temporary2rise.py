import matplotlib.pyplot as plt
from scipy import optimize
from lmfit import Parameters,minimize, fit_report
import numpy as np

#x=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#y=[0,1,4,9,16,25,36,49,64,81,100,121,144,169,196,225,256,289,324,361,400]
x=[186.0105, 186.031, 186.0514]
x1=[]
for i in range(len(x)):
	x1.append(x[i]-x[0])
y=[0.99909, 1.00315, 1.00852]
e=2.7182818


def fit_exp(x,I,a,T):

		#y = parameters['A1']*sin(2*pi*parameters['a1']*x+paramters['b1']) + paramters['A2']*sin(2*pi*paramters['a2']*x+paramters['b2']) + Z
	return I*(e**((-(a/10000)*x)))+T

params, params_covariance = optimize.curve_fit(fit_exp, x1, y,
		                                                       p0=[2, 2, 2], maxfev = 100000)
		                                                     
print(params[0], params[1], params[2])
y_fit = []
x2 = np.arange(x1[0],x1[len(x1)-1],x1[1]/10)
for i in x2:	                                                     
	y_value = fit_exp(i,params[0],params[1],params[2])
	y_fit.append(y_value)
	

plt.plot(x1,y,"o")
plt.plot(x2,y_fit)


plt.show()
