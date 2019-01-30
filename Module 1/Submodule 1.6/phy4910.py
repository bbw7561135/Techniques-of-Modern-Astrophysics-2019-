import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

"""rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size':14})
plt.rc('axes', labelsize=16)
plt.rcParams.update({'figure.autolayout': True})"""


def plot(x, f, xlabel="", ylabel=""):
	"""
	This function will help plot an list of arrays. 
	plot(x, f, xlabel="", ylabel="")
	Where x is the independent variable. f is a list of dependent variables. Using this function you can plot multiple dependent variables  	against x
	"""
	fig = plt.figure(figsize=(6, 4))
	ax = fig.add_subplot(1,1,1)
	ax.grid(True)
	for i in range(len(f)):     
		ax.plot(x, f[i])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	plt.show()
    

def plot_print(x, f, filename, xlabel="", ylabel=""):
	"""
	This function will help plot an list of arrays. 
	plot_print(x, f, filename, xlabel="", ylabel="")
	Where x is the independent variable. f is a list of dependent variables. Using this function you can plot multiple dependent variables  	against x. Filename is the name of the saved pdf file
	"""
	fig = plt.figure(figsize=(6, 4))
	ax = fig.add_subplot(1,1,1)
	for i in range(len(f)):     
		ax.plot(x, f[i])    
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	plt.savefig(filename)
    

def ode_euler(x_start, x_end, h, y0, z0, f, g):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y = np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
            
        y[i+1] = y[i] + k1
        z[i+1] = z[i] + l1
                
    return x, y, z
    
def ode_rk4(x_start, x_end, h, y0, z0, f, g):
    
    x = np.arange(x_start, x_end, h)
    N = len(x)
    y= np.zeros(N)
    y[0] = y0
    z = np.zeros(N)
    z[0] = z0
    
    for i in range(0, N-1):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
        k2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1)
        l2 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k1, z[i] + 0.5 * l1)
        k3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2)
        l3 = h * g(x[i] + 0.5 * h, y[i] + 0.5 * k2, z[i] + 0.5 * l2)
        k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3)
        l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3)
    
        y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        z[i+1] = z[i] + (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0
        
    return x, y, z

if __name__=="__main__":
	x=np.linspace(0,10)
	f=x**2
	#print (x,f)
	#a = np.exp(x)
	#print (f.shape, a.shape)	
	#plot(x,[f,a])	
	#print(plot.__doc__)
	plot_print(x,[f], 'test1.pdf')	
	print(plot_print.__doc__)
