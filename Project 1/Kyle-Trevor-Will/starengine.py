from numpy import *
from matplotlib.pyplot import *
from IPython.display import clear_output



#Input
X = 0.70                                             #Percent mass of Hydrogen
Y = 0.28                                             #Percent mass of Helium
Z = 1-X-Y                                            #Percent mass of Metals

M = 1.00                                             #Mass of star in solar masses

# Constants
G = 6.67408e-11                                      #m^3 * kg^-1 * s^-2  Gravitational constant
a = 7.565767e-16                                     #J * m^-3 * K^-4     Radiation density constant
c = 2.99792458e8                                     #m * s^-1            Speed of light
k = 1.38064852e-23                                   #J * K^-1            Boltzmann Constant
sigma = 5.670367e-8                                  #W * m^-2 * K^-4     Stefan-Boltzmann Constant
mH = 1.673532499e-27                                 #kg                  Mass of Hydrogen
gamma = 5.0 / 3.0                                    #                    Ideal Gas

mu_i  = (2 * X + (3.0/4.0) * Y + (1.0/2.0) * Z) ** -1   #Assuming complete ionization this equation approximates mu_i for standard stellar properties

r0 = 6.957e8                                                   #m                solar radius
M0 = 1.98847e30                                                #kg               solar mass
L0 = 3.828e26                                                  #w                solar luminosity
rho0 = M0 / (4 * pi *r0 ** 3)                                  #kg * m^-3        density
P0 = G * M0 * rho0 / r0                                        #Pa               pressure
eps0 = L0 / (4 * pi * rho0 * r0 ** 3)                          #w * kg^-1        energy from nuclear reactions
T0 = (1.0 - 1.0 / gamma) * (mu_i * mH / k) * (G * M0 / r0)     #K                temperature
kap0 = (16 * pi * a * c * r0 * T0 ** 4) / (3 * rho0 * L0)      #m^2 * kg^-1      opacity

#print(format(rho0,'e'))                                       #used to check a constant's value                  


#Defining Fucntions to be integrated

def f (eta, m, l, t, vrho, kappa, e):                               #Equation 1       d(p)/d(eta)
    return -m * vrho / eta ** 2

def g (eta, m, l, t, vrho, kappa, e):                               #Equation 2       d(m)/d(eta)
    return vrho * eta ** 2

def h (eta, m, l, t, vrho, kappa, e):                               #Equation 3       d(l)/d(eta)
    return vrho * e * eta ** 2

def i (eta, m, l, t, vrho, kappa, e):                               #Equation 4A      d(t)/d(eta) radiation
    return -(kappa * vrho * l) / (t ** 3 * eta ** 2)

def j (eta, m, l, t, vrho, kappa, e):                               #Equation 4B      d(t)/d(eta) convection
    return -m / eta ** 2


#Functions for parameters

#This function should probably removed
#It produces values equivalent to the alternat Vrho function.
#Initially we wanted to use a dimensionless equation but it was much more intuitive to nondimensionalize afterwards
# def getVrho (p,t):                                   #Equation 5       vrho
#     return (gamma / (gamma - 1)) * p / t  - 4.0 * pi / 3  * (1 - 1 / gamma) ** 3  * (mu_i * mH / k) ** 4  * (a * G ** 3 * M0 ** 2) * t ** 3  

def getVrho (p,t):                                                  #Equation 5     vrho
    
    """
    p   - dimensionless pressure
    t   - dimensionless temperature
    
    """
    
    
    T = T0 * t                                        #Dimensionalize the variables
    P = P0 * p
    Prad = 1.0 / 3.0 * a * T ** 4                     #Radiation from pressure
    rho = (P - Prad) * mu_i * mH / (k * T)            #Calculate density with units
    
    return rho / rho0                                 #Nondimensionalize the density
    
    
def getE (vrho, t):                       #Equation 6       energy from nuclear reactions
    
    """
    vrho   - dimensionless density
    t      - dimensionless temperature
    
    """
    
    
    #Calculate energy from Proton-Proton(pp) chains
    T6 = T0 * t / 10 ** 6                                                                #Temperature in units of 10^6K
    psiPP = 1 + 1.412e8 * (1 / X - 1) * exp(-49.98 * T6**(-1 / 3))                       #correction factor for simultaneous reactions of pp chains
    cPP = 1 + 0.0123 * T6 ** (1.0/3.0) + 0.0109 * T6 ** (2.0/3.0) + 0.000938 * T6        #another correction factor
    fPP = 1                                                                              #pp chain screeing factor
    epsPP = 0.241 * rho0 * vrho * X ** 2 * fPP * psiPP * cPP * T6 ** (-2.0/3.0) * exp(-33.80 * T6 ** (-1.0/3.0)) #total energy from pp chains
    
    #Calculate energy from Carbon-Nitrogen-Oxygen(CNO) cycle
    XCNO = Z / 2.0                                                                       #mass fraction of CNO, approximately half of the metals contribute to CNO
    cCNO =  1 + 0.0027 * T6 ** (1.0/3.0) - 0.00778 * T6 ** (2.0/3.0) - 0.000149 * T6     #correction factor for CNO
    epsCNO = 8.67e20 * rho0 * vrho * X * XCNO * cCNO * T6 ** (-2.0/3.0) * exp(-152.28 * T6 ** (-1.0/3.0))  #total energy from CNO cycle
    
    #Calculate energy from triple alpha(3alpha) process 
    T8 = T0 * t / 10 ** 8                                                                #Temperature in units of 10^8K
    f3alpha = 1.0                                                                        #3alpha screening factor
    eps3alpha = 50.9 * (rho0 * vrho) ** 2 * Y ** 3  * T8 ** -3  * f3alpha * exp(-44.027 * T8 ** -1) #total energy from 3alpha process
    
    return (epsPP + epsCNO + eps3alpha) / eps0                                           #return the sum of energies and nondimensionalize


def getKappa (vrho, t):                      #Equation 7       Opacities
    
    """
    vrho  - dimensionless density
    t     - dimensionless temperature
    """
    
    
    gbf_t = (0.708 * (rho0 * vrho * (1+X)) ** (1.0/5.0)) ** -1                           #Gaunt factor for bound-free process divided by guillotine factor
    kbf = 4.34e21 * gbf_t * Z * (1 + X) * rho0 * vrho / (T0 * t) ** 3.5                  #opacity from bound-free process
    
    gff = 1.0                                                                            #Gaunt factor for free-free process
    kff = 3.68e18 * gff * (1 - Z) * (1 + X) * rho0 * vrho / (T0 * t) ** 3.5              #opacity for free-free process
    
    kes = 0.02 * (1 + X)                                                                 #opacity for electron scattering
    
    #Calculate opacity due to Hydrogen ions
    if T0 * t >= 3000 and T0 * t <= 6000 and rho0 * vrho >= 10**-7 and rho0 * vrho <= 10**-2 and Z >= 0.001 and Z <= 0.03:  #Conditions for H- ions to play a role
        kH = 7.9e-34 * ( Z / 0.02) * (rho0 * vrho) ** (1.0/2.0) * (T0 * t) ** 9
    else:
        kH = 0
    
    return (kbf + kff + kes + kH) / (kap0)                                               #return the sum of opacities and nondimensionalize


#Integration  using 4th order runge kutta method
def integrate(eta_start, eta_end, delta, m0, l0, p0, t0, f, g, h, i, j, stop = lambda eta, m, p, l, t, vrho, kappa, e,data:False , data = None):
    """
    eta_start - starting value of radius for integration
    eta_end   - final value of radius for integration
    delta     - step size for integration
    m0        - initial dimensionless mass
    l0        - initial dimensionless luminosity
    p0        - initial dimensionless pressure
    t0        - initial dimensionless temperature
    f         - Equation 1       d(p)/d(eta)
    g         - Equation 2       d(m)/d(eta)
    h         - Equation 3       d(l)/d(eta)
    i         - Equation 4A      d(t)/d(eta) radiation
    j         - Equation 4B      d(t)/d(eta) convection
    stop      - Function to stop integration - False by default
    data      - data to be passed to stop function - none by default
    """



    eta = arange(eta_start, eta_end, delta)                                   #Range of integration
    N = len(eta)                                                              #Number of steps
    
    m = zeros(N)                                                              #Array for mass
    m[0] = m0
    p = zeros(N)                                                              #Array for pressure
    p[0] = p0
    l = zeros(N)                                                              #Array for luminosity
    l[0] = l0
    t = zeros(N)                                                              #Array for temperature
    t[0] = t0
    vrho = zeros(N)                                                           #Array for density
    vrho[0]= getVrho(p[0],t[0])

    dp = 0                                                                  #initialize change in pressure
    dt = 0                                                                  #initialize change in temperature
    
    for num in range(0,N-1):
        
        if num > 0 and t[num] * dp / (p[num]*dt) < gamma / (gamma-1):       #if the condition is true convection occurs
            Radiation = False
            #print("Convective!")
        else:
            Radiation  = True
        
        vrhoTemp = vrho[num]                                                #used in calling the functions
        eTemp = getE(vrhoTemp, t[num])
        kappaTemp = getKappa(vrhoTemp, t[num])
        
        #print(num,1, vrhoTemp)
            
        if stop(eta[num], m[num], p[num], l[num], t[num], vrhoTemp, kappaTemp, eTemp, data) == True:  #stop function
            #print('stopped')
            return eta[0:num], m[0:num], p[0:num], l[0:num], t[0:num], vrho[0:num]
        
        
        
        
        k1 = delta * f(eta[num], m[num], l[num], t[num], vrhoTemp, kappaTemp, eTemp)                        #change in pressure
        l1 = delta * g(eta[num], m[num], l[num], t[num], vrhoTemp, kappaTemp, eTemp)                        #change in mass
        m1 = delta * h(eta[num], m[num], l[num], t[num], vrhoTemp, kappaTemp, eTemp)                        #change in luminosity
        
        
        if Radiation:                                                                        #change in temperature depends on radiation
            n1 = delta * i(eta[num], m[num], l[num], t[num], vrhoTemp, kappaTemp, eTemp)
        else:
            n1 = delta * j(eta[num], m[num], l[num], t[num], vrhoTemp, kappaTemp, eTemp)
           
        vrhoTemp = getVrho(p[num] + 0.5 * k1, t[num] + 0.5 * n1)                             #update these variables before calling the functions
        eTemp = getE(vrhoTemp, t[num] + 0.5 * n1)
        kappaTemp = getKappa(vrhoTemp, t[num] + 0.5 * n1)

        
        k2 = delta * f(eta[num] + 0.5 * delta, m[num] + 0.5 * l1, l[num] + 0.5 * m1, t[num] + 0.5 * n1, vrhoTemp, kappaTemp, eTemp)                        #change in pressure
        l2 = delta * g(eta[num] + 0.5 * delta, m[num] + 0.5 * l1, l[num] + 0.5 * m1, t[num] + 0.5 * n1, vrhoTemp, kappaTemp, eTemp)                        #change in mass
        m2 = delta * h(eta[num] + 0.5 * delta, m[num] + 0.5 * l1, l[num] + 0.5 * m1, t[num] + 0.5 * n1, vrhoTemp, kappaTemp, eTemp)                           #change in luminosity
                      
        if Radiation:                                                                         #change in temperature depends on radiation
            n2 = delta * i(eta[num] + 0.5 * delta, m[num] + 0.5 * l1, l[num] + 0.5 * m1, t[num] + 0.5 * n1, vrhoTemp, kappaTemp, eTemp)
        else:
            n2 = delta * j(eta[num] + 0.5 * delta, m[num] + 0.5 * l1, l[num] + 0.5 * m1, t[num] + 0.5 * n1, vrhoTemp, kappaTemp, eTemp)
  
        vrhoTemp = getVrho(p[num] + 0.5 * k2, t[num] + 0.5 * n2)                              #update these variables before calling the functions
        eTemp = getE(vrhoTemp, t[num] + 0.5 * n2)
        kappaTemp = getKappa(vrhoTemp, t[num] + 0.5 * n2)
        
        
        k3 = delta * f(eta[num] + 0.5 * delta, m[num] + 0.5 * l2, l[num] + 0.5 * m2, t[num] + 0.5 * n2, vrhoTemp, kappaTemp, eTemp)                        #change in pressure
        l3 = delta * g(eta[num] + 0.5 * delta, m[num] + 0.5 * l2, l[num] + 0.5 * m2, t[num] + 0.5 * n2, vrhoTemp, kappaTemp, eTemp)                        #change in mass
        m3 = delta * h(eta[num] + 0.5 * delta, m[num] + 0.5 * l2, l[num] + 0.5 * m2, t[num] + 0.5 * n2, vrhoTemp, kappaTemp, eTemp)                           #change in luminosity
                
        if Radiation:                                                                          #change in temperature depends on radiation
            n3 = delta * i(eta[num] + 0.5 * delta, m[num] + 0.5 * l2, l[num] + 0.5 * m2, t[num] + 0.5 * n2, vrhoTemp, kappaTemp, eTemp)
        else:
            n3 = delta * j(eta[num] + 0.5 * delta, m[num] + 0.5 * l2, l[num] + 0.5 * m2, t[num] + 0.5 * n2, vrhoTemp, kappaTemp, eTemp)

        vrhoTemp = getVrho(p[num] + k3, t[num] + n3)                                           #update these variables before calling the functions
        eTemp = getE(vrhoTemp, t[num] + n3)
        kappaTemp = getKappa(vrhoTemp, t[num] + n3)
        
        k4 = delta * f(eta[num] + delta, m[num] + l3, l[num] + m3, t[num] + n3, vrhoTemp, kappaTemp, eTemp)                        #change in pressure
        l4 = delta * g(eta[num] + delta, m[num] + l3, l[num] + m3, t[num] + n3, vrhoTemp, kappaTemp, eTemp)                        #change in mass
        m4 = delta * h(eta[num] + delta, m[num] + l3, l[num] + m3, t[num] + n3, vrhoTemp, kappaTemp, eTemp)                           #change in luminosity
 
        if Radiation:                                                                           #change in temperature depends on radiation
            n4 = delta * i(eta[num] + delta, m[num] + l3, l[num] + m3, t[num] + n3, vrhoTemp, kappaTemp, eTemp)
        else:
            n4 = delta * i(eta[num] + delta, m[num] + l3, l[num] + m3, t[num] + n3, vrhoTemp, kappaTemp, eTemp)
        
     
        dp = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0                     #change in pressure
        dm = (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0                     #change in mass
        dl = (m1 + 2.0 * m2 + 2.0 * m3 + m4) / 6.0                     #change in luminosity
        dt = (n1 + 2.0 * n2 + 2.0 * n3 + n4) / 6.0                     #change in temperature
        
        
        p[num+1] = p[num] + dp                                         #generate new values
        m[num+1] = m[num] + dm
        l[num+1] = l[num] + dl
        t[num+1] = t[num] + dt
        vrho[num+1] = getVrho(p[num+1], t[num+1])

            
    return eta, m, p, l, t, vrho


def IsNaN(number):  
    """
    this function will determine if an output is not a number
    number - either a number or not a number
    
    """
    return number != number


def STOP (eta, m, p, l, t, vrho, kappa, e,data):
    
    """
    This function will stop integration if any of the values explode and return nan
    
    eta       - dimensionless radius
    m         - dimensionless mass
    l         - dimensionless luminosity
    p         - dimensionless pressure
    t         - dimensionless temperature
    vrho      - dimensionless density
    kappa     - dimensionless opacity
    e         - dimensionless energy
    data      - unused
    """
    
    if IsNaN(vrho) or IsNaN(p) or IsNaN(m) or IsNaN(l) or IsNaN(t):
        print('is nan')
        return True
    else:
        return False

    
    
#Below are two functions to find the specific t0 value that will create the necessary boundary conditions for our
#star model given some p0 value. The integration of these four coupled equations creates a very irregular pattern
#in the the final results and the bisection method (Findt0) often will not work to find the correct value.
#A second method used was to try every possible digit in decreasing decimal places to see if the sign changes
#This method also did not work very often.

def Findt0 (t1,t2, eta_start, eta_end, delta, m0, l0, p0, f,g,h,i,j,stop = STOP):       #find t0 through bisection method
    
    """
    This function will find t0 through the bisection method
    
    t1        - lower bound guess at where t0 really is
    t2        - upper bound guess at where t0 really is
    eta_start - starting value of radius for integration
    eta_end   - final value of radius for integration
    delta     - step size for integration
    m0        - initial dimensionless mass
    l0        - initial dimensionless luminosity
    p0        - initial dimensionless pressure
    t0        - initial dimensionless temperature
    f         - Equation 1       d(p)/d(eta)
    g         - Equation 2       d(m)/d(eta)
    h         - Equation 3       d(l)/d(eta)
    i         - Equation 4A      d(t)/d(eta) radiation
    j         - Equation 4B      d(t)/d(eta) convection
    stop      - Function to stop integration - STOP by default  
    """
    
    
    tmid = t1+(t2-t1)/2                                                                  #pick a value between t1 and t2
    
    print('Processing Temperature value, please wait...')
    
    model = integrate(eta_start, eta_end, delta, m0, l0, p0, tmid,f,g,h,i,j,stop = STOP)
    t_fmid = model[4][-1]                                                                #find final temperature when tmid is used
    
    while  abs(t_fmid) > 10**-2:                                                         #controls the accuracy
        t_f1 = integrate(eta_start, eta_end, delta, m0, l0, p0, t1,f,g,h,i,j,stop = STOP)[4][-1]     #final temp when t1 is used
        t_f2 = integrate(eta_start, eta_end, delta, m0, l0, p0, t2,f,g,h,i,j,stop = STOP)[4][-1]     #final temp when t2 is used
        
        tmid = t1+ (t2-t1)/2                                                             #pick a new value between t1 and t2               
        
        if t_f1 * t_f2 < 0:# or t_f1 == t1 and t_f2 != t2:                               #if the true value is between t1 and t2 the two final temps will have opposite sign
                                                                                         #there is a specific problem however where if t1 is too low or t2 is too high the integration
                                                                                         #will not be carried out correctly, so some knowledge is required for initial guess
            tmid = t1+ (t2-t1)/2                                                         #pick a new value between t1 and t2                    
            clear_output()
            print('t1 = ',t1, 't2 = ',t2)
            
        else:
            print('not in range')                                                        #the true value was not in between these too (or the initial guess was too far off)
            return(tmid)
        
        t_fmid = integrate(eta_start, eta_end, delta, m0, l0, p0, tmid,f,g,h,i,j,stop = STOP)[4][-1]    #find final temp with tmid      
        if t_fmid * t_f1 < 0:                                     #if the true value is on the left half move the right boundary in
            t2 = tmid
        else:                                                     #if the true value is on the right half move the left boundary in
            t1 = tmid
    return(tmid)



def Findt02(N, eta_start, eta_end, delta, m0, l0, p0, t0, f,g,h,i,j,stop = STOP):        
    """
    An alternate procedure to find some information about where the true value of t0 is
    A second method that will try every possible digit in decreasing decimal places to see if the sign of the final
    t value changes.
    
    
    N         - 10^N will be the first decimal place checked
    eta_start - starting value of radius for integration
    eta_end   - final value of radius for integration
    delta     - step size for integration
    m0        - initial dimensionless mass
    l0        - initial dimensionless luminosity
    p0        - initial dimensionless pressure
    t0        - initial dimensionless temperature
    f         - Equation 1       d(p)/d(eta)
    g         - Equation 2       d(m)/d(eta)
    h         - Equation 3       d(l)/d(eta)
    i         - Equation 4A      d(t)/d(eta) radiation
    j         - Equation 4B      d(t)/d(eta) convection
    stop      - Function to stop integration - STOP by default  
    """

    
    
    
    powers10 = zeros(12+N)                                                           #will be a list of powers of 10 decreasing
    digits = list(range(10))                                                         #a list of all digits
    for b in range(12+N):                                                            #generate the powers of 10 
        powers10[b] = 10**(N-b)                                                      #N selects the highest power
#     model = integrate(eta_start, eta_end, delta, m0, l0, p0, t0,f,g,h,i,j,stop = STOP) #find final temp when the initial guess t0 is used
#     t_f = model[4][-1]
    ttemp = t0                                                                       
    for b in powers10:                                                             #for each power of 10
        for dig in digits:                                                         #for each digit
            model = integrate(eta_start, eta_end, delta, m0, l0, p0, ttemp + b*dig,f,g,h,i,j,stop = STOP)  #find the final temp when the current
                                                                                                           #guess + some digit in some decimal place is used
            t_f = model[4][-1]
            if t_f > 0 and t_f != ttemp + b*dig:                         #the sign for t_f should initially be negative, it switches when the value is too high
                print (t_f)                                              #if the guess is extremely low the integration is not carried out the second condition ignores this integration
                ttemp += b*(dig-1)                                       #add the previous used digit in the current decimal place
                print (b,dig,ttemp)
                break                                                    #don't attempt the remaining digits
            elif dig == 9 and t_f != ttemp + b*dig:                      #if none of the digits switched the sign then it must be 9
                print (t_f)
                ttemp += b*(dig)
                print (b,dig,ttemp)
    print('ttemp =',ttemp)
    return ttemp
    




#set p0 = some value, then try to find t0 that will have both p = 0 and t = 0 at the same radius.

eta_start = 0.001                                                      #starting radius for integration
eta_end = 5                                                            #end radius for integration
delta = 0.001                                                          #step size for integration

#for boundary conditions
m0 = 0.001                                                             #initial dimensionless mass for integration
l0 = 0.001                                                             #initial dimensionless luminosity for integration
p0 = 4.586e3                                                           #initial dimensionless pressure for integration, vary this for different models

#Switch between the next three lines to change how t0 is determined

t0 = 4.06363084                                                         #just assign a number and perform integration only once
#t0 =Findt0(4.0,4.1, eta_start, eta_end, delta, m0, l0, p0, f,g,h,i,j)#,stop = STOP)   #find t0 using bisection method#
#t0 =Findt02(0, eta_start, eta_end, delta, m0, l0, p0, 0, f,g,h,i,j)#,stop = STOP)    #find t0 using digit building method


# # From here until line 481 the code is used to loop the input for t0 if you want to repeatedly guess at t0 very quickly (useful for guessing values manually)
# # Uncomment this code if you want to enter values manually.
# print(t0)

# Flag = True                                                      #keep looping if true
# while Flag:

#     t0 = float(input('enter t0: '))                              #input your guess at t0

#     if t0  < 0:                                                  #enter a negative number to escape the loop
#         clear_output()                                           #clear the screen so that it isn't filled up
#         Flag = False
#     else:

#         clear_output()

#         print('t0 = ',t0)
#         eta,m,p,l,t,vrho = integrate(eta_start, eta_end, delta, m0, l0, p0, t0,f,g,h,i,j,stop = STOP)  #integrate using given values

# #Plotting
#         fig = figure(figsize=(9.0,11.0))

#         ax = fig.add_subplot(3,2,1)
#         ax.plot(eta,m, '.', label ='m' ,ls = '-', ms = 0)
#         ax.set_xlabel("$\\eta$")
#         ax.set_ylabel("m")
#         ax.ticklabel_format(axis='y', style='sci')

#         ax = fig.add_subplot(3,2,2)
#         ax.plot(eta,p, '.', label ='p',ls = '-', ms = 0)
#         ax.set_xlabel("$\\eta$")
#         ax.set_ylabel("p")

#         ax = fig.add_subplot(3,2,3)
#         ax.plot(eta,l, '.', label ='l',ls = '-', ms = 0)
#         ax.set_xlabel("$\\eta$")
#         ax.set_ylabel("l")
#         ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

#         ax = fig.add_subplot(3,2,4)
#         ax.plot(eta,t, '.', label ='t',ls = '-', ms = 0)
#         ax.set_xlabel("$\\eta$")
#         ax.set_ylabel("t")
#         ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

#         ax = fig.add_subplot(3,2,5)
#         ax.plot(eta,vrho, '.', label ='vrho',ls = '-', ms = 0)
#         ax.set_xlabel("$\\eta$")
#         ax.set_ylabel("$\\varrho$")

#         savefig('Solar Model PLots.pdf', dpi = 1000)
#         show()

#         print('final p: ' , p[-1], 'final t :', t[-1], 'final eta:', eta[-1], 'final m:', m[-1], 'final l:', l[-1], 'final vrho:', vrho[-1])





#Comment out the rest of this code if you are entering t0 values manually

#Use the t0 value you have found and generate the model

print('t0 = ',t0)

#perform the integration using given values

Data = array(integrate(eta_start, eta_end, delta, m0, l0, p0, t0,f,g,h,i,j,stop = STOP))  #

eta,m,p,l,t,vrho = Data

#plotting
fig = figure(figsize=(9.0,11.0))

ax = fig.add_subplot(3,2,1)
ax.plot(eta,m, '.', label ='m' ,ls = '-', ms = 0)
ax.set_xlabel("$\\eta$")
ax.set_ylabel("m")
ax.ticklabel_format(axis='y', style='sci')

ax = fig.add_subplot(3,2,2)
ax.plot(eta,p, '.', label ='p',ls = '-', ms = 0)
ax.set_xlabel("$\\eta$")
ax.set_ylabel("p")

ax = fig.add_subplot(3,2,3)
ax.plot(eta,l, '.', label ='l',ls = '-', ms = 0)
ax.set_xlabel("$\\eta$")
ax.set_ylabel("l")
ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

ax = fig.add_subplot(3,2,4)
ax.plot(eta,t, '.', label ='t',ls = '-', ms = 0)
ax.set_xlabel("$\\eta$")
ax.set_ylabel("t")
ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

ax = fig.add_subplot(3,2,5)
ax.plot(eta,vrho, '.', label ='vrho',ls = '-', ms = 0)
ax.set_xlabel("$\\eta$")
ax.set_ylabel("$\\varrho$")

ax = fig.add_subplot(3,2,6)

text(0.1,0.8,'initial p: ' + str(p[0]))
text(0.1,0.7,'final t: ' + str(t[0]))

text(0.1,0.1,'final p: ' + str(p[-1]))
text(0.1,0.2,'final t :' + str(t[-1]))
text(0.1,0.3,'final eta:' + str(eta[-1]) + ' not necessarily radius')
text(0.1,0.4,'final m:' + str(m[-1]))
text(0.1,0.5,'final l:' + str(l[-1]))
text(0.1,0.6,'final vrho:' + str(vrho[-1]))

file_name = str(round(m[-1],3)) + 'M0 ' +'p0 ' + str(round(p[0],3)) +' t0 ' + str(round(t[0],3)) #generate an appropriate file name


savefig( file_name +' Model Plots.pdf', dpi = 1000)
show()

savetxt( file_name +' M0 Data.txt', Data.T, header = 'eta' + 21 * ' ' + 'm' + 24 * ' ' + 'p' + 24* ' ' +'l' + 24 * ' ' + 't' + 24* ' ' + 'vrho')

print('final p: ' , p[-1], 'final t :', t[-1], 'final eta:', eta[-1], 'final m:', m[-1], 'final l:', l[-1], 'final vrho:', vrho[-1])
