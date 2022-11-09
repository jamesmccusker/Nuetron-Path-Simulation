#!/usr/bin/env python
# coding: utf-8

# # Project 3 - Monte Carlo Techniques
# # Penetration of Neutrons Through Shielding
# ### James McCusker - 10624751

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit as cf


# This project simulates the movement of thermal nuetrons through different materials to test which material would be best to sheild the core of a nuclear reactor. 

# ### Spectral Test
# All linear congruential generators generate points which lie on hyperplanes. This is a problem for a project which depends on modelling the movement of thermal nuetrons as completely random. The spectral test is is a statistical test for the radomness of pseudorandom number generators. The code below generates random points in a 3D plane whcih when orientated, does not have any obvious pattern or planes 

# In[2]:


def random_3D_points(lim1,lim2,n):
    """
    Function generates random points in 3D plane
    Parameters
 
    ----------
 
    lim1 : integer
 
    Lower limit of random number.
 
    lim2 : integer
 
    Upper limit of random number.
    
    n: integer 
    
    number of random numbers
     
     Returns
 
    ---------
    x,y,z: 1D array of floats
 
 
    """
    x=np.random.uniform(lim1,lim2,size = n)
    y=np.random.uniform(lim1,lim2,size = n)
    z=np.random.uniform(lim1,lim2,size = n)
    return(x,y,z)
get_ipython().run_line_magic('matplotlib', 'notebook')
n = 10
(x,y,z) = random_3D_points(0,10,1000)
fig = plt.figure()
plt.rcParams['figure.figsize'] = (9,5)
ax = plt.axes(projection='3d')
ax.scatter(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()


# However, the point generator below is seen to be not as random . If orientated correctly, the points are seen to lie on uniform planes. So for this project the first random number generator will be used.

# In[3]:


def randssp(p,q):
    
    global m, c, x
        
    try: a
    except NameError:
        m = pow(2, 31)
        a = pow(2, 16) + 3
        c = 0
        x = 123456789
    
    try: p
    except NameError:
        p = 1
    try: q
    except NameError:
        q = p
    
    r = np.zeros([p,q])

    for l in range (0, q):
        for k in range (0, p):
            x = np.mod(a*x + c, m)
            r[k, l] = x/m
    
    return r

random_points = randssp(3,1000)
x = random_points[0,:]
y = random_points[1,:]
z = random_points[2,:]
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()


# In[4]:


# defining physical constants
N_A = sc.N_A # Avagadros constant
MR_water = 18.0153 # g/mol
MR_lead = 207.2 # g/mol
MR_graphite = 12.011 # g/mol
absorption_water = 0.6652*10**-24 # cm^2
scattering_water = 103.0*10**-24 # cm^2
density_water = 1 # g/cm^3
absorption_lead = 0.158*10**-24 # cm^2
scattering_lead = 11.221*10**-24 # cm^2
density_lead = 11.35 # g/cm^3
absorption_graphite = 0.0045*10**-24 # cm^2
scattering_graphite = 4.74 *10**-24 # cm^2
density_graphite = 1.67 # g/cm^3
T = 10 #cm


# To find the mean free path of each function as well as the probability, we need to find the macroscopic cross section for absorption and scattering. The macroscopic cross sections can be found from the microscopic cross sections and the number of absorbing molecules. The density of the number of absorbing or scattering molecules can be found through the equation: 
# $$n = \frac{\rho N_{A}}{M}$$
# where $M$ is the molar mass, $N_{A}$ is Avagadro's constant, and $\rho$ is the density of the material. It is assumes that the density of absorption and density of scattering is the same. The macroscopic cross section can be found through the equation: 
# $$\Sigma = n\sigma$$
# where $\sigma$ is the microscopic cross section.

# In[5]:


def macro_cross_section(density, absorption,scattering, MR):
    """
    Function calculates macroscopic absorption scattering cross section
    
    Parameters
 
    ----------
 
    density : float
 
    density of of different material.
 
    absorption : float
 
    cross sectional area 
    
     Returns
 
    ---------
     macro_abs : float
     
     absorption area for given material
     
     macro_scat : float
     
     scattering area for given material
 
    """
    density_abs = density * N_A / MR
    density_scat = density_abs
    macro_abs = density_abs*absorption
    macro_scat = density_scat*scattering
    return macro_abs, macro_scat

macro_abs_water, macro_scat_water = macro_cross_section(density_water, absorption_water,scattering_water, MR_water)
macro_abs_lead, macro_scat_lead = macro_cross_section(density_lead, absorption_lead,scattering_lead, MR_lead)
macro_abs_graphite, macro_scat_graphite = macro_cross_section(density_graphite, absorption_graphite,scattering_graphite, MR_graphite)


# The mean free path can be calculated through the equation:
# $$\lambda = \frac{1}{\Sigma_{T}}$$
# where $\Sigma_{T}$ is the total macroscopic cross section which can be found from the addition of the two macroscopic cross section for scattering and absorption.

# In[6]:


def mean_free_path(macro_abs, macro_scat):
    
    """
    Function calculates the mean free path
    
    Parameters
 
    ----------
 
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    scattering area for given material
 
    
     Returns
 
    ---------
     mean_free_path : float
     
     
     mean free path of nuetron between collisions   
 
    """  
    total_macro = macro_abs + macro_scat
    mean_free_path = 1/ total_macro
    return mean_free_path
water_mfp = mean_free_path(macro_abs_water, macro_scat_water)
lead_mfp = mean_free_path(macro_abs_lead, macro_scat_lead)
graphite_mfp = mean_free_path(macro_abs_graphite, macro_scat_graphite)


# To find the probabilty for a nuetron to be absorbed at each collision, we have to use the macroscopic cross sections for each process as seen in the equation below:
# $$prob(Absorption) = \frac{\Sigma_{a}}{\Sigma_{a}\Sigma_{s}}$$

# In[7]:


def probability_absorption(macro_abs,macro_scat):
    """
    Function calculates probability of absorption
    
    Parameters
 
    ----------
 
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
 
    cross sectional area 
    
     Returns
 
    ---------
    probability_abs: float
     
    probability that nuetron is absorbed during collision
      
    """
    probability_abs = macro_abs/ (macro_scat+macro_abs)
    return probability_abs


# ### Exponential distribution

# The number of molecules absorbed is distributed exponentially as seen below. 
# $$n.absorbed = e^{\frac{-x}{\lambda}}$$
# We can use random numbers to model this by using the pdf of an exponential distribution as seen below:
# $$pdf(x) = e^{-x}$$
# $$cdf(x) = \int_0^xpdf(x) = 1 - e^{-x}$$
# $$cdf^{-1}(x) = -ln(1-x)$$

# In[8]:


n = 10000
def exponential_distribution(n,mean_free_path):
    """
    Creates exponential distribution for given mean_free_path
    
    Parameters
 
    ----------
 
    n : integer
 
    number of random points
 
    mean_free_path: float
    
    mean free path of nuetron
 
    
     Returns
 
    ---------
    X : 1D numpy array
    
    array of numbers correlating to exponential distribution
     
    """
    
    u=np.random.uniform(0,1,n)
    X=-mean_free_path*np.log(u)
    return X

water_exponential = exponential_distribution(n,water_mfp)

fig = plt.figure()
plt.hist(water_exponential, 50)
plt.xlabel("Thickness of water")
plt.ylabel('number of molecules absorbed')
plt.show()


# Now we can plot the natural log of the height of each bin against the number of values in each bin. The negative reciprocal of the gradient of this graph is equal to mean free path.

# In[9]:


frequency, bin_edges = np.histogram(water_exponential,50) # obtaining frequencies and edges of bins from histogram

bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # calculating centre of each bin using the edges
bin_uncertainty = (bin_edges[1:]-bin_edges[0:len(bin_edges)-1]) / 2 # calculating the uncertainty in each bin

bin_centers = bin_centers[np.where((frequency>0))] # removing all zero values so the natural log can be taken
frequency = frequency[np.where((frequency>0))]
bin_uncertainty = bin_uncertainty[np.where((frequency>0))]


# In[10]:


(coef, covr) = np.polyfit(bin_centers, np.log(frequency), 1, cov=True)
Pfit1 = np.polyval(coef, bin_centers)
err_coef1 = np.sqrt(covr[0][0])
sigma_P1 = np.sqrt(np.sum(np.power(frequency-Pfit1, 2))/(len(frequency)-2))

print('Error on the fit = +/-{:04.3g}'.format(err_coef1))
print('The mean free path is = {:,.3g},'.format(-1/coef[0]), 'b = {:,.3f}'.format(coef[1]))

# plotting resulting graph
fig = plt.figure()
plt.errorbar(bin_centers, np.log(frequency), xerr= bin_uncertainty, fmt = '.')
plt.xlabel('Thickness of water')
plt.ylabel('ln(number of molecules absorbed)')
plt.plot(bin_centers, Pfit1,'--k', color='orange')
plt.show()


# ### Vectors

# We can use use random numbers to generate random vectors in 3D space. To do this we need to use the equations for polar coordinates as seen in the function below. However, a uniformly distributed theta angle will cause the density of points across the sphere to be uneven. This is because at the top of a sphere there is a smaller x radius and smaller cross section causing a bigger density of points. 

# In[11]:


def random_3d_vector(n,r):
    """
    Function generates random points from polar coordinates
    
    Parameters
 
    ----------
 
    n : integer
 
    number of points.
 
    r : float
 
    length of points from origin
    
     Returns
 
    ---------
     x,y,z :  1D array of floats
 
    """
    u = np.random.uniform(0,1,n)
    phi = np.random.uniform(0,2*np.pi,n)
    theta = np.arccos(1-2*u) # this ensures points are evenly distributed across the sphere 

    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    
    return(x,y,z)

(x,y,z) = random_3d_vector(1000,1)


get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure()
plt.rcParams['figure.figsize'] = (9,5)
ax = plt.axes(projection='3d')
ax.scatter(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()


# Using the fucntion above, we can generate an isotropic step generator which chooses its length of step based on the exponential distribution function and the direction based on the random 3D vector function.

# In[12]:


def isotropic_step_generator(n,mfp):
    """
    Fucntion generates random points at different distances from the origin
    
    Parameters
 
    ----------
 
    n : integer
 
    number of points.
 
    mfp : float
 
    mean free path of nuetron
    
     Returns
 
    ---------
     x,y,z :  1D array of floats
 
    """
    
    
    theta = np.arccos(1-2*np.random.uniform(0,1,n))
    phi = np.random.uniform(0, 2*np.pi, n)
    length = mfp * np.log(np.random.uniform(0,1,n))
 
    x = length * np.sin(theta) * np.cos(phi)
    y = length * np.sin(theta) * np.sin(phi)
    z = length * np.cos(theta)
    
    return x, y, z

empty_array = np.empty((n,3))

i=0
while i<n:
    empty_array[i] = isotropic_step_generator(1,1) 
    i+=1

x_array = empty_array[:, 0]
y_array = empty_array[:, 1]
z_array = empty_array[:, 2]

get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure()
plt.rcParams['figure.figsize'] = (9,5)
ax = plt.axes(projection='3d')
ax.scatter(x_array, y_array, z_array)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()


# Next, a fucntion to determine which process will happen at each step is used. It takes the length of each step from the exponential distribution and compares it to the thickness of the material. A random number is also generated and compared to the probability of absorption.

# In[13]:


def determine_process(x,T,a,macro_abs,macro_scat):
    """
    Fucntion determins if nuetron is absorbed scattered or transmitted
    
    Parameters
 
    ----------
 
    x : float
 
    Distance moved by nuetron
 
    T : float
 
    Thickness of material
    
    a : float
    
    Random number
    
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    scattering area for given material
     
     Returns
 
    ---------
     1,2,3,0 :  integer
     
     integer corresponds to process nuetron undergoes
 
    """
    if x > T:
        return 1
    if x < 0:
        return 2 
    if a < probability_absorption(macro_abs,macro_scat):
        return 3 
    else:
        return 0


# The function below creates a random walk which checks which process is occouring at each point using the previous function.

# In[14]:


def random_walk(T,mfp,macro_abs,macro_scat):
    """
    Fucntion generates a random walk for nuetron
    
    Parameters
 
    ----------
 
    T : float
 
    Thickness of material
    
    mfp : float
 
    mean free path of nuetron
    
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    scattering area for given material
     
     Returns
 
    ---------
     x_history, y_history, z_history :  1D array of floats
     
     the path taken by nuetron
     
     outcome: string
     
     process nuetron undergoes
 
    """
    x_history = [0] #array of steps taken
    y_history = [0]
    z_history = [0]
    
    x = exponential_distribution(1,mfp) # intial steo is in the x direction
    y = 0
    z = 0
    a = np.random.uniform(0,1,1)
    while determine_process(x,T,a,macro_abs,macro_scat) == 0:

        x_vector, y_vector, z_vector = isotropic_step_generator(1,mfp)
        x = x + x_vector
        y = y + y_vector
        z = z + z_vector
        
        x_history = np.append(x_history,x)
        y_history = np.append(y_history,y)
        z_history = np.append(z_history,z)
    
        
    if determine_process(x,T,a,macro_abs,macro_scat) == 1:
        outcome = 'Transmitted'
    if determine_process(x,T,a,macro_abs,macro_scat) == 2:
        outcome = 'reflected'
    if determine_process(x,T,a,macro_abs,macro_scat) == 3:
        outcome = 'absorbed'
      
    return x_history, y_history, z_history,outcome


# In[15]:


x_water, y_water, z_water, outcome = random_walk(T,water_mfp,macro_abs_water, macro_scat_water)
get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure()
plt.rcParams['figure.figsize'] = (9,5)
ax = plt.axes(projection='3d')
ax.plot(x_water, y_water, z_water)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
print('The mean free path of water is {:,.3f}'.format(water_mfp),'cm')
print('The nuetron was', outcome)


# In[16]:


x_lead, y_lead, z_lead, outcome = random_walk(T,lead_mfp,macro_abs_lead, macro_scat_lead)
get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure()
plt.rcParams['figure.figsize'] = (9,5)
ax = plt.axes(projection='3d')
ax.plot(x_lead, y_lead, z_lead)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
print('The mean free path of lead is {:,.3f}'.format(lead_mfp),'cm')
print('The nuetron was', outcome)


# In[17]:


x_graphite, y_graphite, z_graphite, outcome = random_walk(T,graphite_mfp,macro_abs_graphite, macro_scat_graphite)
get_ipython().run_line_magic('matplotlib', 'notebook')
fig = plt.figure()
plt.rcParams['figure.figsize'] = (9,5)
ax = plt.axes(projection='3d')
ax.plot(x_graphite, y_graphite, z_graphite)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
print('The mean free path of graphite is {:,.3f}'.format(graphite_mfp),'cm')
print('The nuetron was', outcome)


# Next we want to calculate the number of process occouring in each material to see which is the best. 

# In[18]:


def number_of_process(macro_abs,macro_scat,mfp, n, T):
    """
    Fucntion calculates proportion of each process the nuetron undergoes
    
    Parameters
 
    ----------
    
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    scattering area for given material
    
    mfp : float
 
    mean free path of nuetron
 
    T : float
 
    Thickness of material
     
     Returns
 
    ---------
    num_transmission, num_scattering, num_absorbtion :  integer
     
    number of nuetrons which went through each process

 
    """
    num_transmission = 0
    num_scattering = 0
    num_absorption = 0
    x = -mean_free_path(macro_abs, macro_scat)*np.log(np.random.uniform(size=n))
    prob = probability_absorption(macro_abs,macro_scat)
    while len(x) > 0:
        u=np.random.uniform(size=len(x))
        num_transmission = num_transmission + np.count_nonzero(x>T)
        num_scattering = num_scattering + np.count_nonzero(x<0)
        num_absorption = num_absorption + np.count_nonzero(u[np.argwhere((x>0) & (x<T))]<(prob))
        x = np.delete(x, np.argwhere((x<0) | (x>T) | (u<prob)))
        x = x + isotropic_step_generator(len(x),mfp)[0]
    return num_transmission, num_scattering, num_absorption
# calling the number of processes for each material
water_trans, water_scat, water_abs = number_of_process(macro_abs_water,macro_scat_water,water_mfp, n, T)
lead_trans, lead_scat, lead_abs = number_of_process(macro_abs_lead,macro_scat_lead,lead_mfp, n, T)
graphite_trans, graphite_scat, graphite_abs = number_of_process(macro_abs_graphite,macro_scat_graphite,graphite_mfp, n, T)


# In[19]:


# plotting results below on a pie chart
fig,(ax1,ax2,ax3) = plt.subplots(1,3)
water_pie = water_trans, water_scat, water_abs
lead_pie = lead_trans, lead_scat, lead_abs
graphite_pie = graphite_trans, graphite_scat, graphite_abs
process = ['transmitted','scattered','absorbed']

ax1.pie(water_pie, labels = process,  autopct='%1.1f%%', startangle=90,  textprops={'fontsize': 8})
ax1.axis('equal')
ax1.set_title('Water')

ax2.pie(lead_pie, labels = process, autopct='%1.1f%%', startangle=90,  textprops={'fontsize': 8})
ax2.axis('equal') 
ax2.set_title('Lead')

ax3.pie(graphite_pie, labels = process, autopct='%1.1f%%', startangle=90,  textprops={'fontsize': 8})
ax3.axis('equal')
ax3.set_title('Graphite')
plt.show()


# We now want to see how varying the thickness on each material effects the number of processes taking place. To do this we create a numpy geomspace so we can have more points towards the begging of the range where the exponential distribution is more accurate. Below we collect an array of values for scattered absorbed and transmitted for a range of thickness values.

# In[20]:


T = np.geomspace(0.001,10,num=50) 
def thickness_change(macro_abs, macro_scat, mfp):
    """
    Fucntion calculates number of nuetrons absorbed scattered and transmitted for a range of thickness
    
    Parameters
 
    ----------
    
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    scattering area for given material
    
    mfp : float
 
    mean free path of nuetron
 
     Returns
 
    ---------
    num_transmission, num_scattering, num_absorbtion :  integer
     
    number of nuetrons which went through each process

    """
  
    T = np.geomspace(0.001,10,num=50)
    num_transmission = np.zeros(len(T))
    num_scattering = np.zeros(len(T))
    num_absorption = np.zeros(len(T))
    for i in range(len(T)):
        num_transmission[i] =  number_of_process(macro_abs,macro_scat,mfp, n, T[i])[0]
        num_scattering[i] = number_of_process(macro_abs,macro_scat,mfp, n, T[i])[1]
        num_absorption[i] = number_of_process(macro_abs,macro_scat,mfp, n, T[i])[2]
        
    return num_transmission, num_scattering, num_absorption
    
num_water_trans, num_water_abs, num_water_scat = thickness_change(macro_abs_water, macro_scat_water, water_mfp)
num_lead_trans, num_lead_abs, num_lead_scat = thickness_change(macro_abs_lead, macro_scat_lead, lead_mfp)
num_graphite_trans, num_graphite_abs, num_graphite_scat = thickness_change(macro_abs_graphite, macro_scat_graphite, graphite_mfp)


# Below are the plots of each of the three processes in water, lead and graphite respectively. 

# In[21]:


def plot_function(num_trans,num_abs,num_scat):
    """
    plots graph for transmission, scattering and absorption for a range of thickness
    
    Parameters
 
    ----------
    
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    scattering area for given material
    
    mfp : float
 
    mean free path of nuetron
 
    T : float
 
    Thickness of material
     
     Returns
 
    ---------
    num_transmission, num_scattering, num_absorbtion :  integer
     
    number of nuetrons which went through each process

 
    """
    plt.rcParams['figure.figsize']=(8,2)
    fig,(ax1,ax2,ax3) = plt.subplots(1,3)
    
    ax1.scatter(T,num_trans)
    ax1.set_title('Transmission')
    ax1.set(xlabel='Thickness')

    ax2.scatter(T,num_scat) 
    ax2.set_title('Scatterring')
    ax2.set(xlabel='Thickness')

    ax3.scatter(T,num_abs)
    ax3.set_title('Absorption')
    ax3.set(xlabel='Thickness')
    plt.show()
    
    return None
    
water_plots = plot_function(num_water_trans, num_water_abs, num_water_scat)
lead_plots = plot_function(num_lead_trans, num_lead_abs, num_lead_scat)
graphite_plots = plot_function(num_graphite_trans, num_graphite_abs, num_graphite_scat)


# We create an exponential fit for the transmission curve to find the attenuation length. We apply the fit to the transmission curves for each material and plot the natural log of the number of nuetrons transmitted. The attenuation length is the negative reciprocal of the gradient.

# In[22]:


def fit(T,num_transmitted,title):
    """
    Fucntion does a fit for the thickness function with varying thickness
    
    Parameters
 
    ----------
    
     T : float
 
    Thickness of material
    
    num_transmitted : integer
    
    number of transmitted nuetrons
    
    title : string
    
    title of graph
     
     Returns
 
    ---------
    None

 
    """
    (coef, covr) = np.polyfit(T, np.log(num_transmitted), 1, cov=True)
    Pfit2 = np.polyval(coef, T)
    err_coef1 = np.sqrt(covr[0][0])
    err_coef2 = np.sqrt(covr[1][1])
    sigma_P1 = np.sqrt(np.sum(np.power(num_transmitted-Pfit2, 2))/(len(num_transmitted)-2))

    print('Error on the fit = {:04.3g} +/- {:04.3g}'.format(coef[0], err_coef1))
    print('The attenuation length is {:,.3g},'.format(-1/coef[0]), 'b = {:,.3f}'.format(coef[1]))

    fig = plt.figure()
    plt.rcParams['figure.figsize']=(8,6)
    plt.errorbar(T, np.log(num_transmitted), fmt = '.')
    plt.xlabel('Thickness / cm')
    plt.ylabel('ln(number of molecules transmitted)')
    plt.title(title)
    plt.plot(T, Pfit2,'--k', color='orange')
    plt.show()
    
    return None

water_fit = fit(T,num_water_trans,'Transmission for water')
lead_fit = fit(T,num_lead_trans, 'Transmission for lead')
graphite_fit = fit(T,num_graphite_trans,'Transmission for graphite')


# The function below runs the random walk function multiple times to create an array for the number of transmitted, absorbed and scattered nuetrons at a given thickness. This data is then used in error analysis.

# In[25]:


T = 10

def percent_process(num_trans, num_scat, num_abs, macro_abs,macro_scat,mfp, iterations):
    """
    Function calculates percentage of each process.
    
    Parameters
 
    ----------
 
    num_trans : integer
 
    number of transmitted nuetrons
 
    num_scat : integer
 
    number of scattered nuetrons
    
    num_abs : integer
    
    number of absorbed nuetrons
    
    macro_abs : float
 
    absorption area for given material.
 
    macro_scat : float
    
    mfp : float
    
    mean free path
    
    iterations : integer
    
    number of times function repeated
    
     Returns
 
    ---------
    transmitted, scattered, absorbed : 1D array of float
    
    fraction of each process
 
    """
    transmitted = np.empty((0))
    scattered = np.empty((0))
    absorbed = np.empty((0))
    for iteration in np.arange(iterations):  
        num_trans, num_scat, num_abs =  number_of_process(macro_abs,macro_scat,mfp, n, T)[0:3]
        total_processes = num_trans+num_scat+num_abs
        transmitted_percent = num_trans / total_processes
        scattered_percent = num_scat / total_processes
        absorbed_percent = num_abs / total_processes 
        transmitted = np.append(transmitted, transmitted_percent)
        scattered = np.append(scattered, scattered_percent) 
        absorbed = np.append(absorbed, absorbed_percent)
    return transmitted, scattered, absorbed
transmitted_water, scattered_water, absorbed_water = percent_process(num_water_trans, num_water_scat, num_water_abs, macro_abs_water,macro_scat_water,water_mfp,50)
transmitted_lead, scattered_lead, absorbed_lead = percent_process(num_lead_trans, num_lead_scat, num_lead_abs, macro_abs_lead,macro_scat_lead,lead_mfp,50)
transmitted_graphite, scattered_graphite, absorbed_graphite = percent_process(num_graphite_trans, num_graphite_scat, num_graphite_abs, macro_abs_graphite,macro_scat_graphite,graphite_mfp,50)


# The error is then calculated from the equation:
# $$error = \frac{\sigma}{mean}$$
# The results are plotted below. 

# In[26]:


def error(array):
    """
    Function calculates mean and error for each process.
    
    Parameters
 
    ----------
    
    array: 1D array of floats
    
    fraction of each process
    
     Returns
 
    ---------
    error : float
    
    error of array
    
    mean_array: float
    
    mean of array
 
    """
    mean_array = np.mean(array)
    std_array = np.std(array)
    error = std_array / mean_array
    return error, mean_array

water_trans_error,mean_water_trans = error(transmitted_water)
water_scat_error,mean_water_abs = error(scattered_water)
water_abs_error,mean_water_scat = error(absorbed_water)

lead_trans_error,mean_lead_trans = error(transmitted_lead)
lead_scat_error,mean_lead_scat = error(scattered_lead)
lead_abs_error,mean_lead_absorbed = error(absorbed_lead)

graphite_trans_error,mean_graphite_trans = error(transmitted_graphite)
graphite_scat_error, mean_graphite_scat = error(scattered_graphite)
graphite_abs_error, mean_graphite_abs = error(absorbed_graphite)

def print_errors(material,mean_trans,mean_scat,mean_abs,error_trans,error_scat,error_abs):
    print('For a thickness of 10cm in',material,'the fraction of each process is:')
    print('Transmission:{:,.3f}+/-{:,.3f}'.format(mean_trans,error_trans))
    print('Scattering:{:,.3f}+/-{:,.3f}'.format(mean_scat,error_scat))
    print('Absorption:{:,.3f}+/-{:,.3f}'.format(mean_abs,error_abs))
    
    return None

water_errors = print_errors('Water',mean_water_trans,mean_water_abs,mean_water_scat,water_trans_error,water_scat_error,water_abs_error)
print('')
lead_errors = print_errors('Lead',mean_lead_trans,mean_lead_scat,mean_lead_absorbed,lead_trans_error,lead_scat_error,lead_abs_error)
print('')
graphite_errors = print_errors('Graphite',mean_graphite_trans,mean_graphite_scat,mean_graphite_abs,graphite_trans_error,graphite_scat_error,graphite_abs_error)


# ## Conclusion

# In conclusion, the results show that the best material to use in sheilding a nuclear reactor is water. It has the lowest fraction of transmitted nuetrons. However the error on this value is large leaving some uncertainty in this result. We can see that Graphite would not be a good material to use as it has the highest tranmission rate and a very low absorption.
