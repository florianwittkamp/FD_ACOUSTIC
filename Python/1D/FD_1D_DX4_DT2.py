## FD_1D_DX4_DT2.py 1-D acoustic Finite-Difference modelling
# GNU General Public License v3.0
#
# Author: Florian Wittkamp 2016
#
# Finite-Difference acoustic seismic wave simulation
# Discretization of the first-order acoustic wave equation
#
# Temporal second-order accuracy O(DT^2)
# Spatial fourth-order accuracy  O(DX^4)

##  Initialisation
print(" ")
print("Starting FD_1D_DX4_DT2")

import numpy as np
import time as tm
import matplotlib.pyplot as plt

## Input Parameter

# Discretization
c1=20   # Number of grid points per dominant wavelength
c2=0.5  # CFL-Number
nx=2000 # Number of grid points
T=10     # Total propagation time

# Source Signal
f0=10      # Center frequency Ricker-wavelet
q0=1       # Maximum amplitude Ricker-Wavelet
xscr=100  # Source position (in grid points)

# Receiver
xrec1=400  # Position Reciever 1 (in grid points)
xrec2=800  # Position Reciever 2 (in grid points)
xrec3=1800 # Position Reciever 3 (in grid points)

# Velocity and density
modell_v = np.hstack((1000*np.ones((np.int(nx/2))),1500*np.ones((np.int(nx/2)))))
rho=np.hstack((1*np.ones((np.int(nx/2))),1.5*np.ones((np.int(nx/2)))))

## Preparation

# Init wavefields
vx=np.zeros((nx),float)
p=np.zeros((nx),float)

# Calculate first Lame-Paramter
l=rho * modell_v * modell_v

cmin=min(modell_v.flatten())  # Lowest P-wave velocity
cmax=max(modell_v.flatten())  # Highest P-wave velocity
fmax=2*f0                     # Maximum frequency
dx=cmin/(fmax*c1)             # Spatial discretization (in m)
dt=dx/(cmax)*c2               # Temporal discretization (in s)
lampda_min=cmin/fmax          # Smallest wavelength

# Output model parameter:
print("Model size: x:",dx*nx,"in m")
print("Temporal discretization: ",dt," s")
print("Spatial discretization: ",dx," m")
print("Number of gridpoints per minimum wavelength: ",lampda_min/dx)

# Create space and time vector
x=np.arange(0,dx*nx,dx) # Space vector
t=np.arange(0,T,dt)     # Time vector
nt=np.size(t)           # Number of time steps

# Plotting model
plt.figure(1)
plt.plot(x,modell_v)
plt.ylabel('VP in m/s')
plt.xlabel('Depth in m')
plt.figure(2)
plt.plot(x,rho)
plt.ylabel('Density in g/cm^3')
plt.xlabel('Depth in m')
plt.draw()
plt.pause(0.001)

# Source signal - Ricker-wavelet
tau=np.pi*f0*(t-1.5/f0)
q=q0*(1-2*tau**2)*np.exp(-tau**2)

# Plotting source signal
plt.figure(3)
plt.plot(t,q)
plt.title('Source signal Ricker-Wavelet')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')
plt.draw()
plt.pause(0.001)

# Init Seismograms
Seismogramm=np.zeros((3,nt)); # Three seismograms

# Calculation of some coefficients
i_dx=1.0/(dx)

## Time stepping
print("Starting time stepping...")
for n in range(2,nt):

        # Inject source wavelet
        p[xscr]=p[xscr]+q[n]

        # Update velocity
        for kx in range(5,nx-4):
            # Calculating spatial derivative
            p_x=i_dx*9.0/8.0*(p[kx+1]-p[kx])-i_dx*1.0/24.0*(p[kx+2]-p[kx-1])

            # Update velocity
            vx[kx]=vx[kx]-dt/rho[kx]*p_x

        # Update pressure
        for kx in range(5,nx-4):

            # Calculating spatial derivative
            vx_x= i_dx*9.0/8.0*(vx[kx]-vx[kx-1])-i_dx*1.0/24.0*(vx[kx+1]-vx[kx-2])

            # Update pressure
            p[kx]=p[kx]-l[kx]*dt*(vx_x);

        # Save seismograms
        Seismogramm[0,n]=p[xrec1]
        Seismogramm[1,n]=p[xrec2]
        Seismogramm[2,n]=p[xrec3]

print("Finished time stepping!")

## Plot seismograms
plt.figure(4)
plt.plot(t,Seismogramm[0,:])
plt.title('Seismogram 1')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')

plt.figure(5)
plt.plot(t,Seismogramm[1,:])
plt.title('Seismogram 2')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')

plt.figure(6)
plt.plot(t,Seismogramm[2,:])
plt.title('Seismogram 3')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')
plt.draw()

plt.show()

print(" ")
