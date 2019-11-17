# # FD_2D_DX4_DT2 2-D acoustic Finite-Difference modelling
#
# GNU General Public License v3.0
#
# Author: Florian Wittkamp
#
# Finite-Difference acoustic seismic wave simulation
#
# Discretization of the first-order acoustic wave equation
#
# Temporal second-order accuracy $O(\Delta T^2)$
# Spatial fourth-order accuracy  $O(\Delta X^4)$

# ##  Initialisation
import numpy as np
import matplotlib.pyplot as plt

## Input Parameter

# Discretization
c1=20   # Number of grid points per dominant wavelength
c2=0.5  # CFL-Number
nx=200 # Number of grid points in X
ny=200 # Number of grid points in Y
T=1     # Total propagation time

# Source Signal
f0= 5      # Center frequency Ricker-wavelet
q0= 1       # Maximum amplitude Ricker-Wavelet
xscr = 100  # Source position (in grid points) in X
yscr = 100  # Source position (in grid points) in Y

# Receiver
xrec1=100; yrec1=80;  # Position Reciever 1 (in grid points)
xrec2=100; yrec2=100;  # Position Reciever 2 (in grid points)
xrec3=100; yrec3=120;# Position Reciever 3 (in grid points)

# Velocity and density
modell_v = 3000*np.ones((ny,nx))
rho=2.2*np.ones((ny,nx))

## Preparation

# Init wavefields
vx=np.zeros((ny,nx))
vy=np.zeros((ny,nx))
p=np.zeros((ny,nx))

# Calculate first Lame-Paramter
l=rho * modell_v * modell_v

cmin=min(modell_v.flatten())  # Lowest P-wave velocity
cmax=max(modell_v.flatten())  # Highest P-wave velocity
fmax=2*f0                     # Maximum frequency
dx=cmin/(fmax*c1)             # Spatial discretization (in m)
dy=dx                         # Spatial discretization (in m)
dt=dx/(cmax)*c2               # Temporal discretization (in s)
lampda_min=cmin/fmax          # Smallest wavelength

# Output model parameter:
print("Model size: x:",dx*nx,"in m, y:",dy*ny,"in m")
print("Temporal discretization: ",dt," s")
print("Spatial discretization: ",dx," m")
print("Number of gridpoints per minimum wavelength: ",lampda_min/dx)

# ## Create space and time vector

x=np.arange(0,dx*nx,dx) # Space vector in X
y=np.arange(0,dy*ny,dy) # Space vector in Y
t=np.arange(0,T,dt)     # Time vector
nt=np.size(t)           # Number of time steps

# Plotting model
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.subplots_adjust(wspace=0.4,right=1)
ax1.plot(x,modell_v)
ax1.set_ylabel('VP in m/s')
ax1.set_xlabel('Depth in m')
ax1.set_title('P-wave velocity')
plt.draw()
plt.pause(0.001)

ax2.plot(x,rho)
ax2.set_ylabel('Density in g/cm^3')
ax2.set_xlabel('Depth in m')
ax2.set_title('Density');

plt.draw()
plt.pause(0.001)

# ## Source signal - Ricker-wavelet

tau=np.pi*f0*(t-1.5/f0)
q=q0*(1.0-2.0*tau**2.0)*np.exp(-tau**2)

# Plotting source signal
plt.figure(3)
plt.plot(t,q)
plt.title('Source signal Ricker-Wavelet')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')
plt.draw()
plt.pause(0.001)

# ## Time stepping

# Init Seismograms
Seismogramm=np.zeros((3,nt)); # Three seismograms

# Calculation of some coefficients
i_dx=1.0/(dx)
i_dy=1.0/(dy)
c1=9.0/(8.0*dx)
c2=1.0/(24.0*dx)
c3=9.0/(8.0*dy)
c4=1.0/(24.0*dy)
c5=1.0/np.power(dx,3)
c6=1.0/np.power(dy,3)
c7=1.0/np.power(dx,2)
c8=1.0/np.power(dy,2)
c9=np.power(dt,3)/24.0

## Time stepping
print("Starting time stepping...")
for n in range(2,nt):

        # Update velocity
        for kx in range(5,nx-4):
            for ky in range(5,ny-4):

                p_x=c1*(p[ky,kx+1]-p[ky,kx])-c2*(p[ky,kx+2]-p[ky,kx-1])
                p_y=c3*(p[ky+1,kx]-p[ky,kx])-c4*(p[ky+2,kx]-p[ky-1,kx])

                vx[ky,kx]=vx[ky,kx]-dt/rho[ky,kx]*p_x
                vy[ky,kx]=vy[ky,kx]-dt/rho[ky,kx]*p_y

        # Inject source wavelet
        p[yscr,xscr]=p[yscr,xscr]+q[n]

        # Update pressure
        for kx in range(5,nx-4):
            for ky in range(5,ny-4):

                vx_x=c1*(vx[ky,kx]-vx[ky,kx-1])-c2*(vx[ky,kx+1]-vx[ky,kx-2])
                vy_y=c3*(vy[ky,kx]-vy[ky-1,kx])-c4*(vy[ky+1,kx]-vy[ky-2,kx])

                p[ky,kx]=p[ky,kx]-l[ky,kx]*dt*(vx_x+vy_y)

        # Save seismograms
        Seismogramm[0,n]=p[yrec1,xrec1]
        Seismogramm[1,n]=p[yrec2,xrec2]
        Seismogramm[2,n]=p[yrec3,xrec3]

print("Finished time stepping!")


# ## Save seismograms

## Save seismograms
np.save("Seismograms/FD_2D_DX4_DT2",Seismogramm)

# ## Plotting

## Image plot
fig, ax = plt.subplots(1,1)
img = ax.imshow(p);
ax.set_title('P-Wavefield')
ax.set_xticks(range(0,nx+1,np.int(nx/5)))
ax.set_yticks(range(0,ny+1,np.int(ny/5)))
ax.set_xlabel('Grid-points in X')
ax.set_ylabel('Grid-points in Y')
fig.colorbar(img)
plt.draw()
plt.pause(0.001)

## Plot seismograms
plt.figure()
plt.plot(t,Seismogramm[0,:])
plt.title('Seismogram 1')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')

plt.figure()
plt.plot(t,Seismogramm[1,:])
plt.title('Seismogram 2')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')

plt.figure()
plt.plot(t,Seismogramm[2,:])
plt.title('Seismogram 3')
plt.ylabel('Amplitude')
plt.xlabel('Time in s')
plt.draw()

plt.show()

print(" ")
