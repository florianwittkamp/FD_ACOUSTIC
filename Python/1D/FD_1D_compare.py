## FD_1D_compare.py 1-D acoustic Finite-Difference modelling
# GNU General Public License v3.0
#
# Author: Florian Wittkamp 2016
#
# Finite-Difference acoustic seismic wave simulation
# Discretization of the first-order acoustic wave equation
#
# Compare seismograms, which are calculated with different
# spatial and temporal accuracy.
#
# This script will read in the seismograms to compare them.

## Initialisation
import numpy as np
from matplotlib.pyplot import *

## Execute the following python code
# python FD_1D_DX4_DT3_ABS_fast.py
# python FD_1D_DX4_DT3_ABS_fast.py
# python FD_1D_DX4_DT4_ABS_fast.py
# python FD_1D_DX4_DT4_LW_fast.py
# python FD_1D_DX8_DT2_fast.py

## Read Seismograms
FD_1D_DX4_DT2=np.load("Seismograms/FD_1D_DX4_DT2_fast.npy")
FD_1D_DX4_DT3_ABS=np.load("Seismograms/FD_1D_DX4_DT3_ABS_fast.npy")
FD_1D_DX4_DT4_ABS=np.load("Seismograms/FD_1D_DX4_DT4_ABS_fast.npy")
FD_1D_DX4_DT4_LW=np.load("Seismograms/FD_1D_DX4_DT4_LW_fast.npy")
FD_1D_DX8_DT2=np.load("Seismograms/FD_1D_DX8_DT2_fast.npy")

## Plotting
figure(1)
Line1,=plot(FD_1D_DX4_DT2[0,:])
Line2,=plot(FD_1D_DX4_DT3_ABS[0,:])
Line3,=plot(FD_1D_DX4_DT4_ABS[0,:])
Line4,=plot(FD_1D_DX4_DT4_LW[0,:])
Line5,=plot(FD_1D_DX8_DT2[0,:])
legend((Line1, Line2, Line3, Line4, Line5), ('FD_1D_DX4_DT2', 'FD_1D_DX4_DT3_ABS', 'FD_1D_DX4_DT4_ABS', 'FD_1D_DX4_DT4_LW', 'FD_1D_DX8_DT2'))
title('Seismogram 1')
ylabel('Amplitude')
xlabel('Time in s')

figure(2)
Line1,=plot(FD_1D_DX4_DT2[1,:])
Line2,=plot(FD_1D_DX4_DT3_ABS[1,:])
Line3,=plot(FD_1D_DX4_DT4_ABS[1,:])
Line4,=plot(FD_1D_DX4_DT4_LW[1,:])
Line5,=plot(FD_1D_DX8_DT2[1,:])
legend((Line1, Line2, Line3, Line4, Line5), ('FD_1D_DX4_DT2', 'FD_1D_DX4_DT3_ABS', 'FD_1D_DX4_DT4_ABS', 'FD_1D_DX4_DT4_LW', 'FD_1D_DX8_DT2'))
title('Seismogram 2')
ylabel('Amplitude')
xlabel('Time in s')

figure(3)
Line1,=plot(FD_1D_DX4_DT2[2,:])
Line2,=plot(FD_1D_DX4_DT3_ABS[2,:])
Line3,=plot(FD_1D_DX4_DT4_ABS[2,:])
Line4,=plot(FD_1D_DX4_DT4_LW[2,:])
Line5,=plot(FD_1D_DX8_DT2[2,:])
legend((Line1, Line2, Line3, Line4, Line5), ('FD_1D_DX4_DT2', 'FD_1D_DX4_DT3_ABS', 'FD_1D_DX4_DT4_ABS', 'FD_1D_DX4_DT4_LW', 'FD_1D_DX8_DT2'))
title('Seismogram 3')
ylabel('Amplitude')
xlabel('Time in s')
show()
