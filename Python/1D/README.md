# Python Finite-Difference-Code

The Python Finite-Difference code is tested with **Python 2.7** and **3.7**. The modules **numpy** and **matplotlib** are required.

There are two versions of each file. The one without `_fast` in the file name is identical to the Matlab version. The `_fast` versions use a vectorization of the loops over the grid points, which will speed up the code significantly.
In Matlab this vectorization is done automatically in the background, therefore, no speed up could be archived by adding this to the Matlab code.  

### Compare accuracy

To compare the accuracy of the different temporal and spatial orders you can use the script `FD_1D_compare.py`, that will plot the seismograms of each FD-script which is available in this repository.

Prior to use this script, you have to run the `*_fast.py` scripts, which will produce and save the seismograms:
```
python FD_1D_DX4_DT3_ABS_fast.py
python FD_1D_DX4_DT4_ABS_fast.py
python FD_1D_DX4_DT4_LW_fast.py
python FD_1D_DX4_DT2_fast.py
python FD_1D_DX8_DT2_fast.py
```
If you get an error using `FD_compare.py` double check that the seismograms are stored inside the folder `Seismograms/`.
