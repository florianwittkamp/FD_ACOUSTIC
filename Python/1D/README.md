# Jupyter Notebook Finite-Difference Code

The Jupyter Notebook Python Finite-Difference code has been tested with **Python 3.9**. The required modules are **numpy** and **matplotlib**.

There are two versions of each file. Files without `_fast` in their names are identical to the Matlab versions. The `_fast` versions utilize vectorization of the loops over the grid points, resulting in a significant speedup of the code. In Matlab, this vectorization is automatically performed in the background; therefore, no speedup could be achieved by adding this to the Matlab code.

## Compare Accuracy

To compare the accuracy of the different temporal and spatial orders, you can use the `FD_1D_compare.py` notebook, which plots the seismograms of each FD-script available in this repository.

Before using this notebook, you must run the `*_fast.py` notebooks, which produce and save the seismograms:
```
python FD_1D_DX4_DT3_ABS_fast.py
python FD_1D_DX4_DT4_ABS_fast.py
python FD_1D_DX4_DT4_LW_fast.py
python FD_1D_DX4_DT2_fast.py
python FD_1D_DX8_DT2_fast.py
```
If you get an error using `FD_compare.py` double check that the seismograms are stored inside the folder `Seismograms/`.
