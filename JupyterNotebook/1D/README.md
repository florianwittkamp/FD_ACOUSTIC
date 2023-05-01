# Jupyter Notebook Finite-Difference Code

The Jupyter Notebook Python Finite-Difference code has been tested with **Python 3.9**. The required modules are **numpy** and **matplotlib**.

There are two versions of each file. Files without `_fast` in their names are identical to the Matlab versions. The `_fast` versions utilize vectorization of the loops over the grid points, resulting in a significant speedup of the code. In Matlab, this vectorization is automatically performed in the background; therefore, no speedup could be achieved by adding this to the Matlab code.

## Compare Accuracy

To compare the accuracy of the different temporal and spatial orders, you can use the `FD_1D_compare.ipynb` notebook, which plots the seismograms of each FD-script available in this repository.

Before using this notebook, you must run the `*_fast.ipynb` notebooks, which produce and save the seismograms:

```
FD_1D_DX4_DT3_ABS_fast.ipynb
FD_1D_DX4_DT4_ABS_fast.ipynb
FD_1D_DX4_DT4_LW_fast.ipynb
FD_1D_DX4_DT2_fast.ipynb
FD_1D_DX8_DT2_fast.ipynb
```
If you encounter an error while using `FD_1D_compare.ipynb`, double-check that the seismograms are stored inside the `Seismograms/` folder.

