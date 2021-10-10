# Visualizing Snapshots

The script `plot_snapshots.py` located in the [Postprocessing](https://github.com/PrincetonUniversity/EDIPIC-2D/tree/main/Postprocessing) directory enables rapid visualization of EDIPIC-2D snapshot data output by the code.

## System Requirements

To run the script your system must have a recent installation of **Python 3** with the following packages:
- numpy
- scipy
- matplotlib


## How to Use the Script

Copy `plot_snapshots.py` into the directory where you are running EDIPIC-2D (i.e. where the `edipic2d` executable is located)

Execute the code using `python3 plot_snapshots.py`.

The script will prompt you with four questions. Example responses are given below:

```
Enter the snapshot parameter name (e.g. Ne/EX,TXi):
Ne          <--- Here the user has chose to plot the electron density
Enter the number of steps between plots (default is 1):
10          <--- Here the user has chosen to plot every 10th snapshot
Enter zmin for the surface plot (default is the minimum of each data set):
0.0         <--- Here the user has chosen to set the lower bound of the contour plots to 0.0
Enter zmax for the surface plot (default is the maxmimum of each data set):
1.0e17      <--- Here the user has chosen to set the upper bound of the contour plots to 1.0e17
```
The script will then count off the data files which it is reading. All snapshots are saved into a new directory with the name of the variable appended by `_plots`.


Some notes on the inputs:
1. The plotting variable names correspond exactly to those in the snapshot binary file names. The list of possible names correspond to those listed in the `init_snapshots.dat` file. If the file name is not detected then the script will exit.
2. If the user does not enter lower and upper bounds for the contour plots then each plot adjusts these based on the individual data set.
3. If the directory for plotting is already detected then the script will delete any content within that directory before plotting
