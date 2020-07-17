# Surface Adjusted Distance Measurement

This repository stores programs for the simulation of surface-adjusted distance measurement in digital elevation models. The methodology and results of the case study are documented in the manuscript "_A Systematic Evaluation of Surface-Adjusted Distance Measurements using HPC-enabled Monte Carlo Simulation_" under review by the journal of Geographical Analysis.

## Data collection
This folder contains Python programs for collecting real-time and benchmark travel times from Google Maps API. Specifically:
- **main.py**: The main function. It uses the _map_ function in IPyParallel to calculate distance using different surface-adjustment methods (calls _surfaceAdjusted.py_), in different resolutions of DEM and in 6 different study areas.
- **main.ipynb**: Jupyter Notebook of the main function.

- **surfaceAdjusted.py**: contain the eight surface-adjusted distance measurement methods.

- **neighborhood.py**: get neighborhood pixels around a selected pixel in a raster. Used in the polynomial methods (_local_interp.py_).

- **local_interp.py**: fitting the polynomial functions in the neighborhood pixels around a selected pixel.

- **misc.py**: contains miscellaneous functions.
