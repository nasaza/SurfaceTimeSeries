# SurfaceTimeSeries 
Package to work with surface time series

1) Overview

This MATLAB package implements functional time series reconstruction and forecasting as described in 
"On Forecasting Environmental Data: An example to ground-level ozone concentrations" (by Alexander Gleim and Nazarii Salish) 

2) Prerequisites
- Matlab R2024
- Before running this package, ensure you have installed the MATLAB FDA package by Ramsay and Silverman. You can download it from:

https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/

3) Usage

Follow these steps to use the package:

Run Step1_Reconstruction.m: Reconstruct functional time series from gridded data.

Run Step2_Forecating: Use the reconstructed data to build forecasts.

Run the scripts in the indicated order to ensure proper execution.

4) Output Structure

Outputs/ – Final forecast results and supporting figures are saved in this folder (created after running Step1_Reconstruction.m).

Data/ – Intermediate results, such as reconstructed functional time series, are stored here.

5) Licence/Reference

If you use this package in your research, please cite our paper:

"On Forecasting Environmental Data: An example to ground-level ozone concentrations" (by Alexander Gleim and Nazarii Salish)

This package is provided for research purposes. Please acknowledge its use in any resulting publications.
