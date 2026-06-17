# SurfaceTimeSeries

MATLAB replication package for the empirical ozone-forecasting illustration in the manuscript **“Practical Forecasting of Environmental Maps: A Functional Data Approach.”** by Alexander Gleim and Nazarii Salish

The package reconstructs environmental surfaces observed at irregular spatial locations, extracts static or dynamically relevant low-dimensional components, and compares several forecasting approaches.

## Main workflow

The package implements the following steps:

1. reconstruct surface time series over the irregular geographic domain using finite-element basis functions;
2. estimate static FPCA components and dynamic components based on cumulative autocovariances;
3. forecast ozone surfaces using functional, score-based, multivariate-perspective, nearest-neighbour, and random-forest methods;
4. evaluate forecasts against reconstructed ozone surfaces;
5. produce the main figures and additional diagnostics used in the empirical analysis.

## Repository structure

```text
SurfaceTimeSeries/
├── Step1_CreateSurfaceTimeSeries.m
├── Step2_RunForecastingComparison.m
├── Figures_PollutedDays.m
├── README.md
├── LICENSE
├── .gitignore
├── AddFunc/
│   └── supporting MATLAB functions
├── Data/
│   ├── SeasonAdjData.mat
│   ├── SeasComp.mat
│   ├── cleaned CSV data
│   └── GeoConstraints/
│       └── DE_Constraints.csv
└── Diagnostics/
    ├── AdditionalDiagnostics.m
    ├── GridSearch_OttoSalish_Ozone.m
    └── SelectFPCA_AueFFPE_Ozone.m
```

The generated folders and files

```text
Outputs/
Data/FTSs.mat
```

are intentionally not stored in the repository. They are created locally when the scripts are run.

## Requirements

A recent MATLAB release is recommended. The code uses functionality such as `tiledlayout`, `exportgraphics`, and `clim`.

Required MATLAB products include:

- MATLAB;
- Mapping Toolbox;
- Econometrics Toolbox;
- Statistics and Machine Learning Toolbox.

The package also requires the **Functional Data Analysis MATLAB package by Ramsay and co-authors**. This package provides the `fd` class and related functionality used throughout the code. In particular, functions such as `eval_FEM_fd.m` rely on the Ramsay FDA package.

Before running the replication scripts, install the Ramsay FDA MATLAB package and ensure that its folders are available on the MATLAB search path.

## Data

The `Data` folder contains the cleaned inputs needed for the empirical illustration:

- seasonally adjusted ozone and weather observations;
- seasonal components used to return selected forecasts to the original ozone scale;
- station coordinates and cleaned CSV inputs;
- the geographic boundary of Germany.

The large file `Data/FTSs.mat` is not distributed because it can be reproduced from the included data by running Step 1.

### Step 1: Reconstruct the surface time series

Run:

```matlab
Step1_CreateSurfaceTimeSeries
```

This script:

- reads the cleaned data;
- projects the coordinates to Gauss–Krüger Zone 3 (EPSG:31467);
- constructs and cleans the triangulations;
- reconstructs the functional surface observations;
- saves the required coefficient matrices and basis objects in:

```text
Data/FTSs.mat
```

It also creates the triangulation figure in the `Outputs` folder.

Step 1 must be run before the remaining scripts!!!

### Step 2: Run the forecasting comparison

Run:

```matlab
Step2_RunForecastingComparison
```

The script compares:

1. mean forecast;
2. naive forecast;
3. functional autoregression;
4. static FPCA scores with VARX;
5. static FPCA scores with KNN;
6. dynamic scores with VARX;
7. dynamic scores with KNN;
8. multivariate-perspective factors with VARX;
9. multivariate-perspective factors with KNN;
10. dynamic scores with random forest.

The default specification reproduces the main one-step-ahead forecasting exercise. The current implementation of the multivariate-perspective KNN procedure is one-step ahead and therefore requires `h = 1`.

Numerical results, tables, figures, and spatial MSE surfaces are written to `Outputs/`.

### Polluted-day forecast surfaces

Run:

```matlab
Figures_PollutedDays
```

This script creates a single two-row figure for the two high-pollution days discussed in the paper. It compares:

- the reconstructed reference surface;
- PCA VAR;
- dynamic-score VAR;
- multivariate-perspective VAR.

The seasonal component is added back so that all panels are shown on the original ozone-concentration scale.

## Additional diagnostics

The scripts in `Diagnostics/` can be run after Step 1.

### Static and dynamic score diagnostics

```matlab
run(fullfile('Diagnostics','AdditionalDiagnostics.m'))
```

This produces:

- static and dynamic scree plots;
- ACF and PACF plots for the first dynamic score series;
- the first three static and dynamic loading surfaces with sign alignment.

### Otto–Salish dynamic-score selection

```matlab
run(fullfile('Diagnostics','GridSearch_OttoSalish_Ozone.m'))
```

This jointly examines the number of dynamic scores and the VAR lag order using BIC- and HQC-type criteria.

### Aue et al. fFPE selection

```matlab
run(fullfile('Diagnostics','SelectFPCA_AueFFPE_Ozone.m'))
```

This jointly examines the number of static FPCA scores and the VAR lag order using the functional final prediction error criterion.

All diagnostic results are saved in `Outputs/`.

## User settings

The main scripts collect the principal tuning choices near the beginning of each file. These include:

- forecast horizon and evaluation length;
- maximum number of components;
- static and dynamic score dimensions;
- VAR lag orders;
- the cumulative-autocovariance lag;
- KNN search range;
- triangulation threshold.

Users changing the main forecasting specification should keep the settings in `Figures_PollutedDays.m` consistent with those in `Step2_RunForecastingComparison.m`.

## Citation

When using this software in academic or applied research, please cite the associated manuscript:

> *Practical Forecasting of Environmental Maps: A Functional Data Approach.* 
Working-paper version, 2026. https://arxiv.org/abs/2202.03332

A complete bibliographic citation can be added here once the article is published.

## License

This software is released under the MIT License. It may be used, copied, modified, and redistributed, provided that the original copyright and license notice are retained.

Copyright © 2026 The authors.

Academic and research users are kindly requested to cite the associated paper listed above.
