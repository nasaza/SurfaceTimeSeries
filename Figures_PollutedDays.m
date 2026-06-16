%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures: Polluted Days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Short Info
% This script creates forecast-surface figures for the two high-pollution
% days discussed in the paper.
%
% Each row of the final figure corresponds to one polluted day. The columns
% show the observed surface and forecasts from the three linear methods that
% performed best across the considered forecasting horizons:
%
%   1. True surface
%   2. PCA VAR
%   3. DS VAR
%   4. MP VAR
%
% Forecasts are first produced for the seasonally adjusted data. The
% seasonal component is then added back so that all surfaces are displayed
% on the original ozone-concentration scale.
%
% Outputs are saved directly to the Outputs folder.

%% Add Libraries

addpath AddFunc
addpath Data

%% User settings
% Practitioners may want to change only this block.

H       = 165;      % Number of forecast origins in the evaluation period
h       = 1;        % Forecast horizon
K_max   = 10;       % Maximum number of scores computed
q_dyn   = 2;        % Lag order in cumulative autocovariance operator

% Forecast-period indices of the two polluted days discussed in the paper.
% These indices refer to positions within the H-period forecast evaluation
% sample, not to positions in the full data set.
PollutedDays = [36 38];

% Seasonal component added back to the seasonally adjusted surfaces.
% In the original empirical illustration both selected observations used
% SeasOz(:,8). This should be checked if the selected days are changed.
seasonalComponentIndex = 8;

% Static PCA specification:
% entries correspond to Ozone, Sun, Precipitation, Temperature, and Wind.
Static_L_set = [5 1 1 1 1];
Static_m_set = [2 1 1 1 1];

% Dynamic-score specification:
Dyn_L_set = [5 1 1 1 1];
Dyn_m_set = [2 1 1 1 1];

Threshold = 0.75;   % Triangulation cleaning threshold

outFolder = fullfile(pwd,'Outputs');

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% Step 1: Read Data

load(fullfile('Data','SeasonAdjData'));  % Seasonally adjusted grid data
load(fullfile('Data','FTSs'));           % Functional data created in Step 1
load(fullfile('Data','SeasComp'));       % Seasonal components

ConstrReg = csvread( ...
    fullfile('Data','GeoConstraints','DE_Constraints.csv'));

[~,T] = size(OzoneS);

%% Step 2: Project Coordinates
% Coordinates are projected to Gauss-Krüger Zone 3, EPSG:31467.

proj = projcrs(31467);

[LonOz,LatOz] = projfwd(proj,LatOz,LonOz);

[x,y]    = projfwd(proj,ConstrReg(:,1),ConstrReg(:,2));
ConstrReg = [y,x];

[XPrecip,YPrecip] = projfwd(proj,YPrecip,XPrecip);
[XWind,YWind]     = projfwd(proj,YWind,XWind);
[XSun,YSun]       = projfwd(proj,YSun,XSun);
[XTemp,YTemp]     = projfwd(proj,YTemp,XTemp);

%% Step 3: Triangulation for Plotting

DTOzone = delaunayTriangulation([LonOz,LatOz]);

CleanDTOz = CleanTriangulation( ...
    DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);

Region     = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord = polyshape(Region);

FEMBasOz = OzoneBasis;
SeasComp = SeasOz(:,seasonalComponentIndex);

%% Step 4: Forecast the Polluted Days

Models = {'True Surface','PCA VAR','DS VAR','MP VAR'};

nDays  = length(PollutedDays);
nModels = length(Models);

PollutedResults = struct();

AllSurfaces = cell(nDays,nModels);

tic

for pp = 1:nDays

    i = PollutedDays(pp);

    fprintf('\nForecasting polluted day %d/%d: index %d.\n', ...
            pp,nDays,i);

    %% Forecast origin

    nobs_i = T-H+i-1;

    % True observation on the original scale
    TrueVal = OzoneS(:,nobs_i+h) + SeasComp;

    %% Static and dynamic ozone scores

    pcastrOz_i = pca3D( ...
        {OzoneCoef(:,1:nobs_i),OzoneBasis},K_max,1);

    dynstrOz_i = DynamScoresSurf( ...
        {OzoneCoef(:,1:nobs_i),OzoneBasis},K_max,1,q_dyn);

    %% Static PCA scores for weather covariates

    pcastrS_i = pca3D( ...
        {SunCoef(:,1:nobs_i),SunBasis},K_max,1);

    pcastrP_i = pca3D( ...
        {PrecipCoef(:,1:nobs_i),PrecipBasis},K_max,1);

    pcastrT_i = pca3D( ...
        {TempCoef(:,1:nobs_i),TempBasis},K_max,1);

    pcastrW_i = pca3D( ...
        {WindCoef(:,1:nobs_i),WindBasis},K_max,1);

    PCA_addreg = {pcastrS_i,pcastrP_i,pcastrT_i,pcastrW_i};

    %% Method 1: PCA VAR

    FuncPredCoeff = FFM_VARX( ...
        pcastrOz_i,PCA_addreg,Static_L_set,Static_m_set,h);

    FuncPred   = fd(FuncPredCoeff,FEMBasOz);
    FPointEval = eval_FEM_fd(LonOz,LatOz,FuncPred);

    MSurf_PCA_VAR = FPointEval + SeasComp;

    %% Method 2: DS VAR

    FuncPredCoeff = FFM_VARX( ...
        dynstrOz_i,PCA_addreg,Dyn_L_set,Dyn_m_set,h);

    FuncPred   = fd(FuncPredCoeff,FEMBasOz);
    FPointEval = eval_FEM_fd(LonOz,LatOz,FuncPred);

    MSurf_DS_VAR = FPointEval + SeasComp;

    %% Method 3: MP VAR
    FPredGrid = DFFM_VARX_GRID( ...
        OzoneS(:,1:nobs_i), ...
        WindS(:,1:nobs_i), ...
        SunS(:,1:nobs_i), ...
        PrecipS(:,1:nobs_i), ...
        TempS(:,1:nobs_i), ...
        Static_L_set,Static_m_set,h);

    MSurf_MP_VAR = FPredGrid + SeasComp;

    %% Store surfaces

    SurfArray = {TrueVal,MSurf_PCA_VAR,MSurf_DS_VAR,MSurf_MP_VAR};

    AllSurfaces(pp,:) = SurfArray;

    PollutedResults(pp).ForecastIndex = i;
    PollutedResults(pp).nobs_i        = nobs_i;
    PollutedResults(pp).Models        = Models;
    PollutedResults(pp).SurfArray     = SurfArray;

end

RunTime = toc;

%% Step 5: Plot Both Polluted Days in One Figure
allVals = [];

for pp = 1:nDays
    for mm = 1:nModels
        allVals = [allVals; AllSurfaces{pp,mm}(:)];
    end
end

bottom = min(allVals);
top    = max(allVals);

fig = figure;
fig.Position = [100,100,1500,700];

TL = tiledlayout(2,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for pp = 1:nDays

    for mm = 1:nModels

        nexttile

        SurfaceFD = Grid2Func( ...
            AllSurfaces{pp,mm},DTOzone,CleanDTOz);

        hold on
            plot(RegionBord,'FaceColor','none');
            plot(SurfaceFD,[],[],[],100);
        hold off

        axis equal tight
        view(2)

        colormap(jet)
        colorbar;
        title(Models{mm});

        xlabel('Easting');
        ylabel('Northing');

        % % Display one color bar per row to reduce visual clutter.
        % if mm == nModels
        %     colorbar;
        % end

    end

end

%title(TL,'Polluted Days: Observed Surfaces and Linear Forecasts');

%% Step 6: Save Figure and Forecast Surfaces

pdfName = fullfile(outFolder,'Figures_PollutedDays_LinearForecasts.pdf');
pngName = fullfile(outFolder,'Figures_PollutedDays_LinearForecasts.png');
figName = fullfile(outFolder,'Figures_PollutedDays_LinearForecasts.fig');
matName = fullfile(outFolder,'Figures_PollutedDays_Results.mat');

exportgraphics(fig,pdfName, ...
    'BackgroundColor','none', ...
    'Resolution',300, ...
    'ContentType','vector');

exportgraphics(fig,pngName, ...
    'BackgroundColor','white', ...
    'Resolution',300);

savefig(fig,figName);

save(matName, ...
    'PollutedResults', ...
    'PollutedDays', ...
    'Models', ...
    'RunTime');

fprintf('\nFigures_PollutedDays completed successfully.\n');
fprintf('Outputs saved to: %s\n',outFolder);

