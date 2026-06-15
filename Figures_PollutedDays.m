%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures: Polluted Days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Short Info
% This script creates forecast-surface figures for the two high-pollution
% days in the forecasting period. The script uses the same forecasting
% specifications as the main forecasting comparison and adds back the
% seasonal component in order to plot the surfaces on the original scale.
%
% Output files are saved directly to the Outputs folder.

%% Add Libraries

addpath AddFunc
addpath Data

%% User settings
% Practitioners may want to change only this block.

H        = 165;      % Number of forecast origins in the evaluation period
h        = 1;        % Forecast horizon
K_max    = 10;       % Maximum number of scores computed
KNN_max  = 50;       % Maximum number of neighbours searched in KNN methods
q_dyn    = 2;        % Lag order in cumulative autocovariance operator

% These are the two forecast-period indices corresponding to the polluted
% days discussed in the paper. If the data set is changed, these indices
% should be updated accordingly.
PollutedDays = [36 38];

% Seasonal component added back to seasonally adjusted forecasts.
% In the original empirical illustration this was SeasOz(:,8).
% If the polluted days change, verify that this seasonal component is still
% the appropriate one.
seasonalComponentIndex = 8;

% Number of components and lags used in the forecasting methods.
% The first entry corresponds to ozone. The following entries correspond to
% Sun, Precipitation, Temperature, and Wind, respectively.
Static_L_set = [3 1 1 1 1];
Static_m_set = [2 1 1 1 1];

Dyn_L_set    = [4 1 1 1 1];
Dyn_m_set    = [2 1 1 1 1];

Threshold = 0.75;   % Triangulation cleaning threshold

outFolder = fullfile(pwd,'Outputs');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% Step 1: Read Data

load(fullfile('Data','SeasonAdjData'));  % Seasonally adjusted grid data
load(fullfile('Data','FTSs'));           % Functional/surface data from Step 1
load(fullfile('Data','SeasComp'));       % Seasonal components
ConstrReg = csvread(fullfile('Data','GeoConstraints','DE_Constraints.csv'));
[~,T] = size(OzoneS);

%% Step 2: Project Coordinates
% Coordinates are projected to Gauss-Krüger Zone 3, EPSG:31467.

proj = projcrs(31467);

[LonOz,  LatOz] = projfwd(proj, LatOz, LonOz);
[x, y]          = projfwd(proj, ConstrReg(:,1), ConstrReg(:,2));
ConstrReg       = [y,x];

[XPrecip,YPrecip] = projfwd(proj, YPrecip, XPrecip);
[XWind,  YWind]   = projfwd(proj, YWind, XWind);
[XSun,   YSun]    = projfwd(proj, YSun, XSun);
[XTemp,  YTemp]   = projfwd(proj, YTemp, XTemp);

%% Step 3: Triangulation for Plotting

DTOzone   = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);

Region     = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord = polyshape(Region);

FEMBasOz = OzoneBasis;
SeasComp = SeasOz(:,seasonalComponentIndex);

%% Step 4: Forecast and Plot the Polluted Days

Models = {'True Surface','MF','NF','FAR','PCA VAR','PCA KNN','DS VAR','DS KNN', ...
          'MP VAR','MP KNN','DS RF'};
PollutedResults = struct();

tic

for pp = 1:length(PollutedDays)

    i = PollutedDays(pp);

    fprintf('\nCreating polluted-day figure %d/%d: forecast-period index %d.\n', ...
            pp, length(PollutedDays), i);

    %% Forecast origin

    nobs_i  = T-H+i-1;
    TrueVal = OzoneS(:,nobs_i+h) + SeasComp;

    %% Static and dynamic scores

    dynstrOz_i = DynamScoresSurf({OzoneCoef(:,1:nobs_i),OzoneBasis}, ...
                                  K_max, 1, q_dyn);

    pcastrOz_i = pca3D({OzoneCoef(:,1:nobs_i),OzoneBasis}, K_max, 1);

    pcastrS_i  = pca3D({SunCoef(:,1:nobs_i),SunBasis}, K_max, 1);
    pcastrP_i  = pca3D({PrecipCoef(:,1:nobs_i),PrecipBasis}, K_max, 1);
    pcastrT_i  = pca3D({TempCoef(:,1:nobs_i),TempBasis}, K_max, 1);
    pcastrW_i  = pca3D({WindCoef(:,1:nobs_i),WindBasis}, K_max, 1);

    PCA_addreg = {pcastrS_i,pcastrP_i,pcastrT_i,pcastrW_i};

    %% Method 1: Mean Forecast

    MSurf_MF = mean(OzoneS(:,1:nobs_i),2) + SeasComp;

    %% Method 2: Naive Forecast

    MSurf_NF = OzoneS(:,nobs_i) + SeasComp;

    %% Method 3: FAR

    FuncPredCoeff = SFAR1(pcastrOz_i,Static_L_set(1),h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_FAR     = FPointEval + SeasComp;

    %% Method 4: PCA VAR

    FuncPredCoeff = FFM_VARX(pcastrOz_i,PCA_addreg,Static_L_set,Static_m_set,h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_PCA_VAR = FPointEval + SeasComp;

    %% Method 5: PCA KNN
    K_NNv       = zeros(KNN_max,1);
    pcastrKNN_i = pca3D({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, K_max, 1);

    for knn = 1:KNN_max

        FuncPredCoeff = FFM_KNN(pcastrKNN_i,Static_L_set(1),Static_m_set(1),knn,h);
        FuncPred      = fd(FuncPredCoeff,FEMBasOz);
        FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);

        K_NNv(knn,1) = mean((OzoneS(:,nobs_i)-FPointEval).^2);

    end

    [~,K_min_PCA] = min(K_NNv);

    FuncPredCoeff = FFM_KNN(pcastrOz_i,Static_L_set(1),Static_m_set(1),K_min_PCA,h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_PCA_KNN = FPointEval + SeasComp;

    %% Method 6: DS VAR
    FuncPredCoeff = FFM_VARX(dynstrOz_i,PCA_addreg,Dyn_L_set,Dyn_m_set,h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_DS_VAR  = FPointEval + SeasComp;

    %% Method 7: DS KNN
    K_NNv       = zeros(KNN_max,1);
    dynstrKNN_i = DynamScoresSurf({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, ...
                                   K_max, 1, q_dyn);

    for knn = 1:KNN_max

        FuncPredCoeff = FFM_KNN(dynstrKNN_i,Dyn_L_set(1),Dyn_m_set(1),knn,h);
        FuncPred      = fd(FuncPredCoeff,FEMBasOz);
        FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);

        K_NNv(knn,1) = mean((OzoneS(:,nobs_i)-FPointEval).^2);

    end

    [~,K_min_DS]  = min(K_NNv);
    FuncPredCoeff = FFM_KNN(dynstrOz_i,Dyn_L_set(1),Dyn_m_set(1),K_min_DS,h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_DS_KNN  = FPointEval + SeasComp;

    %% Method 8: MP VAR
    FPredGrid = DFFM_VARX_GRID(OzoneS(:,1:nobs_i), WindS(:,1:nobs_i), ...
                               SunS(:,1:nobs_i), PrecipS(:,1:nobs_i), ...
                               TempS(:,1:nobs_i), Static_L_set,Static_m_set,h);

    MSurf_MP_VAR = FPredGrid + SeasComp;

    %% Method 9: MP KNN
    K_NNv    = zeros(KNN_max,1);
    pcastr_i = FactorDecompFTSGrid(OzoneS(:,1:nobs_i-h), K_max, 1);

    for knn = 1:KNN_max

        FPredGrid = DFFM_KNN_GRID(pcastr_i,Static_L_set(1),Static_m_set(1),knn);
        K_NNv(knn,1) = mean((OzoneS(:,nobs_i)-FPredGrid).^2);

    end

    [~,K_min_MP] = min(K_NNv);

    pcastr_i     = FactorDecompFTSGrid(OzoneS(:,1:nobs_i), K_max, 1);
    FPredGrid    = DFFM_KNN_GRID(pcastr_i,Static_L_set(1),Static_m_set(1),K_min_MP);
    MSurf_MP_KNN = FPredGrid + SeasComp;

    %% Method 10: DS Random Forest

    FuncPredCoeff = FRandomForest(dynstrOz_i,Dyn_L_set(1),h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_RF      = FPointEval + SeasComp;

    %% Store forecast surfaces

    SurfArray = {TrueVal, MSurf_MF, MSurf_NF, MSurf_FAR, ...
                 MSurf_PCA_VAR, MSurf_PCA_KNN, MSurf_DS_VAR, MSurf_DS_KNN, ...
                 MSurf_MP_VAR, MSurf_MP_KNN, MSurf_RF};

    PollutedResults(pp).ForecastIndex = i;
    PollutedResults(pp).nobs_i        = nobs_i;
    PollutedResults(pp).Models        = Models;
    PollutedResults(pp).SurfArray     = SurfArray;
    PollutedResults(pp).K_min_PCA     = K_min_PCA;
    PollutedResults(pp).K_min_DS      = K_min_DS;
    PollutedResults(pp).K_min_MP      = K_min_MP;

    %% Plot forecast surfaces

    allVals = cell2mat(SurfArray);
    bottom  = min(allVals,[],'all');
    top     = max(allVals,[],'all');

    fig = figure;
    fig.Position = [100, 100, 1600, 850];

    for mm = 1:length(Models)

        subplot(3,4,mm)

        Season = Grid2Func(SurfArray{mm},DTOzone,CleanDTOz);

        hold on
            plot(RegionBord,'FaceColor','none');
            plot(Season,[],[],[],100);
        hold off

        axis equal tight
        colormap(jet)
        clim([bottom top])
        colorbar
        view(2)

        title(Models{mm});
        xlabel('Easting');
        ylabel('Northing');

    end

    sgtitle(['Polluted Day Forecasts: forecast-period index ',num2str(i)]);

    pdfName = fullfile(outFolder,['Figures_PollutedDay_',num2str(i),'.pdf']);
    pngName = fullfile(outFolder,['Figures_PollutedDay_',num2str(i),'.png']);
    figName = fullfile(outFolder,['Figures_PollutedDay_',num2str(i),'.fig']);

    exportgraphics(fig,pdfName, ...
        'BackgroundColor','none', ...
        'Resolution',300, ...
        'ContentType','vector');

    exportgraphics(fig,pngName, ...
        'BackgroundColor','white', ...
        'Resolution',300);

    savefig(fig,figName);

    fprintf('Saved polluted-day figures for index %d.\n', i);

end

toc

save(fullfile(outFolder,'Figures_PollutedDays_Results.mat'),'PollutedResults');

fprintf('\nFigures_PollutedDays completed successfully.\n');
fprintf('Outputs saved to: %s\n', outFolder);
