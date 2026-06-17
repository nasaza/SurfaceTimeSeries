%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional Diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Add Libraries and path
scriptFolder = fileparts(mfilename('fullpath'));
rootFolder   = fileparts(scriptFolder);
addpath(fullfile(rootFolder,'AddFunc'));

dataFolder = fullfile(rootFolder,'Data');
outFolder  = fullfile(rootFolder,'Outputs');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% User settings
K_max     = 15;
T_tr      = 200;
q_dyn     = 2;     % cumulative autocovariance lag
nShow     = 6;     % number of scores shown in ACF/PACF
nLoad     = 3;     % number of loading surfaces to plot

%% Step 1: Read Data
% This step loads seasonally adjusted ozone concentration data observed
% at different stations, the surface data (created in Step 1), and the
% geographic borders of Germany.

load(fullfile(dataFolder,'SeasonAdjData.mat'),'LonOz','LatOz');
ftsFile = fullfile(dataFolder,'FTSs.mat');
if ~isfile(ftsFile)
    error(['Data/FTSs.mat was not found. Run ', ...
           'Step1_CreateSurfaceTimeSeries.m first.']);
end
load(ftsFile,'OzoneCoef','OzoneBasis');

ConstrReg = readmatrix(fullfile(dataFolder,'GeoConstraints','DE_Constraints.csv'));

%% Projected Coordinates
% Project data to Gauss-Krüger Zone 3.

proj = projcrs(31467);

[LonOz,LatOz] = projfwd(proj,LatOz,LonOz);
[x,y]         = projfwd(proj,ConstrReg(:,1),ConstrReg(:,2));
ConstrReg     = [y,x];

%% Static and dynamic scores

dyn_scs  = DynamScoresSurf({OzoneCoef(:,1:T_tr),OzoneBasis}, K_max, 1, q_dyn);
stat_scs = pca3D({OzoneCoef(:,1:T_tr),OzoneBasis}, K_max, 1);

%% Figure 1: Scree plots

fg1 = figure(1);
fg1.Position = [100,100,900,600];

subplot(2,2,1)
plot(stat_scs.varprop(1:K_max),'-o')
title('Static Scores: Cumulative Variance Proportion')
xlabel('Component')
ylabel('Cumulative proportion')

subplot(2,2,2)
plot(stat_scs.values(1:K_max),'-o')
title('Static Scores: Eigenvalues')
xlabel('Component')
ylabel('Eigenvalue')

subplot(2,2,3)
plot(dyn_scs.varprop(1:K_max),'-o')
title('Dynamic Scores: Cumulative Dynamic Proportion')
xlabel('Component')
ylabel('Cumulative proportion')

subplot(2,2,4)
plot(dyn_scs.values(1:K_max),'-o')
title('Dynamic Scores: Eigenvalues')
xlabel('Component')
ylabel('Eigenvalue')

exportgraphics(fg1, ...
    fullfile(outFolder,'Diagnostics_ScreePlots.pdf'), ...
    'BackgroundColor','none', ...
    'Resolution',300);

%% Figure 2: ACF of first dynamic scores
fg2 = figure(2);
fg2.Position = [100,100,1100,650];

for k = 1:nShow
    subplot(2,3,k)
    autocorr(dyn_scs.pcascr(:,k));
    title(sprintf('Dynamic Score %d: ACF',k))
end

exportgraphics(fg2, ...
    fullfile(outFolder,'Diagnostics_DynamicScores_ACF.pdf'), ...
    'BackgroundColor','none', ...
    'Resolution',300);

%% Figure 3: PACF of first dynamic scores
fg3 = figure(3);
fg3.Position = [100,100,1100,650];

for k = 1:nShow
    subplot(2,3,k)
    parcorr(dyn_scs.pcascr(:,k));
    title(sprintf('Dynamic Score %d: PACF',k))
end

exportgraphics(fg3, ...
    fullfile(outFolder,'Diagnostics_DynamicScores_PACF.pdf'), ...
    'BackgroundColor','none', ...
    'Resolution',300);

%% Figure 4: Loading surfaces comparison: PCA vs Dynamic Scores
% We plot the first 3 PCA loading surfaces and the first 3 dynamic loading
% surfaces. Since loading signs are arbitrary, we align them for clearer
% comparison:
%
% (i)  each PCA loading is oriented so that its largest absolute value at
%      the ozone locations is positive;
% (ii) each dynamic loading is then oriented to have positive inner product
%      with the corresponding PCA loading, based on station evaluations.

Region     = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord = polyshape(Region);

% Extract loading coefficients
stat_load_coef = getcoef(stat_scs.pcafd);
dyn_load_coef  = getcoef(dyn_scs.pcafd);

% Containers for sign-aligned fd objects and station evaluations
stat_load_fd   = cell(1,nLoad);
dyn_load_fd    = cell(1,nLoad);
stat_vals_cell = cell(1,nLoad);
dyn_vals_cell  = cell(1,nLoad);

for k = 1:nLoad

    %--- Static loading
    coef_stat_k = stat_load_coef(:,k);
    fd_stat_k   = fd(coef_stat_k,OzoneBasis);
    vals_stat_k = eval_FEM_fd(LonOz,LatOz,fd_stat_k);

    % orient so largest absolute station value is positive
    [~,idxMax] = max(abs(vals_stat_k));
    if vals_stat_k(idxMax) < 0
        coef_stat_k = -coef_stat_k;
        vals_stat_k = -vals_stat_k;
        fd_stat_k   = fd(coef_stat_k,OzoneBasis);
    end

    stat_load_fd{k}   = fd_stat_k;
    stat_vals_cell{k} = vals_stat_k;

    %--- Dynamic loading
    coef_dyn_k = dyn_load_coef(:,k);
    fd_dyn_k   = fd(coef_dyn_k,OzoneBasis);
    vals_dyn_k = eval_FEM_fd(LonOz,LatOz,fd_dyn_k);

    % orient to match corresponding PCA loading
    if sum(vals_stat_k .* vals_dyn_k) < 0
        coef_dyn_k = -coef_dyn_k;
        vals_dyn_k = -vals_dyn_k;
        fd_dyn_k   = fd(coef_dyn_k,OzoneBasis);
    end

    dyn_load_fd{k}   = fd_dyn_k;
    dyn_vals_cell{k} = vals_dyn_k;
end

% Common symmetric color scale across all 6 loadings
vals_load = [];
for k = 1:nLoad
    vals_load = [vals_load; stat_vals_cell{k}; dyn_vals_cell{k}];
end
cmax_load = max(abs(vals_load),[],'all');

% Plot
fg4 = figure(4);
fg4.Position = [100,100,1300,750];

for k = 1:nLoad
    subplot(2,3,k)
    hold on
    plot(RegionBord,'FaceColor','none');
    plot(stat_load_fd{k},[],[],[],100);
    hold off
    axis equal tight
    view(2)
    colormap(jet)
    % clim([-cmax_load cmax_load])
    colorbar
    title(sprintf('PCA: Loading %d',k))
    xlabel('Easting')
    ylabel('Northing')
end

for k = 1:nLoad
    subplot(2,3,nLoad+k)
    hold on
    plot(RegionBord,'FaceColor','none');
    plot(dyn_load_fd{k},[],[],[],100);
    hold off
    axis equal tight
    view(2)
    colormap(jet)
    % clim([-cmax_load cmax_load])
    colorbar
    title(sprintf('DS: Loading %d',k))
    xlabel('Easting')
    ylabel('Northing')
end

exportgraphics(fg4, ...
    fullfile(outFolder,'Diagnostics_LoadingSurfaces_PCA_DS.pdf'), ...
    'BackgroundColor','none', ...
    'Resolution',300, ...
    'ContentType','vector');

exportgraphics(fg4, ...
    fullfile(outFolder,'Diagnostics_LoadingSurfaces_PCA_DS.png'), ...
    'BackgroundColor','white', ...
    'Resolution',300);

fprintf('\nAdditional diagnostics completed successfully.\n');
fprintf('Outputs saved to: %s\n', outFolder);