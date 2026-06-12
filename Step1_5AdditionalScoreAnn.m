%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aditional Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Add Libraries
addpath AddFunc
addpath Data

%% Step 1: Read Data
% This step loads seasonally adjusted ozone concentration data observed
% at different stations, the surface data (created in Step 1), and the
% geographic borders of Germany

load('Data\SeasonAdjData'); % Data on the grid seasonaly addjusted
load('Data\FTSs');          % Data given in the surface form (Step 1 has to be executed for that)
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');

% Projected Coordinates
Threshold = 0.75;
% project data on the plain (Gauss-Krüger Zone 3)
wgs84           = geocrs(4326);
proj            = projcrs(31467);
[LonOz,  LatOz] = projfwd(proj, LatOz, LonOz);
[x, y]          = projfwd(proj, ConstrReg(:,1), ConstrReg(:,2));
ConstrReg       = [y,x];

[~,T]           = size(OzoneS);
K_max           = 15;
T_tr            = 200;
OzoneFTS        = fd(OzoneCoef(:,1:T_tr),OzoneBasis);
    
%% Static and dynamic scores
q_dyn           = 3; % cumulative autocovariance    
dyn_scs         = DynamScoresSurf({OzoneCoef(:,1:T_tr),OzoneBasis}, K_max, 1, q_dyn); 
stat_scs        = pca3D({OzoneCoef(:,1:T_tr),OzoneBasis}, K_max, 1); 

% scree plots
fg1 = figure(1);
    plot(dyn_scs.values(1:K_max));
    title('Scree plot');

% comapre with static scores
figure(2);
subplot(2,2,1)
    plot(stat_scs.varprop(1:K_max));
    title('Static Scores: Prop');
subplot(2,2,2)    
    plot(stat_scs.values(1:K_max));
    title('Static Scores: Vals');
subplot(2,2,3)
    plot(dyn_scs.varprop(1:K_max));
    title('Dynamic Scores: Prop');
subplot(2,2,4)    
    plot(dyn_scs.values(1:K_max));
    title('Dynamic Scores: Vals');


% Autocovariances     
fg2=figure(3);
subplot(2,3,1)
    autocorr(dyn_scs.pcascr(:,1));
    title('1st Score Series');
subplot(2,3,2)
    autocorr(dyn_scs.pcascr(:,2));
    title('2nd Score Series');
subplot(2,3,3)
    autocorr(dyn_scs.pcascr(:,3));
    title('3rd Score Series');
subplot(2,3,4)
    autocorr(dyn_scs.pcascr(:,4));
    title('4th Score Series');
subplot(2,3,5)
    autocorr(dyn_scs.pcascr(:,5));
    title('5th Score Series');
subplot(2,3,6)
    autocorr(dyn_scs.pcascr(:,6));
    title('6th Score Series');    

% Partial Autocovariances     
figure(4);
subplot(2,3,1)
    parcorr(dyn_scs.pcascr(:,1));
    title('1st Score Series');
subplot(2,3,2)
    parcorr(dyn_scs.pcascr(:,2));
    title('2nd Score Series');
subplot(2,3,3)
    parcorr(dyn_scs.pcascr(:,3));
    title('3rd Score Series');
subplot(2,3,4)
    parcorr(dyn_scs.pcascr(:,4));
    title('4th Score Series');
subplot(2,3,5)
    parcorr(dyn_scs.pcascr(:,5));
    title('5th Score Series');
subplot(2,3,6)
    parcorr(dyn_scs.pcascr(:,6));
    title('6th Score Series');    



outFolder = fullfile(pwd, 'Outputs');
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end
exportgraphics(fg1, ['Outputs/ScreePlot.pdf'], 'BackgroundColor', 'none', 'Resolution', 300);


%% Loading surfaces comparison: PCA vs Dynamic Scores

% Prepare Germany border
Region     = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord = polyshape(Region);

% Extract loading coefficients
stat_load_coef = getcoef(stat_scs.pcafd);
dyn_load_coef  = getcoef(dyn_scs.pcafd);

% First two PCA loading surfaces
stat_load_1 = fd(stat_load_coef(:,1), OzoneBasis);
stat_load_2 = fd(stat_load_coef(:,2), OzoneBasis);

% First two dynamic loading surfaces
dyn_load_1  = fd(dyn_load_coef(:,1), OzoneBasis);
dyn_load_2  = fd(dyn_load_coef(:,2), OzoneBasis);

% Common color scale based on values at ozone stations
vals_load = [ ...
    eval_FEM_fd(LonOz,LatOz,stat_load_1); ...
    eval_FEM_fd(LonOz,LatOz,stat_load_2); ...
    eval_FEM_fd(LonOz,LatOz,dyn_load_1); ...
    eval_FEM_fd(LonOz,LatOz,dyn_load_2)];

cmax_load = max(abs(vals_load),[],'all');

% Plot
fg5 = figure(5);
fg5.Position = [100, 100, 1100, 700];

subplot(2,2,1)
    hold on
    plot(RegionBord, 'FaceColor', 'none');
    plot(stat_load_1, [], [], [], 100);
    hold off
    axis equal tight
    view(2)
    colormap(jet)
    clim([-cmax_load cmax_load])
    colorbar
    title('PCA: 1st loading')

subplot(2,2,2)
    hold on
    plot(RegionBord, 'FaceColor', 'none');
    plot(stat_load_2, [], [], [], 100);
    hold off
    axis equal tight
    view(2)
    colormap(jet)
    clim([-cmax_load cmax_load])
    colorbar
    title('PCA: 2nd loading')

subplot(2,2,3)
    hold on
    plot(RegionBord, 'FaceColor', 'none');
    plot(dyn_load_1, [], [], [], 100);
    hold off
    axis equal tight
    view(2)
    colormap(jet)
    clim([-cmax_load cmax_load])
    colorbar
    title('DS: 1st loading')

subplot(2,2,4)
    hold on
    plot(RegionBord, 'FaceColor', 'none');
    plot(dyn_load_2, [], [], [], 100);
    hold off
    axis equal tight
    view(2)
    colormap(jet)
    clim([-cmax_load cmax_load])
    colorbar
    title('DS: 2nd loading')

exportgraphics(fg5, ...
    fullfile(outFolder,'LoadingSurfaces_PCA_DS.pdf'), ...
    'BackgroundColor', 'none', ...
    'Resolution', 300, ...
    'ContentType', 'vector');