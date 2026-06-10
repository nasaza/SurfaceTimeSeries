%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forecasts Comparison
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
q_dyn           = 2; % cumulative autocovariance    
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
subplot(2,2,1)
    autocorr(dyn_scs.scr(:,1));
    title('1st Score Series');
subplot(2,2,2)
    autocorr(dyn_scs.scr(:,2));
    title('2nd Score Series');
subplot(2,2,3)
    autocorr(dyn_scs.scr(:,3));
    title('3rd Score Series');
subplot(2,2,4)
    autocorr(dyn_scs.scr(:,4));
    title('4th Score Series');


outFolder = fullfile(pwd, 'Outputs');
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end
exportgraphics(fg1, ['Outputs/ScreePlot.pdf'], 'BackgroundColor', 'none', 'Resolution', 300);