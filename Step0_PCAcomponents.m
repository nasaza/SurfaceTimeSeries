clear all;
clc;

%% Short Info:
% Analyzez strucutre of the PCA

%% Add Libraries and Data
addpath AddFunc
addpath Data

%% Step 1: Read Data and Project
load('Data\SeasonAdjData');
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');
mkdir('Outputs');  % creates a folder to save all of the outputs 
Threshold = 0.75;

% project data on the plain (Gauss-Krüger Zone 3)
wgs84           = geocrs(4326);
proj            = projcrs(31467);
[LonOz,  LatOz] = projfwd(proj, LatOz, LonOz);
[x, y]          = projfwd(proj, ConstrReg(:,1), ConstrReg(:,2));
ConstrReg       = [y,x];

%% Step 2: Triangulation of the selected geographic area
% This step defines the geographical area for our analysis. In particular, 
% it creates a triangulation of the area for each variable and generates Figure 7 in the paper's appendix.

DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);
% [RemP, regP] = TriangStats(DTOzone, [ConstrReg(:,2),ConstrReg(:,1)], Threshold);
[OzoneFTS,OzoneCoef,OzoneBasis]    = Grid2Func(Ozone,DTOzone,CleanDTOz);
%% Plots 
% This block obtains additional information about the characteristics of our surface data.

K_max       = 15;    
pcastrOzone = pca3D(OzoneFTS(1:200), K_max, 1);   
% plot a scree plot of explaind variance for an intuition
h = findobj('type','figure');
n = length(h);

% scree plots

figure(2);
subplot(1,2,1)
    plot(pcastrOzone.varprop(1:K_max));
    title('Ozone');
subplot(1,2,2)    
    plot(pcastrOzone.values(1:K_max));

fh=figure(3);
subplot(2,2,1)
    autocorr(pcastrOzone.pcascr(:,1));
    title('1st Score Series');
subplot(2,2,2)
    autocorr(pcastrOzone.pcascr(:,2));
    title('2nd Score Series');
subplot(2,2,3)
    autocorr(pcastrOzone.pcascr(:,3));
    title('3rd Score Series');
subplot(2,2,4)
    autocorr(pcastrOzone.pcascr(:,4));
    title('4th Score Series');

exportgraphics(fh, ['Outputs/AutocorrScores.pdf'], 'BackgroundColor', 'none', 'Resolution', 300);