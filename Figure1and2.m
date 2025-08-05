clear all;
clc;

%% Short Info:
% This file creates functiona time sereis

%% Add Libraries and Data
addpath AddFunc
addpath Data

%% Step 1: Read Data

load('Data\SeasonAdjData');
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');
Threshold = 0.75;

%% Project lat-lon degrees on a flat surface (in meters, Gauss-Krüger Zone 3)
wgs84           = geocrs(4326);
proj            = projcrs(31467);
[LonOz,  LatOz] = projfwd(proj, LatOz, LonOz);
[x, y]          = projfwd(proj, ConstrReg(:,1), ConstrReg(:,2));
ConstrReg       = [y,x];



%% Build surface
DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);

[OzoneFTS,OzoneCoef,OzoneBasis]    = Grid2Func(Ozone,DTOzone,CleanDTOz);
border_poly = polyshape(ConstrReg(:,2), ConstrReg(:,1));


%% Figure 1

Day         = 151;
bottom      = 40;
top         = 100;
fig1        = figure(1);
fig1.Position = [100, 100, 900, 400]

subplot(1,2,1) 
    hold on
    plot(border_poly,'FaceColor', 'none');
    scatter(LonOz,LatOz,20,Ozone(:,Day),'filled'); % For Figure 1
    axis equal tight
    hold off
    colormap(jet)
    colorbar;
    caxis([bottom top]);  
    title('Grid Data');
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    xlabel('Easting (meters)');
    ylabel('Northing (meters)');

subplot(1,2,2)    
    hold on
    plot(border_poly,'FaceColor', 'none');
    plot(OzoneFTS(Day),[],[],[],100);
    view(2);
    axis equal tight
    hold off
    colorbar
    colormap(jet)    
    caxis([40 top]);  
    title('Surface');
    xlabel('Easting (meters)');
    ylabel('Northing (meters)');


% Adjust the figure properties to minimize white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

exportgraphics(fig1,['Outputs/Figure1_m',num2str(Threshold),'.pdf'],'BackgroundColor','none','Resolution',300)  



fig1        = figure(2);
fig1.Position = [100, 100, 900, 400]

hold on
    plot(border_poly,'FaceColor', 'none');
    % scatter(LonOz,LatOz,20,Ozone(:,Day),'filled'); % For Figure 1
    triplot(DTOzone(logical(CleanDTOz),:),DTOzone.Points(:,1), DTOzone.Points(:,2));   % for Figure 2
    scatter(LonOz,LatOz,"r","."); % for Figure 2
    axis equal tight
    hold off
    % colormap(jet)
    % colorbar;
    % caxis([bottom top]);  
    % title('Grid Data');
    % ax = gca;
    % outerpos = ax.OuterPosition;
    % ti = ax.TightInset; 
    % left = outerpos(1) + ti(1);
    % bottom = outerpos(2) + ti(2);
    % ax_width = outerpos(3) - ti(1) - ti(3);
    % ax_height = outerpos(4) - ti(2) - ti(4);
    % ax.Position = [left bottom ax_width ax_height];
    xlabel('Easting (meters)');
    ylabel('Northing (meters)');
% Adjust the figure properties to minimize white space


exportgraphics(fig1,['Outputs/Figure2_m',num2str(Threshold),'.pdf'],'BackgroundColor','none','Resolution',300)   


 


