%% Add Libraries and Data
addpath Data/
addpath AddFunc
load('Data\SeasonAdjData');
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');

%% Preparation of the FEM and data
DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = ones(326,1);
OutDTs      = [2,5,20,35,69,70,76,80,81,128,131,133,134,135,159,...
                        173,239,248,261,308,312,315,322,323,324];

CleanDTOz(OutDTs)   = zeros(length(OutDTs),1);     

%% PCA
[OzoneFTS,OzoneCoef,OzoneBasis]    = Grid2Func(OzoneS,DTOzone,CleanDTOz);
K_max       = 15;
pcastrOzone = pca3D(OzoneFTS, K_max, 1);    
bottom      = min(OzoneS,[],'all');
top         = max(OzoneS,[],'all');

Region      = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord  = polyshape(Region);

%% Figure
fig = figure;
subplot(1,3,1)
    hold on
    plot(RegionBord,'FaceColor', 'none');
    plot(pcastrOzone.pcafd(1),[],[],[],100);
    view(2)
    hold off
    colormap(jet);
    xlabel('');
    ylabel('');
subplot(1,3,2)
    hold on
    plot(RegionBord,'FaceColor', 'none');
    plot(pcastrOzone.pcafd(2),[],[],[],100);
    view(2)
    hold off
    colormap(jet);
    xlabel('');
    ylabel('');    
subplot(1,3,3)
    hold on
    plot(RegionBord,'FaceColor', 'none');
    plot(pcastrOzone.pcafd(3),[],[],[],100);
    view(2)
    hold off
    colormap(jet);
    %caxis([bottom top])
    %colorbar;
    xlabel('');
    ylabel('');    


fig.Position = [100, 100, 700, 250];    
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


exportgraphics(fig,['Outputs/PCAComponets.pdf'],'BackgroundColor','none','Resolution',300,'ContentType', 'vector')  