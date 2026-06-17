%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select FPCA Dimension and VAR Order by Aue et al. fFPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Short Info
% Jointly selects:
%   d - number of static FPCA scores;
%   p - VAR lag order;
% using the functional final prediction error criterion of
% Aue, Dubart Norinho and Hoermann (2015).

%% Add Libraries and path

scriptFolder = fileparts(mfilename('fullpath'));
rootFolder   = fileparts(scriptFolder);
addpath(fullfile(rootFolder,'AddFunc'));

dataFolder   = fullfile(rootFolder,'Data');
outFolder    = fullfile(rootFolder,'Outputs');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% User settings
H      = 165;
d_max  = 10;
p_max  = 7;

%% Read functional data

ftsFile = fullfile(dataFolder,'FTSs.mat');

if ~isfile(ftsFile)
    error(['Data/FTSs.mat was not found. Run ', ...
           'Step1_CreateSurfaceTimeSeries.m first.']);
end

load(ftsFile,'OzoneCoef','OzoneBasis');

T    = size(OzoneCoef,2);
T_tr = T-H;

if T_tr <= 0
    error('H must be smaller than the number of observations.');
end

if d_max > size(OzoneCoef,1)
    error('d_max exceeds the number of ozone basis coefficients.');
end


%% Static FPCA on initial training sample
stat_scs = pca3D({OzoneCoef(:,1:T_tr),OzoneBasis},d_max,1);

%% fFPE selection
FFPE = AueFFPE(stat_scs,d_max,p_max);

fprintf('\nAue et al. fFPE selection\n');
fprintf('Selected FPCA dimension d: %d\n',FFPE.selectedD);
fprintf('Selected VAR order p:      %d\n',FFPE.selectedP);
fprintf('Minimum fFPE:              %.6f\n',FFPE.minimumFFPE);

%% Save numerical output
Results.T_tr = T_tr;
Results.FFPE = FFPE;

save(fullfile(outFolder,'FPCA_AueFFPE_Selection.mat'),'Results');

Rows = [];
for ip = 1:length(FFPE.pGrid)
    for id = 1:length(FFPE.dGrid)
        Rows = [Rows;
                FFPE.pGrid(ip), ...
                FFPE.dGrid(id), ...
                FFPE.TraceSigma(ip,id), ...
                FFPE.TailMass(id), ...
                FFPE.FFPE(ip,id)];
    end
end

GridTable = array2table(Rows, ...
    'VariableNames',{'VAR_Order','FPCA_Dimension', ...
                     'ResidualTrace','TruncationMass','fFPE'});

writetable(GridTable,fullfile(outFolder,'FPCA_AueFFPE_Grid.csv'));

%% Plot criterion grid
fig = figure;
fig.Position = [100,100,850,550];

imagesc(FFPE.dGrid,FFPE.pGrid,FFPE.FFPE);
set(gca,'YDir','normal');
colorbar

xlabel('Number of FPCA scores, d');
ylabel('VAR order, p');
title('Aue et al. Functional Final Prediction Error');

hold on
plot(FFPE.selectedD,FFPE.selectedP,'kp', ...
     'MarkerSize',14,'MarkerFaceColor','w','LineWidth',1.5);
hold off

exportgraphics(fig,fullfile(outFolder,'FPCA_AueFFPE_Selection.pdf'), ...
    'BackgroundColor','none','Resolution',300,'ContentType','vector');

exportgraphics(fig,fullfile(outFolder,'FPCA_AueFFPE_Selection.png'), ...
    'BackgroundColor','white','Resolution',300);

fprintf('Outputs saved to: %s\n',outFolder);
