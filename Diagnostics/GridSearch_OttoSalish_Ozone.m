%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Otto-Salish Information Criterion: Ozone Dynamic Scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
close all;

%% Short Info
% This script jointly selects:
%
%   J = number of dynamic ozone scores;
%   m = lag order for those scores.
%
% It implements the BIC- and HQC-type information criteria proposed by
% Otto and Salish. The cumulative-autocovariance lag q_dyn is selected by the
% practitioner
%
% The analysis uses the first T_tr ozone surfaces, matching the initial
% training sample used in the forecasting exercise. Results, tables, and
% figures are saved directly to the Outputs folder.

%% Add Libraries and paths

scriptFolder = fileparts(mfilename('fullpath'));
rootFolder   = fileparts(scriptFolder);
addpath(fullfile(rootFolder,'AddFunc'));

dataFolder   = fullfile(rootFolder,'Data');
outFolder    = fullfile(rootFolder,'Outputs');

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end


%% User settings
% Practitioners may want to change only this block.

T_tr   = 200;  % Initial training-sample size
K_max  = 15;   % Maximum number of dynamic scores considered
p_max  = 7;    % Maximum VAR lag order considered
q_dyn  = 2;    % Lag order in cumulative autocovariance operator


%% Load functional ozone data created in Step 1

ftsFile = fullfile(dataFolder,'FTSs.mat');

if ~isfile(ftsFile)
    error(['Data/FTSs.mat was not found. Run ', ...
           'Step1_CreateSurfaceTimeSeries.m first.']);
end

load(ftsFile,'OzoneCoef','OzoneBasis');

if T_tr > size(OzoneCoef,2)
    error('T_tr exceeds the number of available ozone observations.');
end

if K_max > size(OzoneCoef,1)
    error('K_max exceeds the number of ozone basis coefficients.');
end


%% Run Otto-Salish grid search

FTSobj = {OzoneCoef(:,1:T_tr),OzoneBasis};
OS     = OttoSalishIC(FTSobj,K_max,p_max,q_dyn);

%% Display selected specifications

fprintf('\nOtto-Salish selection for ozone\n');
fprintf('Training observations: T = %d\n',T_tr);
fprintf('Dynamic-operator lag: q = %d\n\n',q_dyn);

fprintf('BIC selection: J = %d dynamic scores, m = %d VAR lags.\n', ...
        OS.BICSelection.J,OS.BICSelection.m);
fprintf('BIC cumulative eigenvalue mass at J: %.4f.\n', ...
        OS.EigenMass(OS.BICSelection.J));

fprintf('\nHQC selection: J = %d dynamic scores, m = %d VAR lags.\n', ...
        OS.HQCSelection.J,OS.HQCSelection.m);
fprintf('HQC cumulative eigenvalue mass at J: %.4f.\n', ...
        OS.EigenMass(OS.HQCSelection.J));

%% Save numerical results

fileStem = sprintf('OttoSalishGrid_Ozone_T%d_q%d',T_tr,q_dyn);

save(fullfile(outFolder,[fileStem,'.mat']),'OS','T_tr','K_max','p_max','q_dyn');
writetable(OS.GridTable,fullfile(outFolder,[fileStem,'.csv']));

SelectionTable = table( ...
    ["BIC";"HQC"], ...
    [OS.BICSelection.J;OS.HQCSelection.J], ...
    [OS.BICSelection.m;OS.HQCSelection.m], ...
    [OS.BICSelection.value;OS.HQCSelection.value], ...
    [OS.BICSelection.MSE;OS.HQCSelection.MSE], ...
    [OS.EigenMass(OS.BICSelection.J);OS.EigenMass(OS.HQCSelection.J)], ...
    'VariableNames',{'Criterion','SelectedScores','SelectedLags', ...
                     'CriterionValue','MSE','CumulativeEigenvalueMass'});

writetable(SelectionTable, ...
    fullfile(outFolder,[fileStem,'_Selections.csv']));

%% Plot MSE, BIC, and HQC grid surfaces

fig1 = figure(1);
fig1.Position = [100,100,1350,420];

TL = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% MSE surface
nexttile
imagesc(1:p_max,1:K_max,OS.MSE)
set(gca,'YDir','normal')
colorbar
xlabel('VAR lag order, m')
ylabel('Number of dynamic scores, J')
title('Otto-Salish MSE')

% BIC surface
nexttile
imagesc(1:p_max,1:K_max,OS.BIC)
set(gca,'YDir','normal')
colorbar
hold on
plot(OS.BICSelection.m,OS.BICSelection.J,'kp', ...
    'MarkerSize',13,'MarkerFaceColor','w')
hold off
xlabel('VAR lag order, m')
ylabel('Number of dynamic scores, J')
title(sprintf('BIC: J=%d, m=%d', ...
    OS.BICSelection.J,OS.BICSelection.m))

% HQC surface
nexttile
imagesc(1:p_max,1:K_max,OS.HQC)
set(gca,'YDir','normal')
colorbar
hold on
plot(OS.HQCSelection.m,OS.HQCSelection.J,'kp', ...
    'MarkerSize',13,'MarkerFaceColor','w')
hold off
xlabel('VAR lag order, m')
ylabel('Number of dynamic scores, J')
title(sprintf('HQC: J=%d, m=%d', ...
    OS.HQCSelection.J,OS.HQCSelection.m))

title(TL,sprintf('Ozone Dynamic-Score Selection: T=%d, q=%d',T_tr,q_dyn));

exportgraphics(fig1, ...
    fullfile(outFolder,[fileStem,'.pdf']), ...
    'BackgroundColor','none', ...
    'Resolution',300, ...
    'ContentType','vector');

exportgraphics(fig1, ...
    fullfile(outFolder,[fileStem,'.png']), ...
    'BackgroundColor','white', ...
    'Resolution',300);

savefig(fig1,fullfile(outFolder,[fileStem,'.fig']));

%% Plot cumulative dynamic eigenvalue mass
% This plot is descriptive only. Selection is based on BIC/HQC above.

fig2 = figure(2);
fig2.Position = [100,100,650,420];

plot(1:K_max,OS.EigenMass,'-o','LineWidth',1.2)
hold on
xline(OS.BICSelection.J,'--',sprintf('BIC J=%d',OS.BICSelection.J))
xline(OS.HQCSelection.J,':',sprintf('HQC J=%d',OS.HQCSelection.J))
hold off

grid on
xlabel('Number of dynamic scores, J')
ylabel('Cumulative dynamic eigenvalue mass')
title('Dynamic Eigenvalue Mass')
ylim([0,1.02])

exportgraphics(fig2, ...
    fullfile(outFolder,[fileStem,'_EigenvalueMass.pdf']), ...
    'BackgroundColor','none', ...
    'Resolution',300, ...
    'ContentType','vector');

exportgraphics(fig2, ...
    fullfile(outFolder,[fileStem,'_EigenvalueMass.png']), ...
    'BackgroundColor','white', ...
    'Resolution',300);

fprintf('\nGrid-search results saved to: %s\n',outFolder);
