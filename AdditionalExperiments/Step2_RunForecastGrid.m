%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Forecasting Grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Add Libraries

addpath AddFunc
addpath Data

%% Output folder

outFolder = fullfile(pwd,'Outputs','ForecastGrid');

if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% Specification grid

% d_grid = 1:3;   % lags in cumulative autocovariance operator
% L_grid = 3:5;   % number of scores/components
% m_grid = 1:3;   % number of lags in forecasting step
% 
% h      = 1;     % forecast horizon

d_grid  = 1:3;   % lags in cumulative autocovariance operator
Lm_grid = [3 1;
           4 2;
           5 3]; % paired (L,m) specifications

h = 1;

%% Run all specifications

counter = 0;
Summary = struct();

% for d_dyn = d_grid
% 
%     for L = L_grid
% 
%         % for m = m_grid
% 
%             counter = counter + 1;
% 
%             fprintf('\n====================================================\n');
%             fprintf('Running specification %d: d=%d, L=%d, m=%d, h=%d\n', ...
%                     counter, d_dyn, L, m, h);
%             fprintf('====================================================\n');
% 
%             try
% 
%                 Results = RunForecastSpec(d_dyn, L, m, h, outFolder);
% 
%                 Summary(counter).d_dyn   = d_dyn;
%                 Summary(counter).L       = L;
%                 Summary(counter).m       = m;
%                 Summary(counter).h       = h;
%                 Summary(counter).Models  = Results.Models;
%                 Summary(counter).MeanMSE = Results.MeanMSE;
%                 Summary(counter).RunTime = Results.RunTime;
%                 Summary(counter).status  = 'ok';
%                 Summary(counter).message = '';
% 
%             catch ME
% 
%                 warning('Specification failed: d=%d, L=%d, m=%d. Error: %s', ...
%                          d_dyn, L, m, ME.message);
% 
%                 Summary(counter).d_dyn   = d_dyn;
%                 Summary(counter).L       = L;
%                 Summary(counter).m       = m;
%                 Summary(counter).h       = h;
%                 Summary(counter).Models  = {};
%                 Summary(counter).MeanMSE = NaN;
%                 Summary(counter).RunTime = NaN;
%                 Summary(counter).status  = 'failed';
%                 Summary(counter).message = ME.message;
% 
%             end
% 
%             save(fullfile(outFolder,'ForecastGrid_Summary.mat'),'Summary');
% 
%         % end
% 
%     end
% 
% end
counter = 0;
Summary = struct();

for d_dyn = d_grid

    for lm = 1:size(Lm_grid,1)

        L = Lm_grid(lm,1);
        m = Lm_grid(lm,2);

        counter = counter + 1;

        fprintf('\n====================================================\n');
        fprintf('Running specification %d: d=%d, L=%d, m=%d, h=%d\n', ...
                counter, d_dyn, L, m, h);
        fprintf('====================================================\n');

        try

            Results = RunForecastSpec(d_dyn, L, m, h, outFolder);

            Summary(counter).d_dyn   = d_dyn;
            Summary(counter).L       = L;
            Summary(counter).m       = m;
            Summary(counter).h       = h;
            Summary(counter).Models  = Results.Models;
            Summary(counter).MeanMSE = Results.MeanMSE;
            Summary(counter).RunTime = Results.RunTime;
            Summary(counter).status  = 'ok';
            Summary(counter).message = '';

        catch ME

            warning('Specification failed: d=%d, L=%d, m=%d. Error: %s', ...
                     d_dyn, L, m, ME.message);

            Summary(counter).d_dyn   = d_dyn;
            Summary(counter).L       = L;
            Summary(counter).m       = m;
            Summary(counter).h       = h;
            Summary(counter).Models  = {};
            Summary(counter).MeanMSE = NaN;
            Summary(counter).RunTime = NaN;
            Summary(counter).status  = 'failed';
            Summary(counter).message = ME.message;

        end

        save(fullfile(outFolder,'ForecastGrid_Summary.mat'),'Summary');

    end

end
%% Convert successful results to table

Rows = table();

for ii = 1:length(Summary)

    if strcmp(Summary(ii).status,'ok')

        tmp = array2table(Summary(ii).MeanMSE, ...
            'VariableNames', matlab.lang.makeValidName(Summary(ii).Models));

        tmp.d_dyn  = Summary(ii).d_dyn;
        tmp.L      = Summary(ii).L;
        tmp.m      = Summary(ii).m;
        tmp.h      = Summary(ii).h;
        tmp.RunTime = Summary(ii).RunTime;

        tmp = movevars(tmp, {'d_dyn','L','m','h','RunTime'}, 'Before', 1);

        Rows = [Rows; tmp];

    end

end

writetable(Rows, fullfile(outFolder,'ForecastGrid_MeanMSE.csv'));

disp(Rows)

%% Find best specification by DS VAR and DS KNN

if ~isempty(Rows)

    if ismember('DSVAR', Rows.Properties.VariableNames)

        [~,idxBestDSVAR] = min(Rows.DSVAR);

        fprintf('\nBest DS VAR specification:\n');
        disp(Rows(idxBestDSVAR,:));

    end

    if ismember('DSKNN', Rows.Properties.VariableNames)

        [~,idxBestDSKNN] = min(Rows.DSKNN);

        fprintf('\nBest DS KNN specification:\n');
        disp(Rows(idxBestDSKNN,:));

    end

end