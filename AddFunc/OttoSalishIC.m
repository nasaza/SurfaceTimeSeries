function IC = OttoSalishIC(FTSobj,K_max,p_max,q_dyn)
% OTTOSALISHIC Joint selection of dynamic-score dimension and VAR lag order.
%
% Implements the BIC- and HQC-type criteria proposed by Otto and Salish for
% selecting the number of dynamic factors/scores J and the VAR order m:
%
%   BIC(J,m) = log(MSE(J,m)) + J*m*log(T)/T
%   HQC(J,m) = log(MSE(J,m)) + 2*J*m*log(log(T))/T
%
% where
%
%   MSE(J,m) = trace(Sigma_eta(J,m)) + V_epsilon(J).
%
% Sigma_eta(J,m) is the covariance matrix of the residuals from a VAR(m)
% fitted to the first J dynamic scores, and V_epsilon(J) is the average
% integrated squared norm of the residual curves after projection on the
% first J dynamic loading functions.
%
% IMPORTANT:
%   q_dyn determines the cumulative autocovariance operator used to obtain
%   the dynamic scores. The Otto-Salish criterion selects J and m conditional
%   on q_dyn; it does not provide a criterion for selecting q_dyn itself.
%
% Inputs:
%   FTSobj - FDA fd object or cell {coefficient matrix, basis object}.
%            The coefficient matrix must be nbasis x T.
%   K_max  - maximum number of dynamic scores considered.
%   p_max  - maximum VAR lag order considered.
%   q_dyn  - maximum lag in the cumulative autocovariance operator.
%
% Output structure IC:
%   .MSE                    K_max x p_max MSE surface
%   .BIC                    K_max x p_max BIC surface
%   .HQC                    K_max x p_max HQC surface
%   .ScoreInnovationVar     trace of estimated VAR residual covariance
%   .ResidualCurveVar       functional residual variance for each J
%   .DesignRank             rank of each VAR design matrix
%   .BICSelection           selected J, m, and criterion value
%   .HQCSelection           selected J, m, and criterion value
%   .GridTable              long-form table with all grid results
%   .DynamicDecomposition   output from DynamScoresSurf
%   .EigenMass              cumulative dynamic eigenvalue mass
%   .T, .K_max, .p_max, .q_dyn
%
% The VAR is estimated by conditional least squares without an intercept,
% matching the score VAR formulation in the Otto-Salish manuscript. The
% dynamic scores are centered by construction.

%% Input checks

if nargin < 4
    q_dyn = 1;
end

validateattributes(K_max,{'numeric'},{'scalar','integer','positive'});
validateattributes(p_max,{'numeric'},{'scalar','integer','positive'});
validateattributes(q_dyn,{'numeric'},{'scalar','integer','positive'});

if isa(FTSobj,'fd')
    fdobj = FTSobj;
elseif iscell(FTSobj) && numel(FTSobj)==2
    fdobj = fd(FTSobj{1},FTSobj{2});
else
    error('FTSobj must be an fd object or a cell {coefficients,basis}.');
end

coef_FTS = getcoef(fdobj)';                  % T x nbasis
[T,nbasis] = size(coef_FTS);

if q_dyn >= T
    error('q_dyn must be smaller than the sample size T.');
end

if K_max > nbasis
    error('K_max cannot exceed the number of basis functions.');
end

if p_max >= T
    error('p_max must be smaller than the sample size T.');
end

%% Dynamic decomposition
% The data must be centered for the score VAR and residual-curve variance.

DS = DynamScoresSurf(FTSobj,K_max,1,q_dyn);

Scores = DS.pcascr(:,1:K_max);               % T x K_max
Theta  = getcoef(DS.pcafd);                   % nbasis x K_max
muCoef = getcoef(DS.meanfd)';                 % 1 x nbasis

coefCentered = coef_FTS - muCoef;

%% FEM inner-product matrix

BasisFD = getbasis(fdobj);
G = FEMMassMatrix(BasisFD);
G = (G + G')/2;

%% Functional residual variance V_epsilon(J)
% epsilon_t^(J) = Y_t - mean - sum_{j=1}^J score_{j,t} psi_j.

ResidualCurveVar = NaN(K_max,1);

for J = 1:K_max

    fittedCoef = Scores(:,1:J) * Theta(:,1:J)';
    residCoef  = coefCentered - fittedCoef;

    % Row-wise squared functional norms: r_t' G r_t.
    residNormSq = sum((residCoef*G).*residCoef,2);
    ResidualCurveVar(J) = mean(residNormSq);

end

%% Grid over score dimension J and VAR order m

MSE                = NaN(K_max,p_max);
BIC                = NaN(K_max,p_max);
HQC                = NaN(K_max,p_max);
ScoreInnovationVar = NaN(K_max,p_max);
DesignRank         = NaN(K_max,p_max);
EffectiveSample    = NaN(K_max,p_max);

for J = 1:K_max

    F = Scores(:,1:J);

    for m = 1:p_max

        nEff = T-m;
        nReg = J*m;

        EffectiveSample(J,m) = nEff;

        % Avoid an unidentified or saturated conditional LS problem.
        if nEff <= nReg
            continue
        end

        Y = F(m+1:T,:);
        X = zeros(nEff,nReg);

        for lag = 1:m
            cols = (lag-1)*J + (1:J);
            X(:,cols) = F(m+1-lag:T-lag,:);
        end

        DesignRank(J,m) = rank(X);

        % The manuscript assumes an invertible lag-score covariance matrix.
        % Mark rank-deficient specifications as invalid rather than using a
        % generalized inverse that could artificially reduce the MSE.
        if DesignRank(J,m) < nReg
            continue
        end

        % Conditional least-squares VAR estimate. B is (J*m) x J.
        Bhat = X \ Y;
        Ehat = Y - X*Bhat;

        SigmaEta = (Ehat'*Ehat)/nEff;
        ScoreInnovationVar(J,m) = trace(SigmaEta);

        MSE(J,m) = ScoreInnovationVar(J,m) + ResidualCurveVar(J);

        if isfinite(MSE(J,m)) && MSE(J,m) > 0
            BIC(J,m) = log(MSE(J,m)) + J*m*log(T)/T;
            HQC(J,m) = log(MSE(J,m)) + 2*J*m*log(log(T))/T;
        end

    end
end

%% Select minima

BICtmp = BIC;
BICtmp(~isfinite(BICtmp)) = Inf;
[minBIC,idxBIC] = min(BICtmp(:));

HQCtmp = HQC;
HQCtmp(~isfinite(HQCtmp)) = Inf;
[minHQC,idxHQC] = min(HQCtmp(:));

if isinf(minBIC) || isinf(minHQC)
    error(['No valid (J,m) specification was found. Reduce K_max or p_max, ', ...
           'or inspect the rank of the dynamic-score lag matrices.']);
end

[J_BIC,m_BIC] = ind2sub(size(BIC),idxBIC);
[J_HQC,m_HQC] = ind2sub(size(HQC),idxHQC);

%% Long-form output table

[Jgrid,mgrid] = ndgrid(1:K_max,1:p_max);

GridTable = table( ...
    Jgrid(:), ...
    mgrid(:), ...
    MSE(:), ...
    BIC(:), ...
    HQC(:), ...
    ScoreInnovationVar(:), ...
    ResidualCurveVar(Jgrid(:)), ...
    DesignRank(:), ...
    EffectiveSample(:), ...
    'VariableNames',{'J','m','MSE','BIC','HQC', ...
                     'ScoreInnovationVariance','ResidualCurveVariance', ...
                     'DesignRank','EffectiveSample'});

%% Save output structure

IC.MSE                = MSE;
IC.BIC                = BIC;
IC.HQC                = HQC;
IC.ScoreInnovationVar = ScoreInnovationVar;
IC.ResidualCurveVar   = ResidualCurveVar;
IC.DesignRank         = DesignRank;
IC.EffectiveSample    = EffectiveSample;
IC.GridTable          = GridTable;

IC.BICSelection.J     = J_BIC;
IC.BICSelection.m     = m_BIC;
IC.BICSelection.value = minBIC;
IC.BICSelection.MSE   = MSE(J_BIC,m_BIC);

IC.HQCSelection.J     = J_HQC;
IC.HQCSelection.m     = m_HQC;
IC.HQCSelection.value = minHQC;
IC.HQCSelection.MSE   = MSE(J_HQC,m_HQC);

IC.DynamicDecomposition = DS;
IC.EigenMass = DS.varprop(1:K_max);

IC.T      = T;
IC.K_max  = K_max;
IC.p_max  = p_max;
IC.q_dyn  = q_dyn;

end
