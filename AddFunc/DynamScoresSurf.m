function DS = DynamScoresSurf(FTSobj, npca, centerfd,q)
%  Extract dynamic components from surface time series based on  Bathia et al., 2010
%  and Otto and Salish, 2026.
%  That is the dynamic components are based on lagged autocovariance operators:
%
%     K = sum_{h=1}^q C_h C_h^*
%
% where C_h is the lag-h autocovariance operator.
%
%
% Inputs:
%  - fdobj          functional time series;
%  - npca           number of pca components/factors to be computed;
%  - centerfd       centered or not centered data; 0 for not centered and 1
%                   for centered
%   - q             maximum lag used in the dynamic operator
% Output:
%   - DS            structure containing dynamic directions, scores, eigenvalues,
%                   cumulative dynamic variation, and mean function


%% Step 0: Checking initial values
if nargin < 3
    centerfd = 0;
end

if nargin < 4
    q = 1;
end

if isa(FTSobj,'fd')
    fdobj = FTSobj;
elseif iscell(FTSobj) && length(FTSobj)==2
    fdobj = fd(FTSobj{1},FTSobj{2});
else
    error('Wrong FTS object. Use either fd object or {coefficients,basis}.');
end

%% Step 1: mean and centering

meanfd = mean(fdobj); % calculate the mean

if centerfd == 1     
    fdobj  = center(fdobj);
end


%% Step 2: get coefficients of FTS in basis representation

BasisFD     = getbasis(fdobj); 
coef_FTS    = getcoef(fdobj)'; % getcoef(fdobj) is usually [nbasis x T]
[T, ~]      = size(coef_FTS);

if q >= T
    error('q must be smaller than T.');
end

if npca > size(coef_FTS,2)
    error('npca cannot exceed the number of basis functions.');
end


G   = FEMMassMatrix(BasisFD);
G   = (G + G')/2;
R   = chol(G + 1e-10*speye(size(G,1)), 'upper');
Z   = coef_FTS * R';


%% Step 3: Calculate estimate of the Cumulative AutoCov Operator

DynCov = zeros(size(Z,2));

for h = 1:q
    Gamma_h = (Z(1+h:T,:)' * Z(1:T-h,:)) / (T-h);
    DynCov  = DynCov + Gamma_h * Gamma_h';
end

DynCov = (DynCov + DynCov')/2;

%% Step 6: eigen-decomposition of cumulative Autocovariance

[U,D] = eig(DynCov);

eigvals = diag(D);
[~,inds] = sort(eigvals);
eigvals = eigvals(flip(inds));

if sum(eigvals) > 0
    varprop = cumsum(eigvals) / sum(eigvals);
else
    varprop = NaN(size(eigvals));
end

U       = U(:,flip(inds));
Theta   = R \ U;
Theta   = Theta(:,1:npca);
scores  = coef_FTS * G * Theta;
harmfd  = fd(Theta, BasisFD);

%% Step 5: Saving outputs 

DS.pcafd   = harmfd;
DS.values  = eigvals(1:npca);
DS.scr     = scores(:,1:npca);
DS.varprop = varprop;
DS.meanfd  = meanfd;


end