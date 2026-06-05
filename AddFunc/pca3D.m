function pcastr = pca3D(FTSobj, npca, centerfd)
% PCA3D Functional principal component analysis for surface time series.
%
% Inputs:
%   FTSobj   - fd object or {coefficients, basis}
%   npca     - number of principal components
%   centerfd - 1 if data should be centered, 0 otherwise
%
% Output:
%   pcastr   - structure with loadings, scores, eigenvalues, varprop, meanfd

%% Step 0: defaults and input

if nargin < 3
    centerfd = 0;
end

if isa(FTSobj,'fd')
    fdobj = FTSobj;
elseif iscell(FTSobj) && length(FTSobj)==2
    fdobj = fd(FTSobj{1}, FTSobj{2});
else
    error('Wrong FTS object. Use either fd object or {coefficients,basis}.');
end

%% Step 1: mean and centering

meanfd = mean(fdobj);

if centerfd == 1
    fdobj = center(fdobj);
end

%% Step 2: coefficients and FEM mass matrix

BasisFD  = getbasis(fdobj);
coef_FTS = getcoef(fdobj)';     % T x nbasis
[T, nbasis] = size(coef_FTS);

if npca > nbasis
    error('npca cannot exceed the number of basis functions.');
end

G = FEMMassMatrix(BasisFD);
G = (G + G')/2;

R = chol(G + 1e-10*speye(size(G,1)), 'upper');

% Orthonormalized coordinates
Z = coef_FTS * R';

%% Step 3: covariance in orthonormalized coordinates

CovZ = (Z' * Z) / T;
CovZ = (CovZ + CovZ')/2;

%% Step 4: eigendecomposition

[U,D] = eig(CovZ);

eigvals = diag(D);
[eigvals, inds] = sort(eigvals,'descend');
U = U(:,inds);

eigvals(eigvals < 0 & abs(eigvals) < 1e-10) = 0;

if sum(eigvals) > 0
    varprop = cumsum(eigvals) / sum(eigvals);
else
    varprop = NaN(size(eigvals));
end

%% Step 5: transform loadings back to FEM basis

Theta = R \ U;
Theta = Theta(:,1:npca);

scores = coef_FTS * G * Theta;
harmfd = fd(Theta, BasisFD);

%% Step 6: save output

pcastr.pcafd   = harmfd;
pcastr.values  = eigvals(1:npca);
pcastr.pcascr  = scores;
pcastr.varprop = varprop;
pcastr.meanfd  = meanfd;
pcastr.method  = 'static';

end