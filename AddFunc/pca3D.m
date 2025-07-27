function pcastr = pca3D(FTSobj, npca, centerfd)
%  PCA Functional principal components analysis for surface (functional)
%  time series
% Inputs:
%  - fdobj          functional time series;
%  - npca           number of pca components/factors to be computed;
%  - centerfd       centered or not centered data; 0 for not centered and 1
%                   for centered

%% Step 0: Checking wich variable are passe to the function
n_inputs = length(FTSobj);
if n_inputs==1
    fdobj=FTSobj;
elseif n_inputs==2
    fdobj=fd(FTSobj{1},FTSobj{2});
% elseif n_inputs==3
%     fdobj=FTSobj{1};    
else
    disp('Wring FTS object is used as input. Use either only FTS or FTS and corresponding FTS coeff with basis to speed up code');
    return
end

%% Step 1:
meanfd = mean(fdobj); % calculate the mean

if nargin < 3
    centerfn = 0;   %  subtract mean from data before PCA
end

if centerfd == 1     
    fdobj  = center(fdobj);
end

%% Step 2: get coefficients of FTS in basis representation
BasisFD    = getbasis(fdobj);
coef_FTS   = getcoef(fdobj)';
[T,~]      = size(coef_FTS);

%% Step 3: Calculate estimate of the Cov of fdobj

CovX       = 1/T*(coef_FTS'*coef_FTS);

%% Step 4: calculate  nharm eigenvectors of CovX

[Theta,D]  = eig(CovX); % Theta is a matrix that contains coeff of eigencalues in BasisFD
                        % D if diagonal matrix with eigenvalues
eigvals    = diag(D);
[~,inds]   = sort(eigvals);
eigvals    = eigvals(flip(inds));
Theta      = Theta(:,flip(inds));

%% Step 5: Saving outputs 
scores     = coef_FTS*Theta;
varprop    = 1/sum(eigvals)*cumsum(eigvals);
harmfd     = fd(Theta, BasisFD);


pcastr.pcafd   = harmfd;
pcastr.values  = eigvals(1:npca);
pcastr.pcascr  = scores(:,1:npca);
pcastr.varprop = varprop;
pcastr.meanfd  = meanfd;


