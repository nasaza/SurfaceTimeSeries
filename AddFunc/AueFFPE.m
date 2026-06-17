function Result = AueFFPE(PCA_Str,d_max,p_max)
% AueFFPE
% Implements the functional final prediction error criterion of
% Aue, Dubart Norinho and Hoermann (2015).
%
% For FPCA dimension d and VAR order p:
%
%   fFPE(p,d) = ((n + p*d)/(n - p*d))*trace(SigmaHat_Z)
%               + sum_{ell>d} lambdaHat_ell.
%
% Inputs:
%   PCA_Str - output from pca3D, calculated with centered functional data
%   d_max   - largest candidate FPCA dimension
%   p_max   - largest candidate VAR order
%
% The VARs are estimated by conditional least squares without an intercept,
% because the FPCA scores are centered.

requiredFields = {'pcascr','values','varprop'};
for ii = 1:length(requiredFields)
    if ~isfield(PCA_Str,requiredFields{ii})
        error('PCA_Str is missing field %s.',requiredFields{ii});
    end
end

Scores  = PCA_Str.pcascr;
eigvals = PCA_Str.values(:);
varprop = PCA_Str.varprop(:);
n       = size(Scores,1);

d_max = min([d_max,size(Scores,2),length(eigvals)]);

if d_max < 1
    error('No FPCA components are available.');
end
if p_max < 0 || floor(p_max) ~= p_max
    error('p_max must be a nonnegative integer.');
end

dGrid = 1:d_max;
pGrid = 0:p_max;

dRef = d_max;
if length(varprop) >= dRef && varprop(dRef) > 0
    totalVariance = sum(eigvals(1:dRef))/varprop(dRef);
else
    warning(['Complete eigenvalue mass could not be recovered from varprop. ', ...
             'The available eigenvalue sum is used instead.']);
    totalVariance = sum(eigvals);
end

FFPE       = NaN(length(pGrid),length(dGrid));
TraceSigma = NaN(size(FFPE));
TailMass   = NaN(1,length(dGrid));

for id = 1:length(dGrid)

    d = dGrid(id);
    Y = Scores(:,1:d);

    TailMass(id) = max(totalVariance-sum(eigvals(1:d)),0);

    for ip = 1:length(pGrid)

        p = pGrid(ip);

        if n-p*d <= 0
            continue
        end

        if p == 0
            Resid = Y;
        else
            if n <= p+1
                continue
            end

            Ytarget = Y(p+1:n,:);
            Xlag    = zeros(n-p,p*d);

            for lag = 1:p
                cols = (lag-1)*d + (1:d);
                Xlag(:,cols) = Y(p+1-lag:n-lag,:);
            end

            Bhat  = Xlag\Ytarget;
            Resid = Ytarget-Xlag*Bhat;
        end

        SigmaHat = (Resid'*Resid)/n;
        trSigma  = trace(SigmaHat);

        TraceSigma(ip,id) = trSigma;
        FFPE(ip,id) = ((n+p*d)/(n-p*d))*trSigma + TailMass(id);
    end
end

[minValue,linearIndex] = min(FFPE,[],'all','omitnan');
[idxP,idxD] = ind2sub(size(FFPE),linearIndex);

Result.FFPE          = FFPE;
Result.TraceSigma    = TraceSigma;
Result.TailMass      = TailMass;
Result.dGrid         = dGrid;
Result.pGrid         = pGrid;
Result.selectedD     = dGrid(idxD);
Result.selectedP     = pGrid(idxP);
Result.minimumFFPE   = minValue;
Result.totalVariance = totalVariance;
Result.n             = n;

end
