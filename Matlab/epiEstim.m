function [Rt_est, Rt_CI, Rt_block_est, Rt_block_CI] = epiEstim(nCasesImp, nCasesLoc, GTD, priorShape, priorScale, relInfImp, windowSize, iBreaks, Alpha)

% Code to implement the epiEstim method for estimating the time-varying
% effective reproduction number Reff - see Cori et al for details
%
% USAGE [Rt_est, Rt_est_CI] = epiEstim(nCases, GTD, a, b, windowSize, Alpha)
%
% INPUTS: nCases - a row vector of daily cases (or a matrix whose rows are
% vectors of daily cases to be treated independently)
%         GTD - a vector of PMF values for the generation interval
%         distribution by day (on t = 1, 2, 3, ...)
%         priorShape, priorScale - prior gamma shape and scale parameters
%         for Reff - set shape=1 and scale=inf for a uniform prior
%         windowSize - size of the rolling observation window (in days)
%         iBreaks (optional - set iBreaks = [] for no breaks) - indices of break points to estimate Rt in
%         defined time blocks
%         Alpha - significance level for the confidence interval required
%         (e.g. Alpha=0.05 for a 95% CI)
%
% OUTPUTS: Rt_est - posterior mean estimate for Reff by day
%          Rt_CI - lower and upper bounds of the 95% posterior CrI for Reff
%          Rt_block_est - posterior mean estimate for Reff in specified
%          time blocks
%          Rt_block_CI - lower and upper bounds of the 95% posterior CrI
%          for Reff in specified time blocks

[nSeries, nDays] = size(nCasesLoc);

nCasesTot = relInfImp*nCasesImp + nCasesLoc;        % Effective total cases contributing to force of infections (with imported cases weighted by relInfImp)


% Make sure GTD is a normalised PMF
GTD = GTD/sum(GTD);

c = zeros(nSeries, nDays+length(GTD)-1);
% Analyse each for of nCases one at a time
for iSeries = 1:nSeries
    % Convolution of nCases with generation interval PMF (quantifying
    % contribution to force of infection by day)
    c(iSeries, :) = conv(nCasesTot(iSeries, :), GTD);
end
% Pad with a lead zero and combine each time series of c into a single
% matrix
Gamma = [zeros(nSeries, 1), c(:, 1:nDays-1)];

% Calculate posterior shape and scale parameters using formulae from Cori
% et al
nc = [zeros(nSeries, 1), cumsum(nCasesLoc, 2)];
gc = [zeros(nSeries, 1), cumsum(Gamma, 2)];

% Compute total new cases and aggregated force of infection in a rolling window
nr = nc(:, windowSize+1:end)-nc(:, 1:end-windowSize);
gr = gc(:, windowSize+1:end)-gc(:, 1:end-windowSize);

postShape = priorShape + nr;
postScale = 1 ./ (1/priorScale + gr);

% Calculate mean/median and CI
%Rt_est = [nan(nSeries, windowSize-1), gaminv(1/2, postShape, postScale)]; % posterior median
Rt_est = [nan(nSeries, windowSize-1), postShape.*postScale];              % posterior mean
CI_lower = gaminv(Alpha/2, postShape, postScale);
CI_upper = gaminv(1-Alpha/2, postShape, postScale);
Rt_CI = [nan(2, windowSize-1, nSeries), [shiftdim(CI_lower', -1); shiftdim(CI_upper', -1)]];


% If time blocks are specified, also compute nr and gr in each block
%if ~isempty(iBreaks)                            
    nr = diff(nc(:, [1, iBreaks, end]), [], 2 );
    gr = diff(gc(:, [1, iBreaks, end]), [], 2 );
    postShape = priorShape + nr;
    postScale = 1 ./ (1/priorScale + gr);
    
    % Calculate mean/median and CI
    %Rt_est = [nan(nSeries, windowSize-1), gaminv(1/2, postShape, postScale)]; % posterior median
    Rt_block_est = postShape.*postScale;              % posterior mean
    CI_lower = gaminv(Alpha/2, postShape, postScale);
    CI_upper = gaminv(1-Alpha/2, postShape, postScale);
    Rt_block_CI = [shiftdim(CI_lower', -1); shiftdim(CI_upper', -1)];

%else % return calculation of average over  whole time period
%    Rt_block_est = [];
%    Rt_block_CI = [];
%end

