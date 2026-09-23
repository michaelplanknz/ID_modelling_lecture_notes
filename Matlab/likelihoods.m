clear
close all

rng(55012);

% True parameter values
R0 = 2;
Gamma = 0.1;
pDet = 0.4;

% Define grid of parameters to use
R0_arr = linspace(1.9, 2.1, 100);
Gamma_arr = linspace(0.09, 0.11, 100);

% Initial infected fraction and pop size
i0 = 0.001;
N = 1e4;

% Time to run model for (days)
tMax = 200;

% Significance level for confidence sets
Alpha = 0.05;


% Define ODE function, time span and IC vector
mySIR_ODEs = @(t, y, R0, Gamma, N)( Gamma * [-R0*y(1)*y(2)/N;  R0*y(1)*y(2)/N - y(2) ]  );
tSpan = 0:1:tMax;
IC = N*[1-i0; i0];


% Likelihood ratio threshold 
ls = chi2inv(1-Alpha, 2)/2;

% Array sizes
nRows = length(R0_arr);
nCols = length(Gamma_arr);

h = figure(1);
h.Position = [     68   277   934   685];
tiledlayout(2, 2, "TileSpacing", "compact");


%%
% S and I data

% Generate synthetic data
[t, Y] = ode45(@(t, y)mySIR_ODEs(t, y, R0, Gamma, N), tSpan, IC);
S_obs = poissrnd( max(0, Y(:, 1)));
I_obs = poissrnd( max(0,  Y(:, 2)));

% Set up matrix to store log likelihoods
LL = zeros(nRows, nCols);

% Cycle through parameters and calculate log likelihood for each combination
for iRow = 1:nRows
    for jCol =1:nCols
        LL(iRow, jCol) = logLik_SIR(R0_arr(iRow), Gamma_arr(jCol), 1, S_obs, I_obs, nan(length(t)-1, 1), i0, N, tMax);
    end
end

% Indices of max likelihood
[im, jm] = find(LL == max(max(LL)));

% Normalised log likelihood (with max 0)
LLn = LL - max(max(LL));


% Plotting
nexttile;
imagesc( R0_arr, Gamma_arr, LLn');
h = gca; h.YDir = "normal";
clim([-ls 0])
hold on
plot(R0, Gamma, 'ko')
plot(R0_arr(im), Gamma_arr(jm) , 'ro')
xlabel('R_0')
ylabel('\gamma (day^{-1})')
title('(a)')


%%
% Prevalence data

% Generate synthetic data
[t, Y] = ode45(@(t, y)mySIR_ODEs(t, y, R0, Gamma, N), tSpan, IC);
I_obs = poissrnd( max(0, Y(:, 2)));


% Set up matrix to store log likelihoods
LL = zeros(nRows, nCols);

% Cycle through parameters and calculate log likelihood for each combination
for iRow = 1:nRows
    for jCol =1:nCols
        LL(iRow, jCol) = logLik_SIR(R0_arr(iRow), Gamma_arr(jCol), 1, nan(length(t), 1), I_obs, nan(length(t)-1, 1), i0, N, tMax);
    end
end

% Indices of max likelihood
[im, jm] = find(LL == max(max(LL)));

% Normalised log likelihood (with max 0)
LLn = LL - max(max(LL));


% Plotting
nexttile;
imagesc( R0_arr, Gamma_arr, LLn');
h = gca; h.YDir = "normal";
clim([-ls 0])
hold on
plot(R0, Gamma, 'ko')
plot(R0_arr(im), Gamma_arr(jm) , 'ro')
xlabel('R_0')
ylabel('\gamma (day^{-1})')
title('(b)')







%%
% Incidence data

% Generate synthetic data
[t, Y] = ode45(@(t, y)mySIR_ODEs(t, y, R0, Gamma, N), tSpan, IC);
inc_daily = max(0, -diff(Y(:, 1)));
inc_obs = poissrnd(inc_daily);

% Set up matrix to store log likelihoods
LL = zeros(nRows, nCols);

% Cycle through parameters and calculate log likelihood for each combination
for iRow = 1:nRows
    for jCol =1:nCols
        LL(iRow, jCol) = logLik_SIR(R0_arr(iRow), Gamma_arr(jCol), 1, nan(length(t), 1), nan(length(t), 1), inc_obs, i0, N, tMax);
    end
end

% Indices of max likelihood
[im, jm] = find(LL == max(max(LL)));

% Normalised log likelihood (with max 0)
LLn = LL - max(max(LL));

% Plotting
nexttile;
imagesc( R0_arr, Gamma_arr, LLn');
h = gca; h.YDir = "normal";
hold on
clim([-ls 0])
plot(R0, Gamma, 'ko')
plot(R0_arr(im), Gamma_arr(jm) , 'ro')
xlabel('R_0')
ylabel('\gamma (day^{-1})')
title('(c)')


%%
% Prevalence data with underascertainment

R0_arr = linspace(1.6, 2.8, 100);
Gamma_arr = linspace(0.06, 0.14, 100);

% Set IC so that the initial observable prevalence is always i0 regardless of pDet
IC = N*[1-i0/pDet; i0/pDet];

% Generate synthetic data
[t, Y] = ode45(@(t, y)mySIR_ODEs(t, y, R0, Gamma, N), tSpan, IC);
inc_daily = max(0, -diff(Y(:, 1)));
%I_obs = poissrnd( max(0, Y(:, 2)));
inc_obs = poissrnd( pDet * inc_daily);

% Set up matrix to store log likelihoods
LL = zeros(nRows, nCols);
pDetEst = zeros(nRows, nCols);

% Cycle through parameters and calculate log likelihood for each combination
p0 = 0.5;
opts = optimset('Display','off');
for iRow = 1:nRows
    for jCol =1:nCols
        objFn = @(x)(-logLik_SIR(R0_arr(iRow), Gamma_arr(jCol), x, nan(length(t), 1), nan(length(t), 1), inc_obs, i0, N, tMax));
        [x, f] = fmincon(objFn, p0, [], [], [], [], 0, 1, [], opts );
        LL(iRow, jCol) = -f;
        pDetEst(iRow, jCol) = x;
    end
end


% Indices of max likelihood
[im, jm] = find(LL == max(max(LL)));

% Normalised log likelihood (with max 0)
LLn = LL - max(max(LL));


% Plotting
nexttile;
imagesc( R0_arr, Gamma_arr, LLn');
colorbar;
h = gca; h.YDir = "normal";
clim([-ls 0])
hold on
plot(R0, Gamma, 'ko')
plot(R0_arr(im), Gamma_arr(jm) , 'ro')
plot([1.9 2.1 2.1 1.9 1.9], [0.09 0.09 0.11 0.11 0.09], '-', 'Color', [1 1 1] )
xlabel('R_0')
ylabel('\gamma (day^{-1})')
title('(d)')


Z = pDetEst;
Z(LLn < -ls) = nan;

figure;
imagesc( R0_arr, Gamma_arr, Z');
h = gca; h.YDir = "normal";
colorbar
hold on
plot(R0, Gamma, 'ko')
plot(R0_arr(im), Gamma_arr(jm) , 'ro')
xlabel('R_0')
ylabel('\gamma (day^{-1})')


