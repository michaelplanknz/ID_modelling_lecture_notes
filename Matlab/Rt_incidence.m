clear
close all

% Filename to save data
fOut = 'incidence_data.csv';

% Time period to work in
t = 0:80;
nDays = length(t);

% Max, mean and SD for generation interval
GTmax = 15;
GTmean = 5;
GTsd = 2.5;

% Rt values before and after the ramp, and start/end times for the ramp
R0 = 1.5;
R1 = 0.7;
ramp1 = 30;
ramp2 = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crete a vector of daily Rt values
Rtrue = nan(1, nDays);
Rtrue(t < ramp1) = R0;
Rtrue(t >= ramp1 & t < ramp2) = linspace(R0, R1, ramp2-ramp1);
Rtrue(t >= ramp2) = R1;

% Create a generation interval PMF vector using a discretised gamma
% distribution on [1, 2, ...] with specified mean and SD
[sh, sc] = gamShapeScale(GTmean, GTsd);
pdfFunc = @(x)gampdf(x, sh, sc);
gs = discDist(pdfFunc, 1, GTmax);



% Setup an array for incidence data
inc = nan(1, nDays);

% Initialise incidence data in a burn in period up to t=GTmax using
% exponential growth
r = (R0-1)/GTmean;
X = exp(r*t(1:GTmax+1));
inc = poissrnd(X);

% Simulate daily incidence data using the Poisson renewal model
for iDay = GTmax+2:nDays
    X = Rtrue(iDay) * sum(gs.*inc(iDay-1:-1:iDay-GTmax));
    inc(iDay) = poissrnd(X);
end


% Save data as a table
tbl.t = t';
tbl.Rtrue = Rtrue';
tbl.inc = inc';
tbl = struct2table(tbl);

writetable(tbl, fOut);
