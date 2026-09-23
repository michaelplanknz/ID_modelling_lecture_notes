function l = logLik_SIR(R0, Gamma, pDet, S_obs, I_obs, inc_obs, i0, N, tMax)

tSpan = 0:1:tMax;


eps = 0.001/N;


% Set IC so that the initial observable prevalence is always i0 regardless of pDet
IC = N*[1-i0/pDet; i0/pDet];

% Define ODE function
mySIR_ODEs = @(t, y, R0, Gamma)( Gamma * [-R0*y(1)*y(2)/N;  R0*y(1)*y(2)/N - y(2) ]  );

% Solve ODE
[t, Y] = ode45(@(t, y)mySIR_ODEs(t, y, R0, Gamma), tSpan, IC);

% Define observables
S = max(eps, Y(:, 1));
I = max(eps, pDet * Y(:, 2));
inc =  max(eps, -pDet * diff(Y(:, 1)));

% Calculate likelihood
l =  nansum(S_obs.*log(S) - S) + nansum(I_obs.*log(I) - I) + nansum(inc_obs.*log(inc) - inc);
