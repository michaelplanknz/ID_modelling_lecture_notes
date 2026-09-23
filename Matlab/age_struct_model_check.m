clear
close all

% Settings for NGM
act_lvl = [1.5, 1.25, 1, 0.75, 0.5];
eps = 0.2;

% ODE model settings to check equilibrium results
Gamma = 0.2;
w = 0.01;
tSpan = [0, 10000];
i0 = 1e-6;

% NGM
K = 5 * (eps * (act_lvl'*act_lvl)/sum(act_lvl)  + (1-eps) * diag(act_lvl));


% Values of R0 to run
R0_vec = [1.1, 3];

% Number of age groups
[na, ~] = size(K);

% Number of iterations to do
nIts = 1000;

% Index of dominant eigenvalue (largest absolute value)
[v_dom, lambda_dom] = eigs(K, 1);
v_dom = v_dom'/sum(v_dom)

% Population size vector (na equally sized groups)
N = 1/na * ones(na, 1);


% Define model ODEs and IC for checking equilibrium results
IC = [(1-i0)*N; i0*N];
mySIR = @(t, y, K, Gamma, N)( Gamma * [-y(1:na)./N.*(K*y(na+1:end));  y(1:na)./N.*(K*y(na+1:end)) - y(na+1:end)  ]  );
mySIRS = @(t, y, K, Gamma, N)( Gamma * [-y(1:na)./N.*(K*y(na+1:end)) + w/Gamma*(N-y(1:na)-y(na+1:end));  y(1:na)./N.*(K*y(na+1:end)) - y(na+1:end)  ]  );




for jj = 1:length(R0_vec)

  R0 = R0_vec(jj);
  
  Kn = R0/lambda_dom * K;

  % Initial guess for final size vector z and equilibrium susceptible pop S
  z = (1-1/R0) * ones(na, 1);
  NS = 1/na * (1-1/R0) * ones(na, 1);

  % Use 1/R0 as a damping factor for fixed pt iteration
  Omega = 1/R0;

  % Iterate final size and endemic equilibrium equations
  for ii = 1:nIts
    z = 1 - exp(-Kn' * z);
    NSsav = NS;
    M = (1-NS./N).*Kn;
    NS = (1-Omega)*NS + Omega * M*NS;
  end

  [t, Y] = ode45(@(t, y)mySIR(t, y, Kn, Gamma, N), tSpan, IC  );
  zODE = 1 - (Y(end, 1:na) + Y(end, na+1:end))'./N;

  [t, Y] = ode45(@(t, y)mySIRS(t, y, Kn, Gamma, N), tSpan, IC  );
  nsODE = 1 - Y(end, 1:na)'./N;


  R0
  z = z'
  zODE = zODE'
  ns = (NS./N)'
  nsODE = nsODE'
  ratio = z./ns
end
