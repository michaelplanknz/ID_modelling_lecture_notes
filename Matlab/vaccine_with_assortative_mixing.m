clear
close all

% Array of assortativity parameters to use
eps_arr = 0:0.01:1;

% Array of vaccine effectiveness to use
VE_arr = 0:0.1:1;

% Proportion vaccinated
pv = 0.9;

% R0
R0 = 4;

% Initially infected fraction
I0 = 1e-4;

% Time span, subpopulation sizes, and initial condition for solving SIR model
tSpan = [0 100];
N = [1-pv; pv];
IC = [(1-pv)*(1-I0); pv*(1-I0); (1-pv)*I0; pv*I0];

n1 = length(VE_arr);
n2 = length(eps_arr);

% Initialise arrays for the vaccinated reproduction number and final epidemic size
Rv = zeros(n1, n2);
FS = zeros(n1, n2);

for i1 = 1:n1
    for i2 = 1:n2
        eps = eps_arr(i2);
        VE = VE_arr(i1);
        
        % Contact matrix
        C =  [ (1-eps)*(1-pv) + eps,  (1-eps)*(1-pv) ;
                (1-eps)*pv         ,  (1-eps)*pv + eps];

        % Next generation matrix
        NGM = R0* C.*[1; 1-VE];

        % Dominant eigenvalue of NGM
        Rv(i1, i2) = eigs(NGM, 1);
        
        % Solve SIR ODE to get final size
        [t, Y] = ode45(@(t, y)myVaccineODE(t, y, R0, C, N, VE), tSpan, IC);
        Rinf = N - Y(end, 1:2)';
        FS(i1, i2) = sum(Rinf);
    end
end

colOrd = flipud(hot(n1+2));

% Plot results
h = figure;
h.Position = [ 680   631   880   347];
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
colororder(colOrd(3:end, :));
plot(eps_arr, Rv, 'LineWidth', 2)
ylim([0 R0])
xlabel('assortativity (\epsilon)')
ylabel('R_v')
title('(a)')
lgd = legend(string(VE_arr), 'location', 'southeast');
lgd.Title.String = "v_I";
nexttile;
plot(eps_arr, FS, 'LineWidth', 2)
ylim([0 1])
xlabel('assortativity (\epsilon)')
ylabel('proportion infected')
title('(b)')

