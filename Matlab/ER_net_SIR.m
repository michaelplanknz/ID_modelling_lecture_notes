clear
close all

rng(95026);

% Pop size
nPop = 10000;

% Mean degree, recovery rate and per-edge transmission rate parameters
kMean = 5;
Gamma = 0.2;

% Calculate Beta from specified value of R0
R0 = 5/3;
Beta = Gamma/(kMean/R0-1);

% Initial number infected
I0 = 50;

% Calculate R0 
%R0 = kMean*Beta/(Beta+Gamma);

% Gamma for the SIR MODE model 
% Set GammaSIR = Gamma to compare models with the same R0 and same infectious period
% or  GamamSIR = Gamma+Beta to compare models with the same R0 and same
% epidemic growth rate
GammaSIR = Gamma;      


% Set initial condition and time span for SIR model
IC = [1-I0/nPop ;I0/nPop];
tSpan = 0:100;

% Solve SIR model with equivalent R0
myODE = @(t, y, R0, Gamma)( Gamma*[-R0*y(1)*y(2); R0*y(1)*y(2)-y(2)  ]  );
[tt_fixR0, Y_fixR0] = ode45(@(t, y)myODE(t, y, R0, GammaSIR), tSpan, IC);


% Solve SIR model with equivalent Beta - denoted C1 fpr closure at level of individuals
[tt_C1, Y_C1] = ode45(@(t, y)myODE(t, y, Beta*kMean/Gamma, Gamma), tSpan, IC);

% Initial number of [S, I, SI, SS] for pair approximation:
IC_C2 = [nPop-I0; I0; (nPop-I0)*I0*kMean/nPop; (nPop-I0)^2*kMean/nPop ];

% Solve SIR model pair approximation - denoted C2 for closure at level of pairs
[tt_C2, Y_C2] = ode45(@(t, y)myPairApprox(t, y, Beta, Gamma, kMean), tSpan, IC_C2);


% Initialise node states for network model (0 for S, 1 for I, 2 for R)
state = zeros(nPop, 1);

% Choose a random sample of nodes and make them infecious
seedIDs = randsample(nPop, I0, false);
state(seedIDs) = 1;

% Calculate the network density
pEdge = kMean/(nPop-1);

% Generate the adjacency matrix for an ER network
A = triu(sprand(nPop, nPop, pEdge), 1) > 0;
A = A+A';

% Vector containing the number of neighbours for each node
nNeighs = full(sum(A));

% Maintain a list of nodes in the S and I states
iSus = find(state == 0);
iInf = find(state == 1);

% Maintain a vector of FOIs (only non zero for susceptible nodes)
FOI = Beta * full(sum(A.*(state == 1)))' .* (state == 0);

% Initialise variables for the simulation
finished = false;
evCount = 0;
t = nan(1, 2*nPop+1);
nInf = nan(1, 2*nPop+1);
nSus = nan(1, 2*nPop+1);

t(1) = 0;
nInf(1) = I0;
nSus(1) = nPop-I0;

% Simulate with Gillespie algorithm until there are no infectious nodes
while nInf(evCount+1) > 0

    % Update event counter and time variable
    evCount = evCount + 1;
    t(evCount+1) = t(evCount) + exprnd(1/(sum(FOI)+Gamma*nInf(evCount) ));

    % Decide whether the event was a transmission or a recovery
    u = rand;
    if u < sum(FOI)/(sum(FOI)+Gamma*nInf(evCount))
        % If transmission, choose the newly infected node in proportion to the FOI
        ind = randsample(length(iSus), 1, true, FOI(iSus));
        iEvent = iSus(ind);
        % Update the state of the node the event is happening to and set is
        % FOI to zero
        state(iEvent) = 1;
        FOI(iEvent) = 0;
        % Update the list of susceptible and infecious nodes
        iSus = setdiff(iSus, iEvent);
        iInf = [iInf; iEvent];
        % Get a list of susceptible neighbours
        susNeighs = find(A(:, iEvent) & state == 0);
        % Increase the FOI of each susceptible neighbour by Tau
        FOI(susNeighs) = FOI(susNeighs) + Beta;
        % Record the number of infectious and susceptivle nodes at this event
        nInf(evCount+1) = nInf(evCount)+1;
        nSus(evCount+1) = nSus(evCount)-1;
    else
        % If recovery, choose the recovering node at random from the infectious nodes
        ind = randi(length(iInf));
        iEvent = iInf(ind);
        % Update the state of the node the event is happening to        
        state(iEvent) = 2;
        % Update the list of infecious nodes
        iInf = setdiff(iInf, iEvent);
        % Get a list of susceptible neighbours
        susNeighs = find(A(:, iEvent) & state == 0);
        % Reduce the FOI of each susceptible neighbour by Tau
        FOI(susNeighs) = max(0, FOI(susNeighs) - Beta);
        % Record the number of infectious and susceptivle nodes at this event
        nInf(evCount+1) = nInf(evCount)-1;
        nSus(evCount+1) = nSus(evCount);
    end
end

% Plot results
h = figure(1);
h.Position = [   129         563        1057         418];
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
plot(t, nSus/nPop, tt_C1, Y_C1(:, 1), tt_fixR0, Y_fixR0(:, 1), 'LineWidth', 2)
ylim([0 1])
xlim([0 100])
xlabel('time (days)')
ylabel('fraction susceptible')
legend('network simulation', 'SIR ODE (same beta)', 'SIR ODE (same R_0)');
nexttile;
plot(t, nInf/nPop, tt_C1, Y_C1(:, 2), tt_fixR0, Y_fixR0(:, 2), 'LineWidth', 2)
xlim([0 100])
xlabel('time (days)')
ylabel('fraction infectious')



% Plot pair approximation comparison
h = figure(2);
h.Position = [   129         563        1057         418];
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
plot(t, nSus/nPop, tt_C1, Y_C1(:, 1), 'LineWidth', 2)
hold on
plot(tt_C2, Y_C2(:, 1)/nPop, 'LineWidth', 2)
ylim([0 1])
xlim([0 100])
xlabel('time (days)')
ylabel('fraction susceptible')
legend('network simulation', 'mean-field', 'pair. approx.');
nexttile;
plot(t, nInf/nPop, tt_C1, Y_C1(:, 2), 'LineWidth', 2)
hold on
plot(tt_C2, Y_C2(:, 2)/nPop, 'LineWidth', 2)
xlim([0 100])
xlabel('time (days)')
ylabel('fraction infectious')

