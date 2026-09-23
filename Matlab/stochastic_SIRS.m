clear
close all

rng(53304);

N = 1e4;
Gamma = 0.25;
R0 = 1.5;
Beta = R0 * Gamma;
w = [0.05, 0.005];
I0 = 20;

nReps = 1;
tMax = 1000;
dt = 0.1;

tSav = 0:dt:tMax;

IC = [N-I0; I0];

greyCol = [0.5 0.5 0.5];

h = figure;
h.Position = [  194   338   981   612];
tiledlayout(2, 2, 'TileSpacing', 'compact');
colOrd = colororder;
letters = ["(a)", "(b)"];

nChange = length(w);
for iChange = 1:nChange
    
    myFn = @(t, y, Beta, Gamma, w, N)([-Beta*y(1)*y(2)/N + w*(N-y(1)-y(2)); Beta*y(1)*y(2)/N - Gamma*y(2) ]);
    
    [~, Y] = ode45(@(t, y)myFn(t, y, Beta, Gamma, w(iChange), N), tSav, IC);
    
    Imat = zeros(nReps, length(tSav));
    Smat = N*ones(nReps, length(tSav));
    
    for iRep = 1:nReps
        t = 0;
        I = I0;
        S = N-I0;
        
        Isav = zeros(1, length(tSav));
        Ssav = zeros(1, length(tSav));
    
        while I > 0 & t < tMax
            rInf = Beta*I*S/N;
            rRec = Gamma*I;
            rWane = w(iChange)*(N-S-I);
            tInc = exprnd(1/(rInf+rRec+rWane));
            ind = find(tSav >= t & tSav < t+tInc);
            Isav(ind) = I;
            Ssav(ind) = S;
    
            t = t+tInc;
            u = rand;
            if u < rInf/(rInf+rRec+rWane)
                I = I+1;
                S = S-1;
            elseif u < (rInf+rRec)/(rInf+rRec+rWane)
                I = I-1;
            else
                S = S+1;
            end
        end
        Smat(iRep, :) = Ssav;
        Imat(iRep, :) = Isav;
    end
   
    nexttile(iChange);
    hold on
    plot(tSav, Imat(1, :), 'Color', greyCol)
    plot(tSav, Y(:, 2), 'LineWidth', 2, 'LineStyle', '--', 'Color', colOrd(1, :))
    xlabel('time (days)')
    ylabel('I(t)')
    ylim([0 1000])
    if iChange == 1
        legend('Stochastic realisation', 'SIRS ODE' )
    end
    title(letters(iChange) + sprintf(' w = %.3f days^{-1}', w(iChange)))
    
    nexttile(iChange+2);
    hold on
    plot(Smat(1, :), Imat(1, :), 'Color', greyCol)
    plot(Y(:, 1), Y(:, 2), 'LineWidth', 2, 'LineStyle', '--', 'Color', colOrd(1, :))
    xlim([0 N])
    ylim([0 1000])
    xlabel('S(t)')
    ylabel('I(t)')


end










