clear
close all


N = 1e4;
Gamma = 0.25;
Beta = 1.5 * Gamma;
I0 = 20;

nReps = 1e4;
tMax = 200;
dt = 0.01;

tSav = 0:dt:tMax;

IC = [N-I0; I0];

myFn = @(t, y, Beta, Gamma, N)([-Beta*y(1)*y(2)/N; Beta*y(1)*y(2)/N - Gamma*y(2) ]);

[~, Y] = ode45(@(t, y)myFn(t, y, Beta, Gamma, N), tSav, IC);
pMax_ODE = max(Y(:, 2))/N;
FS_ODE = 1 - Y(end, 1)/N;

FS = zeros(nReps, 1);

Imat = zeros(nReps, length(tSav));
Smat = N*ones(nReps, length(tSav));

parfor iRep = 1:nReps
    t = 0;
    I = I0;
    S = N-I0;
    maxI = I;
    
    Isav = zeros(1, length(tSav));
    Ssav = zeros(1, length(tSav));

    while I > 0
        rInf = Beta*I*S/N;
        rRec = Gamma*I;
        tInc = exprnd(1/(rInf+rRec));
        ind = find(tSav >= t & tSav < t+tInc);
        Isav(ind) = I;
        Ssav(ind) = S;

        t = t+tInc;
        u = rand;
        if u < rInf/(rInf+rRec)
            I = I+1;
            S = S-1;
            maxI = max(maxI, I);
        else
            I = I-1;
        end
    end
    FS(iRep) = 1-S/N;
    Smat(iRep, :) = Ssav;
    Imat(iRep, :) = Isav;
end


greyCol = [0.5 0.5 0.5];

figure;
hold on
plot(tSav, Imat(1:19, :), 'Color', greyCol ,'HandleVisibility', 'off')
plot(tSav, Imat(20, :), 'Color', greyCol)
plot(tSav, mean(Imat), 'LineWidth', 2)
plot(tSav, Y(:, 2), 'LineWidth', 2)
xlabel('time (days)')
ylabel('I(t)')
legend('Stochastic realisations', 'Mean of realisations', 'SIR ODE' )
xlim([0 120])


figure;
histogram(FS);
xline(mean(FS), 'k-')
xline(FS_ODE, 'r-');
xlabel('final size')
ylabel('count')
legend('IBM', sprintf('IBM mean %.3f', mean(FS)), sprintf('ODE %.3f', FS_ODE)   )


