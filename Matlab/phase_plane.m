clear
close all

dx = 0.05;

Beta = 0.6;
Gamma = 0.2;
w = 0.01;

Ii = 1e-5;
Ri = [0 0.1 0.2 0.3 0.4];

tSpan = 0:0.1:2000;

ICs = [1-Ii; Ii] .* (1-Ri);

Sv = dx/2:dx:1;
Iv = dx/2:dx:1;

R0 = Beta/Gamma;


VF = @(Y)([-Beta*Y(1, :, :).*Y(2, :, :); Beta*Y(1, :, :).*Y(2, :, :) - Gamma*Y(2, :, :) ]  ) ;
myODE = @(t, y)(VF(y));

VFw = @(Y)([-Beta*Y(1, :, :).*Y(2, :, :) + w*(1-Y(1, :, :)-Y(2, :, :)); Beta*Y(1, :, :).*Y(2, :, :) - Gamma*Y(2, :, :) ]  ) ;
myODEw = @(t, y)(VFw(y));


[SM, IM] = meshgrid(Sv, Iv);
SM(SM+IM > 1) = nan;
IM(SM+IM > 1) = nan;

Z = VF( [shiftdim(SM, -1); shiftdim(IM, -1)] );
Zmag = sqrt(squeeze(Z(1, :, :)).^2 + squeeze(Z(2, :, :)).^2 );

Zw = VFw( [shiftdim(SM, -1); shiftdim(IM, -1)] );
Zwmag = sqrt(squeeze(Zw(1, :, :)).^2 + squeeze(Zw(2, :, :)).^2 );

nICs = size(ICs, 2);
for iIC = 1:nICs
    [~, Y(:, :, iIC)] = ode45(myODE, tSpan, ICs(:, iIC));
    [~, Yw(:, :, iIC)] = ode45(myODEw, tSpan, ICs(:, iIC));

end

peakPrev = max(Y(:, 2, 1));

greyCol = 0.7*[1 1 1];
redCol = [1 0.9 0.9];
greenCol = [0.9 1 0.9];

colOrd = colororder;



h = figure(1);
 h.Position = [680   429   669   549];
hold on
h1 = fill([1/R0 1 1/R0], [0 0 1-1/R0], redCol, 'DisplayName', 'epidemic growing', 'LineStyle', 'none');
h2 = fill([0 1/R0 1/R0 0], [0 0 1-1/R0, 1], greenCol, 'DisplayName', 'epidemic shrinking', 'LineStyle', 'none');
plot([1 0 0], [0 1 0], 'k-', 'HandleVisibility', 'off')
h3 = quiver(SM, IM, squeeze(Z(1, :, :))./Zmag, squeeze(Z(2, :, :))./Zmag, 0.5, 'DisplayName', 'Vector field' );
h3.Color = greyCol;

for iIC = 1:nICs
    St = Y(:, 1, iIC);
    It = Y(:, 2, iIC);
    pl = plot(St, It, 'Color', 0.5*colOrd(1, :) + 0.5*[1 1 1], 'LineWidth', 2);
    if iIC == 1
        h4 = pl;
        pl.DisplayName = 'Trajectories';
    else
        pl.HandleVisibility = 'off';
    end
    ia = find( St-It < 1/R0, 1, 'first' );
    ax = Y(ia:ia+1, 1, iIC);
    ay = Y(ia:ia+1, 2, iIC);
    an = annotation('arrow' );
    an.Parent = gca;
    an.X = ax;
    an.Y = ay;
    an.Color = 0.5*colOrd(1, :) + 0.5*[1 1 1];


end
h5 = plot(1/R0*[1 1], [0, 1-1/R0], 'r-', 'DisplayName', 'HIT (I nullcline)', 'LineWidth', 2, 'Color', colOrd(2, :));
h6 = plot([0, 1-peakPrev], peakPrev*[1 1], 'k--', 'DisplayName', 'peak I(t) when S(0)=1');
h7 = plot([0, 1/R0], [0 0], 'LineWidth', 2, 'Color', colOrd(3, :), 'DisplayName', 'Neutrally stable equilibria');
h8 = plot([1/R0, 1], [0 0], '--', 'LineWidth', 2, 'Color', colOrd(3, :), 'DisplayName', 'Unstable equilibria');
xlim([0 1])
ylim([-0.005 1])
legend([h3, h4, h5, h6, h7, h8, h1, h2], 'Location', 'northeast')
title(sprintf('R_0=%.1f', R0))
xlabel('susceptible fraction')
ylabel('infectious fraction')

saveas(h, 'SIR_phase_plane.png');



EE = [1/R0;  (1-1/R0)/(1+Gamma/w) ];
Snull = 0:0.01:1;
Inull = w*(1-Snull)./(Beta*Snull+w);

h = figure(2);
 h.Position = [680   429   669   549];
hold on
h1 = fill([1/R0 1 1/R0], [0 0 1-1/R0], redCol, 'DisplayName', 'epidemic growing', 'LineStyle', 'none');
h2 = fill([0 1/R0 1/R0 0], [0 0 1-1/R0, 1], greenCol, 'DisplayName', 'epidemic shrinking', 'LineStyle', 'none');
plot([0 1 0 0], [0 0 1 0], 'k-', 'HandleVisibility', 'off')
h3 = quiver(SM, IM, squeeze(Zw(1, :, :))./Zwmag, squeeze(Zw(2, :, :))./Zwmag, 0.5, 'DisplayName', 'Vector field' );
h3.Color = greyCol;

for iIC = nICs:-1:1
    St = Yw(:, 1, iIC);
    It = Yw(:, 2, iIC);
    if iIC == 1
        h4 = plot(St, It, 'Color', colOrd(1, :), 'LineWidth', 2, 'DisplayName', 'Heteroclinic trajectory');
    elseif iIC == 2
        h5 = plot(St, It, 'Color', 0.5*colOrd(1, :) + 0.5*[1 1 1], 'LineWidth', 2, 'DisplayName', 'Trajectories');
    else
        plot(St, It, 'Color', 0.5*colOrd(1, :) + 0.5*[1 1 1], 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    ia = find( St-It < 1/R0, 1, 'first' );
    ax = Yw(ia:ia+1, 1, iIC);
    ay = Yw(ia:ia+1, 2, iIC);
    an = annotation('arrow' );
    an.Parent = gca;
    an.X = ax;
    an.Y = ay;
    if iIC == 1
        an.Color = colOrd(1, :);
    else
        an.Color = 0.5*colOrd(1, :) + 0.5*[1 1 1];
    end
end
h6 = plot(1/R0*[1 1], [0, 1-1/R0], 'r-', 'DisplayName', 'HIT (I nullcline)', 'LineWidth', 2, 'Color', colOrd(2, :));
h7 = plot(Snull, Inull, 'LineWidth', 2, 'Color', colOrd(4, :), 'DisplayName', 'S nullcline');
h8 = plot(1, 0 , 'ko', 'DisplayName', 'Disease-free equilibrium');
h9 = plot(EE(1), EE(2), 'ko', 'MarkerFaceColor', [0 0 0 ], 'DisplayName', 'Endemic equilibrium');
plot([0 1], [1 0], 'k-', 'HandleVisibility', 'off')
xlim([0 1])
ylim([0 1])
legend([h3, h4, h5, h6, h7, h8, h9, h1, h2], 'Location', 'northeast')
title(sprintf('R_0=%.1f', R0))
xlabel('susceptible fraction')
ylabel('infectious fraction')

saveas(h, 'SIRS_phase_plane.png');


