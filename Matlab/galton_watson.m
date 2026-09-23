clear 
close all

R0 = [0.9, 2];
nSteps = 6;
nToLabel = 3;
TOL = 1e-10;

R0_arr = 0:0.05:5;
k_vals = [inf, 2, 1, 0.5, 0.1];

s = 0:0.01:1;
PUE = min(1, 1./R0);


h = figure;
h.Position = [   680   608   849   370];
tiledlayout(1, 2,"TileSpacing", "compact");
letters = ["(a)", "(b)"];
colOrd = colororder;

nCases = length(R0);
for iCase = 1:nCases
    PGF = @(ss)( 1./(1+R0(iCase)*(1-ss)) );
    G = PGF(s);
    xs = zeros(1, 2*nSteps-1);
    ys = zeros(1, 2*nSteps-1);
    ys(1) = PGF(0);
    for iStep = 1:nSteps-1
        xs(2*iStep) = ys(2*iStep-1);
        ys(2*iStep) = ys(2*iStep-1);
        xs(2*iStep+1) = xs(2*iStep);
        ys(2*iStep+1) = PGF(xs(2*iStep+1));
    end

    nexttile;
    hold on
    plot(s, s, 'k-')
    plot(s, G, 'Color', colOrd(1, :), 'LineWidth', 2)
    plot(xs, ys, 'color', colOrd(2, :))
    for iLbl = 1:nToLabel
        plot(xs(2*iLbl-1)*[1 1], [0 ys(2*iLbl-1)], '--', 'Color', colOrd(2, :))
    end
    plot(PUE(iCase)*[1 1], [0 PUE(iCase)],  '--', 'Color', colOrd(5, :))
    xlim([-0.005 1])
    ylim([0 1])
    xlabel('s')
    ylabel('G_X(s)')
    title(letters(iCase) + sprintf(' R_0 = %.1f', R0(iCase)))
end
saveas(h, "staircase.png");


nk = length(k_vals);
nr = length(R0_arr);
PUE = nan(nk, nr);
for ik = 1:nk
    for ir = 1:nr
        if isinf(k_vals(ik))
           PGF = @(ss)( exp(-R0_arr(ir)*(1-ss)) );
        else
           PGF = @(ss)( (k_vals(ik)./(k_vals(ik) + R0_arr(ir)*(1-ss) )).^k_vals(ik)  );
        end
        q = 0;
        convFlag = false;
        while convFlag == false
            qSav = q;
            q = PGF(q);
            convFlag = abs(qSav-q) < TOL;
        end
        PUE(ik, ir) = q;
    end
end


h = figure;
h.Position = [ 680   647   644   331];
plot(R0_arr, PUE, 'LineWidth', 2)
xlabel('R_0')
ylabel('PUE')
lgd = legend(string(k_vals), 'Location', 'northwest');
lgd.Title.String = '\kappa';
grid on

saveas(h, "PUE.png");
