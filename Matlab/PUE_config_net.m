clear 
close all

R0_arr = 1:0.1:6;
Kappa_arr = [inf, 2, 1, 0.5, 0.1];
P_arr = [0.01, 0.3, 0.7, 0.999];


kMax = 500;
tol = 1e-6;
maxIts = 1000;

n1 = length(R0_arr);
n2 = length(Kappa_arr);
n3 = length(P_arr);

PUE = ones(n1, n2, n3);
PUE_GW = ones(n1, n2);
kk = 0:kMax;

i1_min = find(R0_arr > 1, 1, 'first');

for i3 = 1:n3
    P = P_arr(i3)
    for i2 = 1:n2
        Kappa = Kappa_arr(i2);
        for i1 = i1_min:n1
            R0 = R0_arr(i1);
            % Calculate mean and variance of degree
            if isfinite(Kappa)
                km = R0/P*Kappa/(Kappa+1);
                kvar = km*(km+Kappa)/Kappa;
                negbin_p = Kappa/(Kappa+km);
                [m, v] = nbinstat(Kappa, negbin_p);
                assert(abs(m-km) < tol & abs(v-kvar) < tol );
                pk = nbinpdf(kk, Kappa, negbin_p);
            else
                km = R0/P;
                pk = poisspdf(kk, km);
            end
            pk = pk/sum(pk);
            nuk = kk.*pk;
            nuk = nuk/sum(nuk);
            convFlag = false;
            count = 0;
            piv = 0;
            while ~convFlag && count < maxIts
                piv_sav = piv;
                piv = sum( nuk.*(1-P+P*piv).^(kk-1) );
                convFlag = abs(piv_sav-piv) < tol;
            end
            if ~convFlag
                fprintf('Warning: not converged\n');
            end
            PUE(i1, i2, i3) = sum( pk.*(1-P+P*piv).^kk );
        end
    end
end

% Calculate Galton-Watson PUE
for i2 = 1:n2
    Kappa = Kappa_arr(i2);
    for i1 = i1_min:n1
        R0 = R0_arr(i1);
        if isfinite(Kappa)
            negbin_p = Kappa/(Kappa+R0);
            [m, v] = nbinstat(Kappa, negbin_p);
            pk = nbinpdf(kk, Kappa, negbin_p);
        else
            pk = poisspdf(kk, R0);
        end
        pk = pk/sum(pk);
        convFlag = false;
        count = 0;
        piv = 0;
        while ~convFlag && count < maxIts
            piv_sav = piv;
            piv = sum( pk.*piv.^kk );
            convFlag = abs(piv_sav-piv) < tol;
        end
        if ~convFlag
            fprintf('Warning: GW not converged\n');
        end
        PUE_GW(i1, i2) = piv;
    end
end



h = figure(1);
h.Position = [  45         209        1022         754];
tiledlayout(2, 2, "TileSpacing", "compact");
for i3 = 1:n3
    nexttile;
    plot(R0_arr, PUE(:, :, i3));
    xlabel('R_0')
    ylabel('PUE')
    if i3 == 1
        lgd = legend(string(Kappa_arr));
        lgd.Title.String = '\kappa';
    end
    title(sprintf('P = %.3f', P_arr(i3)));
end


h = figure(2);
h.Position = [   324   639   916   339];
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
plot(R0_arr, PUE(:, :, 2));
xlabel('R_0')
ylabel('PUE')
title('(a) configuration network');
nexttile;
plot(R0_arr, PUE_GW(:, :));
xlabel('R_0')
ylabel('PUE')
title('(b) Galton-Watson branching process');
lgd = legend(string(Kappa_arr));
lgd.Title.String = '\kappa';
lgd.Location = "northeastoutside";

