clear
close all

R0 = 2;
pInf = 0.1;
nReps = 1e5;
pThresh = 0.2;
N_arr = [3, 4, 5, 6, 8, 10, 12, 14, 20, 30, 40, 50, 75, 100, 125, 150];
I0 = 1;

nN = length(N_arr);

FS = nan(nReps, nN);
for iN = 1:nN
    N = N_arr(iN)
    p = R0/N;
    parfor iRep = 1:nReps
        Ik = I0;
        Sk = N;
        while Ik > 0 & Sk > 0
            Ik = binornd(Sk, 1-(1-p)^Ik);
            Sk = Sk-Ik;
        end
        FS(iRep, iN) = N-Sk;
    end
end



FSD = nan(nReps, nN);
for iN = 1:nN
    N = N_arr(iN)
    parfor iRep = 1:nReps
        Ik = I0;
        Sk = N;
        while Ik > 0 & Sk > 0
            Ik = binornd(Sk, 1-(1-pInf)^Ik);
            Sk = Sk-Ik;
        end
        FSD(iRep, iN) = N-Sk;
    end
end






myF = @(z)(1-z-exp(-R0*z));
zs = fzero(myF, [1-1/R0, 1])
sd = sqrt(zs*(1-zs)./N_arr)/(1-(1-zs)*R0);

h = figure(1);
h.Position = [    112         125        1128         818];
tiledlayout(4, 4, "TileSpacing", "compact");
for iN = 1:nN
    N = N_arr(iN);
    nexttile;
    nHist = histcounts(FS(:, iN), 0:N+1);
    bar((0:N)/N, nHist/nReps);
    xlim([0 1])
    if iN >= 10
        xline(zs, 'r-')
        xline(zs+1.96*sd(iN), 'r--')
        xline(zs-1.96*sd(iN), 'r--')
    end
    title(sprintf('N=%i, p_{20}=%.3f', N, mean(FS(:, iN)/N >= pThresh) ));
    if iN >= 13
        xlabel('final size')
    end
    if mod(iN-1, 4) == 0
        ylabel('probability')
    end
end





h = figure(2);
h.Position = [    112         125        1128         818];
tiledlayout(4, 4, "TileSpacing", "compact");
for iN = 1:nN
    N = N_arr(iN);
    R0i = pInf*N;
    myF = @(z)(1-z-exp(-R0i*z));
    zs = fzero(myF, [1-1/R0i, 1]);
    sd = sqrt(zs*(1-zs)/N)/(1-(1-zs)*R0i);
    nexttile;
    nHist = histcounts(FSD(:, iN), 0:N+1);
    bar((0:N)/N, nHist/nReps);
    if iN >= 10
        xline(zs, 'r-')
        xline(zs+1.96*sd, 'r--')
        xline(zs-1.96*sd, 'r--')
    end
    xlim([0 1])
    title(sprintf('N=%i, p_{20}=%.3f', N, mean(FSD(:, iN)/N >= pThresh) ));
    if iN >= 13
        xlabel('final size')
    end
    if mod(iN-1, 4) == 0
        ylabel('probability')
    end
end



