clear
close all

% Set parameter values and ranges
Gamma = 0.2;
cMean = 0.6;
cr_arr = 1:0.1:5;
eps_arr = 0:0.02:1;
N1_arr = [0.5 0.8 0.95];

maxIter = 500;
tol = 1e-7;


% Calculate R0 and final size for homogeneous model with same avg contact rate 
Rh = cMean/Gamma
myf = @(x)(1-x-exp(-Rh*x));
FSh = fzero(myf, [1-1/Rh ,1])

nChange1 = length(N1_arr);
nChange2 = length(cr_arr);
nChange3 = length(eps_arr);

% Initial infected fraction
i0 = 1e-4;
tSpan = [0, 500];

% Initialise arrays
[R0, FS, FS2] = deal(zeros(nChange3, nChange2, nChange1));

for iChange1 = 1:nChange1
    % Pop size vector
    N = [N1_arr(iChange1), 1-N1_arr(iChange1)];

    % Initial condition
    IC = N'*[1-i0, i0];
    IC = IC(:);

    for iChange2 = 1:nChange2
        % Calculate contact rates so their ratio is cr and their mean is cMean
        c1 = cMean/(N(1) + (1-N(1))*cr_arr(iChange2));
        c = c1*[1, cr_arr(iChange2)];
        for iChange3 = 1:nChange3
            eps = eps_arr(iChange3);
            % Contact matrix: Cij is contacts a type j individual has with type i
            C =  (1-eps)*N'.*c'.*c/cMean + eps*diag(c) ;
            % Calculate R0
            R0(iChange3, iChange2, iChange1) = eigs(C/Gamma, 1);
            % Solve ODEs
            [t, Y] = ode45(@(t, y)myTwoTypeODE(t, y, C, N', Gamma), tSpan, IC);
            FS(iChange3, iChange2, iChange1) = 1-sum(Y(end, :));

            % Check final size by Andreasen method
            K = C/Gamma;
            A = K.*N./N';
            convFlag = false;
            iIter = 0;
            z = [0.5; 0.5];
            while ~convFlag & iIter < maxIter
                zSav = z;
                z = 1-exp(-A*z);
                convFlag = norm(z-zSav, inf) < tol;
            end
            FS2(iChange3, iChange2, iChange1) = sum(N'.*z);
        end
    end
end


h = figure(1);
h.Position = [     88         367        1166         596];
tiledlayout(2, nChange1, "TileSpacing", "compact");
for iChange1 = 1:nChange1
    nexttile(iChange1);
    imagesc(cr_arr, eps_arr, R0(:, :, iChange1)/Rh); 
    colormap(hot);
    h = gca; h.YDir = 'normal';
    clim([1 4.5])
    hc = colorbar;
    if iChange1 == 1
        ylabel('\epsilon')
    elseif iChange1 == nChange1
        hc.Label.String = "R_0/R_H";
    end
    title(sprintf('N_2 = %.2f', 1-N1_arr(iChange1)))
    nexttile(iChange1+nChange1);
    imagesc(cr_arr, eps_arr, FS1(:, :, iChange1))
    colormap(hot);
    h = gca; h.YDir = 'normal';
    hc = colorbar;
    if iChange1 == 1
        ylabel('\epsilon')
    elseif iChange1 == nChange1
       hc.Label.String = "final size";
    end
    clim([0.5 0.95])
    xlabel('c_2/c_1')

end


