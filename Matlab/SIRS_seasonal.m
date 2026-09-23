clear 
close all

par.tI = 5;                             % Infectious period (days)
par.eps = 1e-7;                         % External infectoin rate (proportion of population infected per day)
change1 = [0 0.1 0.2 0.3];              % seasonality amplitude
change2 = [100 200 300 1000];           % waning time (days)
change3 = [2 4 6 8 10];                 % R0

% Initial infectious fraction
I0 = 0.0001;

% Time span inlcuding a burin period (t<0) and period to be plotted (t>0)
tSpan = -365.25*80:1:365.25*4; 

nChange1 = length(change1);
nChange2 = length(change2);
nChange3 = length(change3);

nSteps = sum(tSpan >= 0);
tYears = tSpan(tSpan >= 0)/365.25;

% IC
y0 = [1-I0; I0];

% ODE solver options
opts = odeset('NonNegative', [1; 1]);



for iChange1 = 1:nChange1
    par.seasonAmp = change1(iChange1);
    for iChange2 = 1:nChange2
        par.tw = change2(iChange2);
        Isol = zeros(nChange3, nSteps);
        Ssol = zeros(nChange3, nSteps);
        Rt = zeros(nChange3, nSteps);
        for iChange3 = 1:nChange3
            par.R0 = change3(iChange3);
            [t, Y] = ode45(@(t, y)mySIRseas(t, y, par), tSpan, y0, opts);
            Isol(iChange3, :) = Y(t >= 0, 2)';
            Ssol(iChange3, :) = Y(t >= 0, 1)';
            seas = 1 - par.seasonAmp*cos(2*pi*t(t >= 0)/365.25)';
            Rt(iChange3, :) = par.R0*seas.*Ssol(iChange3, :);
        end

        figure(1);
        subtightplot(nChange1, nChange2, nChange2*(iChange1-1)+iChange2)
        plot(tYears, Isol, 'LineWidth', 2)
        if iChange1 > 1 || iChange2 == nChange2
            set(gca, 'ColorOrderIndex', 1);
            Iavg = trapz(tYears, Isol')'./trapz(tYears, ones(size(tYears)));
            hold on
            %plot(tYears, Iavg.*ones(size(tYears)), ':')
        end
        ylim([0 0.08])
        grid on
        if iChange1 == 1
           title(sprintf('immunity = %i days', par.tw));
           if iChange2 == nChange2
              lg = legend( string(change3) );
              title(lg, 'R0');
           end
        end
        if iChange1 == nChange1 
            xlabel('time (years)')
        end
        if iChange1 < nChange1        
            set(gca, 'XTickLabel', '');
        end
        if iChange2 == 1
           ylabel(sprintf('infectious (seasonality = %i%%)', 100*par.seasonAmp))
        else
            set(gca, 'YTickLabel', '');
        end

        
        figure(2);
        subtightplot(nChange1, nChange2, nChange2*(iChange1-1)+iChange2)
        plot(tYears, Rt)
        ylim([0.8 1.3])
        grid on
        if iChange1 == 1
           title(sprintf('immunity = %i days', par.tw));
           if iChange2 == nChange2
              lg = legend( string(change3) );
              title(lg, 'R0');
           end
        end
        if iChange1 == nChange1 
            xlabel('time (years)')
        end
        if iChange1 < nChange1        
            set(gca, 'XTickLabel', '');
        end
        if iChange2 == 1
           ylabel(sprintf('Rt (seasonality = %i%%)', 100*par.seasonAmp))
        else
            set(gca, 'YTickLabel', '');
        end
   
        
        
        
        
        
        figure(3);
        subtightplot(nChange1, nChange2, nChange2*(iChange1-1)+iChange2)
        if iChange1 > 1 || iChange2 == nChange2 
            plot(Ssol', Isol')
        else
             plot(Ssol', Isol', '.')
        end
        xlim([0 0.8])
        ylim([0 0.1])
        grid on
        if iChange1 == 1
           title(sprintf('immunity = %i days', par.tw));
           if iChange2 == nChange2
              lg = legend( string(change3) );
              title(lg, 'R0');
           end
        end
        if iChange1 == nChange1 
            xlabel('susceptible')
        end
        if iChange1 < nChange1        
            set(gca, 'XTickLabel', '');
        end
        if iChange2 == 1
           ylabel(sprintf('infectious (seasonality = %i%%)', 100*par.seasonAmp))
        else
            set(gca, 'YTickLabel', '');
        end
       
        
        
    end
end

