clear
close all

% File name where incidence data is stored
fIn = 'incidence_data.csv';


% Max, mean and SD for generation interval
GTmax = 15;
GTmean = 5;
GTsd = 2.5;


% EpiEstim settings
priorShape = 1;
priorScale = 2;
windowSize = 7;
Alpha = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read incidence data from file
tbl = readtable(fIn);
t = tbl.t';
Rtrue = tbl.Rtrue';
inc = tbl.inc';

nDays = length(t);


% Create a generation interval PMF vector using a discretised gamma
% distribution on [1, 2, ...] with specified mean and SD
[sh, sc] = gamShapeScale(GTmean, GTsd);
pdfFunc = @(x)gampdf(x, sh, sc);
gs = discDist(pdfFunc, 1, GTmax);



% Calculate naive Rt estimate by inverting the renewal equation
inf_press = conv(inc, [0, gs]);
inf_press = inf_press(1:nDays);
Rt_naive =  inc./inf_press;

% Calculate Rt using epiEstim
[Rt_est, Rt_CI] = epiEstim(0*inc, inc, gs, priorShape, priorScale, 1, windowSize, [], Alpha);


% Calculate Rt using epiEstim with a mis-specified GI
[sh, sc] = gamShapeScale(2*GTmean, GTsd);
pdfFunc = @(x)gampdf(x, sh, sc);
gs_mis = discDist(pdfFunc, 1, GTmax);
[Rt_est_mis, Rt_CI_mis] = epiEstim(0*inc, inc, gs_mis, priorShape, priorScale, 1, windowSize, [], Alpha);



% Plot results
lightRed = [1 0.6 0.6];
colOrd = colororder;

h = figure;
h.Position = [   93    52   748   943];
tiledlayout(3, 1, "TileSpacing", "compact")

nexttile;
plot(t, inc, '.')
xlabel('time (days)')
ylabel('daily incidence')
title('(a)')

nexttile;
plot(t, Rtrue, 'LineWidth', 2)
hold on
plot(t, Rt_naive, '.')
xlabel('time (days)')
ylabel('R_t')
ylim([0 2.5])
legend('true value', 'estimated value')
title('(b)')

tClip = t(windowSize:end);
Rt_CI = Rt_CI(:, windowSize:end);
nexttile;
h1 = plot(t, Rtrue, 'LineWidth', 2);
hold on
h2 = fill([tClip, fliplr(tClip)]-windowSize/2, [Rt_CI(1, :), fliplr(Rt_CI(2, :))], lightRed, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h3 = plot(t-windowSize/2, Rt_est, 'color', colOrd(2, :), 'LineWidth', 2);
ylim([0 2.5])
xlabel('time (days)')
ylabel('R_t')
ha = gca;
ha.Children = [h2; h3; h1];
legend('true value', 'estimated value', '95% CI')
title('(c)')
saveas(h, 'Rt_inference.png')

Rt_CI_mis = Rt_CI_mis(:, windowSize:end);
h = figure;
h.Position = [ 205   708   776   354];
h1 = plot(t, Rtrue, 'LineWidth', 2);
hold on
h2 = fill([tClip, fliplr(tClip)]-windowSize/2, [Rt_CI_mis(1, :), fliplr(Rt_CI_mis(2, :))], lightRed, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
h3 = plot(t-windowSize/2, Rt_est_mis, 'color', colOrd(2, :), 'LineWidth', 2);
xlabel('time (days)')
ylabel('R_t')
ylim([0 5])
ha = gca;
ha.Children = [h2; h3; h1];
legend('true value', 'estimated value', '95% CI')
saveas(h, 'Rt_inference2.png')

