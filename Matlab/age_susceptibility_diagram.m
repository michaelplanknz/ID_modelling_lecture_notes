clear
close all

L = 70;
A = 10;
p = 0.6;
Lambda = 1/A;

a = 0:0.01:L;
S = exp(-Lambda*a);

Av = A/(1-p);

colOrd = colororder;
blue = [0.8 0.8 1];
grey = [0.8, 0.8 0.8];

h = figure;
tl = tiledlayout(1, 1);
ax2 = axes(tl);
plot(ax2, a, ones(size(a)))
ax2.XTickLabel = {'', 'A', '', '', '', '', '', 'L'};
ax2.XAxisLocation = 'top';
ax2.YLim = [0 1];
ax1 = axes(tl);
hold on
fill(ax1, [0, A, A, 0], [0, 0, 1, 1], blue, 'EdgeColor', 'none')
fill(ax1, [A, L, L, A], [0, 0, 1, 1], grey, 'EdgeColor', 'none')
plot(ax1, a, a <= A, '-', 'LineWidth', 2, 'color', colOrd(1, :))
plot(ax1, a, S, '--', 'LineWidth', 2, 'Color', colOrd(1, :))
legend('susceptible', 'immune')
ylabel('susceptible fraction')
xlabel('age (years)')

h = figure;
h.Position = [ 808   377   930   338];
tl = tiledlayout(1, 2);
nexttile;
hold on
fill([0, A, A, 0], [0, 0, 1, 1], blue)
fill([A, L, L, A], [0, 0, 1, 1], grey)
legend('susceptible', 'immune')
ylabel('susceptible fraction')
xlabel('age (years)')
title('(a) without immunusation at birth')

nexttile;
hold on
fill([0, Av, Av, 0], [p, p, 1, 1], blue)
fill([0, Av, Av, 0], [0, 0, p, p], grey)
fill([Av, L, L, Av], [p, p, 1, 1], grey)
fill([Av, L, L, Av], [0, 0, p, p], grey)
legend('susceptible', 'immune')
ylabel('susceptible fraction')
xlabel('age (years)')
title('(b) with immunusation at birth')

