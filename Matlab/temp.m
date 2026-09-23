clear 
close all

N = 1e6;
R0 = 2;
knb = 0.1;

Si = zeros(N, 1);

Z = nbinrnd(knb, knb/(knb+R0), 1, N);
inz = find(Z > 0);
nnz = length(inz)
% Count number of "ancestors" Si of each individual in the pop
for i = 1:nnz
   jCol = inz(i);
   ind = randsample(N, Z(jCol), true);
   Si(ind) = Si(ind) + 1;
end

% Ancestor dist is Poisson even if offspring dist is heavy tailed, provided
% offspring are randomly selected 
histogram(Si, 'Normalization', 'probability');
hold on
plot(0:12, poisspdf(0:12, R0), 'o-')
mean(Si)
var(Si)

