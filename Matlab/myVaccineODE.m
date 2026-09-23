function dydt = myVaccineODE(t, y, R0, C, N, VE)

S = y(1:2);
I = y(3:4);

FOI = R0 * (C'*(I./N)) .* [1; 1-VE];

dSdt = -FOI.*S;
dIdt = FOI.*S - I;

dydt = [dSdt; dIdt];

