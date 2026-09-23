function dydt = myTwoTypeODE(t, y, C, N, Gamma)

S = y(1:2);
I = y(3:4);

Lambda = (C')*(I./N);

dSdt = -Lambda.*S;
dIdt = Lambda.*S - Gamma*I;

dydt = [dSdt; dIdt];

