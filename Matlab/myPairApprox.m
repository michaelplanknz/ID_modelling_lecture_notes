function dydt = myPairApprox(t, y, Beta, Gamma, kMean)

S = y(1);
I = y(2);
SI = y(3);
SS = y(4);

dSdt = -Beta*SI;
dIdt = Beta*SI - Gamma*I;
dSIdt = Beta*(kMean-1)/kMean * (SS-SI)*SI/S - Beta*SI - Gamma*SI;
dSSdt = -2*Beta*(kMean-1)/kMean * SI*SS/S;


dydt = [dSdt, dIdt, dSIdt, dSSdt]';