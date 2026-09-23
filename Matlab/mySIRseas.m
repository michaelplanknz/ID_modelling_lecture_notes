function f = mySIRseas(t, y, par)

S = y(1);  I = y(2);

seas = 1 - par.seasonAmp*cos(2*pi*t/365.25);

f = [ -seas*par.R0/par.tI * S*I + 1/par.tw * (1-S-I) - par.eps;
       seas*par.R0/par.tI * S*I - 1/par.tI * I + par.eps];
 