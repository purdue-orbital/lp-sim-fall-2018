syms dP Kt A2 A1 Cc A2A1 q2 v2 rho g mdot

v2 = mdot/rho/g/A2*144; %ft/s
q2 = rho/2*v2^2/144;%psi
Kt = (1/Cc)^2;
A1=1.431388153;%inches
dP=281.0677292;%psi
rho = 1.9364;%slug/ft^3
g=32.2;%ft/s^2
mdot=8.335288406/2;%lbf/s


A1A2 = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
cC = [.63 .635 .65 .665 .69 .715 .75 .8 .87 1];

x1 = (0.1:0.01:1)

p=polyfit(A1A2*A1,cC,5);
cC2 = poly2sym(p,A2)

F = (1/cC2)^2== dP/q2;

