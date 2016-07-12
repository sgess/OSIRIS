function [Ez,Ez0,w0,chi] = holo_wake(n0,a,b,zz)

SI_consts;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);

b0 = B0(a,b);
b3 = B3(a,b);
chi = sqrt(2*b3/(2*b3+a*b0));

Ez = (SI_e*omega_p^2/(pi*SI_eps0*SI_c^2))*(b0/(a*(2*b3+a*b0)))*cos(chi*omega_p*zz/SI_c);
Ez0 = (SI_e*omega_p^2/(pi*SI_eps0*SI_c^2))*(b0/(a*(2*b3+a*b0)));
w0 = chi*omega_p;

end

function x = B0(a,b)

b1 = besselk(0,a);
b2 = besseli(0,b);
b3 = besselk(0,b);
b4 = besseli(0,a);

x = b1.*b2 - b3.*b4;

end

function x = B3(a,b)

b1 = besselk(1,a);
b2 = besseli(0,b);
b3 = besselk(0,b);
b4 = besseli(1,a);

x = b1.*b2 + b3.*b4;

end