function [Ez0,w0,lam0,chi] = holo_wake_pars(n0,r_in,r_out)

SI_consts;
[omega_p,lambda_p, skin_depth] = plasma_parameters(n0);

a = r_in/skin_depth;
b = r_out/skin_depth;

b0 = B0(a,b);
b3 = B3(a,b);
chi = sqrt(2*b3/(2*b3+a*b0));

Ez0 = (SI_e*omega_p^2/(pi*SI_eps0*SI_c^2))*(b0/(a*(2*b3+a*b0)));
w0 = chi*omega_p;
lam0 = lambda_p/chi;

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