function [km,w0,lam0,chi] = holo_wake_scaled(n0,r_in,r_out,order)

SI_consts;
[omega_p,lambda_p, skin_depth] = plasma_parameters(n0);

a = r_in/skin_depth;
b = r_out/skin_depth;

if order == 0
    bess00 = Bij(0,0,a,b);
    bess10 = Bij(1,0,a,b);
    chi = sqrt(2*bess10/(2*bess10-a*bess00));
    km = -2*bess00/(a*(2*bess10-a*bess00));
end

if order == 1
    bess11 = Bij(1,1,a,b);
    bess21 = Bij(2,1,a,b);
    chi = sqrt(2*bess21/(4*bess21-a*bess11));
    km = -4*bess11/(a*(4*bess21-a*bess11));
end

w0 = chi*omega_p;
lam0 = lambda_p/chi;

