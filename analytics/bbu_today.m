%lg

% load SI_constants
SI_consts;

% load plasma parameters
n0 = 8e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
kp = omega_p/SI_c;

% beam parameters
gamma = 39824;
sig_z = 35e-6;
N = 5.34e9;
I = N*SI_e*SI_c/sig_z;

% load geometric parameters
r_in = 240E-6;
r_out = 290E-6;

R_in = kp*r_in;
R_out = kp*r_out;

B11 = Bij(1,1,R_in,R_out);
B21 = Bij(2,1,R_in,R_out);

B_core = B11/(4*B21-R_in*B11);
chi = sqrt(2*B21/(4*B21-R_in*B11))


wp0 = -16*pi*(kp^3/R_in^3)*(B_core/chi);

lg_me = (1/4)*sqrt(gamma*(lambda_p*1e-6)/(N*SI_re*wp0*chi*sig_z))
lg_me_too = 


K1=besselk(1,R_in);
K2=besselk(2,R_in);
kap1 = kp^2*(K1/(R_in*K2))*(1+R_in*K1/(4*K2))^(-1)
kap1_me = kp^2*4*K1/(R_in*(4*K2+R_in*K1))

kap1_meToo = -kp^2*4*B11/(R_in*(4*B21-R_in*B11))


Lg = 2^(-3/2)*(N*SI_re*kap1*sig_z/(gamma*r_in^2))^(-1/2)
Lg0 = 2^(-5/2)*(N*SI_re*kap1*sig_z/(gamma*r_in^2))^(-1/2)