function [EZ_out, rho_b, min_ind] = CompareTheory(zz,EZ_line,n0,N,charge,sz,cent,a,b,shift_ind)

ax_ind = numel(zz);
[~,ind] = min(abs(zz));
zz_cos = zz(ind:end)/1e6;
zz_rho = zz/1e6;
dz = zz_cos(2)-zz_cos(1);
mu = cent/1e6;
sigz = sz/1e6;

rho_b = N/(sqrt(2*pi)*sigz)*exp(-(zz_rho-mu).^2/(2*sigz^2));

Ez = holo_wake(n0,a,b,zz_cos);
Wz = charge*dz*conv(rho_b,Ez,'full')/1e9;

diff_ind = numel(Wz)-ax_ind;
diff_sum = zeros(1,diff_ind);
for i = 1:diff_ind
    
    ind1 = i;
    ind2 = ax_ind+(i-1);
    Wz_comp = Wz(ind1:ind2);
    
    diff_sum(i) = sum((Wz_comp-EZ_line).^2);
end

if nargin < 10
    [~,min_ind] = min(diff_sum);
else
    min_ind = shift_ind;
end

EZ_out = Wz(min_ind:(ax_ind+(min_ind-1)));

