function [R_vec, N_vec] = hollow_channel(param_struct)

n_points = param_struct.hollow.n_points;
type = param_struct.hollow.type;

R_vec = zeros(1,n_points);
N_vec = zeros(1,n_points);

radius = param_struct.hollow.radius/param_struct.plasma.SD;
width = param_struct.hollow.width/param_struct.plasma.SD;
ramp = param_struct.hollow.ramp/param_struct.plasma.SD;

if strcmp(type,'flat')
    
    if n_points ~= 6
        error('More code needed');
    end
    
    R_vec(1) = 0.0;
    N_vec(1) = 0.0;
    R_vec(2) = radius - width/2 - ramp;
    N_vec(2) = 0.0;
    R_vec(3) = radius - width/2;
    N_vec(3) = 1.0;
    R_vec(4) = radius + width/2;
    N_vec(4) = 1.0;
    R_vec(5) = radius + width/2 + ramp;
    N_vec(5) = 0.0;
    R_vec(6) = param_struct.size.Box_R;
    N_vec(6) = 0.0;
    
end

if strcmp(type,'hollow')
    
    if n_points ~= 6
        error('More code needed');
    end
    
    R_vec(1) = 0.0;
    N_vec(1) = 0.0;
    R_vec(2) = radius - ramp;
    N_vec(2) = 0.0;
    R_vec(3) = radius;
    N_vec(3) = 1.0;
    R_vec(4) = radius + width/2;
    N_vec(4) = 1.0;
    R_vec(5) = radius + width/2 + ramp;
    N_vec(5) = 1.0;
    R_vec(6) = param_struct.size.Box_R;
    N_vec(6) = 1.0;
    
end

if strcmp(type,'cdf')
    
    R_vec = linspace(0,param_struct.size.Box_R,n_points);
    dR_vec = R_vec(2)-R_vec(1);
    gauss = (sqrt(2*pi)*width)^(-1)*exp(-(R_vec-radius).^2/(2*width^2));
    N_vals = 1000*dR_vec*cumsum(gauss);
    N_int_vals = floor(N_vals);
    N_vec = N_int_vals/1000;
    N_vec(end) = 0;
    
end